#include <algorithm>
#include <cassert>
#include <jla_boost/graph/PairToRangeAdaptor.hpp>
#include <mod/Chem.hpp>
#include <mod/graph/Graph.hpp>
#include <mod/graph/internal/Graph.hpp>
#include <mod/rule/Rule.hpp>
#include <mod/rule/internal/Rule.hpp>
#include <mod/lib/Chem/MoleculeUtil.hpp>
#include <mod/lib/Graph/Single.hpp>
#include <mod/lib/Graph/Properties/Molecule.hpp>
#include <mod/lib/Graph/LabelledGraph.hpp>
#include <mod/lib/LabelledUnionGraph.hpp>
#include <mod/lib/Rules/Real.hpp>

#include <boost/bimap.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

#include <map>
#include <optional>
#include <set>
#include <iostream>
#include <fstream>

// import asRange which wraps a std::pair of iterators with a proxy class
// to enable ADL of begin() and end(), so it can be used with ranged-based for loops
using mod::asRange;

using GraphType = mod::lib::Graph::GraphType;
using Vertex = mod::lib::Graph::Vertex;
using Edge = mod::lib::Graph::Edge;
using LUG = mod::lib::LabelledUnionGraph<mod::lib::Graph::LabelledGraph>;
using UVertex = boost::graph_traits<LUG::GraphType>::vertex_descriptor;
using UEdge = boost::graph_traits<LUG::GraphType>::edge_descriptor;
using VertexMap = boost::bimap<UVertex, UVertex>;
using AtomData = mod::AtomData;
using BondType = mod::BondType;
using PropString = mod::lib::Graph::LabelledGraph::PropStringType;
using MolView = mod::lib::Graph::PropMolecule;

using mod::lib::Chem::bondToChar;

std::shared_ptr<mod::rule::Rule> createRule(const VertexMap &vertexMap,
                                            const LUG &lgEduct, const LUG &lgProduct, bool doChemistryCheck);

template<typename Graph>
std::size_t getVertexId(const typename boost::graph_traits<Graph>::vertex_descriptor &v, const Graph &g) {
	return get(boost::vertex_index_t(), g, v);
}

template<typename Graph>
auto getVertexFromId(std::size_t id, const Graph &g) {
	assert(id < num_vertices(g));
	return vertex(id, g);
}

int bondValue(BondType bt) {
	switch(bt) {
	case BondType::Single:
		return 1;
	case BondType::Double:
		return 2;
	case BondType::Triple:
		return 3;
	case BondType::Aromatic:
		throw mod::InputError("Aromatic bonds are not supported.");
	}
	__builtin_unreachable();
}

template< typename T >
struct Comb {
    Comb( std::vector<T> data ) 
    : data(data) 
    {
        for ( int i = 0; i < data.size(); ++i ) {
            chosen.push_back(false);
        }
    }

    Comb() {}

    std::vector<T> &get_data() {
        return data;
    }

    std::vector<bool> &get_chosen() {
        return chosen;
    }

    unsigned int count_chosen() const {
        unsigned int c = 0;

        for ( auto it = chosen.rbegin(); it != chosen.rend() && *it; ++it ) {
            c++;
        }

        return c;
    }

    void push_back( T e ) {
        data.push_back( e );
        chosen.push_back( false );
    }

    std::vector<T> &current() {
        if ( data_image ) {
            return *data_image;
        }
        data_image = std::vector<T>{};
        auto it_d = data.begin();
        auto it_c = chosen.begin();
        while ( it_d != data.end() ) {
            if ( *it_c ) {
                data_image->push_back(*it_d);
            }
            it_d++;
            it_c++;
        }
        return *data_image;
    }

    bool next() {
        if ( false == std::next_permutation(chosen.begin(), chosen.end()) ) {
            return false;
        }
        data_image = std::nullopt;
        return true;
    }

    void reset_chosen( unsigned int n ) {
        std::fill(chosen.begin(), chosen.end()-n, false);
        std::fill(chosen.end()-n, chosen.end(), true);
        data_image = std::nullopt;
    }

private:
    std::optional<std::vector<T>> data_image;
    std::vector<T> data;
    std::vector<bool> chosen;
};

template<class SomeGraph>
struct Generator {

    Generator( const SomeGraph &gEduct, const SomeGraph &gProduct, const LUG &lgEduct, const LUG &lgProduct, unsigned int cycle_length )
        : gEduct(gEduct)
        , gProduct(gProduct)
        , lgEduct(lgEduct)
        , lgProduct(lgProduct)
        , cycle_length(cycle_length)
        {
            // We are assuming that we have the same atoms on both sides.
            for ( UVertex v_product : asRange(vertices(gProduct)) ) {
                AtomData molecule = mod::graph::internal::getMolecule(v_product, lgProduct);
                if ( product_molecules.find(molecule) == product_molecules.end() ) {
                    product_molecules.insert({molecule, Comb<UVertex>({v_product})});
                } else {
                    product_molecules[molecule].push_back(v_product);
                }
            }
            
            for ( std::pair<const AtomData, Comb<UVertex>> &pair : product_molecules ) {
                std::sort(pair.second.get_data().begin(), pair.second.get_data().end());
            }
            
            const auto vsEduct = vertices(gEduct);
			std::vector<UVertex> eductVertices(vsEduct.first, vsEduct.second);
            std::sort(eductVertices.begin(), eductVertices.end(), [&lgEduct](UVertex l, UVertex r){
                return mod::graph::internal::getMolecule(l, lgEduct).getAtomId() < mod::graph::internal::getMolecule(r, lgEduct).getAtomId();
            });
            for ( auto v : eductVertices ) {
                chosen_molecules.push_back(v);
            }
            chosen_molecules.reset_chosen(cycle_length);
            {
                // First count how many we have picked of each type.
                std::map<AtomData, unsigned int> molecule_count;
                for ( auto it = chosen_molecules.current().begin(); it != chosen_molecules.current().end(); it++ ) {
                    auto mol = mod::graph::internal::getMolecule(*it, lgEduct);
                    if ( molecule_count.find(mol) == molecule_count.end() ) {
                        molecule_count[mol] = 1;
                    } else {
                        molecule_count[mol]++;
                    }
                    if ( educt_molecules.find(mol) == educt_molecules.end() ) {
                        educt_molecules[mol] = {*it};
                    } else {
                        educt_molecules[mol].push_back(*it);
                    }
                }
                // Specify that we want that many of each type.
                for ( auto &[mol, c] : molecule_count ) {
                    product_molecules[mol].reset_chosen(c);
                    mapping_order.push_back(mol);
                }
            }
            // Now we have seperated the vertices into containers of unique type.
            // And we have selected cycle_length amount of vertices from the educt side to map.

            // Sort so we place the Hydrogens last.
            make_mapping();
        }

    bool operator++() {
        bool permute_next_molecule;
        int current_molecule = 0;
        // std::cout << "Mapping order: ";
        // for ( auto &atom : mapping_order ) {
        //     std::cout << atom << " ";
        // }
        // std::cout << "\n";
        do {
            permute_next_molecule = false;
            AtomData molecule = mapping_order[current_molecule];
            std::vector<UVertex> &educt = educt_molecules[molecule];
            Comb<UVertex> &product = product_molecules[molecule];
            // std::cout << "Current molecule index " << current_molecule << " With molecule " << molecule << "\n";
            
            // std::cout << "Product: ";
            // for ( auto v : product.current() ) { std::cout << v; }
            // std::cout << "\nEduct:   ";
            // for ( auto v : educt ) { std::cout << v; }
            // std::cout << "\n";
            assert(product.current().size() == educt.size());

            // std::cout << "Before: ";
            // for ( auto v : product.current() ) {
            //     std::cout << v << ' ' << mod::graph::internal::getMolecule(v, lgProduct) << ' ';
            // }
            // std::cout << "\nAfter:  "; 

            
            if ( product.next() == false ) {
                // Progress to the next molecule
                // std::cout << "We are false?";
                permute_next_molecule = true;
                current_molecule = (current_molecule+1) % mapping_order.size();
                if (current_molecule == 0) {
                    // Select new educt molecules.
                    if ( false == permute_educts() ) {
                        // There are no more new combinations of molecules.
                        std::cout << "No more permutations\n";
                        return false;
                    }
                }
            }
            // for ( auto v : product_molecules[molecule].current() ) {
            //     std::cout << v << ' ' << mod::graph::internal::getMolecule(v, lgProduct) << ' ';
            // }
            // std::cout << std::endl;
            
            make_mapping();

        } while (permute_next_molecule);
        return true;
    }

    VertexMap operator*() const {
        // std::cout << "Size of VM: " << vm.size() << std::endl;
        return this->vm;
    }
private:

    bool permute_educts() {
        // std::cout << "Call to permute educts!!!!!!!!!!!!!!!!!!!!!\n";
        if ( !chosen_molecules.next() ) {
            return false;
        }
        mapping_order.clear();
        educt_molecules.clear();

        // First count how many we have picked of each type.
        std::map<AtomData, unsigned int> molecule_count;
        for ( auto it = chosen_molecules.current().begin(); it != chosen_molecules.current().end(); it++ ) {
            auto mol = mod::graph::internal::getMolecule(*it, lgEduct);
            // std::cout << *it << " " << mol << "\n";
            if ( molecule_count.find(mol) == molecule_count.end() ) {
                molecule_count[mol] = 1;
            } else {
                molecule_count[mol]++;
            }
            if ( educt_molecules.find(mol) == educt_molecules.end() ) {
                educt_molecules[mol] = {*it};
            } else {
                educt_molecules[mol].push_back(*it);
            }
        }
        // Specify that we want that many of each type.
        for ( auto &[mol, c] : molecule_count ) {
            // std::cout << "We want " << c << " of " << mol << std::endl;
            product_molecules[mol].reset_chosen(c);
            mapping_order.push_back(mol);
        }
        std::sort(mapping_order.begin(), mapping_order.end(), [](auto &l, auto &r) {return l.getAtomId() > r.getAtomId(); });
           

        return true;
    }

    void make_mapping() {
        vm.clear();
        std::sort(mapping_order.begin(), mapping_order.end(), [](auto &l, auto &r) { return l.getAtomId() > r.getAtomId(); });
        for ( AtomData &molecule : mapping_order ) {                
            auto it_educt = educt_molecules[molecule].begin();
            auto it_product = product_molecules[molecule].current().begin();
            while ( it_educt != educt_molecules[molecule].end() ) {
                vm.insert({*it_educt, *it_product});
                it_educt++;
                it_product++;
            }
        }
    }

    VertexMap vm;
    const SomeGraph &gEduct;
    const SomeGraph &gProduct;
    const LUG &lgEduct;
    const LUG &lgProduct;
    unsigned int cycle_length;
    
    Comb<UVertex> chosen_molecules;

    std::vector<AtomData> mapping_order;
    std::map<AtomData, Comb<UVertex>> product_molecules;
    std::map<AtomData, std::vector<UVertex>> educt_molecules;
};

template <class Vertex>
struct cycle_detector : public boost::dfs_visitor<> {
    cycle_detector( bool& has_cycle, std::vector<std::vector<Vertex>> &cycles )
            : _has_cycle(has_cycle)
            , _cycles(cycles) {}

    template<class Edge, class Graph>
    void back_edge(Edge e, Graph &g) {
        auto source = boost::source(e, g);
        auto target = boost::target(e, g);

        if (parent.find(target) != parent.end() && parent[target] == source) {
            return;
        }
//        std::cout << "--[[ Cycle detected ]]--\n";

        _has_cycle = true;

        std::vector<Vertex> _cycle;
        _cycle.push_back(source);

//        std::cout << "Adding vertex to cycle " << source << "\n";
        for (auto current = target; current != source; current = parent[current]) {
            _cycle.push_back(current);
//            std::cout << "Adding vertex to cycle " << current << "\n";
        }
        _cycles.push_back(_cycle);
    }

    template<class Edge, class Graph>
    void examine_edge( Edge e, Graph &g) {
        auto source = boost::source(e, g);
        auto target = boost::target(e, g);
        if (parent.find(target) != parent.end() && parent[target] == source) {
            return;
        }
        parent[source] = target;
    }


protected:
    bool &_has_cycle;
    std::vector<std::vector<Vertex>> &_cycles;
    std::map<Vertex, Vertex> parent;
};

using CycleGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_name_t, std::string>>;
using CycleVDescriptor = boost::graph_traits<CycleGraph>::vertex_descriptor;

/**
 * This functions takes a morphism and two graphs and does an XOR on all the edges.
 * The resulting graph will be used to detect cycles between the two.
 * @tparam SomeGraph Because of the auto keywords everywhere it is quite bothersome to figure out the correct type...
 * @param vm The morphism mapping the educt to the product
 * @param graph_map The mapping to and from the returned CycleGraph and the educt
 * @param gEduct The graph representation of the Educt
 * @param gProduct The graph representation of the Product
 * @param lgEduct The labelled graph representation of the Educt
 * @param lgProduct The labelled graph representation of the Product
 * @return An XOR graph from the educt and product graphs where each edge is labelled with either "ADD" or "REMOVE" depending on if it adds an edge or removes an edge from the Educt to the Product.
 */
template< class SomeGraph >
CycleGraph XOR_graphs( VertexMap &vm, boost::bimap<UVertex, CycleVDescriptor> &graph_map, const SomeGraph &gEduct, const SomeGraph &gProduct, const LUG &lgEduct, const LUG &lgProduct ) {
    /**
     * Implementation details:
     *      Generally there are three cases where we want to XOR the two graphs:
     *      - [1] There exists an edge on the educt side where the morphism maps to a corresponding edge on the product side
     *           a) They share label -> Don't add it
     *           b) They do NOT share the same label -> Add it
     *      - [2] The edge is unique on the educt side -> Add it
     *      - [3] The edge is unique on the product side -> Add it
     *
     *      Because the mød graphs are undirected graphs, when traversing all the edges,
     *      the edge will appear twice in the forms (v, u) and (u, v). To avoid doubling
     *      the amount of edges we do a quick check to see if the edge already exists in
     *      the XOR graph.
     *
     *      We are assuming that the only valid labels for edges are '-' and '='.
     */
    CycleGraph cycle_graph;
    auto edge_types = boost::get(boost::edge_name, cycle_graph);



    // This lambda function returns a vertex descriptor for the corresponding input vertex. If it has already been inserted into the cycle_graph it will return that, otherwise return a new one and insert it into the map.
    // We are only mapping the left graph to make it simple. In order to call this with a right mapping, simply use the vertex map's right side to map it into the left graph.
    auto getVertex = [&graph_map, &cycle_graph]( UVertex vertex_to_find ) -> CycleVDescriptor {
        CycleVDescriptor v;
        auto it = graph_map.left.find(vertex_to_find);
        if ( it != graph_map.left.end() ) {
            v = it->second;
        } else {
            v = boost::add_vertex(cycle_graph);
            graph_map.insert({vertex_to_find, v});
        }
        return v;
    };

    for ( const auto &vertex_pair: vm ) {
        // Go through all the neighbors of the left hand side, and find the matching right hand side from the vertex_map. Check if this edge also exists on the right side, if so, remove.

        const UVertex left_source = vertex_pair.left;
        const UVertex right_source = vertex_pair.right;

        for ( const UEdge &left_edge : asRange(out_edges(left_source, gEduct)) ) {
            
            const UVertex left_target = target(left_edge, gEduct);

            // Do we actually have the target in the VertexMap. If not -> continues;
            if ( vm.left.find(left_target) == vm.left.end() ) { continue; }
            
            // Find the right side
            const UVertex right_target = vm.left.at(left_target);
            
            // Do we have a matching right edge?
            const std::pair<UEdge, bool> &right_edge = edge(right_source, right_target, gProduct);

            if ( right_edge.second ) {
                // We found a matching edge.
                // At this point, we have to check the difference between the left side and the right side bond. If they differ, we want to add this edge.
                const std::string &left_bond = mod::graph::internal::getString(left_edge, lgEduct);
                const std::string &right_bond = mod::graph::internal::getString(right_edge.first, lgProduct);
                if ( left_bond != right_bond ) {
                    // Case [1b] Add the edge. There is a difference.

                    CycleVDescriptor from = getVertex(left_source);
                    CycleVDescriptor to = getVertex(left_target);
                    if ( boost::edge(from, to, cycle_graph).second == false ) {
                        
                        // Assuming only '-' and '=' are possible.
                        if (left_bond == "-" && right_bond == "=") {
                            auto e = boost::add_edge(from, to, cycle_graph).first;
                            edge_types[e] = "ADD";
                        } else if (left_bond == "=" && right_bond == "-"){
                            auto e = boost::add_edge(from, to, cycle_graph).first;
                            edge_types[e] = "REMOVE";
                        }
                    }
                }
            } else {
                // Case [2]
                // There exists an edge on the left hand side but NOT on the right hand side, therefore we want to add this edge.
                const std::string &left_bond = mod::graph::internal::getString(left_edge, lgEduct);

                CycleVDescriptor from = getVertex(left_source);
                CycleVDescriptor to = getVertex(left_target);
                if ( boost::edge(from, to, cycle_graph).second == false ) {
                    
                    if ( left_bond == "-" ) {
                        auto e = boost::add_edge(from, to, cycle_graph).first;
                        edge_types[e] = "REMOVE";
                    }
                    
                }
            }
        }
        // At this point, do the same but for the right side. However, we only want to add the edges that only appear on the right side.
        for ( const UEdge &right_edge : asRange(out_edges(right_source, gProduct)) ) {
            // Find the left side
            const UVertex right_target = target(right_edge, gProduct);

            if ( vm.right.find(right_target) == vm.right.end() ) { continue; }
            const UVertex left_target = vm.right.at(right_target);

            // Do we have a matching left edge?
            const std::pair<UEdge, bool> &left_edge = edge(left_source, left_target, gEduct);

            if ( !left_edge.second ) {
                // Case [3]
                // We want to add this edge.
                const std::string &right_bond = mod::graph::internal::getString(right_edge, lgProduct);
                CycleVDescriptor from = getVertex(left_source);
                CycleVDescriptor to = getVertex(left_target);
                if (boost::edge(from, to, cycle_graph).second == false) {
                    if ( right_bond == "-" ) {
                        auto e = boost::add_edge(from, to, cycle_graph).first;
                        edge_types[e] = "ADD";
                    }
                }

            }
        }
    }
    return cycle_graph;
}

std::vector<std::shared_ptr<mod::rule::Rule> > doStuff(
		const std::vector<std::shared_ptr<mod::graph::Graph>> &educts,
		const std::vector<std::shared_ptr<mod::graph::Graph>> &products,
		bool doChemistryCheck,
        int k,
        int c) {

	// first make objects representing the disjoint union of the educts and products
	const auto makeUnion = [](const auto &gs) {
		LUG lug;
		for(const auto &g: gs)
			mod::graph::internal::push_back(lug, &mod::graph::internal::getLabelledGraph(g->getGraph()));
		return lug;
	};
	const auto lgEduct = makeUnion(educts);
	const auto lgProduct = makeUnion(products);
	// the underlying graphs of the unions
	const auto &gEduct = get_graph(lgEduct);
    const auto &gProduct = get_graph(lgProduct);

	if(num_vertices(gEduct) != num_vertices(gProduct)) {
		throw mod::InputError("Error: the graphs do not have the same number of vertices:\n\teduct: "
		                      + std::to_string(num_vertices(gEduct)) + "\tproduct: "
		                      + std::to_string(num_vertices(gProduct)) + "\n");
	}

	std::vector<VertexMap> vertexMaps;
    CycleGraph cycle_graph;
	{ // START WORKING HERE
		{ // Example THREE --- an actual enumeration algorithm
			// In essence, we iterate over all permutations of the vertices in the educt graph,
			// and for each of them we make a mapping to the vertices of the product graph.
			// Note, this does no 'alternating addition/removal' check for the edges, so you will very likely
			// get non-chemical reactions. This will in essence then be your main task:
			// make sure you will have only chemically valid rules as described on the webpage.
			// While it can be done by skipping permutations, a fast solution would probably not rely on simple
			// permutation generation.
			// For now, we just disable the chemistry check as a hax to get something printed:

			// Copy the educt vertices into a vector so we can permute them.
			const auto vsEduct = vertices(gEduct);
			std::vector<UVertex> eductVertices(vsEduct.first, vsEduct.second);
			const auto productVertices = asRange(vertices(gEduct));

			// std::next_permutation uses sortedness to determine the last permuation, so start by sorting.
			std::sort(eductVertices.begin(), eductVertices.end());

			constexpr int limit = 10;
			// We will limit to 'limit' permutations, as otherwise the PDF might get too large.
			// Note that we also reject all rules where the atom type would change.
			int permutationCount = 0;
            int num = 0;
            std::cout << "0" << std::endl;
            Generator vertexmap_gen(gEduct, gProduct, lgEduct, lgProduct, c);
            bool repeat = false;
			do {
                std::cout << "1" << std::endl;
                VertexMap vm = *vertexmap_gen;
                std::cout << "2" << std::endl;
                for(const auto &vt: vm) {
                    std::cout << "\"" << mod::graph::internal::getMolecule(vt.left, lgEduct) << "\" "
                            << getVertexId(vt.left, gEduct)
                            << "\t<=>   "
                            << getVertexId(vt.right, gProduct) << " \""
                            << mod::graph::internal::getMolecule(vt.right, lgProduct) << "\"\n";
                }
                std::cout << std::endl;
                // std::string a;
                // std::cin >> a;
                
                repeat = ++vertexmap_gen;
                bool valid = true;
                
                // Make sure to copy the VertexMap before modifying it.

//                 for (int i = 0; i != eductVertices.size(); i++) {
//                     // check if the molecule data for the atoms matches, otherwise we can't make a rule
// //                    std::cout << "coutIteration " << i << '\n';
//                     std::cout << mod::graph::internal::getMolecule(eductVertices[i], lgEduct) << " I " << mod::graph::internal::getMolecule(productVertices[i], lgProduct) << "\n";
//                     if (mod::graph::internal::getMolecule(eductVertices[i], lgEduct) !=
//                         mod::graph::internal::getMolecule(productVertices[i], lgProduct)) {
//                         valid = false;
//                         break;
//                     } else {
//                         std::cout << "Valid mapping\n";
//                     }
//                     vm.insert(VertexMap::value_type(eductVertices[i], productVertices[i]));
//                 }
                // if (!valid) continue;

                boost::bimap <UVertex, CycleVDescriptor> cycle_graph_map;
                cycle_graph.clear();
                cycle_graph = XOR_graphs(vm, cycle_graph_map, gEduct, gProduct, lgEduct, lgProduct);
                auto edge_types = boost::get(boost::edge_name, cycle_graph);
                // std::cout << "Size of XOR graph: " << boost::num_vertices(cycle_graph) << "\n";

                bool has_cycle = false;
                std::vector <std::vector<CycleVDescriptor>> cycles;
                cycle_detector vis(has_cycle, cycles);
                boost::depth_first_search(cycle_graph, visitor(vis));

                if (!has_cycle) {
                    std::cout << "No Cycle :(" << std::endl;
                    continue;
                }

                std::cout << "Found cycle" << std::endl;

                // From this point, we know we have a cycle.

                if (cycles.size() != 1) {
                    std::cout << "Found " << cycles.size() << " cycles" << std::endl;
                    // for ( auto &cycle : cycles ) {
                    //     for (auto it = cycle.begin(); it != cycle.end() && valid; it++) {
                    //         std::cout << *it;
                    //         if (it != cycle.begin()) {
                    //             std::cout << ", ";
                    //         }
                    //     }
                    //     std::cout << "\nPrint cycle? ";

                    //     std::string answer = "";
                    //     std::cin >> answer;
                    //     if ( answer == "y" ) {
                    //         std::string filename = "graph" + std::to_string(permutationCount) + ".dot";
                    //         std::ofstream dotFile(filename);
                    //         boost::write_graphviz(dotFile, cycle_graph);
                    //         dotFile.close();
                    //         auto r = createRule(vm, lgEduct, lgProduct, true);
                    //         r->print();
                    //         std::cin >> answer;
                    //     }
                    // }
                    continue;
                }
                
                std::vector <CycleVDescriptor> cycle = cycles.at(0);

                // Check if it is a valid cycle.
                if (cycle.size() % 2 == 1)
                    continue; // We know that we can't possibly alternate between an add and a remove.
                std::cout << "With even length" << std::endl;
                // If the cycle is not of the desired length, discard it.
                if (cycle.size() != c)
                    continue;
                std::cout << "And with the correct cycle length of " << c << std::endl;


                // Cycle is only valid if we are alternating between a deleted and an add.
                std::string prev_operation = "";
                bool first = true;
                CycleVDescriptor prev;
                for (auto it = cycle.begin(); it != cycle.end() && valid; it++) {
                    const CycleVDescriptor v = *it;
                    if (!first) {
                        auto e = boost::edge(prev, v, cycle_graph).first;
                        // Check if we are alternating between ADD and REMOVE
                        std::string this_operation = edge_types[e];
                        if (prev_operation == "") {
                            prev_operation = this_operation;
                        } else {
                            if (prev_operation == this_operation) {
                                valid = false;
                            }
                            prev_operation = this_operation;
                        }
                    }
                    first = false;
                    prev = v;
                }
                if (!valid) continue;
                
                std::cout << "The cycle is correctly alternating" << std::endl;

                // We can also use the XOR graph to detect if more than 2 changes per vertex are happening. We can do this because each edge represents a change in the production.
                // If so, discard that solution.
                auto range = boost::vertices(cycle_graph);
                for (auto v = range.first; v != range.second; v++) {
                    if (boost::degree(*v, cycle_graph) > 2) {
                        valid = false;
                    }
                }

                if (!valid) continue;
                std::cout << "The degree for all vertices are less than 3 and the mapping is VALID!\n";
                // We know there is only a single cycle.

                

                // Now that we have ensured that the graph is chemically valid, we want to remove the vertex-mappings where the shortest distance to the cycle is greater than k.
                
                // We can go k out, for each vertex of the cycle, and add these vertices to the mapping.

                std::vector<std::pair<UVertex, UVertex>> work_set;
                unsigned int current_distance_from_cycle = 1;
                for ( auto &vt : vm ) {
                    work_set.push_back({vt.left, vt.right});
                }
                while ( current_distance_from_cycle <= k ) {
                    std::vector<std::pair<UVertex, UVertex>> subset;
                    
                    for ( auto &[left_source, right_source] : work_set ) {

                         for ( const auto &left_edge : asRange(out_edges(left_source, gEduct)) ) {
                            UVertex left_target = target(left_edge, gEduct);
                            // If we already have this in the mapping, skip
                            if ( vm.left.find(left_target) != vm.left.end() ) continue;
                            AtomData left_mol = mod::graph::internal::getMolecule(left_target, lgEduct);

                            // Find a matching pair on the right side.
                            for ( const auto &right_edge : asRange( out_edges(right_source, gProduct)) ) {
                                UVertex right_target = target(right_edge, gProduct);
                                // If we already have this mapped, skip
                                if ( vm.right.find( right_target ) != vm.right.end() ) continue;
                                AtomData right_mol = mod::graph::internal::getMolecule(right_target, lgProduct);

                                // And if the atoms don't match, skip
                                if ( left_mol != right_mol ) continue;
                                // Check if the bond is the same
                                if ( mod::graph::internal::getMolecule(left_edge, lgEduct) != mod::graph::internal::getMolecule(right_edge, lgProduct) ) continue;
                                std::cout << "Adding context " << left_mol << std::endl;
                                // std::string a;
                                // std::cin >> a;
                                // It is a match. We can potentially have multiple matches here.
                                // But it might not be that interesting to get the same cycle with a different context.
                                // A chemist would know.

                                // We need to check if the new mapping will make either the left or the right side have more than two edges changing for each vertex. If it does, we have to discard it.

                                // If a potential vertex is being added, we must ensure that for both sides no unique edges connect it to the cycle.
                                bool valid = true;
                                for ( auto &[l, r] : work_set ) {
                                    if ( edge(l, left_target, gEduct).second != edge(r, right_target, gProduct).second ) {
                                        valid = false;
                                        break;
                                    }
                                }
                                
                                if ( valid ) {
                                    vm.insert({left_target, right_target});
                                    subset.push_back({left_target, right_target});
                                }
                            }

                        }

                    }

                    work_set = subset;
                    current_distance_from_cycle++;
                }


                // Or we can:
                
                // To do this we do a johnson_all_pairs_shortest_paths which will give us a distance matrix. We can iterate over all the vertices in the cycle, and add the vertices which distances are less than or equal to k.
                // For this we need the entire graph which means Educt OR Product, so we will need to do an OR on all the edges.

                // Make a copy of the graph and give each edge a weight of 1.
#if 0
                using DstGraph = boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,
                /* Vertex prop */ boost::no_property,
                /* Edge prop   */ boost::property<boost::edge_weight_t, int>>;
                using DstVDescriptor = boost::graph_traits<DstGraph>::vertex_descriptor;

                DstGraph dst_graph;
                auto edge_weight = boost::get(boost::edge_weight, dst_graph);
                boost::bimap <UVertex, DstVDescriptor> dst_graph_map;
                for (auto e: asRange(edges(gEduct))) {
                     const UVertex left_s = source(e, gEduct);
                     const UVertex left_t = target(e, gEduct);

                     DstVDescriptor right_s;
                     DstVDescriptor right_t;

                     auto s_it = dst_graph_map.left.find(left_s);
                     auto t_it = dst_graph_map.left.find(left_t);
                     if (s_it != dst_graph_map.left.end()) {
                        right_s = s_it->second;
                     } else {
                        right_s = add_vertex(dst_graph);
                        dst_graph_map.insert({left_s, right_s});
                     }
                     if (t_it != dst_graph_map.left.end()) {
                        right_t = t_it->second;
                     } else {
                        right_t = add_vertex(dst_graph);
                        dst_graph_map.insert({left_t, right_t});
                     }
                     if (boost::edge(right_s, right_t, dst_graph).second == false) {
                        auto e = boost::add_edge(right_s, right_t, dst_graph).first;
                        edge_weight[e] = 1;
                     }
                }
                for (auto e: asRange(edges(gProduct))) {
                    // We don't care about the Product side, so we map them to the educt side.
                    if ( vm.right.find( source(e, gProduct) ) == vm.right.end() ) continue;
                    if ( vm.right.find( target(e, gProduct) ) == vm.right.end() ) continue;
                    const UVertex left_s = vm.right.at(source(e, gProduct));
                    const UVertex left_t = vm.right.at(target(e, gProduct));

                    DstVDescriptor right_s;
                    DstVDescriptor right_t;

                    auto s_it = dst_graph_map.left.find(left_s);
                    auto t_it = dst_graph_map.left.find(left_t);
                    if (s_it != dst_graph_map.left.end()) {
                       right_s = s_it->second;
                    } else {
                        right_s = add_vertex(dst_graph);
                        dst_graph_map.insert({left_s, right_s});
                    }
                    if (t_it != dst_graph_map.left.end()) {
                       right_t = t_it->second;
                    } else {
                       right_t = add_vertex(dst_graph);
                       dst_graph_map.insert({left_t, right_t});
                    }
                    if (boost::edge(right_s, right_t, dst_graph).second == false) {
                       auto e = boost::add_edge(right_s, right_t, dst_graph).first;
                       edge_weight[e] = 1;
                    }
                }
                std::cout << "We have successfully made the OR graph\n";
                unsigned int num_v = boost::num_vertices(dst_graph);
                std::vector <std::vector<int>> distance_matrix(num_v, std::vector<int>(num_v));

                boost::johnson_all_pairs_shortest_paths(dst_graph, distance_matrix);

                // Now that we have a distance matrix, iterate over all the vertices in the cycle, and only add the vertices that are <= k.

                std::set <UVertex> retain_vertices;
                for (const CycleVDescriptor cycle_v : cycle) {
                    // Map the vertex from the cycle graph -> Educt graph -> Dst Graph
                    DstVDescriptor from_vertex = dst_graph_map.left.at(cycle_graph_map.right.at(cycle_v));
                    for (UVertex v: asRange(vertices(gEduct))) {
                        // Map the vertex from the Educt graph -> Dst Graph
                        DstVDescriptor to_vertex = dst_graph_map.left.at(v);
                        if (distance_matrix[from_vertex][to_vertex] <= k) {
                            // Add the vertex v
                            retain_vertices.insert(v);
                        }
                    }
                }
                // Only retain the mappings which are in the set
                auto it = vm.begin();
                while (it != vm.end()) {
                    if (retain_vertices.find(it->left) == retain_vertices.end()) {
                        it = vm.erase(it);
                    } else {
                        ++it;
                    }
                }
#endif
#if 0
                std::string filename = "graph" + std::to_string(permutationCount) + ".dot";
                std::ofstream dotFile(filename);
                boost::write_graphviz(dotFile, cycle_graph);
                dotFile.close();
#endif
                std::cout << "We found a good mapping!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                vertexMaps.push_back(std::move(vm));
				++permutationCount;
				if(permutationCount == limit) break;
			} while(repeat);
		} // end of Example THREE
	} // END WORKING (roughly) HERE.

	{ // debug stuff
		std::cout << "Maps:\n";
		for(const auto &vm: vertexMaps) {
			for(const auto &vt: vm) {
				std::cout << "\"" << mod::graph::internal::getMolecule(vt.left, lgEduct) << "\" "
				          << getVertexId(vt.left, gEduct)
				          << "\t<=>   "
				          << getVertexId(vt.right, gProduct) << " \""
				          << mod::graph::internal::getMolecule(vt.right, lgProduct) << "\"\n";
			}
			std::cout << '\n';
		}
	}

	std::vector<std::shared_ptr<mod::rule::Rule>> rules;
	for(const VertexMap &vm: vertexMaps) {
		auto r = createRule(vm, lgEduct, lgProduct, doChemistryCheck);
		// use r->setName(name); to give the rule a nicer name
		const auto iter = std::find_if(
				rules.begin(), rules.end(),
				[&r](const auto &rOther) {
					auto labelSettings = mod::LabelSettings(mod::LabelType::String,
					                                        mod::LabelRelation::Isomorphism);
					return 1 == r->isomorphism(rOther, 1, labelSettings);
				});
		if(iter != rules.end()) {
			std::cout << "Duplicate rule deleted\n";
		} else {
			rules.push_back(r);
		}
	}
	return rules;
}

#ifdef AS_PYTHON_EXTENSION
#include <boost/python.hpp>
#include <mod/Config.hpp>

namespace py = boost::python;

namespace {
// this can be used to make sure the extension and mod is using the same shared library

uintptr_t magicLibraryValue() {
	return (uintptr_t)&mod::getConfig();
}

} // namespace

BOOST_PYTHON_MODULE(pydoStuff) { // this macro argument is the name of the module, it must be the same as the .so file name.
	py::def("magicLibraryValue", &magicLibraryValue);

	// Change the string in the first argument to give the function another name in Python.
	// The second argument is the function pointer to the function above.
	py::def("doStuff", &doStuff);
}

#else // not AS_PYTHON_EXTENSION

int main(int argc, char **argv) {
	std::vector<std::shared_ptr<mod::graph::Graph> > educts, products;
	{ // do something else, e.g., take SMILES strings from argv
		// std::shared_ptr<mod::graph::Graph> g1, g2;
		// g1 = mod::graph::Graph::fromSMILES("OCC=O");
		// g2 = mod::graph::Graph::fromSMILES("OC=CO");
		// educts.push_back(g1);
		// products.push_back(g2);
        educts.push_back( mod::graph::Graph::fromSMILES("O") );
        educts.push_back( mod::graph::Graph::fromSMILES("Cl") );
        educts.push_back( mod::graph::Graph::fromSMILES("CC(=O)OCC") );
        products.push_back( mod::graph::Graph::fromSMILES("Cl") );
        products.push_back( mod::graph::Graph::fromSMILES("OCC") );
        products.push_back( mod::graph::Graph::fromSMILES("CC(=O)O") );
	}
    constexpr int c = 6;
    constexpr int k = 0;
	auto rules = doStuff(educts, products, true, k, c);
	for(auto r: rules) r->print();
	return 0;
}

#endif // AS_PYTHON_EXTENSION

//------------------------------------------------------------------------------
// Library stuff
//------------------------------------------------------------------------------

#include <mod/graph/Printer.hpp>
#include <mod/rule/Rule.hpp>
#include <mod/lib/Graph/Properties/String.hpp>
#include <mod/lib/Rules/Properties/Molecule.hpp>

std::shared_ptr<mod::rule::Rule> createRule(const VertexMap &vertexMap,
                                            const LUG &lgEduct, const LUG &lgProduct, bool doChemistryCheck) {
	using namespace mod;
	using namespace mod::lib;
	const auto &gEduct = mod::graph::internal::getGraph(lgEduct);
	const auto &gProduct = mod::graph::internal::getGraph(lgProduct);
	const auto n = num_vertices(gEduct);
	{ // error checking
		const auto n2 = num_vertices(gProduct);
		if(n != n2) {
			std::cout << "Different number of vertices: " << n << " and " << n2 << std::endl;
			assert(false);
			std::exit(1);
		}

		for(const auto &vt: vertexMap) {
			const auto vLeft = vt.left;
			const auto vRight = vt.right;
			const auto vLeftId = get(boost::vertex_index_t(), gEduct, vLeft);
			const auto vRightId = get(boost::vertex_index_t(), gProduct, vRight);
			if(vLeftId >= n) {
				std::cout << "Invalid left vertex, index is " << vLeftId << std::endl;
				assert(false);
				std::exit(1);
			}
			if(vRightId >= n) {
				std::cout << "Invalid right vertex, index is " << vRightId << std::endl;
				assert(false);
				std::exit(1);
			}
			const auto &lLabel = mod::graph::internal::getString(vLeft, lgEduct);
			const auto &rLabel = mod::graph::internal::getString(vRight, lgProduct);
			if(lLabel != rLabel) {
				std::cout << "Label mismatch: " << lLabel << " " << vLeftId << " <=> " << vRightId << " " << rLabel
				          << std::endl;
				assert(false);
				std::exit(1);
			}
		}
	} // end of error checking

	using CoreVertex = mod::lib::Rules::Vertex;
	using CoreEdge = mod::lib::Rules::Edge;
	using Membership = mod::lib::Rules::Membership;
	auto dpoRule = mod::rule::internal::makeLabelledRule();
	dpoRule.pString = mod::rule::internal::makePropStringCore(dpoRule);
	auto &pStringCore = *dpoRule.pString;
	auto &core = mod::rule::internal::getGraph(dpoRule);
	std::vector<CoreVertex> leftToCore(n, core.null_vertex()), rightToCore(n, core.null_vertex());

	// copy all matched vertices form the left
	for(const auto v: asRange(vertices(gEduct))) {
		const auto iter = vertexMap.left.find(v);
		if(iter == vertexMap.left.end()) continue;
		const auto vId = get(boost::vertex_index_t(), gEduct, v);
		const auto &label = mod::graph::internal::getString(v, lgEduct);
		const auto vCore = add_vertex(core);
		core[vCore].membership = Membership::K;
		mod::rule::internal::add(pStringCore, vCore, label, label);
		leftToCore[vId] = vCore;
	}

	// link vertices from right side
	for(const auto v: asRange(vertices(gProduct))) {
		const auto iter = vertexMap.right.find(v);
		if(iter == vertexMap.right.end())continue;
		const auto vId = get(boost::vertex_index_t(), gProduct, v);
		const auto left = iter->second;
		const auto leftId = get(boost::vertex_index_t(), gEduct, left);
		const auto vCore = leftToCore[leftId];
		assert(core[vCore].membership == Membership::K);
		rightToCore[vId] = vCore;
	}

	// copy left edges
	for(const auto e: asRange(edges(gEduct))) {
		const auto src = source(e, gEduct);
		const auto tar = target(e, gEduct);
		const auto srcCore = leftToCore[get(boost::vertex_index_t(), gEduct, src)];
		const auto tarCore = leftToCore[get(boost::vertex_index_t(), gEduct, tar)];
		if(srcCore == core.null_vertex()) continue;
		if(tarCore == core.null_vertex()) continue;
		const auto &label = mod::graph::internal::getString(e, lgEduct);
		const auto p = add_edge(srcCore, tarCore, core);
		assert(p.second);
		core[p.first].membership = Membership::L;
		mod::rule::internal::add(pStringCore, p.first, label, "");
	}

	// copy right edges, or promote to context
	for(const auto e: asRange(edges(gProduct))) {
		const auto src = source(e, gProduct);
		const auto tar = target(e, gProduct);
		const auto srcCore = rightToCore[get(boost::vertex_index_t(), gProduct, src)];
		const auto tarCore = rightToCore[get(boost::vertex_index_t(), gProduct, tar)];
		if(srcCore == core.null_vertex()) continue;
		if(tarCore == core.null_vertex()) continue;
		auto p = edge(srcCore, tarCore, core);
		const auto &label = mod::graph::internal::getString(e, lgProduct);
		if(p.second) {
			core[p.first].membership = Membership::K;
			mod::rule::internal::setRight(pStringCore, p.first, label);
		} else {
			p = add_edge(srcCore, tarCore, core);
			assert(p.second);
			core[p.first].membership = Membership::R;
			mod::rule::internal::add(pStringCore, p.first, "", label);
		}
	}

	if(doChemistryCheck) {
		const auto error = [&]() {
			auto r = mod::rule::internal::makeRule(std::move(dpoRule));
			r->print();
			std::cout << "Run 'mod_post' to see rule." << std::endl;
			std::exit(1);
		};

		const auto molView = mod::rule::internal::makePropMoleculeCore(dpoRule, pStringCore);
		for(const auto v: asRange(vertices(core))) {
			unsigned int left = 0, right = 0;
			int bondChange = 0;
			for(const auto e: asRange(out_edges(v, core))) {
				const auto membership = core[e].membership;
				if(membership == Membership::L) {
					left++;
					const auto bt = mod::rule::internal::getMoleculeLeft(e, molView);
					bondChange -= bondValue(bt);
				} else if(membership == Membership::R) {
					right++;
					const auto bt = mod::rule::internal::getMoleculeRight(e, molView);
					bondChange += bondValue(bt);
				} else {
					const auto btLeft = mod::rule::internal::getMoleculeLeft(e, molView);
					const auto btRight = mod::rule::internal::getMoleculeRight(e, molView);
					if(btLeft == btRight) continue;
					++left;
					++right;
					bondChange -= bondValue(btLeft);
					bondChange += bondValue(btRight);
				}
			}
			if(bondChange != 0) {
				std::cout << "Non-zero bond change for vertex " << get(boost::vertex_index_t(), core, v) << ": "
				          << bondChange << "\n";
				std::cout << "Label: (" << mod::rule::internal::getStringLeft(v, pStringCore)
				          << ", " << mod::rule::internal::getStringRight(v, pStringCore) << ")" << std::endl;
				error();
			}
			if(left > 2) {
				std::cout << "Too many left side edges (" << left << ") for vertex "
				          << get(boost::vertex_index_t(), core, v) << ".\n";
				std::cout << "Label: (" << mod::rule::internal::getStringLeft(v, pStringCore)
				          << ", " << mod::rule::internal::getStringRight(v, pStringCore) << ")" << std::endl;
				error();
			}
			if(right > 2) {
				std::cout << "Too many right side edges (" << right << ") for vertex "
				          << get(boost::vertex_index_t(), core, v) << ".\n";
				std::cout << "Label: (" << mod::rule::internal::getStringLeft(v, pStringCore)
				          << ", " << mod::rule::internal::getStringRight(v, pStringCore) << ")" << std::endl;
				error();
			}
		}
	}

	return mod::rule::internal::makeRule(std::move(dpoRule));
}
