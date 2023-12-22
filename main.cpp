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

#include <map>
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
using VDescriptor = boost::adjacency_list<>::vertex_descriptor;

/**
 * This functions takes a morphism and two graphs and does an XOR on all the edges.
 * The resulting graph will be used to detect cycles between the two.
 * @tparam SomeGraph Because of the auto keywords everywhere it is quite bothersome to figure out the correct type...
 * @param vm The morphism mapping the educt to the product
 * @param gEduct The graph representation of the Educt
 * @param gProduct The graph representation of the Product
 * @param lgEduct The labelled graph representation of the Educt
 * @param lgProduct The labelled graph representation of the Product
 * @return An XOR graph from the educt and product graphs where each edge is labelled with either "ADD" or "REMOVE" depending on if it adds an edge or removes an edge from the Educt to the Product.
 */
template< class SomeGraph >
CycleGraph XOR_graphs( VertexMap &vm, const SomeGraph &gEduct, const SomeGraph &gProduct, const LUG &lgEduct, const LUG &lgProduct ) {
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

    std::map<const UVertex, VDescriptor> to_cycle_graph;

    // This lambda function returns a vertex descriptor for the corresponding input vertex. If it has already been inserted into the cycle_graph it will return that, otherwise return a new one and insert it into the map.
    // We are only mapping the left graph to make it simple. In order to call this with a right mapping, simply use the vertex map's right side to map it into the left graph.
    auto getVertex = [&to_cycle_graph, &cycle_graph]( UVertex vertex_to_find ) -> VDescriptor {
        VDescriptor v;
        auto it = to_cycle_graph.find(vertex_to_find);
        if ( it != to_cycle_graph.end() ) {
            v = it->second;
        } else {
            v = boost::add_vertex(cycle_graph);
            to_cycle_graph.insert({vertex_to_find, v});
        }
        return v;
    };

    for ( const auto &vertex_pair: vm ) {
        // Go through all the neighbors of the left hand side, and find the matching right hand side from the vertex_map. Check if this edge also exists on the right side, if so, remove.

        const UVertex left_source = vertex_pair.left;
        const UVertex right_source = vertex_pair.right;

        for ( const UEdge &left_edge : asRange(out_edges(left_source, gEduct)) ) {
            // Find the right side
            const UVertex left_target = target(left_edge, gEduct);
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

                    VDescriptor from = getVertex(left_source);
                    VDescriptor to = getVertex(left_target);
                    if ( boost::edge(from, to, cycle_graph).second == false ) {
                        auto e = boost::add_edge(from, to, cycle_graph).first;
                        // Assuming only '-' and '=' are possible.
                        if (left_bond == "-" && right_bond == "=") {
                            edge_types[e] = "ADD";
                        } else if (left_bond == "=" && right_bond == "-"){
                            edge_types[e] = "REMOVE";
                        }
                    }
                }
            } else {
                // Case [2]
                // There exists an edge on the left hand side but NOT on the right hand side, therefore we want to add this edge.
                VDescriptor from = getVertex(left_source);
                VDescriptor to = getVertex(left_target);
                if ( boost::edge(from, to, cycle_graph).second == false ) {
                    auto e = boost::add_edge(from, to, cycle_graph).first;
                    edge_types[e] = "REMOVE";
                }
            }
        }
        // At this point, do the same but for the right side. However, we only want to add the edges that only appear on the right side.
        for ( const UEdge &right_edge : asRange(out_edges(right_source, gProduct)) ) {
            // Find the left side
            const UVertex right_target = target(right_edge, gProduct);
            const UVertex left_target = vm.right.at(right_target);

            // Do we have a matching left edge?
            const std::pair<UEdge, bool> &left_edge = edge(left_source, left_target, gEduct);

            if ( !left_edge.second ) {
                // Case [3]
                // We want to add this edge.
                VDescriptor from = getVertex(vm.right.at(right_source));
                VDescriptor to = getVertex(left_target);
                if (boost::edge(from, to, cycle_graph).second == false) {
                    auto e = boost::add_edge(from, to, cycle_graph).first;
                    edge_types[e] = "ADD";
                }

            }
        }
    }
    return cycle_graph;
}

std::vector<std::shared_ptr<mod::rule::Rule> > doStuff(
		const std::vector<std::shared_ptr<mod::graph::Graph>> &educts,
		const std::vector<std::shared_ptr<mod::graph::Graph>> &products,
		bool doChemistryCheck) {
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

			constexpr int limit = 20;
			// We will limit to 'limit' permutations, as otherwise the PDF might get too large.
			// Note that we also reject all rules where the atom type would change.
			int permutationCount = 0;
			do {
				VertexMap vm;
				bool valid = true;
				for(int i = 0; i != eductVertices.size(); i++) {
					// check if the molecule data for the atoms matches, otherwise we can't make a rule
					if(mod::graph::internal::getMolecule(eductVertices[i], lgEduct) !=
					   mod::graph::internal::getMolecule(productVertices[i], lgProduct)) {
						valid = false;
						break;
					}
					vm.insert(VertexMap::value_type(eductVertices[i], productVertices[i]));
				}
				if(!valid) continue;

                CycleGraph cycle_graph = XOR_graphs( vm, gEduct, gProduct, lgEduct, lgProduct );
                auto edge_types = boost::get(boost::edge_name, cycle_graph);

//                std::ofstream dotFile("graph.dot");
//                boost::write_graphviz(dotFile, cycle_graph);
//                dotFile.close();

                bool has_cycle = false;
                std::vector<std::vector<Vertex>> cycles;
                cycle_detector vis(has_cycle, cycles);
                boost::depth_first_search(cycle_graph, visitor(vis));
//                std::cout << "The graph has a cycle? " << has_cycle << " Number of cycles: " << cycles.size() << std::endl;

                if ( has_cycle ) {
                    if ( cycles.size() > 1 ) valid = false;
                    // Check if it is a valid cycle.
                    for ( const auto &cycle : cycles ) {
                        if ( cycle.size() % 2 == 1 ) valid = false; // We know that we can't possibly alternate between an add and a remove.

                        std::string prev_operation = "";
                        bool first = true;
                        VDescriptor prev;
                        for ( auto it = cycle.begin(); it != cycle.end(); it++ ) {
                            const VDescriptor v = *it;
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
//                                std::cout << " [ " << this_operation << " ], ";
                            }
//                            std::cout << v;
                            first = false;
                            prev = v;
                        }
                        auto e = boost::edge(prev, *cycle.begin(), cycle_graph).first;
//                        std::cout << " [ " << edge_types[e] << " ], " << *cycle.begin() << "\n";
                    }
//                    std::cout << std::endl;
                }
//                std::cout << (valid ? "Valid cycle detected!\n" : "Invalid cycle detected!\n");

//                std::cout << std::endl;
//                std::cout << std::endl;
                if (!valid) continue;

                // We can also use the XOR graph to detect if more than 2 changes per vertex are happening.
                // If so, discard that solution. TODO

				vertexMaps.push_back(std::move(vm));
				++permutationCount;
				if(permutationCount == limit) break;
			} while(std::next_permutation(eductVertices.begin(), eductVertices.end()));
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
		std::shared_ptr<mod::graph::Graph> g1, g2;
		g1 = mod::graph::Graph::fromSMILES("OCC=O");
		g2 = mod::graph::Graph::fromSMILES("OC=CO");
		educts.push_back(g1);
		products.push_back(g2);
	}
	auto rules = doStuff(educts, products, true);
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
