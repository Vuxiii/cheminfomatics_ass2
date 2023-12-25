#include <iostream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <optional>

int main() {
    int n, r;
    std::cin >> n;
    std::cin >> r;

    std::vector<bool> v(n);
    std::fill(v.end() - r, v.end(), true);

    do {
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                std::cout << (i + 1) << " ";
            }
        }
        std::cout << "\n";
    } while (std::next_permutation(v.begin(), v.end()));
    return 0;
}

template< typename T >
struct Combinatorical_Generator {
    Combinatorical_Generator( std::vector<T> data ) 
    : data(data) 
    {
        for ( int i = 0; i < data.size(); ++i ) {
            chosen.push_back(false);
        }
    }

    std::vector<T > current() const {
        std::vector<T > out;
        for ( T v : data ) {
            out.push_back(v);
        }
        return out;
    }

    std::optional<std::vector<T>> next() {
        if ( false == std::next_permutation(chosen.begin(), chosen.end()) ) {
            // No more permutations
            return std::nullopt;
        }
        return current();
    }

    void reset_chosen( int n ) {
        std::fill(chosen.begin(), chosen.end()-n, false);
        std::fill(chosen.end()-n, chosen.end(), true);
    }

private:
    std::vector<T> data;
    std::vector<bool> chosen;
};