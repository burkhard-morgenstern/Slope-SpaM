#include <iosfwd>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <threadpool/ThreadPool.h>

namespace spam {
	struct pattern {
		std::string bits;
		std::vector<size_t> indices;
		
		pattern(std::string _bits)
			: bits(std::move(_bits))
		{
			for (size_t i = 0; i < bits.size(); ++i) {
				if (bits[i] == '1') {
					indices.push_back(i);
				}
			}
		}

		auto begin() const {
			return indices.begin();
		}

		auto end() const {
			return indices.end();
		}

		size_t size() const {
			return indices.size();
		}
	};
struct sequence {
	std::string name;
	std::string bases;

	static auto from_multi_fasta_file(std::string const& filename)
		-> std::optional<std::vector<sequence>>;
};

auto operator>>(std::istream& is, sequence& seq)
	-> std::istream&;


class distance_matrix {
	std::vector<spam::sequence> sequences;
	spam::pattern pattern;
	size_t kmax;

	std::vector<std::vector<std::string>> wordlists;
	std::vector<std::vector<std::pair<size_t, std::string_view>>> viewlists;
	std::vector<std::vector<double>> matrix;

	ThreadPool threadpool;

public:
	distance_matrix(
		std::vector<spam::sequence> const& sequences,
		spam::pattern const& pattern,
		size_t kmax);

	auto size() const
		-> size_t;

	auto column(size_t i) const
		-> std::pair<spam::sequence const&, std::vector<double> const&>;

private:
	void calculate();

	void initialize_matrix();

	void create_wordlists();
	void create_wordlists_par();

	void create_viewlists();

	void calculate_matrix();
	void calculate_matrix_par();

	auto calculate_element(size_t i, size_t j) const
		-> std::pair<double, double>;

	static auto merge_viewlists(
		std::vector<std::pair<size_t, std::string_view>> const& lhs,
		std::vector<std::pair<size_t, std::string_view>> const& rhs)
		-> std::vector<std::pair<size_t, std::string_view>>;

	static auto calculate_matches(
		std::vector<std::pair<size_t, std::string_view>> wordlist,
		size_t kmin,
		size_t kmax)
		-> std::vector<size_t>;

	static auto calculate_distance(std::vector<size_t> const& matches,
		size_t length1,
		size_t length2,
		size_t kmin)
		-> std::pair<double, double>;
};

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix);

} // namespace spam
