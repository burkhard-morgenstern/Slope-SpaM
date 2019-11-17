#pragma once

#include <threadpool/ThreadPool.h>

#include "pattern.hpp"
#include "sequence.hpp"
#include "wordlist.hpp"

namespace spam {

class insufficient_pattern_exception
	: std::exception
{
	std::string message;

public:
	explicit insufficient_pattern_exception(std::string message)
		: message(std::move(message))
	{}

	auto what() const noexcept
		-> const char *
		override
	{
		return message.c_str();
	}
};

class distance_matrix {
	std::vector<spam::sequence> sequences;
	spam::pattern pattern;
	std::vector<size_t> wordlengths;

	std::vector<wordlist> wordlists;
	std::vector<std::vector<double>> matrix;

	std::optional<std::shared_ptr<ThreadPool>> threadpool;

public:
	distance_matrix(
		std::vector<spam::sequence> sequences,
		spam::pattern pattern,
		std::vector<size_t> wordlengths,
		std::optional<std::shared_ptr<ThreadPool>> threadpool = {});

	auto size() const
		-> size_t;

	auto column(size_t i) const
		-> std::pair<spam::sequence const&, std::vector<double> const&>;

private:
	void calculate();

	void initialize_matrix();
	void create_wordlists();
	void create_wordlists_seq();
	void create_wordlists_par();
	void calculate_matrix();
	void calculate_matrix_seq();
	void calculate_matrix_par();

	auto calculate_element(size_t i, size_t j) const
		-> std::pair<double, double>;
};

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix);

} // namespace spam
