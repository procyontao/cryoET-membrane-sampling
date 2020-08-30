#pragma once
#include <cstddef>
#include <string>

void imod_process(
		const std::string &input,
		const std::string &reference,
		const std::string &output,
		const std::size_t halfboxsize,
		const std::size_t sampling,
		const float distance_factor,
		const std::ptrdiff_t zbridge,
		const std::ptrdiff_t intmode,
		const bool debug = false);
