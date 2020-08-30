#include <iostream>
#include <boost/program_options.hpp>
#include <common.h>
#include "imodinterp.h"

int main(int argc, char* argv[]) {
    
    namespace po = boost::program_options;
		po::options_description desc(
					"TCL-SURF (" __DATE__ " " __TIME__ ").\nAllowed options",
					128);

    desc.add_options()
        ("help,h", "produce help message")
				("input,i", po::value<std::string>(), "[p] path to the input model.")
				("reference,r", po::value<std::string>(), "[p] path to the reference model.")
				("output,o", po::value<std::string>(), "[p] path to the output.")
				("halfboxsize,l", po::value<std::size_t>(), "[ull] the half size of the box side.")
				("sampling,s", po::value<std::size_t>()->default_value(2), "[ull] the oversampling rate.")
				("distancefactor,f", po::value<float>()->default_value(1.1f, "1.1"), "[f] the distance factor.")
				("zbridge,z", po::value<std::ptrdiff_t>()->default_value(32), "[ll] number of slices over which to connect contours for interpolation.")
				("interpmethod,t", po::value<std::ptrdiff_t>()->default_value(std::ptrdiff_t(intmodes::INT_LINEAR)), "[ll] interpolation method.")
				("debug,d", po::bool_switch(), "[b] whether to output debug information.")
				;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    
		imod_process(
				vm["input"].as<std::string>(),
				vm["reference"].as<std::string>(),
				vm["output"].as<std::string>(),
				vm["halfboxsize"].as<std::size_t>(),
				vm["sampling"].as<std::size_t>(),
				vm["distancefactor"].as<float>(),
				vm["zbridge"].as<std::ptrdiff_t>(),
				vm["interpmethod"].as<std::ptrdiff_t>(),
				vm["debug"].as<bool>()
				);

		return 0;
}
