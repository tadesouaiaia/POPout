#include "pop_data.h"
#include <iostream>
#include "argparse.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// helper lives here
static std::string pop_name(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string name = (slash == std::string::npos)
        ? path
        : path.substr(slash + 1);

    size_t dot = name.find('.');
    if (dot != std::string::npos)
        name = name.substr(0, dot);
    return name;
}



int main(int argc, char** argv) {
    argparse::ArgumentParser program("POPout");

    program.add_argument("popfiles")
        .help("One or more popfiles")
        .nargs(1, -1);
    
    program.add_argument("--out")
        .help("Output prefix")
        .default_value(std::string("out"));
    
    program.add_argument("--tailSize")
    .help("Tail Bin Size")
    .default_value(float(1.0))
    .scan<'g', float>();

	program.add_argument("--regRange")
    .help("Regression Size")
    .default_value(int(100))
    .scan<'i', int>();



    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n" << program;
        return 1;
    }

    const std::vector<std::string> popfiles = program.get<std::vector<std::string>>("popfiles");
    const std::string out = program.get<std::string>("--out");
    const float tailSize = program.get<float>("--tailSize");
    const int   regRange = program.get<int>("--regRange");
    
    
    
    const std::string outfile = out + "-popout.txt";
    std::ofstream ofs(outfile);
    if (!ofs) {
          std::cerr << "ERROR: cannot open output file " << outfile << "\n";
            return 1;
    }


    PopData::print_tail_header(ofs);
    for (size_t i = 0; i < popfiles.size(); i++) {
        //const std::string& popfile = popfiles[i];
        std::string name = pop_name(popfiles[i]);
        PopData pd(popfiles[i].c_str(), name, tailSize, regRange);   // exact argv[i] semantics (const char*)
        //pd.run_reg(100);
        //#pd.calculate_means();
        pd.calculate_bin_results(80, 80);
		bool qc_ok = (pd.min_pval() * 80.0 >= 0.05); 
        

 
    	pd.calculate_bin_results(regRange, 100);
    	pd.print_tail_popout(ofs, qc_ok);

        
    }

    return 0;
}









