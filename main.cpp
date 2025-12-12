
// AUTHOR "Xinwei Zhao <zhaoxi@uni-muenster.de>"


#include "src/options.h"
#include "src/run_blat.h"
#include "src/run_minimap2sam.h"
#include "src/run_minimap2paf.h"

#include <iostream>
#include <string>





// Main function for the whole programm, the detailed main function for each progress is divided into own file, This file
// present the rough structure for the programm
int main(int argc, char** argv) {

    // Parse command line options to determine the alignment method to use, which is the basic step for the programm
    options_t options = option_parser(argc, argv);

    // Call the corresponding analysis function according to the user input method, the input type determines the calculation
    // in which form, we support PSL, PAF and SAM input form
    std::string alignment_method = options.input_type;
    if (alignment_method == "blat") {
        run_blat(options);
    } else if (alignment_method == "minimap2sam") {
        run_minimap2sam(options);
    } else if (alignment_method == "minimap2paf") {
        run_minimap2paf(options);
    } else {
        std::cerr << "Invalid alignment method: " << alignment_method << "please check you input with -m, we support blat, minimap2sam and minimap2paf as input" << std::endl;
        return 1;
    }
    return 0;
}


