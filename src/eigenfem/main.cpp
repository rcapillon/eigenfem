//
// main.cpp
//

#include <iostream>
#include <ctime>

#include "solvers.h"
#include "io.h"


int main(int argc, char *argv[])
{
    time_t global_timer_start = time(nullptr); // Starts timer for whole code execution

    // Value used for pi
    const float PI = 3.141592653589793;

    // count and print of arguments fed to main
    // std::cout << "Number of arguments: " << argc << std::endl;
    // for (size_t i = 0; i < argc; i++)
    // {
    //     std::cout << "Argument " << i << ":" << std::endl;
    //     std::cout << argv[i] << std::endl << std::endl;
    // }

    std::string path_to_input_file = argv[1];

    // Parse input file
    InputParser ip(path_to_input_file);
    ip.parse_input_file();
    // for (size_t i = 0; i < ip.lines.size(); i++)
    // {
    //     std::cout << ip.lines[i] << std::endl;        
    // }
    
    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}