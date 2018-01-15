#ifndef UTILS_HEADER_DEF
#define UTILS_HEADER_DEF

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

template <typename T>
struct input_data {

  std::string output_file;         ///< Output file to which data will be written.
  T ReL;                           ///< Reynolds number based on length.
  T L;                             ///< Reference length (m).
  T Me;                            ///< Edge Mach number.
  T gamma;                         ///< Ratio of specific heats.
  T Te;                            ///< Edge temperature.

};

/**
 * Function used to read data from input file.
 *
 * @param[in] input_file std::string which contains input file name.
 * @return input_data struct which contains input data.
 */
template <typename T>
input_data<T> read_input_file(const std::string& input_file) {

  input_data<T> tmp;
  std::string token,line;
  std::istringstream iss;
  std::ifstream infile(input_file.c_str());
  if (!infile.is_open()) {
    std::cerr << "\nError: Can't open " << input_file << " .\nExiting.\n\n";
    exit(1);
  }
  while (std::getline(infile,line)) {

    // Making sure current line isn't a comment
    if (line.find("#")==std::string::npos) {

      if (line.find("ReL")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.ReL;
        iss.clear();
        std::cout << "Reynolds number: " << tmp.ReL << std::endl;
      }
      if (line.find("lref")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.L;
        iss.clear();
        std::cout << "Reference length (m): " << tmp.L << std::endl;
      }
      if (line.find("Me")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.Me;
        iss.clear();
        std::cout << "Edge Mach number: " << tmp.Me << std::endl;
      }
      if (line.find("gamma")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.gamma;
        iss.clear();
        std::cout << "Ratio of specific heats: " << tmp.gamma << std::endl;
      }
      if (line.find("Te")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.Te;
        iss.clear();
        std::cout << "Edge temperature (K): " << tmp.Te << std::endl;
      }
    }
  }
  infile.close();

  return tmp;

}

#endif
