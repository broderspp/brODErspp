#include <iostream>
#include <stdlib.h>

#include "SetupProperties.h"

SetupProperties::SetupProperties(std::string& s_file_input)
        : s_problem_type("empty"),
          s_mixture("empty"),
          s_state_model("empty"),
          s_thermo_db("empty") {

    // Open File
    readInputFile(s_file_input);

    // Check Members not Empty
    errorInputFileInaccurate();

}

// ------------------------------------------------------------------------------------

std::vector< std::string > SetupProperties::getInputFileLines() {

  return v_input_file;

}

// ------------------------------------------------------------------------------------

void SetupProperties::errorInputFileInaccurate(){
   
    if (s_problem_type.compare("empty") == 0) {
        std::cerr << "In the input file a problem type should be provided!" << std::endl;
        exit(1);
    }

    if (s_mixture.compare("empty") == 0) {
        std::cerr << "In the input file a mixture should be provided!" << std::endl;
        exit(1);
    }
    if (s_state_model.compare("empty") == 0) {
        std::cerr << "In the input file a state_model should be provided!" << std::endl;
        exit(1);
    }
    if (s_thermo_db.compare("empty") == 0) {
        std::cerr << "In the input file a thermodynamic database should be provided!" << std::endl;
        exit(1);
    }
    
}

// ------------------------------------------------------------------------------------

std::string SetupProperties::parseline(std::string line)
{
  // Remove comments starting somewhere in the line
  for(size_t ii = 0; ii < line.length(); ++ii) {
    if(line[ii] == '#') { // then the rest is a comment, trim it!
      line = line.substr(0,ii);
    }
  }
  // Removing whitespaces in the beginning of the line
  while( line[0] == ' ' ) {
    line = line.substr(1, line.length());
  }
  // Removing whitespaces in the end of the line
  while( line[line.length()-1] == ' ' ) {
    line = line.substr(0, line.length()-1);
  }
  return line;
};

// ------------------------------------------------------------------------------------

void SetupProperties::readInputFile(std::string& s_file_input)
{

  std::cout << "Loading Input file '" << s_file_input << "'" << std::endl;

  // Open file
  const char * filename = s_file_input.c_str();

  std::ifstream filein;
  filein.open(filename);
  if (!filein) {
      std::cout << "ERROR: couldn't open input file." << std::endl;
      exit(1);
  }

  // Auxiliary variables
  std::string line;

  // read lines
  while(getline(filein, line)){
      line = parseline(line); // Parse the line

      // If the line is not empty, save it in a vector of strings
      if(line.length() > 0) {
         v_input_file.push_back(line);
      }

  }

  // Check for known syntax into the vector of saved lines
  std::string lineNow;
  for(size_t ii = 0; ii < v_input_file.size(); ++ii) {

    lineNow = v_input_file.at(ii);

    if(lineNow.compare("Problem type:")==0) {
        s_problem_type   = v_input_file.at(ii+1);
    }
    if(lineNow.compare("Name of the mixture:")==0) {
        s_mixture        = v_input_file.at(ii+1);
    }
    if(lineNow.compare("State model:")==0) {
        s_state_model    = v_input_file.at(ii+1);
    }
    if(lineNow.compare("Thermodynamic database:")==0) {
        s_thermo_db      = v_input_file.at(ii+1);
    }

  }

  // Variables to save read strings
  std::cout << "\nInput file was read: "                 << std::endl;
  std::cout << "    Problem type: "                      << s_problem_type        << std::endl;
  std::cout << "    Mixture: "                           << s_mixture             << std::endl;
  std::cout << "    State model: "                       << s_state_model         << std::endl;
  std::cout << "    Thermodynamic database: "            << s_thermo_db           << std::endl;
  std::cout << std::endl;

  return;
}

