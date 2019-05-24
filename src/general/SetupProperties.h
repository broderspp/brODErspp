#ifndef SETUPPROPERTIES_H
#define SETUPPROPERTIES_H

#include <string>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

class SetupProperties {
public:
    SetupProperties(std::string& s_file_input);
    ~SetupProperties(){}

    std::string getProblemType(){ return s_problem_type; }
    std::string getMixture(){ return s_mixture; }
    std::string getStateModel(){ return s_state_model; }
    std::string getThermoDB(){ return s_thermo_db; }

    // Function for reading input file
    std::string parseline(std::string line);
    void readInputFile(std::string& s_file_input);
    std::vector< std::string > getInputFileLines();

private:
    std::string s_problem_type;
    std::string s_mixture;
    std::string s_state_model;
    std::string s_thermo_db;

    std::vector< std::string > v_input_file;

    inline void errorInputFileInaccurate();

};

#endif /* SETUPPROPERTIES_H */
