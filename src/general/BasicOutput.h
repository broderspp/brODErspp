#ifndef BASICOUTPUT_H
#define BASICOUTPUT_H

#include <iostream>
#include <boost/numeric/odeint.hpp>

#include "mutation++.h"
#include "Data.h"

typedef boost::numeric::ublas::vector<double> vector_type;

class BasicOutput {
public:
    BasicOutput(const Mutation::Mixture& l_mix, Data& l_data,
                const std::vector< std::string > l_input_file);
    ~BasicOutput();
    void operator()(const vector_type& l_sol, const double x);

    void printSolutionHeader(std::string s_style);

private:
    // Variables
    const size_t n_eq;
    const size_t n_sp;
    size_t n_steps;
    size_t steps_counter;

    const Mutation::Mixture* p_mix;
    Data* p_data;

    std::string s_species_output_type;
 
    std::vector<double> v_species;
    std::vector<double> v_sol;
  
    // Functions
    void inputFileParse(const std::vector< std::string > l_input_file);
};

#endif /* BASICOUTPUT_H */
