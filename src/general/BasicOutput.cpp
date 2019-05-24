#include "BasicOutput.h"

BasicOutput::BasicOutput(const Mutation::Mixture& l_mix, Data& l_data, 
                         const std::vector< std::string > l_input_file)
           : p_mix(&l_mix),
             p_data(&l_data),
             n_eq(l_data.nEquations()),
             n_sp(l_mix.nSpecies()),
             n_steps(1),
             s_species_output_type("empty"),
             v_species(l_mix.nSpecies(), 0.0),
             v_sol(l_data.nEquations(), 0.0)
{   
    steps_counter = 0;     

    // Parse input file and print options
    inputFileParse(l_input_file);
}

// -----------------------------------------------------------------

BasicOutput::~BasicOutput(){}

// -----------------------------------------------------------------

void BasicOutput::operator()(const vector_type& l_sol, const double x)
{
  // Increment the internal counter
  steps_counter++;

  if(steps_counter >= n_steps) { // If it's time to print..

    steps_counter = 0;   // restart the counter

    // Initialize solution vector v_sol (it will be eventually changed)
    for(size_t i_eq = 0; i_eq < n_eq; ++i_eq){
      v_sol.at(i_eq) = l_sol[i_eq];
    }

    // Convert in mole fraction if needed
    if(s_species_output_type.compare("mole_fraction") == 0)  // User wants mole fraction as output
    {
      // Copy species into a vector
      for(size_t i_sp = 0; i_sp < n_sp; ++i_sp) {
        v_species.at(i_sp) = l_sol[i_sp];
      }

      // Convert the vector
      p_mix->convert<Mutation::Thermodynamics::Y_TO_X>(&v_species[0], &v_species[0]);

      // Rewrite solution vector
      for(size_t i_sp = 0; i_sp < n_sp; ++i_sp) {
        v_sol.at(i_sp) = v_species.at(i_sp);
      }
    }

    // Print solution vector

    std::cout.precision(15);
    std::cout << "Sol: " << x << "  ";
    for (int i_eq = 0; i_eq < n_eq; i_eq++){
        std::cout.precision(15);
        std::cout << v_sol.at(i_eq) << "  ";
    }
    std::cout << std::endl;
  }
}

// -----------------------------------------------------------------

void BasicOutput::inputFileParse(const std::vector< std::string > l_input_file)
{
  // ====== Parses the input file for getting infos on output options

  // Get the output type for species (mass or mole fraction)
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("Species output type:") == 0) {
      s_species_output_type = l_input_file.at(id_l+1);
    }

    if( l_input_file.at(id_l).compare("Print sol each N steps:") == 0) {
      n_steps = atoi(l_input_file.at(id_l+1).c_str());
    }
  }

  // ====== Check the options

  // species output type
  if (s_species_output_type.compare("empty") == 0) { // How to output species mass/mole
    std::cerr << " ATTENTION: species output type not specified in the input "
              << "file.\n  Setting to default value: 'mass_fraction'."<< std::endl;
    s_species_output_type = "mass_fraction";
  }

  if (n_steps < 1) {
    std::cerr << " ATTENTION: output-each-N-timesteps not smaller than 1."
              << " Check the input file!" << std::endl;
    std::cerr << " Aborting." << std::endl;
    exit(1);
  }

  // ======= Print the options

  std::cout << "Output settings:" << std::endl;

  std::cout << "    Species output type: " << s_species_output_type << std::endl;

  std::cout << "    Solution is printed each " << n_steps << " steps."
            << std::endl << std::endl;
}

// -----------------------------------------------------------------

void BasicOutput::printSolutionHeader(std::string s_style)
{
  std::string s_spaces = "       "; // spaces..

  // === Position
  std::string s_position = "Position x [m]" + s_spaces;
  std::string s_curv_abs = "Curvilinear Abscissa s [m]" + s_spaces;
  
  // === Chemical species

  // Solution output type
  std::string frac_prefix;
  if(s_species_output_type.compare("mass_fraction") == 0) {
    frac_prefix = "Y_";
  } else if(s_species_output_type.compare("mole_fraction") == 0) {
    frac_prefix = "X_";
  } 

  std::string s_species;
  for(size_t id_s = 0; id_s < n_sp; ++id_s) {
    s_species += frac_prefix + p_mix->speciesName(id_s) + s_spaces;
  }

  // === Velocity
  std::string s_vel = "U [m/s]" + s_spaces;

  // === Temperature(s)
  std::string s_temp;

  // multi-Temperature model
  if (p_mix->nEnergyEqns() > 1) {
    std::ostringstream ss; // temporary

    for(size_t id_T = 0; id_T < p_mix->nEnergyEqns(); ++id_T) {
      ss << "T_" << id_T << " [K]" << s_spaces;
    }
    s_temp = ss.str();
  }
  // 1 Temperature model
  else {  s_temp = "T [K]" + s_spaces;  }

  // === Output header
  std::string s_header;
  // Shocking
  if(s_style.compare("shocking") == 0){
    s_header = s_position + s_species + s_vel + s_temp;
  // LARSEN
  } else if(s_style.compare("larsen") == 0) {
    s_header = s_curv_abs + s_species + s_temp;
  }
  
  std::cout << s_header << std::endl;
}
