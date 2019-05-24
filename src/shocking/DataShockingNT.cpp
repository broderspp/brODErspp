#include <sstream>
#include "DataShockingNT.h"

DataShockingNT::DataShockingNT(Mutation::Mixture& l_mix)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(1), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_meq+n_eneq),
                          pos_V(n_sp),
                          pos_T(n_sp+n_meq),
                          m_P(0.0),
                          m_V(0.0),
                          m_rho(0.0),
                          m_mdot(0.0),
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0)
{}

// -----------------------------------------------------------------

DataShockingNT::DataShockingNT(Mutation::Mixture& l_mix, 
                               const std::vector< std::string > l_input_file)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(1), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_meq+n_eneq),
                          pos_V(n_sp),
                          pos_T(n_sp+n_meq),
                          m_P(0.0),
                          m_V(0.0),
                          m_rho(0.0),
                          m_mdot(0.0),
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0)
{
    // Parse inut file
    inputFileParse(l_input_file);
    buildState();

    // Check the input file and print settings
    checkAndPrintFreeStream(); // Do this AFTER buildState()
}

// -----------------------------------------------------------------

DataShockingNT::~DataShockingNT() {}

// -----------------------------------------------------------------

void DataShockingNT::inputFileParse(const std::vector< std::string > l_input_file){
  // Parses the input file lines

  // working variables
  double T_now_str; 
  std::vector<double> v_T_now;

  // Get Pressure, Temperatures and Velocity
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {

    // PRESSURE
    if( l_input_file.at(id_l).compare("FS Press:") == 0) {
      m_P = atof(l_input_file.at(id_l+1).c_str());
    }

    // TEMPERATURES
    if( l_input_file.at(id_l).compare("FS Temp:") == 0) {

      // Unpack the line in the temperatures
      std::istringstream iss(l_input_file.at(id_l+1));
      while(iss >> T_now_str) {  // Read the values up to N
        v_T_now.push_back(T_now_str);
      }

      // Check if the number of temperatures given is fine
      if(v_T_now.size() != n_eneq) {
        std::cerr << " ATTENTION: " << v_T_now.size() << " temperatures have been"
                  << " specified, while the state model supports " << n_eneq 
                  << ". Check the input file." << std::endl;
        std::cerr << " Aborting." << std::endl;
        exit(1);
      } else { // Assign temperatures
        for(size_t en_id = 0; en_id < n_eneq; ++en_id) {
          v_T[en_id] = v_T_now.at(en_id);
        }
      }
    }

    // VELOCITY
    if( l_input_file.at(id_l).compare("FS Vel:") == 0) {
      m_V = atof(l_input_file.at(id_l+1).c_str());
    }
  }
}

// -----------------------------------------------------------------

void DataShockingNT::buildState(){

    m_mix.equilibriumComposition(v_T[0], m_P, &v_xi[0]);

    m_mix.convert<Mutation::Thermodynamics::X_TO_Y>(&v_xi[0], &v_yi[0]);
    m_rho = m_mix.density(v_T[0], m_P, &v_xi[0]);

    // Convert yi to rhoi
    for (int i_sp = 0; i_sp < n_sp; ++i_sp){
        v_rhoi[i_sp] = m_rho *v_yi[i_sp];
    }

    // compute rho*u
    m_mdot = m_rho * m_V;

    fillStateVectorX();
}

// -----------------------------------------------------------------

void DataShockingNT::buildStatePostShock(){

    m_mix.convert<Mutation::Thermodynamics::Y_TO_X>(&v_yi[0], &v_xi[0]);
    m_rho = m_mix.density(v_T[0], m_P, &v_xi[0]);

    // compute rho*u
    // Important Note:
    // m_mdot is slightly changed after the shock...
    // Taking it from the pre shock and computing rho that way might be more consistent
    m_mdot = m_rho * m_V;

    fillStateVectorX();
}

// -----------------------------------------------------------------

void DataShockingNT::fillStateVectorX(){

    for (int i_sp = 0; i_sp < n_sp; ++i_sp){
        v_X[i_sp] = v_yi[i_sp];
    }
    v_X[pos_V] = m_V;
    for (int i_en = 0; i_en < n_eneq; ++i_en){
        v_X[pos_T+i_en] = v_T[i_en];
    }
}

// -----------------------------------------------------------------

void DataShockingNT::checkAndPrintFreeStream(){

  // =====  Check settings (from input file)  ======

  // Check freestream Pressure
  if (m_P == 0.0){
      std::cerr << " ATTENTION: pressure not correctly set in the input file!"
                << std::endl;
      std::cerr << std::endl << "Aborting." << std::endl;
      exit(1);
  }

  // Check freestream Temperature
  if (v_T[0] == 0.0){
      std::cerr << " ATTENTION: temperature not correctly set in the input file!"
                << std::endl;
      std::cerr << std::endl << "Aborting." << std::endl;
      exit(1);
  }

  // Check freestream Velocity
  if (m_V == 0.0){
      std::cerr << " ATTENTION: velocity not correctly set in the input file!"
                << std::endl;
      std::cerr << std::endl << "Aborting." << std::endl;
      exit(1);
  }

  // ======  Print Settings  ======
  // Prints the freestream conditions, read from file
  std::cout << "Free-stream conditions:"    << std::endl;
  std::cout << "    Press [Pa]: " << m_P    << std::endl;
  std::cout << "    Temp [K]:   " << v_T[0] << std::endl;
  std::cout << "    Vel [m/s]:  " << m_V    << std::endl;
  std::cout << "    Computed mass flux [kg/m2 s]: " << m_mdot << std::endl;
  std::cout << "    Computed density [kg/m3]:     " << m_rho  << std::endl;
  std::cout << std::endl;

}

