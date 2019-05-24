#include "DataLarsenTTv.h"

DataLarsenTTv::DataLarsenTTv(Mutation::Mixture& l_mix, 
                             std::vector< std::string > l_input_file)
                        : m_mix(l_mix),
                          n_sp(l_mix.nSpecies()),
                          n_meq(0), // An 1D solver can have 1 component of velocity
                          n_eneq(l_mix.nEnergyEqns()),
                          n_eq(n_sp+n_eneq),
                          v_rhoi(n_sp, 0.0),
                          v_xi(n_sp, 0.0),
                          v_yi(n_sp, 0.0),
                          v_T(n_eneq, 0.0),
                          v_X(n_eq, 0.0),
                          s_ext_flow_filename("empty"),
                          s_species_input_type("empty"),
                          s_diffusion_Htot("empty"),
                          s_diffusion_J("empty"),
                          s_diffusion_J_scheme("empty")
{
  // ----  Check flags (from input file)  ----
  setFlags(l_input_file);

  // ----  Resize variables  ----
  Y0_ext.resize(m_mix.nSpecies()); // Right dimension to Y0
  Y_ext.resize(n_sp);              // Also to Y_ext
  ud_ext.resize(n_sp);             // Also to ud_ext, diffusion velocities
  J_ext.resize(n_sp);              // Also to J_ext, diffusion mass fluxes
  hin_ext.resize(n_sp);            // Also to hin_ext, internal enthalpies
}

// -----------------------------------------------------------------------------

DataLarsenTTv::~DataLarsenTTv() {}

// -----------------------------------------------------------------------------

vector_type DataLarsenTTv::getInitialState()
{
  vector_type x(n_eq); // Species + 1 temperature
  
  // chemical species
  for(size_t ii = 0; ii < m_mix.nSpecies(); ++ii) {
    x[ii] = Y0_ext[ii];
  }
  // translational temperature
  x[m_mix.nSpecies()] = T_ext.at(0);
  // vibrational temperature
  x[m_mix.nSpecies()+1] = Tv_ext.at(0);

  return x;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::setCurrentStep(int pos)
{
  currentStep = pos;

  return; 
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepValues(double x_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function
  x_now[0]   = x_ext.at(currentStep);
  x_now[1]   = x_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepValues(double x_now[], double u_now[], double rho_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function
  x_now[0]   = x_ext.at(currentStep);
  x_now[1]   = x_ext.at(currentStep+1);
  u_now[0]   = u_ext.at(currentStep);
  u_now[1]   = u_ext.at(currentStep+1);
  rho_now[0] = rho_ext.at(currentStep);
  rho_now[1] = rho_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepEnthalpy(double H_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function
  H_now[0]   = H_ext.at(currentStep);
  H_now[1]   = H_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepInternalEnthalpy(size_t s_id, double hin_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function.
  // s_id is the species id of whom the mass diffusion J is seeked.

  hin_now[0]   = hin_ext[s_id].at(currentStep);
  hin_now[1]   = hin_ext[s_id].at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepMassDiffusion(size_t s_id, double J_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function.
  // s_id is the species id of whom the mass diffusion J is seeked.

  J_now[0]   = J_ext[s_id].at(currentStep);
  J_now[1]   = J_ext[s_id].at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepDiffusionVel(size_t s_id, double ud_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function.
  // s_id is the species id of whom the diffusion velocity ud is seeked.

  ud_now[0]   = ud_ext[s_id].at(currentStep);
  ud_now[1]   = ud_ext[s_id].at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepStaglineR(double stagline_r_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function.

  stagline_r_now[0]   = stagline_r_ext.at(currentStep);
  stagline_r_now[1]   = stagline_r_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::getStepStaglineV(double stagline_v_now[])
{
  // Gets the values at position currentStep (private variable to be set with
  // the setCurrentStep() function.

  stagline_v_now[0]   = stagline_v_ext.at(currentStep);
  stagline_v_now[1]   = stagline_v_ext.at(currentStep+1);

  return;
}

// -----------------------------------------------------------------------------

int DataLarsenTTv::getStrNele()
{
  // Returns the number of elements of the imported streamline
  return x_ext.size();
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::readDataFile()
{
  std::string lineRead;
  std::ifstream filein(s_ext_flow_filename.c_str());
  if(filein.is_open() == 0)
  {
    std::cerr << "I could not open the file '" << s_ext_flow_filename
              << "' Does it exist?\nAborting!\n";
    exit(EXIT_FAILURE);
  }

  getline(filein,lineRead); // read first line
  std::stringstream nowstream(lineRead); // and unpack it in the nowstream variable

  // working variables
  double xNow, TNow, TvNow, rhoNow, uNow, FracNow, udiNow;
  double stagline_rNow, stagline_vNow;
  std::vector<double> FracNowVect(n_sp, 0.0);  // mass or mole fractions

  // -----   Reading first line and saving initial conditions
  nowstream >> xNow >> TNow >> TvNow >> rhoNow >> uNow;

  x_ext.push_back(xNow);
  T_ext.push_back(TNow);
  Tv_ext.push_back(TvNow);
  u_ext.push_back(uNow);
  rho_ext.push_back(rhoNow);

  // -----   then read and save value for species mass fractions
  for(size_t ii = 0; ii < n_sp; ++ii) // first, read all the species and fill vector
  {
    nowstream >> FracNow;
    FracNowVect.at(ii) = FracNow;
  }

  if(s_species_input_type.compare("mole_fraction") == 0) // convert if needed to mass fract
  {
    m_mix.convert<Mutation::Thermodynamics::X_TO_Y>(&FracNowVect[0], &FracNowVect[0]);
  }
  else if(s_species_input_type.compare("mass_fraction") == 0) // nothing to be done
  { /* nothing to be done in this case */ }
  else  // keyword not recognized
  { 
    std::cerr << "ATTENTION: value not recognized for ''Species input type:'' flag."
              << " Only 'mass_fraction' or 'mole_fraction' supported."
              << " Check the input file!" << std::endl;
    exit(1);
  }
  
  for(size_t ii = 0; ii < n_sp; ++ii) // then, save it to data structure
  {
    Y0_ext.at(ii) = FracNowVect.at(ii);      // Also, save the INITIAL value!!
    Y_ext[ii].push_back(FracNowVect.at(ii));
  }

  // -----   Diffusion velocities
  // <<<<<<-----------------   this will be modified.. for now let's keep it like this.
  // <<<<<<-----------------   In the future I may want to have diffusion velocity for
  // <<<<<<-----------------   for just some species, such as the ablated species (this
  // <<<<<<-----------------   is the only way I have to take them into account)
  for(size_t ii = 0; ii < n_sp; ++ii)
  {
    if( (getJdiffFlag().compare("external_streamwise") == 0) ||
        (getJdiffFlag().compare("external_stagline_sphere") == 0) )
    { // If they are externally-given:
      nowstream >> udiNow;
    } else {                                  // otherwise set them to zero
      udiNow = 0.0;
    }
    
    ud_ext[ii].push_back(udiNow); // Save the value into diffusion velocities
  }

  // If I'm working on stagline, I need to get two more parameters from the input file
  if (getJdiffFlag().compare("external_stagline_sphere") == 0)
  {
    nowstream >> stagline_rNow >> stagline_vNow;
    stagline_r_ext.push_back(stagline_rNow);
    stagline_v_ext.push_back(stagline_vNow);
  }

  // Keep reading now until the EOF
  while(getline(filein,lineRead))
  {
    nowstream.clear();          // clear stream variable nowstream
    nowstream.str(lineRead);    // put newly read line into nowstream

    // ...and unpack it into my variables
    nowstream >> xNow >> TNow >> TvNow >> rhoNow >> uNow; 

    x_ext.push_back(xNow);
    T_ext.push_back(TNow);
    Tv_ext.push_back(TvNow);
    u_ext.push_back(uNow);
    rho_ext.push_back(rhoNow);

    // -----   then read and save value for species mass fractions
    for(size_t ii = 0; ii < n_sp; ++ii) // first, read all the species and fill vector
    {
      nowstream >> FracNow;
      FracNowVect.at(ii) = FracNow;
    }
  
    if(s_species_input_type.compare("mole_fraction") == 0) // convert if needed
    {                                                      // to mole fractions
      m_mix.convert<Mutation::Thermodynamics::X_TO_Y>(&FracNowVect[0], &FracNowVect[0]);
    }
    
    for(size_t ii = 0; ii < n_sp; ++ii) // then, save it to data structure
    {
      Y_ext[ii].push_back(FracNowVect.at(ii));
    }

    // Diffusion velocities
    // <<<<<<-----------------   this will be modified.. for now let's keep it like this.
    // <<<<<<-----------------   In the future I may want to have diffusion velocity for
    // <<<<<<-----------------   for just some species, such as the ablated species (this
    // <<<<<<-----------------   is the only way I have to take them into account)
    for(size_t ii = 0; ii < n_sp; ++ii)
    {
      if( (getJdiffFlag().compare("external_streamwise") == 0) ||
          (getJdiffFlag().compare("external_stagline_sphere") == 0) )
      { // If they are externally-given:
        nowstream >> udiNow;
      } else {                                  // otherwise set them to zero
        udiNow = 0.0;
      }

      ud_ext[ii].push_back(udiNow); // Save the value into diffusion velocities
    }

    // If I'm working on stagline, I need to get two more parameters from the input file
    if (getJdiffFlag().compare("external_stagline_sphere") == 0)
    {
      nowstream >> stagline_rNow >> stagline_vNow;
      stagline_r_ext.push_back(stagline_rNow);
      stagline_v_ext.push_back(stagline_vNow);
    }

  }

  filein.close();  // I'm done with this data file

  return;
}

// -----------------------------------------------------------------------------

std::string DataLarsenTTv::getHtotDiffusionFlag()
{
  // Returns the flag telling wether ehthalpy is to be conserved (adiabatic 
  // fluid particle) or diffuses
  return s_diffusion_Htot;
}

// -----------------------------------------------------------------------------

std::string DataLarsenTTv::getJdiffFlag()
{
  // Returns the flag telling wether diffusive mass fluxes are to be imported
  // from the external computation or not.
  return s_diffusion_J;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::computeEnthalpy()
{
  // Working variables
  double T_now[2];
  std::vector<double> v_hi_now;  v_hi_now.resize(n_eneq*n_sp);

  // Assign size to total enthalpy vector and internal enthalpy matrix
  H_ext.resize(x_ext.size());

  for(size_t s_id = 0; s_id < n_sp; ++s_id) {
    hin_ext.at(s_id).resize(x_ext.size());
  }

  // Marching on the streamline
  for(size_t kk = 0; kk < x_ext.size(); ++kk)
  {
    T_now[0] = T_ext.at(kk);
    T_now[1] = Tv_ext.at(kk);

    // Compute densities
    for(size_t ii = 0; ii < n_sp; ++ii)
    {
      v_rhoi[ii] = Y_ext[ii][kk]*rho_ext[kk];
    }

    // Set state
    m_mix.setState(&v_rhoi[0], T_now, 1);

    // Compute total enthalpy
    H_ext.at(kk) = 0.5*pow(u_ext.at(kk), 2) + m_mix.mixtureHMass();

    // Compute internal enthalpies
    m_mix.getEnthalpiesMass(&v_hi_now[0]);

    for(size_t s_id = 0; s_id < n_sp; ++s_id) {
      hin_ext.at(s_id).at(kk) = v_hi_now.at(s_id + n_sp);
    }
  }
 
  return;
}
        
// -----------------------------------------------------------------------------

void DataLarsenTTv::computeMassDiffusion()
{
  // Resize the diffusive mass fluxes vector
  for(size_t i_sp = 0; i_sp < n_sp; ++i_sp)
  {
    J_ext[i_sp].resize(x_ext.size());
  }

  // Marching on the streamline
  for(size_t kk = 0; kk < x_ext.size(); ++kk)
  {
    // Compute densities and species diffusive mass flux
    for(size_t i_sp = 0; i_sp < n_sp; ++i_sp)
    {
      v_rhoi[i_sp]    = Y_ext[i_sp][kk]*rho_ext[kk];
      J_ext[i_sp][kk] = v_rhoi[i_sp]*ud_ext[i_sp][kk];
    }

  }

  return;
}

// -----------------------------------------------------------------------------

void DataLarsenTTv::setFlags(std::vector< std::string > l_input_file)
{
  // This function sets some flags from the input file.
  // Afterwards, checks that important flags were set and alerts the user if it's
  // not the case, setting them to default values.
  // 
  // The input file is given as a vector of strings, prepared by the m_setup_props
  // object, available as variable "l_input_file".

  // ----  Set flags and stuff

  // Get the name for the external flowfield file
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("External flowfield filename:") == 0) {
      s_ext_flow_filename = l_input_file.at(id_l+1);
    }
  }

  // Get the input type for species (mass or mole fraction)
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("Species input type:") == 0) {
      s_species_input_type = l_input_file.at(id_l+1);
    }
  }

  // Get flag for total enthalpy diffusion
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("Total enthalpy diffusion:") == 0 ) {
      s_diffusion_Htot = l_input_file.at(id_l+1);
    }
  }

  // Get flag for mass fluxes
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("Diffusive mass fluxes:") == 0 ) {
      s_diffusion_J = l_input_file.at(id_l+1);
    }
  }

  // Get flag for mass fluxes scheme
  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {
    if( l_input_file.at(id_l).compare("Diffusive mass fluxes scheme:") == 0 ) {
      s_diffusion_J_scheme = l_input_file.at(id_l+1);
    }
  }

  // ----  Perform some checks and eventually alert the user

  if (s_ext_flow_filename.compare("empty") == 0) { // Name of external simulation file
    std::cerr << "ATTENTION: name for the simulation file not specified! ";
    std::cerr << "It may be specified using the flag 'External simulation file:'";
    std::cerr << std::endl << "Aborting." << std::endl;
    exit(1);
  }

  if (s_species_input_type.compare("empty") == 0) { // Species definition in input file
    std::cerr << " ATTENTION: species input type not specified in the input "
              << "file.\n  Setting to default value: 'mass_fraction'."<< std::endl;
    s_species_input_type = "mass_fraction";
  }

  if (s_diffusion_Htot.compare("empty") == 0) { // Total enthalpy diffusion
    std::cerr << " ATTENTION: model for total enthalpy diffusion not specified in the "
              << "input file.\n  Setting to default value: 'adiabatic'." << std::endl;
    s_diffusion_Htot = "adiabatic";
  }

  if (s_diffusion_J.compare("empty") == 0) { // Diffusive mass fluxes
    std::cerr << " ATTENTION: model for diffusive mass fluxes not specified in the "
              << "input file.\n  Setting to default value: 'none'." << std::endl;
    s_diffusion_J = "none";
  }

  // ----  Print what was found
  std::cout << "\nLARSEN settings from input file:" << std::endl;

  std::cout << "    External flowfield filename: "  << s_ext_flow_filename  << std::endl;
  std::cout << "    Species input type: "           << s_species_input_type << std::endl;
  std::cout << "    Total enthalpy diffusion: "     << s_diffusion_Htot     << std::endl;
  std::cout << "    Mass fluxes: "                  << s_diffusion_J        << std::endl;

  std::cout << std::endl;

}

