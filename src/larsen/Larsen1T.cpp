#include "Larsen1T.h"

Larsen1T::Larsen1T(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nen),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_hi(m_ns, 0.),
              v_cpi(m_ns, 0.),
              v_Ji(m_ns, 0.),
              v_divJi(m_ns, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1) {

    // Set the molar weights..
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

// ---------------------------------------------------------------

void Larsen1T::computeHeatFlux()
{
  // Returns the heat flux from the flud particle, according to the flag
  // in the input file.
  // It can be zero (adiabatic particle) or externally given (evaluated
  // from the difference of total enthalpy from one point to the other
  // along the streamline).

  // =====  No heat flux (adiabatic particle)
  if(m_data.getHtotDiffusionFlag().compare("adiabatic") == 0) {
    Q_now = 0.0;
  }
  // =====  Heat flux from enthalpy difference (external)
  else if(m_data.getHtotDiffusionFlag().compare("external") == 0) {
    m_data.getStepEnthalpy(H_ext);
    Q_now = u_now*(H_ext[1] - H_ext[0])/(x_ext[1] - x_ext[0] + 1E-15);
  }
  // =====  Option not recognized
  else {
    std::cerr << "Choice for Heat flux '"<< m_data.getHtotDiffusionFlag()
              << "' not recognized. Aborting." << std::endl;
    exit(1);
  }

}

// ---------------------------------------------------------------

void Larsen1T::computeJdivJ(const double t) {
  // Computes J and divJ according to the requested scheme.
  // Fills the member vectors v_Ji and v_divJi.

  // =====  No diffusion mass fluxes 
  if(m_data.getJdiffFlag().compare("none") == 0) { 
    std::fill(v_Ji.begin(), v_Ji.end(), 0.0); // Set everything to zero.
    std::fill(v_divJi.begin(), v_divJi.end(), 0.0);
  }
  
  // =====  Externally given, compute derivative only streamwise
  else if(m_data.getJdiffFlag().compare("external_streamwise") == 0){

    // For each species, get J and the derivative
    for(size_t id_s = 0; id_s < m_ns; ++id_s) {

      m_data.getStepMassDiffusion(id_s, J_ext);

      v_Ji.at(id_s)    = J_ext[0] 
           + (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);

      v_divJi.at(id_s) = (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15); 
    }
  }

  // =====  Externally given, compute derivatives in stagnation line mode
  else if(m_data.getJdiffFlag().compare("external_stagline_sphere") == 0){

    // Get the radius at current step's integration extremes
    m_data.getStepStaglineR(stagline_r_ext);

    // Compute current radius by linear interpolation among the extremes
    stagline_r_now = stagline_r_ext[0] +
                    (stagline_r_ext[1] - stagline_r_ext[0])/(x_ext[1]
                                    - x_ext[0] + 1E-15)*(t - x_ext[0]);

    // For each species..
    for(size_t id_s = 0; id_s < m_ns; ++id_s) {

      // Get diffusion mass fluxes at current step's integration extremes
      m_data.getStepMassDiffusion(id_s, J_ext);

      // Compute diffusion mass fluxes by (linear) interpolation among extremes 
      v_Ji.at(id_s)    = J_ext[0] 
           + (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);

      v_divJi.at(id_s) = (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15) 
                        - 2/stagline_r_now*v_Ji.at(id_s); 
    }
  }

  // =====  Computed
  else if(m_data.getJdiffFlag().compare("computed_streamwise") == 0){
    std::cerr << "TO BE IMPLEMENTED.." << std::endl;
  } 

  // =====  Request not recognized
  else {
    std::cerr << "DIFFUSION MODEL NOT SUPPORTED YET!" << std::endl;
  }

}

// ---------------------------------------------------------------

vector_type Larsen1T::computedxdt(const vector_type& xState, const double t){

  // This function computes the right hand side of the system (aka the value of the
  // derivatives).
  // 
  // STEP 1) Initialization of variables at the current timestep and extraction
  //         of relevant quantities (rho, h, cp, J..)
  //
  // STEP 2) Computation of the derivatives
  //
  // One note.. Here "x" refers to the state of the system: mass fractions Y and 
  // temperature T. The variable "t" is the coordinate along the streamline (it's
  // not a time), so that "dxdt" is the SPATIAL derivative of the state along the 
  // streamline.
  // 

  // ###########   STEP 1) Preparing the variables   #######################

  // Initializing variables to current state
  m_data.getStepValues(x_ext, u_ext, rho_ext);

  u_now   = u_ext[0] + (u_ext[1] - u_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);
  rho_now = rho_ext[0] + (rho_ext[1] - rho_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);

  // Compute densities..
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_yi[ii] = xState[ii];
    if(xState[ii] < 0.0){   // It might become negative during some kind of integration..
      v_yi[ii] = 0.0;
    }
    v_rhoi[ii] = v_yi[ii] * rho_now;
  }

  T_now = xState[m_ns];

  // Set mixture state
  // NOTE: &rhoi[0] is a tricky way to pass a VECTOR where an ARRAY was asked..
  m_mix.setState(&v_rhoi[0], &T_now, 1);

  // ----  Get properties from the mixture  ----
  m_mix.netProductionRates(&v_omegai[0]);
  m_mix.speciesHOverRT(T_now, &v_hi[0]);
  m_mix.speciesCpOverR(T_now, &v_cpi[0]);

  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_hi[ii]  *= Mutation::RU*T_now/m_mix.speciesMw(ii);
    v_cpi[ii] *= Mutation::RU/m_mix.speciesMw(ii);
  }

  // ----  Get or compute diffusive mass fluxes and heat flux  ----
  computeJdivJ(t); // Just updates the member variables v_Ji, v_divJi
  computeHeatFlux(); // Updates the member variable Q_now


  // ##########   STEP 2) COMPUTING DERIVATIVES   ####################

  // ====> Mass equations
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_dxdt[ii] = (v_omegai[ii] - v_divJi[ii])/rho_now;
  }

  // ====> Translational temperature
  double num_part = 0.0;
  double den      = 0.0;
  for(size_t ii = 0; ii < m_ns; ++ii) {
    num_part += v_hi[ii]*v_dxdt[ii];
    den      += v_yi[ii]*v_cpi[ii];
  }

  // Approximate derivative of u^2
  double Du2Dt = u_now*(pow(u_ext[1],2) 
                        - pow(u_ext[0],2))/(x_ext[1] - x_ext[0] + 1E-15);

  // Temperature equation
  v_dxdt[m_ns] = - (0.5*Du2Dt + num_part - Q_now)/den;

  // Correcting ALL THE derivatives with the velocity 
  // (remember, these are material derivatives!!)
  for(size_t ii = 0; ii < m_ns+ 1; ++ii) {
    v_dxdt[ii] /= u_now;
  }

  return v_dxdt;

};

