#include "LarsenTTv.h"

LarsenTTv::LarsenTTv(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nen),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_T(m_nen, 0.),
              v_hi(m_ns*m_nen, 0.),
//              v_ei(m_ns*m_nen, 0.), // REMOVE ME
              v_cpi(m_ns*m_nen, 0.),
//              v_cvi(m_ns*m_nen, 0.), // REMOVE ME
              v_Ji(m_ns, 0.),
              v_divJi(m_ns, 0.),
              v_divhinJi(m_ns, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1),
              v_omega_int_en(m_nen-1,0.),
              m_energy_transfer_source(0.0){

    // Set the molar weights..
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

// ----------------------------------------------------------------------

void LarsenTTv::computeHeatFlux()
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

// ----------------------------------------------------------------------

void LarsenTTv::computeJdivJ(const double t) {
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
      m_data.getStepInternalEnthalpy(id_s, hin_ext);

      v_Ji.at(id_s)    = J_ext[0] 
           + (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);

      v_divJi.at(id_s) = (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15); 

      v_divhinJi.at(id_s) = (hin_ext[1]*J_ext[1] - hin_ext[0]*J_ext[0])
                           /(x_ext[1] - x_ext[0] + 1E-15); 
    }
  }

  // =====  Externally given, compute derivatives in stagnation line mode
  else if(m_data.getJdiffFlag().compare("external_stagline_sphere") == 0){

    std::cerr << "DIFFUSION MODEL external_stagline_sphere NOT YET IMPLEMENTED "
              << "FOR THE TTv CASE. USE external_streamwise INSTEAD." << std::endl;

    // THE FOLLOWING IMPLEMENTATION SHOULD BE MORE OR LESS FINE.. EXCEPT THE 
    // v_divhinJi THAT I DID RANDOMLY (not in Munafo's thesis). TO BE CHECKED.
    // 
    //      // Get the radius at current step's integration extremes
    //      m_data.getStepStaglineR(stagline_r_ext);
    //  
    //      // Compute current radius by linear interpolation among the extremes
    //      stagline_r_now = stagline_r_ext[0] +
    //                      (stagline_r_ext[1] - stagline_r_ext[0])/(x_ext[1]
    //                                      - x_ext[0] + 1E-15)*(t - x_ext[0]);
    //  
    //      // For each species..
    //      for(size_t id_s = 0; id_s < m_ns; ++id_s) {
    //  
    //        // Get diffusion mass fluxes at current step's integration extremes
    //        m_data.getStepMassDiffusion(id_s, J_ext);
    //        m_data.getStepInternalEnthalpy(id_s, hin_ext);
    //  
    //        // Compute diffusion mass fluxes by (linear) interpolation among extremes 
    //        v_Ji.at(id_s)    = J_ext[0] 
    //             + (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);
    //  
    //        v_divJi.at(id_s) = (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15) 
    //                          - 2/stagline_r_now*v_Ji.at(id_s); 
    //  
    //        // CHECK THIS:
    //        v_divhinJi.at(id_s) = (hin_ext[1]*J_ext[1] - hin_ext[0]*J_ext[0])
    //                             /(x_ext[1] - x_ext[0] + 1E-15)
    //                          - 2/stagline_r_now*hin_ext[0]*v_Ji.at(id_s); // CHECK ME
    //  
    //      }
  }

  // =====  Computed
  else if(m_data.getJdiffFlag().compare("computed_streamwise") == 0){
    std::cerr << "TO BE IMPLEMENTED.." << std::endl;
  } 

  // =====  Request not recognized
  else {
    std::cerr << "DIFFUSION MODEL NOT SUPPORTED YET! Check the input file" << std::endl;
  }

}

// ----------------------------------------------------------------------

vector_type LarsenTTv::computedxdt(const vector_type& xState, const double t){

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
  for(size_t i_sp = 0; i_sp < m_ns; ++i_sp) {
    v_yi[i_sp] = xState[i_sp];
    if(xState[i_sp] < 0.0){  // It might become negative during some kind of integration..
      v_yi[i_sp] = 0.0;
    }
    v_rhoi[i_sp] = v_yi[i_sp] * rho_now;
  }

  // right after the m_ns mass equations, there are temperatures
  for(size_t i_en = 0; i_en < m_nen; ++i_en) {
    v_T[i_en] = xState[m_ns + i_en];
  }

  // ----  Set mixture state and get properties from the mixture  ----
  // NOTE: &rhoi[0] is a tricky way to pass a VECTOR where an ARRAY was asked..
  m_mix.setState(&v_rhoi[0], &v_T[0], set_state_rhoi_T);

  m_mix.getEnthalpiesMass(&v_hi[0]);
//  m_mix.getEnergiesMass(&v_ei[0]); REMOVE ME
  m_mix.getCpsMass(&v_cpi[0]);
//  m_mix.getCvsMass(&v_cvi[0]); REMOVE ME
  m_mix.netProductionRates(&v_omegai[0]);
  m_mix.energyTransferSource(&v_omega_int_en[0]);

  // ----  Get or compute diffusive mass fluxes and heat flux  ----
  computeJdivJ(t);   // Updates the member variables v_Ji, v_divJi, v_divhinJi
  computeHeatFlux(); // Updates the member variable Q_now

  // ##########   STEP 2) COMPUTING DERIVATIVES   ####################

  // ====> Mass equations
  for(size_t ii = 0; ii < m_ns; ++ii) {
    v_dxdt[ii] = (v_omegai[ii] - v_divJi[ii])/rho_now;
  }

  // ====> Vibrational Temperature
  double fac_num = 0.0;
  double fac_den = 0.0;
  for (int i_sp = 0; i_sp < m_ns; ++i_sp ){
///// OLD OLD //      fac_num += v_omegai[i_sp] * v_ei[i_sp + m_ns];
// NO! Cp not Cv! h not e!!    fac_num += rho_now * v_ei[i_sp + m_ns] * v_dxdt[i_sp] - v_divhinJi[i_sp];
// NO! Cp not Cv! h not e!!    fac_den += v_yi[i_sp] * v_cvi[i_sp + m_ns];
    fac_num += rho_now * v_hi[i_sp + m_ns] * v_dxdt[i_sp] - v_divhinJi[i_sp];
    fac_den += v_yi[i_sp] * v_cpi[i_sp + m_ns];
  }

  v_dxdt[m_ns+1] = (v_omega_int_en[0] - fac_num) / (fac_den * rho_now);

  // ====> Translational Temperature
  double num_part1 = 0.0;
  double num_part2 = 0.0;
  double den       = 0.0;
  for(size_t ii = 0; ii < m_ns; ++ii) {
//    num_part1 += v_hi[ii]*v_omegai[ii]/rho_now;
    num_part1 += v_hi[ii]*v_dxdt[ii];
    num_part2 += v_yi[ii]*v_cpi[ii+m_ns]*v_dxdt[m_ns+1];
    den       += v_yi[ii]*v_cpi[ii];
  }

  // Approximate derivative of u^2
  double Du2Dt = u_now*(pow(u_ext[1],2) - pow(u_ext[0],2))/(x_ext[1] - x_ext[0] + 1E-15);

  // Temperature equation
  v_dxdt[m_ns] = - (0.5*Du2Dt + num_part1 + num_part2 - Q_now)/den;

  // Correcting ALL THE derivatives with the velocity 
  // (remember, these are material derivatives!!)
  for(size_t ii = 0; ii < m_ns+m_nen; ++ii) {
    v_dxdt[ii] /= u_now;
  }

  return v_dxdt;

}; 

