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
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1){
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

vector_type Larsen1T::computedxdt(const vector_type& xState, const double t){

  double x_ext[2];
  double u_ext[2];
  double rho_ext[2];
  double H_ext[2];
  double J_ext[2];
  double ud_ext[2];
  double stagline_r_ext[2];
  double stagline_v_ext[2];

  double u_now, rho_now, T_now;

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

  // ----  Computing derivatives  -----
  // ====> Mass equations
  double dJdx;
  double J_now, stagline_r_now, stagline_v_now, geom_source;

std::cout << "Ad x = " << t << " ";

  // Simple implementation for the geometry term in stagline
  m_data.getStepStaglineR(stagline_r_ext);
  m_data.getStepStaglineV(stagline_v_ext);

  stagline_r_now = stagline_r_ext[0] + 
                  (stagline_r_ext[1] - stagline_r_ext[0])/(x_ext[1] 
                                  - x_ext[0] + 1E-15)*(t - x_ext[0]);

  stagline_v_now = stagline_v_ext[0] + 
                  (stagline_v_ext[1] - stagline_v_ext[0])/(x_ext[1]  
                                  - x_ext[0] + 1E-15)*(t - x_ext[0]);

std::cout << "x is: " << x_ext[0] << " and r: " << stagline_r_now << " and v: " 
          << stagline_v_now << std::endl;

  for(size_t ii = 0; ii < m_ns; ++ii) {
    m_data.getStepMassDiffusion(ii, J_ext); // Get mass flux
    dJdx = (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15); // Approximate div_x

std::cout << J_ext[0] << " ";
//std::cout << "UNO: " << v_omegai[ii] << " DUE: " << dJdx << std::endl;

    J_now = J_ext[0] + (J_ext[1] - J_ext[0])/(x_ext[1] - x_ext[0] + 1E-15)*(t - x_ext[0]);
    geom_source = 2/stagline_r_now*(v_rhoi[ii]*(+u_now + stagline_v_now) + J_now);
//BOH BOH BOH HO MESSO SEGNI A CASO!

geom_source=0.0;
dJdx = 0.0;
    v_dxdt[ii] = (v_omegai[ii] - dJdx + geom_source)/rho_now;
  }
std::cout << std::endl;

std::cout << "Alla x = " << x_ext[0] << " ";
for(size_t ii = 0; ii < m_ns; ++ii) {
  m_data.getStepDiffusionVel(ii, ud_ext);
  std::cout << " " << ud_ext[0];
}
std::cout << std::endl;

for(size_t ii = 0; ii < m_ns; ++ii) {
  m_data.getStepMassDiffusion(ii, J_ext);
  m_data.getStepDiffusionVel(ii, ud_ext);
  std::cout << "COMPARA - J_ext: " << J_ext[0] << "  rhoi*ud: " << v_rhoi[ii]*ud_ext[0] << "  rhoi: " << v_rhoi[ii] << "  ud_ext: " << ud_ext[0] << std::endl;
}
std::cout << "------------------------------------\n";

std::cout << "Density at x = " << x_ext[0] << " ";
for(size_t ii = 0; ii < m_ns; ++ii) {
  std::cout << v_rhoi[ii] << " ";
}

std::cout << "\n------------------------------------\n";






  // ====> Translational temperature
  double num_part = 0.0;
  double den      = 0.0;
  for(size_t ii = 0; ii < m_ns; ++ii) {
    num_part += v_hi[ii]*v_omegai[ii]/rho_now;
    den      += v_yi[ii]*v_cpi[ii];
  }
  // Approximate derivative of u^2
  double Du2Dt = u_now*(pow(u_ext[1],2) - pow(u_ext[0],2))/(x_ext[1] - x_ext[0] + 1E-15);

  // Approximate derivative of enthalpy, AKA thermal flux
  m_data.getStepEnthalpy(H_ext);
  double Q_ext = u_now*(H_ext[1] - H_ext[0])/(x_ext[1] - x_ext[0] + 1E-15);

// ..
//std::cout << "TOTAL ENTHALPY at x = " << x_ext[0] << "   " << H_ext[0] << std::endl;
//std::cout << "Internal Enthalpy at x = " << x_ext[0] << "   " << m_mix.mixtureHMass() << std::endl;
Q_ext = 0.0;
  v_dxdt[m_ns] = - (0.5*Du2Dt + num_part - Q_ext)/den;

  // Correcting ALL THE derivatives with the velocity (remember, these are material derivatives!!)
  for(size_t ii = 0; ii < m_ns+ 1; ++ii) {
    v_dxdt[ii] /= u_now;
  }

  return v_dxdt;

};

