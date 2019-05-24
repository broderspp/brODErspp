#include "ShockingTTv.h"

ShockingTTv::ShockingTTv(Mutation::Mixture& l_mix, Data& l_data)
            : Problem(l_mix, l_data),
              m_ns(l_mix.nSpecies()),
              m_nmom(1),
              m_nen(l_mix.nEnergyEqns()),
              m_neq(m_ns+m_nmom+m_nen),
              pos_u(m_ns),
              pos_T(m_ns+m_nmom),
              pos_T_int(pos_T+1),
              m_mdot(0.),
              v_Mw(m_ns, 0.),
              v_yi(m_ns, 0.),
              v_rhoi(m_ns, 0.),
              v_omegai(m_ns, 0.),
              v_omega_int_en(m_nen-1, 0.),
              v_T(m_nen, 0.),
              v_hi(m_ns*m_nen, 0.),
              v_ei(m_ns*m_nen, 0.),
              v_cpi(m_ns*m_nen, 0.),
              v_cvi(m_ns*m_nen, 0.),
              v_dxdt(m_neq, 0.),
              set_state_rhoi_T(1),
              m_energy_transfer_source(0.0){
    for (int i_sp = 0; i_sp < m_ns; ++i_sp){
        v_Mw[i_sp] = m_mix.speciesMw(i_sp);
    }
}

vector_type ShockingTTv::computedxdt(const vector_type& x, const double t){
// This equation computes the RHS of the system, required by the integration routine.
// See equations in the doc folder.

    m_mdot = m_data.mdot();
    double l_u = x[pos_u];
    double l_rho = m_mdot/l_u;

    // Compute species densities
    for (int i_sp = 0; i_sp < m_ns; i_sp++){

        if(x[i_sp] >= 0.0 ) { ////////////////////
          v_yi[i_sp] = x[i_sp];
        } else {              //////////////////// 
          v_yi[i_sp] = 0.0;   ////////////////////
        }                     ////////////////////

        v_rhoi[i_sp] = v_yi[i_sp] * l_rho;
    }
    // Fill array of temperatures
    for (int i_en = 0; i_en < m_nen; ++i_en){
        v_T[i_en] = x[pos_T+i_en];
    }

    // Get mixture quantities at the current thermodynamic state
    m_mix.setState(&v_rhoi[0], &v_T[0], set_state_rhoi_T);
    m_mix.getEnthalpiesMass(&v_hi[0]);
    m_mix.getEnergiesMass(&v_ei[0]);
    m_mix.getCpsMass(&v_cpi[0]);
    m_mix.getCvsMass(&v_cvi[0]);
    m_mix.netProductionRates(&v_omegai[0]);
    m_mix.energyTransferSource(&v_omega_int_en[0]);

    // Working variables
    double y_sum     = 0.0;
    double omega_sum = 0.0;
    double cp_tr_sum = 0.0; // roto-translational cp
    double cp_in_sum = 0.0; // internal cp
    double cv_in_sum = 0.0; // internal cv
    double h_tot_sum = 0.0;
    double e_in_sum  = 0.0;

    for (int i_sp = 0; i_sp < m_ns; i_sp++){
      y_sum     += v_yi[i_sp]/v_Mw[i_sp];
      omega_sum += v_omegai[i_sp]/v_Mw[i_sp];
      cp_tr_sum += v_cpi[i_sp]*v_yi[i_sp];
      cp_in_sum += v_cpi[i_sp+m_ns]*v_yi[i_sp];
      cv_in_sum += v_cvi[i_sp+m_ns]*v_yi[i_sp];
      h_tot_sum += v_hi[i_sp]*v_omegai[i_sp];
      e_in_sum  += v_ei[i_sp+m_ns]*v_omegai[i_sp];
    }

    // Matrix elements for momentum equation (note the correction for free electrons)
    double ele_a = m_mdot*l_u/(Mutation::RU*v_T[0]) - l_rho*y_sum 
                   + m_mix.hasElectrons()*l_rho*(1 - v_T[1]/v_T[0])*v_yi[0]/v_Mw[0];
    double ele_b = m_mdot/v_T[0]*(y_sum - m_mix.hasElectrons()*v_yi[0]/v_Mw[0]);
    double ele_c = m_mdot/v_T[0]*m_mix.hasElectrons()*v_yi[0]/v_Mw[0];
    double ele_d = -omega_sum + m_mix.hasElectrons()*v_omegai[0]/v_Mw[0]*(1 - v_T[1]/v_T[0]);

    // Matrix elements for total energy equation
    double ele_e = l_u;
    double ele_f = cp_tr_sum;
    double ele_g = cp_in_sum;
    double ele_h = m_energy_transfer_source/m_mdot - h_tot_sum/m_mdot; // CHECK RADIATION!!!! 

    // Matrix elements for internal energy equation
    double ele_i = m_mix.hasElectrons()*v_rhoi[0]*Mutation::RU/v_Mw[0]*v_T[1]; // electrons pressure
    double ele_j = 0.0;
    double ele_k = m_mdot*cv_in_sum;
    double ele_l = v_omega_int_en[0] - e_in_sum;

    // Assemble matrix and solve system
    Eigen::Matrix3d A;
    A << ele_a, ele_b, ele_c,
         ele_e, ele_f, ele_g,
         ele_i, ele_j, ele_k;

    Eigen::Vector3d b;
    b << ele_d, ele_h, ele_l;

    Eigen::Matrix3d iA = invert_3by3_matrix(A);
    Eigen::Vector3d RHS = iA*b;

    // Filling the derivatives vector
    for (int i_sp = 0; i_sp < m_ns; i_sp++){
      v_dxdt[i_sp] = v_omegai[i_sp]/m_mdot;
    }

    v_dxdt[pos_u]   = RHS(0);
    v_dxdt[pos_T]   = RHS(1);
    v_dxdt[pos_T+1] = RHS(2);


// Old implementation, missing electrons //    // Sum Variables
// Old implementation, missing electrons //    double wrk1, wrk2, wrk3;
// Old implementation, missing electrons //    double h_sum =0.0;
// Old implementation, missing electrons //    double cp_sum = 0.0;
// Old implementation, missing electrons //    double y_sum = 0.0;
// Old implementation, missing electrons //    double omega_sum = 0.0;
// Old implementation, missing electrons //
// Old implementation, missing electrons //    for (int i_sp = 0; i_sp < m_ns; i_sp++){
// Old implementation, missing electrons //        wrk1 = 1.0/v_Mw[i_sp];
// Old implementation, missing electrons //        wrk2 = v_yi[i_sp];
// Old implementation, missing electrons //        wrk3 = v_omegai[i_sp];
// Old implementation, missing electrons //
// Old implementation, missing electrons //        h_sum += v_hi[i_sp]*wrk3;
// Old implementation, missing electrons //        cp_sum += wrk2 * v_cpi[i_sp];
// Old implementation, missing electrons //        y_sum += wrk2 * wrk1;
// Old implementation, missing electrons //        omega_sum += wrk3 * wrk1;
// Old implementation, missing electrons //    }
// Old implementation, missing electrons //
// Old implementation, missing electrons //    // Helper Variables;
// Old implementation, missing electrons //    double a, b, c, d, e, f;
// Old implementation, missing electrons //    a = m_mdot * l_u / (Mutation::RU * v_T[0]) - l_rho * y_sum;
// Old implementation, missing electrons //    b = m_mdot / v_T[0] * y_sum;
// Old implementation, missing electrons //    c = - omega_sum;
// Old implementation, missing electrons //    d = m_mdot * l_u;
// Old implementation, missing electrons //    e = m_mdot * cp_sum;
// Old implementation, missing electrons //    f = -h_sum + m_mdot * m_energy_transfer_source; // This is for radiation
// Old implementation, missing electrons //
// Old implementation, missing electrons //    for (int i_en = 0; i_en < m_nen-1; i_en++){
// Old implementation, missing electrons //        f = f - v_omega_int_en[i_en];
// Old implementation, missing electrons //    }
// Old implementation, missing electrons //
// Old implementation, missing electrons //    wrk1 = 1.0/m_mdot;
// Old implementation, missing electrons //    wrk2 = 1.0/(a*e - b*d);
// Old implementation, missing electrons //
// Old implementation, missing electrons //    // Filling the derivatives vector
// Old implementation, missing electrons //    for (int i_sp = 0; i_sp < m_ns; i_sp++){
// Old implementation, missing electrons //        v_dxdt[i_sp] = v_omegai[i_sp] * wrk1;
// Old implementation, missing electrons //    }
// Old implementation, missing electrons //
// Old implementation, missing electrons //    v_dxdt[pos_u] = (c*e - b*f) * wrk2;
// Old implementation, missing electrons //    v_dxdt[pos_T] = (a*f - c*d) * wrk2;
// Old implementation, missing electrons //
// Old implementation, missing electrons //    // This is the part that need to be generalized in order to be more than TTv
// Old implementation, missing electrons //    // Probably a function in Thermodynamics should be added that computes w*h
// Old implementation, missing electrons //    // Fix it for electrons
// Old implementation, missing electrons //    double fac_num = 0.;
// Old implementation, missing electrons //    double fac_den = 0.;
// Old implementation, missing electrons //    for (int i_sp = 0; i_sp < m_ns; ++i_sp ){
// Old implementation, missing electrons //        fac_num += v_omegai[i_sp] * v_ei[i_sp + m_ns];
// Old implementation, missing electrons //        fac_den += v_yi[i_sp] * v_cvi[i_sp + m_ns];
// Old implementation, missing electrons //    }
// Old implementation, missing electrons //
// Old implementation, missing electrons //    for (int i_en = 0; i_en < m_nen - 1; ++i_en){
// Old implementation, missing electrons //        v_dxdt[pos_T_int + i_en] = (v_omega_int_en[i_en] - fac_num) / fac_den * wrk1;
// Old implementation, missing electrons //    }

    return v_dxdt;
}

// -------------------------------------------------------------

Eigen::Matrix3d ShockingTTv::invert_3by3_matrix(Eigen::Matrix3d A)
{
  // This function inverts analytically the 3x3 matrix A, into the matrix iA
  double det;

  det = A(0,0)*(A(2,2)*A(1,1) - A(2,1)*A(1,2))
      - A(1,0)*(A(2,2)*A(0,1) - A(2,1)*A(0,2))
      + A(2,0)*(A(1,2)*A(0,1) - A(1,1)*A(0,2));

  det = 1.0/det;

  Eigen::Matrix3d iA; // Inverse matrix

  // First column 
  iA(0,0) = det*(A(2,2)*A(1,1) - A(1,2)*A(2,1));
  iA(1,0) = - det*(A(2,2)*A(1,0) - A(2,0)*A(1,2));
  iA(2,0) = det*(A(2,1)*A(1,0) - A(2,0)*A(1,1));

  // Second column 
  iA(0,1) = - det*(A(2,2)*A(0,1) - A(2,1)*A(0,2));
  iA(1,1) = det*(A(2,2)*A(0,0) - A(2,0)*A(0,2));
  iA(2,1) = - det*(A(2,1)*A(0,0) - A(2,0)*A(0,1));

  // Third column
  iA(0,2) = det*(A(1,2)*A(0,1) - A(1,1)*A(0,2));
  iA(1,2) = - det*(A(1,2)*A(0,0) - A(1,0)*A(0,2));
  iA(2,2) = det*(A(1,1)*A(0,0) - A(1,0)*A(0,1));

// std::cerr<< "DBDBDB   INVERTING: " << std::endl 
//          << iA(0,0) << "  " << iA(0,1) << "  " << iA(0,2) << std::endl
//          << iA(1,0) << "  " << iA(1,1) << "  " << iA(1,2) << std::endl
//          << iA(2,0) << "  " << iA(2,1) << "  " << iA(2,2) << std::endl;
//  
  return iA;
}

