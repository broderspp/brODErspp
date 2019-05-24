#ifndef LARSENTTv_H
#define LARSENTTv_H

#include <boost/numeric/odeint.hpp>
#include "mutation++.h"

#include "Problem.h"

class LarsenTTv: public Problem {
public:
    LarsenTTv(Mutation::Mixture& l_mix, Data& l_data);
    ~LarsenTTv(){}

    vector_type computedxdt(const vector_type& x, const double t);

private:
    size_t m_ns, m_nen, m_neq;

    mutable std::vector<double> v_Mw;
    mutable std::vector<double> v_yi, v_rhoi, v_omegai, v_hi, v_cpi, 
                                v_ei, v_cvi, v_omega_int_en, v_T;
    mutable std::vector<double> v_Ji, v_divJi, v_divhinJi;
    mutable vector_type v_dxdt;

    const size_t set_state_rhoi_T;
    const double m_energy_transfer_source;

    // Auxiliary variables.. Initialized here once for all.
    double x_ext[2], u_ext[2], rho_ext[2], H_ext[2], J_ext[2], ud_ext[2], hin_ext[2];
    double u_now, rho_now, T_now, Tv_now, Q_now;

    double stagline_r_ext[2], stagline_v_ext[2]; // just for Stagline
    double stagline_r_now, stagline_v_now;       // just for Stagline

    // Functions
    void computeJdivJ(const double t);
    void computeHeatFlux();
};

#endif /* LARSENTTv_H */
