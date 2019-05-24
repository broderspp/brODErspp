#include "RankineHugoniot1T.h"

RankineHugoniot1T::RankineHugoniot1T(Mutation::Mixture& l_mix)
                                        : m_mix(l_mix),
                                          set_state_rhoi_T(1),
                                          v_y(m_mix.nSpecies(),0.0),
                                          v_x(m_mix.nSpecies(),0.0),
                                          v_rhoi1(m_mix.nSpecies(),0.0),
                                          v_rhoi2(m_mix.nSpecies(),0.0),
                                          v_hi1(m_mix.nSpecies(),0.0),
                                          v_hi2(m_mix.nSpecies(),0.0),
                                          v_cpi1(m_mix.nSpecies(),0.0),
                                          v_cpi2(m_mix.nSpecies(),0.0),
                                          m_tol(1.0E-6){}

RankineHugoniot1T::~RankineHugoniot1T(){}

void RankineHugoniot1T::applyShockRelations(const Data& l_data_before, Data& l_data_after){

    double l_u1 = l_data_before.getVelocity();
    double l_p1 = l_data_before.getPressure();
    double l_T1 = l_data_before.getTTrans();

    double l_u2;
    double l_p2;
    double l_T2;

    double l_rho1, l_rho2;
    double m_dot;

    double h1, h2;

    for(size_t s_id = 0; s_id < m_mix.nSpecies(); ++s_id) {
      v_rhoi1(s_id) = l_data_before.getPartialDensities().at(s_id);
    }

    m_mix.setState(v_rhoi1.data(), &l_T1, set_state_rhoi_T);

    l_rho1 = v_rhoi1.sum();
    v_y = v_rhoi1 / l_rho1;
    // v_x = l_data_before.getMoleFrac();
    m_mix.convert<Mutation::Thermodynamics::Y_TO_X>(v_y.data(), v_x.data());

    m_dot = l_data_before.mdot();

    m_mix.getEnthalpiesMass(v_hi1.data());
    m_mix.getCpsMass(v_cpi1.data());

    h1 = v_y.dot(v_hi1);

    double ratio, resR;
    ratio = 0.1;
    resR  = 1.0;
    l_T2  = l_T1;

    while (resR > m_tol){

        l_p2 = l_p1 + m_dot * l_u1 * (1.0 - ratio);
        l_rho2 = l_rho1 / ratio;

        double rhs = h1 + 0.5 * l_u1 * l_u1 * (1.0 - ratio * ratio);
        double resT = 1.;

        while (resT > m_tol){
            double f, fp;

            m_mix.setState(v_rhoi2.data(), &l_T2, set_state_rhoi_T);

            m_mix.getEnthalpiesMass(v_hi2.data());
            m_mix.getCpsMass(v_cpi2.data());

            f = v_y.dot(v_hi2);
            fp = v_y.dot(v_cpi2);

            double T_old = l_T2;
            f -= rhs;
            l_T2 -= f/fp;
            resT = std::abs( (l_T2 - T_old) / T_old );

        }

        double ratio_old = ratio;

        l_rho2 = m_mix.density(l_T2, l_p2, v_x.data());

        ratio = l_rho1/l_rho2;
        resR = std::( (ratio -ratio_old)/ratio );

    }

    l_u2 = l_u1 * ratio;

    l_data_after.setPressure(l_p2);
    l_data_after.setVelocity(l_u2);
    l_data_after.setTTrans(l_T2);
    l_data_after.setMassFractions(l_data_before.getMassFractions());

    l_data_after.buildStatePostShock();

}

// void RankineHugoniot1T::applyShockRelations(const Data& l_data_before, Data& l_data_after){
//
//     double l_u1 = l_data_before.getVelocity();
//     double l_p1 = l_data_before.getPressure();
//     double l_T1 = l_data_before.getTTrans();
//
//     // Set the state with rhoi, T
//     m_mix.setState(&l_data_before.getPartialDensities()[0],
//                    &l_T1,
//                    set_state_rhoi_T);
//
//     double l_gamma =  m_mix.mixtureFrozenCpMass()/m_mix.mixtureFrozenCvMass();
//     double l_c1 = m_mix.frozenSoundSpeed();
//
//     double l_gp1 = l_gamma+1.;
//     double l_gm1 = l_gamma-1.;
//     double l_M1 = l_u1/l_c1;
//     double l_M1s = l_M1 * l_M1;
//
//     l_data_after.setPressure(l_p1 * (2.0 * l_gamma * l_M1s - l_gm1) / l_gp1);
//     l_data_after.setVelocity(l_u1 - l_c1 * 2.0 / l_gp1 * (l_M1 - 1.0 / l_M1));
//     l_data_after.setTTrans(l_T1 * (2.0 * l_gamma * l_M1s - l_gamma + 1) * (l_gm1 + 2.0 / l_M1s) / (l_gp1 * l_gp1) );
//     l_data_after.setMassFractions(l_data_before.getMassFractions());
//
//     l_data_after.buildStatePostShock();
//
// }
