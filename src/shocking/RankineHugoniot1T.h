#ifndef RANKINEHUGONIOT1T_H
#define RANKINEHUGONIOT1T_H

#include "mutation++.h"

#include "ShockRelations.h"

class RankineHugoniot1T: public ShockRelations {
public:
    RankineHugoniot1T(Mutation::Mixture& l_mix);
    ~RankineHugoniot1T();

    void applyShockRelations(const Data& l_data_pre, Data& l_data_post);

private:
    Mutation::Mixture& m_mix;
    const size_t set_state_rhoi_T;
    const size_t m_ns, m_nen;

    Eigen::VectorXd v_rhoi1;
    Eigen::VectorXd v_y;
    Eigen::VectorXd v_x;

    void findPostShockState(double Us, double P1, double T1,double RHO1,
                            double Yi[] ,Mutation::Mixture* gasMix, double sol[]);

};

#endif /* RANKINEHUGONIOT1T_H */
