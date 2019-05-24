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

    Eigen::VectorXd v_y;
    Eigen::VectorXd v_x;
    Eigen::VectorXd v_rhoi1;
    Eigen::VectorXd v_rhoi2;
    Eigen::VectorXd v_hi1;
    Eigen::VectorXd v_hi2;
    Eigen::VectorXd v_cpi1;
    Eigen::VectorXd v_cpi2;

    const double m_tol;
};

#endif /* RANKINEHUGONIOT1T_H */
