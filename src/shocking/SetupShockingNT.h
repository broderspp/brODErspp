#ifndef SETUPSHOCKINGNT_H
#define SETUPSHOCKINGNT_H

#include "SetupProblem.h"

class SetupShockingNT: public SetupProblem {
public:
    SetupShockingNT();
    ~SetupShockingNT(){}

    Data* getDataPreShock(Mutation::Mixture& l_mix, 
                          const std::vector< std::string > l_input_file);

    Data* getDataPostShock(Mutation::Mixture& l_mix);

    ShockRelations* getShockRelations(Mutation::Mixture& l_mix);

    Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data);
};

#endif /* SETUPSHOCKINGNT_H */
