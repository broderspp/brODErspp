#ifndef SETUPSHOCKING1T_H
#define SETUPSHOCKING1T_H

#include "SetupProblem.h"

class SetupShocking1T: public SetupProblem {
public:
    SetupShocking1T();
    ~SetupShocking1T(){}

    Data* getDataPreShock(Mutation::Mixture& l_mix,
                          const std::vector< std::string > l_input_file);

    Data* getDataPostShock(Mutation::Mixture& l_mix);

    ShockRelations* getShockRelations(Mutation::Mixture& l_mix);

    Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data);
};

#endif /* SETUPSHOCKING1T_H */
