#ifndef SETUPPROBLEM_H
#define SETUPPROBLEM_H

#include <string>
#include "mutation++.h"

#include "Data.h"
#include "SetupProperties.h"
#include "ShockRelations.h"
#include "Problem.h"
#include "Solution.h"

class SetupProblem {
public:
  SetupProblem(){}
  virtual ~SetupProblem(){}

  // Common for all the solvers
  static SetupProblem* createFactory(std::string l_problem_type, 
                                     std::string l_state_model);

  virtual Problem* getProblem(Mutation::Mixture& l_mix, Data& l_data){};

  virtual Solution* getSolution(Mutation::Mixture& l_mix, Data& l_data){};


  // S H O C K I N G 
  virtual Data* getDataPreShock(Mutation::Mixture& l_mix, 
                        const std::vector< std::string > l_input_parameters){};

  virtual Data* getDataPostShock(Mutation::Mixture& l_mix){};

  virtual ShockRelations* getShockRelations(Mutation::Mixture& l_mix){};

  // L A R S E N    and   M U L T I L A R S E N 
  virtual Data* getData(Mutation::Mixture& l_mix, 
                        const std::vector< std::string > l_input_parameters){};
};

#endif /* SETUPPROBLEM_H */
