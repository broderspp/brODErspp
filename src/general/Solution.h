#ifndef SOLUTION_H 
#define SOLUTION_H

#include "mutation++.h"
#include "Data.h"

class Solution {
// This class is used to store the solution. Also provides some exporting functions

public:
    Solution(Mutation::Mixture& l_mix, Data& l_data)
            : m_mix(l_mix),
              m_data(l_data){}
    virtual ~Solution(){}


    // For MULTILLARSEN
    virtual void initSol() { errorFunctionCall("initSol()"); };

    virtual void saveStepSolution(size_t step_id, const vector_type& x_sol)
                               { errorFunctionCall("saveStepSolution()"); };

    virtual void exportBaselineSolVTK(const char* filename) 
                               { errorFunctionCall("exportBaselineSolVTK()"); };

    virtual void exportSolVTK(const char* filename) 
                               { errorFunctionCall("exportSolVTK()"); };

    virtual void exportAllVTK(const char* filename) 
                               { errorFunctionCall("exportAllVTK()"); };

protected:
    Mutation::Mixture& m_mix;
    Data& m_data;

private:
    void errorFunctionCall(const char * funcName){
      std::cerr << " In: Solution.h - " 
                << funcName << " was called, however it is not defined for ";
      std::cerr << "the current problem_type!\n"; 
      exit(1); 
    }

};

#endif /* SOLUTION_H */
