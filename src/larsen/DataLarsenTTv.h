#ifndef DATALARSENTTv_H
#define DATALARSENTTv_H

#include <string>
#include "mutation++.h"

#include "Data.h"
#include "SetupProperties.h"

class DataLarsenTTv: public Data {
public:
    DataLarsenTTv(Mutation::Mixture& l_mix, std::vector< std::string > l_input_file);
    ~DataLarsenTTv();
    
    vector_type getInitialState();
    int nEquations() const { return n_eq; }
   
    void readDataFile();     // reads external data file and fills internal vectors
    void computeEnthalpy();       // computes the enthalpy along the streamline
    void computeMassDiffusion();  // computes diffusive mass fluxes along streamline
    void setCurrentStep(int pos); // sets an internal variable to the value 'pos'

    int  getStrNele(); // Gets number of elements of a streamline

    // The following functions return values at the set step.
    // s_id stands for "species_id"
    void getStepValues(double x_now[]);
    void getStepValues(double x_now[], double u_now[], double rho_now[]);
    void getStepEnthalpy(double H_now[]);
    void getStepInternalEnthalpy(size_t s_id, double hin_now[]);
    void getStepMassDiffusion(size_t s_id, double J_now[]);
    void getStepDiffusionVel(size_t s_id, double ud_now[]);
    void getStepStaglineR(double stagline_r_now[]);
    void getStepStaglineV(double stagline_v_now[]);

    // Those functions return flags (coming from input file)
    void setFlags(std::vector< std::string > l_input_file);
    std::string getHtotDiffusionFlag();
    std::string getJdiffFlag();

private:
    // Functions:
    void errorStateNotSetProperly();
    inline void fillStateVectorX();

    // Variables 
    Mutation::Mixture& m_mix;          // The MIXTURE

    std::string s_ext_flow_filename;   // External flowfield file
    std::string s_species_input_type;  // Flag (from input file)
    std::string s_diffusion_Htot;      // Flag (from input file)
    std::string s_diffusion_J;         // Flag (from input file)
    std::string s_diffusion_J_scheme;  // Flag (from input file)

    const size_t n_sp;
    const size_t n_meq;
    const size_t n_eneq;
    const size_t n_eq;

    std::vector<double> v_rhoi;
    std::vector<double> v_xi;
    std::vector<double> v_yi;
    std::vector<double> v_T;
    std::vector<double> v_Tv;

    vector_type v_X;

    // Those vectors store the whole streamline
    std::vector<double> x_ext, u_ext, rho_ext, T_ext, Tv_ext, H_ext;
    std::vector< std::vector<double> > Y_ext;
    std::vector< std::vector<double> > ud_ext;  // diffusion velocities
    std::vector< std::vector<double> > J_ext;   // diffusive mass fluxes
    std::vector< std::vector<double> > hin_ext; // internal enthalpy of species


    std::vector<double> stagline_r_ext; // r coordinate for stagline approach
    std::vector<double> stagline_v_ext; // v parameter for stagline approach

    // Initial condition
    std::vector<double> Y0_ext;

    // Current step
    int currentStep;

};

#endif /* DATALARSENTTv_H */
