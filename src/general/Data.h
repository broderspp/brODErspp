#ifndef DATA_H
#define DATA_H

#include <boost/numeric/odeint.hpp>
#include <vector>
#include <string>

typedef boost::numeric::ublas::vector<double> vector_type;

class Data {
public:
    virtual ~Data(){}
    
    virtual vector_type getInitialState() = 0;
    virtual double mdot() const { }
    virtual int nEquations() const = 0;
    
    virtual std::vector<double> getPartialDensities() const { };
    virtual double getPressure() const {  };
    virtual double getVelocity() const {  };
    virtual double getTTrans() const {  };
    virtual std::vector<double> getTemperatures() const { };
    virtual std::vector<double> getMassFractions() const { };
    
    virtual void setPressure(const double& l_P){  };
    virtual void setTTrans(const double& l_T)  {  };
    virtual void setTemperatures(const std::vector<double>& l_T)   {  };
    virtual void setMassFractions(const std::vector<double>& l_yi) {  };

    virtual std::string getSpeciesOutputFlag() 
                           { errorFunctionCall("getSpeciesOutputFlag()"); };

    // S H O C K I N G
    virtual void setVelocity(const double& l_V) { errorFunctionCall("setVelocity()"); };
    virtual void buildState() { errorFunctionCall("buildState()"); };
    virtual void buildStatePostShock() { errorFunctionCall("buildStatePostShock()"); };

    // L A R S E N
    virtual void readDataFile() 
                           { errorFunctionCall("readDataFile()"); };

    virtual void computeEnthalpy() 
                           { errorFunctionCall("computeEnthalpy()"); };

    virtual void computeMassDiffusion() 
                           { errorFunctionCall("computeMassDiffusion()"); };

    virtual void setCurrentStep(int pos) 
                           { errorFunctionCall("setCurrentStep()"); };

    virtual void getStepValues(double x[], double u[], double rho[]) 
                           { errorFunctionCall("getStepValues()"); };

    virtual void getStepEnthalpy(double H[]) 
                           { errorFunctionCall("getStepEnthalpy()"); };

    virtual void getStepInternalEnthalpy(size_t s_id, double hin[]) 
                           { errorFunctionCall("getStepInternalEnthalpy()"); };

    virtual void getStepMassDiffusion(size_t s_id, double J[]) 
                           { errorFunctionCall("getStepMassDiffusion()"); };

    virtual void getStepDiffusionVel(size_t s_id, double ud[]) 
                           { errorFunctionCall("getStepDiffusionVel()"); };

    virtual void getStepStaglineR(double r[]) 
                           { errorFunctionCall("getStepStaglineR()"); };

    virtual void getStepStaglineV(double v[]) 
                           { errorFunctionCall("getStepStaglineV()"); };

    virtual int  getStrNele() 
                           { errorFunctionCall("getStrNele()"); };

    virtual std::string getHtotDiffusionFlag()
                           { errorFunctionCall("getHtotDiffusionFlag()"); };

    virtual std::string getJdiffFlag() 
                           { errorFunctionCall("getJdiffFlag()"); };

    // M U L T I L A R S E N   (also uses some of the previous)
    virtual int getStreamlinesNumber()
                   { errorFunctionCall("getStreamlinesNumber()"); };

    virtual void getStepXY(size_t str_id, double x[], double y[])
                   { errorFunctionCall("getStepXY()"); };

    virtual void getStepX(double x[]) 
                   { errorFunctionCall("getStepX(double[])");};

    virtual double getStepSlope(size_t str_id)
                   { errorFunctionCall("getStepSlope(size_t)");};

    virtual double getStepDiffU2(size_t str_id)
                   { errorFunctionCall("getStepDiffU2(size_t)");};

    virtual void getStepValues(size_t str_id, double modU[], double rho[]) 
                   { errorFunctionCall("getStepValues(size_t, double[], double[]");};

    virtual void getStepValues(size_t str_id, double x[], double y[], 
                                  double u[], double v[], double rho[])
                   { errorFunctionCall("getStepValues(..with five double[])"); }
  
    virtual void getStepBaselineYi(size_t str_id, size_t s_id, double Yi_now[])
                   { errorFunctionCall("getStepBaselineYi(size_t,size_t,double[][])"); }

    virtual void getStepBaselineT(size_t str_id, double T_now[])
                   { errorFunctionCall("getStepBaselineT(size_t,double[])"); }

    virtual void getStepBaselineTv(size_t str_id, double Tv_now[])
                   { errorFunctionCall("getStepBaselineTv(size_t,double[])"); }

    virtual void getStepEnthalpy(size_t str_id, double H[]) 
                   { errorFunctionCall("getStepEnthalpy(size_t, double[])"); };

    virtual double getPointDensity(size_t str_id, size_t x_id)
                   { errorFunctionCall("getPointDensity(size_t, size_t)"); };

    virtual std::string getGeometryFlag() 
                           { errorFunctionCall("getGeometryFlag()"); };

    // setFlags is PROBABLY USELESS!
    virtual void setFlags(std::vector< std::string > v_str)
                   { errorFunctionCall("setFlags()"); }; // Used to set simulation flags

    virtual void readDataFile(std::string l_filename)
                   { errorFunctionCall("readDataFile(std::string)"); };

    virtual void computeEnthalpy(size_t str_id) 
                   { errorFunctionCall("computeEnthalpy(size_t)"); };

    virtual void computeCurvilinearAbscissa(size_t str_id)
                   { errorFunctionCall("computeCurvilinearAbscissa(size_t)"); };

    virtual void computeSlope(size_t str_id)
                   { errorFunctionCall("computeSlope(size_t)"); };

    virtual void computeVelocityModule(size_t str_id)
                   { errorFunctionCall("computeVelocityModule(size_t)"); };



private:
    void errorFunctionCall(const char * funcName){
      std::cerr << " In: Data.h - " 
                << funcName << " was called, however it is not defined for ";
      std::cerr << "the current problem_type!\n"; 
      exit(1); 
    }
};  

#endif /* DATA_H */
