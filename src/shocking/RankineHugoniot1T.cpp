#include "RankineHugoniot1T.h"

typedef boost::numeric::ublas::vector< double > vec_type;
typedef boost::numeric::ublas::matrix< double > mat_type;
typedef boost::numeric::ublas::permutation_matrix<size_t> pmatrix;

RankineHugoniot1T::RankineHugoniot1T(Mutation::Mixture& l_mix)
                                        : m_mix(l_mix),
                                          m_ns(m_mix.nSpecies()),
                                          m_nen(m_mix.nEnergyEqns()),
                                          set_state_rhoi_T(1) {
  // Resize Eigen vectors
  v_rhoi1.resize(m_ns);
  v_y.resize(m_ns);
  v_x.resize(m_ns);
}

RankineHugoniot1T::~RankineHugoniot1T(){}

void RankineHugoniot1T::applyShockRelations(const Data& l_data_before, Data& l_data_after){

    double l_u1 = l_data_before.getVelocity();
    double l_p1 = l_data_before.getPressure();
    double l_T1 = l_data_before.getTTrans();

    double mdot;
    double l_rho1, l_rho2;
    double h1, h2;

    for(size_t s_id = 0; s_id < m_ns; ++s_id) { // fill v_rhoi1
      v_rhoi1(s_id) = l_data_before.getPartialDensities().at(s_id);
    }

    // Set mixture state, pre-shock
    m_mix.setState(v_rhoi1.data(), &l_T1, set_state_rhoi_T);

    // Compute mass fractions vector v_y
    l_rho1 = v_rhoi1.sum();
    v_y = v_rhoi1/l_rho1;

    // Iteratively solve the problem
    double sol[3]; // iterative solution
    findPostShockState(l_u1, l_p1, l_T1, l_rho1, &v_y[0], &m_mix, sol);

    // Save solution into post-shock data structure
    l_data_after.setVelocity(sol[0]);
    l_data_after.setPressure(sol[1]);
    l_data_after.setTTrans(sol[2]);
    l_data_after.setMassFractions(l_data_before.getMassFractions());

    l_data_after.buildStatePostShock();
}

// ---------------------------------------------------------------------------


void RankineHugoniot1T::findPostShockState(double Us, double P1, double T1,double RHO1,
                              double Yi[] ,Mutation::Mixture* gasMix, double sol[])
{  
    // This function is from the C++ shocking version, by Pierre Schrooyen, 2014.
    //
    // The function Takes the free-stream parameters and computes the post-shock state,
    // iterating in order to find the right tempertaure, pressure and velocity.
    // In fact, as the gas mixture gets hot across the shock, internal modes get excited.

    double U1 = Us;
    size_t ns = gasMix->nSpecies();
    double rho_i [ns];
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = RHO1*Yi[iSpecies];
    }
    gasMix->setState(rho_i,&T1,1);
    
    std::cout<<"Compute jump conditions..."<<std::endl;
    std::cout<<"---------------------------------------------------------------------------"<<std::endl;
    std::cout <<" Pre-shock conditions (given)"<<std::endl;
    std::cout <<"P1 : "<<std::setw(15)<< P1<<" [Pa]"<<std::endl;
    std::cout <<"T1 : "<<std::setw(15)<< T1<<" [K]"<<std::endl;
    std::cout <<"US : "<<std::setw(15)<< U1<<" [m/s]"<<std::endl;
    std::cout<<"---------------------------------------------------------------------------"<<std::endl;
    double Mm1 = gasMix->mixtureMwMass(Yi);
    double Rstar = Mutation::RU/Mm1;
    double Gamma = gasMix->mixtureFrozenCpMass()/gasMix->mixtureFrozenCvMass();
    
    double GM1 = Gamma-1.0;
    double GP1 = Gamma+1.0;
    
    double C1 = std::sqrt(Gamma*Rstar*T1);
    double Ms = U1/C1;
    
    double P2 = P1 *(1.0+2.0*Gamma/GP1*(Ms*Ms-1.0));
    double T2 = T1 *(1.0+((2.0*GM1/(GP1*GP1))*((Gamma*Ms*Ms+1.0)/(Ms*Ms))*(Ms*Ms-1.0)));
    double U2 = C1 *2.0 /GP1 *(Ms - 1.0/Ms);
    double RHO2 = P2/(Rstar*T2);
    
    std::cout <<" Post-shock conditions (cold gas approximation)"<<std::endl;
    std::cout <<"P2 :    "<<std::setw(15)<< P2<<" [Pa]"<<std::endl;
    std::cout <<"T2 :    "<<std::setw(15)<< T2<<" [K]"<<std::endl;
    std::cout <<"US-U2 : "<<std::setw(15)<< Us-U2<<" [m/s]"<<std::endl;
    std::cout <<"RHO2  : "<<std::setw(15)<< RHO2<<" [kg/m^3]"<<std::endl;

    std::cout<<"---------------------------------------------------------------------------"<<std::endl;
    
    //Hot gas approximation
    double rho1  = RHO1;
    double rhou1 = RHO1*U1;
    double h1    = gasMix->mixtureHMass();
    double Ekin1 = 0.5*U1*U1;
    int niter = 0;
    double TP_temp[2];
    
    vec_type F (3,0.0);
    mat_type J(3,3);
    double error = 1.0;
    vec_type xi(3,0.0);
    
    //Initial guess :
    xi[0] = Us-U2;
    xi[1] = P2;
    xi[2] = T2;
    double Ui = xi[0];
    double Pi = xi[1];
    double Ti = xi[2];
    double hi = 0.0;
    hi = gasMix->mixtureHMass();
    for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
        rho_i[iSpecies] = Yi[iSpecies]*Pi*Mm1/(Mutation::RU*Ti);
    }
    gasMix->setState(rho_i,&Ti,1);
    F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
    F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
    F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
    double Resini = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
    
    while (std::abs(error)>1.0e-12 && niter<100) {
        Ui = xi[0];
        Pi = xi[1];
        Ti = xi[2];
        TP_temp[0] = Ti;
        TP_temp[1] = Pi;
        double rho = Pi*Mm1/(Mutation::RU*Ti);
        for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
            rho_i[iSpecies] = Yi[iSpecies]*rho;
        }
        gasMix->setState(rho_i,&Ti,1);
        hi = gasMix->mixtureHMass();
        
        F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
        F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
        F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
        
        J(0,0) = -Pi/(Rstar*Ti);
        J(0,1) = -Ui/(Rstar*Ti);
        J(0,2) = (Ui*Pi)/(Rstar*Ti*Ti);
        J(1,0) = -(2.0*Ui*Pi)/(Rstar*Ti);
        J(1,1) = -1.0-(Ui*Ui)/(Rstar*Ti);
        J(1,2) = (Ui*Ui*Pi)/(Rstar*Ti*Ti);
        J(2,0) = - Ui;
        J(2,1) = 0.0;
        J(2,2) = - gasMix->mixtureFrozenCpMass();
        
        pmatrix piv(3);
        lu_factorize(J, piv);
        lu_substitute(J, piv, F);
        xi[0] += F[0];
        xi[1] += F[1];
        xi[2] += F[2];
        
        //Compute Error
        Ui = xi[0];
        Pi = xi[1];
        Ti = xi[2];
        TP_temp[0] = Ti;
        TP_temp[1] = Pi;
        rho = Pi*Mm1/(Mutation::RU*Ti);
        for (size_t iSpecies=0; iSpecies<ns; iSpecies++) {
            rho_i[iSpecies] = Yi[iSpecies]*rho;
        }
        gasMix->setState(rho_i,&Ti,1);
        hi = gasMix->mixtureHMass();
        
        F[0] = -(rhou1 - Pi/(Rstar*Ti)*Ui);
        F[1] = -(P1 - Pi - Pi/(Rstar*Ti)*(Ui*Ui)+rhou1*U1);
        F[2] = -(h1 + Ekin1 - hi - 0.5*Ui*Ui);
        
        error = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2])/Resini;
        niter++;
    }
    if (niter>=99) {
        std::cout<<"Newton-Raphson for the post-shock value did not converge, the program will terminate :"<<error<<std::endl;
        exit(0);
    }
    else{
        std::cout<<"Newton-Raphson converged in "<< niter<< " iterations and the residual are "<<error<<std::endl;
    }
    sol[0] = xi[0];
    sol[1] = xi[1];
    sol[2] = xi[2];
    double c2 = gasMix->frozenSoundSpeed();
    double Mach2 =(Us-U2)/c2;
    std::cout <<" Post-shock conditions (hot gas approximation)"<<std::endl;
    std::cout <<"P2 :    "<<std::setw(15)<< sol[1]<<" [Pa]"<<std::endl;
    std::cout <<"T2 :    "<<std::setw(15)<< sol[2]<<" [K]"<<std::endl;
    std::cout <<"US-U2 : "<<std::setw(15)<< sol[0]<<" [m/s]"<<std::endl;
        std::cout<<"Mach number"<<std::setw(15)<< Mach2<<" [-]"<<std::endl;
    std::cout<<"---------------------------------------------------------------------------"<<std::endl;
}
