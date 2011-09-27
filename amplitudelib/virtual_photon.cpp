/*
 * Virtual photon wave function
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "virtual_photon.hpp"
#include "../tools/tools.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

const double ZINTACCURACY=0.001;
const int MAXITER_ZINT=1000;
const double MINZ=0.00001;  // Integration limits
const double MAXZ=0.9999;


VirtualPhoton::VirtualPhoton()
{
    // Quark charges, 0=u, 1=d, 2=s
    // Parameters same as in Ref. 0902.1112
    e_f[0]=2.0/3.0; e_f[1]=-1.0/3.0; e_f[2]=2.0/3.0;
    m_f[0]=0.14; m_f[1]=0.14; m_f[2]=0.14;
}

/*
 * Transversially polarized component of the overlap
 */

double VirtualPhoton::PsiSqr_T(double Qsqr, double r, double z)
{
    double result=0;
    for (int f=0; f<=2; f++)     // Sum quar flavors
    {
        double epstmp=Epsilon(Qsqr,z,f);
        result += SQR(e_f[f])*(
            (SQR(z)+SQR(1.0-z))*SQR(epstmp)
                    *SQR(gsl_sf_bessel_K1(epstmp*r))
                + SQR(m_f[f]*gsl_sf_bessel_K0(epstmp*r))
            );
    }
    result *= Nc/(2.0*SQR(M_PI))*ALPHA_e;

    return result;
}

/*
 * Longitudinally polarized component of the overlap
 */

double VirtualPhoton::PsiSqr_L(double Qsqr, double r, double z)
{
    double result=0;
    for (int f=0; f<=2; f++)     // Sum quar flavors
    {
        double epstmp=Epsilon(Qsqr,z,f);
        result += SQR(e_f[f])* SQR( gsl_sf_bessel_K0(epstmp*r) );
    }
    result *= 2.0*Nc/SQR(M_PI)*ALPHA_e*Qsqr*SQR(z)*SQR(1.0-z);

    return result;
}

/* 
 * Wave function overlap integrated over z=[0,1]
 * PsiSqr_T/L is quite a smooth function so there is 
 * nothing tricky to do here, just use gsl_integration_qng
 */
 
/* As we have to integrate a member function of this class by GSL,
 * we need some helper structures
 */
struct zinthelper{
    VirtualPhoton *vm_p;
    double  Qsqr;
    double  r;
};

double  zhelperfuncT(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_T(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

double  zhelperfuncL(double z, void * p){
  return ((zinthelper*)p)->vm_p->PsiSqr_L(
							((zinthelper*)p)->Qsqr,
							((zinthelper*)p)->r,
							z);
}

double VirtualPhoton::PsiSqr_T_intz(double Qsqr, double r)
{
    double result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncT;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, MINZ, MAXZ,  0, ZINTACCURACY, 
        &result, &abserr, &eval);
    //gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    //status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
    //    MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    //gsl_integration_workspace_free(ws);
/*
    if(status){ std::cerr<< "z integral in Photon failed with code " 
        << status << " (transverse, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}
*/  
    return result;
}

double VirtualPhoton::PsiSqr_L_intz(double Qsqr, double r)
{
    double result,abserr;
    size_t eval;
    struct zinthelper zintpar;
    zintpar.vm_p=this;
    zintpar.Qsqr=Qsqr;
    zintpar.r=r;
    gsl_function int_helper;
    int_helper.function=&zhelperfuncL;
    int_helper.params=&zintpar;
    
    int status = gsl_integration_qng(&int_helper, MINZ, MAXZ, 0, ZINTACCURACY, 
        &result, &abserr, &eval);
    //gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_ZINT);
    //int status = gsl_integration_qag(&int_helper, 0, 1, 0, ZINTACCURACY,
    //    MAXITER_ZINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    //gsl_integration_workspace_free(ws);
  /*  
    if(status){ std::cerr<< "z integral in VirtualPhoton failed: code " 
        << status << " (longitudinal, Qsqr=" << Qsqr << ", r=" << r 
        << "relerr=" << abserr/result << ") at " << LINEINFO << std::endl;}
*/
    return result;
}



double VirtualPhoton::Epsilon(double Qsqr, double z, int f)
{
    return std::sqrt( z*(1.0-z)*Qsqr+SQR(m_f[f]) );
}

std::string VirtualPhoton::GetParamString()
{
    std::stringstream str;
    str << "e_f[0]=" << e_f[0] << ", e_f[1]=" << e_f[1] << ", e_f[2]=" << e_f[2]
    << ", m_f[0]=" << m_f[0] << ", m_f[1]=" << m_f[1] << ", m_f[2]=" << m_f[2]
    << ", m_f[3]=" << m_f[3] << endl;
    return str.str();
}

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic)
{
    return os << " Virtual photon wave function. Params: "
        << ic.GetParamString() << " .";
        
}

