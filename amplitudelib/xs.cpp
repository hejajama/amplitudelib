/*
 * AmplitudeLib cross section calculation methods
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

#include "amplitudelib.hpp"
#include "../tools/tools.hpp"
#include <cmath>
#include "../pdf/pdf.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../fragmentation/kkp.hpp"
#include "../fragmentation/pkhff.hpp"
#include "../pdf/cteq.hpp"
#include <gsl/gsl_integration.h>

/* Differential forward hadron production multiplicity
 * dN_h / (dy_h d^2 p_T)
 * Ref 1001.1378 eq (1), normalization is arbitrary
 */
struct Inthelper_hadronprod
{   
    AmplitudeLib *N;
    double y, pt, xf;
    double sqrts;
    double miny,maxy;
    double minpt,maxpt;
    PDF* pdf;
    FragmentationFunction* frag;
    bool deuteron;
    Hadron final;
    bool ptint;     // if false, don't integrate over pt
    double scale;
};

double Inthelperf_hadronprod(double z, void *p)
{
    // x1: fraction of parent hadron momentum carried by quark,
    // x2: Bjorken x for the dense system
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    double x1 = par->xf/z; double x2 = x1*std::exp(-2.0*par->y);
    if (x1>1)
    {
		cerr << "Entering kinematically fobidden region at y=" << par->y <<", pt=" << par->pt << " " << LINEINFO << endl;
		return 0;
	}
	
	double scale=par->scale;
	if (scale<0) scale = par->pt;

    double y_A = std::log(par->N->X0()/x2);

    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", z " << z << " xf " <<
            par->xf << " x1 " << x1 << " x2 " << x2 << " y_A " << y_A 
            << " y " << par->y << " sqrts " 
            << par->sqrts << " pt " << par->pt << endl ;
        return 0;
    }

    bool deuteron = par->deuteron;
    par->N->InitializeInterpolation(y_A);

    double result = 0;
    
    // Quark from proton:
    // UGD in fundamental representation
    double nf = par->N->S_k(par->pt/z, y_A);
    // PDF and fragmentation
    double xqf = par->pdf->xq(x1, scale, U)*par->frag->Evaluate(U, par->final, z, scale)
        + par->pdf->xq(x1, scale, D)*par->frag->Evaluate(D, par->final, z, scale);
        //+ par->pdf->xq(x1, scale, S)*par->frag->Evaluate(S, par->final, z, scale);

    if (deuteron)
    {
        // isospin symmetry, u in p -> d in n
        xqf += par->pdf->xq(x1, scale, U)*par->frag->Evaluate(D, par->final, z, scale)
        + par->pdf->xq(x1, scale, D)*par->frag->Evaluate(U, par->final, z, scale)
        ;//+ par->pdf->xq(x1, scale, S)*par->frag->Evaluate(S, par->final, z, scale);
    }
        
    result = nf*xqf;
	
    // Adjoint representation, gluon scatters
    double na = par->N->S_k(par->pt/z, y_A, true);
    double xgf = par->pdf->xq(x1, scale, G)*par->frag->Evaluate(G, par->final, z, scale);
    if (deuteron) xgf *= 2.0;   // gluon pdf gets multiplied by 2
    result += na*xgf;

    return result/SQR(z);
}

/*
 * Calculate dN / dyd^2 pt
 * if deuteron is true (default: false), the probe is deuteron, not proton
 * => use isospin symmetry
 * Default final state particle is PI0
 * If optional parameter scale is given (default: -1), use it when evaluating PDF/FF
 * (used e.g. when calculating DPS contribution (b))
 */
double AmplitudeLib::dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
    FragmentationFunction* fragfun, PDF* pdf, bool deuteron, Hadron final, double scale )
{
    const double K = 1.0; // normalization factor
    // We assume light hadrons
    double xf = pt/sqrts*std::exp(y);
    
    SetOutOfRangeErrors(false);
    
    if (xf > 1 or sqrts < 10)
    {
        cerr << "Parameters don't make sense, xf=" << xf << ", sqrts="
            << sqrts << ", y=" << y << " pt=" << pt << " " << LINEINFO << endl;
        return 0;
    }
    Inthelper_hadronprod helper;
    helper.N=this; helper.y=y; helper.pt=pt; helper.xf=xf;
    helper.deuteron=deuteron;
    helper.final=final;
    helper.sqrts=sqrts;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.y=y; helper.scale=scale;

    double result=0; double abserr=0;
    const int MULTIPLICITYXINTPOINTS=4;

    
    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_hadronprod;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MULTIPLICITYXINTPOINTS);

    int status;
    status=gsl_integration_qag(&fun, xf, 1.0,
            0, 0.05, MULTIPLICITYXINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "z integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << " y " << y << " pt " << pt << endl;


    result *= K / SQR(2.0*M_PI);

    return result;
}

/*
 * Multiplicity integrated over rapidity and p_T range
 */
double Inthelperf_hadronprod_yint(double y, void *p);
double Inthelperf_hadronprod_ptint(double pt, void *p);
const int HADRONPROD_YINTPOINTS=3;
const int HADRONPROD_PTINTPOINTS=3;
const double HADRONPROD_INTACCURACY=0.01;

double AmplitudeLib::HadronMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    Inthelper_hadronprod helper;
    helper.miny=miny; helper.maxy=maxy;

    helper.N=this; 
    helper.deuteron=deuteron;
    helper.final=final;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.sqrts=sqrts;
    helper.maxpt=maxpt; helper.minpt=minpt;
    helper.ptint=true;

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_yint;
    fun.params=&helper;
    
    // If miny=maxy, don't integrate over y, only over p_T
    if (std::abs(maxy-miny)<0.001)
		return Inthelperf_hadronprod_yint(maxy, &helper);

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(HADRONPROD_YINTPOINTS);

    int status=0; double abserr, result;
    status=gsl_integration_qag(&fun, miny, maxy,
            0, HADRONPROD_INTACCURACY, HADRONPROD_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "Yint failed at " << LINEINFO <<", result " << result <<", relerr "
            << std::abs(abserr/result) << endl;
    }
    
    gsl_integration_workspace_free(workspace);

    return result;
}

double Inthelperf_hadronprod_yint(double y, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    

    // Don't integrate over pt
    if (par->ptint == false)
        return par->N->dHadronMultiplicity_dyd2pt(y, par->pt, 
            par->sqrts, par->frag, par->pdf, par->deuteron, par->final); 
    
    cout << "# yint y=" << y << endl;
    gsl_function fun;
    fun.function=Inthelperf_hadronprod_ptint;
    par->y = y;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(HADRONPROD_PTINTPOINTS);
    status=gsl_integration_qag(&fun, par->minpt, par->maxpt,
            0, HADRONPROD_INTACCURACY, HADRONPROD_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "yint failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << " y " << y << endl;
    }
    
    gsl_integration_workspace_free(workspace);

    return result*2.0*M_PI; // 2\pi from angular integral

}

double Inthelperf_hadronprod_ptint(double pt, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    return pt*par->N->dHadronMultiplicity_dyd2pt(par->y, pt, par->sqrts, par->frag,
        par->pdf, par->deuteron, par->final);   
}

/*
 * Average hadron yield in rapidity range
 * Integrate HadronMultiplicity_dyd2pt over y region and averge
 */
double AmplitudeLib::AverageHadronMultiplicity(double miny, double maxy, double pt, double sqrts, 
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    Inthelper_hadronprod helper;
    helper.N=this; 
    helper.deuteron=deuteron;
    helper.final=final;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.sqrts=sqrts;
    helper.ptint=false;
    helper.pt=pt;

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_yint;
    fun.params=&helper;
    const int YINT_HADRONAVERAGE = 2;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(YINT_HADRONAVERAGE);
    status=gsl_integration_qag(&fun, miny, maxy,
            0, 0.01, YINT_HADRONAVERAGE,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "yint failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << " pt " << pt << endl;
    }
    
    gsl_integration_workspace_free(workspace);
    
    return result / (maxy-miny);
}

/*
 * Douple parton scattering
 */
struct Inthelper_dps
{
    double sqrts;
    FragmentationFunction* fragfun;
    bool deuteron;
    Hadron final;
    AmplitudeLib* N;
    double y1,y2,pt1,pt2;
    double xf1,xf2;
    double z1;
    PDF *pdf;
    char dps_mode;
    
    double cache_sk1_fund;	// cached N(x_A, pt1/z1)
    double cache_sk1_adj;	// cached N_A(x_A, pt1/z1) (adjoint rep)
};


const int DPS_ZINTPOINTS = 1;

double Inthelperf_dps_z1(double z1, void* p);
double Inthelperf_dps_z2(double z2, void* p);

/*
 * Calculate dN/d^2 p_1 d^2 p_2 dy_1 dy_2 for DPS contribution a+c
 * Note: deuteron is not supported, so one should multiply the
 * final result by 2 or implement isospin-symmetrized deuteron
 */

double AmplitudeLib::DPS(double y1, double y2, double pt1, double pt2, double sqrts,
              FragmentationFunction* fragfun, PDF *pdf, bool deuteron, Hadron final, char dps_mode)
{	const double ncoll = 15.1;	// central d+Au
	const double C_p = 1;
	Inthelper_dps par;
	if (dps_mode != 'c')
	{
		cerr << "dps_mode " << dps_mode << " is not supported " << LINEINFO << endl;
		exit(1);
	}
	par.xf1 = pt1/sqrts*std::exp(y1);
	par.xf2 = pt2/sqrts*std::exp(y2);
	
	if (par.xf1>1 or par.xf2>1)
		return 0;	// kinematically forbidden
	
	par.y1=y1; par.y2=y2; par.pt1=pt1; par.pt2=pt2; par.N=this;
	par.pdf= pdf; par.sqrts=sqrts; par.deuteron=deuteron; par.final=final;
	par.fragfun=fragfun;
	par.dps_mode=dps_mode;
	
	gsl_function fun;
    fun.function=Inthelperf_dps_z1;
    fun.params=&par;
    
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_ZINTPOINTS);
    status=gsl_integration_qag(&fun, par.xf1, 1.0,
            0, 0.1, DPS_ZINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "z1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    
    gsl_integration_workspace_free(workspace);
	
	result *= C_p;
	
	
	result /= std::pow(2.0*M_PI, 4);
	return result;
	
}

double Inthelperf_dps_z1(double z1, void* p)
{
	Inthelper_dps* par = (Inthelper_dps*)p;
	par->z1=z1;
	gsl_function fun;
    fun.function=Inthelperf_dps_z2;
    fun.params=par;
    
	// Calculate scattering amplitudes with momenta p_1/z_1 as from now on 
	// z_1 and x_{A1} are constant
	double xa1 = par->pt1/z1*std::exp(-par->y1)/par->sqrts;
	double ya1 = std::log(par->N->X0() / xa1);
	par->N->InitializeInterpolation(ya1);
	par->cache_sk1_fund = par->N->S_k(par->pt1/z1, ya1, false);	// fundamental
	par->cache_sk1_adj = par->N->S_k(par->pt1/z1, ya1, true); // adjoint
	
	double minz2 = z1*par->xf2 / (z1-par->xf1);
	if (minz2>1) return 0;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_ZINTPOINTS);
    status=gsl_integration_qag(&fun, minz2, 1.0,
            0, 0.1, DPS_ZINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "z2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
    return result/SQR(z1);
    
}

double Inthelperf_dps_z2(double z2, void* p)
{
	Inthelper_dps* par = (Inthelper_dps*)p;
	par->N->SetOutOfRangeErrors(false);
	
	double xp1 = par->xf1 / par->z1;
	double xp2 = par->xf2 / z2;
	double scale = std::max(par->pt1, par->pt2);
	
	if (xp1 + xp2 >=1)
		return 0;	// Kinematical constraint, two quarks from the same proton
		// probably never true, as lower limit of z2 integral is chosen s.t.
		// the kinematical constraint is fulfilled
	

	//double xa1 = xp1 * std::exp(-2.0*par->y1);
	double xa2 = xp2 * std::exp(-2.0*par->y2);
	//double ya1 = std::log( par->N->X0() / xa1);
	double ya2 = std::log( par->N->X0() / xa2);
	par->N->InitializeInterpolation(ya2);
	double nf2 = par->N->S_k(par->pt2/z2, ya2);
	double na2 = par->N->S_k(par->pt2/z2, ya2, true);	// adjoint rep
	
	std::vector<Parton> partons; partons.push_back(U); partons.push_back(D); partons.push_back(G);
	double result=0;
	for (int p1ind=0; p1ind<partons.size(); p1ind++)
	{
		for (int p2ind=0; p2ind<=p1ind; p2ind++)
		{
			// Ddpf() is symmetrized DPDF with kinematical constraint
			// Combinatorics: hadron 1/2 can be produced by both partons 1/2
			double cont_1 = par->pdf->Dpdf(xp1, xp2, scale, partons[p1ind], partons[p2ind])
			 * par->fragfun->Evaluate(partons[p1ind], par->final, par->z1, scale)
			 * par->fragfun->Evaluate(partons[p2ind], par->final, z2, scale);
			
			// Gluon scattering: S() in adjoint rep
			if (partons[p1ind]!=G) cont_1*= par->cache_sk1_fund;
			else cont_1 *= par->cache_sk1_adj;
			if (partons[p2ind]!=G) cont_1*= nf2;
			else cont_1 *= na2;
			 
			double cont_2 = par->pdf->Dpdf(xp2, xp1, scale, partons[p1ind], partons[p2ind])
			 * par->fragfun->Evaluate(partons[p2ind], par->final, par->z1, scale)
			 * par->fragfun->Evaluate(partons[p1ind], par->final, z2, scale);	
			 
			 double sk1=0, sk2=0;
			 // Gluon scattering: S() in adjoint rep
			 if (partons[p1ind]!=G) cont_2 *= nf2;
			 else cont_2 *= na2;
			 if (partons[p2ind]!=G) cont_2 *= par->cache_sk1_fund;
			 else cont_2 *= par->cache_sk1_adj;
			 
			 result += cont_1 + cont_2;
		}
	}		
	
	
	return result/SQR(z2);		
}

/*
 * DPS in parton level
 * Includes sum over quark flavors&gluons, no fragmentation
 * 
 * Code is very similar than the code above, but as we are not here
 * integrating over z1 and z2, we keep this separate as above we
 * use some caching techiques to optimize code
 */
double AmplitudeLib::DPS_partonlevel(double y1, double y2, double pt1, double pt2, double sqrts,
              PDF* pdf, bool deuteron, char dps_mode, double scale)
{
	if (scale<0) scale=std::max(pt1, pt2);
	
	SetOutOfRangeErrors(false);
	
	double xp1 = pt1*std::exp(y1) / sqrts; double xp2 = pt2*std::exp(y2) / sqrts;
	double xa1 = pt1*std::exp(-y1) / sqrts; double xa2 = pt2*std::exp(-y2) / sqrts;
	
	double ya1 = std::log(X0() / xa1); double ya2 = std::log(X0()/xa2);
	InitializeInterpolation(ya1);
	double nf1 = S_k(pt1, ya1); double na1 = S_k(pt1, ya1, true);
	InitializeInterpolation(ya2);
	double nf2= S_k(pt2, ya2); double na2 = S_k(pt2, ya2, true);
	
	cout << "nf1 " << nf1 << " na1 " << na1 << " nf2 " << nf2 << " na2 " << na2 << endl;
	
	
	std::vector<Parton> partons; partons.push_back(U); partons.push_back(D); partons.push_back(G);
	double result=0;

	for (int p1ind=0; p1ind<partons.size(); p1ind++)
	{
		for (int p2ind=0; p2ind<partons.size(); p2ind++)
		{
			// Combinations uu, ud, ug, du, dd, dg, gu, gd, gg
			double xf=0;
			if (dps_mode=='c')
			{
				double tmp = pdf->Dpdf(xp1, xp2, scale, partons[p1ind], partons[p2ind]);
				if (partons[p1ind]==G) tmp *= na1; else tmp*=nf1;
				if (partons[p2ind]==G) tmp*=na2; else tmp*= nf2;
				if (deuteron) tmp*=2;	// neglect isospin
				result += tmp;
			}
			else
			{
				double tmp = 2.0*pdf->xq(xp1, scale, partons[p1ind])*pdf->xq(xp2, scale, partons[p2ind]);	// 2.0: neglect isospin symmetry, hadron 1/2 can come from p/n
				if (partons[p1ind]==G) tmp *= na1; else tmp*=nf1;
				if (partons[p2ind]==G) tmp*=na2; else tmp*= nf2;
				result += tmp;
			}
		}
	}
	
	result /= std::pow(2.0*M_PI, 4.0);
	return result;
}	

/*
 * DPS integrated over pt and y range
 * Integrates over d^2 p_1 d^2 p_2 dy1 dy2
 * Notice that for CP() one must calculate dN/d\phi
 * That is, this must be divided by 2\pi when calculating
 * DPS
 * 
 * Also Deuteron is not supported: if isospin is neglected the final
 * result is obtained by multiplying this by 2
 */
struct Inthelper_dpsint
{
	double miny,maxy,minpt,maxpt,sqrts;
	double pt1,y1,y2;
	AmplitudeLib* N;
	FragmentationFunction* fragfun;
	PDF* pdf;
	bool deuteron;
	Hadron final;
	char dps_mode;
};
const int DPS_YINTPOINTS = 1;
const int DPS_PTINTPOINTS = 1;
double Inthelperf_dpsint_y1(double y1, void* p);
double Inthelperf_dpsint_y2(double y2, void* p);
double Inthelperf_dpsint_pt1(double pt1, void* p);
double Inthelperf_dpsint_pt2(double pt2, void* p);
double AmplitudeLib::DPSMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
			FragmentationFunction* fragfun, PDF* pdf, bool deuteron, Hadron final, char dps_mode)
{
	Inthelper_dpsint par;
	par.dps_mode = dps_mode;
	par.N=this;
	par.pdf= pdf; par.sqrts=sqrts; par.deuteron=deuteron; par.final=final;
	par.fragfun=fragfun; par.miny=miny; par.maxy=maxy; par.minpt=minpt;
	par.maxpt=maxpt;
	
	if (deuteron)
		cerr << "Deuteron DPS is not supported, falling back to proton" << endl;
	
	cout <<"# DPS contribution "
	/*integrating yrange " << miny << " - " << maxy 
		<<", ptrange " << minpt << " " << maxpt */
		<< " contrib: " << par.dps_mode << endl;
	
	if (dps_mode != 'b' and dps_mode != 'c')
	{
		cerr << "DPSmode " << dps_mode << " is invalid!" << endl;
		return 0;
	}

	gsl_function fun;
    fun.function=Inthelperf_dpsint_y1;
    fun.params=&par;

	if (std::abs(miny-maxy)<0.1) 
	{
		cout <<"# Not integrating over rapidity range... setting y=" << miny << endl;
		par.y1=miny; par.y2=miny;
		return SQR(2.0*M_PI)*Inthelperf_dpsint_y2(miny, &par); //(2\pi)^2 from angural integrals
	}
	std::vector<double> yvals; 
    yvals.push_back(2.4);  yvals.push_back(3.2); yvals.push_back(4);
    //yvals.push_back(3); yvals.push_back(3.4); yvals.push_back(3.8);
    //yvals.push_back(3); yvals.push_back(3.266); yvals.push_back(3.533); yvals.push_back(3.8);
    double result=0;
    for (int y1ind=0; y1ind<yvals.size(); y1ind++)
    {
		double tmpres=0;
		for (int y2ind=0; y2ind<yvals.size(); y2ind++)
		{
			cout << "# y1 " << yvals[y1ind] << " y2 " << yvals[y2ind] << endl;
			par.y1=yvals[y1ind]; par.y2=yvals[y2ind];
			double res = Inthelperf_dpsint_y2(par.y2, &par);
			/*if (y2ind==0 or y2ind==4) tmpres += res;
			else if (y2ind==1 or y2ind==3) tmpres += 4.0*res;
			else if (y2ind==2) tmpres += 2.0*res;*/
			if (yvals.size()==3)
			{
				if (y2ind==1) tmpres += 4.0*res;
				else tmpres += res;
			} 
			else if (yvals.size()==4)
			{
				if (y2ind==0 or y2ind==3) tmpres += res;
				else tmpres += 3.0*res;
			}
						
		}
		if (yvals.size()==3)
		{
			tmpres *=  (yvals[yvals.size()-1] - yvals[0])/6.0;
			if (y1ind==1) result += 4.0*tmpres;
			else result += tmpres;
		}
		else if (yvals.size()==4)
		{
			tmpres *= ( (yvals[yvals.size()-1] - yvals[0])/8.0);
			if (y1ind==0 or y1ind==3) result += tmpres;
			else result += 3.0*tmpres;
		}
		else
			cerr << "WTF! " << LINEINFO << endl;
	}
	if (yvals.size()==3)
		result *= (yvals[yvals.size()-1] - yvals[0])/6.0;
	else if (yvals.size()==4)
		result *= ( (yvals[yvals.size()-1] - yvals[0])/8.0);
    
    /*int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_YINTPOINTS);
    status=gsl_integration_qag(&fun, miny, maxy,
            0, 0.1, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "y1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	*/
	return result*SQR(2.0*M_PI);	// 4\pi^2 from angular integrals
	
}
double Inthelperf_dpsint_y1(double y1, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	par->y1=y1;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_y2;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_YINTPOINTS);
    status=gsl_integration_qag(&fun, par->miny, par->maxy,
            0, 0.1, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "y2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}
double Inthelperf_dpsint_y2(double y2, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	//cout <<"# y1 " << par->y1 << " y2 "  << y2 << endl;
	par->y2=y2;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_pt1;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);
    status=gsl_integration_qag(&fun, 2, 4,
            0, 0.1, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "pt1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}

double Inthelperf_dpsint_pt1(double pt1, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	//cout <<"# pt1 " << pt1 << endl;
	par->pt1=pt1;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_pt2;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);
    status=gsl_integration_qag(&fun,1, pt1,
            0, 0.1, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "pt2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}

double Inthelperf_dpsint_pt2(double pt2, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	//cout <<"# pt1 " << par->pt1 << " pt2 " << pt2 << endl;
	
	
	double result = 0;
	if (par->dps_mode=='a' or par->dps_mode=='c')
		result = par->N->DPS(par->y1, par->y2, par->pt1, pt2, 
			par->sqrts, par->fragfun, par->pdf, par->deuteron, par->final, par->dps_mode);
	// b: isospin symmetry is neglected, here we take one quark/gluon from p and an another one from d
	// but there is also similar contribution where p and n are exchanged, thus this must be multiplied
	// by 2 when calculating pedestal (this code does not multiply the result by 2)
	else 
	{
		
		double scale = std::max(pt2, par->pt1);
		result = par->N->dHadronMultiplicity_dyd2pt(par->y1, par->pt1, par->sqrts, par->fragfun, par->pdf, false, par->final, scale)
			*  par->N->dHadronMultiplicity_dyd2pt(par->y2, pt2, par->sqrts, par->fragfun, par->pdf, false, par->final, scale);
	}
	return pt2 * par->pt1 * result;
}
