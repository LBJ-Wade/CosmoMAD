///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CosmoMad.                                    //
//                                                                   //
// CosmoMad is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published //
// by the Free Software Foundation, either version 3 of the License, //
// or (at your option) any later version.                            //
//                                                                   //
// CosmoMad is distributed in the hope that it will be useful, but   //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CosmoMad.  If not, see <http://www.gnu.org/licenses/>. //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cosmo_mad.h"

//Number of points for the a(t) relation
#define CSM_NINTERP_A 1000
//For x>CSM_SPH_BESSEL_XMAX, j_l(x) ~ cos(x-(l+1)*pi/2)/x
#define CSM_SPH_BESSEL_XMAX 100.
//For x>CSM_CYL_BESSEL_XMAX, J_0(x) ~ sqrt(2/(pi*x))*cos(x-pi/4)
#define CSM_CYL_BESSEL_XMAX 10.

//Scales for splines in Mpc/h
//Logarithmic splines are used for CSM_RMIN_LOG < r < CSM_RMIN_LIN
//Linear splines are used for CSM_RMIN_LIN < r < CSM_RMAX_LIN
#define CSM_RMIN_LOG 0.1
#define CSM_RMIN_LIN 15.
#define CSM_RMAX_LIN 500.

//Logarithmic splines exist for
//  CSM_RMIN_LOG_SPLINE < r < CSM_RMAX_LOG_SPLINE
//Linear splines exist for
//  CSM_RMIN_LIN_SPLINE < r < CSM_RMAX_LIN_SPLINE
#define CSM_RMIN_LOG_SPLINE 0.05
#define CSM_RMAX_LOG_SPLINE 20.
#define CSM_RMIN_LIN_SPLINE 10.
#define CSM_RMAX_LIN_SPLINE 550.

//Interval in log(r) for logarithmic splines
#define CSM_DLOG 0.01
//Interval in r for linear splines (in Mpc/h)
#define CSM_DR 1.

/**** Error handler ****/
static gsl_error_handler_t *csm_gsl_error_handler_old;


/****************************/
/*     General routines     */
/****************************/
static void error_mem_out(void)
{
  //////
  // Memory shortage handler
  fprintf(stderr,"CosmoMad: Out of memory \n");
  exit(1);
}

static int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}
  
static void int_error_handle(int status,double result,
                             double error)
{
  //////
  // Error handler for gsl
  if(isnan(result)) {
    fprintf(stderr,"CosmoMad: NAN found \n");
  }
  else{
    if(status==GSL_EROUND) {
      fprintf(stderr,"CosmoMad: Roundoff error: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EMAXITER) {
      fprintf(stderr,"CosmoMad: Ran out of iterations: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ESING) {
      fprintf(stderr,"CosmoMad: Singularity found: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EDIVERGE) {
      fprintf(stderr,"CosmoMad: Integral seems to diverge: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ETOL) {
      fprintf(stderr,"CosmoMad: Can't reach tolerance: %lE %lE\n",
	      result,error);
    }
    else if(status==GSL_EUNDRFLW)
      fprintf(stderr,"CosmoMad: Underflow: %lE %lE \n",result,error);
    else if(status==GSL_EDOM) {
      fprintf(stderr,"CosmoMad: Outside interpolation range!! %lE %lE\n",
	      result,error);
      exit(1);
    }
    else if(status) {
      fprintf(stderr,"CosmoMad: Unknown error code %d %lf %lf \n",
	      status,result,error);
      exit(1);
    }
  }
}

void csm_unset_gsl_eh(void)
{
  //////
  // Disables GSL default error handler
  csm_gsl_error_handler_old=gsl_set_error_handler_off();
}

void csm_set_verbosity(Csm_params *par,int verb)
{
  //////
  // Sets verbosity level
  par->flag_verbose=verb;
}

static void csm_bg_params_free(Csm_bg_params *par)
{
  //////
  // bg_params destructor
  free(par->eh);
  if(par->at_spline_set) {
    gsl_interp_accel_free(par->intacc_at);
    gsl_spline_free(par->spline_at);
  }
  par->at_spline_set=0;
  free(par);
}

static Csm_bg_params *csm_bg_params_new(void)
{
  //////
  // bg_params creator
  // Default Planck parameters
  Csm_bg_params *bg=(Csm_bg_params *)malloc(sizeof(Csm_bg_params));
  if(bg==NULL) error_mem_out();
  bg->OM=0.315;
  bg->OL=0.685;
  bg->OB=0.049;
  bg->h=0.673;
  bg->w0=-1;
  bg->wa=0;
  bg->TCMB=2.725;
  bg->OK=0;
  bg->ksign=0;
  bg->normalDE=1;
  bg->constantw=1;

  bg->eh=(Csm_EH_params *)malloc(sizeof(Csm_EH_params));
  bg->at_spline_set=0;
  
  return bg;
}

static void csm_pk_params_free(Csm_pk_params *par)
{
  //////
  // Pk params destructor
  if(par->karr_set)
    free(par->logkarr);
  par->karr_set=0;

  if(par->pkL_spline_set) {
    gsl_interp_accel_free(par->intacc_pkL);
    gsl_spline_free(par->spline_pkL);
  }
  par->pkL_spline_set=0;

  if(par->pkNL_spline_set) {
    gsl_interp_accel_free(par->intacc_pkNL);
    gsl_spline_free(par->spline_pkNL);
  }
  par->pkNL_spline_set=0;

  if(par->pkmulti_spline_set) {
    int ii;
    for(ii=0;ii<par->l_multi_max/2+1;ii++) {
      gsl_interp_accel_free((par->intacc_pkmulti)[ii]);
      gsl_spline_free((par->spline_pkmulti)[ii]);
    }
    free(par->intacc_pkmulti);
    free(par->spline_pkmulti);
  }
  par->pkmulti_spline_set=0;

  free(par);
} 

static Csm_pk_params *csm_pk_params_new(void)
{
  //////
  // pk_params creator
  Csm_pk_params *pk=(Csm_pk_params *)malloc(sizeof(Csm_pk_params));
  if(pk==NULL) error_mem_out();
  
  pk->s8=0.8;
  pk->ns=0.96;
  pk->fg=0;
  pk->bias=1;
  pk->sigv=0;
  pk->gf2=1;
  pk->norm=1;
  pk->l_multi_max=0;
  pk->use_RPT=0;
  pk->use_RPT_ss=0;
  pk->pkNL_set=0;

  pk->numk=0;
  pk->logkcamb=3;
  pk->logkmin=-3;
  pk->logkcambNL=3;
  pk->logkminNL=-3;
  pk->karr_set=0;
  pk->logkarr=NULL;
  pk->pkL_spline_set=0;
  pk->intacc_pkL=NULL;
  pk->spline_pkL=NULL;
  pk->pkNL_spline_set=0;
  pk->intacc_pkNL=NULL;
  pk->spline_pkNL=NULL;
  pk->pkmulti_spline_set=0;
  pk->intacc_pkmulti=NULL;
  pk->spline_pkmulti=NULL;

  return pk;
}

static void csm_xi_params_free(Csm_xi_params *par,int nl)
{
  //////
  // xi_params creator
  int ii;
  for(ii=0;ii<nl;ii++) {
    gsl_interp_accel_free(par->intacc_ximulti_log[ii]);
    gsl_spline_free(par->spline_ximulti_log[ii]);
    gsl_interp_accel_free(par->intacc_ximulti_lin[ii]);
    gsl_spline_free(par->spline_ximulti_lin[ii]);
  }
  free(par->intacc_ximulti_log);
  free(par->intacc_ximulti_lin);
  free(par->spline_ximulti_log);
  free(par->spline_ximulti_lin);
  free(par->ximultimin);
  free(par);
}

static Csm_xi_params *csm_xi_params_new(void)
{
  //////
  // xi_params creator
  Csm_xi_params *xi=(Csm_xi_params *)malloc(sizeof(Csm_xi_params));
  if(xi==NULL) error_mem_out();

  xi->intacc_ximulti_log=NULL;
  xi->intacc_ximulti_lin=NULL;
  xi->spline_ximulti_log=NULL;
  xi->spline_ximulti_lin=NULL;

  return xi;
}

static void csm_mf_params_free(Csm_mf_params *par)
{
  if(par->sM_spline_set) {
    gsl_interp_accel_free(par->intacc_sM);
    gsl_spline_free(par->spline_sM);
  }
  if(par->dsM_spline_set) {
    gsl_interp_accel_free(par->intacc_dsM);
    gsl_spline_free(par->spline_dsM);
  }
  free(par);
}

static Csm_mf_params *csm_mf_params_new(void)
{
  //////
  // mf_params creator
  Csm_mf_params *mf=(Csm_mf_params *)malloc(sizeof(Csm_mf_params));
  if(mf==NULL) error_mem_out();
  mf->sM_spline_set=0;
  mf->dsM_spline_set=0;

  return mf;
}

void csm_params_free(Csm_params *par)
{
  //////
  // csm_params destructor
  if(par->bg_params_set)
    csm_bg_params_free(par->bg);
  par->bg_params_set=0;
  if(par->pk_params_set)
    csm_pk_params_free(par->pk);
  par->pk_params_set=0;
  if(par->xi_spline_set)
    csm_xi_params_free(par->xi,par->pk->l_multi_max/2+1);
  par->xi_spline_set=0;
  if(par->mf_params_set)
    csm_mf_params_free(par->mf);
  par->mf_params_set=0;
  free(par);
}

Csm_params *csm_params_new(void)
{
  //////
  // csm_params destructor
  Csm_params *par=(Csm_params *)malloc(sizeof(Csm_params));
  if(par==NULL) error_mem_out();
  par->flag_verbose=0;
  par->bg_params_set=0;
  par->bg=NULL;
  par->pk_params_set=0;
  par->pk=NULL;
  par->xi_spline_set=0;
  par->xi=NULL;
  par->mf_params_set=0;
  par->mf=NULL;

  return par;
}

/****************************/
/*   Background cosmology   */
/****************************/
static void set_EH_params(Csm_bg_params *par)
{
  //////
  // Computes Eisenstein & Hu parameters for
  // P_k and r_sound
  double OMh2,OBh2;
  double th2p7;
  OMh2=par->OM*par->h*par->h;
  OBh2=par->OB*par->h*par->h;
  th2p7=par->TCMB/2.7;
  par->eh->th2p7=th2p7;
  par->eh->zeq=25000*OMh2/pow(th2p7,4);
  par->eh->keq=0.0746*OMh2/(par->h*th2p7*th2p7);

  double b1,b2;
  b1=0.313*pow(OMh2,-0.419)*(1+0.607*pow(OMh2,0.674));
  b2=0.238*pow(OMh2,0.223);
  par->eh->zdrag=1291*pow(OMh2,0.251)*(1+b1*pow(OBh2,b2))/
    (1+0.659*pow(OMh2,0.828));

  double Req,Rd;
  //EH actually says 1100 instead of 1000...
  Req=31.5*OBh2*1100/(par->eh->zeq*pow(th2p7,4));
  Rd=31.5*OBh2*1100/(par->eh->zdrag*pow(th2p7,4));
  par->eh->rsound=2/(3*par->eh->keq)*sqrt(6/Req)*
    log((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)));

  par->eh->kSilk=1.6*pow(OBh2,0.52)*pow(OMh2,0.73)*
    (1+pow(10.4*OMh2,-0.95))/par->h;

  double a1,a2,b_frac;
  a1=pow(46.9*OMh2,0.670)*(1+pow(32.1*OMh2,-0.532));
  a2=pow(12.0*OMh2,0.424)*(1+pow(45.0*OMh2,-0.582));
  b_frac=OBh2/OMh2;
  par->eh->alphac=pow(a1,-b_frac)*pow(a2,-b_frac*b_frac*b_frac);

  double bb1,bb2;
  bb1=0.944/(1+pow(458*OMh2,-0.708));
  bb2=pow(0.395*OMh2,-0.0266);
  par->eh->betac=1/(1+bb1*(pow(1-b_frac,bb2)-1));

  double y=(1+par->eh->zeq)/(1+par->eh->zdrag);
  double gy=y*(-6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
  par->eh->alphab=2.07*par->eh->keq*par->eh->rsound*pow(1+Rd,-0.75)*gy;

  par->eh->betab=0.5+b_frac+(3-2*b_frac)*sqrt(pow(17.2*OMh2,2)+1);

  par->eh->bnode=8.41*pow(OMh2,0.435);

  par->eh->rsound_approx=par->h*44.5*log(9.83/OMh2)/
    sqrt(1+10*pow(OBh2,0.75));
}

static double aeqmin(Csm_bg_params *par)
{
  //////
  // Returns MIN(aeq_k,aeq_L), where aeq_i is the scale
  // factor at equality of M with k or L.
  double aeqK=1;
  double aeqL=1;
  
  if(par->ksign!=0)
    aeqK=par->OM/fabs(par->OK);

  if (par->OL!=0) {
    if(par->normalDE)
      aeqL=pow(par->OM/par->OL,0.333);
    else
      aeqL=pow(par->OM/par->OL,-1/(3*par->w0));
  }

  return CSM_MIN(aeqK,aeqL);
}

static double sinn(double x,int sign)
{
  //////
  //         { sin(x)  , if k==1
  // sinn(x)={  x      , if k==0
  //         { sinh(x) , if k==-1
  double dum;

  if(sign==-1)
    dum=sinh(x);
  else if(sign==1)
    dum=sin(x);
  else
    dum=x;
  
  return dum;
}

static double naHm1(double a,void *params)
{
  //////
  // H0/(a*H[a])
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) 
    dum=sqrt(a/(par->OM+par->OL*a*a*a+par->OK*a));
  else if(par->constantw) {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum;
}

static double na2Hm1(double a,void *params)
{
  //////
  // H0/(a^2*H[a])
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) 
    dum=1/sqrt(a*(par->OM+par->OL*a*a*a+par->OK*a));
  else if(par->constantw) {
    dum=1/sqrt(a*(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=1/sqrt(a*(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		  exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum;
}

static double naHm3(double a,void *params)
{
  //////
  // (H0/(a*H[a]))^3
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) {
    dum=sqrt(a/(par->OM+par->OL*a*a*a+par->OK*a));
  }
  else if(par->constantw) {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		  exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum*dum*dum;
}

static void costime(Csm_bg_params *par,
		    double aa,double *t,double *delt)
{
  //////
  // Cosmic time. The returned value is in Gyr/h.
  double alim;
  
  alim=0.01*par->a_equality;
  if(aa<=alim) {
    *t=0.666666666666666666*sqrt(aa*aa*aa/par->OM)*CSM_HGYR;
    *delt=0;
  }
  else {
    double relerrt=1E-6;
    double integral,errintegral,int0;
    size_t sdum;
    gsl_function integrand;
    int stat;

    int0=0.666666666666666666*sqrt(alim*alim*alim/par->OM);
    integrand.function=&naHm1;
    integrand.params=par;

    stat=gsl_integration_qng(&integrand,alim,aa,0,relerrt,
			     &integral,&errintegral,&sdum);
    int_error_handle(stat,integral,errintegral);
    *t=(int0+integral)*CSM_HGYR;
    //    *delt=errintegral*CSM_HGYR;
  }
}

static void parthor(Csm_bg_params *par,
		    double aa,double *ph,double *delph)
{
  //////
  // Particle horizon. The returned value is in Mpc/h.
  double alim;

  alim=0.01*par->a_equality;
  if(aa<=alim) {
    *ph=2*sqrt(aa/par->OM)*CSM_HMPC;
    *delph=0;
  }
  else {
    double relerrt=1E-6;
    double integral,errintegral,int0;
    size_t sdum;
    gsl_function integrand;
    int stat;

    int0=2*sqrt(alim/par->OM);
    integrand.function=&na2Hm1;
    integrand.params=par;

    stat=gsl_integration_qng(&integrand,alim,aa,0,relerrt,
			     &integral,&errintegral,&sdum);
    int_error_handle(stat,integral,errintegral);
    *ph=(int0+integral)*CSM_HMPC;
    //    *delph=errintegral*CSM_HMPC;
  }
}

static double h_normalized(Csm_params *par,double aa)
{
  //////
  // Normalized hubble rate
  if(par->bg->normalDE) {
    return sqrt((par->bg->OM+par->bg->OL*aa*aa*aa+par->bg->OK*aa)/
		(aa*aa*aa));
  }
  else if(par->bg->constantw) {
    return sqrt((par->bg->OM+par->bg->OL*pow(aa,-3*par->bg->w0)+
		 par->bg->OK*aa)/(aa*aa*aa));
  }
  else {
    return sqrt((par->bg->OM+par->bg->OL*pow(aa,-3*(par->bg->w0+par->bg->wa))*
		 exp(3*par->bg->wa*(aa-1))+par->bg->OK*aa)/(aa*aa*aa));
  }
}

double csm_omega_matter(Csm_params *par,double aa)
{
  //////
  // Normalized hubble rate
  if(par->bg->normalDE) {
    return par->bg->OM/(par->bg->OM+par->bg->OL*aa*aa*aa+par->bg->OK*aa);
  }
  else if(par->bg->constantw) {
    return par->bg->OM/(par->bg->OM+par->bg->OL*pow(aa,-3*par->bg->w0)+
			par->bg->OK*aa);
  }
  else {
    return par->bg->OM/(par->bg->OM+par->bg->OL*pow(aa,-3*(par->bg->w0+par->bg->wa))*
			exp(3*par->bg->wa*(aa-1))+par->bg->OK*aa);
  }
}

double csm_hubble(Csm_params *par,double aa)
{
  //////
  // Hubble rate at aa in h/Mpc
  return h_normalized(par,aa)/CSM_HMPC;
}

double csm_cosmic_time(Csm_params *par,double aa)
{
  //////
  // Cosmic time in Gyr/h.
  double ct, ect;

  costime(par->bg,aa,&ct,&ect);
  return ct;
}

double csm_particle_horizon(Csm_params *par,double aa)
{
  //////
  // Particle horizon
  double ph, eph;

  parthor(par->bg,aa,&ph,&eph);
  return ph;
}

double csm_radial_comoving_distance(Csm_params *par,double aa)
{
  //////
  // chi(a)
  double rcd;

  rcd=csm_particle_horizon(par,aa);
  return par->bg->phorizon-rcd;
}

double csm_curvature_comoving_distance(Csm_params *par,double aa)
{
  //////
  // r(a)
  if(par->bg->ksign==0)
    return csm_radial_comoving_distance(par,aa);
  else {
    double dum;
    double ksq=sqrt(fabs(par->bg->OK));
    dum=csm_radial_comoving_distance(par,aa)/CSM_HMPC;
    dum=sinn(ksq*dum,par->bg->ksign)/ksq;
    return dum*CSM_HMPC;
  }
}

double csm_angular_diameter_distance(Csm_params *par,double aa)
{
  //////
  // d_A(a)
  double dum;

  dum=aa*csm_curvature_comoving_distance(par,aa);
  return dum;
}

double csm_luminosity_distance(Csm_params *par,double aa)
{
  //////
  // d_L(a)
  double dum;

  dum=csm_curvature_comoving_distance(par,aa)/aa;
  return dum;
}

static int growth_ode_system(double a,const double y[],double dydt[],void *params)
{
  Csm_params *par=(Csm_params *)params;
  double hnorm=h_normalized(par,a);
  
  dydt[0]=y[1]/(a*a*a*hnorm);
  dydt[1]=1.5*hnorm*a*csm_omega_matter(par,a)*y[0];

  return GSL_SUCCESS;
}

#define CSM_GROWTH_AINIT 1E-6
void csm_growth_factor_and_growth_rate(Csm_params *par,double aa,double *gf,double *fg)
{
  if(aa<CSM_GROWTH_AINIT) {
    *gf=aa;
    *fg=1;
  }
  else {
    double y[2];
    double ainit=CSM_GROWTH_AINIT;
    gsl_odeiv2_system sys={growth_ode_system,NULL,2,par};
    gsl_odeiv2_driver *d=
      gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkck,0.1*CSM_GROWTH_AINIT,0,1E-10);

    y[0]=CSM_GROWTH_AINIT;
    y[1]=CSM_GROWTH_AINIT*CSM_GROWTH_AINIT*CSM_GROWTH_AINIT*h_normalized(par,CSM_GROWTH_AINIT);

    int status=gsl_odeiv2_driver_apply(d,&ainit,aa,y);
    if(status!=GSL_SUCCESS) {
      fprintf(stderr,"CosmoMad: ODE didn't converge\n");
      exit(1);
    }
    
    *gf=y[0];
    *fg=y[1]/(aa*aa*h_normalized(par,aa)*y[0]);
  }
}      

double csm_growth_factor(Csm_params *par,double aa)
{
  //////
  // D(a)
  double gf,dum;
  csm_growth_factor_and_growth_rate(par,aa,&gf,&dum);

  return gf;
}

double csm_f_growth(Csm_params *par,double aa)
{
  double fg,dum;
  csm_growth_factor_and_growth_rate(par,aa,&dum,&fg);

  return fg;
}

static void set_at_spline(Csm_bg_params *par)
{
  //////
  // Sets spline for a(t)
  double aarr[CSM_NINTERP_A],tarr[CSM_NINTERP_A];
  int ii;

  par->intacc_at=gsl_interp_accel_alloc();
  par->spline_at=gsl_spline_alloc(gsl_interp_cspline,CSM_NINTERP_A);
  for(ii=0;ii<CSM_NINTERP_A;ii++) {
    double ct,ect;
    double aa=ii/(CSM_NINTERP_A-1.);
    costime(par,aa,&ct,&ect);
    aarr[ii]=aa;
    tarr[ii]=ct;
  }
  gsl_spline_init(par->spline_at,tarr,aarr,CSM_NINTERP_A);
  par->at_spline_set=1;
}

double csm_scale_factor(Csm_params *par,double t)
{
  //////
  // a(t) using spline
  double result;
  int stat;

  if(!(par->bg->at_spline_set))
    set_at_spline(par->bg);

  stat=gsl_spline_eval_e(par->bg->spline_at,t,
			 par->bg->intacc_at,&result);
  int_error_handle(stat,t,result);
  
  return result;
}

void csm_background_set(Csm_params *par,
			double OmegaM,double OmegaL,double OmegaB,
			double ww,double wwa,double hh,double T_CMB)
{
  //////
  // This initializes the cosmological parameters.
  if(par->bg_params_set)
    csm_bg_params_free(par->bg);
  par->bg_params_set=0;

  par->bg=csm_bg_params_new();
  par->bg_params_set=1;
  par->bg->h=hh;
  par->bg->w0=ww;
  par->bg->wa=wwa;
  par->bg->OM=OmegaM;
  par->bg->OL=OmegaL;
  par->bg->OK=1-par->bg->OM-par->bg->OL;
  par->bg->OB=OmegaB;
  par->bg->TCMB=T_CMB;

  //Check parameters
  if(fabs(par->bg->wa)<1E-6) {
    if(fabs(par->bg->w0+1)<1E-6) {
      par->bg->constantw=1;
      par->bg->normalDE=1;
    }
    else {
      par->bg->constantw=1;
      par->bg->normalDE=0;
    }
  }
  else {
    par->bg->constantw=0;
    par->bg->normalDE=0;
  }

  if(fabs(par->bg->OK)<1E-6)
    par->bg->ksign=0;
  else if(par->bg->OK>0)
    par->bg->ksign=-1;
  else
    par->bg->ksign=1; 

  if(par->bg->OM<=0) {
    fprintf(stderr,"CosmoMad: Wrong matter parameter %.3lf \n",
	    par->bg->OM);
    exit(1);
  }
  if(par->bg->OM<par->bg->OB) {
    fprintf(stderr,"CosmoMad: Wrong M/B parameter %.3lf > %.3lf \n",
	    par->bg->OB,par->bg->OM);
    exit(1);
  }
  if(par->bg->w0>-0.333333) {
    fprintf(stderr,"CosmoMad: DE is too exotic (w=%.3lf \n",par->bg->w0);
    exit(1);
  }
  if(par->bg->TCMB<0) {
    fprintf(stderr,"CosmoMad: Wrong CMB temperature %.3lf \n",
	    par->bg->TCMB);
    exit(1);
  }

  set_EH_params(par->bg);

  if(par->flag_verbose) {
    printf("The cosmological model is:\n");
    printf(" O_M=%.3f O_L=%.3f O_K=%.3f\n",
	   par->bg->OM,par->bg->OL,par->bg->OK);
    printf(" O_B=%.3f w=%.3f h=%.3f\n",
	   par->bg->OB,par->bg->w0,par->bg->h);
    if(par->bg->ksign==0)
      printf(" Flat universe, ");
    else if(par->bg->ksign==1)
      printf(" Closed universe, ");
    else if(par->bg->ksign==-1)
      printf(" Open universe, ");
    if(par->bg->normalDE)
      printf("standard cosmological constant\n");
    else {
      printf("non-standard dark energy");
      if(par->bg->constantw)
	printf("\n");
      else
	printf(": w(a) = %.3lf + %.3lf*(1-a) \n",par->bg->w0,par->bg->wa);
    }
  }
  par->bg->a_equality=aeqmin(par->bg);
  par->bg->phorizon=csm_particle_horizon(par,1);
  par->bg->bang_time=csm_cosmic_time(par,1);
  par->bg->growth0=csm_growth_factor(par,1);
  if(par->flag_verbose) {
    printf("\n Time of equality: a_eq=%.5lf\n",par->bg->a_equality);
    printf(" Particle horizon: ");
    printf("chi_H(0)=%.3lE Mpc/h\n",par->bg->phorizon);
    printf(" Sound horizon at drag epoch: ");
    printf("r_s(z_d)=%.3lE Mpc/h\n",par->bg->eh->rsound);
    printf(" Age of the Universe: ");
    printf("t_BB=%.3lE Gyr/h\n",par->bg->bang_time);
    printf(" Present growth factor: ");
    printf("D_0=%.3lf\n\n",par->bg->growth0);
  }
}

double csm_theta_BAO(Csm_params *par,double aa)
{
  //////
  // Position of the angular BAO peak in deg.
  return CSM_RTOD*par->bg->eh->rsound/
    csm_curvature_comoving_distance(par,aa);
}

double csm_Dz_BAO(Csm_params *par,double aa)
{
  //////
  // Position of the radial BAO peak in Dz
  double rs,dz;

  rs=par->bg->eh->rsound;
  dz=csm_hubble(par,aa)*rs;

  return dz;
}


/****************************/
/*      Power Spectrum      */
/****************************/
static double EH_tk0(double keq,double k,double a,double b)
{
  //////
  // Eisentstein & Hu's Tk_0
  double q=k/(13.41*keq);
  double c=14.2/a+386./(1+69.9*pow(q,1.08));
  double l=log(M_E+1.8*b*q);
  return l/(l+c*q*q);
}

static double EH_tkc(Csm_EH_params *eh,double k)
{
  //////
  // Eisenstein & Hu's Tk_c
  double f=1/(1+pow(k*eh->rsound*0.185185185185,4));
  return f*EH_tk0(eh->keq,k,1,eh->betac)+
    (1-f)*EH_tk0(eh->keq,k,eh->alphac,eh->betac);
}


static double EH_tkb(Csm_EH_params *eh,double k)
{
  //////
  // Eisenstein & Hu's Tk_b
  double x=k*eh->rsound;
  double x_bessel;
  double part1;
  double part2;

  if(k==0) x_bessel=0;
  else {
    x_bessel=x*pow(1+eh->bnode*eh->bnode*eh->bnode
		   /(x*x*x),-0.3333333333333333333);
  }

  part1=EH_tk0(eh->keq,k,1,1)/(1+0.03698224852*x*x);
  
  if(k==0) part2=0;
  else part2=eh->alphab/(1+pow(eh->betab/x,3))*
	 exp(-pow(k/eh->kSilk,1.4));
  
  return csm_j_bessel(0,x_bessel)*(part1+part2);
}

static double EH_Pk(Csm_bg_params *bg,double k,
		    double nsc,int with_wiggles)
{
  //////
  // Eisenstein & Hu's power spectrum
  double Tk;
  double b_frac=bg->OB/bg->OM;
  
  if(with_wiggles)
    Tk=b_frac*EH_tkb(bg->eh,k)+(1-b_frac)*EH_tkc(bg->eh,k);
  else {
    double OMh2=bg->OM*bg->h*bg->h;
    double alpha_gamma=1-0.328*log(431*OMh2)*b_frac+0.38*log(22.3*OMh2)*
      b_frac*b_frac;
    double gamma_eff=bg->OM*bg->h*
      (alpha_gamma+(1-alpha_gamma)/
       (1+pow(0.43*k*bg->eh->rsound_approx,4)));
    double q=k*bg->eh->th2p7*bg->eh->th2p7/gamma_eff;
    double l0=log(2*M_E+1.8*q);
    double c0=14.2+731/(1+62.5*q);
    Tk=l0/(l0+c0*q*q);
  }
  
  return pow(k,nsc)*Tk*Tk;
}

static double BBKSPk(Csm_bg_params *bg,double k,double nsc)
{
  //////
  // BBKS's power spectrum.
  double q, Gamma, Tk;
  
  Gamma=bg->OM*bg->h*exp(-bg->OB*(1+sqrt(2*bg->h)/bg->OM));
  q=k/Gamma;
  Tk=log(1.+2.34*q)/
     (2.34*q*pow((1.+3.89*q+259.21*q*q+
                  162.771336*q*q*q+
	          2027.16958*q*q*q*q),0.25));
  return Tk*Tk*pow(k,nsc);
}

void csm_set_linear_pk(Csm_params *par,char *fname,
		       double lkmn,double lkmx,double dlk,
		       double nns,double s8)
{
  //////
  // Reads (or generates) linear matter power spectrum
  int ii;
  double kk,ppk;
  double *pkarr;

  if(!(par->bg_params_set)) {
    fprintf(stderr,"You must set background before you set pk\n");
    exit(1);
  }
  if(par->pk_params_set) { 
    csm_pk_params_free(par->pk);
    par->pk_params_set=0;
  }
  par->pk=csm_pk_params_new();
  par->pk_params_set=1;
  par->pk->ns=nns;

  if(!strcmp(fname,"BBKS")) {
    if(par->flag_verbose)
      printf("P_k set using the BBKS transfer function\n");
    par->pk->logkcamb=lkmx;
    par->pk->logkmin=lkmn;
    par->pk->numk=(int)((lkmx-lkmn)/dlk)+1;
    par->pk->logkarr=malloc(par->pk->numk*sizeof(double));
    if(par->pk->logkarr==NULL) error_mem_out();
    pkarr=malloc(par->pk->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    
    for(ii=0;ii<par->pk->numk;ii++) { 
      par->pk->logkarr[ii]=lkmn+ii*dlk; //log(k) in h Mpc^-1
      kk=pow(10.,par->pk->logkarr[ii]);
      pkarr[ii]=BBKSPk(par->bg,kk,par->pk->ns);
    }
    par->pk->karr_set=1;
  }
  else if(!strcmp(fname,"EH")) {
    if(par->flag_verbose)
      printf("P_k set using the Eisenstein & Hu transfer function\n");
    par->pk->logkcamb=lkmx;
    par->pk->logkmin=lkmn;
    par->pk->numk=(int)((lkmx-lkmn)/dlk)+1;
    par->pk->logkarr=malloc(par->pk->numk*sizeof(double));
    if(par->pk->logkarr==NULL) error_mem_out();
    pkarr=malloc(par->pk->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    
    for(ii=0;ii<par->pk->numk;ii++) { 
      par->pk->logkarr[ii]=lkmn+ii*dlk; //log(k) in h Mpc^-1
      kk=pow(10.,par->pk->logkarr[ii]);
      pkarr[ii]=EH_Pk(par->bg,kk,par->pk->ns,1);
    }
    par->pk->karr_set=1;
  }
  else if(!strcmp(fname,"EH_smooth")) {
    if(par->flag_verbose) {
      printf("P_k set using the smooth Eisenstein & Hu");
      printf(" transfer function\n");
    }
    par->pk->logkcamb=lkmx;
    par->pk->logkmin=lkmn;
    par->pk->numk=(int)((lkmx-lkmn)/dlk)+1;
    par->pk->logkarr=malloc(par->pk->numk*sizeof(double));
    if(par->pk->logkarr==NULL) error_mem_out();
    pkarr=malloc(par->pk->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    
    for(ii=0;ii<par->pk->numk;ii++) { 
      par->pk->logkarr[ii]=lkmn+ii*dlk; //log(k) in h Mpc^-1
      kk=pow(10.,par->pk->logkarr[ii]);
      pkarr[ii]=EH_Pk(par->bg,kk,par->pk->ns,0);
    }
    par->pk->karr_set=1;
  }
  else {
    FILE *fpk;
    if(par->flag_verbose)
      printf("Reading P_k from file: %s\n",fname);
    fpk=fopen(fname,"r");
    if(fpk==NULL) {
      fprintf(stderr,"CosmoMad: Error opening file %s \n",fname);
      exit(1);
    }
    par->pk->numk=linecount(fpk);
    par->pk->logkarr=malloc(par->pk->numk*sizeof(double));
    if(par->pk->logkarr==NULL) error_mem_out();
    pkarr=malloc(par->pk->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    rewind(fpk);
    for(ii=0;ii<par->pk->numk;ii++) {
      int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
      if(stat!=2) {
	fprintf(stderr,"CosmoMad: Error reading file %s, line %d \n",
		fname,ii+1);
	exit(1);
      }
      pkarr[ii]=ppk;
      par->pk->logkarr[ii]=log10(kk); //log(k) in h Mpc^-1
    }
    par->pk->karr_set=1;
    par->pk->logkmin=par->pk->logkarr[0];
    par->pk->logkcamb=par->pk->logkarr[par->pk->numk-1];
    fclose(fpk);
  }
  
  par->pk->pkL_spline_set=1;
  par->pk->intacc_pkL=gsl_interp_accel_alloc();
  par->pk->spline_pkL=gsl_spline_alloc(gsl_interp_cspline,
				       par->pk->numk);
  gsl_spline_init(par->pk->spline_pkL,par->pk->logkarr,
		  pkarr,par->pk->numk);
    
  // normalize
  double sigma_8_ini=sqrt(csm_sig0_L(par,8,8,"TopHat","TopHat"));
  if(par->flag_verbose) {
    printf(" Original sigma8=%lf, ",sigma_8_ini);
  }
  if(s8>0) {
    par->pk->s8=s8;
    if(par->flag_verbose) {
      printf("normalized to sigma8=%lf\n\n",s8);
    }
  }
  else {
    if(par->flag_verbose) {
      printf("not modified\n\n");
    }
    par->pk->s8=sigma_8_ini;
  }

  par->pk->norm=par->pk->s8*par->pk->s8/
    (sigma_8_ini*sigma_8_ini);
  for(ii=0;ii<par->pk->numk;ii++)
    pkarr[ii]*=par->pk->norm;

  gsl_spline_free(par->pk->spline_pkL);
  par->pk->spline_pkL=gsl_spline_alloc(gsl_interp_cspline,
				       par->pk->numk);
  gsl_spline_init(par->pk->spline_pkL,par->pk->logkarr,
		  pkarr,par->pk->numk);
  free(pkarr);
}

static double pk_linear(Csm_pk_params *pkpar,double kk)
{
  //////
  // Returns linear P_k at z=0
  double lgk=log10(kk);
  double pk;

  if(!pkpar->pkL_spline_set) {
    fprintf(stderr,"CosmoMad: Pk is not set \n");
    exit(1);
  }

  if(lgk<pkpar->logkmin) {
    double pk0;
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,pkpar->logkmin,
			       pkpar->intacc_pkL,&pk0);
    int_error_handle(stat,pkpar->logkmin,pk0);
    pk=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkmin));
  }
  else if(lgk<pkpar->logkcamb) {
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,lgk,
			       pkpar->intacc_pkL,&pk);
    int_error_handle(stat,lgk,pk);
  }
  else {
    double pkf;
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,pkpar->logkcamb,
			       pkpar->intacc_pkL,&pkf);
    int_error_handle(stat,pkpar->logkcamb,pkf);
    pk=pkf*pow(10,-3*(lgk-pkpar->logkcamb));
  }

  return pk;
}

double csm_Pk_linear_0(Csm_params *par,double kk)
{
  //////
  // Returns linear P_k at z=0
  return pk_linear(par->pk,kk);
}

static double wind(double x,int setwf)
{
  //////
  // Window function. setwf=0 -> top hat, setwf=1 -> gaussian
  if(setwf==1) {
    return exp(-0.5*x*x);
  }
  else if(setwf==0) {
    if(x<0.1) {
      return 1.-0.1*x*x+0.003571429*x*x*x*x
	-6.61376E-5*x*x*x*x*x*x
	+7.51563E-7*x*x*x*x*x*x*x*x;
    }
    else
      return 3*(sin(x)-x*cos(x))/(x*x*x);
  }
  else if(setwf==2) {
    if(x>1)
      return 0;
    else
      return 1;
  }
  else {
    fprintf(stderr,"CosmoMad: Wrong window function \n");
    exit(1);
  }
}

typedef struct {
  double r;
  double R1;
  double R2;
  int wf1;
  int wf2;
  Csm_pk_params *pkpar;
} xiparam; //Parameters for xi(r) integrals

static double integxiL_O(double kk,void *params)
{
  //////
  // Integrand for xi_L (for oscillatory integration)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  par=(xiparam *)params;
  
  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);
  
  dum=pk_linear(par->pkpar,kk)*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)/xr;

  return dum;
}

static double integxiL_NO(double logk,void *params)
{
  //////
  // integrand for xi_L (including oscillatory term in j0)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  par=(xiparam *)params;
  
  double kk=pow(10,logk);
  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);

  dum=pk_linear(par->pkpar,kk)*kk*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)*csm_j_bessel(0,xr);

  return dum;
}

double csm_xi2p_L(Csm_params *par,double r,
		  double R1,double R2,
		  char *wf1,char *wf2,double errfac)
{
  //////
  // Correlation function between the linear density contrast smoothed
  // with window function (wf1,R1) and with window function (wf2,R2)
  // at two points separated by a distance r:
  //              <delta_(R1,wf1)(x)*delta_(R2,wf2)(x+r)>
  gsl_function integrand;
  double relerrt=1E-6;
  double integral1=0,integral2=0,errintegral;
  xiparam pars;
  double klim,logklim;

  if(r==0) {
    klim=pow(10.,par->pk->logkcamb);
    logklim=par->pk->logkcamb;
  }
  else {
    klim=CSM_SPH_BESSEL_XMAX/r;
    logklim=log10(klim);
  }

  pars.pkpar=par->pk;
  pars.r=r;
  pars.R1=R1;
  pars.R2=R2;
  if(!strcmp(wf1,"Gauss"))
    pars.wf1=1;
  else if(!strcmp(wf1,"TopHat"))
    pars.wf1=0;
  else if(!strcmp(wf1,"SharpK"))
    pars.wf1=2;
  else {
    fprintf(stderr,"CosmoMad: Unknown window function %s \n",wf1);
    exit(1);
  }
  if(!strcmp(wf2,"Gauss"))
    pars.wf2=1;
  else if(!strcmp(wf2,"TopHat"))
    pars.wf2=0;
  else if(!strcmp(wf2,"SharpK"))
    pars.wf2=2;
  else {
    fprintf(stderr,"CosmoMad: Unknown window function %s \n",wf2);
    exit(1);
  }

  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(10000);
  integrand.params=&pars;
  integrand.function=&integxiL_NO;
  
  int stat=gsl_integration_qagil(&integrand,logklim,0,relerrt,10000,w,
				 &integral1,&errintegral);
  int_error_handle(stat,integral1,errintegral);
  integral1*=CSM_TWOPIPIINVLOGTEN;

  if(r>0) {
    gsl_integration_workspace *cw
      =gsl_integration_workspace_alloc(1000);
    gsl_integration_qawo_table *wf
      =gsl_integration_qawo_table_alloc(r,0.1,GSL_INTEG_SINE,100);
    
    integrand.function=&integxiL_O;
    stat=gsl_integration_qawf(&integrand,klim,relerrt,1000,
			      w,cw,wf,&integral2,&errintegral);
    int_error_handle(stat,integral2,errintegral);
    integral2*=CSM_TWOPIPIINV;

    gsl_integration_qawo_table_free(wf);
    gsl_integration_workspace_free(cw);
  }
  gsl_integration_workspace_free(w);

  return integral1+integral2;
}

double csm_sig0_L(Csm_params *par,
		  double R1,double R2,char *wf1,char *wf2)
{
  //////
  // Covariance between the linear density contrast smoothed with 
  // window function (wf1,R1) and with window function (wf2,R2) at
  // the same point:  <delta_(R1,wf1)(x)*delta_(R2,wf2)(x)>  
  return csm_xi2p_L(par,0,R1,R2,wf1,wf2,1);
}

static double integsigv(double logk,void *params)
{
  //////
  // Integrand for sigma_v
  double dum;
  double kk=pow(10,logk);
  Csm_pk_params *pk=(Csm_pk_params *)params;

  dum=0.3333333333333333*CSM_TWOPIPIINVLOGTEN*
    pk_linear(pk,kk)*kk;
  
  return dum;
}

void csm_set_nonlinear_pk(Csm_params *par,char *fnamePkHFIT)
{
  //////
  // Sets non-linearities for Pk
  if(!(par->pk_params_set)) {
    fprintf(stderr,"Linear Pk not set\n");
    exit(1);
  }

  if(!strcmp(fnamePkHFIT,"RPT")) {
    gsl_function integrand;
    double relerrt=1E-4;
    double sigma_v,errsigma_v;

    if(par->flag_verbose) {
      printf("Setting nonlinear contribution using RPT");
      printf(" gaussian smoothing \n");
    }

    gsl_integration_workspace *w
      =gsl_integration_workspace_alloc(1000);
    integrand.params=par->pk;
    
    integrand.function=&integsigv;
    int stat=gsl_integration_qagil(&integrand,par->pk->logkcamb,0,
				   relerrt,1000,w,
				   &sigma_v,&errsigma_v);
    int_error_handle(stat,sigma_v,errsigma_v);
    
    gsl_integration_workspace_free(w);

    par->pk->sigv=sqrt(sigma_v);

    if(par->flag_verbose)
      printf(" sigma_v = %.3lE Mpc/h\n\n",par->pk->sigv);

    par->pk->use_RPT=1;
    par->pk->use_RPT_ss=0;
    if(par->pk->pkNL_spline_set) {
      gsl_interp_accel_free(par->pk->intacc_pkNL);
      gsl_spline_free(par->pk->spline_pkNL);
      par->pk->pkNL_spline_set=0;
    }
  }
  else if(!strcmp(fnamePkHFIT,"RPT_ss")) {
    gsl_function integrand;
    double relerrt=1E-4;
    double sigma_v,errsigma_v;
    
    if(par->flag_verbose) {
      printf("Setting nonlinear contribution using");
      printf(" RPT gaussian smoothing corrected on small scales\n");
    }

    gsl_integration_workspace *w
      =gsl_integration_workspace_alloc(1000);
    integrand.params=par->pk;
    
    integrand.function=&integsigv;
    int stat=gsl_integration_qagil(&integrand,par->pk->logkcamb,0,
				   relerrt,1000,w,
				   &sigma_v,&errsigma_v);
    int_error_handle(stat,sigma_v,errsigma_v);
    
    gsl_integration_workspace_free(w);
    par->pk->sigv=sqrt(sigma_v);
    if(par->flag_verbose)
      printf(" sigma_v = %.3lE Mpc/h\n\n",par->pk->sigv);

    double *pkarr,*pkarr_lin,norm;
    int ii;

    //Create arrays
    pkarr=malloc(par->pk->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    pkarr_lin=malloc(par->pk->numk*sizeof(double));
    if(pkarr_lin==NULL) error_mem_out();
    for(ii=0;ii<par->pk->numk;ii++) { 
      double kk=pow(10.,par->pk->logkarr[ii]);
      pkarr[ii]=EH_Pk(par->bg,kk,par->pk->ns,0);
      pkarr_lin[ii]=pk_linear(par->pk,kk);
    }

    //Normalize smooth pk
    if(par->pk->pkL_spline_set) {
      gsl_spline_free(par->pk->spline_pkL);
      gsl_interp_accel_free(par->pk->intacc_pkL);
      par->pk->pkL_spline_set=0;
    }
    par->pk->intacc_pkL=gsl_interp_accel_alloc();
    par->pk->spline_pkL=gsl_spline_alloc(gsl_interp_cspline,
					 par->pk->numk);
    gsl_spline_init(par->pk->spline_pkL,par->pk->logkarr,
		    pkarr,par->pk->numk);
    par->pk->pkL_spline_set=1;
    norm=par->pk->s8*par->pk->s8/csm_sig0_L(par,8,8,"TopHat","TopHat");
    for(ii=0;ii<par->pk->numk;ii++)
      pkarr[ii]*=norm;
    if(par->pk->pkL_spline_set) {
      gsl_spline_free(par->pk->spline_pkL);
      gsl_interp_accel_free(par->pk->intacc_pkL);
      par->pk->pkL_spline_set=0;
    }
    par->pk->intacc_pkL=gsl_interp_accel_alloc();
    par->pk->spline_pkL=gsl_spline_alloc(gsl_interp_cspline,
					 par->pk->numk);
    gsl_spline_init(par->pk->spline_pkL,par->pk->logkarr,
		    pkarr_lin,par->pk->numk);
    par->pk->pkL_spline_set=1;

    //Initialize smooth pk spline
    if(par->pk->pkNL_spline_set) {
      gsl_interp_accel_free(par->pk->intacc_pkNL);
      gsl_spline_free(par->pk->spline_pkNL);
      par->pk->pkNL_spline_set=0;
    }
    par->pk->logkminNL=par->pk->logkarr[0];
    par->pk->logkcambNL=par->pk->logkarr[par->pk->numk-1];
    par->pk->intacc_pkNL=gsl_interp_accel_alloc();
    par->pk->spline_pkNL=gsl_spline_alloc(gsl_interp_cspline,
					  par->pk->numk);
    gsl_spline_init(par->pk->spline_pkNL,par->pk->logkarr,
		    pkarr,par->pk->numk);
    par->pk->pkNL_spline_set=1;

    free(pkarr);
    free(pkarr_lin);

    par->pk->use_RPT=0;
    par->pk->use_RPT_ss=1;
  }
  else {
    int ii,nk;
    double kk,ppk;
    double *lkarr,*pkarr;
    FILE *fpk;
    
    if(par->flag_verbose)
      printf("Reading non-linear P_k from file: %s\n",fnamePkHFIT);
    fpk=fopen(fnamePkHFIT,"r");
    if(fpk==NULL) {
      fprintf(stderr,"CosmoMad: Error opening file %s \n",fnamePkHFIT);
      exit(1);
    }
    nk=linecount(fpk);
    lkarr=malloc(nk*sizeof(double));
    if(lkarr==NULL) error_mem_out();
    pkarr=malloc(nk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    rewind(fpk);
    for(ii=0;ii<nk;ii++) {
      int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
      if(stat!=2) {
	fprintf(stderr,"CosmoMad: Error reading file %s, line %d \n",
		fnamePkHFIT,ii+1);
	exit(1);
      }
      pkarr[ii]=ppk*par->pk->norm;
      lkarr[ii]=log10(kk); //log(k) in h Mpc^-1
    }
    par->pk->logkminNL=lkarr[0];
    par->pk->logkcambNL=lkarr[nk-1];
    fclose(fpk);
    
    if(par->pk->pkNL_spline_set) {
      gsl_interp_accel_free(par->pk->intacc_pkNL);
      gsl_spline_free(par->pk->spline_pkNL);
      par->pk->pkNL_spline_set=0;
    }
    par->pk->intacc_pkNL=gsl_interp_accel_alloc();
    par->pk->spline_pkNL=gsl_spline_alloc(gsl_interp_cspline,nk);
    gsl_spline_init(par->pk->spline_pkNL,lkarr,pkarr,nk);
    par->pk->pkNL_spline_set=1;

    free(pkarr);
    free(lkarr);

    par->pk->use_RPT=0;
    par->pk->use_RPT_ss=0;
  }

  par->pk->pkNL_set=1;
}

void csm_set_Pk_params(Csm_params *par,double fg,double gf,
		       double bias,int l_max)
{
  //////
  // Sets beta, gf==D(a) and bias for Pk  
  if(!(par->pk_params_set)) {
    fprintf(stderr,"Pk is not set \n");
    exit(1);
  }
  if(l_max%2!=0) {
    fprintf(stderr,"CosmoMad: multipole number must be even\n");
    exit(1);
  }

  par->pk->fg=fg;
  par->pk->gf2=gf*gf;
  par->pk->bias=bias;
  par->pk->l_multi_max=l_max;
}

static double pk_nonlinear(Csm_pk_params *pkpar,double kk)
{
  //////
  // Real-space, non-linear matter power spectrum
  // normalized with gf given in csm_set_Pk_params
  // but without bias
  double pkNL;
  
  if(!pkpar->pkNL_set) {
    pkNL=pk_linear(pkpar,kk);
  }
  else {
    double lgk=log10(kk);
    
    if(pkpar->use_RPT) {
      double pkL=pk_linear(pkpar,kk);
      pkNL=pkL*exp(-kk*kk*pkpar->sigv*pkpar->sigv*pkpar->gf2);
    }
    else if(pkpar->use_RPT_ss) {
      double pkL=pk_linear(pkpar,kk);
      double pkS;
      if(lgk<pkpar->logkminNL) {
	double pk0;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkminNL,
				   pkpar->intacc_pkNL,&pk0);
	int_error_handle(stat,pkpar->logkminNL,pk0);
	pkS=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkminNL));
      }
      else if(lgk<pkpar->logkcambNL) {
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,lgk,
				   pkpar->intacc_pkNL,&pkS);
	int_error_handle(stat,lgk,pkS);
      }
      else {
	double pkf;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkcambNL,
				   pkpar->intacc_pkNL,&pkf);
	int_error_handle(stat,pkpar->logkcambNL,pkf);
	pkS=pkf*pow(10,-3*(lgk-pkpar->logkcambNL));
      }
      
      pkNL=(pkL-pkS)*exp(-kk*kk*pkpar->sigv*pkpar->sigv*pkpar->gf2)+pkS;
    }
    else {
      if(lgk<pkpar->logkminNL) {
	double pk0;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkminNL,
				   pkpar->intacc_pkNL,&pk0);
	int_error_handle(stat,pkpar->logkminNL,pk0);
	pkNL=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkminNL));
      }
      else if(lgk<pkpar->logkcambNL) {
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,lgk,
				   pkpar->intacc_pkNL,&pkNL);
	int_error_handle(stat,lgk,pkNL);
      }
      else {
	double pkf;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkcambNL,
				   pkpar->intacc_pkNL,&pkf);
	int_error_handle(stat,pkpar->logkcambNL,pkf);
	pkNL=pkf*pow(10,-3*(lgk-pkpar->logkcambNL));
      }
    }
  }

  return pkpar->gf2*pkNL;
}

double csm_Pk_nonlinear(Csm_params *par,double kk)
{
  //////
  // Real-space, non-linear matter power spectrum
  // normalized with gf given in csm_set_Pk_params
  // but without b  
  return pk_nonlinear(par->pk,kk);
}

static double pk_full(Csm_pk_params *pkpar,double kk,double muk)
{
  //////
  // Full power spectrum (with RSDs, bias and non-linearities).
  double pkNL,pk;
  double kaiser;
  double b=pkpar->bias;

  pkNL=pk_nonlinear(pkpar,kk);
  kaiser=(1+pkpar->fg/b*muk*muk);
  pk=kaiser*kaiser*pkNL;

  return b*b*pk;
}

double csm_Pk_full(Csm_params *par,double kk,double muk)
{
  //////
  // Full power spectrum (with RSDs, bias and non-linearities).
  return pk_full(par->pk,kk,muk);
}

double csm_p_leg(int l,double x)
{
  //////
  // Legendre polynomials
  gsl_sf_result res;

  int stat=gsl_sf_legendre_Pl_e(l,x,&res);
  int_error_handle(stat,res.val,res.err);

  return res.val;
}

typedef struct {
  double kk;
  int l;
  Csm_pk_params *pkpar;
} pkl_par; //Parameters for Pk multipoles integrals

static double integ_pkl(double muk,void *params)
{
  //////
  // Integrand for P^l(k)
  pkl_par *par=(pkl_par *)params;
  double kk=par->kk;
  int l=par->l;

  return (l+0.5)*csm_p_leg(l,muk)*pk_full(par->pkpar,kk,muk);
}

static double pk_multipole(Csm_pk_params *pkpar,double kk,int l)
{
  //////
  // Returns multipole l of Pk (RSDs, gf and bias included)
  double pk;
  if(pkpar->pkmulti_spline_set) {
    int il=l/2;
    double logk=log10(kk);
    if((l%2!=0)||(l>pkpar->l_multi_max)) {
      fprintf(stderr,"CosmoMad: multipoles must be even and < %d\n",
	      pkpar->l_multi_max);
      exit(1);
    }

    if(logk<pkpar->logkarr[0]) {
      double pk0;
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],
				 pkpar->logkarr[0],
				 pkpar->intacc_pkmulti[il],&pk0);
      int_error_handle(stat,pkpar->logkarr[0],pk0);
      pk=pk0*pow(10,pkpar->ns*(logk-pkpar->logkarr[0]));
    }
    else if(logk<pkpar->logkarr[pkpar->numk-1]) {
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],logk,
				 pkpar->intacc_pkmulti[il],&pk);
      int_error_handle(stat,logk,pk);
    }
    else {
      double pkf;
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],
				 pkpar->logkarr[pkpar->numk-1],
				 pkpar->intacc_pkmulti[il],&pkf);
      int_error_handle(stat,pkpar->logkarr[pkpar->numk-1],pkf);
      pk=pkf*pow(10,-3*(logk-pkpar->logkarr[pkpar->numk-1]));
    }
  }
  else {
    gsl_function integrand;
    pkl_par pars;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    int status;
    double relerr=1E-4;
    double integral,errintegral;

    pars.kk=kk;
    pars.l=l;
    pars.pkpar=pkpar;
    integrand.function=&integ_pkl;
    integrand.params=&pars;
    
    status=gsl_integration_qag(&integrand,-1,1,0,relerr,1000,
			       GSL_INTEG_GAUSS41,w,&integral,
			       &errintegral);
    
    int_error_handle(status,integral,errintegral);
    gsl_integration_workspace_free(w);

    pk=integral;
  }
  
  return pk;
}

double csm_Pk_multipole(Csm_params *par,double kk,int l)
{
  return pk_multipole(par->pk,kk,l);
}

static void set_pk_multipoles(Csm_pk_params *pkpar)
{
  //////
  // Initializes splines for fast P_l(k) calculation
  // This function is only used if the xi_l(r) are
  // calculated
  int ii;
  int nl=1+pkpar->l_multi_max/2;
  
  if(pkpar->pkmulti_spline_set) {
    for(ii=0;ii<nl;ii++) {
      gsl_interp_accel_free((pkpar->intacc_pkmulti)[ii]);
      gsl_spline_free((pkpar->spline_pkmulti)[ii]);
    }
    free(pkpar->intacc_pkmulti);
    free(pkpar->spline_pkmulti);
    pkpar->pkmulti_spline_set=0;
  }

  pkpar->intacc_pkmulti=(gsl_interp_accel **)
    malloc(nl*sizeof(gsl_interp_accel *));
  if(pkpar->intacc_pkmulti==NULL) error_mem_out();
  pkpar->spline_pkmulti=(gsl_spline **)
    malloc(nl*sizeof(gsl_spline *));
  if(pkpar->spline_pkmulti==NULL) error_mem_out();

  for(ii=0;ii<nl;ii++) {
    int jj;

    double *pkarr=(double *)malloc(pkpar->numk*sizeof(double));
    if(pkarr==NULL) error_mem_out();
    for(jj=0;jj<pkpar->numk;jj++) {
      double kk=pow(10,pkpar->logkarr[jj]);
      pkarr[jj]=pk_multipole(pkpar,kk,2*ii);
    }
    pkpar->intacc_pkmulti[ii]=gsl_interp_accel_alloc();
    pkpar->spline_pkmulti[ii]=gsl_spline_alloc(gsl_interp_cspline,
					       pkpar->numk);
    gsl_spline_init(pkpar->spline_pkmulti[ii],pkpar->logkarr,
		    pkarr,pkpar->numk);
    free(pkarr);
  }

  pkpar->pkmulti_spline_set=1;
}


/****************************/
/*   Correlation Function   */
/****************************/
typedef struct {
  double rr;
  int l;
  Csm_pk_params *pkpar;
} xil_par; //Params for xi_l(r) integrals

static double integxi_l_NO(double logk,void *params)
{
  //////
  // Integrand for xi^l(r) including all oscillatory
  // terms for logarithmic integration
  double pk,dum;
  xil_par *par=(xil_par *)params;
  double kk=pow(10,logk);
  double xr=kk*par->rr;
  int l=par->l;
  double jbes=csm_j_bessel(l,xr);

  pk=pk_multipole(par->pkpar,kk,l);
  dum=pk*kk*kk*kk*jbes;

  return dum;
}

static double integxi_l_O(double kk,void *params)
{
  //////
  // Integrand for xi^l(r) without oscillatory terms
  // for oscillatory integration
  double pk,dum;
  xil_par *par=(xil_par *)params;
  double xr=kk*par->rr;
  int sign=1;

  if(par->l%4==2)
    sign=-1;

  pk=pk_multipole(par->pkpar,kk,par->l);
  dum=sign*pk*kk*kk/xr;

  return dum;
}

double csm_xi_multipole(Csm_params *par,double rr,int l)
{
  //////
  // l-th multipole of the correlation function
  if(par->xi_spline_set) {
    if(l>par->pk->l_multi_max) {
      fprintf(stderr,"CosmoMad: Multipole l=%d not calculated \n",l);
      exit(1);
    }
    if(l%2!=0) {
      fprintf(stderr,"CosmoMad: multipoles must be even \n");
      exit(1);
    }
    int il=l/2;
    if(rr<CSM_RMIN_LOG) return par->xi->ximultimin[il];
    else if(rr<CSM_RMIN_LIN) {
      double result;
      int stat;
      double logr=log10(rr);
      stat=gsl_spline_eval_e(par->xi->spline_ximulti_log[il],logr,
			     par->xi->intacc_ximulti_log[il],&result);
      int_error_handle(stat,CSM_RMIN_LOG,logr);
      return result;
    }
    else if(rr<CSM_RMAX_LIN) {
      double result;
      int stat;
      stat=gsl_spline_eval_e(par->xi->spline_ximulti_lin[il],rr,
			     par->xi->intacc_ximulti_lin[il],&result);
      int_error_handle(stat,CSM_RMAX_LIN,rr);

      return result;
    }
    else {
      return 0;
      //      fprintf(stderr,"CosmoMad: r=%lE is out of range \n",rr);
      //      exit(1);
    }
  }
  else {
    if(!(par->pk_params_set)) {
      fprintf(stderr,"CosmoMad: must set Pk params ");
      fprintf(stderr,"before calculating xi!\n");
      exit(1);
    }
    if(!par->pk->pkmulti_spline_set) {
      if(par->flag_verbose)
	printf("Setting Pk multipoles\n");
      set_pk_multipoles(par->pk);
    }
    
    gsl_function integrand;
    double relerrt=1E-4;
    double integral1=0,integral2=0,errintegral;
    double klim,logklim;
    int status;
    xil_par pars;
    int sign;

    if(rr==0) {
      klim=pow(10.,par->pk->logkcamb);
      logklim=par->pk->logkcamb;
    }
    else {
      klim=CSM_SPH_BESSEL_XMAX/rr;
      logklim=log10(klim);
    }

    if(l%4==0) sign=1;
    else if(l%4==2) sign=-1;
    else {
      fprintf(stderr,"CosmoMad: multipoles must be even \n");
      exit(1);
    }

    pars.rr=rr;
    pars.l=l;
    pars.pkpar=par->pk;

    gsl_integration_workspace *w
      =gsl_integration_workspace_alloc(1000);
    integrand.params=&pars;
    integrand.function=integxi_l_NO;

    status=gsl_integration_qagil(&integrand,logklim,0,relerrt,1000,
    				 w,&integral1,&errintegral);
    int_error_handle(status,integral1,errintegral);
    integral1*=CSM_TWOPIPIINVLOGTEN;

    if(rr>0) {
      gsl_integration_workspace *cw
	=gsl_integration_workspace_alloc(1000);
      gsl_integration_qawo_table *wf
	=gsl_integration_qawo_table_alloc(rr,0.1,GSL_INTEG_SINE,100);
      
      integrand.function=&integxi_l_O;
      status=gsl_integration_qawf(&integrand,klim,relerrt,1000,
				  w,cw,wf,&integral2,&errintegral);
      int_error_handle(status,integral2,errintegral);
      integral2*=CSM_TWOPIPIINV;
      
      gsl_integration_qawo_table_free(wf);
      gsl_integration_workspace_free(cw);
    }
    gsl_integration_workspace_free(w);
    
    return sign*(integral1+integral2);
  }
}

void csm_set_xi_multipole_splines(Csm_params *par)
{
  //////
  // Initializes all splines for correlation function multipoles.
  if(par->xi_spline_set)
    csm_unset_xi_multipole_splines(par);

  int ii;
  int nl=1+par->pk->l_multi_max/2;
  double r_min_log=CSM_RMIN_LOG_SPLINE;
  double r_max_log=CSM_RMAX_LOG_SPLINE;
  double logr_min_log=log10(r_min_log);
  double logr_max_log=log10(r_max_log);
  double dlogr_log=CSM_DLOG;
  double r_min_lin=CSM_RMIN_LIN_SPLINE;
  double r_max_lin=CSM_RMAX_LIN_SPLINE;
  double dr_lin=CSM_DR;
  int nr_log=(int)((logr_max_log-logr_min_log)/dlogr_log);
  int nr_lin=(int)((r_max_lin-r_min_lin)/dr_lin);

  par->xi=csm_xi_params_new();

  if(par->flag_verbose)
    printf("Setting splines for xi \n\n");

  par->xi->intacc_ximulti_lin=(gsl_interp_accel **)
    malloc(nl*sizeof(gsl_interp_accel *));
  if(par->xi->intacc_ximulti_lin==NULL) error_mem_out();
  par->xi->intacc_ximulti_log=(gsl_interp_accel **)
    malloc(nl*sizeof(gsl_interp_accel *));
  if(par->xi->intacc_ximulti_log==NULL) error_mem_out();
  par->xi->spline_ximulti_lin=(gsl_spline **)
    malloc(nl*sizeof(gsl_spline *));
  if(par->xi->spline_ximulti_lin==NULL) error_mem_out();
  par->xi->spline_ximulti_log=(gsl_spline **)
    malloc(nl*sizeof(gsl_spline *));
  if(par->xi->spline_ximulti_log==NULL) error_mem_out();
  par->xi->ximultimin=(double *)malloc(nl*sizeof(double));
  if(par->xi->ximultimin==NULL) error_mem_out();

  for(ii=0;ii<nl;ii++) {
    int jj;
    
    par->xi->ximultimin[ii]=csm_xi_multipole(par,CSM_RMIN_LOG,2*ii);

    //Set array for xi_lin
    double *xi_arr_lin=(double *)malloc(nr_lin*sizeof(double));
    double *rarr_lin=(double *)malloc(nr_lin*sizeof(double));
    if(rarr_lin==NULL) error_mem_out();
    if(xi_arr_lin==NULL) error_mem_out();
    for(jj=0;jj<nr_lin;jj++) {
      double rr;
      rarr_lin[jj]=r_min_lin+jj*(r_max_lin-r_min_lin)/nr_lin;
      rr=rarr_lin[jj];
      xi_arr_lin[jj]=csm_xi_multipole(par,rr,2*ii);
    }
    par->xi->intacc_ximulti_lin[ii]=gsl_interp_accel_alloc();
    par->xi->spline_ximulti_lin[ii]=gsl_spline_alloc(gsl_interp_cspline,
						     nr_lin);
    gsl_spline_init(par->xi->spline_ximulti_lin[ii],
		    rarr_lin,xi_arr_lin,nr_lin);
    free(xi_arr_lin);
    free(rarr_lin);
    
    //Set array for xi_log
    double *xi_arr_log=(double *)malloc(nr_log*sizeof(double));
    double *lograrr_log=(double *)malloc(nr_log*sizeof(double));
    if(lograrr_log==NULL) error_mem_out();
    if(xi_arr_log==NULL) error_mem_out();
    for(jj=0;jj<nr_log;jj++) {
      double logr,rr;
      logr=logr_min_log+jj*(logr_max_log-logr_min_log)/nr_log;
      lograrr_log[jj]=logr;
      rr=pow(10.,logr);
      //      printf("%lE \n",logr);
      xi_arr_log[jj]=csm_xi_multipole(par,rr,2*ii);
    }
    par->xi->intacc_ximulti_log[ii]=gsl_interp_accel_alloc();
    par->xi->spline_ximulti_log[ii]=gsl_spline_alloc(gsl_interp_cspline,
						     nr_log);
    gsl_spline_init(par->xi->spline_ximulti_log[ii],
		    lograrr_log,xi_arr_log,nr_log);
    free(xi_arr_log);
    free(lograrr_log);
  }
 
  par->xi_spline_set=1;
}

void csm_unset_xi_multipole_splines(Csm_params *par)
{
  //////
  // Frees al memory associated with the xi_l(r) splines
  // (undoes csm_set_xi_multipole_splines
  if(par->xi_spline_set)
    csm_xi_params_free(par->xi,par->pk->l_multi_max/2+1);
  par->xi_spline_set=0;
}

double csm_xi_3D(Csm_params *par,double rr,double mu)
{
  //////
  // Returns the 3-D correlation function at radius rr
  // and cos(theta)==mu
  double result=0;
  int ii;

  for(ii=0;ii<par->pk->l_multi_max/2+1;ii++) {
    double xi_l=csm_xi_multipole(par,rr,2*ii);
    double p_l=csm_p_leg(2*ii,mu);
    result+=p_l*xi_l;
  }

  return result;
}

typedef struct {
  double k_tra;
  double r_par;
  double r_tra;
  double r_mod;
  Csm_pk_params *pkpar;
} xiPS_par; //Params for the integrals of xi(pi,sigma)

static double integxi_PS_O_level0(double k_par,void *params)
{
  //////
  // Integrand for xi(pi,sigma) (for oscillatory integration)  
  double pk;
  xiPS_par *par=(xiPS_par *)params;
  double k_tra=par->k_tra;
  double k_mod=sqrt(k_par*k_par+k_tra*k_tra);

  if(k_mod==0)
    pk=0;
  else
    pk=pk_full(par->pkpar,k_mod,k_par/k_mod);
  
  return pk;
}

static double integxi_PS_NO_level1(double log_k_tra,void *params)
{
  //////
  // Integrand for xi(pi,sigma) (for logarithmic integration)
  gsl_function integrand;
  double relerrt=1E-4;
  xiPS_par *par=(xiPS_par *)params;
  double integral,errintegral;
  double k_tra=pow(10.,log_k_tra);
  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *cw
    =gsl_integration_workspace_alloc(1000);
  gsl_integration_qawo_table *wf
    =gsl_integration_qawo_table_alloc(par->r_par,0.1,
				      GSL_INTEG_COSINE,100);
  double x_tra=k_tra*par->r_tra;
  double J0=gsl_sf_bessel_J0(x_tra);

  par->k_tra=k_tra;

  integrand.params=par;
  integrand.function=&integxi_PS_O_level0;
  int stat=gsl_integration_qawf(&integrand,0,relerrt,1000,
				w,cw,wf,&integral,&errintegral);
  int_error_handle(stat,integral,errintegral);

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(cw);
  gsl_integration_workspace_free(w);

  return integral*J0*k_tra*k_tra;
}

static double integxi_PS_O_level1(double k_tra,void *params)
{
  //////
  // Integrand for xi(pi,sigma)  
  gsl_function integrand;
  double relerrt=1E-4;
  xiPS_par *par=(xiPS_par *)params;
  double integral,errintegral;
  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *cw
    =gsl_integration_workspace_alloc(1000);
  gsl_integration_qawo_table *wf
    =gsl_integration_qawo_table_alloc(par->r_par,0.1,
				      GSL_INTEG_COSINE,100);
  par->k_tra=k_tra;

  integrand.params=par;
  integrand.function=&integxi_PS_O_level0;
  int stat=gsl_integration_qawf(&integrand,0,relerrt,1000,
				w,cw,wf,&integral,&errintegral);
  int_error_handle(stat,integral,errintegral);

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(cw);
  gsl_integration_workspace_free(w);

  return integral*k_tra/sqrt(M_PI*k_tra*par->r_tra);
}

double csm_xi_pi_sigma(Csm_params *par,double pi,double sigma,
		       int use_multipoles)
{
  //////
  // Redshift-space 3D correlation function. If use_multipoles!=0
  // multipoles will be used. Otherwise the slow double integral
  // defined above will be used.
  if(use_multipoles) {
    double r=sqrt(pi*pi+sigma*sigma);
    double mu;
    if(r>0) mu=pi/r;
    else mu=1;

    return csm_xi_3D(par,r,mu);
  }
  else {
    gsl_function integrand;
    double relerrt=1E-4;
    xiPS_par pars;
    double integral1=0,integral2=0,errintegral;
    double klim,logklim;
    pars.r_par=pi;
    pars.r_tra=sigma;
    pars.r_mod=sqrt(pi*pi+sigma*sigma);
    pars.pkpar=par->pk;

    if(sigma==0) {
      klim=pow(10.,par->pk->logkcamb);
      logklim=par->pk->logkcamb;
    }
    else {
      klim=CSM_CYL_BESSEL_XMAX/sigma;
      logklim=log10(klim);
    }

    gsl_integration_workspace *w
      =gsl_integration_workspace_alloc(1000);
    integrand.params=&pars;
    integrand.function=&integxi_PS_NO_level1;
    int stat=gsl_integration_qagil(&integrand,logklim,0,relerrt,1000,w,
    				   &integral1,&errintegral);
    int_error_handle(stat,integral1,errintegral);
    integral1*=CSM_TWOPIPIINVLOGTEN;

    if(sigma>0) {
      gsl_integration_workspace *cw
	=gsl_integration_workspace_alloc(1000);
      gsl_integration_qawo_table *wf;
      double res_sin,res_cos;

      integrand.function=&integxi_PS_O_level1;

      wf=gsl_integration_qawo_table_alloc(sigma,0.1,GSL_INTEG_SINE,100);
      stat=gsl_integration_qawf(&integrand,klim,relerrt,1000,
				  w,cw,wf,&res_sin,&errintegral);
      int_error_handle(stat,res_sin,errintegral);
      gsl_integration_qawo_table_free(wf);

      wf=gsl_integration_qawo_table_alloc(sigma,0.1,
					  GSL_INTEG_COSINE,100);
      stat=gsl_integration_qawf(&integrand,klim,relerrt,1000,
				  w,cw,wf,&res_cos,&errintegral);
      int_error_handle(stat,res_cos,errintegral);
      gsl_integration_qawo_table_free(wf);
      
      integral2=CSM_TWOPIPIINV*(res_cos+res_sin);
      gsl_integration_workspace_free(cw);
    }
    gsl_integration_workspace_free(w);

    return integral1+integral2;
  }
}

double csm_M2R(Csm_params *par,double mass)
{
  return pow(mass/(CSM_FOURPITHIRD*(par->bg->OM)*CSM_RHOCRIT),
	     0.3333333333333333);
}

double csm_R2M(Csm_params *par,double radius)
{
  return CSM_FOURPITHIRD*(par->bg->OM)*CSM_RHOCRIT*radius*radius*radius;
}

void csm_set_mf_params(Csm_params *par,double lmmn,double lmmx,double dlm)
{
  int im,nm;
  double *lmarr,*sigmarr,*dsigmarr;

  if(par->mf_params_set)
    csm_mf_params_free(par->mf);
  par->mf_params_set=0;

  if(!(par->pk_params_set)) {
    fprintf(stderr,"You must set pk before you set mf\n");
    exit(1);
  }

  par->mf=csm_mf_params_new();

  //Set mass/sigma arrays
  nm=(int)((lmmx-lmmn)/dlm)+1;
  lmarr=malloc(nm*sizeof(double));
  if(lmarr==NULL) error_mem_out();
  sigmarr=malloc(nm*sizeof(double));
  if(sigmarr==NULL) error_mem_out();
  dsigmarr=malloc(nm*sizeof(double));
  if(dsigmarr==NULL) error_mem_out();
  
  for(im=0;im<nm;im++) {
    double mass,rad;
    lmarr[im]=lmmn+dlm*im;
    mass=pow(10.,lmarr[im]);
    rad=csm_M2R(par,mass);
    sigmarr[im]=sqrt(csm_sig0_L(par,rad,rad,"TopHat","TopHat"));
  }

  //Numerically differentiate sigma(M)
  dsigmarr[0]=-log10(sigmarr[1]/sigmarr[0])/dlm;
  for(im=1;im<nm-1;im++)
    dsigmarr[im]=-log10(sigmarr[im+1]/sigmarr[im-1])/(2*dlm);
  dsigmarr[nm-1]=-log10(sigmarr[nm-1]/sigmarr[nm-2])/dlm;

  //Initialize all parameters
  par->mf->logmmin=lmarr[0];
  par->mf->logmmax=lmarr[nm-1];

  par->mf->intacc_sM=gsl_interp_accel_alloc();
  par->mf->spline_sM=gsl_spline_alloc(gsl_interp_cspline,nm);
  gsl_spline_init(par->mf->spline_sM,lmarr,sigmarr,nm);
  par->mf->sM_spline_set=1;

  par->mf->intacc_dsM=gsl_interp_accel_alloc();
  par->mf->spline_dsM=gsl_spline_alloc(gsl_interp_cspline,nm);
  gsl_spline_init(par->mf->spline_dsM,lmarr,dsigmarr,nm);
  par->mf->dsM_spline_set=1;

  par->mf_params_set=1;
}

double csm_sigmaM(Csm_params *par,double m)
{
  int stat;
  double sM;
  double lm=log10(m);

  if(!(par->mf_params_set)) {
    fprintf(stderr,"You have to set the mass-function params\n");
    exit(1);
  }    

  if((lm<par->mf->logmmin) || (lm>par->mf->logmmax)) {
    fprintf(stderr,"Mass %lE is out of range\n",m);
    exit(1);
  }

  stat=gsl_spline_eval_e(par->mf->spline_sM,lm,par->mf->intacc_sM,&sM);
  int_error_handle(stat,lm,sM);

  return sM;
}

double csm_dlsigMdlM(Csm_params *par,double m)
{
  int stat;
  double dsM;
  double lm=log10(m);

  if(!(par->mf_params_set)) {
    fprintf(stderr,"You have to set the mass-function params\n");
    exit(1);
  }    

  if((lm<par->mf->logmmin) || (lm>par->mf->logmmax)) {
    fprintf(stderr,"Mass %lE is out of range\n",m);
    exit(1);
  }

  stat=gsl_spline_eval_e(par->mf->spline_dsM,lm,par->mf->intacc_dsM,&dsM);
  int_error_handle(stat,lm,dsM);

  return dsM;
}

double csm_multiplicity_function(Csm_params *par,double mass,double z,char *mftype)
{
  double fM,nu;
  double gf=csm_growth_factor(par,1./(1+z))/csm_growth_factor(par,1);
  double sigma=gf*csm_sigmaM(par,mass);

  if(!strcmp(mftype,"PS")) {
    nu=CSM_DELTA_C/sigma;
    fM=sqrt(2/M_PI)*nu*exp(-0.5*nu*nu);
  }
  else if(!strcmp(mftype,"ST")) {
    nu=CSM_DELTA_C/sigma;
    fM=CSM_ST_A*sqrt(2*CSM_ST_a/M_PI)*(1+pow(1./(CSM_ST_a*nu*nu),CSM_ST_p))*nu*exp(-CSM_ST_a*nu*nu*0.5);
  }
  else if(!strcmp(mftype,"JP")) {
    nu=CSM_DELTA_C/sigma;
    fM=(CSM_JP_a*CSM_JP_b*pow(nu,CSM_JP_b)+2*CSM_JP_c*nu*nu*(1+CSM_JP_a*pow(nu,CSM_JP_b)))*
      exp(-CSM_JP_c*nu*nu)/pow(1+CSM_JP_a*pow(nu,CSM_JP_b),2);
  }
  else if(!strcmp(mftype,"Jenkins")) {
    double x=pow(fabs(log(1./sigma)+CSM_JEN_b),CSM_JEN_q);
    fM=CSM_JEN_B*exp(-x);
  }
  else if(!strcmp(mftype,"Warren")) {
    fM=CSM_WAR_A*(1./pow(sigma,CSM_WAR_a)+CSM_WAR_b)*exp(-CSM_WAR_c/(sigma*sigma));
  }
  else if(!strcmp(mftype,"Tinker200")) {
    double A=CSM_TINKER_A_200*pow(1+z,CSM_TINKER_Aexp_200);
    double a=CSM_TINKER_a_200*pow(1+z,CSM_TINKER_aexp_200);
    double b=CSM_TINKER_b_200*pow(1+z,CSM_TINKER_bexp_200);
    double c=CSM_TINKER_c_200;
    fM=A*(1+pow(b/sigma,a))*exp(-c/(sigma*sigma));
  }
  else if(!strcmp(mftype,"Tinker500")) {
    double A=CSM_TINKER_A_500*pow(1+z,CSM_TINKER_Aexp_500);
    double a=CSM_TINKER_a_500*pow(1+z,CSM_TINKER_aexp_500);
    double b=CSM_TINKER_b_500*pow(1+z,CSM_TINKER_bexp_500);
    double c=CSM_TINKER_c_500;
    fM=A*(1+pow(b/sigma,a))*exp(-c/(sigma*sigma));
  }
  else if(!strcmp(mftype,"Watson")) {
    fM=CSM_WATSON_A*(1+pow(CSM_WATSON_BETA/sigma,CSM_WATSON_ALPHA))*exp(-CSM_WATSON_GAMMA/(sigma*sigma));
  }
  else {
    fprintf(stderr,"Wrong mass function type %s\n",mftype);
    exit(1);
  }

  fM*=csm_dlsigMdlM(par,mass);
  
  return fM;
}

double csm_mass_function_logarithmic(Csm_params *par,double mass,double z,char *mftype)
{
  double fM=csm_multiplicity_function(par,mass,z,mftype);
  return M_LN10*fM*par->bg->OM*CSM_RHOCRIT/mass;
}
