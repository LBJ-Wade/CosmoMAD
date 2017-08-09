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

#include "cosmo_mad.h"

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
    csm_int_error_handle(stat,integral,errintegral);
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
    csm_int_error_handle(stat,integral,errintegral);
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
    if(status!=GSL_SUCCESS)
      csm_report_error(1,"CosmoMad: ODE didn't converge\n");
    
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
  csm_int_error_handle(stat,t,result);
  
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

  if(par->bg->OM<=0)
    csm_report_error(1,"Wrong matter parameter %.3lf \n",par->bg->OM);
  if(par->bg->OM<par->bg->OB)
    csm_report_error(1,"Wrong M/B parameter %.3lf > %.3lf \n",par->bg->OB,par->bg->OM);
  if(par->bg->w0>-0.333333)
    csm_report_error(1,"DE is too exotic (w=%.3lf \n",par->bg->w0);
  if(par->bg->TCMB<0)
    csm_report_error(1,"Wrong CMB temperature %.3lf \n",par->bg->TCMB);

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
