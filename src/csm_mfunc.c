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

  if(!(par->pk_params_set))
    csm_report_error(1,"You must set pk before you set mf\n");

  par->mf=csm_mf_params_new();

  //Set mass/sigma arrays
  nm=(int)((lmmx-lmmn)/dlm)+1;
  lmarr=csm_malloc(nm*sizeof(double));
  sigmarr=csm_malloc(nm*sizeof(double));
  dsigmarr=csm_malloc(nm*sizeof(double));
  
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

  if(!(par->mf_params_set))
    csm_report_error(1,"You have to set the mass-function params\n");

  if((lm<par->mf->logmmin) || (lm>par->mf->logmmax))
    csm_report_error(1,"Mass %lE is out of range\n",m);

  stat=gsl_spline_eval_e(par->mf->spline_sM,lm,par->mf->intacc_sM,&sM);
  csm_int_error_handle(stat,lm,sM);

  return sM;
}

double csm_dlsigMdlM(Csm_params *par,double m)
{
  int stat;
  double dsM;
  double lm=log10(m);

  if(!(par->mf_params_set))
    csm_report_error(1,"You have to set the mass-function params\n");

  if((lm<par->mf->logmmin) || (lm>par->mf->logmmax))
    csm_report_error(1,"Mass %lE is out of range\n",m);

  stat=gsl_spline_eval_e(par->mf->spline_dsM,lm,par->mf->intacc_dsM,&dsM);
  csm_int_error_handle(stat,lm,dsM);

  return dsM;
}

double csm_halo_bias(Csm_params *par,double mass,double z,char *mftype)
{
  double nu,bM;
  double gf=csm_growth_factor(par,1./(1+z))/csm_growth_factor(par,1);
  double sigma=gf*csm_sigmaM(par,mass);

  if(!strcmp(mftype,"PS")) {
    nu=CSM_DELTA_C/sigma;
    bM=1+(nu*nu-1)/CSM_DELTA_C;
  }
  else if(!strcmp(mftype,"ST")) {
    double a=0.707;
    double sqa=sqrt(a);
    double an2=a*nu*nu;
    double b=0.5;
    double c=0.6;
    bM=1.+(sqa*an2+sqa*b*pow(an2,1-c)-
	   pow(an2,c)/(pow(an2,c)+b*(1-c)*(1-0.5*c)))/(sqa*CSM_DELTA_C);
  }
  else if((!strcmp(mftype,"Tinker10_200")) ||
	  (!strcmp(mftype,"Tinker10_500"))) {
    double y,ey,fit_A,fit_a,fit_B,fit_b,fit_C,fit_c;
    if(!strcmp(mftype,"Tinker10_200"))
      y=2.3010299956639813; //log10(200)
    else
      y=2.6989700043360187; //log10(500)
    ey=exp(-pow(4./y,4.));
    fit_A=1.0+0.24*y*ey;
    fit_a=0.44*y-0.88;
    fit_B=0.183;
    fit_b=1.5;
    fit_C=0.019+0.107*y+0.19*ey;
    fit_c=2.4;
    nu=CSM_DELTA_C/sigma;
    
    bM=1.-fit_A*pow(nu,fit_a)/(pow(nu,fit_a)+pow(CSM_DELTA_C,fit_a))+
      fit_B*pow(nu,fit_b)+fit_C*pow(nu,fit_c);
  }    
  else
    csm_report_error(1,"Wrong mass function type %s\n",mftype);

  return bM;
}

double csm_multiplicity_function(Csm_params *par,double mass,double z,
				 char *mftype)
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
  else if(!strcmp(mftype,"Tinker10_200")) {
    double alpha=CSM_TINKER10_ALPHA_200;
    double beta=CSM_TINKER10_BETA_200*pow(1+z,CSM_TINKER10_BETAexp);
    double gamma=CSM_TINKER10_GAMMA_200*pow(1+z,CSM_TINKER10_GAMMAexp);
    double phi=CSM_TINKER10_PHI_200*pow(1+z,CSM_TINKER10_PHIexp);
    double eta=CSM_TINKER10_ETA_200*pow(1+z,CSM_TINKER10_ETAexp);
    nu=CSM_DELTA_C/sigma;
    fM=alpha*(1+pow(beta*nu,-2*phi))*pow(nu,2*eta)*exp(-0.5*gamma*nu*nu);
  }
  else if(!strcmp(mftype,"Tinker10_500")) {
    double alpha=CSM_TINKER10_ALPHA_500;
    double beta=CSM_TINKER10_BETA_500*pow(1+z,CSM_TINKER10_BETAexp);
    double gamma=CSM_TINKER10_GAMMA_500*pow(1+z,CSM_TINKER10_GAMMAexp);
    double phi=CSM_TINKER10_PHI_500*pow(1+z,CSM_TINKER10_PHIexp);
    double eta=CSM_TINKER10_ETA_500*pow(1+z,CSM_TINKER10_ETAexp);
    nu=CSM_DELTA_C/sigma;
    fM=alpha*(1+pow(beta*nu,-2*phi))*pow(nu,2*eta)*exp(-0.5*gamma*nu*nu);
  }
  else if(!strcmp(mftype,"Watson")) {
    fM=CSM_WATSON_A*(1+pow(CSM_WATSON_BETA/sigma,CSM_WATSON_ALPHA))*
      exp(-CSM_WATSON_GAMMA/(sigma*sigma));
  }
  else
    csm_report_error(1,"Wrong mass function type %s\n",mftype);

  fM*=csm_dlsigMdlM(par,mass);
  
  return fM;
}

double csm_mass_function_logarithmic(Csm_params *par,double mass,double z,char *mftype)
{
  double fM=csm_multiplicity_function(par,mass,z,mftype);
  return M_LN10*fM*par->bg->OM*CSM_RHOCRIT/mass;
}

