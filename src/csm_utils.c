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

int csm_linecount(FILE *f)
{
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void csm_report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr," CosmoMad, Fatal error: %s",msg);
    exit(level);
  }
  else
    fprintf(stderr," CosmoMad, Warning: %s",msg);
}

void csm_int_error_handle(int status,double result,double error)
{
  //////
  // Error handler for gsl
  if(isnan(result))
    csm_report_error(1,"NAN found \n");
  else {
    if(status==GSL_EROUND)
      csm_report_error(0,"Roundoff error: %lE %lE \n",result,error);
    else if(status==GSL_EMAXITER)
      csm_report_error(0,"Ran out of iterations: %lE %lE \n",result,error);
    else if(status==GSL_ESING)
      csm_report_error(0,"Singularity found: %lE %lE \n",result,error);
    else if(status==GSL_EDIVERGE)
      csm_report_error(0,"Integral seems to diverge: %lE %lE \n",result,error);
    else if(status==GSL_ETOL)
      csm_report_error(0,"Can't reach tolerance: %lE %lE\n",result,error);
    else if(status==GSL_EUNDRFLW)
      csm_report_error(0,"Underflow: %lE %lE \n",result,error);
    else if(status==GSL_EDOM)
      csm_report_error(1,"Outside interpolation range!! %lE %lE\n",result,error);
    else if(status)
      csm_report_error(1,"Unknown error code %d %lf %lf \n",status,result,error);
  }
}

void *csm_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) csm_report_error(1,"Out of memory\n");

  return outptr;
}

void *csm_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL) csm_report_error(1,"Out of memory\n");

  return outptr;
}

FILE *csm_fopen(const char *path,const char *mode)
{
  FILE *fout=fopen(path,mode);
  if(fout==NULL)
    csm_report_error(1,"Couldn't open file %s\n",path);

  return fout;
}

void csm_unset_gsl_eh(Csm_params *par)
{
  //////
  // Disables GSL default error handler
  par->gsl_error_handler_old=gsl_set_error_handler_off();
}

void csm_set_verbosity(Csm_params *par,int verb)
{
  //////
  // Sets verbosity level
  par->flag_verbose=verb;
}

void csm_bg_params_free(Csm_bg_params *par)
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

Csm_bg_params *csm_bg_params_new(void)
{
  //////
  // bg_params creator
  // Default Planck parameters
  Csm_bg_params *bg=csm_malloc(sizeof(Csm_bg_params));
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

  bg->eh=csm_malloc(sizeof(Csm_EH_params));
  bg->at_spline_set=0;

  return bg;
}

void csm_pk_params_free(Csm_pk_params *par)
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

Csm_pk_params *csm_pk_params_new(void)
{
  //////
  // pk_params creator
  Csm_pk_params *pk=csm_malloc(sizeof(Csm_pk_params));

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

void csm_xi_params_free(Csm_xi_params *par,int nl)
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

Csm_xi_params *csm_xi_params_new(void)
{
  //////
  // xi_params creator
  Csm_xi_params *xi=csm_malloc(sizeof(Csm_xi_params));

  xi->intacc_ximulti_log=NULL;
  xi->intacc_ximulti_lin=NULL;
  xi->spline_ximulti_log=NULL;
  xi->spline_ximulti_lin=NULL;

  return xi;
}

void csm_mf_params_free(Csm_mf_params *par)
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

Csm_mf_params *csm_mf_params_new(void)
{
  //////
  // mf_params creator
  Csm_mf_params *mf=csm_malloc(sizeof(Csm_mf_params));
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
  Csm_params *par=csm_malloc(sizeof(Csm_params));
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

double csm_p_leg(int l,double x)
{
  //////
  // Legendre polynomials
  gsl_sf_result res;

  int stat=gsl_sf_legendre_Pl_e(l,x,&res);
  csm_int_error_handle(stat,res.val,res.err);

  return res.val;
}

#define CSM_GAMMA1 2.6789385347077476336556 //Gamma(1/3)
#define CSM_GAMMA2 1.3541179394264004169452 //Gamma(2/3)
#define CSM_ROOTPI12 21.269446210866192327578 //12*sqrt(pi)

static double sph_bessel(int l,double x)
{
  double jl;
  double ax=fabs(x);
  double ax2=x*x;
  if(l<0) {
    fprintf(stderr,"CosmoMas: l>0 for Bessel function");
    exit(1);
  }

  if(l<7) {
    if(l==0) {
      if(ax<0.1) jl=1-ax2*(1-ax2/20.)/6.;
      else jl=sin(x)/x;
    }
    else if(l==1) {
      if(ax<0.2) jl=ax*(1-ax2*(1-ax2/28)/10)/3;
      else jl=(sin(x)/ax-cos(x))/ax;
    }
    else if(l==2) {
      if(ax<0.3) jl=ax2*(1-ax2*(1-ax2/36)/14)/15;
      else jl=(-3*cos(x)/ax-sin(x)*(1-3/ax2))/ax;
    }
    else if(l==3) {
      if(ax<0.4)
	jl=ax*ax2*(1-ax2*(1-ax2/44)/18)/105;
      else 
	jl=(cos(x)*(1-15/ax2)-sin(x)*(6-15/ax2)/ax)/ax;
    }
    else if(l==4) {
      if(ax<0.6)
	jl=ax2*ax2*(1-ax2*(1-ax2/52)/22)/945;
      else
	jl=(sin(x)*(1-(45-105/ax2)/ax2)+cos(x)*(10-105/ax2)/ax)/ax;
    }
    else if(l==5) {
      if(ax<1.0)
	jl=ax2*ax2*ax*(1-ax2*(1-ax2/60)/26)/10395;
      else {
	jl=(sin(x)*(15-(420-945/ax2)/ax2)/ax-
	    cos(x)*(1-(105-945/ax2)/ax2))/ax;
      }
    }
    else {
      if(ax<1.0)
	jl=ax2*ax2*ax2*(1-ax2*(1-ax2/68)/30)/135135;
      else {
	jl=(sin(x)*(-1+(210-(4725-10395/ax2)/ax2)/ax2)+
	    cos(x)*(-21+(1260-10395/ax2)/ax2)/ax)/ax;
      }
    }
  }
  else {
    double nu=l+0.5;
    double nu2=nu*nu;
    
    if(ax<1.0E-40) jl=0;
    else if((ax2/l)<0.5) {
      jl=(exp(l*log(ax/nu)-M_LN2+nu*(1-M_LN2)-(1-(1-3.5/nu2)/(30*nu2))/(12*nu))/nu)*
	(1-ax2/(4*nu+4)*(1-ax2/(8*nu+16)*(1-ax2/(12*nu+36))));
    }
    else if((l*l/ax)<0.5) {
      double beta=ax-0.5*M_PI*(l+1);
      jl=(cos(beta)*(1-(nu2-0.25)*(nu2-2.25)/(8*ax2)*(1-(nu2-6.25)*(nu2-12.25)/(48*ax2)))-
	  sin(beta)*(nu2-0.25)/(2*ax)*(1-(nu2-2.25)*(nu2-6.25)/(24*ax2)*
				       (1-(nu2-12.25)*(nu2-20.25)/(80*ax2))))/ax;
    }
    else {
      double l3=pow(nu,0.325);
      if(ax<nu-1.31*l3) {
	double cosb=nu/ax;
	double sx=sqrt(nu2-ax2);
	double cotb=nu/sx;
	double secb=ax/nu;
	double beta=log(cosb+sx/ax);
	double cot3b=cotb*cotb*cotb;
	double cot6b=cot3b*cot3b;
	double sec2b=secb*secb;
	double expterm=((2+3*sec2b)*cot3b/24
			-((4+sec2b)*sec2b*cot6b/16
			  +((16-(1512+(3654+375*sec2b)*sec2b)*sec2b)*cot3b/5760
			    +(32+(288+(232+13*sec2b)*sec2b)*sec2b)*sec2b*cot6b/(128*nu))*
			  cot6b/nu)/nu)/nu;
	jl=sqrt(cotb*cosb)/(2*nu)*exp(-nu*beta+nu/cotb-expterm);
      }
      else if(ax>nu+1.48*l3) {
	double cosb=nu/ax;
	double sx=sqrt(ax2-nu2);
	double cotb=nu/sx;
	double secb=ax/nu;
	double beta=acos(cosb);
	double cot3b=cotb*cotb*cotb;
	double cot6b=cot3b*cot3b;
	double sec2b=secb*secb;
	double trigarg=nu/cotb-nu*beta-0.25*M_PI-
	  ((2+3*sec2b)*cot3b/24+(16-(1512+(3654+375*sec2b)*sec2b)*sec2b)*
	   cot3b*cot6b/(5760*nu2))/nu;
	double expterm=((4+sec2b)*sec2b*cot6b/16-
			(32+(288+(232+13*sec2b)*sec2b)*sec2b)*
			sec2b*cot6b*cot6b/(128*nu2))/nu2;
	jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);
      }
      else {
	double beta=ax-nu;
	double beta2=beta*beta;
	double sx=6/ax;
	double sx2=sx*sx;
	double secb=pow(sx,0.3333333333333333333333);
	double sec2b=secb*secb;
	
	jl=(CSM_GAMMA1*secb+beta*CSM_GAMMA2*sec2b
	    -(beta2/18-1.0/45.0)*beta*sx*secb*CSM_GAMMA1
	    -((beta2-1)*beta2/36+1.0/420.0)*sx*sec2b*CSM_GAMMA2
	    +(((beta2/1620-7.0/3240.0)*beta2+1.0/648.0)*beta2-1.0/8100.0)*sx2*secb*CSM_GAMMA1
	    +(((beta2/4536-1.0/810.0)*beta2+19.0/11340.0)*beta2-13.0/28350.0)*beta*sx2*sec2b*CSM_GAMMA2
	    -((((beta2/349920-1.0/29160.0)*beta2+71.0/583200.0)*beta2-121.0/874800.0)*
	      beta2+7939.0/224532000.0)*beta*sx2*sx*secb*CSM_GAMMA1)*sqrt(sx)/CSM_ROOTPI12;
      }
    }
  }
  if((x<0)&&(l%2!=0)) jl=-jl;
  
  return jl;
}

double csm_j_bessel(int l,double x)
{
  //////
  // Spherical Bessel functions j_l(x)  
#ifdef _USE_GSL_BESSEL
  gsl_sf_result res;
  int stat=gsl_sf_bessel_jl_e(l,x,&res);
  if(stat==GSL_EUNDRFLW)
    return 0;
  else {
    int_error_handle(stat,res.val,res.err);
    return res.val;
  }
#else //_USE_GSL_BESSEL
  return sph_bessel(l,x);
#endif //_USE_GSL_BESSEL
}
