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

  if(!(par->bg_params_set))
    csm_report_error(1,"You must set background before you set pk\n");
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
    par->pk->logkarr=csm_malloc(par->pk->numk*sizeof(double));
    pkarr=csm_malloc(par->pk->numk*sizeof(double));
    
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
    par->pk->logkarr=csm_malloc(par->pk->numk*sizeof(double));
    pkarr=csm_malloc(par->pk->numk*sizeof(double));
    
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
    par->pk->logkarr=csm_malloc(par->pk->numk*sizeof(double));
    pkarr=csm_malloc(par->pk->numk*sizeof(double));
    
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
    fpk=csm_fopen(fname,"r");
    par->pk->numk=csm_linecount(fpk);
    par->pk->logkarr=csm_malloc(par->pk->numk*sizeof(double));
    pkarr=csm_malloc(par->pk->numk*sizeof(double));
    rewind(fpk);
    for(ii=0;ii<par->pk->numk;ii++) {
      int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
      if(stat!=2)
	csm_report_error(1,"CosmoMad: Error reading file %s, line %d \n",fname,ii+1);
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

  if(!pkpar->pkL_spline_set)
    csm_report_error(1,"CosmoMad: Pk is not set \n");

  if(lgk<pkpar->logkmin) {
    double pk0;
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,pkpar->logkmin,
			       pkpar->intacc_pkL,&pk0);
    csm_int_error_handle(stat,pkpar->logkmin,pk0);
    pk=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkmin));
  }
  else if(lgk<pkpar->logkcamb) {
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,lgk,
			       pkpar->intacc_pkL,&pk);
    csm_int_error_handle(stat,lgk,pk);
  }
  else {
    double pkf;
    int stat=gsl_spline_eval_e(pkpar->spline_pkL,pkpar->logkcamb,
			       pkpar->intacc_pkL,&pkf);
    csm_int_error_handle(stat,pkpar->logkcamb,pkf);
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
  else
    csm_report_error(1,"CosmoMad: Wrong window function \n");
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
  else
    csm_report_error(1,"CosmoMad: Unknown window function %s \n",wf1);
  if(!strcmp(wf2,"Gauss"))
    pars.wf2=1;
  else if(!strcmp(wf2,"TopHat"))
    pars.wf2=0;
  else if(!strcmp(wf2,"SharpK"))
    pars.wf2=2;
  else
    csm_report_error(1,"CosmoMad: Unknown window function %s \n",wf2);

  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(10000);
  integrand.params=&pars;
  integrand.function=&integxiL_NO;
  
  int stat=gsl_integration_qagil(&integrand,logklim,0,relerrt,10000,w,
				 &integral1,&errintegral);
  csm_int_error_handle(stat,integral1,errintegral);
  integral1*=CSM_TWOPIPIINVLOGTEN;

  if(r>0) {
    gsl_integration_workspace *cw
      =gsl_integration_workspace_alloc(1000);
    gsl_integration_qawo_table *wf
      =gsl_integration_qawo_table_alloc(r,0.1,GSL_INTEG_SINE,100);
    
    integrand.function=&integxiL_O;
    stat=gsl_integration_qawf(&integrand,klim,relerrt,1000,
			      w,cw,wf,&integral2,&errintegral);
    csm_int_error_handle(stat,integral2,errintegral);
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
  if(!(par->pk_params_set))
    csm_report_error(1,"Linear Pk not set\n");

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
    csm_int_error_handle(stat,sigma_v,errsigma_v);
    
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
    csm_int_error_handle(stat,sigma_v,errsigma_v);
    
    gsl_integration_workspace_free(w);
    par->pk->sigv=sqrt(sigma_v);
    if(par->flag_verbose)
      printf(" sigma_v = %.3lE Mpc/h\n\n",par->pk->sigv);

    double *pkarr,*pkarr_lin,norm;
    int ii;

    //Create arrays
    pkarr=csm_malloc(par->pk->numk*sizeof(double));
    pkarr_lin=csm_malloc(par->pk->numk*sizeof(double));
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
    fpk=csm_fopen(fnamePkHFIT,"r");
    nk=csm_linecount(fpk);
    lkarr=csm_malloc(nk*sizeof(double));
    pkarr=csm_malloc(nk*sizeof(double));
    rewind(fpk);
    for(ii=0;ii<nk;ii++) {
      int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
      if(stat!=2)
	csm_report_error(1,"CosmoMad: Error reading file %s, line %d \n",fnamePkHFIT,ii+1);
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
  if(!(par->pk_params_set))
    csm_report_error(1,"Pk is not set \n");
  if(l_max%2!=0)
    csm_report_error(1,"CosmoMad: multipole number must be even\n");

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
	csm_int_error_handle(stat,pkpar->logkminNL,pk0);
	pkS=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkminNL));
      }
      else if(lgk<pkpar->logkcambNL) {
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,lgk,
				   pkpar->intacc_pkNL,&pkS);
	csm_int_error_handle(stat,lgk,pkS);
      }
      else {
	double pkf;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkcambNL,
				   pkpar->intacc_pkNL,&pkf);
	csm_int_error_handle(stat,pkpar->logkcambNL,pkf);
	pkS=pkf*pow(10,-3*(lgk-pkpar->logkcambNL));
      }
      
      pkNL=(pkL-pkS)*exp(-kk*kk*pkpar->sigv*pkpar->sigv*pkpar->gf2)+pkS;
    }
    else {
      if(lgk<pkpar->logkminNL) {
	double pk0;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkminNL,
				   pkpar->intacc_pkNL,&pk0);
	csm_int_error_handle(stat,pkpar->logkminNL,pk0);
	pkNL=pk0*pow(10,pkpar->ns*(lgk-pkpar->logkminNL));
      }
      else if(lgk<pkpar->logkcambNL) {
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,lgk,
				   pkpar->intacc_pkNL,&pkNL);
	csm_int_error_handle(stat,lgk,pkNL);
      }
      else {
	double pkf;
	int stat=gsl_spline_eval_e(pkpar->spline_pkNL,pkpar->logkcambNL,
				   pkpar->intacc_pkNL,&pkf);
	csm_int_error_handle(stat,pkpar->logkcambNL,pkf);
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
    if((l%2!=0)||(l>pkpar->l_multi_max))
      csm_report_error(1,"CosmoMad: multipoles must be even and < %d\n",pkpar->l_multi_max);

    if(logk<pkpar->logkarr[0]) {
      double pk0;
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],
				 pkpar->logkarr[0],
				 pkpar->intacc_pkmulti[il],&pk0);
      csm_int_error_handle(stat,pkpar->logkarr[0],pk0);
      pk=pk0*pow(10,pkpar->ns*(logk-pkpar->logkarr[0]));
    }
    else if(logk<pkpar->logkarr[pkpar->numk-1]) {
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],logk,
				 pkpar->intacc_pkmulti[il],&pk);
      csm_int_error_handle(stat,logk,pk);
    }
    else {
      double pkf;
      int stat=gsl_spline_eval_e(pkpar->spline_pkmulti[il],
				 pkpar->logkarr[pkpar->numk-1],
				 pkpar->intacc_pkmulti[il],&pkf);
      csm_int_error_handle(stat,pkpar->logkarr[pkpar->numk-1],pkf);
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
    
    csm_int_error_handle(status,integral,errintegral);
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

  pkpar->intacc_pkmulti=csm_malloc(nl*sizeof(gsl_interp_accel *));
  pkpar->spline_pkmulti=csm_malloc(nl*sizeof(gsl_spline *));

  for(ii=0;ii<nl;ii++) {
    int jj;

    double *pkarr=csm_malloc(pkpar->numk*sizeof(double));
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
    if(l>par->pk->l_multi_max)
      csm_report_error(1,"CosmoMad: Multipole l=%d not calculated \n",l);
    if(l%2!=0)
      csm_report_error(1,"CosmoMad: multipoles must be even \n");
    int il=l/2;
    if(rr<CSM_RMIN_LOG) return par->xi->ximultimin[il];
    else if(rr<CSM_RMIN_LIN) {
      double result;
      int stat;
      double logr=log10(rr);
      stat=gsl_spline_eval_e(par->xi->spline_ximulti_log[il],logr,
			     par->xi->intacc_ximulti_log[il],&result);
      csm_int_error_handle(stat,CSM_RMIN_LOG,logr);
      return result;
    }
    else if(rr<CSM_RMAX_LIN) {
      double result;
      int stat;
      stat=gsl_spline_eval_e(par->xi->spline_ximulti_lin[il],rr,
			     par->xi->intacc_ximulti_lin[il],&result);
      csm_int_error_handle(stat,CSM_RMAX_LIN,rr);

      return result;
    }
    else
      return 0;
  }
  else {
    if(!(par->pk_params_set))
      csm_report_error(1,"CosmoMad: must set Pk params before calculating xi!\n");
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
    else
      csm_report_error(1,"CosmoMad: multipoles must be even \n");

    pars.rr=rr;
    pars.l=l;
    pars.pkpar=par->pk;

    gsl_integration_workspace *w
      =gsl_integration_workspace_alloc(1000);
    integrand.params=&pars;
    integrand.function=integxi_l_NO;

    status=gsl_integration_qagil(&integrand,logklim,0,relerrt,1000,
    				 w,&integral1,&errintegral);
    csm_int_error_handle(status,integral1,errintegral);
    integral1*=CSM_TWOPIPIINVLOGTEN;

    if(rr>0) {
      gsl_integration_workspace *cw
	=gsl_integration_workspace_alloc(1000);
      gsl_integration_qawo_table *wf
	=gsl_integration_qawo_table_alloc(rr,0.1,GSL_INTEG_SINE,100);
      
      integrand.function=&integxi_l_O;
      status=gsl_integration_qawf(&integrand,klim,relerrt,1000,
				  w,cw,wf,&integral2,&errintegral);
      csm_int_error_handle(status,integral2,errintegral);
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

  par->xi->intacc_ximulti_lin=csm_malloc(nl*sizeof(gsl_interp_accel *));
  par->xi->intacc_ximulti_log=csm_malloc(nl*sizeof(gsl_interp_accel *));
  par->xi->spline_ximulti_lin=csm_malloc(nl*sizeof(gsl_spline *));
  par->xi->spline_ximulti_log=csm_malloc(nl*sizeof(gsl_spline *));
  par->xi->ximultimin=csm_malloc(nl*sizeof(double));

  for(ii=0;ii<nl;ii++) {
    int jj;
    
    par->xi->ximultimin[ii]=csm_xi_multipole(par,CSM_RMIN_LOG,2*ii);

    //Set array for xi_lin
    double *xi_arr_lin=csm_malloc(nr_lin*sizeof(double));
    double *rarr_lin=csm_malloc(nr_lin*sizeof(double));
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
    double *xi_arr_log=csm_malloc(nr_log*sizeof(double));
    double *lograrr_log=csm_malloc(nr_log*sizeof(double));
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
  csm_int_error_handle(stat,integral,errintegral);

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
  csm_int_error_handle(stat,integral,errintegral);

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
    csm_int_error_handle(stat,integral1,errintegral);
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
      csm_int_error_handle(stat,res_sin,errintegral);
      gsl_integration_qawo_table_free(wf);

      wf=gsl_integration_qawo_table_alloc(sigma,0.1,
					  GSL_INTEG_COSINE,100);
      stat=gsl_integration_qawf(&integrand,klim,relerrt,1000,
				  w,cw,wf,&res_cos,&errintegral);
      csm_int_error_handle(stat,res_cos,errintegral);
      gsl_integration_qawo_table_free(wf);
      
      integral2=CSM_TWOPIPIINV*(res_cos+res_sin);
      gsl_integration_workspace_free(cw);
    }
    gsl_integration_workspace_free(w);

    return integral1+integral2;
  }
}
