///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CosmoMAD.                                    //
//                                                                   //
// CosmoMAD is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published //
// by the Free Software Foundation, either version 3 of the License, //
// or (at your option) any later version.                            //
//                                                                   //
// CosmoMAD is distributed in the hope that it will be useful, but   //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CosmoMAD.  If not, see <http://www.gnu.org/licenses/>. //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/**
 * @file cosmo_mad.h
 * @author David Alonso
 * @date 13 Oct 2013
 * @brief Header file for the CosmoMAD library
 *
 * This header file contains the definitions of all
 * the available CosmoMAD functions and macros.
 */

#ifndef _COSMO_MAD_H
#define _COSMO_MAD_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_odeiv2.h>

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

#define CSM_FOURPITHIRD 4.18879020478639 //4*pi/3
#define CSM_LOGTEN 2.302585092994 //ln(10)
#define CSM_TWOPIPIINVLOGTEN  0.1166503235296796 //ln(10)/(2*pi^2)
#define CSM_TWOPIPIINV  0.05066059182116889 //1/(2*pi^2)
#define CSM_RTOD 57.29577951308232 //180/pi
#define CSM_DTOR 0.017453292519943295 //pi/180
#define CSM_HGYR 9.777751486751187 //H0^-1 in Gyrs/h
#define CSM_HMPC 2997.92458 //H0^-1 in Mpc/h

//Mass function params
#define CSM_DELTA_C 1.686470199841 //Collapse threshold
#define CSM_RHOCRIT 2.77536627E11
#define CSM_ST_A 0.322
#define CSM_ST_a 0.707
#define CSM_ST_p 0.3
#define CSM_JP_a 1.529
#define CSM_JP_b 0.704
#define CSM_JP_c 0.412
#define CSM_JEN_B 0.315
#define CSM_JEN_b 0.61
#define CSM_JEN_q 3.8
#define CSM_WAR_A 0.7234
#define CSM_WAR_a 1.625
#define CSM_WAR_b 0.2538
#define CSM_WAR_c 1.1982
#define CSM_TINKER_A_500 0.215
#define CSM_TINKER_Aexp_500 -0.14
#define CSM_TINKER_a_500 1.585
#define CSM_TINKER_aexp_500 -0.06
#define CSM_TINKER_b_500 1.960
#define CSM_TINKER_bexp_500 -0.128
#define CSM_TINKER_c_500 1.395
#define CSM_TINKER_A_200 0.186
#define CSM_TINKER_Aexp_200 -0.14
#define CSM_TINKER_a_200 1.47
#define CSM_TINKER_aexp_200 -0.06
#define CSM_TINKER_b_200 2.57
#define CSM_TINKER_bexp_200 -0.01067562865
#define CSM_TINKER_c_200 1.19
#define CSM_TINKER10_ALPHA_200 0.368
#define CSM_TINKER10_BETA_200 0.589
#define CSM_TINKER10_GAMMA_200 0.864
#define CSM_TINKER10_PHI_200 -0.729
#define CSM_TINKER10_ETA_200 -0.243
#define CSM_TINKER10_ALPHA_500 0.387
#define CSM_TINKER10_BETA_500 0.5435
#define CSM_TINKER10_GAMMA_500 0.998
#define CSM_TINKER10_PHI_500 -0.98
#define CSM_TINKER10_ETA_500 -0.267
#define CSM_TINKER10_BETAexp 0.2
#define CSM_TINKER10_GAMMAexp -0.01
#define CSM_TINKER10_PHIexp -0.08
#define CSM_TINKER10_ETAexp 0.27
#define CSM_WATSON_A 0.282
#define CSM_WATSON_ALPHA 2.163
#define CSM_WATSON_BETA 1.406
#define CSM_WATSON_GAMMA 1.210

#define CSM_MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define CSM_MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CSM_ABS(a)   (((a) < 0) ? -(a) : (a))
#define CSM_CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
// min(max(x,low),high)


//csm_background
typedef struct {
  double rsound;
  double zeq;
  double keq;
  double zdrag;
  double kSilk;
  double rsound_approx;
  double th2p7;
  double alphac;
  double alphab;
  double betac;
  double betab;
  double bnode;
} Csm_EH_params;

typedef struct {
  //Tuneable parameters
  double OM;
  double OB;
  double OL;
  double h;
  double w0;
  double wa;
  double TCMB;

  // Derived parameters
  double OK;
  int ksign;
  int normalDE;
  int constantw;
  double bang_time;
  double phorizon;
  double a_equality;
  double growth0;
  Csm_EH_params *eh;

  // a(t) spline
  int at_spline_set;
  gsl_interp_accel *intacc_at;
  gsl_spline *spline_at;
} Csm_bg_params;
void csm_bg_params_free(Csm_bg_params *par);
Csm_bg_params *csm_bg_params_new(void);

typedef struct {
  //Parameters
  double s8;
  double ns;
  double fg;
  double bias;
  double sigv;
  double gf2;
  double norm;
  int l_multi_max;
  int use_RPT;
  int use_RPT_ss;
  int pkNL_set;

  //Pk arrays
  int numk;
  double logkcamb;
  double logkmin;
  double logkcambNL;
  double logkminNL;
  int karr_set;
  double *logkarr;
  int pkL_spline_set;
  gsl_interp_accel *intacc_pkL;
  gsl_spline *spline_pkL;
  int pkNL_spline_set;
  gsl_interp_accel *intacc_pkNL;
  gsl_spline *spline_pkNL;
  int pkmulti_spline_set;
  gsl_interp_accel **intacc_pkmulti;
  gsl_spline **spline_pkmulti;
} Csm_pk_params;
void csm_pk_params_free(Csm_pk_params *par);
Csm_pk_params *csm_pk_params_new(void);

typedef struct {
  gsl_interp_accel **intacc_ximulti_log;
  gsl_spline **spline_ximulti_log;
  gsl_interp_accel **intacc_ximulti_lin;
  gsl_spline **spline_ximulti_lin;
  double *ximultimin;
} Csm_xi_params;
void csm_xi_params_free(Csm_xi_params *par,int nl);
Csm_xi_params *csm_xi_params_new(void);

typedef struct {
  //Mass arrays
  double logmmin;
  double logmmax;
  int sM_spline_set;
  gsl_interp_accel *intacc_sM;
  gsl_spline *spline_sM;
  int dsM_spline_set;
  gsl_interp_accel *intacc_dsM;
  gsl_spline *spline_dsM;
} Csm_mf_params;
void csm_mf_params_free(Csm_mf_params *par);
Csm_mf_params *csm_mf_params_new(void);

/**
 * @brief Basic parameter structure in CosmoMAD.
 *
 * A Csm_params structure contains all the information
 * loaded by the user regarding a given cosmological model.
 */
typedef struct {
  /* Background parameters */
  int flag_verbose;
  int bg_params_set;
  Csm_bg_params *bg;
  int pk_params_set;
  Csm_pk_params *pk;
  int xi_spline_set;
  Csm_xi_params *xi;
  int mf_params_set;
  Csm_mf_params *mf;
  gsl_error_handler_t *gsl_error_handler_old;
} Csm_params;


/*****************/
/*   Utilities   */
/*****************/
int csm_linecount(FILE *f);

void csm_report_error(int level,char *fmt,...);

void *csm_malloc(size_t size);

void *csm_calloc(size_t nmemb,size_t size);

FILE *csm_fopen(const char *path,const char *mode);

void csm_int_error_handle(int status,double result,double error);

/**
 * @brief Unset GSL error handler
 *
 * Disables the default GSL error handler. The user will then be
 * notified if any errors or warnings are found regarding the
 * GSL functions (e.g.: some integral was not able to reach the
 * required precision). These warnings will often be unimportant,
 * and the code will not exit.
 */
void csm_unset_gsl_eh(Csm_params *par);

/**
 * @brief Sets verbosity level
 *
 * Determines the amount of information to be output. Only two
 * values for @p verb are supported: 0 (nothing) and 1 (everything).
 * The default verbosity level is 1 (all messages are shown).
 */
void csm_set_verbosity(Csm_params *par,int verb);

/**
 * @brief Csm_params destructor.
 *
 * Frees up all memory associated with @p par.
 */
void csm_params_free(Csm_params *par);

/**
 * @brief Csm_params creator.
 *
 * Returns an initialized Csm_params structure, without any
 * associated cosmological information (see #csm_background_set).
 */
Csm_params *csm_params_new(void);

/**
 * @brief Legendre polynomial
 *
 * Returns the Legendre polynomial of order @p l at @p x.
 */
double csm_p_leg(int l,double x);

/**
 * @brief Spherical Bessel function
 *
 * Returns the spherical Bessel function of order @p l at @p x.
 */
double csm_j_bessel(int l,double x);


/****************************/
/*   Background cosmology   */
/****************************/
/**
 * @brief Matter parameter
 *
 * Returns the matter parameter \f$\Omega_M(a)\f$ at 
 * \f$a = {\tt aa}\f$.
 */
double csm_omega_matter(Csm_params *par,double aa);

/**
 * @brief Hubble parameter
 *
 * Returns the Hubble parameter \f$H(a)\f$ at 
 * \f$a = {\tt aa}\f$ in \f$h/\f$Mpc.
 */
double csm_hubble(Csm_params *par,double aa);

/**
 * @brief Cosmic time
 *
 * Returns the cosmic time \f$t(a)\f$ at \f$a={\tt aa}\f$
 * by calculating the integral
 * \f[
 *    t(a) = \int_0^a \frac{da'}{a'\,H(a')}.
 * \f]
 */
double csm_cosmic_time(Csm_params *par,double aa);

/**
 * @brief Particle horizon
 *
 * Returns the particle horizon \f$\chi_p(a)\f$ at
 * \f$a={\tt aa}\f$ by calculating the integral 
 * \f[
 *    \chi_p(a) \int_0^a\frac{da'}{a'^2\,H(a')}.
 * \f]
 */
double csm_particle_horizon(Csm_params *par,double aa);

/**
 * @brief Radial comoving distance
 *
 * Returns the radial comoving distance \f$\chi(a)\f$ 
 * at \f$a={\tt aa}\f$ as
 * \f$\chi(a)=\chi_p(1)-\chi_p(a)\f$, where \f$\chi_p\f$ is
 * the particle horizon calculated with
 * #csm_particle_horizon.
 */
double csm_radial_comoving_distance(Csm_params *par,double aa);

/**
 * @brief Curvature comoving distance
 *
 * Returns the curvature comoving distance \f$r(a)\f$ 
 * at \f$a={\tt aa}\f$ as
 * \f[
 *   r(a) =  \frac{1}{H_0\sqrt{|\Omega_k|}}{\rm sinn}\left(H_0\sqrt{|\Omega_k|}\chi(a)\right),
 * \f]
 * where \f$\chi(a)\f$ is the radial comoving distance calculated
 * with #csm_radial_comoving_distance.
 */
double csm_curvature_comoving_distance(Csm_params *par,double aa);

/**
 * @brief Angular diameter distance
 *
 * Returns the angular diameter distance \f$d_A(a)\f$ 
 * at \f$a={\tt aa}\f$ as \f$d_A(a) = a\,r(a)\f$, where
 * \f$r(a)\f$ is the curvature comoving distance calculated
 * with #csm_curvature_comoving_distance.
 */
double csm_angular_diameter_distance(Csm_params *par,double aa);

/**
 * @brief Luminosity distance
 *
 * Returns the luminosity distance \f$d_L(a)\f$ 
 * at \f$a={\tt aa}\f$ as \f$d_L(a) = r(a)/a\f$, where
 * \f$r(a)\f$ is the curvature comoving distance calculated
 * with #csm_curvature_comoving_distance.
 */
double csm_luminosity_distance(Csm_params *par,double aa);

/**
 * @brief Growth factor AND growth rate
 *
 * Returns the linear growth factor \f$D(a)\f$ and
 * growth rate \f$f(a)\f$ simultaneously at
 * \f$a={\tt aa}\f$ in the variables \f${\tt gf}\f$
 * and \f${\tt fg}\f$. If both quantities are needed
 * at the same time, calling this function once is
 * more efficient than calling #csm_growth_factor
 * and #csm_f_growth separately, since both are simultaneously
 * estimated when solving the evolution equation 
 * for matter perturbations. 
 */
void csm_growth_factor_and_growth_rate(Csm_params *par,double aa,double *gf,double *fg);

/**
 * @brief Growth factor
 *
 * Returns the linear growth factor \f$D(a)\f$ at
 * \f$a={\tt aa}\f$ (normalized to \f$D(a\ll1)\simeq a\f$)
 * by solving the equation for matter perturbations:
 * \f[
 *   \frac{d}{da}\left(a^3\,H(a)\,\frac{dD}{da}\right)(a)=
 *          \frac{3}{2}H(a)\,a\,\Omega_M(a)D.
 * \f]
 */
double csm_growth_factor(Csm_params *par,double aa);

/**
 * @brief Growth rate
 *
 * Returns the linear growth rate \f$f(a)\equiv 
 * d\,\log D(a)/d\,\log a\f$ at \f$a={\tt aa}\f$.
 */
double csm_f_growth(Csm_params *par,double aa);

/**
 * @brief Scale factor
 *
 * Returns the value of the scale factor \f$a(t)\f$
 * at cosmic time \f$t={\tt t}\f$. The first time this
 * function is called, the relation \f$t(a)\f$ is calculated
 * for several values of \f$a\in[0,1]\f$ using #csm_cosmic_time,
 * and a spline is used to invert this relation in all subsequent
 * calls.
 */
double csm_scale_factor(Csm_params *par,double t);

/**
 * @brief Set background cosmology
 *
 * Sets background cosmology for @p par: \f$\Omega_M={\tt OM}\f$,
 * \f$\Omega_{\Lambda}={\tt OL}\f$, \f$\Omega_b={\tt OB}\f$,
 * \f$w_0={\tt w0}\f$, \f$w_a={\tt wa}\f$, \f$h={\tt hh}\f$
 * and \f$T_{\rm CMB}={\tt T\_CMB}\f$. This function must be
 * called before calculating any \f$a\f$-dependent quantity.
 */
void csm_background_set(Csm_params *par,
			double OmegaM,double OmegaL,double OmegaB,
			double ww,double wwa,double hh,double T_CMB);

/**
 * @brief Angular BAO scale
 * 
 * Returns the angular scale of the BAO (in deg) at
 * \f$a = {\tt aa}\f$ as 
 * \f[
 *   \theta_{\rm BAO}(a) = \frac{r_s}{r(a)},
 * \f]
 * where \f$r_s\f$ is the sound horizon scale estimated
 * through the Eiseinstein \& Hu fitting formula and
 * \f$r(a)\f$ is the curvature comoving distance calculated
 * with #csm_curvature_comoving_distance.
 */
double csm_theta_BAO(Csm_params *par,double aa);

/**
 * @brief Radial BAO scale
 * 
 * Returns the radial scale of the BAO as a redshift separation
 * \f$a = {\tt aa}\f$ as \f$\Delta z_{\rm BAO}(a) = H(a)\,r_s\f$,
 * where \f$r_s\f$ is the sound horizon scale estimated
 * through the Eiseinstein \& Hu fitting formula and
 * \f$H(a)\f$ is the Hubble paramteter calculated with #csm_hubble.
 */
double csm_Dz_BAO(Csm_params *par,double aa);


/****************************/
/*      Power Spectrum      */
/****************************/
/**
 * @brief Set linear power spectrum
 *
 * Initializes the linear power spectrum at \f$a=1\f$. This
 * function does different things depending on the value of
 * @p fname:
 *
 * - If @p fname is "BBKS", the power spectrum will be calculated
 *   using the BBKS transfer function in the interval
 *   \f${\tt lkmn} < \log_10(k) < {\tt lkmx}\f$ in intervals of
 *   \f$\Delta\log_10(k)={\tt dlk}\f$.
 *
 * - If @p fname is "EH", the power spectrum will be calculated
 *   using the Eisenstein \& Hu transfer function.
 *
 * - If @p fname is "EH_smooth" the power spectrum will be calculated
 *   from the Eisenstein \& Hu transfer function <b> without
 *   acoustic oscillations </b>.
 *
 * - Finally, @p fname can be set to the path to an ASCII file
 *   containing the power spectrum. This file must be in CAMB
 *   format, i.e.: two columns \f$(k,P(k))\f$ with \f$k\f$ in units
 *   of \f$h/\f$Mpc and evenly spaced in \f$\log_10(k)\f$.
 *
 * Once the \f$P(k)\f$ is read (or calculated) it is normalized to
 * \f$\sigma_8={\tt s8}\f$. If @p s8 is a negative number, the
 * normalization of the power spectrum is preserved (note that this
 * only makes sense for power spectra read from a CAMB file, since
 * we only use the BBKS of EH transfer functions with the other
 * options). After this is done, a spline is used for fast
 * interpolation. The normalization for \f$P(k)\f$ used here is
 * such that
 * \f[
 *   \langle\delta({\bf x})\delta({\bf x}+{\bf r})\rangle = 
 *   \frac{1}{2\pi^2}\int_0^{\infty}P(k)\,\frac{\sin(kr)}{kr}k^2\,dk.
 * \f]
 */
void csm_set_linear_pk(Csm_params *par,char *fname,
		       double lkmn,double lkmx,double dlk,
		       double nns,double s8);

/**
 * @brief Linear power spectrum.
 *
 * Returns the linear power spectrum at \f$a=1\f$ and wave number
 * \f$k={\tt kk}\f$. If \f${\tt kk}\f$ lies outside the
 * interpolation limits, \f$P(k)\f$ is approximated by
 * \f$P(k)\propto k^{n_s}\f$ for small \f$k\f$ and by \f$P(k)\propto
 * k^{-3}\f$ for large \f$k\f$.
 */
double csm_Pk_linear_0(Csm_params *par,double kk);

/**
 * @brief Linear correlation function
 *
 * Let \f$\delta({\bf x},R,T)\f$ be the density contrast smoothed
 * over a scale \f$R\f$ with window function \f$T\f$. This function
 * returns the correlation function
 * \f[
 *   \langle\delta({\bf x},{\tt R1},{\tt wf1})
 *   \delta({\bf x}+{\bf r},{\tt R2},{\tt wf2})\rangle.
 * \f]
 * The possible values for \f${\tt wf1}\f$ and \f${\tt wf2}\f$
 * are "TopHat" and "Gauss":
 * \f[
 *   W_{\rm TH}(x)=3\,\frac{\sin x - x\,\cos x}{x^3},\hspace{12pt}
 *   W_{\rm G}(x)=\exp(-x^2/2).
 * \f]
 * For some values of the parameters it may be impossible for the
 * GSL integrator to obtain the required accuracy, in which case
 * the error tolerance can be scaled by the argument @p errfac. For
 * most cases the default tolerance (@p errfac = 1) is OK.
 */
double csm_xi2p_L(Csm_params *par,double r,
		  double R1,double R2,
		  char *wf1,char *wf2,double errfac);

/**
 * @brief Linear covariance
 *
 * Returns the covariance of the linear density field smoothed
 * over scales @p R1 and @p R2 (i.e.: it is equivalent to calling
 * #csm_xi2p_L with r = 0).
 */
double csm_sig0_L(Csm_params *par,
		  double R1,double R2,
		  char *wf1,char *wf2);

/**
 * @brief Set non-linear power spectrum
 *
 * Initializes the non-linear power spectrum at \f$a=1\f$.
 * Three options are available depending on the value of 
 * @p fnamePkHFIT function does different things depending on the
 * value of @p fname:
 *
 * - If @p fname is "RPT", the mildly non-linear power spectrum is
 *   approximated by including a Gaussian damping term arising in
 *   renormalized perturbation theory:
 *   \f[
 *     P(k,z)=P_L(k,z)\,e^{-k^2\sigma_v^2(z)},
 *   \f]
 *   where
 *   \f[
 *     \sigma_v^2(z) = \frac{1}{6\pi^2}\int_0^{\infty}P_L(k,z)\,dk.
 *   \f]
 *   Here \f$P_L(k,z)\f$ is the linear power spectrum (see
 *   #csm_Pk_linear_0). Thus, in this case \f$\sigma_v^2(z=0)\f$
 *   is calculated and used as above whenever #csm_Pk_nonlinear.
 *
 * - If @p fnamePkHFIT is "RPT_ss", the Gaussian damping factor
 *   described above is used only on the oscillatory part of the
 *   power spectrum:
 *   \f[
 *     P(k,z)=\left[P_L(k,z)-P^{\rm no\,BAO}_L(k,z)\right]
 *     e^{-k^2\sigma_v^2(z)}+ P^{\rm no\,BAO}_L(k,z),
 *   \f]
 *   where the no-BAO power spectrum is obtained using the
 *   Eisenstein \& Hu  fitting formula.
 *
 * - Finally, @p fnamePkHFIT may be set to the path to a file
 *   containing the \f$z=0\f$ non-linear power spectrum (e.g.:
 *   using HALOFIT). The format for this file must be the same as
 *   the one used in #csm_set_linear_pk. Note that in this case
 *   there is no way to normalize the power spectrum to the value of
 *   \f$\sigma_8\f$ used for the linear case. Therefore the user
 *   must make sure that the files for both the linear and
 *   non-linear spectra were generated using the  same normalization.
 */
void csm_set_nonlinear_pk(Csm_params *par,char *fnamePkHFIT);

/**
 * @brief Set power spectrum params
 *
 * Sets the parameters necessary to calculate the full redshift-space
 * power spectrum: \f$\beta(a)={\tt beta}\f$, \f$D(a)={\tt gf}\f$ and
 * \f$b={\tt bias}\f$ (see #csm_Pk_full).
 * \f${\tt l\_max}\f$ is the maximum multipole that will be used in the
 * calculation of the power spectrum and 3D correlation function (e.g. 4
 * for the Kaiser approximation).
 */
void csm_set_Pk_params(Csm_params *par,double beta,
		       double gf,double bias,int l_max);

/**
 * @brief Non-linear power spectrum
 *
 * Returns the non-linear power spectrum for \f$k={\tt kk}\f$.
 * If \f${\tt kk}\f$ lies outside the interpolation limits,
 * \f$P(k)\f$ is approximated by \f$P(k)\propto k^{n_s}\f$ for
 * small \f$k\f$ and by \f$P(k)\propto k^{-3}\f$ for large \f$k\f$.
 * This function returns the power spectrum normalized with the
 * growth factor supplied by #csm_set_Pk_params, but without
 * bias or RSDs.
 */
double csm_Pk_nonlinear(Csm_params *par,double kk);

/**
 * @brief Full power spectrum
 *
 * Returns the full redshift-space power spectrum for
 * \f$k={\tt kk}\f$ and \f$\mu_k={\tt muk}\f$. The current
 * implementation uses the Kaiser approximation for RSDs:
 * \f[
 *   P_s(a,k,\mu_k)=b^2(1+\beta(a)\mu_k^2)^2P_{\rm NL}(a,k),
 * \f]
 * where \f$P_{\rm NL}(a,k)\f$ is the non-linear power spectrum
 * calculated using #csm_Pk_nonlinear, which is normalized using
 * the supplied growth factor.
 */
double csm_Pk_full(Csm_params *par,double kk,double muk);

/**
 * @brief Power spectrum multipole
 *
 * Returns the l-th power spectrum multipole:
 * \f[
 *   P_l(k)=\frac{2l+1}{2}\int_{-1}^1L_l(\mu_k)\,P(k,\mu_k),
 * \f]
 * where \f$L_l(x)\f$ is the Legendre polynomial of order \f$l\f$.
 */
double csm_Pk_multipole(Csm_params *par,double kk,int l);

/**
 * @brief Correlation function multipole
 *
 * Returns the l-th multipole of the redshift-space correlation
 * function through the integral
 * \f[
 *   \xi_l(r)=\frac{i^l}{2\pi^2}\int_0^{\infty}P_l(k)\,j_l(kr),
 * \f]
 * where \f$P_l(k)\f$ is the l-th power spectrum multipole
 * (as returned by #csm_Pk_multipole). The first time this function
 * is called a spline is created for each power spectrum multipole
 * in order to accelerate the calculation of the integral above.
 */
double csm_xi_multipole(Csm_params *par,double rr,int l);

/**
 * @brief Set multipole splines
 *
 * If the correlation function must be calculated repeatedly, it
 * may be faster to calculate first the multipole once for a set of
 * values of \f$r\f$ and interpolate between these afterwards. This
 * function initializes a set of spline objects that are used
 * thereafter when calling #csm_xi_multipole
 * (directly or indirectly). Specifically, a logarithmic-spaced
 * spline is used for \f$0.1 < r\,h/{\rm Mpc} < 15\f$, and a
 * linear-spaced spline is used for \f$15 < r\,h/{\rm Mpc} < 500\f$.
 * Hence subsequent calls to this function will not calculate
 * the integral in #csm_xi_multipole, but will perform a much faster
 * interpolation. If this function is call for
 * \f$r>500\,{\rm Mpc}/h\f$, #csm_xi_multipole will return 0 and
 * for \f$r<0.1\,{\rm Mpc}/h\f$ it will return the value at 
 * \f$0.1\,{\rm Mpc}/h\f$.
 */
void csm_set_xi_multipole_splines(Csm_params *par);

/**
 * @brief Unset multipole splines
 *
 * Undoes all the operations in #csm_set_xi_multipole_splines,
 * freeing up the allocated memory. It is not necessary to call
 * this function unless the splines need to be reinitialized, since
 * it is implicitly called by #csm_params_free.
 */
void csm_unset_xi_multipole_splines(Csm_params *par);

/**
 * @brief 3D correlation function (polar coords)
 *
 * Returns the anisotropic 3D correlation function as a sum over
 * multipoles
 * \f[
 *   \xi(r,\mu) = \sum_{l=0}^{\infty}\xi_l(r)\,L_l(\mu).
 * \f]
 * Note that under the Kaiser approximation used by CosmoMAD in
 * its present version only the first three multipoles
 * (\f$l=0,\,2,\,4\f$) are used. When may calls to this function
 * are necessary it may be wise to call #csm_set_xi_multipole_splines
 * first for a better performance.
 */
double csm_xi_3D(Csm_params *par,double rr,double mu);

/**
 * @brief 3D correlation function (orthogonal coords)
 *
 * Returns the anisotropic 3D correlation function using longitudinal
 * (\f$\pi=r\,\mu\f$) and transverse (\f$\sigma=r\,\sqrt{1-\mu^2}\f$)
 * coordinates. If @p use_multipoles is set to 1, the sum over
 * multipoles described ind #csm_xi_3D is used. If set to 0 the
 * following double integral is performed:
 * \f[
 *   \xi(\pi,\sigma)=\frac{1}{2\pi^2}
 *   \int_0^{\infty}dk_{\parallel}\,\cos(k_{\parallel}\pi)
 *   \int_0^{\infty}dk_{\perp}\,k_{\perp}\,J_0(k_{\perp}\sigma)
 *   P(k_{\parallel},k_{\perp}),
 * \f]
 * where \f$J_0(x)\f$ is the 0-th order cylindrical Bessel function.
 * Note that the latter approach, although exact, will be much slower
 * than the former, unless a large number of multipoles is needed.
 */
double csm_xi_pi_sigma(Csm_params *par,double pi,double sigma,
		       int use_multipoles);


/***************************/
/*      Mass function      */
/***************************/
/**
 * @brief Radius of a comoving sphere of mass M
 *
 * Returns the comoving radius of a sphere of mass M (in
 * \f$M_{\odot}/h\f$).
 */
double csm_M2R(Csm_params *par,double mass);

/**
 * @brief Mass of a sphere of comoving radius R
 *
 * Returns the comoving mass (in \f$M_{\odot}/h\f$) inside
 * a sphere of comoving radius @p radius.
 */
double csm_R2M(Csm_params *par,double radius);

/**
 * @brief Sets arrays for mass-functions calculations
 *
 * Sets up mass arrays for mass-function estimates
 */
void csm_set_mf_params(Csm_params *par,double lmmn,double lmmx,double dlm);

/**
 * @brief sigma(M)
 * 
 * Returns variance on a top-hat sphere of radius R(M)
 */
double csm_sigmaM(Csm_params *par,double m);

/**
 * @brief d[log(sigma(M))]/d[log(M)]
 * 
 * Returns variance on a top-hat sphere of radius R(M)
 */
double csm_dlsigMdlM(Csm_params *par,double m);

/**
 * @brief -dF(>M)/d(log(M))
 *
 * Returns multiplicity function
 */
double csm_multiplicity_function(Csm_params *par,double mass,double z,char *mftype);

/**
 * @brief b(M,z)
 *
 * Halo bias
 */
double csm_halo_bias(Csm_params *par,double mass,double z,char *mftype);

/**
 * @brief -dn(M)/d(log10(M))
 *
 * Returns mass function in logarithmic bins of mass
 */
double csm_mass_function_logarithmic(Csm_params *par,double mass,double z,char *mftype);

#endif //_COSMO_MAD_
