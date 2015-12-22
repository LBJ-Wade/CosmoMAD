#include <Python.h>
#include "cosmo_mad.h"

static Csm_params *glob_pars;

static PyObject *hubble(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_hubble(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *omega_matter(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_omega_matter(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *cosmic_time(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_cosmic_time(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *particle_horizon(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_particle_horizon(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *radial_comoving_distance(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_radial_comoving_distance(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *curvature_comoving_distance(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_curvature_comoving_distance(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *angular_diameter_distance(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_angular_diameter_distance(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *luminosity_distance(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_luminosity_distance(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *growth_factor(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_growth_factor(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *f_growth(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_f_growth(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *growth_factor_and_growth_rate(PyObject *self,PyObject *args)
{
  double gf_res,fg_res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  csm_growth_factor_and_growth_rate(glob_pars,aa,&gf_res,&fg_res);

  return Py_BuildValue("dd",gf_res,fg_res);
}

static PyObject *theta_BAO(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_theta_BAO(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *Dz_BAO(PyObject *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_Dz_BAO(glob_pars,aa);

  return Py_BuildValue("d",res);
}

static PyObject *scale_factor(PyObject *self,PyObject *args)
{
  double res;
  double tt;

  if(!PyArg_ParseTuple(args,"d",&tt)) return NULL;
  
  res=csm_scale_factor(glob_pars,tt);

  return Py_BuildValue("d",res);
}

static PyObject *background_set(PyObject *self,PyObject *args)
{
  double om,ol,ob,w0,wa,hh,tcmb;

  if(!PyArg_ParseTuple(args,"ddddddd",&om,&ol,&ob,&w0,&wa,&hh,&tcmb)) return NULL;
  
  csm_background_set(glob_pars,om,ol,ob,w0,wa,hh,tcmb);

  Py_RETURN_NONE;
}

static PyObject *set_linear_pk(PyObject *self,PyObject *args)
{
  char *fname;
  double lkmn,lkmx,dlk,nns,s8;

  if(!PyArg_ParseTuple(args,"sddddd",&fname,&lkmn,&lkmx,&dlk,&nns,&s8))
    return NULL;
  
  csm_set_linear_pk(glob_pars,fname,lkmn,lkmx,dlk,nns,s8);

  Py_RETURN_NONE;
}

static PyObject *set_nonlinear_pk(PyObject *self,PyObject *args)
{
  char *fname;

  if(!PyArg_ParseTuple(args,"s",&fname)) return NULL;
  
  csm_set_nonlinear_pk(glob_pars,fname);

  Py_RETURN_NONE;
}

static PyObject *Pk_linear_0(PyObject *self,PyObject *args)
{
  double res;
  double kk;

  if(!PyArg_ParseTuple(args,"d",&kk)) return NULL;
  
  res=csm_Pk_linear_0(glob_pars,kk);

  return Py_BuildValue("d",res);
}

static PyObject *Pk_nonlinear(PyObject *self,PyObject *args)
{
  double res;
  double kk;

  if(!PyArg_ParseTuple(args,"d",&kk)) return NULL;
  
  res=csm_Pk_nonlinear(glob_pars,kk);

  return Py_BuildValue("d",res);
}

static PyObject *xi2p_L(PyObject *self,PyObject *args)
{
  double res;
  double r,R1,R2;
  char *wf1,*wf2;

  if(!PyArg_ParseTuple(args,"dddss",&r,&R1,&R2,&wf1,&wf2)) return NULL;
  
  res=csm_xi2p_L(glob_pars,r,R1,R2,wf1,wf2,1);

  return Py_BuildValue("d",res);
}

static PyObject *sig0_L(PyObject *self,PyObject *args)
{
  double res;
  double R1,R2;
  char *wf1,*wf2;

  if(!PyArg_ParseTuple(args,"ddss",&R1,&R2,&wf1,&wf2)) return NULL;
  
  res=csm_sig0_L(glob_pars,R1,R2,wf1,wf2);

  return Py_BuildValue("d",res);
}

static PyObject *set_Pk_params(PyObject *self,PyObject *args)
{
  double beta,gf,bias;
  int l_max;

  if(!PyArg_ParseTuple(args,"dddi",&beta,&gf,&bias,&l_max)) return NULL;
  
  csm_set_Pk_params(glob_pars,beta,gf,bias,l_max);

  Py_RETURN_NONE;
}

static PyObject *Pk_full(PyObject *self,PyObject *args)
{
  double res;
  double kk,muk;

  if(!PyArg_ParseTuple(args,"dd",&kk,&muk)) return NULL;
  
  res=csm_Pk_full(glob_pars,kk,muk);

  return Py_BuildValue("d",res);
}

static PyObject *Pk_multipole(PyObject *self,PyObject *args)
{
  double res;
  double kk;
  int l;

  if(!PyArg_ParseTuple(args,"di",&kk,&l)) return NULL;
  
  res=csm_Pk_multipole(glob_pars,kk,l);

  return Py_BuildValue("d",res);
}

static PyObject *p_leg(PyObject *self,PyObject *args)
{
  double res;
  double x;
  int l;

  if(!PyArg_ParseTuple(args,"id",&l,&x)) return NULL;
  
  res=csm_p_leg(l,x);

  return Py_BuildValue("d",res);
}

static PyObject *j_bessel(PyObject *self,PyObject *args)
{
  double res;
  double x;
  int l;

  if(!PyArg_ParseTuple(args,"id",&l,&x)) return NULL;
  
  res=csm_j_bessel(l,x);

  return Py_BuildValue("d",res);
}

static PyObject *xi_multipole(PyObject *self,PyObject *args)
{
  double res;
  double r;
  int l;

  if(!PyArg_ParseTuple(args,"di",&r,&l)) return NULL;
  
  res=csm_xi_multipole(glob_pars,r,l);

  return Py_BuildValue("d",res);
}

static PyObject *set_xi_multipole_splines(PyObject *self,PyObject *args)
{
  csm_set_xi_multipole_splines(glob_pars);

  Py_RETURN_NONE;
}

static PyObject *unset_xi_multipole_splines(PyObject *self,PyObject *args)
{
  csm_unset_xi_multipole_splines(glob_pars);

  Py_RETURN_NONE;
}

static PyObject *xi_3D(PyObject *self,PyObject *args)
{
  double res;
  double r,mu;

  if(!PyArg_ParseTuple(args,"dd",&r,&mu)) return NULL;
  
  res=csm_xi_3D(glob_pars,r,mu);

  return Py_BuildValue("d",res);
}

static PyObject *xi_pi_sigma(PyObject *self,PyObject *args)
{
  double res;
  double pi,sigma;
  int use_multi;

  if(!PyArg_ParseTuple(args,"ddi",&pi,&sigma,&use_multi)) return NULL;
  
  res=csm_xi_pi_sigma(glob_pars,pi,sigma,use_multi);

  return Py_BuildValue("d",res);
}

static PyObject *M2R(PyObject *self,PyObject *args)
{
  double res;
  double mass;

  if(!PyArg_ParseTuple(args,"d",&mass)) return NULL;
  
  res=csm_M2R(glob_pars,mass);

  return Py_BuildValue("d",res);
}

static PyObject *R2M(PyObject *self,PyObject *args)
{
  double res;
  double radius;

  if(!PyArg_ParseTuple(args,"d",&radius)) return NULL;
  
  res=csm_R2M(glob_pars,radius);

  return Py_BuildValue("d",res);
}

static PyObject *cfrac(PyObject *self,PyObject *args)
{
  double res;
  double nu;
  char *mfmodel;

  if(!PyArg_ParseTuple(args,"ds",&nu,&mfmodel)) return NULL;
  
  res=csm_cfrac(nu,mfmodel);

  return Py_BuildValue("d",res);
}

static PyObject *collapsed_fraction(PyObject *self,PyObject *args)
{
  double res;
  double mass;
  char *mfmodel;

  if(!PyArg_ParseTuple(args,"ds",&mass,&mfmodel)) return NULL;
  
  res=csm_collapsed_fraction(glob_pars,mass,mfmodel);

  return Py_BuildValue("d",res);
}

static PyMethodDef PyCsmMethods[] = {
  {"hubble",hubble,METH_VARARGS,
   "Hubble parameter"},
  {"omega_matter",omega_matter,METH_VARARGS,
   "Matter parameter"},
  {"cosmic_time",cosmic_time,METH_VARARGS,
   "Cosmic time"},
  {"particle_horizon",particle_horizon,METH_VARARGS,
   "Particle horizon"},
  {"radial_comoving_distance",radial_comoving_distance,METH_VARARGS,
   "Radial comoving distance"},
  {"curvature_comoving_distance",curvature_comoving_distance,METH_VARARGS,
   "Curvature comoving distance"},
  {"angular_diameter_distance",angular_diameter_distance,METH_VARARGS,
   "Angular diameter distance"},
  {"luminosity_distance",luminosity_distance,METH_VARARGS,
   "Luminosity distance"},
  {"growth_factor_and_growth_rate",growth_factor_and_growth_rate,METH_VARARGS,
   "Growth factor (normalized to a at early times) and growth rate"},
  {"growth_factor",growth_factor,METH_VARARGS,
   "Growth factor (normalized to a at early times)"},
  {"f_growth",f_growth,METH_VARARGS,
   "Growth rate"},
  {"theta_BAO",theta_BAO,METH_VARARGS,
   "Angular BAO scale in degrees"},
  {"Dz_BAO",Dz_BAO,METH_VARARGS,
   "Radial BAO scale as a redshift interval"},
  {"scale_factor",scale_factor,METH_VARARGS,
   "Scale factor"},
  {"background_set",background_set,METH_VARARGS,
   "Sets background cosmological parameters"},
  {"set_linear_pk",set_linear_pk,METH_VARARGS,
   "Sets / reads linear power spectrum"},
  {"set_nonlinear_pk",set_nonlinear_pk,METH_VARARGS,
   "Sets / reads non-linear power spectrum"},
  {"Pk_linear_0",Pk_linear_0,METH_VARARGS,
   "Returns linear power spectrum at z=0"},
  {"Pk_nonlinear",Pk_nonlinear,METH_VARARGS,
   "Returns non-linear normalized to the growth factor "
   "supplied through set_Pk_params"},
  {"xi2p_L",xi2p_L,METH_VARARGS,
   "Returns linear correlation function for two density fields "
   "smoothed with different window functions"},
  {"sig0_L",sig0_L,METH_VARARGS,
   "Returns linear covariance for two density fields "
   "smoothed with different window functions"},
  {"set_Pk_params",set_Pk_params,METH_VARARGS,
   "Sets parameters for redshift-space power spectrum"},
  {"Pk_full",Pk_full,METH_VARARGS,
   "Returns redshift-space power spectrum"},
  {"Pk_multipole",Pk_multipole,METH_VARARGS,
   "Returns redshift-space power spectrum multipole"},
  {"p_leg",p_leg,METH_VARARGS,
   "Returns Legendre polynomial"},
  {"j_bessel",j_bessel,METH_VARARGS,
   "Returns spherical Bessel function"},
  {"xi_multipole",xi_multipole,METH_VARARGS,
   "Returns correlation function multipole"},
  {"set_xi_multipole_splines",set_xi_multipole_splines,METH_VARARGS,
   "Sets splines for the correlation function multipoles"},
  {"unset_xi_multipole_splines",unset_xi_multipole_splines,METH_VARARGS,
   "Unsets splines for the correlation function multipoles"},
  {"xi_3D",xi_3D,METH_VARARGS,
   "Returns 3D correlation function by summing over multipoles"},
  {"xi_pi_sigma",xi_pi_sigma,METH_VARARGS,
   "Returns 3D correlation function"},
  {"R2M",R2M,METH_VARARGS,
   "Returns mass for a given Lagrangian radius"},
  {"M2R",M2R,METH_VARARGS,
   "Returns Lagrangian radius for a given mass"},
  {"cfrac",cfrac,METH_VARARGS,
   "Returns the collapsed fraction for universal variable nu"},
  {"collapsed_fraction",collapsed_fraction,METH_VARARGS,
   "Returns the collapsed fraction as a function of mass"},
  {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initpy_cosmo_mad(void)
{
  (void)Py_InitModule("py_cosmo_mad",PyCsmMethods);

  glob_pars=csm_params_new();
  csm_unset_gsl_eh();
  
}
