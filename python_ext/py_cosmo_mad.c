#include <Python.h>
#include <structmember.h>
#include "cosmo_mad.h"

#ifndef PyVarObject_HEAD_INIT
#define PyVarObject_HEAD_INIT(type, size) \
  PyObject_HEAD_INIT(type) size,
#endif
#ifndef Py_TYPE
#define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif

typedef struct {
  PyObject_HEAD
  Csm_params *p;
  /* Type-specific fields go here. */
} PcsPar;

static int PcsPar_init(PcsPar *self,PyObject *args,PyObject *kwds)
{
  self->p=csm_params_new();

  return 0;
}

static void PcsPar_dealloc(PcsPar *self)
{
  csm_params_free(self->p);
  Py_TYPE(self)->tp_free((PyObject *)self);
  //  self->ob_type->tp_free((PyObject *)self);
}

static PyObject *set_verbosity(PcsPar *self,PyObject *args)
{
  int flag_verbose;
  
  if(!PyArg_ParseTuple(args,"i",&flag_verbose)) return NULL;

  csm_set_verbosity(self->p,flag_verbose);

  Py_RETURN_NONE;
}

static PyObject *hubble(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_hubble(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *omega_matter(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_omega_matter(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *cosmic_time(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_cosmic_time(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *particle_horizon(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_particle_horizon(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *radial_comoving_distance(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_radial_comoving_distance(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *curvature_comoving_distance(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_curvature_comoving_distance(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *angular_diameter_distance(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_angular_diameter_distance(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *luminosity_distance(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_luminosity_distance(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *growth_factor(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_growth_factor(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *f_growth(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_f_growth(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *growth_factor_and_growth_rate(PcsPar *self,PyObject *args)
{
  double gf_res,fg_res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  csm_growth_factor_and_growth_rate(self->p,aa,&gf_res,&fg_res);

  return Py_BuildValue("dd",gf_res,fg_res);
}

static PyObject *theta_BAO(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_theta_BAO(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *Dz_BAO(PcsPar *self,PyObject *args)
{
  double res;
  double aa;

  if(!PyArg_ParseTuple(args,"d",&aa)) return NULL;
  
  res=csm_Dz_BAO(self->p,aa);

  return Py_BuildValue("d",res);
}

static PyObject *scale_factor(PcsPar *self,PyObject *args)
{
  double res;
  double tt;

  if(!PyArg_ParseTuple(args,"d",&tt)) return NULL;
  
  res=csm_scale_factor(self->p,tt);

  return Py_BuildValue("d",res);
}

static PyObject *background_set(PcsPar *self,PyObject *args,PyObject *kwds)
{
  double om,ol,ob,w0,wa,hh,tcmb;

  if(!PyArg_ParseTuple(args,"ddddddd",&om,&ol,&ob,&w0,&wa,&hh,&tcmb)) return NULL;
  
  csm_background_set(self->p,om,ol,ob,w0,wa,hh,tcmb);

  Py_RETURN_NONE;
}

static PyObject *set_linear_pk(PcsPar *self,PyObject *args)
{
  char *fname;
  double lkmn,lkmx,dlk,nns,s8;

  if(!PyArg_ParseTuple(args,"sddddd",&fname,&lkmn,&lkmx,&dlk,&nns,&s8))
    return NULL;
  
  csm_set_linear_pk(self->p,fname,lkmn,lkmx,dlk,nns,s8);

  Py_RETURN_NONE;
}

static PyObject *set_nonlinear_pk(PcsPar *self,PyObject *args)
{
  char *fname;

  if(!PyArg_ParseTuple(args,"s",&fname)) return NULL;
  
  csm_set_nonlinear_pk(self->p,fname);

  Py_RETURN_NONE;
}

static PyObject *Pk_linear_0(PcsPar *self,PyObject *args)
{
  double res;
  double kk;

  if(!PyArg_ParseTuple(args,"d",&kk)) return NULL;
  
  res=csm_Pk_linear_0(self->p,kk);

  return Py_BuildValue("d",res);
}

static PyObject *Pk_nonlinear(PcsPar *self,PyObject *args)
{
  double res;
  double kk;

  if(!PyArg_ParseTuple(args,"d",&kk)) return NULL;
  
  res=csm_Pk_nonlinear(self->p,kk);

  return Py_BuildValue("d",res);
}

static PyObject *xi2p_L(PcsPar *self,PyObject *args)
{
  double res;
  double r,R1,R2;
  char *wf1,*wf2;

  if(!PyArg_ParseTuple(args,"dddss",&r,&R1,&R2,&wf1,&wf2)) return NULL;
  
  res=csm_xi2p_L(self->p,r,R1,R2,wf1,wf2,1);

  return Py_BuildValue("d",res);
}

static PyObject *sig0_L(PcsPar *self,PyObject *args)
{
  double res;
  double R1,R2;
  char *wf1,*wf2;

  if(!PyArg_ParseTuple(args,"ddss",&R1,&R2,&wf1,&wf2)) return NULL;
  
  res=csm_sig0_L(self->p,R1,R2,wf1,wf2);

  return Py_BuildValue("d",res);
}

static PyObject *set_Pk_params(PcsPar *self,PyObject *args)
{
  double beta,gf,bias;
  int l_max;

  if(!PyArg_ParseTuple(args,"dddi",&beta,&gf,&bias,&l_max)) return NULL;
  
  csm_set_Pk_params(self->p,beta,gf,bias,l_max);

  Py_RETURN_NONE;
}

static PyObject *Pk_full(PcsPar *self,PyObject *args)
{
  double res;
  double kk,muk;

  if(!PyArg_ParseTuple(args,"dd",&kk,&muk)) return NULL;
  
  res=csm_Pk_full(self->p,kk,muk);

  return Py_BuildValue("d",res);
}

static PyObject *Pk_multipole(PcsPar *self,PyObject *args)
{
  double res;
  double kk;
  int l;

  if(!PyArg_ParseTuple(args,"di",&kk,&l)) return NULL;
  
  res=csm_Pk_multipole(self->p,kk,l);

  return Py_BuildValue("d",res);
}

static PyObject *p_leg(PcsPar *self,PyObject *args)
{
  double res;
  double x;
  int l;

  if(!PyArg_ParseTuple(args,"id",&l,&x)) return NULL;
  
  res=csm_p_leg(l,x);

  return Py_BuildValue("d",res);
}

static PyObject *j_bessel(PcsPar *self,PyObject *args)
{
  double res;
  double x;
  int l;

  if(!PyArg_ParseTuple(args,"id",&l,&x)) return NULL;
  
  res=csm_j_bessel(l,x);

  return Py_BuildValue("d",res);
}

static PyObject *xi_multipole(PcsPar *self,PyObject *args)
{
  double res;
  double r;
  int l;

  if(!PyArg_ParseTuple(args,"di",&r,&l)) return NULL;
  
  res=csm_xi_multipole(self->p,r,l);

  return Py_BuildValue("d",res);
}

static PyObject *set_xi_multipole_splines(PcsPar *self,PyObject *args)
{
  csm_set_xi_multipole_splines(self->p);

  Py_RETURN_NONE;
}

static PyObject *unset_xi_multipole_splines(PcsPar *self,PyObject *args)
{
  csm_unset_xi_multipole_splines(self->p);

  Py_RETURN_NONE;
}

static PyObject *xi_3D(PcsPar *self,PyObject *args)
{
  double res;
  double r,mu;

  if(!PyArg_ParseTuple(args,"dd",&r,&mu)) return NULL;
  
  res=csm_xi_3D(self->p,r,mu);

  return Py_BuildValue("d",res);
}

static PyObject *xi_pi_sigma(PcsPar *self,PyObject *args)
{
  double res;
  double pi,sigma;
  int use_multi;

  if(!PyArg_ParseTuple(args,"ddi",&pi,&sigma,&use_multi)) return NULL;
  
  res=csm_xi_pi_sigma(self->p,pi,sigma,use_multi);

  return Py_BuildValue("d",res);
}

static PyObject *M2R(PcsPar *self,PyObject *args)
{
  double res;
  double mass;

  if(!PyArg_ParseTuple(args,"d",&mass)) return NULL;
  
  res=csm_M2R(self->p,mass);

  return Py_BuildValue("d",res);
}

static PyObject *R2M(PcsPar *self,PyObject *args)
{
  double res;
  double radius;

  if(!PyArg_ParseTuple(args,"d",&radius)) return NULL;
  
  res=csm_R2M(self->p,radius);

  return Py_BuildValue("d",res);
}

static PyObject *set_mf_params(PcsPar *self,PyObject *args,PyObject *kwds)
{
  double lmmn,lmmx,dlm;

  if(!PyArg_ParseTuple(args,"ddd",&lmmn,&lmmx,&dlm)) return NULL;
  
  csm_set_mf_params(self->p,lmmn,lmmx,dlm);

  Py_RETURN_NONE;
}

static PyObject *sigmaM(PcsPar *self,PyObject *args)
{
  double res,m;

  if(!PyArg_ParseTuple(args,"d",&m)) return NULL;

  res=csm_sigmaM(self->p,m);

  return Py_BuildValue("d",res);
}

static PyObject *dlsigMdlM(PcsPar *self,PyObject *args)
{
  double res,m;

  if(!PyArg_ParseTuple(args,"d",&m)) return NULL;
  
  res=csm_dlsigMdlM(self->p,m);

  return Py_BuildValue("d",res);
}

static PyObject *multiplicity_function(PcsPar *self,PyObject *args)
{
  double res,m,z;
  char *mftype;

  if(!PyArg_ParseTuple(args,"dds",&m,&z,&mftype)) return NULL;
  
  res=csm_multiplicity_function(self->p,m,z,mftype);

  return Py_BuildValue("d",res);
}

static PyObject *mass_function_logarithmic(PcsPar *self,PyObject *args)
{
  double res,m,z;
  char *mftype;

  if(!PyArg_ParseTuple(args,"dds",&m,&z,&mftype)) return NULL;
  
  res=csm_mass_function_logarithmic(self->p,m,z,mftype);

  return Py_BuildValue("d",res);
}

/*
static PyObject *cfrac(PcsPar *self,PyObject *args)
{
  double res;
  double nu;
  char *mfmodel;

  if(!PyArg_ParseTuple(args,"ds",&nu,&mfmodel)) return NULL;
  
  res=csm_cfrac(nu,mfmodel);

  return Py_BuildValue("d",res);
}
*/

static PyMethodDef PcsParMethods[] = {
  {"set_verbosity",(PyCFunction)set_verbosity,METH_VARARGS,
   "Control verbosity"},
  {"hubble",(PyCFunction)hubble,METH_VARARGS,
   "Hubble parameter"},
  {"omega_matter",(PyCFunction)omega_matter,METH_VARARGS,
   "Matter parameter"},
  {"cosmic_time",(PyCFunction)cosmic_time,METH_VARARGS,
   "Cosmic time"},
  {"particle_horizon",(PyCFunction)particle_horizon,METH_VARARGS,
   "Particle horizon"},
  {"radial_comoving_distance",(PyCFunction)radial_comoving_distance,METH_VARARGS,
   "Radial comoving distance"},
  {"curvature_comoving_distance",(PyCFunction)curvature_comoving_distance,METH_VARARGS,
   "Curvature comoving distance"},
  {"angular_diameter_distance",(PyCFunction)angular_diameter_distance,METH_VARARGS,
   "Angular diameter distance"},
  {"luminosity_distance",(PyCFunction)luminosity_distance,METH_VARARGS,
   "Luminosity distance"},
  {"growth_factor_and_growth_rate",(PyCFunction)growth_factor_and_growth_rate,METH_VARARGS,
   "Growth factor (normalized to a at early times) and growth rate"},
  {"growth_factor",(PyCFunction)growth_factor,METH_VARARGS,
   "Growth factor (normalized to a at early times)"},
  {"f_growth",(PyCFunction)f_growth,METH_VARARGS,
   "Growth rate"},
  {"theta_BAO",(PyCFunction)theta_BAO,METH_VARARGS,
   "Angular BAO scale in degrees"},
  {"Dz_BAO",(PyCFunction)Dz_BAO,METH_VARARGS,
   "Radial BAO scale as a redshift interval"},
  {"scale_factor",(PyCFunction)scale_factor,METH_VARARGS,
   "Scale factor"},
  {"background_set",(PyCFunction)background_set,METH_VARARGS,
   "Sets background cosmological parameters"},
  {"set_linear_pk",(PyCFunction)set_linear_pk,METH_VARARGS,
   "Sets / reads linear power spectrum"},
  {"set_nonlinear_pk",(PyCFunction)set_nonlinear_pk,METH_VARARGS,
   "Sets / reads non-linear power spectrum"},
  {"Pk_linear_0",(PyCFunction)Pk_linear_0,METH_VARARGS,
   "Returns linear power spectrum at z=0"},
  {"Pk_nonlinear",(PyCFunction)Pk_nonlinear,METH_VARARGS,
   "Returns non-linear normalized to the growth factor "
   "supplied through set_Pk_params"},
  {"xi2p_L",(PyCFunction)xi2p_L,METH_VARARGS,
   "Returns linear correlation function for two density fields "
   "smoothed with different window functions"},
  {"sig0_L",(PyCFunction)sig0_L,METH_VARARGS,
   "Returns linear covariance for two density fields "
   "smoothed with different window functions"},
  {"set_Pk_params",(PyCFunction)set_Pk_params,METH_VARARGS,
   "Sets parameters for redshift-space power spectrum"},
  {"Pk_full",(PyCFunction)Pk_full,METH_VARARGS,
   "Returns redshift-space power spectrum"},
  {"Pk_multipole",(PyCFunction)Pk_multipole,METH_VARARGS,
   "Returns redshift-space power spectrum multipole"},
  {"p_leg",(PyCFunction)p_leg,METH_VARARGS,
   "Returns Legendre polynomial"},
  {"j_bessel",(PyCFunction)j_bessel,METH_VARARGS,
   "Returns spherical Bessel function"},
  {"xi_multipole",(PyCFunction)xi_multipole,METH_VARARGS,
   "Returns correlation function multipole"},
  {"set_xi_multipole_splines",(PyCFunction)set_xi_multipole_splines,METH_VARARGS,
   "Sets splines for the correlation function multipoles"},
  {"unset_xi_multipole_splines",(PyCFunction)unset_xi_multipole_splines,METH_VARARGS,
   "Unsets splines for the correlation function multipoles"},
  {"xi_3D",(PyCFunction)xi_3D,METH_VARARGS,
   "Returns 3D correlation function by summing over multipoles"},
  {"xi_pi_sigma",(PyCFunction)xi_pi_sigma,METH_VARARGS,
   "Returns 3D correlation function"},
  {"R2M",(PyCFunction)R2M,METH_VARARGS,
   "Returns mass for a given Lagrangian radius"},
  {"M2R",(PyCFunction)M2R,METH_VARARGS,
   "Returns Lagrangian radius for a given mass"},
  {"set_mf_params",(PyCFunction)set_mf_params,METH_VARARGS,
   "Sets mass function parameters"},
  {"sigmaM",(PyCFunction)sigmaM,METH_VARARGS,
   "sigma(M)"},
  {"dlsigMdlM",(PyCFunction)dlsigMdlM,METH_VARARGS,
   "-d[log(sigma(M))]/d[log(M)]"},
  {"multiplicity_function",(PyCFunction)multiplicity_function,METH_VARARGS,
   "-d[F(>M)]/d[log(M)]"},
  {"mass_function_logarithmic",(PyCFunction)mass_function_logarithmic,METH_VARARGS,
   "-d[n(M)]/d[log10(M)]"},
  {NULL,NULL,0,NULL}
};

static PyTypeObject pcs_PcsParType = {
  PyVarObject_HEAD_INIT(NULL,0)
  //  0,                          /*ob_size*/
  "py_cosmo_mad.PcsPar",      /*tp_name*/
  sizeof(PcsPar),             /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PcsPar_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,         /*tp_flags*/
  "Py_cosmo_mad params object", /* tp_doc */
  0,                          /* tp_traverse */
  0,                          /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  PcsParMethods,               /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tpi_dictoffset */
  (initproc)PcsPar_init,                /* tp_init */
};

static PyMethodDef PyCsmMethods[] = {
  {NULL}  /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef py_cosmo_mad_module = {
  PyModuleDef_HEAD_INIT,
  "py_cosmo_mad",   /* name of module */
  "Cosmology library",
  -1,
  PyCsmMethods,
  NULL,
  NULL,
  NULL,
  NULL,
};
#endif

//#ifndef PyMODINIT_FUNC/* declarations for DLL import/export */
//#define PyMODINIT_FUNC void
//#endif
#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_py_cosmo_mad(void)
#else
PyMODINIT_FUNC initpy_cosmo_mad(void)
#endif
{
  PyObject *m;
  
  pcs_PcsParType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&pcs_PcsParType)<0)
    return NULL;

#if PY_MAJOR_VERSION >= 3
  m=PyModule_Create(&py_cosmo_mad_module);
#else
  m=Py_InitModule3("py_cosmo_mad",PyCsmMethods,"Cosmology library");
#endif
  if(m==NULL)
    return NULL;

  csm_unset_gsl_eh();
  
  Py_INCREF(&pcs_PcsParType);
  PyModule_AddObject(m,"PcsPar",(PyObject *)&pcs_PcsParType);

  return m;
}
