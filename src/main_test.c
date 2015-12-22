#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmo_mad.h"

int main(int argc,char **argv)
{
  double zz,aa;
  Csm_params *pars=csm_params_new();

  csm_set_verbosity(1);
  csm_background_set(pars,0.25,0.75,0.044,-1,0,0.7,2.75);
  
  printf("Enter redshift: ");
  if(scanf("%lf",&zz)!=1) {
    fprintf(stderr,"error: Didn't get that!\n");
    exit(1);
  }
  aa=1/(1+zz);
  printf(" \n");
  printf(" - Scale factor a = %.3lf \n",aa);
  printf(" - Comoving particle horizon r_p(z) = %.3lf Mpc/h \n",
	 csm_particle_horizon(pars,aa));
  printf(" - Radial comoving distance chi(z) = %.3lf Mpc/h \n",
	 csm_radial_comoving_distance(pars,aa));
  printf(" - Cosmic time t_BB(z) = %.3lf Gyr/h \n",
	 csm_cosmic_time(pars,aa));
  printf(" - Growth factor D(z) = %.3lf \n",
	 csm_growth_factor(pars,aa)/csm_growth_factor(pars,1));
  printf(" - Angular BAO theta_BAO(z) = %.3lf deg\n",
	 csm_theta_BAO(pars,aa));
  printf(" - Radial BAO Dz_BAO(z) = %.3lf \n",
	 csm_Dz_BAO(pars,aa));
  csm_params_free(pars);

  printf(" \n Test done! \n");
  
  return 0;
}
