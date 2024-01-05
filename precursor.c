/** 
   This is the file to generate a precursor simulation for the two-phase TBL. */

#include "embed.h"
#include "navier-stokes/centered.h"   
#include "view.h"
#include "sandbox/perfs.h"
#include "sandbox/profile6.h"   // From Antoon
#include "sandbox/maxruntime.h"   

/**
   The input parameters are 
   */

double RE_tau; // Friction Reynolds number
double ak;     // initial steepness
int MAXLEVEL;  // max level of the simulation
int LEVEL_IN;  // initial level of the simulation
int MINLEVEL;  // min level of the simulation
double uemax = 0.30*0.25;  // refinement criteria for the velocity
double femax = 1.00e-4;    // refinement criteria for the indicator function
int do_en;     // dump or not the field for ensemble
int st_wave;   // initialize or not a Stokes wave
//int RANDOM;    // random number in spectrum.h

/**
   We define these values: the wave number, fluid depth, gravity acceleration, 
   the friction velocity, the expected bulk velocity, the thermophyiscal property and wave period. */

double k_  = 1.0; // we change later 
double h_  = 1.0; // we change later 
double g_  = 1.0; // we change later
double Ustar = 1.0; // we change later
double mu1 = 1.0, mu2 = 1.0, rho1 = 1.0, rho2 = 1.0; // we change later
double T0_ = 1.0; // we change later

/**
   For the restarting step. */

int counter_max = 2;
int counter = 0;
int counter_max_ens = 20;
int counter_ens = 0;

/**
   The density and viscosity ratios are those of air and water. 
   note: the ratio is air properties over water properties */

double RHO_RATIO = 1.0;
double MU_RATIO  = 1.0;

/**
   We need to store some variables. */

face vector muv[], av[], alphav[];
scalar rhov[];

/** 
   Specify the interface shape using a third-order Stokes wave. */

double WaveProfile (double x, double y) {

  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;

}

/** 
   Add some phyiscal noise to promote faster transition. */

double fy (double y) {
  double fy_f = (sq(1.0-sq(y)));
  return fy_f;
}

double dfy (double y) {
  double dfy_f = -4.0*y*(sq(1.0-sq(y)));
  return dfy_f;
}

double gxz (double x, double z) {
  double gxz_f = z*exp(-4.*(4.*sq(x)+sq(z)));
  return gxz_f;
}

double dgxz (double x, double z) {
  double dgxz_f = exp(-4.*(4.*sq(x)+sq(z)))*(1.-8.*sq(z));
  return dgxz_f;
}

/** 
   Set profile. */

void set_profile (scalar cs, face vector fs, vertex scalar phi) {
  foreach_vertex() {
    phi[] = y-WaveProfile(x,z);
  } 
  fractions (phi,cs,fs);
  fractions_cleanup (cs,fs); // remove inconsistencies
}

#include "./spectrum.h" // Used for new input method (the spectrum info)
void import_profile (scalar cs, face vector fs, vertex scalar phi) {
  power_input();
  dkx_ = kx_[1] - kx_[0];
  dky_ = ky_[1] - ky_[0];
  fprintf (stderr, "dkx = %g, dky = %g\n", dkx_, dky_), fflush (stderr);
  foreach_vertex() {
    phi[] = y-(wave(x,z)+h_);
  }
  fractions (phi,cs,fs);
  fractions_cleanup (cs,fs); // remove inconsistencies
}


/**
   Input parameters are passed in from the command line. */

int main(int argc, char *argv[]) {

  /**
     Provide the inputs */

  maxruntime (&argc, argv);

  if (argc > 1)
    RE_tau  = atof (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    LEVEL_IN = atoi(argv[4]);
  if (argc > 5)
    MINLEVEL = atoi(argv[5]);
  if (argc > 6)
    uemax = atof(argv[6]);
  if (argc > 7)
    femax = atof(argv[7]);
  if (argc > 8)
    do_en = atoi(argv[8]);
  if (argc > 9)
    st_wave = atoi(argv[9]);
  if (argc > 10)
    RANDOM = atoi(argv[10]);
 
  if (argc < 11) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 11-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %8E seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " RE_tau = %8E\n ", RE_tau), fflush (stderr);
  fprintf(stderr, " a_0k  = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " MAXLEVEL  = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " Initial LEVEL = %d\n ", LEVEL_IN), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " ue_max = %8E\n ", uemax), fflush (stderr);
  fprintf(stderr, " fe_max = %8E\n ", femax), fflush (stderr);
  fprintf(stderr, " do_en  = %d\n ", do_en), fflush (stderr);
  fprintf(stderr, " st_wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " RANDOM = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     The domain is a cubic box centered on the origin and of length
     $L_0=2*\pi$, periodic in the x- and z-directions. */

  L0 = 2*pi;
  origin (-L0/2., 0, -L0/2.);
  init_grid (1 << LEVEL_IN);
  
  h_ = 1; // Water depth set L0/(2*pi) following Wu, Popinet & Deike JFM2022
  k_ = 4; // Four waves per box

  /**
     Set the boundary conditions */

  periodic (right);
  periodic (front);

  u.r[top] = neumann(0.);
  u.n[top] = dirichlet(0.);  
  u.t[top] = neumann(0.);
  p[top]   = neumann(0.);

  u.r[bottom] = neumann(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = neumann(0.);
  p[bottom]   = neumann(0.);

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
  u.r[embed] = dirichlet(0.);
  p[embed]   = neumann(0.);
 
  /**
     We define the proper values of the variables */

  Ustar = 0.25; // Pick a fixed value
  g_    = 9.950309E-01;
  T0_   = 3.141593;  
  rho1  = 1.225/1000.0;
  rho2  = rho1*RHO_RATIO;
  mu2   = rho2*Ustar*(L0-h_)/RE_tau;
  mu1   = mu2/(MU_RATIO);
  
  /**
     We reduce the tolerance of the solver */

  //TOLERANCE = 1e-6;
  /*
  if(st_wave != 1) {
    //TOLERANCE = 1e-6;
    TOLERANCE = HUGE;
    NITERMIN = 2;
    NITERMAX = 50;
  }
  */

  /**
     Give the address of some variables */

  a     = av;
  alpha = alphav;
  rho   = rhov;
  mu    = muv;
  
  /**
     Run! */

  run();

}

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " RE_tau = %8E\n ", RE_tau), fflush (stderr);
  fprintf(stderr, " a_0k  = %8E\n ", ak), fflush (stderr);
  fprintf(stderr, " MAXLEVEL  = %d\n ", MAXLEVEL), fflush (stderr);
  fprintf(stderr, " Initial LEVEL = %d\n ", LEVEL_IN), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MINLEVEL), fflush (stderr);
  fprintf(stderr, " ue_max  = %8E\n ", uemax), fflush (stderr);
  fprintf(stderr, " fe_max  = %8E\n ", femax), fflush (stderr);
  fprintf(stderr, " do_en  = %d\n ", do_en), fflush (stderr);
  fprintf(stderr, " st_wave = %d\n ", st_wave), fflush (stderr);
  fprintf(stderr, " RANDOM = %d\n ", RANDOM), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Create directories for better organization of the output files */

  fprintf(stderr, "Create directories (if needed) \n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p profiles");
  system(comm);
  sprintf (comm, "mkdir -p restart_bk");
  system(comm);
  sprintf (comm, "mkdir -p restart_ens");
  system(comm);
  /*
  sprintf (comm, "mkdir -p phase");
  system(comm);
  */
  
  /**
     Restart from the latest dump files or initialize a new simulation */

  if (restore ("restart.bin")) {
    fprintf(stderr, "Simulation restarts from a dumped file\n"), fflush (stderr);
    for (scalar s in {u, g}) { 
      s.prolongation = refine_embed_linear;
    }
    ///*
    cs.prolongation = refine_injection;
    adapt_wavelet ({cs,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
    cs.prolongation = fraction_refine;
    //*/
    vertex scalar phi[]; // must be here after restore!
    if(st_wave == 1) {
      fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
      set_profile (cs,fs,phi);
    }
    else {
      fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
      import_profile (cs,fs,phi);
    }
  }
  else {
    fprintf(stderr, "Initialization of all the variables\n"), fflush (stderr);
    refine (fabs(y-h_) < 0.2 && level < MAXLEVEL);
    vertex scalar phi[]; // must be here after the initial grid is set
    if(st_wave == 1) {
      fprintf(stderr, "We use a 3rd-order Stokes wave!\n"), fflush (stderr);
      set_profile (cs,fs,phi);
    }
    else {
      fprintf(stderr, "We import a profile previously generated!\n"), fflush (stderr);
      import_profile (cs,fs,phi);
    }
    double ytau = (L0-h_)/RE_tau;
    double Ubulk_ex = (mu2/(rho2*(L0-h_)))*(pow(0.5*RE_tau/0.09,(1.0/0.88))); // estimation based on Pope's relation 
    fprintf(stderr, "Ubulk_ex = %8E\n", Ubulk_ex), fflush (stderr);
    do {
      foreach() {
	if (phi[] > 0.05) {
	  u.x[] = cs[]*(log(phi[]/ytau)*Ustar/0.41);
          double x_n = 2.0*(x-0.0*L0)/L0;
          double y_n = 2.0*y/L0-1.2;
          double z_n = 2.0*(z-0.0*L0)/L0;
          u.y[] = cs[]*(-1.0*gxz(z_n,x_n)*dfy(y_n)*Ubulk_ex*1.5);
          u.z[] = cs[]*(+1.0*fy(y_n)*dgxz(z_n,x_n)*Ubulk_ex*1.5);
	}
	else {
	  u.x[] = 0;
	  u.y[] = 0;
	  u.z[] = 0;
	}
      }
    }
#if TREE
    while (0);
#else
    while (0);
#endif
  }

  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.y", linear = true, n = {1,0,0}, alpha = -L0/2.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("3D_0", "a");
    save (fp = fp);
    fflush (fp);
  }

  clear();
  view (fov = 22.5, camera = "front", ty = -0.50,
        //theta = 0.0, phi = -pi/240.0, psi = 0.0, ty = -0.5,
        width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("inter_0", "a");
    save (fp = fp);
    fflush (fp);
  }

  stats stat1 = statsf(cs); // log some stats
  double rms_cs = normf(cs).rms;

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"./log_stat_cs.out");
    FILE * log_c1 = fopen(name_1,"a");
    fprintf (log_c1, "%8E %8E %8E %8E %8E %8E\n", 
      	             t, 1.0*i, stat1.min, stat1.sum, stat1.max, rms_cs);
    fclose(log_c1);

  }

  //return 1;

}

/**
   We log the physical time, the number of time-steps and dt of the simulation. */

event log_simulation (i += 10) {
  
  int min_level = +100, max_level = -100;
  foreach(reduction(min:min_level),reduction(max:max_level)) {
    max_level = max(max_level,level);
    min_level = min(min_level,level);
  }

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %d %d\n", t, 1.0*i, dt, max_level, min_level);
    fclose(log_sim);

  }

  double u_air_mean = 0.;
  double u_air_diff = 0.;
  double u_wat_mean = 0.;
  double voll_mean  = 0.;
  double volg_mean  = 0.;
  double volg_diff  = 0.;
  foreach(reduction(+:u_air_mean) reduction(+:u_wat_mean)
          reduction(+:volg_mean ) reduction(+:voll_mean )
	  reduction(+:u_air_diff) reduction(+:volg_diff )) {
    u_air_mean += u.x[]*(0.0+cs[])*dv();
    u_wat_mean += u.x[]*(1.0-cs[])*dv();
    u_air_diff += u.x[]*dv();
    volg_mean  += (0.0+cs[])*dv();
    voll_mean  += (1.0-cs[])*dv();
    volg_diff  += dv();
  }

  /* // it is not compatible with 3D yet
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  */

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_forcing.out");
    FILE * log_sim = fopen(name_1,"a");
    /*
    fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
		      t, 1.0*i, volg_mean, u_air_mean/volg_mean,
		      voll_mean, u_wat_mean/voll_mean,
		      Fp.x, Fp.y, 
		      Fmu.x, Fmu.y);
    */
    fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %8E %8E\n", 
		      t, 1.0*i, volg_mean, u_air_mean/volg_mean,
		      voll_mean, u_wat_mean/voll_mean,
		      u_air_diff/volg_diff,volg_diff);
    fclose(log_sim);

  }

}

void profile_output (int istep) {

  /**
     Compute the vorticity. */

  scalar omega[];
  vorticity (u, omega);

  char file[99];
  sprintf (file, "./profiles/prof_%09d.out", istep);
  scalar uxuy[],uxux[],uyuy[],uzuz[];
  foreach () {
    uxuy[] = u.x[]*u.y[];
    uxux[] = u.x[]*u.x[];
    uyuy[] = u.y[]*u.y[];
    uzuz[] = u.z[]*u.z[];
  }

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = y;
  }
  profiles ({u.x, u.y, u.z, cs, omega, uxuy, uxux, uyuy, uzuz, p}, phi, rf = 1.0, fname = file, n = 1 << MAXLEVEL);

}

event out_pro (t += 2.0*T0_) {

  profile_output (i);
   
  if (pid() == 0) {
    char name_1[80];
    sprintf (name_1,"./profiles/log_pro.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E\n", t, 1.0*i);
    fclose(log_sim);
  }

}

/**
   We update the density and the dynamic viscosity to match the desire RE_tau. */

event properties (i++) {
  foreach() {
    rhov[] = cm[]*rho2;
  }
  foreach_face() {
    alphav.x[] = fm.x[]/rho2;
    muv.x[]    = fm.x[]*mu2;
  }
}

/**
   We prescribe a uniform pressure gradient to drive the flow. */

event acceleration (i++) {
  
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x) {
    av.x[] += fm.x[]*ampl*cs[];
  }

}

/** 
   Output video and field in uncompressed format 
   (we can compress later using convert <FILE>.ppm to <FILE>.mpg and open with mplayer) */

//# define POPEN(name, mode) fopen (name ".ppm", mode)

event movies (t += 5.0*T0_) {
  
  clear();
  view (fov = 27.5,  
        theta = pi/4.0, phi = pi/50.0, psi = 0.0, ty = -0.5,
	width = 1300, height = 1300, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof ("cs", color = "u.x");
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -L0/2.0);
  squares ("u.y", linear = true, n = {0,0,1}, alpha = 0.0);
  squares ("u.z", linear = true, n = {1,0,0}, alpha = 0.0);
  cells (n = {1,0,0}, alpha = -L0/2.0);
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
    fflush (fp);
  }

}

/** 
   Dump every tout_res period */

event dumpstep (t += 2.0*T0_) {

  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  char dname[100];
  sprintf (dname, "dump_%d.bin", counter);
  dump (dname);

  if (pid() == 0) {
    char comm_2[80];
    sprintf (comm_2, "ln -sf dump_%d.bin restart.bin", counter);
    system(comm_2);

    fflush(stderr);
    char name_0[80];
    sprintf (name_0, "dump_%d.bin", counter);
    FILE * fp = fopen(name_0, "r");
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, 1.0*counter, 1.0*res);
    fclose(log_sim);
  }

  if(do_en == 1) {

    if(counter_ens < counter_max_ens) {
      counter_ens++;
    }
    else {
      return 1;
    }
    
    double u_air_mean = 0.;
    double u_wat_mean = 0.;
    double voll_mean  = 0.;
    double volg_mean  = 0.;
    foreach(reduction(+:u_air_mean) reduction(+:u_wat_mean)
            reduction(+:volg_mean ) reduction(+:voll_mean )) {
      u_air_mean += u.x[]*(0.0+cs[])*dv();
      u_wat_mean += u.x[]*(1.0-cs[])*dv();
      volg_mean  += (0.0+cs[])*dv();
      voll_mean  += (1.0-cs[])*dv();
    }
    
    if (pid() == 0) {
    
      fflush(stderr);
      char name_1[80];
      sprintf (name_1,"./restart_ens/log_ens_forcing.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %02d\n", 
    		      t, 1.0*i, volg_mean, u_air_mean/volg_mean,
          	      voll_mean, u_wat_mean/voll_mean, counter_ens);
      fclose(log_sim);
    
    }
    
    sprintf (dname, "./restart_ens/dump_s%02d_ak%.2e_t%.10e.bin", counter_ens, ak, t);
    dump (dname);

  }

}

event dump_backup (t += 10.0*T0_) {
 
  char dname[100];
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = 100000.0) {

  fprintf (fout, "i = %d t = %8E\n", i, t); fflush(fout);
  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}

/** 
   Adaptive function */ 
   
#if TREE
event adapt (i++) {
  cs.prolongation = refine_injection;
  adapt_wavelet ({cs,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL, MINLEVEL);
  cs.prolongation = fraction_refine;
}
#endif
