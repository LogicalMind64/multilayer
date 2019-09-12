/* global constants parsed from ENV variables */
double GEOM_S1;         /* distance from source to center of bottom layer */
double GEOM_S2;         /* distance from center of bottom layer to focus  */
double ENERGY;          /* photon energy in keV                           */
double MIRRORLENGTH;    /* length of the mirror                           */
double THETA;           /* angle of incidence at center of bottom layer   */
int    NUMBERLAYERS;    /* number of bilayers                             */
double LAYERRATIO;      /* ratio of the single layers in each bilayer     */
int    GRIDPOINTS_S;    /* number of grid points along mirror surface     */
int    GRIDPOINTS_T;    /* number of grid points inside the ML structure  */
char   MATERIAL1[32];
char   MATERIAL2[32];

/* derived constants */
double DOMAINSIZE_S;
double DOMAINSIZE_T;

#define GRIDSPACING_S (DOMAINSIZE_S / GRIDPOINTS_S)
#define GRIDSPACING_T (DOMAINSIZE_T / GRIDPOINTS_T)

#define PLOTEVERY (GRIDPOINTS_S / 1000)

illum_t ILLUMINATION;
double ILLUM_CENTER;
double ILLUM_WIDTH ;
double ILLUM_INTENS;
int    ILLUM_NUMBER;
double ILLUM_FACTOR;
char*  ILLUM_FILE;
int    ILLUM_SHM;

double lambda;
double LambdaBragg;


double Coh_V;
double Coh_H;
double Coh_dist;
int Coh_step;
int DivCorr;
double subsampling;
char* fnMLtopo;
int NxCAM;
int NyCAM;
double SxCAM;
double SyCAM;
double FWHMpsf;
double DistanceMIN;
double DistanceMAX;
int DistStep;

char* fnMLsimulation;
char* fnTTfield;
char* fnTTfeature;

//double averageDelta = 1.137e-5;
double averageDelta = 1.137e-5;

#define CORRECTOR (2*sin(local_theta_arr[gridpoint_s])*(averageDelta/sin(local_theta_arr[gridpoint_s])*2/3+thetadiff)*wavenumber)
//#define CORRECTOR (2*sin(local_theta_arr[gridpoint_s])*(thetadiff+averageDelta/sin(local_theta_arr[gridpoint_s]))*wavenumber)
//#define CORRECTOR (2*sin(local_theta_arr[gridpoint_s])*(thetadiff+averageDelta/local_theta_arr[GRIDPOINTS_S/2])*wavenumber)
//#define CORRECTOR (2*sin(local_theta_arr[gridpoint_s])*thetadiff*wavenumber-2*creal(U0))

//#define CORRECTOR (2*wavenumber*(  averageDelta + thetadiff*( sin(local_theta_arr[gridpoint_s])+(s-s0)*local_thetaprime_arr[gridpoint_s]  )  ))

