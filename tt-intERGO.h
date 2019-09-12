#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <neon/ne_session.h>
#include <neon/ne_request.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <sys/types.h>
#include <sys/ipc.h> 
#include <sys/shm.h> 

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef enum
{
    FIRST_EQUATION  = 0,
    SECOND_EQUATION = 1
}
WHICH_EQUATION;

typedef enum
{
    LOWER_ROW       = 0,
    MAIN_AREA       = 1,
    UPPER_ROW       = 2
}
WHICH_BOUNDARY;
typedef struct
{
    complex double leftOfMe;
    complex double coupling;
    complex double topLeft;
    complex double bottomLeft;
    complex double derivative_t;
    complex double RHS;
}
COEFFICIENTS;

typedef enum
{
    SELF     = 0,
    fromFILE = 1,
    fromSHM  = 2
}
illum_t;

typedef struct
{
    int ascii;
    int colour;
    int nowrite;
    int progressbar;
    int smallfiles;
    int flat;
}
FLAGSt;
typedef struct
{
    const char* fieldoutfname;
    const char* featurefname;
    const char* fieldinfname;
    const char* fieldpsi0fname;
    const char* fieldpsi1fname;
    double      delta1;
    double      delta2;
    double      beta1;
    double      beta2;
    size_t      maxmem;
}
OPTIONSt;

FLAGSt FLAGS;
OPTIONSt OPTIONS;

typedef enum
{
    colour_default       =  0,
    colour_black         =  1,
    colour_red           =  2,
    colour_green         =  3,
    colour_brown         =  4,
    colour_blue          =  5,
    colour_magenta       =  6,
    colour_cyan          =  7,
    colour_gray          =  8,
    colour_dark_gray     =  9,
    colour_light_red     = 10,
    colour_light_green   = 11,
    colour_yellow        = 12,
    colour_light_blue    = 13,
    colour_light_magenta = 14,
    colour_light_cyan    = 15,
    colour_white         = 16
}
colour_t;
typedef enum
{
    highlight_default    = colour_default,
    highlight_userinput  = colour_light_blue,
    highlight_emphasize  = colour_white,
    highlight_important  = colour_yellow,
    highlight_ok         = colour_green,
    highlight_notice     = colour_light_green,
    highlight_warning    = colour_light_magenta,
    highlight_error      = colour_light_red
}
highlight_t;

int parseopt(int argc, char* const* argv);

int allocate_memory();
int free_memory();
void resource_usage();
void time_needed();
complex double chi(const char*);
char* read_req(ne_request* req);

int simple_integration();

int initialize_optical_constants();

int calculate_U(int gridpoint_s);
int surface_height();
int interp_Hdev(int NP);

int ttsolver();
int copy_y_to_psi(int gridpoint_s, const double y[]);
int func_P(__attribute__((unused)) double x, const double y[], double f[], __attribute__((unused)) void* params);
int get_first_slice(double y[], const complex float psi0cut[], const complex float psi1cut[]);
int get_boundary_from_psi0(double y[], complex float psi0cut[]);
int get_boundary_from_psi1(double y[], complex float psi1cut[]);

int get_coefficients_for_Dpsi(COEFFICIENTS* coefficients, const complex double value0[], const complex double value1[], int gridpoint_s, int gridpoint_t, WHICH_EQUATION eqn, WHICH_BOUNDARY bnd);
int setDpsi0(complex float* Dpsi0, const complex double value0[], const complex double value1[], WHICH_BOUNDARY bnd);
int setDpsi1(complex float* Dpsi1, const complex double value0[], const complex double value1[], WHICH_BOUNDARY bnd);
int copy_amplitudes_from_y(const double y[], complex double value0[], complex double value1[]);
int copy_derivatives_to_f(double f[], const complex float Dpsi0[], const complex float Dpsi1[]);

double get_s_from_index(int gridpoint_s);
double get_t_from_index(int gridpoint_t);

int convert_incoming_field();
int convert_outgoing_field();
int boundary_conditions();
int write_file_in();
int write_file_out();
int write_To_propag();
int write_SimuFile();
int write_file_psi0();
int write_file_psi1();
int print_efficiency();
double efficiency_integrate_incoming();
double efficiency_integrate_reflected();
double efficiency_integrate_transmitted();

double ALPHA2(__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t);
double BETA2 (__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t);
double SPHERICALWAVE(__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t, __attribute((unused)) WHICH_EQUATION eqn);
int precalculate_alphabeta(int gridpoint_s);
int precalculate_local_theta();
int precalculate_sphericalwave(int gridpoint_s);

void put_to_shm();
void put_plotdata_to_shm();
void put_fieldout_to_shm();
complex float ** psi0;
complex float ** psi1;
complex float *  Dpsi0;
complex float *  Dpsi1;
complex float *  field_in;
complex float *  field_out;
int gridpoint_s, gridpoint_t;
int Prof_num;
int Hdev_col;
double* ALPHA2_arr;
double* BETA2_arr ;
double* ALPHA2_rez;
double* local_theta_arr ;
double* local_thetaprime_arr ;
double* SPHERICALWAVE_arr_FE;
double* SPHERICALWAVE_arr_SE;

double wavenumber;

complex float  U0;
complex float  U1;
complex float  Um1;
complex float U1_init;
complex float Um1_init;
complex float Shf;
complex float Shf_s;
float ** Hdev;
double ** Hdev2;
//double LambdaBragg;

//Hdev matrix size
int l1,l2,Rate_Column; 
float f_l1,f_l2;
float pix_sizeX, pix_sizeY;

double thetadiff = 0.0;

int fatal_error(const char* message);
int fatal_error_maxmem();
int check_maxmem(size_t* allocated, size_t allocating);

const char* get_exponent_utf8(int e);
int parse_exponent_utf8(int expo, char* exponent);
const char* get_utf8(double value);

int print_geometry_information();
double local_theta(int gridpoint_s);
int print_simulation_information();
int print_simulation_information_after();
int print_optical_constants_information();
int progress_meter_init();
int progress_meter(int gridpoint_s, int GRIDPOINTS_S);
int progress_meter_finish();
int print_illumination_information();
int colour(colour_t col);
int highlight(highlight_t high);

int online_analysis();

int getENV(const char* key, const char** value);
int get_config_from_env();



