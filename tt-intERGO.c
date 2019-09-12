/*
This is a script for simulation of X-ray multilayer with heigh defect
by Pierre Piault for ESRF, pierre-etienne.piault@esrf.fr

 * based on Takagi Taupin solver by
 * Markus Osterhoff, ESRF / Uni Göttingen <mosterh1@gwdg.de>, 2009-2011
 * see http://dx.doi.org/10.1364/OL.37.003705
 * and http://dx.doi.org/10.1364/OL.38.005126
 */



#include "tt-intERGO.h"
#include "paramsERGO.h"

/* put in setup
 python ./Propag/TF.py "field_out.dat" "To_propagate.dat"
*/ 

int main(int argc, char* argv[])
{
	
    parseopt(argc, argv);
    get_config_from_env();
    
    wavenumber  = 2*M_PI / 12.4 * ENERGY * 1e10; // en pm-1  (1e-12 m)

    initialize_optical_constants();
    
    printf(" gridspacing s : %f \n", GRIDSPACING_S*1e6);
	printf(" gridspacing t : %f \n", GRIDSPACING_T*1e6);

	printf(" Domaine size S : %f \n", DOMAINSIZE_S);
	printf(" Domaine size T : %f \n", DOMAINSIZE_T);
	printf(" POint s : %d \n", GRIDPOINTS_S);
	printf(" POint t : %d \n", GRIDPOINTS_T);
    /*	set surface height deviation  and allocation of memory for keep data	*/ 
	surface_height();
    print_optical_constants_information();
    print_geometry_information();
    /*  allocate_memory();	*/
	print_simulation_information();
    print_illumination_information();
    /*	setup boundary conditions  */
    boundary_conditions();
	write_SimuFile();

	for (Prof_num=0;Prof_num<l1;Prof_num++)
	{
		printf("# Profil [%d/%d] \n",Prof_num+1,l1);
		Hdev_col=0;
		/* solve system of differential equations */
		ttsolver();

		/* derive outgoing field by applying phase-factor to psi1(bottom) */
		convert_incoming_field();
		convert_outgoing_field();

		/* output data to disk and memory */
		if (FLAGS.nowrite == 0)
		{
		write_file_in();
		write_file_out();

	//	put_to_shm();

		if (FLAGS.smallfiles == 0)
		{
			write_file_psi0();
			write_file_psi1();
		}
		}
		/* if field0 or field1 given, disrespect -n and -s only */
		if (FLAGS.nowrite == 1 || FLAGS.smallfiles == 1)
		{
		if (OPTIONS.fieldpsi0fname != NULL)
			write_file_psi0();
		if (OPTIONS.fieldpsi1fname != NULL)
			write_file_psi1();
		}
		if (FLAGS.nowrite == 1)
		{
		if (OPTIONS.fieldinfname != NULL)
			write_file_in();
		if (OPTIONS.fieldoutfname != NULL)
			write_file_out();
		}
	}
	write_To_propag();

    online_analysis();

    print_simulation_information_after();
    
    free_memory();

    return 0;
}

int getENV(const char* key, const char** value)
{
    const char* s = getenv(key);
    if (s == NULL)
    {
	fprintf(stderr, "could not get ENV variable for key \"%s\"\n", key);
	exit(-1);
    }
    /*
    fprintf(stderr, "\t%s -> %s\n", key, s);
    */

    *value = s;
    return 0;
}

int get_config_from_env()
{
    const char* s = NULL;

    getENV("config_geometry_s1", &s);
    GEOM_S1 = 1e-6 * atof(s);
    getENV("config_geometry_s2", &s);
    GEOM_S2 = 1e-6 * atof(s);

    getENV("config_geometry_theta", &s);
    THETA = atof(s);

    getENV("config_physic_wavelength", &s);
    lambda = atof(s);
    ENERGY = 12.4/lambda * 1e-4;

    getENV("config_mirror_length", &s);
    MIRRORLENGTH = 1e-6 * atof(s);
    DOMAINSIZE_S = MIRRORLENGTH;

    getENV("config_layers_number", &s);
    NUMBERLAYERS = atoi(s);
    LambdaBragg = lambda / (2*sin(THETA));
    double thickness   = LambdaBragg * NUMBERLAYERS;
    DOMAINSIZE_T       = thickness * sin(THETA) * 1e-6;

    getENV("config_simulation_gridpointss", &s);
    GRIDPOINTS_S = atof(s);
    getENV("config_simulation_gridpointst", &s);
    GRIDPOINTS_T = atof(s);

    getENV("config_layers_material1", &s);
    strncpy(MATERIAL1, s, 16);
    getENV("config_layers_material2", &s);
    strncpy(MATERIAL2, s, 16);
    getENV("config_layers_ratio", &s);
    LAYERRATIO = atof(s);

    getENV("config_illumination_method", &s);
    if (strcmp(s, "self") == 0 || strcmp(s, "SELF") == 0)
    {
	ILLUMINATION = SELF;
	getENV("config_illumination_center", &s);
	ILLUM_CENTER = DOMAINSIZE_S / atof(s);
	getENV("config_illumination_width", &s);
	ILLUM_WIDTH = DOMAINSIZE_S / atof(s);
	getENV("config_illumination_intens", &s);
	ILLUM_INTENS = atof(s);
	getENV("config_illumination_number", &s);
	ILLUM_NUMBER = atoi(s);
	getENV("config_illumination_factor", &s);
	ILLUM_FACTOR = atof(s);
    }
    else
	fprintf(stderr, "config_illumination_method %s not implemented.\n", s);
	
	getENV("config_vertical_coherence", &s);
	Coh_V = atof(s);
	getENV("config_horiz_coherence", &s);
	Coh_H = atof(s);
	getENV("config_distance_coherence", &s);
	Coh_dist = atof(s);
	getENV("config_step_coherence", &s);
	Coh_step = atoi(s);
	getENV("config_divergence_correction", &s);
	DivCorr = atoi(s);
	getENV("config_interpolation_rate", &s);
	subsampling = atof(s);
	fnMLtopo = getenv("config_Fname_MLtopo");
		
	getENV("config_detector_NxCAM", &s);
	NxCAM = atoi(s);
	getENV("config_detector_NyCAM", &s);
	NyCAM = atoi(s);
	getENV("config_detector_SxCAM", &s);
	SxCAM = atof(s);
	getENV("config_detector_SyCAM", &s);
	SyCAM = atof(s);
	getENV("config_FWHMpsf", &s);
	FWHMpsf = atof(s);
	getENV("config_detector_DistanceMIN", &s);
	DistanceMIN = atof(s);
	getENV("config_detector_DistanceMAX", &s);
	DistanceMAX = atof(s);
	getENV("config_distance_step", &s);
	DistStep = atof(s);
	
	fnMLsimulation = getenv("config_SimulationFile");
	fnTTfield = getenv("config_TTfieldfile");
	fnTTfeature = getenv("config_TTfeaturefile");
	

    return 0;
}

int parseopt(int argc, char* const* argv)
{
    int c = 0;

    static struct option long_options[] =
    {
	{ "ascii",             no_argument, &FLAGS.ascii,       1    },
	{ "colour",            no_argument, &FLAGS.colour,      1    },
	{ "flat",              no_argument, &FLAGS.flat,        1    },
	{ "nowrite",           no_argument, &FLAGS.nowrite,     1    },
	{ "progressbar",       no_argument, &FLAGS.progressbar, 1    },
	{ "smallfiles",        no_argument, &FLAGS.smallfiles,  1    },

	{ "maxmem",      required_argument, NULL,               1001 },
	{ "delta1",      required_argument, NULL,               1002 },
	{ "delta2",      required_argument, NULL,               1003 },
	{ "beta1",       required_argument, NULL,               1004 },
	{ "beta2",       required_argument, NULL,               1005 },
	{ "field0",      required_argument, NULL,               1006 },
	{ "field1",      required_argument, NULL,               1007 },
	{ "fieldin",     required_argument, NULL,               1008 },
	{ "fieldout",    required_argument, NULL,               1009 },

	{ "thetadiff",   required_argument, NULL,               2001 },

	{0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    while (1)
    {
	c = getopt_long (argc, argv, "acfnps", long_options, &option_index);
	if (c == -1)
	    break;
	//printf("c: %3d (%c)\n", c, c);

	switch (c)
	{
	    case    0: continue;
	    case  '?': continue;

	    case  'a': FLAGS.ascii            = 1;
		       break;
	    case  'c': FLAGS.colour           = 1;
		       break;
	    case  'f': FLAGS.flat             = 1;
		       break;
	    case  'n': FLAGS.nowrite          = 1;
		       break;
	    case  'p': FLAGS.progressbar      = 1;
		       break;
	    case  's': FLAGS.smallfiles       = 1;
		       break;

	    case 1001: OPTIONS.maxmem         = atoi(optarg)*1048576;
		       break;
	    case 1002: OPTIONS.delta1         = atof(optarg);
		       break;
	    case 1003: OPTIONS.delta2         = atof(optarg);
		       break;
	    case 1004: OPTIONS.beta1          = atof(optarg);
		       break;
	    case 1005: OPTIONS.beta2          = atof(optarg);
		       break;
	    case 1006: OPTIONS.fieldpsi0fname = optarg;
		       break;
	    case 1007: OPTIONS.fieldpsi1fname = optarg;
		       break;
	    case 1008: OPTIONS.fieldinfname   = optarg;
		       break;
	    case 1009: OPTIONS.fieldoutfname  = optarg;
		       break;

	    case 2001: thetadiff              = atof(optarg)*1e-6;
		       break;

	    default:  fprintf(stderr, "unrecognized option: %c\n", c);
	}

    }

    /*
    printf("FLAGS.ascii:       %d\n", FLAGS.ascii);
    printf("FLAGS.colour:      %d\n", FLAGS.colour);
    printf("FLAGS.nowrite:     %d\n", FLAGS.nowrite);
    printf("FLAGS.progressbar: %d\n", FLAGS.progressbar);
    printf("FLAGS.smallfiles:  %d\n", FLAGS.smallfiles);

    printf("OPTIONS.maxmem:    %d MB\n", OPTIONS.maxmem);
    */

    //exit(-1);
    return 0;
}

int surface_height()
{	

	FILE *height_file;
 	height_file = fopen(fnMLtopo,"r");
	fscanf (height_file, "%f", &f_l1);
	fscanf (height_file, "%f", &f_l2);
	fscanf (height_file, "%f", &pix_sizeX); // en metre
	fscanf (height_file, "%f", &pix_sizeY); // en metre
	pix_sizeX=pix_sizeX*1e6; //en micron
	pix_sizeY=pix_sizeY*1e6; //en micron
	l1=f_l1;
	l2=f_l2;
	Rate_Column = GRIDPOINTS_S/l2/subsampling;
	
	//--Allocation of memory possible by the fact that the size of data know, so memory size to book known
	allocate_memory();
// be carefull in setup.sh that export config_simulation_gridpointss= is l2*Rate_column in the function int surface_height()
	int l,m;
	for (l=0;l<l1;l++) {
		for (m=0;m<l2;m++) {	
			fscanf (height_file,  "%f", &Hdev[l][m]);
			Hdev[l][m]=Hdev[l][m]*1e9; // Hdev is now in nanometer for everywhere
			//printf(" Hdev () : %f \n", Hdev[l][m]);
		}	
	}
 
	if (subsampling == 0) {	
		printf("Could not accept interpolation rate value : 0 \n");
	}
	else {
		for (l=0; l<l1; l++) {
			interp_Hdev(l);
		}
	}
	free(Hdev);
	fclose(height_file);
	return 0;
}

int interp_Hdev(int NP)
{
	double* x = NULL;
	double* y = NULL;
	x = (double*) malloc((l2+1)*sizeof(double));
	y = (double*) malloc((l2+1)*sizeof(double));
	int m, p;
	
	for (m=0; m<=l2; m++) {
			x[m] = (m-f_l2/2.0)*pix_sizeX;
			y[m] = Hdev[NP][m];
	}
		
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();

    	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, l2+1);

		gsl_spline_init(spline, x, y, l2+1); 
	
	for (p=0; p<=subsampling*l2;p++) {
		double del2 = (p-subsampling*f_l2/2.0)*pix_sizeX/subsampling;
		/*
		printf ("   p: %d \n", p);
		printf ("   NP: %d \n", NP);

		printf("x[0] : %f \n", x[0]);
		printf("x[end] : %f \n", x[l2]);
		printf(" del2 : %f \n", del2);
	//	printf(" del2 : %f \n", (subsampling*(l2-1)-subsampling*l2/2)*pix_sizeX/subsampling);
		*/
		Hdev2[NP][p] = gsl_spline_eval(spline, del2, acc);
		
	}
	fflush(stdout);   		
	gsl_spline_free (spline);
	fflush(stdout);
	gsl_interp_accel_free(acc);
	
	free(x);
    free(y);
	
	return 0;
}
		
	

/* solve the differential equation */

int ttsolver()
{
    /* GSL: choose integration scheme & initialize */
    const gsl_odeiv_step_type*  T = gsl_odeiv_step_rk2;
    gsl_odeiv_step*             s = gsl_odeiv_step_alloc(T, 4*GRIDPOINTS_T);
    gsl_odeiv_system          sys = {func_P, NULL, 4*GRIDPOINTS_T, NULL};
    int status;

    double y    [4*GRIDPOINTS_T];
    double y_err[4*GRIDPOINTS_T];
    progress_meter_init();

    /* reset ode library */
    gsl_odeiv_step_reset(s);

    /* boundary condition: we know all fields at the left boundary (normally they are zero) */
    get_first_slice(y, psi0[0], psi1[0]);

    /* pre-calculate array for local_theta */
    precalculate_local_theta();

    /* loop over gridpoint_s */
    for (gridpoint_s=1; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	double x = gridpoint_s * GRIDSPACING_S;
		
	if ( gridpoint_s >= (Hdev_col+1)*Rate_Column) {
		Hdev_col++;
	}
	calculate_U(Hdev_col);

	/* pre-calculate arrays for ALPHA2, BETA2, SPHERICALWAVE */
	precalculate_alphabeta(gridpoint_s);
	precalculate_sphericalwave(gridpoint_s);

	/* do the next step */
	status = gsl_odeiv_step_apply(s, x, GRIDSPACING_S, y, y_err, NULL, NULL, &sys);

	if (status != GSL_SUCCESS)
	    fatal_error("gsl_odeiv_step_apply() != GSL_SUCCESS!\n");

	/* get boundary values at top and bottom points */
	get_boundary_from_psi0(y, psi0[gridpoint_s]);
	get_boundary_from_psi1(y, psi1[gridpoint_s]);

	/* copy final y[] to psi1/2[gridpoint_s][] */
	copy_y_to_psi(gridpoint_s, y);

	progress_meter(gridpoint_s, GRIDPOINTS_S);
    }
    progress_meter_finish();

    return 0;
}

int copy_y_to_psi(int gridpoint_s, const double y[])
{
    int gridpoint_t;
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double real, imag;
	complex double value;

	real = y[4*gridpoint_t+0];
	imag = y[4*gridpoint_t+1];
	value = real + I*imag;

	psi0[gridpoint_s][gridpoint_t] = value;

	real = y[4*gridpoint_t+2];
	imag = y[4*gridpoint_t+3];
	value = real + I*imag;

	psi1[gridpoint_s][gridpoint_t] = value;
    }
    return 0;
}

int get_first_slice(double y[], const complex float psi0cut[], const complex float psi1cut[])
{
    /* set initial values from boundary psi0,1[0][] */
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double real, imag;
	complex double value;

	value = psi0cut[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	y[4*gridpoint_t+0] = real;
	y[4*gridpoint_t+1] = imag;

	value = psi1cut[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	y[4*gridpoint_t+2] = real;
	y[4*gridpoint_t+3] = imag;
    }
    return 0;
}

int get_boundary_from_psi0(double y[], complex float psi0cut[])
{
    /* get boundary value from psi0(bottom) */
    int gridpoint_t=0;
    {
	double real, imag;
	complex double value;

	value = psi0cut[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	y[4*gridpoint_t+0] = real;
	y[4*gridpoint_t+1] = imag;
    }
    return 0;
}
int get_boundary_from_psi1(double y[], complex float psi1cut[])
{
    /* and those from psi1(top) */
    gridpoint_t=GRIDPOINTS_T-1;
    {
	double real, imag;
	complex double value;

	value = psi1cut[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	y[4*gridpoint_t+2] = real;
	y[4*gridpoint_t+3] = imag;
    }
    return 0;
}

int func_P(__attribute__((unused)) double x, const double y[], double f[], __attribute__((unused)) void* params)
{
    /* array of complex amplitudes, to be copied from y */
    complex double value0[GRIDPOINTS_T];
    complex double value1[GRIDPOINTS_T];

    /* copy aplitudes from y */
    copy_amplitudes_from_y(y, value0, value1);

    /* set Dpsi0 for first equation */
    setDpsi0(Dpsi0, value0, value1, LOWER_ROW);
    setDpsi0(Dpsi0, value0, value1, MAIN_AREA);
    setDpsi0(Dpsi0, value0, value1, UPPER_ROW);

    /* set Dpsi1 for second equation */
    setDpsi1(Dpsi1, value0, value1, UPPER_ROW);
    setDpsi1(Dpsi1, value0, value1, MAIN_AREA);
    setDpsi1(Dpsi1, value0, value1, LOWER_ROW);

    /* copy calculated derivatives back to f[] */
    copy_derivatives_to_f(f, Dpsi0, Dpsi1);

    return GSL_SUCCESS;
}

/* /solve the differential equation */



/* auxiliary functions for solution */

int boundary_conditions()
{
    /* boundary condition: psi0(left border) */
    gridpoint_s = 0;
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	psi0[gridpoint_s][gridpoint_t] = 0.0;
    }

    /* boundary condition: psi0(lower border) */
    double amplitude = 1.0;
    gridpoint_t = 0;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	psi0[gridpoint_s][gridpoint_t] = 0.0;
	int illumcounter;
	for (illumcounter=0; illumcounter<ILLUM_NUMBER; illumcounter++)
	{
	    double x, x0, sigma2;
	    complex double value = 0.0;

	    double where = 1 + illumcounter*(ILLUM_FACTOR-1);

	    x      = gridpoint_s * GRIDSPACING_S;
	    x0     = where * ILLUM_CENTER;
	    sigma2 = pow(ILLUM_WIDTH,2);
	    value  = ILLUM_INTENS * exp(-(x-x0)*(x-x0)/2/sigma2);

	    double s = get_s_from_index(gridpoint_s);
	    double t = get_t_from_index(0);
	    double distance = t+s;

	    psi0[gridpoint_s][gridpoint_t] += value * amplitude / sqrt(distance);
	}
    }

    /* boundary condition: psi1(left border) */
    gridpoint_s = 0;
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	psi1[gridpoint_s][gridpoint_t] = 0.0;
    }

    /* boundary condition: psi1(upper border) */
    gridpoint_t = GRIDPOINTS_T-1;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	/*
	const double x = gridpoint_s * GRIDSPACING_S;
	const double x0 = ILLUM_CENTER;
	const double sigma2 = pow(ILLUM_WIDTH,2);
	double value = ILLUM_INTENS * exp(-(x-x0)*(x-x0)/2/sigma2);
	*/
	double value = 0.0;
	psi1[gridpoint_s][gridpoint_t] = value;
    }

    return 0;
}

int copy_amplitudes_from_y(const double y[], complex double value0[], complex double value1[])
{
    int gridpoint_t;
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double real, imag;
	complex double value;

	real = y[4*gridpoint_t+0];
	imag = y[4*gridpoint_t+1];
	value = real + I*imag;

	value0[gridpoint_t] = value;

	real = y[4*gridpoint_t+2];
	imag = y[4*gridpoint_t+3];
	value = real + I*imag;

	value1[gridpoint_t] = value;
    }
    return 0;
}
int copy_derivatives_to_f(double f[], const complex float Dpsi0[], const complex float Dpsi1[])
{
    int gridpoint_t;
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double real, imag;
	complex double value;

	value = Dpsi0[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	f[4*gridpoint_t+0] = real;
	f[4*gridpoint_t+1] = imag;

	value = Dpsi1[gridpoint_t];
	real = creal(value);
	imag = cimag(value);

	f[4*gridpoint_t+2] = real;
	f[4*gridpoint_t+3] = imag;
    }

    return 0;
}

int get_coefficients_for_Dpsi(COEFFICIENTS* coefficients, const complex double value0[], const complex double value1[], int gridpoint_s, int gridpoint_t, WHICH_EQUATION eqn, WHICH_BOUNDARY bnd)
{
    complex double propagator = 0.0;
    complex double coupler    = 0.0;
    complex double deviator   = 0.0;

    if (eqn == FIRST_EQUATION)
    {
	if (bnd == LOWER_ROW)
	    return -1;

	coefficients->leftOfMe       = value0[gridpoint_t  ];
	coefficients->coupling       = value1[gridpoint_t  ];
	coefficients->bottomLeft     = value0[gridpoint_t-1];

	if (bnd == MAIN_AREA)
	    coefficients->topLeft    = value0[gridpoint_t+1];
	if (bnd == UPPER_ROW)
	    coefficients->topLeft    = value0[gridpoint_t  ];
	coefficients->derivative_t   = coefficients->topLeft - coefficients->bottomLeft;

	propagator  = I*U0 - SPHERICALWAVE_arr_FE[gridpoint_t];
	deviator    = 0.0;
	coupler     = I*Um1;
    }
    if (eqn == SECOND_EQUATION)
    {
	if (bnd == UPPER_ROW)
	    return -1;

	coefficients->leftOfMe       = value1[gridpoint_t  ];
	coefficients->coupling       = value0[gridpoint_t  ];
	coefficients->topLeft        = value1[gridpoint_t+1];

	if (bnd == MAIN_AREA)
	    coefficients->bottomLeft = value1[gridpoint_t-1];
	if (bnd == LOWER_ROW)
	    coefficients->bottomLeft = value1[gridpoint_t  ];
	coefficients->derivative_t   = coefficients->bottomLeft - coefficients->topLeft;

	propagator  = I*U0 + SPHERICALWAVE_arr_SE[gridpoint_t];
	//double s    = get_s_from_index(gridpoint_s);
	//double s0   = get_s_from_index(GRIDSPACING_S/2);
	deviator    = I*CORRECTOR;
	coupler     = I*U1;
    }

    coefficients->RHS          = (propagator+deviator)*coefficients->leftOfMe + coupler*coefficients->coupling;

    return 0;
}
int setDpsi0(complex float* Dpsi0, const complex double value0[], const complex double value1[], WHICH_BOUNDARY bnd)
{
    COEFFICIENTS coefficients;

    switch (bnd)
    {
	case LOWER_ROW:
	    gridpoint_t=0;
	    Dpsi0[gridpoint_t] = 0.0;
	    return 0;

	case MAIN_AREA:
	    for (gridpoint_t=1; gridpoint_t<GRIDPOINTS_T-1; gridpoint_t++)
	    {
		get_coefficients_for_Dpsi(&coefficients, value0, value1, gridpoint_s, gridpoint_t, FIRST_EQUATION, MAIN_AREA);
		complex double derivative_t = coefficients.derivative_t;
		complex double RHS          = coefficients.RHS;

		Dpsi0[gridpoint_t] =
		    1.0*ALPHA2_rez[gridpoint_t] *
			RHS
		    - 0.5*(BETA2_arr[gridpoint_t]*ALPHA2_rez[gridpoint_t]/GRIDSPACING_T) *
			derivative_t;
	    }
	    return 0;

	case UPPER_ROW:
	    gridpoint_t = GRIDPOINTS_T-1;
	    {
		get_coefficients_for_Dpsi(&coefficients, value0, value1, gridpoint_s, gridpoint_t, FIRST_EQUATION, UPPER_ROW);
		complex double derivative_t = coefficients.derivative_t;
		complex double RHS          = coefficients.RHS;

		Dpsi0[gridpoint_t] =
		    1.0*ALPHA2_rez[gridpoint_t] *
			RHS
		    - 1.0*(BETA2_arr[gridpoint_t]*ALPHA2_rez[gridpoint_t]/GRIDSPACING_T) *
			derivative_t;
	    }
	    return 0;
    }

    return -1;
}
int setDpsi1(complex float* Dpsi1, const complex double value0[], const complex double value1[], WHICH_BOUNDARY bnd)
{
    COEFFICIENTS coefficients;

    switch (bnd)
    {
	case LOWER_ROW:
	    gridpoint_t = 0;
	    {
		get_coefficients_for_Dpsi(&coefficients, value0, value1, gridpoint_s, gridpoint_t, SECOND_EQUATION, LOWER_ROW);
		complex double derivative_t = coefficients.derivative_t;
		complex double RHS          = coefficients.RHS;

		Dpsi1[gridpoint_t] =
		    1.0*ALPHA2_rez[gridpoint_t] *
			RHS
		    - 1.0*(BETA2_arr[gridpoint_t]*ALPHA2_rez[gridpoint_t]/GRIDSPACING_T) *
			derivative_t;
	    }
	    return 0;

	case MAIN_AREA:
	    for (gridpoint_t=GRIDPOINTS_T-2; gridpoint_t>0; gridpoint_t--)
	    {
		get_coefficients_for_Dpsi(&coefficients, value0, value1, gridpoint_s, gridpoint_t, SECOND_EQUATION, MAIN_AREA);
		complex double derivative_t = coefficients.derivative_t;
		complex double RHS          = coefficients.RHS;

		Dpsi1[gridpoint_t] =
		    1.0*ALPHA2_rez[gridpoint_t] *
			RHS
		    - 0.5*(BETA2_arr[gridpoint_t]*ALPHA2_rez[gridpoint_t]/GRIDSPACING_T) *
			derivative_t;
	    }
	    return 0;

	case UPPER_ROW:
	    gridpoint_t=GRIDPOINTS_T-1;
	    Dpsi1[gridpoint_t] = 0.0;
	    return 0;
    }

    return -1;
}

/* /auxiliary functions for solution */



/* coefficients for differential equation */

int precalculate_alphabeta(int gridpoint_s)
{
    static int initialized = 0;
    if (FLAGS.flat == 1)
    {
	if (initialized == 1)
	    return 0;
	gridpoint_s = GRIDPOINTS_S/2;
    }

    static double c = -1.0;
    if (c < 0)
	c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

    const double s = get_s_from_index(gridpoint_s);
    const double c2 = c*c;
    int gridpoint_t;

    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double t = get_t_from_index(gridpoint_t);

	double s2 = s*s;
	double t2 = t*t;
	double t2s2 = 1.0/(t2-s2);

	double ALPHA2 = (c2-s2) * t2s2;
	double BETA2  = (t2-c2) * t2s2;

	ALPHA2_arr[gridpoint_t] =     ALPHA2;
	BETA2_arr [gridpoint_t] =      BETA2;
	ALPHA2_rez[gridpoint_t] = 1.0/ALPHA2;
	
	//printf(" loop t Alpha, Beta ... : %f %f %f \n", ALPHA2_arr[gridpoint_t], BETA2_arr[gridpoint_t], ALPHA2_rez[gridpoint_t]);
    }
    initialized = 1;

    return 0;
}

int precalculate_local_theta()
{
    int gridpoint_s;

    if (FLAGS.flat == 1)
    {
	for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
	{
	    local_theta_arr     [gridpoint_s] = local_theta(GRIDPOINTS_S/2);
	    local_thetaprime_arr[gridpoint_s] = 0.0;
	}
    }
    else
    {
	for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
	    local_theta_arr[gridpoint_s] = local_theta(gridpoint_s);

	for (gridpoint_s=1; gridpoint_s<GRIDPOINTS_S-1; gridpoint_s++)
	    local_thetaprime_arr[gridpoint_s] = 
		( local_theta(gridpoint_s+1) - local_theta(gridpoint_s-1) )
		/ (get_s_from_index(gridpoint_s+1)-get_s_from_index(gridpoint_s-1) );

	    local_thetaprime_arr[GRIDPOINTS_S-1] = 0.0;
	    local_thetaprime_arr[       0      ] = 0.0;
    }

    FILE* f = fopen("localtheta.dat", "w");

    const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );
    const double t = get_t_from_index(0);

    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=500)
    {
	double s = get_s_from_index(gridpoint_s);
	double A = s*t/c;
	double Arel = A - t + GEOM_S2;
	fprintf(f, "% 6.0f %8.6f\n", Arel*1e6, local_theta_arr[gridpoint_s]);
    }
    fclose(f);

    return 0;
}

int precalculate_sphericalwave(int gridpoint_s)
{
    static int initialized = 0;
    if (FLAGS.flat == 1)
    {
        if (initialized == 1)
            return 0;
        gridpoint_s = GRIDPOINTS_S/2;
    
        for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
        {   
            SPHERICALWAVE_arr_FE[gridpoint_t] = 0.0;
            SPHERICALWAVE_arr_SE[gridpoint_t] = 0.0;
        }
        initialized = 1;
    }
    
    const double s = get_s_from_index(gridpoint_s);
    int gridpoint_t;
    
    for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {   
	double t = get_t_from_index(gridpoint_t);

	double tpluss  = t+s;
	double tminuss = t-s;

	/* FIRST_EQUATION, SECOND_EQUATION */
	SPHERICALWAVE_arr_FE[gridpoint_t] = 0.5 / tpluss ;
	SPHERICALWAVE_arr_SE[gridpoint_t] = 0.5 / tminuss;
	//printf(" loop t sphericalwave FE , SE : %f %f \n", SPHERICALWAVE_arr_FE[gridpoint_t], SPHERICALWAVE_arr_SE[gridpoint_t]);
    }
    initialized = 1;
    
    return 0;
}

double ALPHA2(__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t)
{
#ifdef _CONSTANT_FACTORS
    return cos(THETA)*cos(THETA);
#else

    const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

    double s = get_s_from_index(gridpoint_s);
    double t = get_t_from_index(gridpoint_t);

    return (c*c-s*s)/(t*t-s*s);

#endif
}
double BETA2(__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t)
{
#ifdef _CONSTANT_FACTORS
    return sin(THETA)*sin(THETA);
#else

    const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

    double s = get_s_from_index(gridpoint_s);
    double t = get_t_from_index(gridpoint_t);

    return  (t*t-c*c)/(t*t-s*s);

#endif
}
double SPHERICALWAVE(__attribute((unused)) int gridpoint_s, __attribute((unused)) int gridpoint_t, __attribute((unused)) WHICH_EQUATION eqn)
{
#ifdef _CONSTANT_FACTORS
    return 0.0;
#else

    double s = get_s_from_index(gridpoint_s);
    double t = get_t_from_index(gridpoint_t);

    switch (eqn)
    {
	case FIRST_EQUATION:
	    return 0.5/(t+s);
	case SECOND_EQUATION:
	    return 0.5/(t-s);
    }
    return 0.0;

#endif
}

double get_s_from_index(int gridpoint_s)
{
    double        x = gridpoint_s * GRIDSPACING_S;
    double    s_cen = .5*(GEOM_S1-GEOM_S2);
    double    s_min = s_cen - .5*DOMAINSIZE_S;
    double        s = s_min + x;

    return s;
}
double get_t_from_index(int gridpoint_t)
{
    double        y = gridpoint_t * GRIDSPACING_T;
    double    t_bot = .5*(GEOM_S1+GEOM_S2);
    double        t = t_bot + y;

    return t;
}

/* /coefficients for differential equation */



/* input, output */

int initialize_optical_constants()
{
    const char* material1    = MATERIAL1;
    const char* material2    = MATERIAL2;

    complex double chi_mat1;
    complex double chi_mat2;

    if (OPTIONS.delta1 != 0)
    {
	complex double n = 1 - OPTIONS.delta1 + I*OPTIONS.beta1;
	complex double chi = n*n-1;

	chi_mat1 = chi;
    }
    else
	chi_mat1 = chi(material1);

    if (OPTIONS.delta2 != 0)
    {
	complex double n = 1 - OPTIONS.delta2 + I*OPTIONS.beta2;
	complex double chi = n*n-1;

	chi_mat2 = chi;
    }
    else
	chi_mat2 = chi(material2);

    complex double u_mat1   = chi_mat1 * wavenumber / 2;
    complex double u_mat2   = chi_mat2 * wavenumber / 2;

    U0 = LAYERRATIO * u_mat1 + (1-LAYERRATIO) * u_mat2;
    U1_init  = 1.0 / M_PI * (u_mat1 - u_mat2);
    //Um1_init = conj(U1_init);
    Um1_init = U1_init;

    return 0;
}

// Function reading heigh matrice to add phase shift in layered medium
int calculate_U(int Hdev_col)
{
//printf(" Prof_Num : %d ",Prof_num);

	double lambda1      = 12.4/ENERGY; //amstrom
	double LambdaBragg1 = lambda1 / (2*sin(THETA)); //amstrom this is the d-spacing
//	complex float Shf1 = -I * (complex float)2.0 * M_PI * (complex float)Hdev2[Prof_num][Hdev_col]*10; // with Hdev2[][] in nanometer so *10 to be homogenous with lambda1
//	complex float Shf2 = LambdaBragg1;
//complex float Shf2 = creal(Shf1)+1.0 + cimag(Shf1)+ I;
//Shf = I * (complex float)4.0 * M_PI * (complex float)Hdev[Prof_num][gridpoint_s] * ( (complex float)LambdaBragg * (complex float)1000000.0);

	float Shf1 = 2.0 * M_PI * Hdev2[Prof_num][Hdev_col]*10;
	float Shf2 = LambdaBragg1;
	float Shfi = Shf1/Shf2;
	
	//float Shfr = creal(Shf1)/Shf2; 
	//float Shfi = cimag(Shf1/Shf2); 
	
	Shf = cos(Shfi)+I*sin(Shfi);
	Shf_s = conj(Shf);

	U1 = U1_init * Shf;
	Um1 = Um1_init * Shf_s;
	
//	printf(" %f ::::: ", Hdev2[Prof_num][Hdev_col]);
//	printf(" 1 : %f \n ", creal((complex float)Shfi));
//	printf(" 2 : %f \n ", cimag(Shf));
//	printf(" 3 : %f \n ", creal(Shf));
//	printf(" 4 : %f \n ", cimag(U1));
//	printf(" 5 : %f \n ", creal(U1));

return 0;
}

int convert_incoming_field()
{
    int gridpoint_s = 0;
    int gridpoint_t = 0;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	complex double value = psi0[gridpoint_s][gridpoint_t];
	field_in[gridpoint_s] = value;
    }

    return 0;
}
int convert_outgoing_field()
{
    int gridpoint_s = 0;
    int gridpoint_t = 0;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	complex double value = psi1[gridpoint_s][gridpoint_t];
	field_out[gridpoint_s] = value;
    }

    return 0;
}
complex double chi(const char* material)
{
    char SERVER[64];
    unsigned int PORT;
    {
	const char* s;
	getENV("config_xocd_hostname", &s);
	strncpy(SERVER, s, 64);
	getENV("config_xocd_port"    , &s);
	PORT = atoi(s);
    }

    const char*        UAGENT = "tt-integrate/0.0x";

    char q[128];
    char* s = (char*) malloc(sizeof(char) * 128);

    snprintf(q, 128, "/%s/%fkeV", material, ENERGY);

    ne_session* sess;
    sess = ne_session_create("http", SERVER, PORT);
    ne_session_proxy(sess, "proxy.esrf.fr", 3128); /* use proxy server to connect to k-raum.org */
    ne_set_useragent(sess, UAGENT);
    ne_set_read_timeout(sess, 3);

    ne_request* req = ne_request_create(sess, "GET", q);

    int ret = ne_begin_request(req);
    const ne_status* st = ne_get_status(req);
    if (st->code != 200 || ret != NE_OK)
	fatal_error("GET of optical constants from server failed!\n");
    s = read_req(req);

    double one, delta, beta;
    sscanf(s, "%lf %lf %lf i", &one, &delta, &beta);

    ne_end_request(req);
    ne_request_destroy(req);
    ne_session_destroy(sess);
    free(s);

    complex double n = one + delta + I*beta;
    complex double chi = n*n-1;

    return chi;
}
char* read_req(ne_request* req)
{
    ne_buffer* buf = ne_buffer_create();
    ssize_t len;
    char data[16384];

    while ( (len = ne_read_response_block(req, data, sizeof data)) > 0 )
        ne_buffer_append(buf, data, len);

    if (len == 0)
        return ne_buffer_finish(buf);

    ne_buffer_destroy(buf);
    return NULL;
}
int write_file_in()
{
	const char* fname;
	if (Prof_num==0) {	

		if (OPTIONS.fieldinfname == NULL)
			fname = "field_in.dat";
		else
			fname = OPTIONS.fieldinfname;
			FILE* fout = fopen(fname, "w");
		if (fout == NULL)
			fatal_error("Could not open file for output!\n");
	}
	   
	if (OPTIONS.fieldinfname == NULL)
		fname = "field_in.dat";
	else
		fname = OPTIONS.fieldinfname;
		FILE* fout = fopen(fname, "a+");
	if (fout == NULL)
		fatal_error("Could not open file for output!\n");
	for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=10)
	{
		double real       = creal(field_in[gridpoint_s]);
		double imag       = cimag(field_in[gridpoint_s]);

		double intensity  = real*real + imag*imag;
		double phase      = atan2(imag,real);

		double s = get_s_from_index(gridpoint_s);
		double t = get_t_from_index(0);

	/* distance focus->origin */
		const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

		double A = s*t/c;
		double Arel = A - t + GEOM_S2;

		fprintf(fout, " %d  % .8e %14.8g %+.6f\n", Prof_num, Arel*1e6, intensity, phase);
	}
    fclose(fout);
	
    return 0;
}
int write_file_out()
{
    const char* fname;
	if (Prof_num==0) {	

		if (OPTIONS.fieldoutfname == NULL)
			//fname = "field_out.dat";
			fname = fnTTfield;
		else
			fname = OPTIONS.fieldoutfname;
			FILE* fout = fopen(fname, "w");
		if (fout == NULL)
			fatal_error("Could not open file for output!\n");
	} 
	if (OPTIONS.fieldoutfname == NULL)
		//fname = "field_out.dat";
		fname = fnTTfield;
	else
		fname = OPTIONS.fieldoutfname;
		FILE* fout = fopen(fname, "a+");
	if (fout == NULL)
		fatal_error("Could not open file for output!\n");

	for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=Rate_Column)

	{
		double real       = creal(field_out[gridpoint_s]);
		double imag       = cimag(field_out[gridpoint_s]);

		double intensity  = real*real + imag*imag;
		double phase      = atan2(imag,real);

		double s = get_s_from_index(gridpoint_s);
		double t = get_t_from_index(0);

	/* distance focus->origin */
		const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

		double A = s*t/c;
		double Arel = A - t + GEOM_S2;

		fprintf(fout, "%d  % .8e %14.8g %+.6f\n", Prof_num, Arel*1e6, intensity, phase);
	}
    fclose(fout);


    return 0;
}
int write_To_propag()
{
    const char* fname;

	if (OPTIONS.fieldoutfname == NULL)
		fname = fnTTfeature;
		//fname = "To_propagate.dat";
	else
		fname = OPTIONS.featurefname;
		FILE* fprop = fopen(fname, "w");
	if (fprop == NULL) {
	
		fatal_error("Could not open file for output!\n");
	}
		//fprintf(fprop, " %f \n", LAMBDA);
		fprintf(fprop, "%f \n", wavenumber*1e-6); // en um-1
		fprintf(fprop, "%d \n", l1);	// nb profil in transverse direction
		fprintf(fprop, "%f \n", l2 * subsampling); 	// nb profil in longitudinal direction
		fprintf(fprop, "%f \n", pix_sizeX / subsampling); //en micron  in longitudinal direction
		fprintf(fprop, "%f \n", THETA);
		fprintf(fprop, "%f \n", pix_sizeY); // en micron in transverse direction
		fprintf(fprop, "%d \n", Coh_step);
		fprintf(fprop, "%f \n", DistanceMIN);
		fprintf(fprop, "%f \n", DistanceMAX);
		fprintf(fprop, "%d \n", DistStep);
		
	
    fclose(fprop);


    return 0;
}
int write_SimuFile()
{
	const char* fname;

	if (OPTIONS.fieldoutfname == NULL)
		fname = fnMLsimulation;
	else
		fname = OPTIONS.featurefname;
		FILE* fprop = fopen(fname, "w");
	if (fprop == NULL) {
	
		fatal_error("Could not open file for output!\n");
	}
		// classic Computation feature 
		fprintf(fprop, "%f \n", GEOM_S1*1e3); // um
		fprintf(fprop, "%f \n", GEOM_S2*1e3); // um
		fprintf(fprop, "%f \n", THETA);	// rad
		fprintf(fprop, "%f \n", DOMAINSIZE_S / ILLUM_WIDTH);
		fprintf(fprop, "%f \n", LambdaBragg); // um
		fprintf(fprop, "%d \n", NUMBERLAYERS); 
		fprintf(fprop, "%f \n", LAYERRATIO);
		fprintf(fprop, "%f \n", MIRRORLENGTH);	// um
		fprintf(fprop, "%f \n", lambda); 	// um
		fprintf(fprop, "%d \n", GRIDPOINTS_S);
		fprintf(fprop, "%d \n", GRIDPOINTS_T);
			
	
		// info for coherence 
		fprintf(fprop, "%f \n", Coh_V);	// um
		fprintf(fprop, "%f \n", Coh_H);	// um
		fprintf(fprop, "%f \n", Coh_dist);	// um
		fprintf(fprop, "%d \n", Coh_step);
		
		// computation parameters
		fprintf(fprop, "%d \n", NxCAM);
		fprintf(fprop, "%d \n", NyCAM);
		fprintf(fprop, "%f \n", SxCAM);	// um
		fprintf(fprop, "%f \n", SyCAM);	// um
		fprintf(fprop, "%f \n", FWHMpsf);	// um
		fprintf(fprop, "%f \n", DistanceMIN);	// um
		fprintf(fprop, "%f \n", DistanceMAX);	// um
		fprintf(fprop, "%d \n", DistStep);
		
		
		
		// info of files
		//fprintf(fprop, "%s \n", fnTTfield);
		//fprintf(fprop, "%s \n", fnTTfeature); 
		

    fclose(fprop);


    return 0;
}
int write_file_psi0()
{
    const char* fname;
    if (OPTIONS.fieldpsi0fname == NULL)
	fname = "field_in2d.dat";
    else
	fname = OPTIONS.fieldpsi0fname;
    FILE* fout0 = fopen(fname, "w");
    if (fout0 == NULL)
	fatal_error("Could not open file for output!\n");
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=Rate_Column)
    {
	for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t+=1)
	{
	    double real      = creal(psi0[gridpoint_s][gridpoint_t]);
	    double imag      = cimag(psi0[gridpoint_s][gridpoint_t]);
	    double intensity = real*real + imag*imag;
	    double phase     = atan2(imag,real);

	    fprintf(fout0, "%8.5f %8.5f %g %.4f\n", gridpoint_s*GRIDSPACING_S*1e6, gridpoint_t*GRIDSPACING_T*1e6, intensity, phase);
	}
	fprintf(fout0, "\n");
    }
    fclose(fout0);

    return 0;
}
int write_file_psi1()
{
    const char* fname;
    if (OPTIONS.fieldpsi1fname == NULL)
	fname = "field_out2d.dat";
    else
	fname = OPTIONS.fieldpsi1fname;
    FILE* fout1 = fopen(fname, "w");
    if (fout1 == NULL)
	fatal_error("Could not open file for output!\n");
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=Rate_Column)
    {
	for (gridpoint_t=0; gridpoint_t<GRIDPOINTS_T; gridpoint_t+=1)
	{
	    double real      = creal(psi1[gridpoint_s][gridpoint_t]);
	    double imag      = cimag(psi1[gridpoint_s][gridpoint_t]);
	    double intensity = real*real + imag*imag;
	    double phase     = atan2(imag,real);

	    fprintf(fout1, "%8.5f %8.5f %g %.4f\n", gridpoint_s*GRIDSPACING_S*1e6, gridpoint_t*GRIDSPACING_T*1e6, intensity, phase);
	}
	fprintf(fout1, "\n");
    }
    fclose(fout1);

    return 0;
}
void put_to_shm()
{
    put_plotdata_to_shm();
    put_fieldout_to_shm();
}
void put_plotdata_to_shm()
{
    const int xmax     = GRIDPOINTS_S / PLOTEVERY;
    const int ymax     = GRIDPOINTS_T;
    const int ymax2    = 50 + 2*GRIDPOINTS_T;
    int x, y;

    key_t key;
    int shmid;
    int size;

    /* shm-key is stored in env variable */
    char key_buf[16];
    const char* s;

    getENV("config_shm_plotdata", &s);
    strncpy(key_buf, s, 16);

    key = atoi(key_buf);

    int ret;
    char ipcrm[32];
    snprintf(ipcrm, 32, "ipcrm -M 0x%x 2>/dev/null", key);
    ret = system(ipcrm);
    if (ret == 0) {;}

    size = xmax * ymax2 * sizeof(float);
    float* shm;

    if ((shmid = shmget(key, size, IPC_CREAT | 0640)) < 0)
	fatal_error("could not shmget()\n");

    if ((shm = (float*)shmat(shmid, NULL, 0)) == (float*) -1)
	fatal_error("could not shmat()\n");

    for (x=0; x<xmax; x++)
	for (y=0; y < ymax; y++)
	{
	    float real, imag, int_psi0, int_psi1;

	    real = creal(psi0[x*PLOTEVERY][y]);
	    imag = cimag(psi0[x*PLOTEVERY][y]);
	    int_psi0 = real*real+imag*imag;

	    real = creal(psi1[x*PLOTEVERY][y]);
	    imag = cimag(psi1[x*PLOTEVERY][y]);
	    int_psi1 = real*real+imag*imag;

	    shm[x*ymax2 + y + ymax+50] = int_psi0;
	    shm[x*ymax2 + y          ] = int_psi1;
	}

    shmdt(shm);
}
void put_fieldout_to_shm()
{
    const int EVERY = 10;
    const int smax  = GRIDPOINTS_S / EVERY;

    key_t key;
    int shmid;
    int size;

    /* shm-key is stored in env variable */
    char key_buf[16];
    const char* s;

    getENV("config_shm_fieldout", &s);
    strncpy(key_buf, s, 16);

    key = atoi(key_buf);

    int ret;
    char ipcrm[32];
    snprintf(ipcrm, 32, "ipcrm -M 0x%x 2>/dev/null", key);
    ret = system(ipcrm);
    if (ret == 0) {;}

    size = 4*smax * sizeof(float);
    float* shm;

    if ((shmid = shmget(key, size, IPC_CREAT | 0640)) < 0)
	fatal_error("could not shmget()\n");

    if ((shm = (float*)shmat(shmid, NULL, 0)) == (float*) -1)
	fatal_error("could not shmat()\n");


    /* distance focus->origin */
    const double c = 0.5*sqrt( GEOM_S1*GEOM_S1 + GEOM_S2*GEOM_S2 - 2*GEOM_S1*GEOM_S2 * cos(M_PI-2*THETA) );

    /* layout in shm:
     * shm[4*index+0] = relative coordinate A, centered at mirror's center
     * shm[4*index+1] = real part of reflected amplitude
     * shm[4*index+2] = imag part of reflected amplitude
     * shm[4*index+3] = 0.0 (unused)
     */

    int gridpoint_s;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s+=EVERY)
    {
	double real       = creal(field_out[gridpoint_s]);
	double imag       = cimag(field_out[gridpoint_s]);

	double s = get_s_from_index(gridpoint_s);
	double t = get_t_from_index(0);

	double A = s*t/c;
	double Arel = A - t + GEOM_S2;

	int index = gridpoint_s / EVERY;

	shm[4*index+0] = Arel*1e6;
	shm[4*index+1] = real;
	shm[4*index+2] = imag;
	shm[4*index+3] = 0.0;
    }

    shmdt(shm);
}

/* /input, output */



/* memory management and stuff */

int free_memory()
{
    int gridpoint_s;

    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	free(psi0[gridpoint_s]);
	free(psi1[gridpoint_s]);
	}

    free(local_theta_arr);
    free(local_thetaprime_arr);
    free(ALPHA2_rez);
    free(ALPHA2_arr);
    free(BETA2_arr);
    free(SPHERICALWAVE_arr_FE);
    free(SPHERICALWAVE_arr_SE);
    free(field_in);
    free(field_out);
    free(Dpsi0);
    free(Dpsi1);
    free(psi0);
    free(psi1);
	free(Hdev2);
    return 0;
}

int allocate_memory()
{
    int gridpoint_s, gridpoint_t;
    size_t allocated = 0;
    size_t allocating;

    allocating = sizeof(complex float*) * GRIDPOINTS_S;

    check_maxmem(&allocated, allocating);
    psi0 = (complex float **) malloc(allocating);
    check_maxmem(&allocated, allocating);
    psi1 = (complex float **) malloc(allocating);

     check_maxmem(&allocated, sizeof(float*) * GRIDPOINTS_S);
    Hdev = (float **) malloc(sizeof(float*) * GRIDPOINTS_S);

    if (psi0 == NULL || psi1 == NULL || Hdev == NULL )
	fatal_error("malloc failed!\n");
		
    // memory allocation of HDEV 
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	allocating = sizeof(complex float ) * GRIDPOINTS_T;

	check_maxmem(&allocated, allocating);
	psi0[gridpoint_s] = (complex float *) malloc(allocating);
	check_maxmem(&allocated, allocating);
	psi1[gridpoint_s] = (complex float *) malloc(allocating);


	// TODO / error Here it should be 
        // Hdev[gridpoint_s] = (float *) malloc( sizeof(float) * GRIDPOINTS_T);
        // But line ~350 is iterated on HDEV until i, j == l2
   	check_maxmem(&allocated, allocating);
    Hdev[gridpoint_s] = (float *) malloc( sizeof(float) * GRIDPOINTS_S);


	if (psi0[gridpoint_s] == NULL || psi1[gridpoint_s] == NULL || Hdev[gridpoint_s] == NULL ) 
	    fatal_error("malloc failed!\n");
    }

	//	memory allocation of  Hdev2
	allocating = sizeof(double *)*l1;
	/*printf (" subsampling : %f \n",subsampling);
	printf (" allocating1 : %li \n", allocating/sizeof(double *));
	printf("\n \n"); */
	Hdev2 = (double **) malloc(allocating); 
    check_maxmem(&allocated, allocating);
    
    
    allocating = sizeof(double)*(l2+1)*subsampling;
//	printf (" allocating2 : %li \n ", allocating/sizeof(double));
//	printf( " l2 *sub : %d \n",l2*subsampling);
	
    for (gridpoint_t=0; gridpoint_t<l1; gridpoint_t++)
    {

    Hdev2[gridpoint_t] = (double *) malloc(allocating);
	check_maxmem(&allocated, allocating);
	
	if(Hdev2[gridpoint_t] == NULL)
	    fatal_error("malloc failed!\n");
    }

    allocating = sizeof(complex float ) * GRIDPOINTS_T;

    check_maxmem(&allocated, allocating);
    Dpsi0 = (complex float *) malloc(allocating);
    check_maxmem(&allocated, allocating);
    Dpsi1 = (complex float *) malloc(allocating);
    if (Dpsi0 == NULL || Dpsi1 == NULL)
	fatal_error("malloc failed!\n");


    allocating = sizeof(complex float ) * GRIDPOINTS_S;

    check_maxmem(&allocated, allocating);
    field_out  = (complex float *) malloc(allocating);
    check_maxmem(&allocated, allocating);
    field_in   = (complex float *) malloc(allocating);

    if (field_out == NULL || field_in == NULL)
	fatal_error("malloc failed!\n");

    allocating = sizeof(double) * GRIDPOINTS_T;

    check_maxmem(&allocated, allocating);
    ALPHA2_arr = (double *) malloc(allocating);
    check_maxmem(&allocated, allocating);
    BETA2_arr  = (double *) malloc(allocating);
    check_maxmem(&allocated, allocating);
    ALPHA2_rez = (double *) malloc(allocating);

    if (ALPHA2_arr == NULL || BETA2_arr == NULL || ALPHA2_rez == NULL)
	fatal_error("malloc failed!\n");

    allocating = sizeof(double) * GRIDPOINTS_S;

    check_maxmem(&allocated, allocating);
    local_theta_arr = (double *) malloc(allocating);
    local_thetaprime_arr = (double *) malloc(allocating);

    if (local_theta_arr == NULL || local_thetaprime_arr == NULL)
	fatal_error("malloc failed!\n");

    allocating = sizeof(double) * GRIDPOINTS_T;

    check_maxmem(&allocated, allocating);
    SPHERICALWAVE_arr_FE = (double *) malloc(allocating);
    check_maxmem(&allocated, allocating);
    SPHERICALWAVE_arr_SE = (double *) malloc(allocating);

    if (SPHERICALWAVE_arr_FE == NULL || SPHERICALWAVE_arr_SE == NULL)
	fatal_error("malloc failed!\n");

    return 0;
}
void resource_usage()
{
    FILE* f = fopen("/proc/self/statm", "r");
    if (f == NULL)
	return;

    int numbers[7];
    int i;
    for (i=0; i<7; i++)
    {
	int ret;
	ret = fscanf(f, "%d", &numbers[i]);
	if (ret != 1)
	    continue;
    }

	if (FLAGS.colour == 1)
	{
	    highlight(highlight_emphasize);
	    fprintf(stderr, "virtual memory usage : %6.1f MB\n", 1.0*numbers[0]*4096/1048576);
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "virtual memory usage : %6.1f MB\n", 1.0*numbers[0]*4096/1048576);

    fclose(f);
    return;

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    fprintf(stderr, "resident size: %ld\n", usage.ru_isrss);
}
void time_needed()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_emphasize);
	    fprintf(stderr, "user time            : %7.2f s\n", usage.ru_utime.tv_sec + 1e-6*usage.ru_utime.tv_usec);
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "user time            : %7.2f s\n", usage.ru_utime.tv_sec + 1e-6*usage.ru_utime.tv_usec);
    fprintf(stderr, "system time          : %7.2f s\n", usage.ru_stime.tv_sec + 1e-6*usage.ru_stime.tv_usec);
}

int print_efficiency()
{

    double integral_incoming    = efficiency_integrate_incoming()    * GRIDSPACING_S;
    double integral_reflected   = efficiency_integrate_reflected()   * GRIDSPACING_S;
    double integral_transmitted = efficiency_integrate_transmitted() * GRIDSPACING_S;

    double reflectivity  = integral_reflected   / integral_incoming;
    double transmittance = integral_transmitted / integral_incoming;

    if (FLAGS.ascii == 0)
    {
	const char* buf1;

	buf1 = get_utf8(integral_incoming);
	fprintf(stderr, "total incoming       :   %s\n",     buf1);
	buf1 = get_utf8(integral_reflected);
	fprintf(stderr, "total reflected      :   %s\n",     buf1);
	buf1 = get_utf8(integral_transmitted);
	fprintf(stderr, "total transmitted    :   %s\n",     buf1);
    } else {
	fprintf(stderr, "total incoming       : %13.4e\n", integral_incoming    );
	fprintf(stderr, "total reflected      : %13.4e\n", integral_reflected   );
	fprintf(stderr, "total transmitted    : %13.4e\n", integral_transmitted );
    }
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_important);
	    fprintf(stderr, "reflectivity         : %11.6f\n", reflectivity         );
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "reflectivity         : %11.6f\n", reflectivity         );
	fprintf(stderr, "transmittance        : %11.6f\n", transmittance        );

    return 0;
}
double efficiency_integrate_incoming()
{
    double integral_incoming = 0.0;
    int gridpoint_t = 0;
    int gridpoint_s;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	double real      = creal(psi0[gridpoint_s][gridpoint_t]);
	double imag      = cimag(psi0[gridpoint_s][gridpoint_t]);
	double intensity = real*real + imag*imag;

	integral_incoming += intensity;
    }
    return integral_incoming;
}
double efficiency_integrate_reflected()
{
    double integral_reflected = 0.0;
    int gridpoint_t = 0;
    int gridpoint_s;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	double real      = creal(psi1[gridpoint_s][gridpoint_t]);
	double imag      = cimag(psi1[gridpoint_s][gridpoint_t]);
	double intensity = real*real + imag*imag;

	integral_reflected += intensity;
    }
    return integral_reflected;
}
double efficiency_integrate_transmitted()
{
    double integral_transmitted = 0.0;
    int gridpoint_t = GRIDPOINTS_T-1;
    int gridpoint_s;
    for (gridpoint_s=0; gridpoint_s<GRIDPOINTS_S; gridpoint_s++)
    {
	double real      = creal(psi0[gridpoint_s][gridpoint_t]);
	double imag      = cimag(psi0[gridpoint_s][gridpoint_t]);
	double intensity = real*real + imag*imag;

	integral_transmitted += intensity;
    }
    gridpoint_s = GRIDPOINTS_S-1;
    for (gridpoint_t=1; gridpoint_t<GRIDPOINTS_T; gridpoint_t++)
    {
	double real      = creal(psi0[gridpoint_s][gridpoint_t]);
	double imag      = cimag(psi0[gridpoint_s][gridpoint_t]);
	double intensity = real*real + imag*imag;

	integral_transmitted += intensity;
    }
    return integral_transmitted;
}

/* /memory management and stuff */



/* information output */

int fatal_error(const char* message)
{
    fprintf(stderr, "%s", message);
    exit(-1);
}

int fatal_error_maxmem()
{
    fprintf(stderr, "error: maxmem reached.\n");
    fprintf(stderr, "aborting.\n");
    exit(-1);
}

int check_maxmem(size_t* allocated, size_t allocating)
{
    if (OPTIONS.maxmem <= 0)
	return 0;

    if (*allocated+allocating > OPTIONS.maxmem)
	fatal_error_maxmem();
    else
	*allocated += allocating;
    return 0;
}

int print_geometry_information()
{
    double lambda      = 12.4e-10/ENERGY;
    double LambdaBragg = lambda / (2*sin(THETA));
    double thickness   = DOMAINSIZE_T / sin(THETA);
    
    fprintf(stderr, "--- geometry information ---\n");
    if (FLAGS.ascii == 0)
    {
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "S₁                   : %9.4f m   \n", GEOM_S1);
	    fprintf(stderr, "S₂                   : %9.4f m   \n", GEOM_S2);
	    highlight(highlight_default);
	}
	else
	{
	    fprintf(stderr, "S₁                   : %9.4f m   \n", GEOM_S1);
	    fprintf(stderr, "S₂                   : %9.4f m   \n", GEOM_S2);
	}
    } else {
	fprintf(stderr, "distance from source : %9.4f m   \n", GEOM_S1);
	fprintf(stderr, "distance to   focus  : %9.4f m   \n", GEOM_S2);
    }
    fprintf(stderr, "s-value at center    : %9.4f     \n", .5*(GEOM_S1-GEOM_S2));
    fprintf(stderr, "t-value at bottom    : %9.4f     \n", .5*(GEOM_S1+GEOM_S2));
    if (FLAGS.ascii == 0)
    {
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "Θ at center          : %9.4f mrad\n", 1e3*local_theta(GRIDPOINTS_S/2));
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "Θ at center          : %9.4f mrad\n", 1e3*local_theta(GRIDPOINTS_S/2));
	fprintf(stderr, "θ at left edge       : %9.4f mrad\n", 1e3*local_theta(             0));
	fprintf(stderr, "θ at right edge      : %9.4f mrad\n", 1e3*local_theta(GRIDPOINTS_S-1));
    } else {
	fprintf(stderr, "angle of incidence   : %9.4f mrad\n", 1e3*local_theta(GRIDPOINTS_S/2));
	fprintf(stderr, "   at left edge      : %9.4f mrad\n", 1e3*local_theta(             0));
	fprintf(stderr, "   at right edge     : %9.4f mrad\n", 1e3*local_theta(GRIDPOINTS_S-1));
    }
    if (FLAGS.ascii == 0)
    {
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "ΔΘ (Bragg-deviation) : %9.4f mrad\n", thetadiff*1e3);
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "ΔΘ (Bragg-deviation) : %9.4f mrad\n", thetadiff*1e3);
    } else {
	fprintf(stderr, "deviation from Bragg : %9.4f mrad\n", thetadiff*1e3);
    }
    fprintf(stderr, "peak expected near   : %9.4f mrad\n", -creal(U0)/wavenumber*1e5);
    fprintf(stderr, "\n");

    fprintf(stderr, "--- layer information ---\n");
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "mirror length        : %9.4f mm\n",  DOMAINSIZE_S*1e3);
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "mirror length        : %9.4f mm\n",  DOMAINSIZE_S*1e3);
    fprintf(stderr, "layer thickness      : %9.4f nm\n",  LambdaBragg*1e9);
    if (FLAGS.ascii == 0)
    {
	fprintf(stderr, "ML structure         : %9.4f µm\n",  thickness*1e6);
    } else {
	fprintf(stderr, "ML structure         : %9.4f microm\n",  thickness*1e6);
    }
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "number of layers     : %4d     \n",  (int)ceil(thickness / LambdaBragg));
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "number of layers     : %4d     \n",  (int)ceil(thickness / LambdaBragg));
    fprintf(stderr, "\n");

    return 0;
}
double local_theta(int gridpoint_s)
{
    double beta2 = BETA2(gridpoint_s, 0);
    double beta  = sqrt(beta2);
    double theta = asin(beta);
    return theta;
}

int print_simulation_information()
{
    double GRIDRATIO = DOMAINSIZE_S/GRIDPOINTS_S * GRIDPOINTS_T/DOMAINSIZE_T * THETA;

    fprintf(stderr, "--- simulation information ---\n");
    if (FLAGS.ascii == 0)
    {
	fprintf(stderr, "grid points, s       : %4d×10³ \n", GRIDPOINTS_S/1000);
    } else {
	fprintf(stderr, "grid points, s       : %4d ths \n", GRIDPOINTS_S/1000);
    }
    fprintf(stderr, "grid points, t       : %4d     \n", GRIDPOINTS_T     );
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_emphasize);
	    fprintf(stderr, "grid ratio           : %9.4f   \n", GRIDRATIO        );
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "grid ratio           : %9.4f   \n", GRIDRATIO        );

    resource_usage();
    fprintf(stderr, "\n");

    return 0;
}
int print_simulation_information_after()
{
    fprintf(stderr, "--- simulation information ---\n");
    time_needed();
    resource_usage();
    fprintf(stderr, "\n");
    return 0;
}

const char* get_exponent_utf8(int e)
{
    switch(e)
    {
	case 0: return "⁰";
	case 1: return "¹";
	case 2: return "²";
	case 3: return "³";
	case 4: return "⁴";
	case 5: return "⁵";
	case 6: return "⁶";
	case 7: return "⁷";
	case 8: return "⁸";
	case 9: return "⁹";
    }
    return "";
}

int parse_exponent_utf8(int expo, char* exponent)
{
    if (expo < 0)
    {
	sprintf(exponent, "⁻");
	expo = abs(expo);
    }
    if (expo > 0)
    {
	if (expo > 9)
	{
	    int e1 = expo / 10;
	    strcat(exponent, get_exponent_utf8(e1));
	}
	int e2 = expo % 10;
	strcat(exponent, get_exponent_utf8(e2));
    }

    return 0;
}

const char* get_utf8(double value)
{
    char buf1[32]     = "";
    char mantissa[32] = "";
    char exponent[32] = "";
    int  expo         = 0;
    static char buf2[32];
    strcpy(buf2, "");

    snprintf(buf1, 32, "% .4e",   value);
    sscanf(  buf1, "%7c",         mantissa);
    expo = atoi(buf1 + 8);
    parse_exponent_utf8(expo, exponent);
    snprintf(buf2, 32, "%s×10%s", mantissa, exponent );

    return buf2;
}

int print_optical_constants_information()
{
    fprintf(stderr, "--- optical constants information ---\n");
    fprintf(stderr, "material 1           : %s\n",              MATERIAL1             );
    fprintf(stderr, "material 2           : %s\n",              MATERIAL2             );
    fprintf(stderr, "E                    : %9.4f keV\n",       ENERGY                );

    if (FLAGS.ascii == 0)
    {
	if (FLAGS.colour == 1)
	{
	    highlight(highlight_userinput);
	    fprintf(stderr, "λ                    : %9.4f Å\n",         12.4/ENERGY           );
	    highlight(highlight_default);
	}
	else
	    fprintf(stderr, "λ                    : %9.4f Å\n",         12.4/ENERGY           );
    } else {
	fprintf(stderr, "lambda               : %9.4f Ang\n",       12.4/ENERGY           );
    }

    if (FLAGS.ascii == 0)
    {
	char buf1[32];
	char buf2[32];

	strncpy(buf1, get_utf8(creal(U0)), 32);
	strncpy(buf2, get_utf8(cimag(U0)), 32);
	fprintf(stderr, "U0                   :  (%s,%s)\n",  buf1, buf2);
	strncpy(buf1, get_utf8(creal(U1_init)), 32);
	strncpy(buf2, get_utf8(cimag(U1_init)), 32);
	fprintf(stderr, "initial U1                   :  (%s,%s)\n",  buf1, buf2);
	strncpy(buf1, get_utf8(creal(Um1_init)), 32);
	strncpy(buf2, get_utf8(cimag(Um1_init)), 32);
	fprintf(stderr, "initial Um1                  :  (%s,%s)\n",  buf1, buf2);
    } else {
	fprintf(stderr, "U0                   :  (% .3e,% .3e)\n",  creal(U0),  cimag(U0) );
	fprintf(stderr, "initial U1                   :  (% .3e,% .3e)\n",  creal(U1_init),  cimag(U1_init) );
	fprintf(stderr, "initial Um1                  :  (% .3e,% .3e)\n",  creal(Um1_init), cimag(Um1_init));
    }
/*
    fprintf(stderr, "Zachariasen g        :   %7.5f       \n",  cimag(U0) / creal(U1) );
#ifdef UTF8SYMBOLS
    fprintf(stderr, "Zachariasen κ        :   %7.5f       \n",  cimag(U1) / creal(U1) );
    fprintf(stderr, "           κ/g       :   %7.5f       \n",  cimag(U1) / cimag(U0) );
#else
    fprintf(stderr, "Zachariasen kappa    :   %7.5f       \n",  cimag(U1) / creal(U1) );
    fprintf(stderr, "          kappa/g    :   %7.5f       \n",  cimag(U1) / cimag(U0) );
#endif
*/
    fprintf(stderr, "\n");

    return 0;
}

int online_analysis()
{
    /* print incoming , reflected intensity, ratio */
    fprintf(stderr, "--- online analysis ---\n");
    print_efficiency();
    fprintf(stderr, "\n");
    return 0;
}

int progress_meter_init()
{
    if (FLAGS.progressbar == 0)
    {
	fprintf(stdout, "[ progress bar disabled ]\n");
	return -1;
    }

    const int number_points = 20;
    int i;

    fprintf(stdout, "[");
    for (i=0; i<number_points; i++)
	fprintf(stdout, " ");

    fprintf(stdout, "]\r[");
    fflush(stdout);
    fflush(stderr);
    return 0;
}
int progress_meter(int gridpoint_s, int gridpoints_s)
{
    if (FLAGS.progressbar == 0)
	return -1;

    const int number_points = 20;
    int every = gridpoints_s / number_points;

    if ( (gridpoint_s%every) == 0)
    {
	fprintf(stdout, "-");
	fflush(stdout);
    }

    return 0;
}
int progress_meter_finish()
{
    if (FLAGS.progressbar == 0)
    {
	fprintf(stdout, "\n");
	return -1;
    }

    fprintf(stdout, "-]\n\n"); 
    fflush(stdout);
    fflush(stderr);
    return 0;
}
int print_illumination_information()
{
    fprintf(stderr, "--- illumination information ---\n");
    switch (ILLUMINATION)
    {
	case SELF:
	    fprintf(stderr, "calculated by        : me\n");
	    fprintf(stderr, "center               :   %5.2f\n", DOMAINSIZE_S/ILLUM_CENTER);
	    fprintf(stderr, "width                :   %5.2f\n", DOMAINSIZE_S/ILLUM_WIDTH );
	    fprintf(stderr, "intensity            :   %5.2f\n", ILLUM_INTENS);
	    fprintf(stderr, "number               :   %2d  \n", ILLUM_NUMBER);
	    fprintf(stderr, "factor               :   %5.2f\n", ILLUM_FACTOR);
	    break;
	case fromFILE:
	    fprintf(stderr, "read in from file    : %s\n",      ILLUM_FILE);
	    break;
	case fromSHM:
	    fprintf(stderr, "read in from shm     : 0x%x\n",    ILLUM_SHM);
	    break;
    }
    fprintf(stderr, "\n");

    return 0;
}
int colour(colour_t col)
{
    switch (col)
    {
	case colour_black         : fprintf(stdout, "\033[22;30m"); fprintf(stderr, "\033[22;30m"); break;
	case colour_red           : fprintf(stdout, "\033[22;31m"); fprintf(stderr, "\033[22;31m"); break;
	case colour_green         : fprintf(stdout, "\033[22;32m"); fprintf(stderr, "\033[22;32m"); break;
	case colour_brown         : fprintf(stdout, "\033[22;33m"); fprintf(stderr, "\033[22;33m"); break;
	case colour_blue          : fprintf(stdout, "\033[22;34m"); fprintf(stderr, "\033[22;34m"); break;
	case colour_magenta       : fprintf(stdout, "\033[22;35m"); fprintf(stderr, "\033[22;35m"); break;
	case colour_cyan          : fprintf(stdout, "\033[22;36m"); fprintf(stderr, "\033[22;36m"); break;
	case colour_gray          : fprintf(stdout, "\033[22;37m"); fprintf(stderr, "\033[22;37m"); break;
	case colour_dark_gray     : fprintf(stdout, "\033[01;30m"); fprintf(stderr, "\033[01;30m"); break;
	case colour_light_red     : fprintf(stdout, "\033[01;31m"); fprintf(stderr, "\033[01;31m"); break;
	case colour_light_green   : fprintf(stdout, "\033[01;32m"); fprintf(stderr, "\033[01;32m"); break;
	case colour_yellow        : fprintf(stdout, "\033[01;33m"); fprintf(stderr, "\033[01;33m"); break;
	case colour_light_blue    : fprintf(stdout, "\033[01;34m"); fprintf(stderr, "\033[01;34m"); break;
	case colour_light_magenta : fprintf(stdout, "\033[01;35m"); fprintf(stderr, "\033[01;35m"); break;
	case colour_light_cyan    : fprintf(stdout, "\033[01;36m"); fprintf(stderr, "\033[01;36m"); break;
	case colour_white         : fprintf(stdout, "\033[01;37m"); fprintf(stderr, "\033[01;37m"); break;
	default                   : fprintf(stdout, "\033[0m"    ); fprintf(stderr, "\033[0m"    );
    }
    fflush(stdout);
    fflush(stderr);

    return col;
}
int highlight(highlight_t high)
{
    return colour(high);
}


/* /information output */

