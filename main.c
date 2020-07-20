#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <stdbool.h>

#include "alloc.h"
#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"
#include "libBitmapWag.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

   Extended by G. Weaver (2020) to use bitmap images as sim input. 
*/

// APP_NAME is used for stderr outputs in this program
// it should be provided by the makefile, if not, it is defined here. 
#ifndef APP_NAME
    #define APP_NAME "navier"
#endif

/**
 * @brief This structure contains a list of options a user can control on the 
 *        command line
 */
typedef struct {
    bool help;           /**< Indicates whether to print the list of command 
                              line options to the user */
    char* inputName;     /**< Indicates path to input file */
    bool hasInput;       /**< Indicates whether and input has been given */
    char* outputName;    /**< Indicates the output name */
    int outputFrequency; /**< Indicates the output frequency */
} argsList;

/**
 * parseArgs implements command line input parsing
 *
 * @param argc number of input arguments
 * @param argv array of arguments 
 * @return struct argsList containing the input options 
 */
argsList parseArgs(const int argc, const char** argv){
    int c;
    int option_index = 0;
    const static struct option long_options[] =
    {
        {"help",        no_argument,       0, 'h'},
        {"input-file",  required_argument, 0, 'i'},
        {"output-name", required_argument, 0, 'o'},
        {"frequency",   required_argument, 0, 'f'},
        {0, 0, 0, 0}
    };
    argsList output;
    
    output.help = false;
    output.inputName = "";
    output.hasInput = false;
    output.outputName = "output";
    output.outputFrequency = 1;

    
    while ((c = getopt_long(argc, argv, "hi:o:f:", long_options, 
        &option_index)) != -1)
    {
        switch (c)
        {
            case 'h':
                output.help = true;
                break;
            case 'i':
                output.inputName = optarg;
                output.hasInput = true;
                break;
            case 'o':
                output.outputName = optarg;
                break;
            case 'f':
                output.outputFrequency = atoi(optarg);
                break;
            default:
                abort ();
                break;
        }
    }
    return output;
}

int main(int argc, char *argv[])
{
    argsList inputArgs = parseArgs(argc, argv);

    int verbose = 1;          /* Verbosity level */
    double xlength = 22.0;    /* Width of simulated domain */
    double ylength = 4.1;     /* Height of simulated domain */
    int imax = 660;           /* Number of cells horizontally */
    int jmax = 120;           /* Number of cells vertically */

    char *outname;
    int output = 0;
    int output_frequency = 0;

    double t_end = 40; //2.1       /* Simulation runtime */
    double del_t = 0.003;      /* Duration of each timestep */
    double tau = 0.5;          /* Safety factor for timestep control */

    int itermax = 100;        /* Maximum number of iterations in SOR */
    double eps = 0.001;        /* Stopping error threshold for SOR */
    double omega = 1.7;        /* Relaxation parameter for SOR */
    double gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    double Re = 150.0;         /* Reynolds number */
    double ui = 1.0;           /* Initial X velocity */
    double vi = 0.0;           /* Initial Y velocity */

    double t, delx, dely;
    int  i, j, itersor = 0, ifluid = 0, ibound = 0;
    double res;
    double **u, **v, **p, **rhs, **f, **g;
    char  **flag;
    int init_case, iters = 0;
    int show_help = 0, show_usage = 0, show_version = 0;

    // Declare two bitmap image structs, img for writing, img2 for reading 
    BitmapWagImg img;

    // BitmapWagError is an enum containing error codes
    BitmapWagError error;

    if(inputArgs.help) {
		printf("Command line flags:\n");
		printf("\t-i, --input-file  [Filename]: Input bitmap file of airfoil\n");
		printf("\t-o, --output-name [name]    : Output folder\n");
		printf("\t-f, --frequency   [int]     : Output Frequency\n");
		printf("\t-h, --help                 : Shows this help option\n");
		return EXIT_SUCCESS;
	}


    outname = inputArgs.outputName;
    output_frequency = inputArgs.outputFrequency;
    output = 1;

    if(inputArgs.hasInput)
    {
        // Try reading the bitmap 
        error = ReadBitmapWag(&img, inputArgs.inputName);
        if(error)
        {
            fprintf(stderr, "%s: error: ReadBitmapWag: %s.\n", APP_NAME, 
                ErrorsToStringBitmapWag(error));
            return EXIT_FAILURE;
        }
     
        fprintf(stderr, "%s: info: Bitmap: %s opened.\n", APP_NAME, 
            inputArgs.outputName);

        imax = GetBitmapWagWidth(&img); 
        jmax = GetBitmapWagHeight(&img); 
    }
    

    delx = xlength/imax;
    dely = ylength/jmax;

    /* Allocate arrays */
    u    = alloc_doublematrix(imax+2, jmax+2);
    v    = alloc_doublematrix(imax+2, jmax+2);
    f    = alloc_doublematrix(imax+2, jmax+2);
    g    = alloc_doublematrix(imax+2, jmax+2);
    p    = alloc_doublematrix(imax+2, jmax+2);
    rhs  = alloc_doublematrix(imax+2, jmax+2); 
    flag = alloc_charmatrix(imax+2, jmax+2);                    

    if (!u || !v || !f || !g || !p || !rhs || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
        return 1;
    }

    unsigned long checker = 0;
    double checker1 = 0.0;

    // Set up initial values
    for (i=0;i<=imax+1;i++) {
         for (j=0;j<=jmax+1;j++) {
            checker += (i*jmax)+ j + 1;
         checker1 += (i*jmax) + j + 1.0;
             u[i][j] = ui;
             v[i][j] = vi;
             p[i][j] = 0.0;
         }
     }

    
    // Populate init_flag
    if(inputArgs.hasInput)
    {
        BitmapWagRgbQuad color;

        // Fill the array with fluids
        for (unsigned long i = 0; i < (imax+2) * (jmax+2); i++)
        {
            flag[0][i] = C_F;
        }
    
        // get and set the color
        for(unsigned i = 0; i < imax; i++)
        {
            for(unsigned j = 0; j < jmax; j++)
            {
                error = GetBitmapWagPixel(&img, i, j, &color);

                if(error)
                {
                    fprintf(stderr, "%s: error: GetBitmapWagPixel %d %d: %s.\n", 
                        APP_NAME, i, j, ErrorsToStringBitmapWag(error));
                    return -1;
                }

                // Anywhere that's white should be filled with fluid
                if(color.rgbBlue == 0xFF && color.rgbRed == 0xFF 
                    && color.rgbBlue == 0xFF)
                {
                    flag[i][jmax - j] = C_F; 
                }
                // Any pixel that is not white should be a boundary
                else
                {
                    flag[i][jmax - j] = C_B; 
                }
            }
        }

        /* Mark the north & south boundary cells */
        for (i=0; i<=imax+1; i++) {
            flag[i][0]      = C_B;
            flag[i][jmax+1] = C_B;
        }
        /* Mark the east and west boundary cells */
        for (j=1; j<=jmax; j++) {
            flag[0][j]      = C_B;
            flag[imax+1][j] = C_B;
        }

        /* flags for boundary cells */
        ibound = 0;
        for (i=1; i<=imax; i++) {
            for (j=1; j<=jmax; j++) {
                if (!(flag[i][j] & C_F)) {
                    (ibound)++;
                    if (flag[i-1][j] & C_F) flag[i][j] |= B_W;
                    if (flag[i+1][j] & C_F) flag[i][j] |= B_E;
                    if (flag[i][j-1] & C_F) flag[i][j] |= B_S;
                    if (flag[i][j+1] & C_F) flag[i][j] |= B_N;
                }
            }
        }

        error = FreeBitmapWag(&img);
        if(error)
        {
            fprintf(stderr, "%s: warning: FreeBitmapWag: %s.\n", APP_NAME, 
                ErrorsToStringBitmapWag(error));
        }
    }
    else
    {
        init_flag(flag, imax, jmax, delx, dely, &ibound);
    }

    

    apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
    
    // Main loop

    for (t = 0.0; t < t_end; t += del_t, iters++) {
        set_timestep_interval(&del_t, imax, jmax, delx, dely, u, v, Re, tau);

        ifluid = (imax * jmax) - ibound;

        compute_tentative_velocity(u, v, f, g, flag, imax, jmax,
            del_t, delx, dely, gamma, Re);

        compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);

        if (ifluid > 0) {
            itersor = poisson(p, rhs, flag, imax, jmax, delx, dely,
                        eps, itermax, omega, &res, ifluid);
        } else {
            itersor = 0;
        }

         printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
                iters, t+del_t, del_t, itersor, res, ibound);

    
        update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);

        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);

        if (output && (iters % output_frequency == 0)) 
        {
          write_ppm(u, v, p, flag, imax, jmax, xlength, ylength, outname,
                iters, output_frequency);
        
        }
    }

    free_matrix(u);
    free_matrix(v);
    free_matrix(f);
    free_matrix(g);
    free_matrix(p);
    free_matrix(rhs);
    free_matrix(flag);

    return EXIT_SUCCESS;
}

// Used for comparing computations when debugging other implementations

unsigned int simplest_checksum_char(char** in, int imax, int jmax) {
  unsigned int checksum = 0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*(i);
    }
  }
  return checksum;
}

double simplest_checksum(double** in, int imax, int jmax) {
  double checksum = 0.0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*((double)(i*jmax)+j);
    }
  }
  return checksum;
}
