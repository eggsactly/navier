#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/**
 * u is the u function for the joukowski equations
 * @param r the parametric variables
 * @param a the x offset
 * @param b the y offset
 * @param polarity + or - for the solution of the sqrt
 */
float u(float r, float a, float b, float polarity) {
	return a + (polarity * sqrt(1 - pow(r - b, 2)));
}

/**
 * xHat is the x component of the paramatric equation for Joukowski
 * @param t the parametric variable
 * @param a the x offset
 * @param b the y offset
 * @param polarity + or - for the solution of the sqrt
 */
float xHat(float t, float a, float b, float polarity) {
	return (pow(u(t, a, b, polarity), 3) + (u(t, a, b, polarity) * pow(t, 2)) + u(t, a, b, polarity)) / (pow(u(t, a, b, polarity), 2) + pow(t, 2));
}

/**
 * yHat is the y component of the paramatric equation for Joukowski
 * @param t the parametric variable
 * @param a the x offset
 * @param b the y offset
 * @param polarity + or - for the solution of the sqrt
 */
float yHat(float t, float a, float b, float polarity) {
	return ((pow(u(t, a, b, polarity), 2) * t) + pow(t, 3) - t) / (pow(u(t, a, b, polarity), 2) + pow(t, 2));
}

/**
 * paintBucket is a recursive algorithm which fills in a group of pixels with
 * the same color
 *
 * @param img pointer to the image to fill in
 * @param height height of the image img
 * @param width width of the image img
 * @param x x coordinate of the pixel to fill in
 * @param y y coordinate of the pixel to fill in
 */
void paintBucket(char ** img, int height, int width, int x, int y)
{
    if(x > 0 && x < width && y > 0 && y < height && img[x][y] == C_F)
    {
        img[x][y] = C_B;
        paintBucket(img, height, width, x+1, y);
        paintBucket(img, height, width, x-1, y);
        paintBucket(img, height, width, x, y+1);
        paintBucket(img, height, width, x, y-1);
    }
}

typedef struct {
    int x;
    int y;
} pixel;

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void init_flag(char **flag, int imax, int jmax, double delx, double dely,
    int *ibound)
{
    const float a = 0.2;
    const float b = -0.1;
    float scale = 100;
    const float angle = 12;

    float stpSze = 0.001;

    int i, j;
    double mx, my, x, y, rad1;

    /* Mask of a circular obstacle */
    mx = 20.0/41.0*jmax*dely;
    my = mx;
    rad1 = 5.0/41.0*jmax*dely;

    // Center of rotation
    const float px = 0;
    const float py = jmax/2;

    // Offset after rotation
    const int ox = -20;
    const int oy = -20;

    float rad;
    float thet;

    // Fill the array
    for (unsigned long i = 0; i < imax * jmax; i++)
    {
        flag[0][i] = C_F;
    }

    // Build up the size of the airfoil dots
    unsigned samples = 0;
    for(float t = -2.0; t < 2.0; t+=stpSze)
    {
        stpSze = 0.01/(2.1*2.1) * (pow(2.1, 2) - pow(t, 2));
        samples++;
    }

    // Allocate the dot array to use to fill the spaces between
    pixel * dotRecord = (pixel *) malloc (samples * sizeof(pixel));
    
    // Initialize dotRecord
    for(unsigned i = 0; i < samples; i++)
    {
        dotRecord[i].x = 0;
        dotRecord[i].y = 0;
    }
    
    samples = 0;
    // populate the dots
    for(float t = -2.0; t < 2.0; t+=stpSze)
    {
        stpSze = 0.01/(2.1*2.1) * (pow(2.1, 2) - pow(t, 2));

        x = (scale * xHat(t, a, b, -1)) + imax/2;
        y = (scale * yHat(t, a, b, -1)) + jmax/2;

        rad = sqrt(pow(x - px, 2) + pow(y - py, 2));
        thet = atan((py - y)/(px - x));

        y = (rad * sin(thet+(angle * M_PI / 180))) + py + oy;
        x = (rad * cos(thet+(angle * M_PI / 180))) + px + ox;

        if(x >= 0 && x < imax && y >= 0 && y < jmax)
        {
            flag[(int)x][(int)y] = C_B;
            dotRecord[samples].x = (int)x;
            dotRecord[samples].y = (int)y;
            samples++;
        }
    }

    // Fill in the space between the dots 
    for(int i = 0; i < samples; i++)
    {
        int nextIndex = (i + 1) % samples;

        int xDiff = dotRecord[i].x - dotRecord[nextIndex].x;
        int yDiff = dotRecord[i].y - dotRecord[nextIndex].y;

        
        // If it's the same pixel or a neighbor, do nothing
        if(abs(xDiff) < 2 && abs(yDiff) < 2)
        {

        }
        // If the slope is infinite or near infinite
        else if (abs(xDiff) <= 1)
        {
            for(int j = 0; j < abs(yDiff); j++)
            {
                if(yDiff < 0)
                {
                    flag[dotRecord[i].x][dotRecord[i].y + j] = C_B;
                }
                else
                {
                    flag[dotRecord[nextIndex].x][dotRecord[nextIndex].y + j] = C_B;
                }
            }
        }
        // In the general case
        else
        {
            // Straddle along the x axis
            for(int j = 0; j < abs(xDiff); j++)
            {
                int xs, ys;
                if(xDiff < 0)
                {
                    xs = dotRecord[i].x + j;
                    ys = dotRecord[i].y + ((j * yDiff)/xDiff);
                }
                else
                {
                    xs = dotRecord[i+1].x + j;
                    ys = dotRecord[i+1].y + ((j * yDiff)/xDiff);
                }
                if(xs > 0 && xs < imax && ys > 0 && ys < jmax)
                {
                    flag[xs][ys] = C_B;
                }
            }
        }
    }

    // Fill in the cavity
    x = (scale * xHat(0, a, b, -1)) + imax/2;
    y = (scale * yHat(0, a, b, -1)) + jmax/2;

    rad = sqrt(pow(x - px, 2) + pow(y - py, 2));
    thet = atan((py - y)/(px - x));

    y = (rad * sin(thet+(angle * M_PI / 180))) + py + oy;
    x = (rad * cos(thet+(angle * M_PI / 180))) + px + ox;
    paintBucket(flag, jmax, imax, (int) x + 2, (int) y);
    
    // fix strange issue where right hand side is filled in
    for(i = imax-100; i <= imax+1; i++)
    {
        for(j = 0; j <= jmax; j++)
        {
            flag[i][j] = C_F;
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
    *ibound = 0;
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (!(flag[i][j] & C_F)) {
                (*ibound)++;
                if (flag[i-1][j] & C_F) flag[i][j] |= B_W;
                if (flag[i+1][j] & C_F) flag[i][j] |= B_E;
                if (flag[i][j-1] & C_F) flag[i][j] |= B_S;
                if (flag[i][j+1] & C_F) flag[i][j] |= B_N;
            }
        }
    }
}
