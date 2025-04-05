// Iterative Reconstruction
// Edit By: Ziyang Li
// Date: 4/5/25
// Description: This file using Iterative Reconstruction algorithm reconstructing

#include <iostream>
#include "Display.h"
#include <cmath>
#include <algorithm>
#include "line_functions.h"
using namespace std;

const int n = 180; // Image Size N x N
int Pixels[n][n];

// Save the image after rotation
int tracy[n][n];


void square(int row1, int row2, int col1, int col2)
{
    for (int i = row1 - 1; i < row2; i++)
    {
        for (int j = col1 - 1; j < col2; j++)
            Pixels[i][j] = 250;
    }
}

void circle(int x1, int y1, int r)
{
    for (int x = x1 - r; x < r + x1; x++)
    {
        int dx = x - x1;
        int right = r * r - dx * dx;
        if (right >= 0)
        {
            int y = sqrt(right);
            line3(x, y + y1, x, y1 - y, Pixels);
        }
    }
}

void trangle(int h1, int h2, int b1, int b2, int bRow)
{
    for (int i = b1; i < b2; i++)
    {
        line3(h1, h2, bRow, i, Pixels);
    }
}

void rotation(int angle, int (&Pixels)[n][n])
{
    int width = n; // using 'constant' rather than literals within the code
    int height = n;
    int background = 2; // this is the background color - use a suitable value here

    float rads = angle * 3.14159265359 / 180.0; // fixed constant PI
    float cs = cos(-rads);                  // precalculate these values
    float ss = sin(-rads);
    float xcenter = (float)(width) / 2.0;
    float ycenter = (float)(height) / 2.0;
    for (int r = 0; r < height; r++)
    {
        for (int c = 0; c < width; c++)
        {
            // now find the pixel of the original image that is rotated to (r,c)
            // rotation formula assumes that origin = top-left and y points down
            int rorig = ycenter + ((float)(r)-ycenter) * cs - ((float)(c)-xcenter) * ss;
            int corig = xcenter + ((float)(r)-ycenter) * ss + ((float)(c)-xcenter) * cs;
            // now get the pixel value if you can
            int p = background; // in case there is no original pixel
            if (rorig >= 0 && rorig < height && corig >= 0 && corig < width)
            {
                p = Pixels[rorig][corig];
            }
            tracy[r][c] = p;
        }
    }
}


// void re(int angle, int (&s)[n][n]){
//     int temp[n][n];
//     for(int i = 0; i < n; i++){
//         for(int j = 0; j < n; j++){
//             temp[j][i] = s[angle][i];
//         }
//         rotation(-angle, temp);
//     }

//     for(int i = 0; i < n; i++){
//         for(int j = 0; j < n; j++){
//             final[i][j] += tracy[i][j];
//         }
//     }
// }

void reduce(int (&s)[n][n])
{
    int max1 = 0;
    int min_val = s[0][0];

    for (int i = 0; i < n; i++)
    {
        int max2 = 0;
        for (int j = 0; j < n; j++)
        {
            min_val = min(min_val, s[i][j]); 
            max2 = max(max2, s[i][j]);
        max1 = max(max1, max2);
        }
        
    }

    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < n; j++)
            s[i][j] = s[i][j] * 255 / max1;
    }
}




void backProject(int (&sinogram)[n][n], int (&f)[n][n]) {
    // // Initialize the reconstructed image to zero
    // for (int row = 0; row < SIZE; row++) {
    //     for (int col = 0; col < SIZE; col++) {
    //         ReconstructedImage[row][col] = 0;
    //     }
    // }
    double PI = 3.1415926;
    // Backprojection loop: iterate over angles (0° to 180° counterclockwise)
    for (int angleIndex = 0; angleIndex < n; angleIndex++) {
        //double angle = (angleIndex * PI / n); // Convert angle index to radians
        double angle = -(angleIndex * PI / n) + PI / 2;
        double cosTheta = cos(angle);
        double sinTheta = sin(angle);

        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                // Compute the sinogram position corresponding to this pixel
                 int t = static_cast<int>(round((x - n / 2) * cosTheta + (y - n / 2) * sinTheta)) + n / 2;


                // Accumulate the value from the sinogram if t is valid
                if (t >= 0 && t < n) {
                    f[x][y] += sinogram[angleIndex][t];
                }
            }
        }
    }

}

int main()
{
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            Pixels[row][col] = 2;
        }
    }

    // Creat three shape
    square(23, 65, 33, 65);
    circle(88, 88, 20);
    trangle(100, 123, 100, 150, 150);

    Display image(Pixels, "image");

    int sinogram[n][n];

    for (int i = 0; i < n; i++)
    {
        rotation(i, Pixels);
        for (int row = 0; row < n; row++)
        {
            int res = 0;

            for (int col = 0; col < n; col++)
            {
                res += tracy[col][row];
            }
            sinogram[i][row] = res;
        }
    }

    reduce(sinogram);
    
    Display image3(sinogram, "sinogram");

    //for(int i = 0; i < n; i++) re(i, sinogram);
    int final[n][n] = {0};

    backProject(sinogram, final);

    reduce(final);

    Display image2(final, "First_Re");

    //------------------------------------------------------------------------------------

    int num_iterations = 180;
double prev_error = 1e6;

for(int iter = 0; iter < num_iterations; iter++){
    int sinogram2[n][n] = {0};
    // Get the sinogram from new images
    for (int i = 0; i < n; i++) {
        rotation(i, final);
        for (int row = 0; row < n; row++) {
            int res = 0;
            for (int col = 0; col < n; col++)
                res += tracy[col][row];
            sinogram2[i][row] = res;
        }
    }
    reduce(sinogram2);
    Display i2(sinogram2, "Sinogram2_150");

    // // cal the error
     int error[n][n] = {0};
    // for(int i = 0; i < n; i++)
    //     for(int j = 0; j < n; j++)
    //         error[i][j] = sinogram[i][j] - sinogram2[i][j];
    // Display i3(sinogram2, "Sinogram_error_300");

    // Calculate the error using division instead of subtraction
        // We need to handle division by zero and very small values
        for(int x = 0; x < n; x++){ // Only process where x < 180
            for(int y = 0; y < n; y++){
                // To prevent division by zero, add a small constant
                const double epsilon = 0.0001; // Small constant to prevent division by zero
                
                if (sinogram2[x][y] < epsilon) {
                    // Avoid division by very small values
                    error[x][y] = 0;
                } else {
                    // Division-based error calculation
                    // We can use Pixels[x][y] / sinogram2[x][y] - 1.0
                    // This gives 0 when the values are equal, positive when Pixels > sinogram2,
                    // and negative when Pixels < sinogram2
                    double ratio = (double)sinogram[x][y] / sinogram2[x][y];
                    // Scale the ratio to a suitable range for your algorithm
                    error[x][y] = (int)((ratio - 1.0) * 100); // Adjust scaling factor as needed
                }
            }
        }
        Display i3(error, "Sinogram_error_150");
    // Checking the global error
    double current_error = 0.0;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            current_error += error[i][j] * error[i][j];
    current_error = sqrt(current_error) / (n * n);
    if(abs(current_error - prev_error) < 1e-5){
        cout << "Converged at iteration " << iter << endl;
        break;
    }
    prev_error = current_error;

    // backprojection error image
    int delta[n][n] = {0};
    backProject(error, delta);

    // delta
    double max_delta = 0.0;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if(abs(delta[i][j]) > max_delta)
                max_delta = abs(delta[i][j]);
    double alpha = 1.0 / max_delta;

    // Update
    double lambda = 0.1; // regularization
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++){
            double reg_term = 0;
            if(i > 0) reg_term += final[i][j] - final[i-1][j];
            if(j > 0) reg_term += final[i][j] - final[i][j-1];
            final[i][j] += alpha * (delta[i][j] - lambda * reg_term);
            final[i][j] = max(0, final[i][j]);
        }
}


reduce(final);
Display image24(final, "Final_output150");
   
     return 0;
}
