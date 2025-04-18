// Back-projection
// Edit By: Ziyang Li
// Date: 12/31/24
// Description: This file using back-projection reconstructing from sinogram
// re function is reconstructed based on the information in each row of the sinogram
// backprojection function is an optimized version of the bp algorithm. The main idea is based on the reverse derivation of the rotation formula.

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

int final[n][n] = {0};



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

void re(int angle, int (&s)[n][n]){
    int temp[n][n];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            temp[j][i] = s[angle][i];
        }
        rotation2(-angle, temp);
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            final[i][j] += tracy[i][j];
        }
    }
}

void reduce(int (&s)[n][n]){
    int max1 = 0;

    for(int i = 0; i < n; i++){
        int max2 = 0;
        for(int j = 0; j < n; j++) max2 = max(max2, s[i][j]);
        max1 = max(max1, max2);
    }
    
    for(int i = 0; i < n; i++){
        
        for(int j = 0; j < n; j++) s[i][j] = s[i][j] * 255 / max1;
        
    }

}

void backProject(int (&sinogram)[n][n]) {
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
                    final[x][y] += sinogram[angleIndex][t];
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

    Display image(Pixels, "output1");

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
    
    Display image3(sinogram, "output2");


    //for(int i = 0; i < n; i++) re(i, sinogram);

    backProject(sinogram);

    reduce(final);

    Display image2(final, "output");

    return 0;
}
