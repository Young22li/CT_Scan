//Rotation image
//Edit By: Ziyang Li
//Date: 10/7/24
//Description: This file create three shape in one pics and rotate this from 0 - 180

#include <iostream>
#include "Display.h"
#include <cmath>
#include "line_functions.h"

using namespace std;

const int n = 180; // Image Size N x N


//Save the image after rotation
int tracy[n][n];

void square(int row1, int row2, int col1, int col2, int (&Pixels)[n][n])
{
    for (int i = row1 - 1; i < row2; i++)
    {
        for (int j = col1 - 1; j < col2; j++)
            Pixels[i][j] = 2;
    }
}

void circle(int x1, int y1, int r, int (&Pixels)[n][n])
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

void trangle(int h1, int h2, int b1, int b2, int bRow, int (&Pixels)[n][n])
{
    for (int i = b1; i < b2; i++)
    {
        line3(h1, h2, bRow, i, Pixels);
    }
}

void rotation(int angle, int (&Pixels)[n][n])
{
    //The original three variables would have been other number
    int width = n; // using 'constant' rather than literals within the code
    int height = n;
    int background = 188; // this is the background color - use a suitable value here


    //The original, angle shoule be theta
    float rads = angle * 3.1415926 / 180.0; // fixed constant PI
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
            tracy[r][c] = p; // Pixels and tracy is the name of the array we created ourselves.
                             //The original used a different name, but I forget.
        }
    }
}



int main()
{
    int Pixels[n][n];

    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            Pixels[row][col] = 188;
        }
    }

    // Creat three shape
    square(23, 65, 33, 65, Pixels);
    circle(88, 88, 20, Pixels);
    trangle(100, 123, 100, 150, 150, Pixels);

    Display image(Pixels, "output1");

    //Rotation
    rotation(15, Pixels);
    //rotation(30);
    //rotation(45);
    //rotation(60);
    //rotation(75);
    //rotation(90);
    //rotation(105);
    //rotation(120);
    //rotation(145);
    //rotation(165);
    //rotation(180);

    Display image2(tracy, "output2");

    return 0;
}
