// FBP
// Edit By: Ziyang Li
// Date: 12/31/24
// Description: Reconstruction of a brain ct image using the fbp algorithm

#include <iostream>
#include <fstream>
#include "Display.h"
using namespace std;
#include <vector>
#include <cmath>

const int n = 400;

int tracy[n][n];

int final[n][n] = {0};

int sino[n][n];

int si[n][n];

void rotation(int angle, int** Pixels)
{
    int width = n; // using 'constant' rather than literals within the code
    int height = n;
    int background = 2; // this is the background color - use a suitable value here

    float rads = angle * 3.14159265359 / 180.0; // fixed constant PI
    float cs = cos(-rads);                      // precalculate these values
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


void reduce(int (&s)[n][n])
{
    int max1 = 0;

    for (int i = 0; i < n; i++)
    {
        int max2 = 0;
        for (int j = 0; j < n; j++)
            max2 = max(max2, s[i][j]);
        max1 = max(max1, max2);
    }

    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < n; j++)
            s[i][j] = s[i][j] * 255 / max1;
    }
}

void dft(int row, vector<double> &real, vector<double> &imag, int (&sinogram)[n][n])
{
    // for X[0]
    int sum = 0;
    for (int i = 0; i < n; i++)
        sum += sinogram[row][i];
    real.push_back(sum);
    imag.push_back(0);
    // cout << nums[0] << ":  "<< accumulate(nums.begin(), nums.end(), 0) << " + 0j" << endl;

    double pi = 3.1415926;

    for (int i = 1; i < n; i++)
    {
        double nonJ = 0.00, jPart = 0.00;
        for (int j = 0; j < n; j++)
        {
            if (j == 0)
            {
                nonJ += sinogram[row][j];
                continue;
            }
            double exponent = (2 * pi * j * i) / n;
            nonJ += sinogram[row][j] * cos(exponent);
            jPart += -(sinogram[row][j] * sin(exponent));
        }
        real.push_back(nonJ);
        imag.push_back(jPart);
    }
}

void idft(int row, vector<double> &real, vector<double> &imag, int (&sinogram)[n][n])
{
    double pi = 3.1415926;
    for (int k = 0; k < n; k++)
    {
        double sumReal = 0.0;
        double sumImag = 0.0;
        for (int j = 0; j < n; j++)
        {
            double angle = (2 * pi * j * k) / n;
            sumReal += real[j] * cos(angle) - imag[j] * sin(angle);
            sumImag += real[j] * sin(angle) + imag[j] * cos(angle);
        }
        // Divide by n to complete the inverse transform
        double res = sumReal / n;
         int resInt = static_cast<int>(ceil(res));
        //int resInt = static_cast<int>(floor(res));
        // cout<<res<<" ";

        sino[row][k] = resInt > 0 ? resInt : 0;
        //sino[row][k] = abs(resInt);
    }
    // return ans;
}

void createFilter(double (&rampFilter)[n])
{
    double c;
    for (int i = 0; i < n; i++) {
        c = abs(i - n / 2);        
        rampFilter[i] = (n - 2 * c) / n;
    }
}


void backProject(int (&sinogram)[n][n]) {
    double PI = 3.1415926;
    int sinogramRows = 180; 
    int sinogramCols = 400; 

    // Backprojection loop: iterate over angles (0° to 180° counterclockwise)
    for (int angleIndex = 0; angleIndex < sinogramRows; angleIndex++) {
        
        double angle = -(angleIndex * PI / sinogramRows) + PI / 2;
        double cosTheta = cos(angle);
        double sinTheta = sin(angle);

        // 遍历重建图片的每个像素点
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                
                int t = static_cast<int>(round((x - n / 2) * cosTheta + (y - n / 2) * sinTheta)) + sinogramCols / 2;

                
                if (t >= 0 && t < sinogramCols) {
                    final[x][y] += sinogram[angleIndex][t];
                }
            }
        }
    }
}


int main() {
    int Pixels[n][n];
	ifstream inFile;

	inFile.open("sinogram_400x400.txt");
		
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			inFile >> Pixels[row][col];
	
	inFile.close();


     for (int i = 0; i < n; i++)
    {
        vector<double> real, imag;
        dft(i, real, imag, Pixels);
        double rampFilter[n];
        createFilter5(rampFilter);
        
        
        for(int j = 0; j < n; j++){
            real[j] *= rampFilter[j];
            imag[j] *= rampFilter[j];
        }

        idft(i, real, imag, Pixels);
    }

    reduce(sino);


   backProject(sino);
    reduce(final);
    

    Display image2(final, "output3");

	return 0;
}
