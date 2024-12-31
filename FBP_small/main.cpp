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

int sinogram[n][n];

int sino[n][n];



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
        rotation(-angle, temp);
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

void sp(int row) {
    //Window[-1, 2, -1]
    for (int col = 0; col < n; col++) {
        int ans = 0, count = 0;
        
        for (int i = col - 1; i <= col + 1; i++) {

            int k = count == 1 ? 2 : -1;
            count++;
            
            if (i >= 0 && i < n) {
                int temp = sinogram[row][i] * k;  
                ans += temp;
            }
        }
        float weightOriginal = 0.08;  // The weight of the original image
        float weightEdge = 0.92;      // Weighting of edge

        // Weighted Sum
        int combinedResult = static_cast<int>(weightOriginal * sinogram[row][col] + weightEdge * ans);
        sino[row][col] = abs(combinedResult/3);
    }
}

void sp2(){
        for(int row = 0; row < n; row++){
            for(int col = 0; col < n; col++){
                int sum = 0, cnt = 0;
                for(int i = row - 1; i <= row + 1; i++){
                    for(int j = col - 1; j <= col + 1; j++){
                        if((i >= 0 && i < n) && (j >= 0 && j < n)){
                            cnt++;
                            if(i == row && j == col){
                                int s = (int)sinogram[i][j] * 8.3;
                                sum += s;
                                //sum += sinogram[i][j] * 8;

                            } 
                            else sum += -sinogram[i][j];
                        }
                    }
                }
                int res = static_cast<int>(ceil(sum));
                sino[row][col] = abs(res);
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
    
    //for(int i = 0; i < n; i++)  sp(i);

    sp2();
    
    reduce(sino);
    Display imagesss(sino, "output6");
    /*
    for(int i = 0; i < n; i++){
        cout<<sinogram[0][i]<<" ";
    }
    cout<<endl<<endl<<endl;
    for(int i = 0; i < n; i++){
        cout<<sino[120][i]<<" ";
    }
    */

    for(int i = 0; i < n; i++) re(i, sino);


    reduce(final);

    Display image2(final, "output");

    return 0;
}