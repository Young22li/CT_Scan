#include <iostream>
#include "Display.h"
#include <cmath>
#include <algorithm>
#include "line_functions.h"
#include <numeric>
using namespace std;

const int n = 180; // Image Size N x N
int Pixels[n][n];

// Save the image after rotation
int tracy[n][n];

int final[n][n] = {0};

int sinogram[n][n];

int sino[n][n];
int si[n][n];

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

void rotation2(int angle, int (&Pixels)[n][n]){
    double PI = 3.1415926;
    double rad = angle * PI / 180.0; 
    int centerX = n / 2;
    int centerY = n / 2;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tracy[i][j] = 0;
        }
    }

    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {

            int transX = x - centerX;
            int transY = y - centerY;
            int originalX = static_cast<int>(cos(-rad) * transX - sin(-rad) * transY + centerX);
            int originalY = static_cast<int>(sin(-rad) * transX + cos(-rad) * transY + centerY);

            if (originalX >= 0 && originalX < n && originalY >= 0 && originalY < n) {
                tracy[y][x] = Pixels[originalY][originalX];
            }
        }
    }

}

void re(int angle, int (&s)[n][n])
{
    int temp[n][n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp[j][i] = s[angle][i];
        }
        rotation2(-angle, temp);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            final[i][j] += tracy[i][j];
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


void gray2(int (&s)[n][n]) {
    double alpha = 1.5;
    int beta = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            s[i][j] = static_cast<int>(alpha * s[i][j] + beta);

            
             if (s[i][j] >= 0 && s[i][j] < 30) {
                s[i][j] += 0;
            }
            else if (s[i][j] >= 30 && s[i][j] < 60) s[i][j] += 35;
            else if (s[i][j] >= 60 && s[i][j] < 100) s[i][j] += 65;
            else if (s[i][j] >= 100 && s[i][j] < 120) s[i][j] += 30;
        }
    }
}

void dft(int row, vector<double> &real, vector<double> &imag)
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

void idft(int row, vector<double> &real, vector<double> &imag)
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
    for (int i = 0; i < n; i++)
    {
        c = abs(i - n / 2);
        rampFilter[i] = (n - 2 * c) / n;
    }
}

void createFilter2(double (&rampFilter)[n])
{
    double c;
    for (int i = 0; i < n; i++)
    {
        // 
        c = abs(i - n / 2);
        double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
        //double window = 0.42 - 0.5 * cos(2 * M_PI * i / (n - 1)) + 0.08 * cos(4 * M_PI * i / (n - 1));
        //double window = 0.3 * (1 - cos(2 * M_PI * i / (n - 1)));

        //
        rampFilter[i] = (2.0 * c / n) * window; // 正U型公式

        //rampFilter[i] = n - 2 * c;
    }
}

void generateGaussianKernel(vector<vector<double>>& kernel, double sigma, int size) {
    int halfSize = size / 2;
    double sum = 0.0;
    for (int i = -halfSize; i <= halfSize; ++i) {
        for (int j = -halfSize; j <= halfSize; ++j) {
            kernel[i + halfSize][j + halfSize] = exp(-(i * i + j * j) / (2 * sigma * sigma));
            sum += kernel[i + halfSize][j + halfSize];
        }
    }
    // 归一化核
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

// 高斯平滑滤波
void sp2Gaussian() {
    double sigma = 0.1;
    int size = 3;
    vector<vector<double>> kernel(size, vector<double>(size, 0.0));
    generateGaussianKernel(kernel, sigma, size);

    int halfSize = size / 2;

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            double sum = 0.0;
            for (int i = -halfSize; i <= halfSize; i++) {
                for (int j = -halfSize; j <= halfSize; j++) {
                    int x = row + i;
                    int y = col + j;
                    if (x >= 0 && x < n && y >= 0 && y < n) {
                        sum += sino[x][y] * kernel[i + halfSize][j + halfSize];
                    }
                }
            }
            si[row][col] = static_cast<int>(ceil(sum));
        }
    }
}

void sp2() {
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            std::vector<int> windowValues;

            // 遍历窗口中的像素
            for (int i = row - 1; i <= row + 1; i++) {
                for (int j = col - 1; j <= col + 1; j++) {
                    // 检查是否在有效范围内
                    if (i >= 0 && i < n && j >= 0 && j < n) {
                        windowValues.push_back(sinogram[i][j]);
                    }
                }
            }

            // 对窗口中的像素值排序
            int sum = accumulate(windowValues.begin(), windowValues.end(), 0);

            // 取中值
            int median = sum / 9;

            // 设置中值到输出图像
            sino[row][col] = median;
        }
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

void createFilter4(double (&rampFilter)[n]) {
    double fc = 350;  // 截止频率，默认值为 0.5
    int butterworthOrder = 3; // Butterworth 滤波器阶数，默认值为 3

    for (int i = 0; i < n; i++)
    {
        // 计算归一化频率 u
        //double u = static_cast<double>(i - n / 2) / n;
        //double u = i / n;
        double c = abs(i - n / 2);

        // 计算 Ramp 部分
        //double ramp = 2.0 * abs(i - n / 2) / n;

        // 计算 Sinc 部分，避免 u = 0 时除以零
        double sinc = sin((M_PI * c)/n);

        // 计算 Butterworth 部分
        double butterworth = M_PI * sqrt(1.0 + pow(c / (2 * fc), 2 * butterworthOrder));

        // 组合滤波器
        rampFilter[i] = sinc / butterworth;
        //rampFilter[i] = sinc;
        //rampFilter[i] = 1.0;
    }
}

void createFilter5(double (&rampFilter)[n])
{
    double c;
    double fc = 350; // 截止频率，默认为 0.5
    for (int i = 0; i < n; i++)
    {
        // 计算 Ramp 滤波器部分
        c = abs(i - n / 2);
        double ramp = 2.0 * c / n;

        double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));

        // 计算归一化频率 u
        double u = static_cast<double>(i / n);
        //double u = static_cast<double>(i - n / 2) / n;

        // 计算 Sinc 函数部分，避免除以零
        double sinc = c == 0 ? 1.0 : sin((M_PI * c)/n) / ((M_PI * c)/n);

        // 组合 Ramp 和 Sinc 滤波器
        rampFilter[i] = ramp * sinc * window;
       //rampFilter[i] = 2 * sinc;
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

    // Display image(Pixels, "output1");

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

    Display image(sinogram, "output10");

    for (int i = 0; i < n; i++)
    {
        vector<double> real, imag;
        dft(i, real, imag);
        double rampFilter[n];
        createFilter4(rampFilter);
        
        
        // for(int j = 0; j < n; j++){
        //     real[j] *= rampFilter[j];
        //     imag[j] *= rampFilter[j];
        // }
        
        
        


      // vector<double> temp1 = real;
       //vector<double> temp2 = imag;
        
        int mid = n / 2;
        for (int j = 0; j < mid; j++)
        {
             real[mid + j] *= rampFilter[j];
             imag[mid + j] *= rampFilter[j];
        }
        for (int j = 0; j < mid; j++)
        {
            real[j] *= rampFilter[mid + j];
            imag[j] *= rampFilter[mid + j];
        }
        
        
        
        

        idft(i, real, imag);
    }

    reduce(sino);

    Display image1(sino, "output11");

    //sp2Gaussian();
    //sp2();

    //gray(sino);
    //gray2(sino);

    //Display image3(si, "output13");

    //for (int i = 0; i < n; i++)
    //    re(i, sino);

    backProject(sino);
    reduce(final);
    // black();

    Display image2(final, "output12");

    return 0;
}