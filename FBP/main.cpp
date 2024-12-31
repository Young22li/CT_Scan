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

void re(int angle, int (&s)[n][n])
{
    int** temp = new int*[n];
    for (int i = 0; i < n; ++i) temp[i] = new int[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp[j][i] = s[angle][i];
        }
        rotation(-angle, temp);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            final[i][j] += tracy[i][j];
        }
    }
    for (int i = 0; i < n; ++i) delete[] temp[i];

    delete[] temp;
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

void gray(int (&s)[n][n])
{
    int grayValue = 148; 
    int threshold = 15;   

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (s[i][j] <= threshold)
            {
                
                s[i][j] = grayValue;
            }    
            
        }
    }
}

void gray2(int (&s)[n][n]) {
    double alpha = 1.2;
    int beta = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            
            s[i][j] = static_cast<int>(alpha * s[i][j] + beta);

            
            if (s[i][j] >= 30 && s[i][j] < 60) s[i][j] += 0;
            else if (s[i][j] >= 60 && s[i][j] < 100) s[i][j] += 55;
            else if (s[i][j] >= 100 && s[i][j] < 120) s[i][j] += 30;
        }
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

void createFilter2(double (&rampFilter)[n])
{
    double c;
    double beta = 1.5;
    double I0_beta = std::cyl_bessel_i(0, beta);
    for (int i = 0; i < n; i++)
    {
        // 
        c = abs(i - n / 2);
        double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
        //double window = 0.42 - 0.5 * cos(2 * M_PI * i / (n - 1)) + 0.08 * cos(4 * M_PI * i / (n - 1));
        //double window = 0.3 * (1 - cos(2 * M_PI * i / (n - 1)));
        double x = 2.0 * i / (n - 1) - 1.0; // 归一化索引 [-1, 1]
        double kaiser_window = std::cyl_bessel_i(0, beta * sqrt(1 - x * x)) / I0_beta;

        //
        //rampFilter[i] = (5.0 * c / n) * window; 
        rampFilter[i] = (5.0 * c / n) * kaiser_window;
        //rampFilter[i] = (5.0 * c / n);

        
    }
}

void createFilter3(double (&rampFilter)[n])
{
     double c; 
    
    for (int i = 0; i < n; ++i) {
        c = abs(i - n / 2); 

        // Hanning 
        double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));

        // Ram-Lak 
        if (c < n / 4){
            //rampFilter[i] = (2.0 * c / n); 
           rampFilter[i] = (7.0 * c / n) * window; 
        } else {
            rampFilter[i] = 1.0; 
        }
    }

    // 
    double maxValue = 0.0;
    for(int i = 0; i < n; i++) maxValue = max(rampFilter[i],maxValue);
    if (maxValue > 0) {
        for (int i = 0; i < n; ++i) {
            rampFilter[i] /= maxValue; // 
        }
    }
}

void createFilter4(double (&rampFilter)[n]) {
    double fc = 500;  // 截止频率，默认值为 0.5
    int butterworthOrder = 3; // Butterworth 滤波器阶数，默认值为 3
    double beta = 5.0;
    double I0_beta = std::cyl_bessel_i(0, beta);

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
        //double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
        //double window = 0.42 - 0.5 * cos(2 * M_PI * i / (n - 1)) + 0.08 * cos(4 * M_PI * i / (n - 1));
        //double window = 0.3 * (1 - cos(2 * M_PI * i / (n - 1)));
        double x = 2.0 * i / (n - 1) - 1.0; // 归一化索引 [-1, 1]
        double kaiser_window = std::cyl_bessel_i(0, beta * sqrt(1 - x * x)) / I0_beta;

        rampFilter[i] = sinc / butterworth * kaiser_window;
        //rampFilter[i] = sinc;
        //rampFilter[i] = butterworthOrder;
        //rampFilter[i] = -1.1;
    }
}

void createFilter5(double (&rampFilter)[n])
{
    double c;
    double beta = 0.5;
    double I0_beta = std::cyl_bessel_i(0, beta);
    
    for (int i = 0; i < n; i++)
    {
        // 计算 Ramp 滤波器部分
        c = abs(i - n / 2);
        double ramp = 5.0 * c / n;

        //double window = 0.54 - 0.46 * cos(2 * M_PI * i / (n - 1));
         double x = 2.0 * i / (n - 1) - 1.0; // 归一化索引 [-1, 1]
        double kaiser_window = std::cyl_bessel_i(0, beta * sqrt(1 - x * x)) / I0_beta;

        // 计算归一化频率 u
        //double u = static_cast<double>(i / n);
        //double u = static_cast<double>(i - n / 2) / n;

        // 计算 Sinc 函数部分，避免除以零
        double sinc = c == 0 ? 0.0 : sin((M_PI * c)/n) / ((M_PI * c)/n);

        // 组合 Ramp 和 Sinc 滤波器
        rampFilter[i] = ramp * sinc * kaiser_window;
       //rampFilter[i] = 2 * sinc;
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
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

// Gaussian
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

void backProject(int (&sinogram)[n][n]) {
    double PI = 3.1415926;
    int sinogramRows = 180; // Sinogram 的行数
    int sinogramCols = 400; // Sinogram 的列数（宽度）

    // Backprojection loop: iterate over angles (0° to 180° counterclockwise)
    for (int angleIndex = 0; angleIndex < sinogramRows; angleIndex++) {
        // 转换角度到弧度，负数表示顺时针，+ PI / 2 矫正
        double angle = -(angleIndex * PI / sinogramRows) + PI / 2;
        double cosTheta = cos(angle);
        double sinTheta = sin(angle);

        // 遍历重建图片的每个像素点
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                // 计算 sinogram 中对应的位置 t
                int t = static_cast<int>(round((x - n / 2) * cosTheta + (y - n / 2) * sinTheta)) + sinogramCols / 2;

                // 如果 t 在 sinogram 范围内，则累加
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

    Display image(Pixels, "output1");
    /*

     for (int i = 0; i < n; i++)
    {
        vector<double> real, imag;
        dft(i, real, imag, Pixels);
        double rampFilter[n];
        createFilter5(rampFilter);
        
        /*
        for(int j = 0; j < n; j++){
            real[j] *= rampFilter[j];
            imag[j] *= rampFilter[j];
        }
        */

        /*
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

        idft(i, real, imag, Pixels);
    }

    reduce(sino);

    Display image1(sino, "output2");
    */

    //sp2Gaussian();

    //gray(sino);
    //gray2(sino);

    //Display image3(sino, "output4");

    for (int i = 0; i < n; i++)
        re(i, Pixels); 

   // backProject(sino);
   //backProject(Pixels);
    reduce(final);
    // black();

    Display image2(final, "output3");

	return 0;
}