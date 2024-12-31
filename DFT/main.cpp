// Discrete Fourier transform
// Edit By: Ziyang Li
// Date: 12/31/24
// Description: The dft of a set of numbers is computed and then the set is rederived by idft 

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
using namespace std;

vector<double> real, imag;

void dft(vector<int> &nums)
{
    int n = nums.size();
    // for X[0]
    real.push_back(accumulate(nums.begin(), nums.end(), 0));
    imag.push_back(0);
    //cout << nums[0] << ":  "<< accumulate(nums.begin(), nums.end(), 0) << " + 0j" << endl;

    double pi = 3.1415926;

    for (int i = 1; i < n; i++)
    {
        double nonJ = 0.00, jPart = 0.00;
        for (int j = 0; j < n; j++)
        {
            if (j == 0)
            {
                nonJ += nums[j];
                continue;
            }
            double exponent = (2 * pi * j * i) / n;
            nonJ += nums[j] * cos(exponent);
            jPart += -(nums[j] * sin(exponent));
        }
        real.push_back(nonJ);
        imag.push_back(jPart);
    }
}

vector<double> idft()
{
    vector<double> ans;
    int n = real.size();
    double pi = 3.1415926;
    ans.resize(n);

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
        ans[k] = sumReal / n;
    }
    return ans;
}

int main()
{
    vector<int> nums = {2,2,2,2,2,2,2,2,2,2};
    //vector<int> nums2 = {22,3,92,41,36,13,42,23,3,65,6};

    dft(nums);
    //dft(nums2);
    auto ans = idft();
    for(int i = 0; i < ans.size(); i++) cout<<real[i]<<" "<<imag[i]<<" "<< ans[i]<<" ";

    // for(int i = 0; i < real.size(); i++){
    //     cout<<real[i] <<" "<<imag[i]<<endl;
    // }

    return 0;
}
