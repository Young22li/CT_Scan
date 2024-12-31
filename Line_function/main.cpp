// Creat image
// Edit By: Ziyang Li
// Date: 10/10/24
// Description: This file will draw three graphs using the line algorithm
// Then using display.h and html.dat, display as a HTML file.

#include <iostream>
#include "Display.h"
#include "line_functions.h"
using namespace std;

const int n = 180; // Image Size N x N
int Pixels[n][n];

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
    return 0;
}
