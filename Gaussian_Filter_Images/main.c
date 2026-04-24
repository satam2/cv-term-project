#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netpbm.h"
#include <stdio.h>
#include <time.h>

//--------------------------------Gaussian Filter-----------------------------------------------------------
Matrix gaussian_filter(Matrix m1, Matrix m2)
{

    Matrix output = createMatrix(m1.height, m1.width);

    int x_center = m2.width / 2;
    int y_center = m2.height / 2;

    for (int x = 0; x < output.width; x++)
    {
        for (int y = 0; y < output.height; y++)
        {
            double sum = 0.0;
            // apply filter
            for (int fx = 0; fx < m2.width; fx++)
            {
                for (int fy = 0; fy < m2.height; fy++)
                {

                    int pointX = x + (fx - x_center);
                    int pointY = y + (fy - y_center);

                    if (pointY >= 0 && pointY < m1.height && pointX >= 0 && pointX < m1.width)
                    {
                        sum += m1.map[pointY][pointX] * m2.map[fy][fx];
                    }
                }
            }
            output.map[y][x] = sum;
        }
    }

    return output;
}


int main(int argc, const char *argv[])
{
    Image inputImage = readImage("DW.ppm");

    for (int x = 0; x < inputImage.width; x++){
        for (int y = 0; y < inputImage.height; y++){
            int r = inputImage.map[y][x].r;
            int g = inputImage.map[y][x].g;
            int b = inputImage.map[y][x].b;

            int gray = (r + g + b) / 3;

            inputImage.map[y][x].r = gray;
            inputImage.map[y][x].g = gray;
            inputImage.map[y][x].b = gray;
        }
    }

    // create gaussian image
    Matrix inputMatrix = image2Matrix(inputImage);

    double gaussianData[5][5] = {
    {1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256},
    {4.0/256,16.0/256,24.0/256,16.0/256, 4.0/256},
    {6.0/256,24.0/256,36.0/256,24.0/256, 6.0/256},
    {4.0/256,16.0/256,24.0/256,16.0/256, 4.0/256},
    {1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256}
    };

    Matrix gaussianKernel = createMatrixFromArray(&gaussianData[0][0], 3, 3);

    Matrix gaussianOutput = gaussian_filter(inputMatrix, gaussianKernel);
    Image outputImage = matrix2Image(gaussianOutput, 1, 1.0);

    writeImage(outputImage, "gaussian.ppm");

    deleteMatrix(gaussianKernel);
    deleteMatrix(gaussianOutput);
    deleteMatrix(inputMatrix);
    deleteImage(inputImage);
    deleteImage(outputImage);

    return 0;

}