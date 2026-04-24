//
//  main.c
//  CS136
//
//  Created by nha2 on 8/27/24.
// Test and demo program for netpbm. Reads a sample image and creates several output images.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netpbm.h"
#include <stdio.h>
#include <time.h>
//--------------------------------Smooting Function-----------------------------------------------------------
Matrix smoothing_filter(Matrix m1, Matrix m2)
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

//--------------------------------Canny Filter----------------------------------------------------------------
Image canny(Image img){
    
    Matrix image_matrix = image2Matrix(img);

    double gaussian[5][5] = {
    {1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256},
    {4.0/256,16.0/256,24.0/256,16.0/256, 4.0/256},
    {6.0/256,24.0/256,36.0/256,24.0/256, 6.0/256},
    {4.0/256,16.0/256,24.0/256,16.0/256, 4.0/256},
    {1.0/256, 4.0/256, 6.0/256, 4.0/256, 1.0/256}
};

    Matrix gaussianFilter = createMatrixFromArray(&gaussian[0][0], 5, 5);
    Matrix smoothFilter = smoothing_filter(image_matrix, gaussianFilter);

    double sobel_x[3][3] = {
        {-1.0, 0.0, 1.0},
        {-2.0, 0.0, 2.0},
        {-1.0, 0.0, 1.0}
    };

    double sobel_y[3][3] = {
        { 1.0, 2.0, 1.0},
        { 0.0, 0.0, 0.0},
        {-1.0,-2.0,-1.0}
    };

    Matrix sobel_x_matrix = createMatrixFromArray(&sobel_x[0][0], 3, 3);  
    Matrix sobel_y_matrix = createMatrixFromArray(&sobel_y[0][0], 3, 3);

    //convolve
    Matrix grad_x = smoothing_filter(smoothFilter, sobel_x_matrix);
    Matrix grad_y = smoothing_filter(smoothFilter, sobel_y_matrix);

    Matrix magnitude = createMatrix(img.height, img.width);
    Matrix direction = createMatrix(img.height, img.width);

    //get direction and magnitude
     for (int x = 0; x < img.width; x++){
        for (int y = 0; y < img.height; y++){
            double gx = grad_x.map[y][x];
            double gy = grad_y.map[y][x]; 

            magnitude.map[y][x] = sqrt((gx * gx) + (gy * gy));
            direction.map[y][x] = atan2(gy, gx); 
        }
    }
    
    Matrix nms = createMatrix(img.height, img.width);

    for (int x = 1; x < img.width-1; x++){
        for (int y = 1; y < img.height-1; y++){

            double angle = direction.map[y][x] * 180.0 / PI;
            if (angle < 0){
                angle += 180.0;
            }

            double current = magnitude.map[y][x];
            double neighbor1 = 0.0;
            double neighbor2 = 0.0;

            if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180)) {
                neighbor1 = magnitude.map[y][x-1];
                neighbor2 = magnitude.map[y][x+1];
            }
            else if (angle >= 22.5 && angle < 67.5) {
                neighbor1 = magnitude.map[y-1][x+1];
                neighbor2 = magnitude.map[y+1][x-1];
            }
            else if (angle >= 67.5 && angle < 112.5) {
                neighbor1 = magnitude.map[y-1][x];
                neighbor2 = magnitude.map[y+1][x];
            }
            else if (angle >= 112.5 && angle < 157.5) {
                neighbor1 = magnitude.map[y-1][x-1];
                neighbor2 = magnitude.map[y+1][x+1];
            }

            //maxima
            if (current >= neighbor1 && current >= neighbor2){
                nms.map[y][x] = current;
            }
            else{
                nms.map[y][x]= 0.0;
            }

        }
    }

    //double thresholding
    double maxMag = 0.0;
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            if (nms.map[y][x] > maxMag){
                maxMag = nms.map[y][x];
            }
        }
    }

    double highThresh = 0.4 * maxMag;
    double lowThresh = 0.2 * maxMag;

    Matrix thresholded = createMatrix(img.height, img.width);

    for (int x = 0; x < img.width; x++) {
        for (int y = 0; y < img.height; y++) {
            double value = nms.map[y][x];


            if (value >= highThresh){
                thresholded.map[y][x] = 255.0;
            }
            else{
                thresholded.map[y][x] = 0.0;
            }

                 
        }
    }

    Matrix finalEdges = createMatrix(img.height, img.width);

    for (int x = 1; x < img.width - 1; x++) {
        for (int y = 1; y < img.height - 1; y++) {

            if (thresholded.map[y][x] == 255.0) {
                finalEdges.map[y][x] = 255.0;
            }
            else if (thresholded.map[y][x] == 100.0) {
                int neighbor = 0;

                for (int j = -1; j <= 1; j++) {
                    for (int i = -1; i <= 1; i++) {
                        if (thresholded.map[y + j][x + i] == 255.0) {
                            neighbor = 1;
                        }
                    }
                }

                if (neighbor){
                    finalEdges.map[y][x] = 255.0;
                }
                    
                else{
                    finalEdges.map[y][x] = 0.0;
                }    
            }
            else {
                finalEdges.map[y][x] = 0.0;
            }
        }
    }

    Image output = matrix2Image(finalEdges, 1, 1.0);

    deleteMatrix(image_matrix);
    deleteMatrix(gaussianFilter);
    deleteMatrix(smoothFilter);

    deleteMatrix(sobel_x_matrix);
    deleteMatrix(sobel_y_matrix);
    deleteMatrix(grad_x);
    deleteMatrix(grad_y);

    deleteMatrix(magnitude);
    deleteMatrix(direction);
    deleteMatrix(nms);
    deleteMatrix(thresholded);
    deleteMatrix(finalEdges);

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

   Image canny_img = canny(inputImage);
   writeImage(canny_img, "canny.ppm");

   deleteImage(canny_img);

    return 0;

}