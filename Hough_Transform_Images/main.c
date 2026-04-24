// netpbm_hough.c
// gcc -Wall -o hough main.c netpbm.c -lm 2>&1
// Test and demo program for the Hough transform using the netpbm.c library.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))

// Build a Hough parameter map for matrix mxSpatial for detecting straight lines.
// Rows in this map represent the normal alpha and columns represent the distance d from the origin.
// Increasing the size of the map in each dimension improves the resolution of the corresponding parameter.
Matrix houghTransformLines(Matrix mxSpatial, int mapHeight, int mapWidth)
{
	int m, n, angle, dist;
	double alpha, maxD = sqrt((double) (SQR(mxSpatial.height) + SQR(mxSpatial.width)));
	Matrix mxParam = createMatrix(mapHeight, mapWidth);
	Matrix sincos = createMatrix(mapHeight, 2);

	// Generate lookup table for sin and cos values to speed up subsequent computation.
	for (angle = 0; angle < mapHeight; angle++)
	{
		alpha = -0.5*PI + PI*(double)angle/(double)mapHeight;
		sincos.map[angle][0] = sin(alpha);
		sincos.map[angle][1] = cos(alpha);
	}

	for (m = 0; m < mxSpatial.height; m++)
		for (n = 0; n < mxSpatial.width; n++)
			if (mxSpatial.map[m][n] > 0.0)
				for (angle = 0; angle < mapHeight; angle++)
				{
					dist = (int) ((m*sincos.map[angle][0] + n*sincos.map[angle][1]) * (double)mapWidth / maxD + 0.5);
					if (dist >= 0 && dist < mapWidth)
						mxParam.map[angle][dist] += mxSpatial.map[m][n];
				}

	deleteMatrix(sincos);
	return mxParam; 
}


// Test whether entry (m, n) in matrix mx is a local maximum, i.e., is not exceeded by any of its 
// maximally 8 neighbors. Return 1 if true, 0 otherwise.
int isLocalMaximum(Matrix mx, int m, int n)
{
	double strength = mx.map[m][n];
	int i, j; 
	int iMin = (m == 0)? 0:(m - 1);
	int iMax = (m == mx.height -1)? m:(m + 1);
	int jMin = (n == 0)? 0:(n - 1);
	int jMax = (n == mx.width -1)? n:(n + 1);

	for (i = iMin; i <= iMax; i++)
		for (j = jMin; j <= jMax; j++)
			if (mx.map[i][j] > strength)
				return 0;
	return 1;
}


// Insert a new entry, consisting of vPos, hPos, and strength, into the list of maxima mx.
void insertMaxEntry(Matrix mx, int vPos, int hPos, double strength)
{
	int m, n = mx.width - 1;

	while (n > 0 && mx.map[2][n - 1] < strength)
	{
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n - 1];
		n--;
	}
	mx.map[0][n] = (double) vPos;
	mx.map[1][n] = (double) hPos;
	mx.map[2][n] = strength;
}


// Delete entry number i from the list of maxima mx.
void deleteMaxEntry(Matrix mx, int i)
{
	int m, n;

	for (n = i; n < mx.width - 1; n++)
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n + 1];

	mx.map[2][mx.width - 1] = -1.0;
}


// Find the <number> highest maxima in a Hough parameter map that are separated by a Euclidean distance 
// of at least <minSeparation> in the map. The result is a three-row matrix with each column representing
// the row, the column, and the strength of one maximum, in descending order of strength.
Matrix findHoughMaxima(Matrix mx, int number, double minSeparation)
{
	int m, n, k, r, index, do_not_insert;
	double minSepSquare = SQR(minSeparation);
	double strength;
	Matrix maxima = createMatrix(3, number);

	for (m = 0; m < number; m++)
		maxima.map[2][m] = -1.0;

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
		{
			strength = mx.map[m][n];
			if (strength <= maxima.map[2][maxima.width - 1])
				continue;
			if (!isLocalMaximum(mx, m, n))
				continue;

			do_not_insert = 0;
			for (k = 0; k < maxima.width; k++)
			{
				if (maxima.map[2][k] < 0.0)
					break;
				r = (int) maxima.map[0][k];
				index = (int) maxima.map[1][k];
				if (SQR(m - r) + SQR(n - index) < minSepSquare)
				{
					if (strength > maxima.map[2][k])
					{
						deleteMaxEntry(maxima, k);
					}
					else
					{
						do_not_insert = 1;
						break;
					}
				}
			}

			if (!do_not_insert)
				insertMaxEntry(maxima, m, n, strength);
		}

	return maxima;
}


// Read image and write Hough transform related output images. 
int main()
{
	int m, n, m1, n1, m2, n2;
	double maxLength, alpha, dist;
	double gaussFilter[3][3] = {{1.0, 2.0, 1.0}, {2.0, 4.0, 2.0}, {1.0, 2.0, 1.0}};
	Matrix gauss = createMatrixFromArray(&gaussFilter[0][0], 3, 3);
	Image inputImage = readImage("DW.ppm");
	Matrix inputMatrix = image2Matrix(inputImage);
	Matrix edgeMatrix = createMatrix(inputImage.height, inputImage.width);
	Matrix houghMatrix;

	// Add code for generating edge matrix here!!!
	// Normalize the Gaussian filter (sum = 16)
	Matrix gaussNorm = createMatrix(3, 3);
	int gi, gj;
	for (gi = 0; gi < 3; gi++)
		for (gj = 0; gj < 3; gj++)
			gaussNorm.map[gi][gj] = gauss.map[gi][gj] / 16.0;

	// Smooth the input
	Matrix smoothed = createMatrix(inputImage.height, inputImage.width);
	int gcx = 1, gcy = 1;
	for (m = 0; m < inputImage.height; m++)
		for (n = 0; n < inputImage.width; n++)
		{
			double sum = 0.0;
			for (gi = 0; gi < 3; gi++)
				for (gj = 0; gj < 3; gj++)
				{
					int pm = m + gi - gcy, pn = n + gj - gcx;
					if (pm >= 0 && pm < inputImage.height && pn >= 0 && pn < inputImage.width)
						sum += inputMatrix.map[pm][pn] * gaussNorm.map[gi][gj];
				}
			smoothed.map[m][n] = sum;
		}

	// Sobel gradients
	double sobelX[3][3] = {{-1.0, 0.0, 1.0}, {-2.0, 0.0, 2.0}, {-1.0, 0.0, 1.0}};
	double sobelY[3][3] = {{ 1.0, 2.0, 1.0}, { 0.0, 0.0, 0.0}, {-1.0,-2.0,-1.0}};
	Matrix gradX = createMatrix(inputImage.height, inputImage.width);
	Matrix gradY = createMatrix(inputImage.height, inputImage.width);
	for (m = 1; m < inputImage.height - 1; m++)
		for (n = 1; n < inputImage.width - 1; n++)
		{
			double gx = 0.0, gy = 0.0;
			for (gi = 0; gi < 3; gi++)
				for (gj = 0; gj < 3; gj++)
				{
					gx += smoothed.map[m + gi - 1][n + gj - 1] * sobelX[gi][gj];
					gy += smoothed.map[m + gi - 1][n + gj - 1] * sobelY[gi][gj];
				}
			gradX.map[m][n] = gx;
			gradY.map[m][n] = gy;
		}

	// Magnitude and direction
	Matrix magnitude = createMatrix(inputImage.height, inputImage.width);
	Matrix direction = createMatrix(inputImage.height, inputImage.width);
	for (m = 0; m < inputImage.height; m++)
		for (n = 0; n < inputImage.width; n++)
		{
			double gx = gradX.map[m][n], gy = gradY.map[m][n];
			magnitude.map[m][n] = sqrt(gx*gx + gy*gy);
			direction.map[m][n] = atan2(gy, gx);
		}

	// Non-maximum suppression
	Matrix nms = createMatrix(inputImage.height, inputImage.width);
	for (m = 1; m < inputImage.height - 1; m++)
		for (n = 1; n < inputImage.width - 1; n++)
		{
			double angle = direction.map[m][n] * 180.0 / PI;
			if (angle < 0) angle += 180.0;
			double cur = magnitude.map[m][n];
			double nb1 = 0.0, nb2 = 0.0;
			if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180)) {
				nb1 = magnitude.map[m][n-1]; nb2 = magnitude.map[m][n+1];
			} else if (angle >= 22.5 && angle < 67.5) {
				nb1 = magnitude.map[m-1][n+1]; nb2 = magnitude.map[m+1][n-1];
			} else if (angle >= 67.5 && angle < 112.5) {
				nb1 = magnitude.map[m-1][n]; nb2 = magnitude.map[m+1][n];
			} else {
				nb1 = magnitude.map[m-1][n-1]; nb2 = magnitude.map[m+1][n+1];
			}
			nms.map[m][n] = (cur >= nb1 && cur >= nb2) ? cur : 0.0;
		}

	// Double thresholding
	double maxMag = 0.0;
	for (m = 0; m < inputImage.height; m++)
		for (n = 0; n < inputImage.width; n++)
			if (nms.map[m][n] > maxMag) maxMag = nms.map[m][n];
	double highThresh = 0.4 * maxMag;

	for (m = 0; m < inputImage.height; m++)
		for (n = 0; n < inputImage.width; n++)
			edgeMatrix.map[m][n] = (nms.map[m][n] >= highThresh) ? 255.0 : 0.0;

	deleteMatrix(gaussNorm);
	deleteMatrix(smoothed);
	deleteMatrix(gradX);
	deleteMatrix(gradY);
	deleteMatrix(magnitude);
	deleteMatrix(direction);
	deleteMatrix(nms);


	houghMatrix = houghTransformLines(edgeMatrix, 360, 500);

	maxLength = sqrt((double) (SQR(inputImage.height) + SQR(inputImage.width)));

	double houghMax = 0.0;
	for (m = 0; m < houghMatrix.height; m++)
		for (n = 0; n < houghMatrix.width; n++)
			if (houghMatrix.map[m][n] > houghMax) houghMax = houghMatrix.map[m][n];
	double houghThresh = 0.5 * houghMax;

	for (m = 0; m < houghMatrix.height; m++)
		for (n = 0; n < houghMatrix.width; n++)
			if (houghMatrix.map[m][n] >= houghThresh && isLocalMaximum(houghMatrix, m, n))
			{
				alpha = -0.5*PI + PI*(double)m/(double) houghMatrix.height;
				dist = maxLength*(double)n/(double) houghMatrix.width;
				m1 = (int) (dist*sin(alpha) - maxLength*cos(alpha) + 0.5);
				n1 = (int) (dist*cos(alpha) + maxLength*sin(alpha) + 0.5);
				m2 = (int) (dist*sin(alpha) + maxLength*cos(alpha) + 0.5);
				n2 = (int) (dist*cos(alpha) - maxLength*sin(alpha) + 0.5);
				line(inputImage, m1, n1, m2, n2, 2, 18, 10, 30, 30, 30, 0);
			}
	writeImage(inputImage, "hough_lines.ppm");

	deleteMatrix(edgeMatrix);
	deleteMatrix(houghMatrix);
	deleteImage(inputImage);
	return 0;
}
