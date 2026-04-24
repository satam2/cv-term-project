// netpbm_hough.c
// Test and demo program for the Hough transform using the netpbm.c library.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define PI 3.141592653589793


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
	//Fill Here
    

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
	int j, m, n, k, l, r, index, do_not_insert;
	double minSepSquare = SQR(minSeparation);
	double strength;
	Matrix maxima = createMatrix(3, number);

	//Fill Here
    
	return maxima;
}


// Read image "desk.ppm" and write Hough transform related output images. 
void main()
{
	int i, m, n, m1, n1, m2, n2;
	double gaussFilter[3][3] = {{1.0, 2.0, 1.0}, {2.0, 4.0, 2.0}, {1.0, 2.0, 1.0}};
	Matrix gauss = createMatrixFromArray(&gaussFilter[0][0], 3, 3); 
	Image inputImage = readImage("desk.ppm");
	Matrix inputMatrix = image2Matrix(inputImage);
	Image edgeImage, houghImage; 
	Matrix edgeMatrix = createMatrix(inputImage.height, inputImage.width);
	Matrix houghMatrix, maxMatrix;
	double maxLength, alpha, dist;

	// Add code for generating edge matrix here!!!

	edgeImage = matrix2Image(edgeMatrix, 1, 1.0);
	writeImage(edgeImage, "desk_edges.pgm");
	
	houghMatrix = houghTransformLines(edgeMatrix, 360, 500);
	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, "desk_hough.pgm");

	maxMatrix = findHoughMaxima(houghMatrix, 5, 50.0);
	for (i = 0; i < 5; i++)
		ellipse(houghImage, maxMatrix.map[0][i], maxMatrix.map[1][i], 20, 20, 2, 10, 7, 255, 255, 255, 0); 
	writeImage(houghImage, "desk_hough_max.ppm");

	maxLength = sqrt((double) (SQR(inputImage.height) + SQR(inputImage.width)));

	for (i = 0; i < 5; i++)
	{
		alpha = -0.5*PI + 1.5*PI*maxMatrix.map[0][i]/(double) houghMatrix.height;
		dist = maxLength*maxMatrix.map[1][i]/(double) houghMatrix.width;
		m1 = (int) (dist*sin(alpha) - maxLength*cos(alpha) + 0.5);
		n1 = (int) (dist*cos(alpha) + maxLength*sin(alpha) + 0.5);
		m2 = (int) (dist*sin(alpha) + maxLength*cos(alpha) + 0.5);
		n2 = (int) (dist*cos(alpha) - maxLength*sin(alpha) + 0.5);
		line(inputImage, m1, n1, m2, n2, 2, 18, 10, 30, 30, 30, 0); 
	}	writeImage(inputImage, "desk_hough_lines.ppm");

	deleteMatrix(edgeMatrix);
	deleteMatrix(houghMatrix);
	deleteImage(inputImage);
	deleteImage(edgeImage);
	deleteImage(houghImage);
}
