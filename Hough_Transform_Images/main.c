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


// Calculate histogram of non-zero NMS values (256 bins)
void calculateHistogram(Matrix mx, int *histogram)
{
	int i, j;
	for (i = 0; i < 256; i++) histogram[i] = 0;
	for (i = 0; i < mx.height; i++)
		for (j = 0; j < mx.width; j++)
			if (mx.map[i][j] > 0.0)
			{
				int val = (int)(mx.map[i][j]);
				if (val > 255) val = 255;
				histogram[val]++;
			}
}

// Otsu's method: find optimal threshold from histogram
double otsuThreshold(int *histogram, int numPixels)
{
	int t;
	double sum = 0.0, sumB = 0.0, wB = 0.0, wF = 0.0;
	double maxVariance = 0.0, threshold = 0.0;

	// Total sum of weighted gray levels
	for (t = 0; t < 256; t++)
		sum += t * histogram[t];

	// Find optimal threshold
	for (t = 0; t < 256; t++)
	{
		wB += histogram[t];
		if (wB == 0) continue;
		wF = numPixels - wB;
		if (wF == 0) break;

		sumB += t * histogram[t];
		double mB = sumB / wB;
		double mF = (sum - sumB) / wF;
		double variance = wB * wF * (mB - mF) * (mB - mF);

		if (variance > maxVariance)
		{
			maxVariance = variance;
			threshold = t;
		}
	}
	return threshold;
}

// Save edge matrix as PPM image (binary: 0=black, 255=white)
void saveEdgeImage(Matrix mx, const char *filename)
{
	int m, n;
	Image img = createImage(mx.height, mx.width);

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
		{
			int val = (mx.map[m][n] > 0.0) ? 255 : 0;
			img.map[m][n].r = val;
			img.map[m][n].g = val;
			img.map[m][n].b = val;
		}

	writeImage(img, (char *)filename);
	deleteImage(img);
	printf("Saved edge image: %s\n", filename);
}

// Save Hough accumulator as PPM (normalized 0-255, angle=y, distance=x)
void saveHoughSpace(Matrix mx, const char *filename)
{
	int m, n;
	double maxVal = 0.0;

	// Find max value for normalization
	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
			if (mx.map[m][n] > maxVal) maxVal = mx.map[m][n];

	if (maxVal == 0.0) maxVal = 1.0;

	Image img = createImage(mx.height, mx.width);

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
		{
			int val = (int)((mx.map[m][n] / maxVal) * 255.0);
			if (val > 255) val = 255;
			img.map[m][n].r = val;
			img.map[m][n].g = val;
			img.map[m][n].b = val;
		}

	writeImage(img, (char *)filename);
	deleteImage(img);
	printf("Saved Hough space: %s\n", filename);
}

// Export detected lines to CSV with angle, distance, strength, length
void exportLinesToCSV(Matrix maxima, Matrix houghMatrix, int imgHeight, int imgWidth, const char *filename)
{
	FILE *fp = fopen(filename, "w");
	if (!fp)
	{
		fprintf(stderr, "Warning: Could not open %s for writing\n", filename);
		return;
	}

	fprintf(fp, "angle_degrees,distance,strength,length_pixels\n");

	double maxLength = sqrt((double)(SQR(imgHeight) + SQR(imgWidth)));

	for (int i = 0; i < maxima.width; i++)
	{
		if (maxima.map[2][i] < 0.0) break;

		int row = (int)maxima.map[0][i];
		int col = (int)maxima.map[1][i];
		double strength = maxima.map[2][i];

		double alpha = -0.5*PI + PI*row/(double)houghMatrix.height;
		double dist = maxLength*col/(double)houghMatrix.width;

		// Calculate line length from image corners
		int m1 = (int)(dist*sin(alpha) - maxLength*cos(alpha) + 0.5);
		int n1 = (int)(dist*cos(alpha) + maxLength*sin(alpha) + 0.5);
		int m2 = (int)(dist*sin(alpha) + maxLength*cos(alpha) + 0.5);
		int n2 = (int)(dist*cos(alpha) - maxLength*sin(alpha) + 0.5);

		double length = sqrt((double)(SQR(m2-m1) + SQR(n2-n1)));
		double angleDeg = alpha * 180.0 / PI;

		fprintf(fp, "%.2f,%.2f,%.2f,%.2f\n", angleDeg, dist, strength, length);
	}

	fclose(fp);
	printf("Exported lines to: %s\n", filename);
}

// Read image and write Hough transform related output images.
int main()
{
	int i, m, n, m1, n1, m2, n2;
	double maxLength, alpha, dist;
	Matrix maxMatrix;
	double gaussFilter[3][3] = {{1.0, 2.0, 1.0}, {2.0, 4.0, 2.0}, {1.0, 2.0, 1.0}};
	Matrix gauss = createMatrixFromArray(&gaussFilter[0][0], 3, 3);
	Image inputImage = readImage("DW.ppm");
	Matrix inputMatrix = image2Matrix(inputImage);
	Matrix edgeMatrix = createMatrix(inputImage.height, inputImage.width);
	Matrix houghMatrix;

	// Edge detection via adaptive thresholding
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

	// Adaptive thresholding using Otsu's method
	int histogram[256];
	calculateHistogram(nms, histogram);
	int numPixels = inputImage.height * inputImage.width;
	double highThresh = otsuThreshold(histogram, numPixels);

	for (m = 0; m < inputImage.height; m++)
		for (n = 0; n < inputImage.width; n++)
			edgeMatrix.map[m][n] = (nms.map[m][n] >= highThresh) ? 255.0 : 0.0;

	printf("Otsu threshold: %.1f\n", highThresh);
	saveEdgeImage(edgeMatrix, "edges_debug.ppm");

	deleteMatrix(gaussNorm);
	deleteMatrix(smoothed);
	deleteMatrix(gradX);
	deleteMatrix(gradY);
	deleteMatrix(magnitude);
	deleteMatrix(direction);
	deleteMatrix(nms);


	houghMatrix = houghTransformLines(edgeMatrix, 360, 500);
	saveHoughSpace(houghMatrix, "hough_space.ppm");

	maxLength = sqrt((double) (SQR(inputImage.height) + SQR(inputImage.width)));

	maxMatrix = findHoughMaxima(houghMatrix, 100, 50.0);
	exportLinesToCSV(maxMatrix, houghMatrix, inputImage.height, inputImage.width, "detected_lines.csv");

	for (i = 0; i < 100; i++)
	{
		if (maxMatrix.map[2][i] < 0.0) break;
		alpha = -0.5*PI + PI*maxMatrix.map[0][i]/(double) houghMatrix.height;
		dist = maxLength*maxMatrix.map[1][i]/(double) houghMatrix.width;
		m1 = (int) (dist*sin(alpha) - maxLength*cos(alpha) + 0.5);
		n1 = (int) (dist*cos(alpha) + maxLength*sin(alpha) + 0.5);
		m2 = (int) (dist*sin(alpha) + maxLength*cos(alpha) + 0.5);
		n2 = (int) (dist*cos(alpha) - maxLength*sin(alpha) + 0.5);
		line(inputImage, m1, n1, m2, n2, 2, 18, 10, 30, 30, 30, 0);
	}
	writeImage(inputImage, "hough_lines.ppm");

	deleteMatrix(edgeMatrix);
	deleteMatrix(houghMatrix);
	deleteMatrix(maxMatrix);
	deleteImage(inputImage);
	return 0;
}
