#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <corecrt_math.h>
#include <windows.h>

#include "tga.h"
#include "model.h"

#define SIZE 1000
#define sinphi 0,5
#define cosphi 0,8660254038
double point[3];

int *findCenter(Model *model, tgaImage *image, double a, double b, int dx, int dy) {
	float meanx; float meany; float meanz;
	meanx = 0.0;
	meany = 0.0;
	meanz = 0.0;
	for (int i = 0; i < model->nvert; i++) {
		meanx += model->vertices[i][0];
		meany += model->vertices[i][1];
		meanz += model->vertices[i][2];
	}
	meanx = meanx / model->nvert;
	meany = meany / model->nvert;
	meanz = meanz / model->nvert;
	point[0] = meanx; point[1] = meany; point[2] = meanz;
	int center[] = { 0, 0 };
	center[0] = (point[0] * a + b + 1) * image->width / 2 + dx;
	center[1] = (1 - (point[1] * a + b)) * image->height / 2 + dy;
	return center;
}

void swap(int *a, int *b) {
	int z = *a;
	*a = *b;
	*b = z;
}

int sign(int a) {
	if (a > 0) {
		return 1;
	}
	else if (a < 0) {
		return -1;
	}
	else {
		return 0;
	}
}

void line(tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color) {
	char steep;
	if(abs(y1 - y0) > abs(x1 - x0)){
		steep = 1;
		swap(&x0, &y0);
		swap(&x1, &y1);
	}
	else
	{
		steep = 0;
	}
	if (x0 > x1) {
		swap(&x0, &x1);
		swap(&y0, &y1);
	}
	int dx = x1 - x0;
	int dy = y1 - y0;
	int de = 2 * abs(dy);
	int e = 0;
	int y = y0;
	for (int x = x0; x < x1; ++x) {
		if (steep) {
			tgaSetPixel(image, y, x, color);
		}
		else
		{
			tgaSetPixel(image, x, y, color);
		}
		e += de;
		if (e > dx) {
			y += sign(dy);
			e -= 2*dx;
		}
	}
}

void triangle(tgaImage *image, int x0, int y0, double z0, int x1, int y1, double z1, int x2, int y2, double z2, tgaColor color, float zbuffer[SIZE][SIZE]) {
	int dots[3][3] = { {x0, x1, x2}, {y0, y1, y2}, {z0, z1, z2} };
	for (int i = 2; i > 0; i--) {
		for (int j = 0; j < i; j++) {
			if (dots[1][j] > dots[1][j + 1]) {
				swap(&dots[0][j], &dots[0][j + 1]);
				swap(&dots[1][j], &dots[1][j + 1]);
				swap(&dots[2][j], &dots[2][j + 1]);
			}
		}
	}
	int start = dots[1][0]; int stop = dots[1][1];
	for (int i = 0; i < 2; i++) {
		for (int y = start; y < stop; y++) {
			int dotx0 = ((float)(y - dots[1][0]) / (dots[1][1] - dots[1][0]))*(dots[0][1] - dots[0][0]) + dots[0][0];
			int dotx1 = ((float)(y - dots[1][0]) / (dots[1][2] - dots[1][0]))*(dots[0][2] - dots[0][0]) + dots[0][0];
			if (dotx0 > dotx1) {
				swap(&dotx0, &dotx1);
			}
			int dotz0 = ((float)(y - dots[1][0]) / (dots[1][1] - dots[1][0]))*(dots[2][1] - dots[2][0]) + dots[2][0];
			int dotz1 = ((float)(y - dots[1][0]) / (dots[1][2] - dots[1][0]))*(dots[2][2] - dots[2][0]) + dots[2][0];
			for (int x = dotx0; x < dotx1; x++) {
				double z;
				if ((dots[1][2] - dots[1][1]) == 0){
					z = (x - dotx0)*(dotz1 - dotz0) / (dotx1 - dotx0) + dotz0;
				}
				else {
					z = (y - start)*(dotz1 - dotz0) / (stop - start) + dotz0;
				}
				if (z > zbuffer[x][y]) {
					zbuffer[x][y] = z;
					tgaSetPixel(image, x, y, color);
				}
			}
		}
		start = dots[1][1]; stop = dots[1][2];
		swap(&dots[0][0], &dots[0][2]);
		swap(&dots[1][0], &dots[1][2]);
		swap(&dots[2][0], &dots[2][2]);

	}
}

int findIntensivity(double p[3][3], float spectator[3]) {
	Vec3 v1, v2;
	for (int i = 0; i < 3; i++) {
		v1[i] = p[1][i] - p[0][i];
		v2[i] = p[2][i] - p[0][i];
	}
	Vec3 n;
	n[0] = v1[1] * v2[2] - v1[2] * v2[1];
	n[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
	n[2] = v1[0] * v2[1] - v1[1] * v2[0];
	double n_len = sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2));
	double costheta = (n[0] * spectator[0] + n[1] * spectator[1] + n[2] * spectator[2]) / n_len;
	int intensivity = -255 * costheta;
	return intensivity;
}

void meshgrid(tgaImage *image, Model *model) {
	int i, j;
	static float zbuffer[SIZE][SIZE];
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			zbuffer[i][j] = -10000;
		}
	}
	for (i = 0; i < model->nface; ++i) {
		int test;
		double minx = (model->vertices[0])[0]; double maxx = minx;
		double miny = (model->vertices[0])[1]; double maxy = miny;
		for (j = 0; j < model->nvert; ++j) {
			Vec3 *vert = &(model->vertices[j]);
			float x = (*vert)[0] * 0.8660254038 - 0.5*(*vert)[2];
			if (x < minx) { minx = x; }
			if ((*vert)[1] < miny) { miny = (*vert)[1]; }
			if (x > maxx) { maxx = x; }
			if ((*vert)[1] > maxy) { maxy = (*vert)[1]; }
		}
		double a; double b; double dx = 0.f; double dy = 0;
		if ((maxx - minx) > (maxy - miny)) {
			a = 2 / (maxx - minx);
			b = (minx + maxx) / (minx - maxx);
			dy = (int)(image->height / 2) - (1 - ((maxy + miny) / 2 * a + b)) * (int)(image->height / 2);
		}
		else
		{
			a = 2 / (maxy - miny);
			b = (miny + maxy) / (miny - maxy);
			dx = (int)(image->width / 2) - ((maxx + minx) / 2 * a + b + 1) * (int)(image->width / 2);
		}
		int screen_coords[3][2];
		double z[3];
		double p[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);
			//static float v1[3] = { 0,0,0 };
			/*float v1; float v2; float v3;
			v1 = 0.8660254038 * ((float)((*v)[0])) - 0.5 * ((float)((*v)[2]));
			v2 = (*v)[1];
			v3 = 0.5 * (*v)[0] + 0.8660254038 * (*v)[2];*/
			/*(*v)[0] = (float)v1[0];
			(*v)[1] = (float)v1[1]; 
			(*v)[2] = (float)v1[2]; */
			p[j][0] = 0.8660254038 * ((float)((*v)[0])) - 0.5 * ((float)((*v)[2]));;
			p[j][1] = (*v)[1];
			p[j][2] = 0.5 * (*v)[0] + 0.8660254038 * (*v)[2];
			screen_coords[j][0] = (p[j][0] * a + b + 1) * image->width / 2 + (int)dx;
			screen_coords[j][1] = (1 - (p[j][1] * a + b)) * image->height / 2 + (int)dy;
			z[j] = 255* p[j][2];
		}

		float spectator[3] = { 0.0, 0.0,-1.0 };
		int intensivity = findIntensivity(p, spectator);
		if (intensivity > 0) {
			triangle(image, screen_coords[0][0], screen_coords[0][1], z[0], screen_coords[1][0], screen_coords[1][1], z[1],
			screen_coords[2][0], screen_coords[2][1], z[2], tgaRGB(intensivity, intensivity, intensivity), zbuffer);
		}

		/*for (j = 0; j < 3; ++j) {
			line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j + 1) % 3][0], screen_coords[(j + 1) % 3][1], tgaRGB(255, 255, 255));
		} */
	}
}

void main() {
	Model *model = loadFromObj("C:\\Users\\sedov\\source\\repos\\Project1\\Project1\\obj\\cat.obj");

	int height = SIZE;
	int width = SIZE;
	tgaImage *image = tgaNewImage(height, width, RGB);

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			//tgaSetPixel(image, i, j, tgaRGB(0,255,0));
		}
	}
	int a[3] = {1, 2, 3};
	a[0] = 2;

	/* line(image, 0, 0, width, 0, tgaRGB(255, 255, 255));
	line(image, 0, 0, 0, height, tgaRGB(255, 255, 255));
	line(image, 0, height-1, width, height-1, tgaRGB(255, 255, 255));
	line(image, width-1, 0, width-1, height, tgaRGB(255, 255, 255)); */
	meshgrid(image, model);

	tgaSaveToFile(image, "C:\\Users\\sedov\\Desktop\\cat.tga");
	tgaFreeImage(image);
	freeModel(model);
	return EXIT_SUCCESS;
}