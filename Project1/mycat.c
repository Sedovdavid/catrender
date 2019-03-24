#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <corecrt_math.h>
#include <windows.h>

#include "tga.h"
#include "model.h"

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

void triangle(tgaImage *image, int x0, int y0, int x1, int y1, int x2, int y2, tgaColor color) {
	int dots[2][3] = { {x0, x1, x2}, {y0, y1, y2} };
	for (int i = 2; i > 0; i--) {
		for (int j = 0; j < i; j++) {
			if (dots[1][j] > dots[1][j + 1]) {
				swap(&dots[0][j], &dots[0][j + 1]);
				swap(&dots[1][j], &dots[1][j + 1]);
			}
		}
	}
	//int a = rand(255);
	//int b = rand(255);
	//int c = rand(255);
	//color = tgaRGB(a,b,c);
	int start = dots[1][0]; int stop = dots[1][1];
	for (int i = 0; i < 2; i++) {
		for (int y = start; y < stop; y++) {
			int dot0 = ((float)(y - dots[1][0]) / (dots[1][1] - dots[1][0]))*(dots[0][1] - dots[0][0]) + dots[0][0];
			int dot1 = ((float)(y - dots[1][0]) / (dots[1][2] - dots[1][0]))*(dots[0][2] - dots[0][0]) + dots[0][0];
			line(image, dot0, y, dot1, y, color);
		}
		start = dots[1][1]; stop = dots[1][2];
		swap(&dots[0][0], &dots[0][2]);
		swap(&dots[1][0], &dots[1][2]);

	}
}

void meshgrid(tgaImage *image, Model *model) {
	int i, j;
	for (i = 0; i < model->nface; ++i) {
		int screen_coords[3][2];
		int test;
		double minx = (model->vertices[0])[0]; double maxx = minx;
		double miny = (model->vertices[0])[1]; double maxy = miny;
		for (j = 0; j < model->nvert; ++j) {
			Vec3 *vert = &(model->vertices[j]);
			if ((*vert)[0] < minx) {				minx = (*vert)[0];			}
			if ((*vert)[1] < miny) {				miny = (*vert)[1];			}
			if ((*vert)[0] > maxx) {				maxx = (*vert)[0];			}
			if ((*vert)[1] > maxy) {				maxy = (*vert)[1];			}
		}
		double a; double b; double dx = 0.f; double dy = 0;
		if ((maxx - minx) > (maxy - miny)) {
			a = 2 / (maxx - minx);
			b = (minx + maxx) / (minx - maxx);
			dy = (int)(image->height / 2) - (1  - ((maxy + miny) / 2 * a + b)) * (int)(image->height / 2);
		}
		else
		{
			a = 2 / (maxy - miny);
			b = (miny + maxy) / (miny - maxy);
			dx = (int)(image->width / 2) - ((maxx + minx) / 2 * a + b + 1) * (int)(image->width / 2);
		}
		double p[3][3];
		for (j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3 * j]]);
			p[j][0] = (*v)[0];
			p[j][1] = (*v)[1];
			p[j][2] = (*v)[2];
			screen_coords[j][0] = ((*v)[0] * a + b + 1) * image->width / 2 + (int)dx;
			screen_coords[j][1] = (1 - ((*v)[1] * a + b)) * image->height / 2 + (int)dy;
		}
		
		/*Vec3 *p0 = &(model->vertices[model->faces[i][0]]);
		Vec3 *p1 = &(model->vertices[model->faces[i][3]]);
		Vec3 *p2 = &(model->vertices[model->faces[i][6]]);*/
		Vec3 v1, v2;
		for (int i = 0; i < 3; i++) {
			v1[i] = p[1][i] - p[0][i];
			v2[i] = p[2][i] - p[0][i];
		}
		Vec3 n;
		n[0] = v1[1] * v2[2] - v1[2] * v2[1];
		n[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
		n[2] = v1[0] * v2[1] - v1[1] * v2[0];
		double n_len = sqrt(pow(n[0],2) + pow(n[1],2) + pow(n[2],2)); 
		double spectator[3] = { 0.0, 0.0,-1.0 };
		/*findCenter(model, image, a, b, dx, dy);
		Vec3 ctt = { p[0][0] - point[0], p[0][1] - point[1], p[0][2] - point[2] };
		double mul = n[0] * ctt[0] + n[1] * ctt[1] + n[2] * ctt[2];
		if (mul < 0) {
			n_len = n_len;
		}*/
		double costheta = (n[0] * spectator[0]+ n[1] * spectator[1]+ n[2] * spectator[2])/n_len;
		int intensivity = -255 * costheta;
		if (intensivity > 0) {
			triangle(image, screen_coords[0][0], screen_coords[0][1], screen_coords[1][0], screen_coords[1][1],
				screen_coords[2][0], screen_coords[2][1], tgaRGB(intensivity, intensivity, intensivity));
		}

		for (j = 0; j < 3; ++j) {
			//line(image, screen_coords[j][0], screen_coords[j][1], screen_coords[(j + 1) % 3][0], screen_coords[(j + 1) % 3][1], tgaRGB(255, 255, 255));
		}
		//tgaSaveToFile(image, "C:\\Users\\sedov\\Desktop\\cat.tga");

		//»щу центр, рисую линию от (0,0) до центра, адрес center сбрасываетс€ при вызове tgaRGB 
		/*int *center = (int *)findCenter(model, image, a, b, dx, dy);
		int qwe = *(center); int rty = *(center+1);
		line(image, 0, 0, qwe, rty, tgaRGB(255, 0, 0)); */
		}
	}

void main(){
	Model *model = loadFromObj("C:\\Users\\sedov\\source\\repos\\Project1\\Project1\\obj\\tree.obj");
	
	int size = 2000;
	int height = size;
	int width = size;
	tgaImage *image = tgaNewImage(height, width, RGB);
	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			//tgaSetPixel(image, i, j, tgaRGB(0,255,0));
		}
	}

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