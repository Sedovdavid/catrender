#include "complex.h"
#include <stdio.h>
#include <stdlib.h>

#include "tga.h"

int calculate_z(int maxiter, Complex z, Complex c) {
  while (maxiter-- && (cabs_(z) < 2)) {
    z = cadd(cmul(z, z), c);
  }
  return maxiter;
}

void draw_julia_set(tgaImage *image, double x1, double x2, double y1, double y2,
                    Complex c, int maxiter) {
  double x_step = (double)(x2 - x1) / image->width;
  double y_step = (double)(y2 - y1) / image->height;

  unsigned i, j;
  double x, y;
  for (i = 0, x = x1; i < image->width; ++i, x += x_step) {
    for (j = 0, y = y1; j < image->height; ++j, y += y_step) {
      Complex z0 = {x, y};
      double intensity =
          (1.0 -
           (double)calculate_z(maxiter, z0, c) / maxiter);
      tgaColor color = tgaRGB((unsigned char)(255 * intensity),
                              (unsigned char)(255 * intensity),
                              0);
      tgaSetPixel(image, i, j, color);
    }
  }
}

int main() {
  /*if (argc < 2) {
    fprintf(stderr, "Usage: %s outfile\n", argv[0]);
    return EXIT_FAILURE;
  }*/

  double x1 = -1.8, x2 = 1.8, y1 = -1.8, y2 = 1.8;
  Complex c = {-0.4, 0.6};

  int width = 800;
  int height = 800;
  int max_iterations = 300;

  tgaImage *image = tgaNewImage(width, height, RGB);
  if (image == NULL) {
    fprintf(stderr, "%s\n", "Failed to create image");
    return EXIT_FAILURE;
  }

  draw_julia_set(image, x1, x2, y1, y2, c, max_iterations);

  if (-1 == tgaSaveToFile(image, "C:\\Users\\sedov\\Desktop\\image1.tga")) {
    perror("tgaSateToFile");
    return EXIT_FAILURE;
  }

  tgaFreeImage(image);
  return EXIT_SUCCESS;
}
