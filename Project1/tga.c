#include "tga.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#pragma pack(push, 1)
struct tgaHeader {
    unsigned char id_len;
    unsigned char color_map_type;
    unsigned char image_type;
    unsigned short color_map_idx;
    unsigned short color_map_len;
    unsigned char color_map_bpp;
    unsigned short image_pos_x;
    unsigned short image_pos_y;
    unsigned short image_width;
    unsigned short image_height;
    unsigned char image_bpp;
    unsigned char image_descriptor;
};
#pragma pack(pop)

tgaColor tgaRGB(unsigned char r, unsigned char g, unsigned char b)
{
    return 0 | (r << 16) | (g << 8) | (b << 0);
}

unsigned char Red(tgaColor c)
{
    return (unsigned char)(c >> 16);
}

unsigned char Blue(tgaColor c)
{
    return (unsigned char)(c >> 0);
}

unsigned char Green(tgaColor c)
{
    return (unsigned char)(c >> 8);
}

static int loadRLE(tgaImage *, FILE *);

tgaImage * tgaNewImage(unsigned int height, unsigned int width, int format)
{
    assert(height && width); /* both must be greater then zero */

    tgaImage *image = (tgaImage *)malloc(sizeof(tgaImage));
    if (!image) {
        return NULL;
    }

    image->height = height;
    image->width = width;
    image->bpp = (unsigned char)format;

    unsigned int data_size = image->height * image->width * image->bpp;
    image->data = (unsigned char *)malloc(data_size);
    if (!image->data) {
        free(image);
        return NULL;
    }

    memset(image->data, 0, data_size); /* TODO: explicit tgaColor? */

    return image;
}

void tgaFreeImage(tgaImage *image)
{
    assert(image);
    assert(image->data);

    free(image->data);
    free(image);
}

int tgaSetPixel(tgaImage *image, unsigned int x, unsigned int y, tgaColor c)
{
    assert(image);
    if (x >= image->width ||
          y >= image->height)
        return -1;

    unsigned char *pixel_pos = (unsigned char *)(image->data + (x + y * image->width) * image->bpp);
    char *b = (char*)&c;

    switch (image->bpp) {
        case RGBA:
            *pixel_pos++ = *b++;
        case RGB:
            *pixel_pos++ = *b++;
            *pixel_pos++ = *b++;
        case GRAYSCALE:
            *pixel_pos++ = *b++;
    }
    // memcpy(pixel_pos, &c, image->bpp);
    return 0;
}

tgaColor tgaGetPixel(tgaImage *image, unsigned int x, unsigned int y)
{
    assert(image);
    assert(x < image->width);
    assert(y < image->height);

    unsigned char *pixel_pos = image->data + (x + y * image->width) * image->bpp;

    tgaColor color;
    memcpy(&color, pixel_pos, image->bpp);
    return color; 
}

int tgaSaveToFile(tgaImage *image, const char *filename)
{
    assert(image);
    assert(filename);

    FILE *fd = fopen(filename, "wb");
    if (!fd) {
        return -1;
    }
    char new_tga_format_signature[] = {'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'};
    unsigned char extension_offset[] = { 0, 0, 0, 0};
    unsigned char developer_offset[] = { 0, 0, 0, 0};

    struct tgaHeader header;
    header.id_len = 0;
    header.color_map_type = 0; /* without ColorMap */ 
    header.image_type = (image->bpp == GRAYSCALE) ? 3 : 2;
    header.color_map_idx = 0;
    header.color_map_len = 0;
    header.color_map_bpp = 0;
    header.image_pos_x = 0;
    header.image_pos_y = 0;
    header.image_width = (unsigned short)image->width;
    header.image_height = (unsigned short)image->height;
    header.image_bpp = image->bpp << 3;
    header.image_descriptor = 0x20; /* top-left origin */

    int rv = 0;
    do {
        if (-1 == fwrite(&header, sizeof(header), 1, fd)) {
            rv = -1;
            break;
        }

        unsigned int data_size = image->height * image-> width * image->bpp;
        if (-1 == fwrite(image->data, data_size, 1, fd)) {
            rv = -1;
            break;
        }

        if (-1 == fwrite(extension_offset, sizeof(extension_offset), 1, fd)) {
            rv = -1;
            break;
        }

        if (-1 == fwrite(developer_offset, sizeof(developer_offset), 1, fd)) {
            rv = -1;
            break;
        }

        if (-1 == fwrite(new_tga_format_signature, sizeof(new_tga_format_signature), 1, fd)) {
            rv = -1;
            break;
        }

    } while (0);

    fclose(fd);
    return rv;;
}

tgaImage * tgaLoadFromFile(const char *filename)
{
    assert(filename);
    FILE *fd = fopen(filename, "r");
    if (!fd) {
        return NULL;
    }

    struct tgaHeader header;
    if (1 != fread(&header, sizeof(header), 1, fd)) {
        fclose(fd);
        return NULL;
    }

    if (header.color_map_type == 1) {
        fprintf(stderr, "TGA files with color map unsupported\n");
        fclose(fd);
        return NULL;
    }
    tgaImage *image = tgaNewImage(header.image_height,
                                  header.image_width,
                                  header.image_bpp >> 3);
    if (!image) {
        fclose(fd);
        return NULL;
    }
    if (header.image_type == 3 || header.image_type == 2) {
        unsigned int size = image->height * image->width * image->bpp;
        if (!fread(image->data, size, 1, fd)) {
            tgaFreeImage(image);
            fclose(fd);
            return NULL;
        }
    } else if (header.image_type == 11 || header.image_type == 10) {
        if (-1 == loadRLE(image, fd)) {
            tgaFreeImage(image);
            fclose(fd);
            return NULL;
        }

    } else {
        fprintf(stderr, "Unknown image type: %u\n", header.image_type);
        tgaFreeImage(image);
    }

    if (!(header.image_descriptor & 0x20)) {
        tgaFlipVertically(image);
    }
    if (!(header.image_descriptor & 0x10)) {
        tgaFlipHorizontally(image);
    }

    fclose(fd);
    return image;
}

void tgaFlipVertically(tgaImage *image)
{
    assert(image);
    unsigned int bytes_per_line = image->width * image->bpp;
    unsigned char *line = (unsigned char *)malloc(bytes_per_line);
    assert(line);
    unsigned int half = image->height / 2;
    unsigned int j;
    for (j = 0; j < half; ++j) {
        // swap lines
        unsigned int l1 = j * bytes_per_line;
        unsigned int l2 = (image->height - 1 - j) * bytes_per_line;
        memcpy(line, image->data + l1, bytes_per_line);
        memmove(image->data + l1, image->data + l2, bytes_per_line);
        memcpy(image->data + l2, line, bytes_per_line);
    }
    free(line);
}

void tgaFlipHorizontally(tgaImage *image)
{
    assert(image);
    unsigned int half = image->width / 2;
    unsigned int i, j;
    for (i = 0; i < half; ++i) {
        for (j = 0; j < image->height; ++j) {
            tgaColor c1 = tgaGetPixel(image, i, j);
            tgaColor c2 = tgaGetPixel(image, image->width - 1 - i, j);
            tgaSetPixel(image, i, j, c2);
            tgaSetPixel(image, image->width - 1 - i , j, c1);
        }
    }
}

int loadRLE(tgaImage *image, FILE *stream)
{
    unsigned int size = image->height * image->width * image->bpp;
    unsigned char *write_buf = image->data;
    while (write_buf < image->data + size) {
        unsigned char chunk_size = 0;
        if (!fread(&chunk_size, sizeof(chunk_size), 1, stream)) {
            return -1;
        }
        if (chunk_size < 128) {
            ++chunk_size;
            if (write_buf + chunk_size * image->bpp > image->data + size) {
                fprintf(stderr, "Chunk size is greater then image data\n");
                return -1;
            }
            if (!fread(write_buf, chunk_size * image->bpp, 1, stream)) {
                return -1;
            }         
            write_buf += chunk_size * image->bpp;
        } else {
            chunk_size -= 127;
            if (write_buf + chunk_size * image->bpp > image->data + size) {
                fprintf(stderr, "Chunk size is greater then image data\n");
                return -1;
            }
            tgaColor color;
            if (!fread(&color, image->bpp, 1, stream)) {
                return -1;
            }
            int i;
            for (i = 0; i < chunk_size; ++i) {
                memcpy(write_buf, &color, image->bpp);
                write_buf += image->bpp;
            }
        }
    }
    return 0;
}



