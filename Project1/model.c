#include "model.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

static ssize_t unigetline(char **lineptr, size_t *n, FILE *stream)
{
    if (!*n) *lineptr = NULL;
    if (!*lineptr) *n = 0;

    size_t i = 0;
    int ch;
    while ((ch = getc(stream)) != EOF) {
        if (ch == '\n') break;

        if (i == *n) {
            *n += 1024;
            *lineptr = realloc(*lineptr, *n);
            assert(*lineptr);
        }
        (*lineptr)[i++] = (char)ch;
    }
    (*lineptr)[i] = '\0';
    return i;
}

Model * loadFromObj(const char *filename)
{
    assert(filename);

    FILE *fd = fopen(filename, "r");
    if (!fd) {
        return NULL;
    }

    Model *model = (Model *)malloc(sizeof(Model));
    model->nvert = 0;
    model->ntext = 0;
    model->nnorm = 0;
    model->nface = 0;
    model->diffuse_map = NULL;
    model->normal_map = NULL;
    model->specular_map = NULL;

    size_t vertcap = 1;
    size_t textcap = 1;
    size_t normcap = 1;
    size_t facecap = 1;

    model->vertices = (Vec3 *)malloc(vertcap * sizeof(Vec3));
    model->textures = (Vec3 *)malloc(textcap * sizeof(Vec3));
    model->normals = (Vec3 *)malloc(normcap * sizeof(Vec3));
    model->faces = (Face *)malloc(facecap * sizeof(Face));


    char *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    while ((linelen = unigetline(&line, &linecap, fd)) > 0) {
        if (!strncmp(line, "vn", 2)) {
            if (model->nnorm >= normcap) { // realloc
                normcap *= 2;
                model->normals = (Vec3 *)realloc(model->normals, normcap * sizeof(Vec3));
                assert(model->normals);
            }
            Vec3 *vn = &model->normals[model->nnorm];
            if(3 != sscanf(line + 2, "%lg %lg %lg\n", &(*vn)[0], &(*vn)[1], &(*vn)[2])) {
                fprintf(stderr, "error while parsing: %s", line + 2);
                exit(EXIT_FAILURE);
            }
            model->nnorm += 1;
        } else if (!strncmp(line, "vt", 2)) {
            if (model->ntext >= textcap) { // realloc
                textcap *= 2;
                model->textures = (Vec3 *)realloc(model->textures, textcap * sizeof(Vec3));
                assert(model->textures);
            }
            Vec3 *vt = &model->textures[model->ntext];
            (*vt)[0] = 0.0;
            (*vt)[1] = 0.0;
            (*vt)[2] = 0.0;
            if (2 != sscanf(line + 2, "%lg %lg\n", &(*vt)[0], &(*vt)[1])) {
                fprintf(stderr, "error while parsing: %s", line + 2);
                exit(EXIT_FAILURE);
            }
            model->ntext += 1;
        } else if (!strncmp(line, "v", 1)) {
            if (model->nvert >= vertcap) { // realloc
                vertcap *= 2;
                model->vertices = (Vec3 *)realloc(model->vertices, vertcap * sizeof(Vec3));
                assert(model->vertices);
            }
            Vec3 *v = &model->vertices[model->nvert];
            if (3 != sscanf(line + 1, "%lg %lg %lg\n", &(*v)[0], &(*v)[1], &(*v)[2])) {
                fprintf(stderr, "error while parsing: %s", line + 2);
                exit(EXIT_FAILURE);                
            }
            model->nvert += 1;
        } else if (!strncmp(line, "f", 1)) {
            if (model->nface >= facecap) { // realloc
                facecap *= 2;
                model->faces = (Face *)realloc(model->faces, facecap * sizeof(Face));
                assert(model->faces);
            }
            Face *f = &model->faces[model->nface];
            assert(9 == sscanf(line + 1, "%u/%u/%u %u/%u/%u %u/%u/%u", &(*f)[0], &(*f)[1], &(*f)[2],
                                                                       &(*f)[3], &(*f)[4], &(*f)[5],
                                                                       &(*f)[6], &(*f)[7], &(*f)[8]));
            assert((*f)[0] <= model->nvert &&
                   (*f)[1] <= model->ntext &&
                   (*f)[2] <= model->nnorm);
            assert((*f)[3] <= model->nvert &&
                   (*f)[4] <= model->ntext &&
                   (*f)[5] <= model->nnorm);
            assert((*f)[6] <= model->nvert &&
                   (*f)[7] <= model->ntext &&
                   (*f)[8] <= model->nnorm);
            int i;
            for (i = 0; i < sizeof(Face)/sizeof(unsigned int); ++i) {
                (*f)[i] -= 1;
            }
            model->nface += 1;
        } else if (!strncmp(line, "#", 1) ||
                   !strncmp(line, "\n", 1)) {
            // skip comments and empty lines
        } else {
            fprintf(stderr, "Warning! Unsupported obj format: %s", line);
        }
    }

    if (line) {
        free(line);
    }
    fclose(fd);
    return model;
}

int loadDiffuseMap(Model *model, const char *filename)
{
    assert(model);
    assert(filename);
    model->diffuse_map = tgaLoadFromFile(filename);
    if (model->diffuse_map) {
        tgaFlipVertically(model->diffuse_map);
        tgaFlipHorizontally(model->diffuse_map);
    }
    return model->diffuse_map != NULL;
}

int loadNormalMap(Model *model, const char *filename)
{
    assert(model);
    assert(filename);
    model->normal_map = tgaLoadFromFile(filename);
    if (model->normal_map) {
        tgaFlipVertically(model->normal_map);
        tgaFlipHorizontally(model->normal_map);
    }
    return model->normal_map != NULL;
}

int loadSpecularMap(Model *model, const char *filename)
{
    assert(model);
    assert(filename);
    model->specular_map = tgaLoadFromFile(filename);
    if (model->specular_map) {
        tgaFlipVertically(model->specular_map);
        tgaFlipHorizontally(model->specular_map);
    }
    return model->specular_map != NULL;
}

Vec3 *getVertex(Model *model, unsigned int nface, unsigned int nvert)
{
    assert(model);
    assert(nface < model->nface);
    assert(nvert < 3);

    return &model->vertices[model->faces[nface][0 + nvert * 3]];  
}

Vec3 *getDiffuseUV(Model *model, unsigned int nface, unsigned int nvert)
{
    assert(model);
    assert(nface < model->nface);
    assert(nvert < 3);

    if (!model->textures) {
        fprintf(stderr, "uv not loaded\n");
        return NULL;
    }    

    return &model->textures[model->faces[nface][1 + nvert * 3]];
}

Vec3 * getNorm(Model *model, unsigned int nface, unsigned int nvert)
{
    assert(model);
    assert(nface < model->nface);
    assert(nvert < 3);

    if (!model->normals) {
        fprintf(stderr, "normals not loaded\n");
        return NULL;
    }    

    return &model->normals[model->faces[nface][2 + nvert * 3]];
}

tgaColor getDiffuseColor(Model *model, Vec3 *uv)
{
    assert(model);
    assert(uv);

    if (!model->diffuse_map) {
        fprintf(stderr, "Diffuse map not loaded\n");
        return tgaRGB(255, 255, 255);
    }
    unsigned int h = model->diffuse_map->height;
    unsigned int w = model->diffuse_map->width;
    return tgaGetPixel(model->diffuse_map,
                       (unsigned int)(w * (*uv)[0]),
                       (unsigned int)(h * (*uv)[1]));
}

int getNormal(Model *model, Vec3 *n, Vec3 *uv)
{
    assert(model);
    assert(uv);

    if (!model->normal_map) {
        fprintf(stderr, "Diffuse map not loaded\n");
        return -1;
    }
    unsigned int h = model->normal_map->height;
    unsigned int w = model->normal_map->width;
    tgaColor c = tgaGetPixel(model->normal_map,
                             (unsigned int)(w * (*uv)[0]),
                             (unsigned int)(h * (*uv)[1]));
    (*n)[2] = (double)Red(c)/255.0 * 2.0 - 1.0;
    (*n)[1] = (double)Green(c)/255.0 * 2.0 - 1.0;
    (*n)[0] = (double)Blue(c)/255.0 * 2.0 - 1.0;
    return 0;
}


void freeModel(Model *model)
{
    assert(model);

    if (model->vertices) 
        free(model->vertices);
    if (model->textures)
        free(model->textures);
    if (model->normals)
        free(model->normals);
    if (model->faces)
        free(model->faces);
    if (model->diffuse_map)
        tgaFreeImage(model->diffuse_map);
    if (model->normal_map)
        tgaFreeImage(model->normal_map);
    if (model->specular_map)
        tgaFreeImage(model->specular_map);
    free(model);
}
