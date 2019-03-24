#ifndef MODEL_H_
#define MODEL_H_

#include "tga.h"

typedef unsigned int Face[9];
typedef double Vec3[3];

typedef struct Model {
    unsigned int nvert; // number of vertices
    unsigned int ntext; // number of texture coords
    unsigned int nnorm; // number of normals
    unsigned int nface; // number of faces
    Vec3 *vertices;
    Vec3 *textures;
    Vec3 *normals;
    Face *faces;
    tgaImage *diffuse_map;
    tgaImage *normal_map;
    tgaImage *specular_map;
} Model;

Model * loadFromObj(const char *filename);

int loadDiffuseMap(Model *model, const char *filename);
int loadNormalMap(Model *model, const char *filename);
int loadSpecularMap(Model *model, const char *filename);

Vec3 * getDiffuseUV(Model *model, unsigned int nface, unsigned int nvert);

Vec3 * getVertex(Model *model, unsigned int nface, unsigned int nvert);

Vec3 * getNorm(Model *model, unsigned int nface, unsigned int nvert);

tgaColor getDiffuseColor(Model *model, Vec3 *uv);

int getNormal(Model *model, Vec3 *n, Vec3 *uv);

void freeModel(Model *);

#endif // MODEL_H_
