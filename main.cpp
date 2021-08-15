#define MC_IMPLEM_ENABLE
#include "SETTINGS.h"

#include "MC.h"
#include "mesh.h"
#include "field.h"

#include <iostream>
#include <fstream>
#include <chrono>

Real sphereFunction(VEC3F pos) {
    return (pos.norm() - 10);
}

Real cubeFunction(VEC3F pos) {
    for (int i = 0; i < 3; ++i) {
        if (abs(pos[i]) > 5)
            return 1;
    }
    return -1;
}

void testMC()
{
    Field3D sphereField(VEC3I(20, 20, 20), VEC3F(-10,-10,-10), VEC3F(10,10,10), &sphereFunction);

    Mesh mesh;
    MC::marching_cube(sphereField, mesh);

    mesh.writeOBJ("object.obj");
}

int main()
{
    testMC();
    return 0;
}
