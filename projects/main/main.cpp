#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdio>
#include <limits>
#include <stdio.h>
#include <cmath>
#include <random>

#include "SETTINGS.h"

#include "MC.h"
#include "mesh.h"
#include "field.h"
#include "Quaternion/QUATERNION.h"
#include "Quaternion/POLYNOMIAL_4D.h"

#include "julia.h"

#include <chrono>
#include <thread>

using namespace std;

POLYNOMIAL_4D randPolynomialInBox(AABB box, Real minPower, Real maxPower, uint nRoots, bool allowNonIntegerPowers = true) {
    random_device rd;
    default_random_engine eng(rd());

    uniform_real_distribution<> spaceXD(box.min()[0], box.max()[0]);
    uniform_real_distribution<> spaceYD(box.min()[1], box.max()[1]);
    uniform_real_distribution<> spaceZD(box.min()[2], box.max()[2]);

    uniform_real_distribution<> powerD(minPower, maxPower);

    vector<QUATERNION> roots{};
    vector<Real> powers{};
    for (uint i = 0; i < nRoots; ++i) {
        QUATERNION rootPos = QUATERNION(spaceXD(eng), spaceYD(eng), spaceZD(eng), 0);
        roots.push_back(rootPos);
        powers.push_back( allowNonIntegerPowers? powerD(eng) : int(powerD(eng)) );
    }

    POLYNOMIAL_4D poly(roots, powers);
    return poly;
}

POLYNOMIAL_4D randPolynomialInObject(Grid3D* sdf, Real minPower, Real maxPower, uint nRoots, bool allowNonIntegerPowers = true) {
    random_device rd;
    default_random_engine eng(rd());

    Real minX = 0, maxX = 1, minY = 0, maxY = 1, minZ = 0, maxZ = 1;
    if (sdf->hasMapBox) {
        minX = sdf->mapBox.min()[0];
        minY = sdf->mapBox.min()[1];
        minZ = sdf->mapBox.min()[2];

        maxX = sdf->mapBox.max()[0];
        maxY = sdf->mapBox.max()[1];
        maxZ = sdf->mapBox.max()[2];
    }

    uniform_real_distribution<> spaceXD(minX, maxX);
    uniform_real_distribution<> spaceYD(minY, maxY);
    uniform_real_distribution<> spaceZD(minZ, maxZ);

    uniform_real_distribution<> powerD(minPower, maxPower);

    vector<QUATERNION> roots{};
    vector<Real> powers{};
    for (uint i = 0; i < nRoots; ++i) {
        QUATERNION rootPos;
        bool insideTarget = false;
        while (not insideTarget) {
            rootPos = QUATERNION(spaceXD(eng), spaceYD(eng), spaceZD(eng), 0);
            insideTarget = (*sdf).getFieldValue(VEC3F(rootPos[0], rootPos[1], rootPos[2])) < 0;
        }
        roots.push_back(rootPos);
        powers.push_back( allowNonIntegerPowers? powerD(eng) : int(powerD(eng)) );
    }

    POLYNOMIAL_4D poly(roots, powers);
    return poly;
}

POLYNOMIAL_4D randPolynomialNearSurface(Grid3D* sdf, Real minPower, Real maxPower, uint nRoots, Real maxDist, bool allowNonIntegerPowers = true) {
    random_device rd;
    default_random_engine eng(rd());

    Real minX = 0, maxX = 1, minY = 0, maxY = 1, minZ = 0, maxZ = 1;
    if (sdf->hasMapBox) {
        minX = sdf->mapBox.min()[0];
        minY = sdf->mapBox.min()[1];
        minZ = sdf->mapBox.min()[2];

        maxX = sdf->mapBox.max()[0];
        maxY = sdf->mapBox.max()[1];
        maxZ = sdf->mapBox.max()[2];
    }

    uniform_real_distribution<> spaceXD(minX, maxX);
    uniform_real_distribution<> spaceYD(minY, maxY);
    uniform_real_distribution<> spaceZD(minZ, maxZ);

    uniform_real_distribution<> powerD(minPower, maxPower);

    vector<QUATERNION> roots{};
    vector<Real> powers{};
    for (uint i = 0; i < nRoots; ++i) {
        QUATERNION rootPos;
        bool insideTarget = false;
        while (not insideTarget) {
            rootPos = QUATERNION(spaceXD(eng), spaceYD(eng), spaceZD(eng), 0);
            insideTarget = abs((*sdf).getFieldValue(VEC3F(rootPos[0], rootPos[1], rootPos[2]))) < maxDist;
        }
        roots.push_back(rootPos);
        powers.push_back( allowNonIntegerPowers? powerD(eng) : int(powerD(eng)) );
    }

    POLYNOMIAL_4D poly(roots, powers);
    return poly;
}

POLYNOMIAL_4D randPolynomialNearSurfaceBlueNoise(Grid3D* sdf, Real minPower, Real maxPower, uint nRoots, Real maxDist, bool allowNonIntegerPowers = true) {
    random_device rd;
    default_random_engine eng(rd());

    Real minX = 0, maxX = 1, minY = 0, maxY = 1, minZ = 0, maxZ = 1;
    if (sdf->hasMapBox) {
        minX = sdf->mapBox.min()[0];
        minY = sdf->mapBox.min()[1];
        minZ = sdf->mapBox.min()[2];

        maxX = sdf->mapBox.max()[0];
        maxY = sdf->mapBox.max()[1];
        maxZ = sdf->mapBox.max()[2];
    }

    uniform_real_distribution<> spaceXD(minX, maxX);
    uniform_real_distribution<> spaceYD(minY, maxY);
    uniform_real_distribution<> spaceZD(minZ, maxZ);

    uniform_real_distribution<> powerD(minPower, maxPower);

    vector<QUATERNION> roots{};
    vector<Real> powers{};
    for (uint i = 0; i < nRoots; ++i) {
        QUATERNION rootPos;
        bool insideTarget = false;
        while (not insideTarget) {
            rootPos = QUATERNION(spaceXD(eng), spaceYD(eng), spaceZD(eng), 0);
            insideTarget = abs((*sdf).getFieldValue(VEC3F(rootPos[0], rootPos[1], rootPos[2]))) < maxDist;
        }
        roots.push_back(rootPos);
        powers.push_back( allowNonIntegerPowers? powerD(eng) : int(powerD(eng)) );
    }

    POLYNOMIAL_4D poly(roots, powers);
    return poly;
}

void dumpVersorModulus(QuatToQuatFn* p, AABB range, int res) {
    QuatQuatRotField x(p, 0);
    QuatQuatRotField y(p, 1);
    QuatQuatRotField z(p, 2);
    QuatQuatMagField mag(p);

    VirtualGrid3D gx(res, res, res, range.min(), range.max(), &x);
    gx.writeF3D("temp/x.f3d", true);

    VirtualGrid3D gy(res, res, res, range.min(), range.max(), &y);
    gy.writeF3D("temp/y.f3d", true);

    VirtualGrid3D gz(res, res, res, range.min(), range.max(), &z);
    gz.writeF3D("temp/z.f3d", true);

    VirtualGrid3D gm(res, res, res, range.min(), range.max(), &mag);
    gm.writeCSV("temp/modulus.csv");
}

Real GLOBAL_fillLevelGradOffset = 0; // XXX This is kinda hacky because otherwise
                                     // we can't capture the local variable in a
                                     // lambda.

int main(int argc, char *argv[]) {
    if(argc != 10 && argc != 11) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <SDF *.f3d> <P(q) *.poly4d, or the word RANDOM> <output resolution> <a> <b> <offset x> <offset y> <offset z> <output *.obj> <optional: octree specifier string>" << endl << endl;

        cout << "    This will compute the 3D Julia quaternion Julia set of the function:" << endl;
        cout << "        f(q) = r(q) * d(q)" << endl;
        cout << "    Where:" << endl;
        cout << "        r(q) = e^( a * SDF(q) + b )" << endl;
        cout << "        d(q) = Normalized( P(q) )" << endl << endl;

        cout << "    P(q) is the quaternion polynomial function, which determines the character of the fractal detail." << endl;
        cout << "    a is a parameter which controls the thickness of the shell in which the chaotic effect has significant influence" << endl;
        cout << "    on set membership, and b is a parameter which controls the position along the surface of the SDF of that shell." << endl << endl;

        cout << "    The offset X, Y, and Z parameters move the origin of the dynamical system around in space, which causes" << endl;
        cout << "    the Julia set to dissolve in interesting ways." << endl;

        cout << "    The octree specifier string is an optional parameter useful for computing large Julia sets in parallel." << endl;
        cout << "    You can select a box in an evenly-subdivided octree of arbitrary depth specified by a string of digits 0-7," << endl;
        cout << "    laid out as follows:" << endl;

        cout << "                 +---+" << endl;
        cout << "               / |4|5|" << endl;
        cout << "              /  +-+-+" << endl;
        cout << "             /   |7|6|" << endl;
        cout << "            /    +---+" << endl;
        cout << "           +---+    / " << endl;
        cout << "           |0|1|   /  " << endl;
        cout << "           +-+-+  /   " << endl;
        cout << "           |3|2| /    " << endl;
        cout << "        ▲  +---+      " << endl;
        cout << "        |             " << endl;
        cout << "        Y X--▶        " << endl;
        cout << "        Z ●           " << endl;
        cout << "        (into page)   " << endl;

        cout << "    Each character of the string will go one level deeper, so the string '5555' specifies the 1/16-edge length box at the far back corner." << endl;

        exit(0);
    }

    // Read distfield
    ArrayGrid3D distFieldCoarse(argv[1]);
    PRINTF("Got distance field with res %dx%dx%d\n", distFieldCoarse.xRes, distFieldCoarse.yRes, distFieldCoarse.zRes);

    // Create interpolation grid (smooth it out)
    InterpolationGrid distField(&distFieldCoarse, InterpolationGrid::LINEAR);
    distField.mapBox.setCenter(VEC3F(0,0,0));

    PRINT("NOTE: Setting simulation bounds to hard-coded values (not from distance field)");
    distField.mapBox.min() = VEC3F(-0.5, -0.5, -0.5);
    distField.mapBox.max() = VEC3F(0.5, 0.5, 0.5);

    RationalQuatPoly *p;

    if (string(argv[2]) == "RANDOM") {
        // Create random polynomial
        PRINT("Generating random polynomial...");
        POLYNOMIAL_4D polynomial = randPolynomialInBox(distField.mapBox, 2, 13, 10, false);

        p = new RationalQuatPoly(polynomial);
    } else {
        // Load from file
        PRINT("Loading polynomial from file...");
        POLYNOMIAL_4D polynomial(argv[2]);
        p = new RationalQuatPoly(polynomial);
    }

    // Save polynomial info for inspection later (helpful if it was randomly generated)
    FILE* out = fopen("temp/p.poly4d", "w");
    p->topPolynomial.write(out);
    PRINT("Wrote polynomial to temp/p.poly4d");
    fclose(out);

    // Finally we actually compute the Julia set
    Real fillLevel = atof(argv[4]); // This is param alpha in the equation
    Real fillOffset = atof(argv[5]); // This is param beta in the equation

    // Offset roots and distance field to reproduce QUIJIBO dissolution
    // effect - this is optional, but produces effects like in Fig. 5
    VEC3F offset3D(atof(argv[6]), atof(argv[7]), atof(argv[8]));

    distField.mapBox.setCenter(offset3D);
    for (unsigned int i = 0; i < p->topPolynomial.roots().size(); i++) {
        p->topPolynomial.rootsMutable()[i] = p->topPolynomial.roots()[i] + offset3D;
    }

    int res = atoi(argv[3]);

    // Set up simulation bounds, taking octree zoom into account
    AABB boundsBox(distField.mapBox.min(), distField.mapBox.max());
    if (argc == 11) { // If an octree specifier string was given, we zoom in on just one box
        char* octreeStr = argv[10];

        for (size_t i = 0; i < strlen(octreeStr); ++i) {
            int oIdx = octreeStr[i] - '0';
            if (oIdx < 0 || oIdx > 7) {
                PRINTF("Found character '%c' in octree specifier string. Valid characters are numbers 0-7, inclusive.\n", octreeStr[i]);
                exit(1);
            }

            boundsBox = boundsBox.subdivideOctree()[oIdx];
        }

        // Pad by one grid cell to avoid gaps
        VEC3F delta = boundsBox.span() / res;
        boundsBox.max() += delta;
        boundsBox.min() -= delta;
    }

    PRINTF("Computing Julia set with resolution %d, a=%f, b=%f, offset=(%f, %f, %f)\n", res, fillLevel, fillOffset, offset3D.x(), offset3D.y(), offset3D.z());

    // Make fields for a and b
    ConstantFunction3D b(fillOffset);
    VirtualGrid3D(100, 100, 100, VEC3F(0,0,0), VEC3F(1,1,1), &b).writeF3D("field.f3d");

    DistanceGuidedQuatFn map(&distField, p, fillLevel, fillOffset);
    QuaternionJuliaSet julia(&map, 3, 20);
    VirtualGrid3DLimitedCache vg(res, res, res, boundsBox.min(), boundsBox.max(), &julia);

    Mesh m;
    MC::march_cubes(&vg, m, true);

    // Currently march_cubes doesn't take the grid's mapBox into account; all vertices are
    // placed in [ (0, xRes), (0, yRes), (0, zRes) ] space. TODO fix march_cubes to account for
    // the mapBox, but for now we'll just manually transform it. Normals should be okay as they are.
    for (uint i = 0; i < m.vertices.size(); ++i) {
        VEC3F v = m.vertices[i];
        m.vertices[i] = vg.gridToFieldCoords(v);
    }

    m.writeOBJ(argv[9]);

    return 0;
}

