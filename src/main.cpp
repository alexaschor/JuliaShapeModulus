#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdio>
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

POLYNOMIAL_4D randPolynomialInBox(Real radius, Real minPower, Real maxPower, uint nRoots, bool allowNonIntegerPowers = true) {
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<> spaceD(-radius, radius);
    uniform_real_distribution<> powerD(minPower, maxPower);

    vector<QUATERNION> roots{};
    vector<Real> powers{};
    for (uint i = 0; i < nRoots; ++i) {
        roots.push_back(QUATERNION(spaceD(eng), spaceD(eng), spaceD(eng), spaceD(eng)));
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

void dumpRotInfo(QuatToQuatFn* p, AABB range, int res) {

    QuatQuatRotField r(p, 0);
    QuatQuatRotField g(p, 1);
    QuatQuatRotField b(p, 2);
    QuatQuatMagField mag(p);

    VirtualGrid3D gr(res, res, res, range.min(), range.max(), &r);
    gr.writeF3D("temp/r.f3d", true);

    VirtualGrid3D gg(res, res, res, range.min(), range.max(), &g);
    gg.writeF3D("temp/g.f3d", true);

    VirtualGrid3D gb(res, res, res, range.min(), range.max(), &b);
    gb.writeF3D("temp/b.f3d", true);

    VirtualGrid3D gm(res, res, res, range.min(), range.max(), &mag);
    gm.writeF3D("temp/mag.f3d", true);
}

// ========================= SPHERE INTERSECTION SAMPLING ============================

// Return percent of sample points inside SDF on surface of sphere
Real sampleSpherePoints(Grid3D* sdf, Real radius, uint numSamples) {
    uint numInside = 0;

    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<> dist(-1, 1);

    for (uint i = 0; i < numSamples; ++i) {
        VEC3F sp(dist(eng), dist(eng), dist(eng));
        sp.normalize();
        sp *= radius;

        if (sdf->getFieldValue(sp) < 0) numInside++;

    }
    return ((Real) numInside) / numSamples;
}

void sampleSphereIntersection(Grid3D* sdf, Real minRadius, Real maxRadius, uint numSpheres, uint numSamples, string outputFile) {
    ofstream out;
    out.open(outputFile);

    PB_START("Sampling sphere intersection curve");

    for (uint i = 0; i < numSpheres; ++i) {
        Real sampleRadius = minRadius + ((maxRadius - minRadius) * ((Real) i / numSpheres));
        out << sampleRadius << ", " << sampleSpherePoints(sdf, sampleRadius, numSamples) << endl;
        PB_PROGRESS((float) i / numSpheres);
    }

    PB_END();

    out.close();
}

// ===================================================================================

int main(int argc, char *argv[]) {

    if(argc != 6) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <4D top polynomial, or RANDOM> <output resolution> <fill level> <output OBJ>" << endl;
        exit(0);
    }

    // Read distfield
    ArrayGrid3D distFieldCoarse(argv[1]);
    PRINTF("Got distance field with res %dx%dx%d\n", distFieldCoarse.xRes, distFieldCoarse.yRes, distFieldCoarse.zRes);

    // Create interpolation grid (smooth it out)
    InterpolationGrid distField(&distFieldCoarse, InterpolationGrid::LINEAR);
    distField.mapBox.setCenter(VEC3F(0,0,0));


    // Save radius-membership relationship
    // sampleSphereIntersection(&distField, 0, 0.6, 2500, 100000, "inside.csv");
    // exit(0);

    PRINT("NOTE: Setting simulation bounds to hard-coded values (not from distance field)");
    distField.mapBox.max() = VEC3F(0.5, 0.5, 0.5);
    distField.mapBox.min() = VEC3F(-0.5, -0.5, -0.5);

    QuatToQuatFn *p;

    if (string(argv[2]) == "RANDOM") {
        // Create random polynomial
        // POLYNOMIAL_4D polyTop = randPolynomialInBox(0.7, 8, 13, 25, false);
        POLYNOMIAL_4D polyTop = randPolynomialInObject(&distField, 2, 13, 10, false);
        p = new RationalQuatPoly(polyTop);
    } else {
        POLYNOMIAL_4D polyTop(argv[2]);
        p = new RationalQuatPoly(polyTop);
    }

    // Save polynomial info for inspection later
    // FILE* out = fopen("temp/p.poly4d", "w");
    //     polyTop.write(out);
    // fclose(out);
    // dumpRotInfo(&p, distField.mapBox, 100);

    // Finally we actually compute the Julia set
    Real fillLevel = atof(argv[4]); // This is C in the equation
    DistanceGuidedQuatMap map(&distField, p, fillLevel, 3);

    int res = atoi(argv[3]);
    VirtualGrid3DLimitedCache vg(res, res, res, distField.mapBox.min(), distField.mapBox.max(), &map);
    Mesh m;
    MC::march_cubes(&vg, m, true);
    m.writeOBJ(argv[5]);

    return 0;
}

