#include "SETTINGS.h"
#include "field.h"
#include "mesh.h"
#include "Quaternion/POLYNOMIAL_4D.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

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

int main(int argc, char* argv[]) {

    if(argc < 6) {
        cout << "USAGE: " << endl;
        cout << "To generate a random polynomial in a bounding box:" << endl;
        cout << " " << argv[0] << " <numRoots> <min power> <max power> <*.poly output> BOX <x min> <y min> <z min> <x max> <y max> <z max>" << endl;;

        cout << "To generate a random polynomial in an OBJ's bounding box:" << endl;
        cout << " " << argv[0] << " <numRoots> <min power> <max power> <*.poly output> OBJBOX <*.obj model>" << endl;;

        cout << "To generate a random polynomial near an SDF's surface:" << endl;
        cout << " " << argv[0] << " <numRoots> <min power> <max power> <*.poly output> SDF <*.f3d sdf> <max distance>" << endl;;
        exit(-1);
    }

    POLYNOMIAL_4D poly;

    int numRoots = atoi(argv[1]);
    int minPower = atoi(argv[2]);
    int maxPower = atoi(argv[3]);

    string outputFilename(argv[4]);

    string mode(argv[5]);

    if (mode == "BOX") {
        AABB box(VEC3F(atof(argv[6]), atof(argv[7]), atof(argv[8])), VEC3F(atof(argv[9]), atof(argv[10]), atof(argv[11])));
        poly = randPolynomialInBox(box, minPower, maxPower, numRoots, false);
    }

    if (mode == "OBJBOX") {
        Mesh m;
        m.readOBJ(string(argv[6]));

        AABB box(m.vertices[0], m.vertices[1]);
        for (const VEC3F& v : m.vertices) {
            box.max() = box.max().cwiseMax(v);
            box.min() = box.min().cwiseMin(v);
        }

        poly = randPolynomialInBox(box, minPower, maxPower, numRoots, false);
    }

    if (mode == "SDF") {
        ArrayGrid3D sdf(argv[6]);
        sdf.mapBox.max() = VEC3F(0.5, 0.5, 0.5);
        sdf.mapBox.min() = VEC3F(-0.5, -0.5, -0.5);

        Real maxDist = atof(argv[7]);

        poly = randPolynomialNearSurface(&sdf, minPower, maxPower, numRoots, maxDist);
    }

    FILE* out = fopen(outputFilename.c_str(), "w");
    poly.write(out);
    fclose(out);

    PRINTF("Wrote polynomial with %d roots to %s", poly.totalRoots(), outputFilename.c_str())

        return 0;
}
