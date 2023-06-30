#include "Quaternion/QUATERNION.h"
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

    if(argc == 1) {
        cout << "USAGE: " << endl;
        cout << "-----------------------------------------------" << endl;

        cout << "To inspect a *.p4d polynomial:" << endl;
        cout << " " << argv[0] << " INS <*.p4d file> " << endl;;

        cout << "-----------------------------------------------" << endl;

        cout << "To generate a random polynomial in a bounding box:" << endl;
        cout << " " << argv[0] << " GEN <numRoots> <min power> <max power> <*.p4d output> BOX <x min> <y min> <z min> <x max> <y max> <z max>" << endl;;

        cout << "To generate a random polynomial in an OBJ's bounding box:" << endl;
        cout << " " << argv[0] << " GEN <numRoots> <min power> <max power> <*.p4d output> OBJBOX <*.obj model>" << endl;;

        cout << "To generate a random polynomial near an SDF's surface:" << endl;
        cout << " " << argv[0] << " GEN <numRoots> <min power> <max power> <*.p4d output> SDF <*.f3d sdf> <max distance>" << endl;;

        cout << "-----------------------------------------------" << endl;

        cout << "To generate a polynomial whose roots are the union of two input polynomials:" << endl;
        cout << " " << argv[0] << " ADD <first *.p4d input> <second *.p4d input> <*.p4d output> " << endl;;

        cout << "-----------------------------------------------" << endl;

        cout << "To scale the positions of a polynomial's roots in space:" << endl;
        cout << " " << argv[0] << " SCALE <*.p4d input> <factor> <*.p4d output>" << endl;;

        cout << "-----------------------------------------------" << endl;

        exit(-1);
    }
    if (string(argv[1]) == "ADD") {
        PRINTF("Reading %s...\n", argv[2]);
        POLYNOMIAL_4D p1(argv[2]);
        PRINTF("Got polynomial with %d roots.\n", p1.totalRoots());

        PRINTF("Reading %s...\n", argv[3]);
        POLYNOMIAL_4D p2(argv[3]);
        PRINTF("Got polynomial with %d roots.\n", p2.totalRoots());

        for (auto q : p2.rootsMutable()) {
            p1.addRoot(q);
        }

        PRINTF("Writing %d roots to %s...\n", p1.totalRoots(), argv[4]);

        FILE* out = fopen(argv[4], "w");
        p1.write(out);
        fclose(out);
    } else if (string(argv[1]) == "SCALE") {
        PRINTF("Reading %s...\n", argv[2]);
        POLYNOMIAL_4D p1(argv[2]);
        PRINTF("Got polynomial with %d roots.\n", p1.totalRoots());


        Real factor = atof(argv[3]);
        PRINTF("Scaling all roots by a factor f=%f\n", factor);

        for (auto& q : p1.rootsMutable()) {
            q *= factor;
        }

        PRINTF("Writing %d roots to %s...\n", p1.totalRoots(), argv[4]);

        FILE* out = fopen(argv[4], "w");
        p1.write(out);
        fclose(out);
    } else if (string(argv[1]) == "GEN") {
        POLYNOMIAL_4D poly;
        int numRoots = atoi(argv[2]);
        int minPower = atoi(argv[3]);
        int maxPower = atoi(argv[4]);

        string outputFilename(argv[5]);

        string mode(argv[6]);

        if (mode == "BOX") {
            AABB box(VEC3F(atof(argv[7]), atof(argv[8]), atof(argv[9])), VEC3F(atof(argv[10]), atof(argv[11]), atof(argv[12])));
            poly = randPolynomialInBox(box, minPower, maxPower, numRoots, false);
        }

        if (mode == "OBJBOX") {
            Mesh m;
            m.readOBJ(string(argv[7]));

            AABB box(m.vertices[0], m.vertices[1]);
            for (const VEC3F& v : m.vertices) {
                box.max() = box.max().cwiseMax(v);
                box.min() = box.min().cwiseMin(v);
            }

            poly = randPolynomialInBox(box, minPower, maxPower, numRoots, false);
        }

        if (mode == "SDF") {
            ArrayGrid3D sdf(argv[7]);
            sdf.mapBox.max() = VEC3F(0.5, 0.5, 0.5);
            sdf.mapBox.min() = VEC3F(-0.5, -0.5, -0.5);

            Real maxDist = atof(argv[8]);

            poly = randPolynomialNearSurface(&sdf, minPower, maxPower, numRoots, maxDist);
        }

        FILE* out = fopen(outputFilename.c_str(), "w");
        poly.write(out);
        fclose(out);

        PRINTF("Wrote polynomial with %d roots to %s", poly.totalRoots(), outputFilename.c_str());

        return 0;
    } else if (string(argv[1]) == "INS") {
        PRINT("Reading polynomial...");
        POLYNOMIAL_4D poly(argv[2]);

        PRINTF("Got polynomial with %d roots:\n", poly.totalRoots());

        AABB bounds = AABB::insideOut();
        Real maxW = numeric_limits<Real>::min(), minW = numeric_limits<Real>::max();
        Real minRad = minW, maxRad = maxW;
        Real minPower = minW, maxPower = maxW;

        int i = 0;
        for (QUATERNION root : poly.roots()) {
            VEC3F rootV3(root.x(), root.y(), root.z());

            bounds.include(rootV3);
            if (root.w() > maxW) maxW = root.w();
            if (root.w() < minW) minW = root.w();

            if (root.magnitude() > maxRad) maxRad = root.magnitude();
            if (root.magnitude() < minRad) minRad = root.magnitude();

            if (poly.powers()[i] > maxPower) maxPower = poly.powers()[i];
            if (poly.powers()[i] < minPower) minPower = poly.powers()[i];

            i++;
        }

        PRINTF("Root power min: %f\n", minPower);
        PRINTF("Root power max: %f\n", maxPower);

        PRINTF("Radius min: %f\n", minRad);
        PRINTF("Radius max: %f\n", maxRad);

        PRINTF("Position min: (%.5f, %.5f, %.5f, %.5f)\n", minW, bounds.min().x(), bounds.min().y(), bounds.min().z()) ;
        PRINTF("Position max: (%.5f, %.5f, %.5f, %.5f)\n", maxW, bounds.max().x(), bounds.max().y(), bounds.max().z()) ;

        PRINTF("Bounding box edges: (%.5f, %.5f, %.5f, %.5f)\n", maxW - minW, bounds.span().x(), bounds.span().y(), bounds.span().z());

    }

}
