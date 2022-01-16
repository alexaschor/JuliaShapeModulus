#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdio>
#include <cmath>

#include "SETTINGS.h"

#include "MC.h"
#include "mesh.h"
#include "field.h"
#include "Quaternion/QUATERNION.h"
#include "Quaternion/POLYNOMIAL_4D.h"

#include <chrono>
#include <thread>

using namespace std;

const Real ESCAPE = 20.0;
const int MAX_ITERATIONS = 3;

ArrayGrid3D* debugField;

class ConstrainedRationalMap: public FieldFunction3D {
private:
    Grid3D* distanceField;
    POLYNOMIAL_4D polynomial;
    Real a, b;
    int maxIterations;
    Real escape;

public:
    ConstrainedRationalMap(Grid3D* distanceField, POLYNOMIAL_4D polynomial, Real a = 2, Real b = 300, int maxIterations = 3, Real escape = 20):
        distanceField(distanceField), polynomial(polynomial), a(a), b(b), maxIterations(maxIterations), escape(escape) {}

    Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION iterate(pos[0], pos[1], pos[2], 0);
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        while (magnitude < ESCAPE && totalIterations < MAX_ITERATIONS) {
            // First lookup the radius we're going to project to and save the original
            const Real distance = (*distanceField)(VEC3F(iterate[0], iterate[1], iterate[2]));
            Real radius = pow(a, b * distance);

            QUATERNION original = iterate;

            // Arbitrary: we define NAN as escaping
            if (isnan(radius)) {
                magnitude = ESCAPE;
                break;
            }

            // If we know it'll escape we can skip quaternion multiplication evaluation
            if (radius > ESCAPE) {
                magnitude = radius;
                totalIterations++;
                break;
            }

            iterate = polynomial.evaluate(iterate);

            // If quaternion multiplication fails, revert back to original
            if (iterate.anyNans()) {
                iterate = original;
            }

            // Try to normalize iterate, might produce infs or nans if magnitude is too small
            QUATERNION normedIterate = iterate;
            normedIterate.normalize();

            bool tooSmall = normedIterate.anyNans();
            if (tooSmall) {
                original.normalize();
                if (original.anyNans()) {
                    // Need to put it at radius but have no direction info
                    iterate = QUATERNION(1,0,0,0) * radius;
                } else {
                    iterate = original * radius;
                }
            } else {
                iterate = normedIterate;
                iterate *= radius;
            }

            magnitude = iterate.magnitude();

            totalIterations++;
        }

        return log(magnitude);
    }

};

Real sphereFunction(VEC3F pos) {
    return (pos.norm() - 1);
}

int main(int argc, char *argv[]) {

    if(argc != 4) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <4D polynomial> <output OBJ>" << endl;
        exit(0);
    }

    ArrayGrid3D distFieldCoarse(argv[1]);
    InterpolationGrid distField(&distFieldCoarse, InterpolationGrid::LINEAR);

    // distField.writeF3D("dfi.field3d", true);

    distField.mapBox.setCenter(VEC3F(0,0,0));

    POLYNOMIAL_4D poly(argv[2]);

    PRINTF("Got distance field with res %dx%dx%d\n", distField.xRes, distField.yRes, distField.zRes);

    ConstrainedRationalMap map(&distField, poly, 2, 300);
    // FieldFunction3D map(sphereFunction);

    VirtualGrid3DLimitedCache vg(100, 100, 100, distField.mapBox.min(), distField.mapBox.max(), &map);
    // ArrayGrid3D vg(300, 300, 300, distField.mapBox.min(), distField.mapBox.max(), &map);

    // FieldFunction3D sf(&sphereFunction);
    // VirtualGrid3DLimitedCache vg(300, 300, 300, VEC3F(-10, -10, -10), VEC3F(10,10,10), &sf);
    // vg.writeF3D("sphere.field3d");

    // vg.writeF3D("fieldOut.field3d");
    Mesh m;
    MC::march_cubes(&vg, m, true);
    m.writeOBJ(argv[3]);

    return 0;
}
