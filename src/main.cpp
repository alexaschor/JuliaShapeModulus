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

ArrayGrid3D* debugField;

class ConstrainedRationalMap: public FieldFunction3D {
//private:
public:
    Grid3D* distanceField;
    POLYNOMIAL_4D topPolynomial;
    bool hasBottomPolynomial;
    POLYNOMIAL_4D bottomPolynomial;
    Real c;
    int maxIterations;
    Real escape;

public:
    ConstrainedRationalMap(Grid3D* distanceField, POLYNOMIAL_4D topPolynomial, Real c = 300, int maxIterations = 3, Real escape = 20):
        distanceField(distanceField), topPolynomial(topPolynomial), hasBottomPolynomial(false), c(c), maxIterations(maxIterations), escape(escape) {}

    Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION iterate(pos[0], pos[1], pos[2], 0);
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        while (magnitude < escape && totalIterations < maxIterations) {
            // First lookup the radius we're going to project to and save the original
            const Real distance = (*distanceField)(VEC3F(iterate[0], iterate[1], iterate[2]));
            Real radius = exp(c * distance);

            QUATERNION original = iterate;

            // Arbitrary: we define NAN as escaping
            if (isnan(radius)) {
                magnitude = escape;
                break;
            }

            // If we know it'll escape we can skip quaternion multiplication evaluation
            if (radius > escape) {
                magnitude = radius;
                totalIterations++;
                break;
            }

            // iterate = topPolynomial.evaluate(iterate);
            iterate = topPolynomial.evaluateScaledPowerFactored(iterate);
            if (hasBottomPolynomial) {
                QUATERNION bottomEval = bottomPolynomial.evaluateScaledPowerFactored(iterate);
                iterate = (iterate / bottomEval);
            }

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

    if(argc != 6) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped Julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <distance field> <4D top polynomial> <4D bottom polynomial, or NONE> <output resolution> <output OBJ>" << endl;
        exit(0);
    }

    ArrayGrid3D distFieldCoarse(argv[1]);

    InterpolationGrid distField(&distFieldCoarse, InterpolationGrid::LINEAR);
    distField.mapBox.setCenter(VEC3F(0,0,0));
    PRINTF("Got distance field with res %dx%dx%d\n", distField.xRes, distField.yRes, distField.zRes);

    PRINT("NOTE: Setting simulation bounds to hard-coded values");
    distField.mapBox.max() = VEC3F(1.25, 1.25, 1.25);
    distField.mapBox.min() = VEC3F(-1.25, -1.25, -1.25);


    if (string(argv[3]) != "NONE") {
        POLYNOMIAL_4D polyTop(argv[2]);
        ConstrainedRationalMap map(&distField, polyTop, 100);

        POLYNOMIAL_4D polyBottom(argv[3]);
        map.hasBottomPolynomial = true;
        map.bottomPolynomial = polyBottom;

        int res = atoi(argv[4]);
        VirtualGrid3DLimitedCache vg(res, res, res, distField.mapBox.min(), distField.mapBox.max(), &map);

        Mesh m;
        MC::march_cubes(&vg, m, true);
        m.writeOBJ(argv[5]);
    } else {
        POLYNOMIAL_4D polyTop(argv[2]);
        ConstrainedRationalMap map(&distField, polyTop, 100);

        int res = atoi(argv[4]);
        VirtualGrid3DLimitedCache vg(res, res, res, distField.mapBox.min(), distField.mapBox.max(), &map);

        Mesh m;
        MC::march_cubes(&vg, m, true);
        m.writeOBJ(argv[5]);
    }

    return 0;
}
