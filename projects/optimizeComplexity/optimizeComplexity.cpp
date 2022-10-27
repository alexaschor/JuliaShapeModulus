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
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/neldermeadsolver.h"

#include "julia.h"

#include <chrono>
#include <thread>

using namespace std;

enum COMPLEXITY_METRIC {
    surfArea, triCount, meshCurvature
};

// Global vars for params, not ideal but expedient here


Real g_a;
Real g_b;
RationalQuatPoly* g_poly;
Grid3D* g_distField;
uint g_iterations = 3;
Real g_escape = 20;

uint g_res;
COMPLEXITY_METRIC metric;

///////////////////////////////////////////////////////////////////////
// Offset roots and march cubes
///////////////////////////////////////////////////////////////////////
Real computeFractal(VEC3F translation, bool save = false, string filename = "")
{
    // Set translation
    VEC3F oldCenter = g_distField->mapBox.center();
    g_distField->mapBox.setCenter(translation);
    for (unsigned int i = 0; i < g_poly->topPolynomial.roots().size(); i++)
        g_poly->topPolynomial.rootsMutable()[i] = g_poly->topPolynomial.roots()[i] + translation;

    PRINTF("Computing Julia set with resolution %d, a=%f, b=%f, offset=(%f, %f, %f)\n", g_res, g_a, g_b, translation.x(), translation.y(), translation.z());

    DistanceGuidedQuatFn map(g_distField, g_poly, g_a, g_b);
    QuaternionJuliaSet julia(&map, g_iterations, g_escape);

    AABB boundsBox(g_distField->mapBox.min(), g_distField->mapBox.max());

    PRINTV3(boundsBox.min());
    PRINTV3(boundsBox.max());

    VirtualGrid3DLimitedCache vg(g_res, g_res, g_res, boundsBox.min(), boundsBox.max(), &julia);
    Mesh m;
    MC::march_cubes(&vg, m, true);

    if (save) m.writeOBJ(filename);

    // Undo translation
    g_distField->mapBox.setCenter(oldCenter);
    for (unsigned int i = 0; i < g_poly->topPolynomial.roots().size(); i++)
        g_poly->topPolynomial.rootsMutable()[i] = g_poly->topPolynomial.roots()[i] - translation;

    return m.computeSurfaceArea();
}


namespace cppoptlib {

    class FractalComplexity : public Problem<double, 3> {
    public:
        double value(const TVector &x) {
            VEC3F translation(x[0], x[1], x[2]);
            // Grab fractal and res from global vars
            Real value = -computeFractal(translation, true, "monitor.obj");
            return value;
        }

        bool detailed_callback(const cppoptlib::Criteria<double> &state, SimplexOp op, int index, const MatrixType &x, std::vector<Scalar> f) {
            TVector xp = x.col(index).transpose();

            printf("* %zu: [%.17g, %.17g, %.17g] => %.17g\n", state.iterations, xp[0], xp[1], xp[2], f[index]);

            return true;
        }
    };

}


int main(int argc, char *argv[]) {

    if(argc != 11) {
        cout << "USAGE: " << endl;
        cout << "To create a shaped julia set from a distance field:" << endl;
        cout << " " << argv[0] << " <SDF *.f3d> <P(q) *.poly4d> <output resolution> <a> <b> <offset x> <offset y> <offset z> <monitor *.obj> <num iterations>" << endl << endl;

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

        exit(0);
    }

    // Read distfield
    ArrayGrid3D distFieldCoarse(argv[1]);
    PRINTF("Got distance field with res %dx%dx%d\n", distFieldCoarse.xRes, distFieldCoarse.yRes, distFieldCoarse.zRes);
    InterpolationGrid distField(&distFieldCoarse, InterpolationGrid::LINEAR);
    distField.mapBox.setCenter(VEC3F(0,0,0));

    PRINT("NOTE: Setting simulation bounds to hard-coded values (not from distance field)");
    distField.mapBox.min() = VEC3F(-0.5, -0.5, -0.5);
    distField.mapBox.max() = VEC3F(0.5, 0.5, 0.5);

    g_distField = &distField;

    PRINT("Loading polynomial from file...");
    POLYNOMIAL_4D polyTop(argv[2]);
    RationalQuatPoly rqp(polyTop);
    g_poly = &rqp;

    g_res = atoi(argv[3]);
    g_a = atof(argv[4]); // This is param a in the equation
    g_b = atof(argv[5]); // This is param b in the equation

    // Now optimize

    // Get args
    cppoptlib::FractalComplexity::TVector x{atof(argv[6]), atof(argv[7]), atof(argv[8])};
    int optIterations = atoi(argv[10]);

    // Construct problem and solver
    cppoptlib::FractalComplexity prob;

    cppoptlib::NelderMeadSolver<cppoptlib::FractalComplexity> solver;
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.iterations = optIterations;
    // crit.fDelta = 0.0005;
    solver.setStopCriteria(crit);

    // Run solver
    solver.minimize(prob, x); // stores argmin to x

    std::cout << "Offset with maximum complexity: " << x.transpose() << std::endl;



    return 0;
}

