#ifndef SETTINGS_H
#define SETTINGS_H

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include <fstream>

using namespace Eigen;

typedef double Real;
typedef unsigned int uint;
typedef Matrix<Real, 2, 1 > VEC2F;
typedef Matrix<Real, 3, 1 > VEC3F;
typedef Matrix<int, 3, 1 > VEC3I;
typedef VectorXd VECTOR;

#define MC_MAX_ROOTFINDING_ITERATIONS 100
#define MC_ROOTFINDING_THRESH 1e-8

// DEBUGGING MACROS
#define DEBUGBOOL true

#define PRINTV3(v) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s={%.2e, %.2e, %.2e};\n", __FILE__, __LINE__, __func__, #v, (v)[0], (v)[1], (v)[2])
#define PRINTV4(v) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s={%.2e, %.2e, %.2e, %.2e};\n", __FILE__, __LINE__, __func__, #v, (v)[0], (v)[1], (v)[2], (v)[3])
#define PRINTD(x) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s=%f;\n", __FILE__, __LINE__, __func__, #x, x)
#define PRINTI(x) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s=%d;\n", __FILE__, __LINE__, __func__, #x, x)
#define PRINTF(format, ...) if (DEBUGBOOL) {fprintf(stderr, "%s:%d (%s)\t| ", __FILE__, __LINE__, __func__); fprintf(stderr, format, ## __VA_ARGS__);}
#define PRINT(x) if (DEBUGBOOL) {fprintf(stderr, "%s:%d (%s)\t| ", __FILE__, __LINE__, __func__); fprintf(stderr, x); fprintf(stderr, "\n");}

#define HERE() PRINT("here!")

#define SOMETIMES() float SOMETIMES_RAND = rand() % 1000; bool SOMETIMES_TRACKER_VAR = SOMETIMES_RAND == 1; if (SOMETIMES_TRACKER_VAR)
#define SCONT() if (SOMETIMES_TRACKER_VAR)

// TIMING MACROS
#define TIMER_INIT() std::chrono::time_point<std::chrono::system_clock> TIMER_START_TIME, TIMER_END_TIME; std::chrono::duration<double> TIMER_DIFF; double TIMER_DURATION;

#define TIMER_START() TIMER_START_TIME = std::chrono::system_clock::now()
#define TIMER_END()   TIMER_END_TIME   = std::chrono::system_clock::now(); TIMER_DIFF = TIMER_END_TIME - TIMER_START_TIME; TIMER_DURATION = TIMER_DIFF.count()

#define TIME(exp) TIMER_START(); exp TIMER_END(); PRINTD(TIMER_DURATION);

// Read and write VEC3F from file
namespace MyEigen {
    inline void write_vec3f(FILE* file, const VEC3F& vec) {
        if (sizeof(Real) == sizeof(double))
        {
            fwrite((void*)&vec[0], sizeof(double), 1, file);
            fwrite((void*)&vec[1], sizeof(double), 1, file);
            fwrite((void*)&vec[2], sizeof(double), 1, file);
        }
        else
        {
            double elementD[3];
            elementD[0] = vec[0];
            elementD[1] = vec[1];
            elementD[2] = vec[2];
            fwrite((void*)&elementD[0], sizeof(double), 1, file);
            fwrite((void*)&elementD[1], sizeof(double), 1, file);
            fwrite((void*)&elementD[2], sizeof(double), 1, file);
        }
    }


    inline void read_vec3f(FILE* file, VEC3F& vec) {
        // make sure it reads into a double
        double elementD[3];
        fread((void*)&elementD[0], sizeof(double), 1, file);
        fread((void*)&elementD[1], sizeof(double), 1, file);
        fread((void*)&elementD[2], sizeof(double), 1, file);

        // let the assign resolve the type issues
        vec[0] = elementD[0];
        vec[1] = elementD[1];
        vec[2] = elementD[2];
    }

}



#endif
