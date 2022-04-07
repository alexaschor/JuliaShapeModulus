#ifndef SETTINGS_H
#define SETTINGS_H

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include <fstream>
#include <chrono>

using namespace Eigen;

typedef double Real;
typedef unsigned int uint;
typedef Matrix<Real, 2, 1 > VEC2F;
typedef Matrix<Real, 3, 1 > VEC3F;
typedef Matrix<int, 3, 1 > VEC3I;
typedef Matrix<int, 1, 1 > VEC2I;
typedef VectorXd VECTOR;

#define MC_MAX_ROOTFINDING_ITERATIONS 100
#define MC_ROOTFINDING_THRESH 1e-8

// DEBUGGING MACROS
#define DEBUGBOOL true

// Printing various data types along with the variable name, filename and line numbet
#define PRINTV2(v) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s={%.2e, %.2e};\n", __FILE__, __LINE__, __func__, #v, (v)[0], (v)[1])
#define PRINTV3(v) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s={%.2e, %.2e, %.2e};\n", __FILE__, __LINE__, __func__, #v, (v)[0], (v)[1], (v)[2])
#define PRINTV4(v) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s={%.2e, %.2e, %.2e, %.2e};\n", __FILE__, __LINE__, __func__, #v, (v)[0], (v)[1], (v)[2], (v)[3])
#define PRINTD(x) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s=%f;\n", __FILE__, __LINE__, __func__, #x, x)
#define PRINTE(x) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s=%e;\n", __FILE__, __LINE__, __func__, #x, x)
#define PRINTI(x) if (DEBUGBOOL) fprintf(stderr, "%s:%d (%s)\t| %s=%d;\n", __FILE__, __LINE__, __func__, #x, x)
#define PRINTF(format, ...) if (DEBUGBOOL) {fprintf(stderr, "%s:%d (%s)\t| ", __FILE__, __LINE__, __func__); fprintf(stderr, format, ## __VA_ARGS__);}
#define PRINT(x) if (DEBUGBOOL) {fprintf(stderr, "%s:%d (%s)\t| ", __FILE__, __LINE__, __func__); fprintf(stderr, x); fprintf(stderr, "\n");}

// Here!
#define HERE() PRINT("here!")

// Useful for inspecting a loop running thousands of times, will sometimes execute the following statement
// SCONT will run only if a SOMETIMES above in the same scope also ran
#define SOMETIMES() float SOMETIMES_RAND = rand() % 1000; bool SOMETIMES_TRACKER_VAR = SOMETIMES_RAND == 1; if (SOMETIMES_TRACKER_VAR)
#define SCONT() if (SOMETIMES_TRACKER_VAR)

// TIMING MACROS
#define TIMER_INIT() std::chrono::time_point<std::chrono::system_clock> TIMER_START_TIME, TIMER_END_TIME; std::chrono::duration<double> TIMER_DIFF; double TIMER_DURATION;

#define TIMER_START() TIMER_START_TIME = std::chrono::system_clock::now()
#define TIMER_END()   TIMER_END_TIME = std::chrono::system_clock::now(); TIMER_DIFF = TIMER_END_TIME - TIMER_START_TIME; TIMER_DURATION = TIMER_DIFF.count()

#define TIME(exp) TIMER_START(); exp TIMER_END(); PRINTD(TIMER_DURATION);

// PROGRESS BAR MACROS

// From https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
#define PB_FILLSTR  "============================================================"
#define PB_EMPTYSTR "                                                            "
#define PBWIDTH 60

namespace progressBar {
    inline void printProgress(float progress) {
        int lpad = (int) (progress * PBWIDTH);
        int rpad = PBWIDTH  - lpad;
        if (progress < ((float) 1 / PBWIDTH)) {
            printf("[%s]", PB_EMPTYSTR);
        } else if (progress > ((float) (PBWIDTH-1) / PBWIDTH)) {
            printf("[%s]", PB_FILLSTR);
        } else {
            printf("[%.*s>%.*s]", lpad - 1, PB_FILLSTR, rpad, PB_EMPTYSTR);
        }
        fflush(stdout);
    }

    inline void printDuration(int seconds) {
        int hours = seconds / 3600;
        int mins = (seconds % 3600) / 60;
        int secs = (seconds % 60);
        printf("%02d:%02d:%02d", hours, mins, secs);
    }
}

#define PB_DECL() char PB_DESCRIPTION[256]; std::chrono::time_point<std::chrono::system_clock> PB_START_TIME, PB_CUR_TIME; std::chrono::duration<double> PB_DIFF; double PB_DURATION = 0
#define PB_STARTD(description_fmt, ...) sprintf(PB_DESCRIPTION, description_fmt, ## __VA_ARGS__); PB_START_TIME = std::chrono::system_clock::now(); PB_PROGRESS(0)

#define PB_START(description_fmt, ...) char PB_DESCRIPTION[256]; sprintf(PB_DESCRIPTION, description_fmt, ## __VA_ARGS__); std::chrono::time_point<std::chrono::system_clock> PB_START_TIME = std::chrono::system_clock::now(), PB_CUR_TIME; std::chrono::duration<double> PB_DIFF; double PB_DURATION = 0; PB_PROGRESS(0)
#define PB_PROGRESS(progress) PB_CUR_TIME = std::chrono::system_clock::now(); PB_DIFF = PB_CUR_TIME - PB_START_TIME; PB_DURATION = PB_DIFF.count(); \
                              printf("\33[2K\r%s: %.2f%% ", PB_DESCRIPTION, (float) progress * 100); \
                              progressBar::printProgress(progress); \
                              printf(" Elapsed: "); progressBar::printDuration(PB_DURATION); printf(", ETA: "); progressBar::printDuration( ((PB_DURATION / (float) (progress))) - PB_DURATION ); fflush(stdout)
#define PB_END() printf("\33[2K\r%s: %.2f%% ", PB_DESCRIPTION, 100.0); progressBar::printProgress(1); printf(" Took: "); progressBar::printDuration(PB_DURATION); printf("\n")

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

    inline void write_vec2f(FILE* file, const VEC2F& vec) {
        if (sizeof(Real) == sizeof(double))
        {
            fwrite((void*)&vec[0], sizeof(double), 1, file);
            fwrite((void*)&vec[1], sizeof(double), 1, file);
        }
        else
        {
            double elementD[2];
            elementD[0] = vec[0];
            elementD[1] = vec[1];
            fwrite((void*)&elementD[0], sizeof(double), 1, file);
            fwrite((void*)&elementD[1], sizeof(double), 1, file);
        }
    }

    inline void read_vec2f(FILE* file, VEC2F& vec) {
        // make sure it reads into a double
        double elementD[2];
        fread((void*)&elementD[0], sizeof(double), 1, file);
        fread((void*)&elementD[1], sizeof(double), 1, file);

        // let the assign resolve the type issues
        vec[0] = elementD[0];
        vec[1] = elementD[1];
    }

}



#endif
