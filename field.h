#ifndef FIELD_H
#define FIELD_H

#include <fstream>

#include "SETTINGS.h"

class Field3D {
public:
    Real* values;
    int xRes, yRes, zRes;

    // Create empty (not zeroed) field with given resolution
    Field3D(int xRes, int yRes, int zRes): xRes(xRes), yRes(yRes), zRes(zRes) {
        values = new Real[xRes * yRes * zRes];
    }

    // Create empty (not zeroed) field with given resolution
    Field3D(VEC3I resolution): Field3D(resolution[0], resolution[1], resolution[2]) {}

    // Destructor
    ~Field3D() {
        delete values;
    }

    // Access value based on integer indices
    Real& at(int x, int y, int z) {
        return values[(z * yRes + y) * xRes + x];
    }

    // Access value from VEC3I position
    Real& at(VEC3I pos) {
        return at(pos[0], pos[1], pos[2]);
    }

    // Access values in their C-contiguous order
    Real& at(int i) {
        return values[i];
    }

    // Create field from scalar function by sampling it on a regular grid
    Field3D(VEC3I resolution, VEC3F functionMin, VEC3F functionMax, Real (*fieldFunction)(VEC3F)): Field3D(resolution) {

        VEC3F gridResF(xRes, yRes, zRes);

        for (int i = 0; i < xRes; i++) {
            for (int j = 0; j < yRes; j++) {
                for (int k = 0; k < zRes; k++) {
                    VEC3F gridPointF(i, j, k);
                    VEC3F fieldDelta = functionMax - functionMin;

                    VEC3F samplePoint = functionMin + (gridPointF.cwiseQuotient(gridResF - VEC3F(1,1,1)).cwiseProduct(fieldDelta));

                    this->at(i, j, k) = fieldFunction(samplePoint);
                }
            }
        }

    }


    void writeCSV(std::string filename){
        std::ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;

        for (int i = 0; i < xRes; ++i) {
            for (int j = 0; j < yRes; ++j) {
                for (int k = 0; k < zRes; ++k) {
                    out << i << ", " << j << ", " << k << ", " << at(i, j, k) << std::endl;
                }
            }
        }

        out.close();

    }

};

#endif
