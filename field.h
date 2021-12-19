#ifndef FIELD_H
#define FIELD_H

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <queue>

#include "SETTINGS.h"

using namespace std;

class AABB: public AlignedBox<Real, 3> {
public:
    using AlignedBox<Real, 3>::AlignedBox;
    AABB(VEC3F c1, VEC3F c2): AlignedBox<Real, 3>(c1.cwiseMin(c2), c1.cwiseMax(c2)) {}

    VEC3F span() const {
        return max()-min();
    }

    VEC3F clamp(VEC3F pos) const {
        return pos.cwiseMax(min()).cwiseMin(max());
    }

    void setCenter(VEC3F newCenter) {
        VEC3F offset = newCenter - this->center();
        max() += offset;
        min() += offset;
    }
};

class FieldFunction3D {
private:
    Real (*fieldFunction)(VEC3F pos);
public:
    FieldFunction3D(Real (*fieldFunction)(VEC3F pos)):fieldFunction(fieldFunction) {}

    FieldFunction3D():fieldFunction(nullptr) {} // Only to be used by subclasses

    virtual Real getFieldValue(const VEC3F& pos) const {
        return fieldFunction(pos);
    }

    virtual Real operator()(const VEC3F& pos) const {
        return getFieldValue(pos);
    }
};

class Grid3D: public FieldFunction3D {
public:
    uint xRes, yRes, zRes;
    bool supportsNonIntegerIndices = false;

    // Bounds (aka center + lengths) for mapping into this grid
    // as a field function
    AABB mapBox;
    bool hasMapBox = false;

    virtual uint totalCells() const {
        return xRes * yRes * zRes;
    }

    virtual Real get(uint x, uint y, uint z) const = 0;

    virtual Real getf(Real x, Real y, Real z) const {
        (void) x; (void) y; (void) z; // Suppress unused argument warning
        printf("This grid doesn't support non-integer indices!\n");
        exit(1);
    }

    virtual Real get(VEC3I pos) const {
        return get(pos[0], pos[1], pos[2]);
    }

    virtual Real getf(VEC3F pos) const {
        return getf(pos[0], pos[1], pos[2]);
    }

    virtual void setMapBox(AABB box) {
        mapBox = box;
        hasMapBox = true;
    }

    virtual Real getFieldValue(const VEC3F& pos) const override {
        if (!hasMapBox) {
            printf("Attempting getFieldValue on a Grid without a mapBox!\n");
            exit(1);
        }

        VEC3F samplePoint = (pos - mapBox.min()).cwiseQuotient(mapBox.span());
        samplePoint = samplePoint.cwiseMax(VEC3F(0,0,0)).cwiseMin(VEC3F(1,1,1));

        const VEC3F indices = samplePoint.cwiseProduct(VEC3F(xRes-1, yRes-1, zRes-1));

        if (supportsNonIntegerIndices) {
            return getf(indices);
        } else {
            return get(indices.cast<int>());
        }
    }

    virtual void writeCSV(string filename) {
        ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;

        for (uint i = 0; i < xRes; ++i) {
            for (uint j = 0; j < yRes; ++j) {
                for (uint k = 0; k < zRes; ++k) {
                    out << i << ", " << j << ", " << k << ", " << get(i, j, k) << endl;
                }
            }
        }

        printf("Wrote %d x %d x % d field (%d values) to %s\n", xRes, yRes, zRes, (xRes * yRes * zRes), filename.c_str());

        out.close();

    }

    // Writes to F3D file using the field bounds if the field has them, otherwise using
    // the resolution ofthe grid.
    void writeF3D(string filename, bool verbose = false) const {
        if (hasMapBox) {
            writeF3D(filename, mapBox, verbose);
        } else {
            writeF3D(filename, AABB(VEC3F(0,0,0), VEC3F(xRes, yRes, zRes)), verbose);
        }
    }

    void writeF3D(string filename, AABB bounds, bool verbose = false) const {
        FILE* file = fopen(filename.c_str(), "wb");

        if (file == NULL) {
            PRINT("Failed to write F3D: file open failed!");
            exit(0);
        }

        if (verbose) {
            printf("Writing %d x %d x %d field to %s... ", xRes, yRes, zRes, filename.c_str());
            fflush(stdout);
        }

        // write dimensions
        fwrite((void*)&xRes, sizeof(int), 1, file);
        fwrite((void*)&yRes, sizeof(int), 1, file);
        fwrite((void*)&zRes, sizeof(int), 1, file);

        MyEigen::write_vec3f(file, bounds.center());
        MyEigen::write_vec3f(file, bounds.span());

        const int totalCells = xRes*yRes*zRes;

        if (totalCells <= 0)
            return;

        // write data
        for (uint i = 0; i < xRes; ++i) {
            for (uint j = 0; j < yRes; ++j) {
                for (uint k = 0; k < zRes; ++k) {
                    const Real val = get(i, j, k);
                    double out;
                    if (sizeof(Real) != sizeof(double)) {
                        out = (double) val;
                    } else {
                        out = val;
                    }

                    fwrite((void*) (&out), sizeof(double), 1, file);

                }
            }

            if (verbose && i % 10 == 0) {
                printf("\rWriting %d x %d x %d field to %s... %.2f%%", xRes, yRes, zRes, filename.c_str(), (Real) (100*i) / xRes);
                fflush(stdout);
            }
        }

        if (verbose) printf("\n");
    }
};

class ArrayGrid3D: public Grid3D {
private:
    Real* values;
public:

    // Create empty (not zeroed) field with given resolution
    ArrayGrid3D(uint xRes, uint yRes, uint zRes) {
        this->xRes = xRes;
        this->yRes = yRes;
        this->zRes = zRes;
        values = new Real[xRes * yRes * zRes];
    }

    // Create empty (not zeroed) field with given resolution
    ArrayGrid3D(VEC3I resolution): ArrayGrid3D(resolution[0], resolution[1], resolution[2]) {}

    // Read ArrayGrid3D from F3D
    ArrayGrid3D(string filename, string format = "f3d", bool verbose = false) {

        if (format == "f3d") {
            FILE* file = fopen(filename.c_str(), "rb");
            if (file == NULL) {
                PRINT("Failed to read F3D: file open failed!");
                exit(0);
            }

            int xRes, yRes, zRes;
            VEC3F center, lengths;

            // read dimensions
            fread((void*)&xRes, sizeof(int), 1, file);
            fread((void*)&yRes, sizeof(int), 1, file);
            fread((void*)&zRes, sizeof(int), 1, file);

            MyEigen::read_vec3f(file, center);
            MyEigen::read_vec3f(file, lengths);

            this->xRes = xRes;
            this->yRes = yRes;
            this->zRes = zRes;

            try {
                values = new Real[xRes * yRes * zRes];
            } catch(bad_alloc& exc) {
                printf("Failed to allocate %.2f MB for ArrayGrid3D read from file!\n", (xRes * yRes * zRes * sizeof(Real)) / pow(2.0,20.0));
                exit(0);
            }

            if (verbose) {
                printf("Reading %d x %d x %d field from %s... ", xRes, yRes, zRes, filename.c_str());
                fflush(stdout);
            }

            AABB bounds((center - lengths/2), (center + lengths/2));
            setMapBox(bounds);

            const int totalCells = xRes * yRes * zRes;
            // always read in as a double
            if (sizeof(Real) != sizeof(double)) {
                double* dataDouble = new double[totalCells];
                fread((void*)dataDouble, sizeof(double), totalCells, file);

                for (int x = 0; x < totalCells; x++)
                    values[x] = dataDouble[x];

                if (verbose) printf("\n");

                delete[] dataDouble;
            } else fread((void*)values, sizeof(Real), totalCells, file);

            if (verbose) {
                printf("done.\n");
            }

        } else {
            PRINT("CSV import not implemented yet!");
            exit(1);
        }

    }
    // Destructor
    ~ArrayGrid3D() {
        delete values;
    }

    // Access value based on integer indices
    Real get(uint x, uint y, uint z) const override {
        return values[(z * yRes + y) * xRes + x];
    }

    // Access value directly (allows setting)
    Real& at(uint x, uint y, uint z) {
        return values[(z * yRes + y) * xRes + x];
    }

    // Access value directly in C-style array (allows setting)
    Real& operator[](size_t x) {
        return values[x];
    }

    // Create field from scalar function by sampling it on a regular grid
    ArrayGrid3D(uint xRes, uint yRes, uint zRes, VEC3F functionMin, VEC3F functionMax, FieldFunction3D *fieldFunction):ArrayGrid3D(xRes, yRes, zRes){

        VEC3F gridResF(xRes, yRes, zRes);

        for (uint i = 0; i < xRes; i++) {
            for (uint j = 0; j < yRes; j++) {
                for (uint k = 0; k < zRes; k++) {
                    VEC3F gridPointF(i, j, k);
                    VEC3F fieldDelta = functionMax - functionMin;

                    VEC3F samplePoint = functionMin + (gridPointF.cwiseQuotient(gridResF - VEC3F(1,1,1)).cwiseProduct(fieldDelta));

                    this->at(i, j, k) = fieldFunction->getFieldValue(samplePoint);
                }
            }
        }

    }


};

class VirtualGrid3D: public Grid3D {
private:
    FieldFunction3D *fieldFunction;
    VEC3F functionMin, functionMax;

    VEC3F getSamplePoint(Real x, Real y, Real z) const {
        VEC3F gridPointF(x, y, z);
        VEC3F gridResF(xRes, yRes, zRes);
        VEC3F fieldDelta = functionMax - functionMin;
        VEC3F samplePoint = functionMin + (gridPointF.cwiseQuotient(gridResF - VEC3F(1,1,1)).cwiseProduct(fieldDelta));

        return samplePoint;
    }
public:

    VirtualGrid3D(uint xRes, uint yRes, uint zRes, VEC3F functionMin, VEC3F functionMax,  FieldFunction3D *fieldFunction):
        fieldFunction(fieldFunction),
        functionMin(functionMin),
        functionMax(functionMax) {
            this->xRes = xRes;
            this->yRes = yRes;
            this->zRes = zRes;

            this->setMapBox(AABB(functionMin, functionMax));

            this->supportsNonIntegerIndices = true;
        }

    // Virtually (shallow) downsample an existing grid
    // VirtualGrid3D(uint xRes, uint yRes, uint zRes, Grid3D *other) {
    // TODO implement
    // }

    virtual Real get(uint x, uint y, uint z) const override {
        return getf(x, y, z);
    }

    virtual Real getf(Real x, Real y, Real z) const override {
        return fieldFunction->getFieldValue(getSamplePoint(x, y, z));
    }
};

// Hash function for Eigen matrix and vector.
// From https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
template<typename T> struct matrix_hash : unary_function<T, size_t> {
    size_t operator()(T const& matrix) const {
        // Note that it is oblivious to the storage order of Eigen matrix (column- or
        // row-major). It will give you the same hash value for two different matrices if they
        // are the transpose of each other in different storage order.
        size_t seed = 0;
        for (long i = 0; i < matrix.size(); ++i) { // For some reason Eigen size is not unsigned?!
            auto elem = *(matrix.data() + i);
            seed ^= hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


class VirtualGrid3DCached: public VirtualGrid3D {
protected:
    unordered_map<VEC3F, Real, matrix_hash<VEC3F>> map;
public:
    using VirtualGrid3D::VirtualGrid3D;

    virtual Real get(uint x, uint y, uint z) {
        return getf(x,y,z);
    }

    virtual Real getf(Real x, Real y, Real z) {
        VEC3F key(x,y,z);

        auto search = map.find(key);
        if (search != map.end()) {
            return search->second;
        }

        Real result = VirtualGrid3D::get(x,y,z);
        map[key] = result;
        return result;
    }

};

class VirtualGrid3DLimitedCache: public VirtualGrid3DCached {
private:
    size_t maxSize;
    queue<VEC3F> cacheQueue;

public:
    // Instantiates a VirtualGrid3D with a limited-size cache. When additional
    // items are inserted into the cache (beyond the capacity), the cache will
    // forget the item that was least recently inserted. If capacity -1 is
    // specified (default), it defaults to a size equal to three XY slices
    // through the field, which is suited for marching cubes.
    VirtualGrid3DLimitedCache(uint xRes, uint yRes, uint zRes, VEC3F functionMin, VEC3F functionMax,  FieldFunction3D *fieldFunction, int capacity = -1):
        VirtualGrid3DCached(xRes, yRes, zRes, functionMin, functionMax, fieldFunction) {
            PRINTV3(functionMin);
            PRINTV3(functionMax);
            if (capacity == -1) {
                maxSize = xRes * yRes * 3;
            } else {
                maxSize = capacity;
            }
        }

    virtual Real getf(Real x, Real y, Real z) override {
        VEC3F key(x,y,z);

        auto search = map.find(key);
        if (search != map.end()) {
            PRINT("CACHE HIT\n");
            return search->second;
        }

        Real result = VirtualGrid3D::get(x,y,z);
        map[key] = result;

        if (cacheQueue.size() >= maxSize) {
            map.erase(cacheQueue.front());
            cacheQueue.pop();
        }

        cacheQueue.push(key);


        return result;
    }
};


class InterpolationGrid: public Grid3D {
private:
    Real interpolate(Real x0, Real x1, Real d) const {
        switch (mode) {
        case LINEAR:
            return ((1 - d) * x0) + (d * x1);
        case SMOOTHSTEP:
            d = (3 * d * d) - (2 * d * d * d);
            return ((1 - d) * x0) + (d * x1);
        }
        assert(false);
        return -1;
    }


public:
    Grid3D* baseGrid;

    enum INTERPOLATION_MODE {
        LINEAR,
        SMOOTHSTEP
    };

    INTERPOLATION_MODE mode;

    InterpolationGrid(Grid3D* baseGrid, INTERPOLATION_MODE mode = LINEAR) {
        this->baseGrid = baseGrid;
        xRes = baseGrid->xRes;
        yRes = baseGrid->yRes;
        zRes = baseGrid->zRes;
        this->mode = mode;
        this->supportsNonIntegerIndices = true;

        if (baseGrid->hasMapBox) this->setMapBox(baseGrid->mapBox);

        if (baseGrid->supportsNonIntegerIndices) {
            PRINT("Warning: laying an InterpolationGrid over a grid which already supports non-integer indices!");
        }
    }

    virtual Real get(uint x, uint y, uint z) const override {
        return baseGrid->get(x, y, z);
    }

    virtual Real getf(Real x, Real y, Real z) const override {
        // "Trilinear" interpolation with whatever technique we select

        const uint x0 = floor(x);
        const uint y0 = floor(y);
        const uint z0 = floor(z);

        const uint x1 = ceil(x);
        const uint y1 = ceil(y);
        const uint z1 = ceil(z);

        const Real xd = (x - x0) / ((Real) x1 / x0);
        const Real yd = (y - y0) / ((Real) y1 / y0);
        const Real zd = (z - z0) / ((Real) z1 / z0);

        // First grab 3D surroundings...
        const Real c000 = baseGrid->get(x0, y0, z0);
        const Real c001 = baseGrid->get(x0, y0, z1);
        const Real c010 = baseGrid->get(x0, y1, z0);
        const Real c011 = baseGrid->get(x0, y1, z1);
        const Real c100 = baseGrid->get(x1, y0, z0);
        const Real c101 = baseGrid->get(x1, y0, z1);
        const Real c110 = baseGrid->get(x1, y1, z0);
        const Real c111 = baseGrid->get(x1, y1, z1);

        // Now create 2D interpolated slice...
        const Real c00 = interpolate(c000, c100, xd);
        const Real c01 = interpolate(c001, c101, xd);
        const Real c10 = interpolate(c010, c110, xd);
        const Real c11 = interpolate(c011, c111, xd);

        // Extract 1D interpolated line...
        const Real c0 = interpolate(c00, c10, yd);
        const Real c1 = interpolate(c01, c11, yd);

        // And grab 0D point.
        return interpolate(c0, c1, zd);
    }



};



#endif
