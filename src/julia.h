#ifndef JULIA_H
#define JULIA_H

#include "mesh.h"
#include "field.h"
#include "Quaternion/QUATERNION.h"
#include "Quaternion/POLYNOMIAL_4D.h"

class QuatToQuatFn {
public:
    virtual QUATERNION getFieldValue(QUATERNION q) const = 0;

    virtual QUATERNION operator()(QUATERNION q) const {
        return getFieldValue(q);
    }
};

class SimpleJuliaQuat: public QuatToQuatFn {
public:
    QUATERNION c;
    SimpleJuliaQuat(QUATERNION c): c(c) {}

    virtual QUATERNION getFieldValue(QUATERNION q) const override {
        return (q * q) + c;
    }
};

class RationalQuatPoly: public QuatToQuatFn {
public:
    POLYNOMIAL_4D topPolynomial;
    bool hasBottomPolynomial;
    POLYNOMIAL_4D bottomPolynomial;

    RationalQuatPoly(POLYNOMIAL_4D topPolynomial) : topPolynomial(topPolynomial), hasBottomPolynomial(false) {}
    RationalQuatPoly(POLYNOMIAL_4D topPolynomial, POLYNOMIAL_4D bottomPolynomial) : topPolynomial(topPolynomial), hasBottomPolynomial(true), bottomPolynomial(bottomPolynomial) {}

    virtual QUATERNION getFieldValue(QUATERNION q) const override {
        QUATERNION out = topPolynomial.evaluateScaledPowerFactored(q);
        if (hasBottomPolynomial) {
            QUATERNION bottomEval = bottomPolynomial.evaluateScaledPowerFactored(q);
            out = (out / bottomEval);
        }
        return out;
    }

};

class DistanceGuidedQuatMap: public FieldFunction3D {
public:
    Grid3D* distanceField;
    QuatToQuatFn* p;

    Real c;
    int maxIterations;
    Real escape;
    Real fitScale;

    ArrayGrid3D* debugGrid;

public:
    DistanceGuidedQuatMap(Grid3D* distanceField,  QuatToQuatFn* p, Real c = 300, int maxIterations = 3, Real escape = 20, Real fitScale = 1):
        distanceField(distanceField), p(p), c(c), maxIterations(maxIterations), escape(escape), fitScale(fitScale) {}

    Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION iterate(pos[0], pos[1], pos[2], 0);
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        while (magnitude < escape && totalIterations < maxIterations) {
            // First lookup the radius we're going to project to and save the original
            const Real distance = (*distanceField)(VEC3F(iterate[0], iterate[1], iterate[2]));
            Real radius = exp(c * distance);

            QUATERNION original = iterate;
            VEC3F iterateV3(iterate[0], iterate[1], iterate[2]);

            // If we know it'll escape we can skip quaternion multiplication evaluation (unless we're scaling)
            if (fitScale == 1 && radius > escape) {
                magnitude = radius;
                totalIterations++;
                break;
            }

            iterate = p->getFieldValue(iterate);

            // If quaternion multiplication fails, revert back to original
            if (iterate.anyNans()) {
                if (debugGrid) (*debugGrid)(pos) = 2;
                iterate = original;
            }

            // Try to normalize iterate, might produce infs or nans if magnitude
            // is too small or zero if magnitude is too large
            QUATERNION normedIterate = iterate;
            normedIterate.normalize();

            // Scale radius (accounting for fitScale)
            // Fitscale is a scalar parameter to smoothly (ish) interpolate between the original P and
            // the mapped version, unless you're doing that it should be always set to 1
            radius = exp((1 - fitScale) * log(iterate.magnitude()) + (fitScale * log(radius)));

            // Now that we know the radius we can break if it's too big:
            if (radius > escape) {
                if (debugGrid) (*debugGrid)(pos) = -1;
                magnitude = radius;
                totalIterations++;
                break;
            }

            bool tooSmall = normedIterate.anyNans();
            if (tooSmall) {
                QUATERNION origNorm = original;
                origNorm.normalize();
                if (origNorm.anyNans()) {
                    // Need to put it at radius but have no direction info
                    iterate = QUATERNION(1,0,0,0) * radius;
                    if (debugGrid) (*debugGrid)(pos) = 4;

                    // This should never happen. If it does, your polynomial is
                    // probably returing wayyy too large values that double
                    // precision can't hold their inverse
                } else {
                    if (debugGrid) (*debugGrid)(pos) = 3;
                    iterate = origNorm * radius;
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

    void sampleFit(string filename, int res, Real fitScale) {
        ofstream out(filename);

        VEC3F min = distanceField->mapBox.min();
        VEC3F max = distanceField->mapBox.max();
        VEC3F inc = distanceField->mapBox.span()/res;

        PB_START("Sampling fit on %d x %d x %d grid...", res, res, res);
        PB_PROGRESS(0);

        for (Real x = min[0]; x < max[0]; x+=inc[0]) {
            for (Real y = min[1]; y < max[1]; y+=inc[1]) {
                for (Real z = min[2]; z < max[2]; z+=inc[2]) {
                    QUATERNION samplePoint(x, y, z, 0);

                    const Real distance = (*distanceField)(VEC3F(samplePoint[0], samplePoint[1], samplePoint[2]));
                    const Real targetMag = exp( c * distance);

                    bool nanOut = false;
                    Real actualMag = p->getFieldValue(samplePoint).magnitude();

                    if (::isnan(actualMag)) {
                        nanOut = true;
                        actualMag = samplePoint.magnitude();
                    }

                    // const Real adjustedMag = ((1 - fitScale) * actualMag) + (fitScale * targetMag);
                    const Real adjustedMag = exp((1 - fitScale) * log(actualMag)) + (fitScale * log(targetMag));

                    if (!::isnan(distance)) out << targetMag << ", " << adjustedMag << ", " << (nanOut? 1 : 0) << endl;
                }
            }
            PB_PROGRESS((x - min[0]) / (max[0] - min[0]));
        }

        PB_END();

        out.close();

    }

};

// =============== INSPECTION FIELDS =======================

class QuatQuatRotField: public FieldFunction3D {
public:
    QuatToQuatFn *func;
    int i;
    QuatQuatRotField(QuatToQuatFn* func, int i): func(func), i(i) {}

    virtual Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION input(pos[0], pos[1], pos[2], 0);
        QUATERNION transformed = func->getFieldValue(input);

        // PRINTV4(transformed);
        // PRINTE(transformed.magnitude());

        QUATERNION normalized = transformed;

        normalized.normalize();
        // PRINTV4(normalized);
        //
        // PRINTE(normalized.magnitude());
        //
        // if(normalized.magnitude() == 0) exit(1);

        return normalized[i];
    }

};

class QuatQuatMagField: public FieldFunction3D {
public:
    QuatToQuatFn *func;
    QuatQuatMagField(QuatToQuatFn* func): func(func) {}

    virtual Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION input(pos[0], pos[1], pos[2], 0);
        QUATERNION transformed = func->getFieldValue(input);
        return transformed.magnitude();
    }

};

#endif
