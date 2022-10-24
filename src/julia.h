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

    virtual void writeCSVPairs(string filename, uint wRes, uint xRes, uint yRes, uint zRes, QUATERNION fieldMin, QUATERNION fieldMax) {
        ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;

        for (uint i = 0; i < wRes; ++i) {
            for (uint j = 0; j < xRes; ++j) {
                for (uint k = 0; k < yRes; ++k) {
                    for (uint l = 0; l < zRes; ++l) {
                        QUATERNION sampleOffset = fieldMax - fieldMin;
                        sampleOffset.w() *= ((float) i / wRes);
                        sampleOffset.x() *= ((float) j / xRes);
                        sampleOffset.y() *= ((float) k / yRes);
                        sampleOffset.z() *= ((float) l / zRes);

                        QUATERNION samplePoint = fieldMin + sampleOffset;
                        QUATERNION value = getFieldValue(samplePoint);

                        out <<
                            samplePoint.w() << "," <<
                            samplePoint.x() << "," <<
                            samplePoint.y() << "," <<
                            samplePoint.z() << "," <<
                            value.w()       << "," <<
                            value.x()       << "," <<
                            value.y()       << "," <<
                            value.z()       << "," <<
                        endl;


                    }
                }
            }
        }

        printf("Wrote %d x %d x %d x %d field (%d values) to %s\n", wRes, xRes, yRes, zRes, (wRes * xRes * yRes * zRes), filename.c_str());

        out.close();

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

class QuaternionJuliaSet: public FieldFunction3D {
public:
    QuatToQuatFn* p;
    int maxIterations;
    Real escape;

public:
    QuaternionJuliaSet(QuatToQuatFn* p, int maxIterations = 3, Real escape = 20):
        p(p), maxIterations(maxIterations), escape(escape) {}

    Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION iterate(pos[0], pos[1], pos[2], 0);
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        while (magnitude < escape && totalIterations < maxIterations) {
            iterate = p->getFieldValue(iterate);
            magnitude = iterate.magnitude();

            totalIterations++;
        }

        return log(magnitude);
    }

};

class DistanceGuidedQuatFn: public QuatToQuatFn {
public:
    Grid3D* distanceField;
    QuatToQuatFn* p;

    FieldFunction3D* a;
    FieldFunction3D* b;

    // When a or b are constant fields, we could use a
    // ConstantFunction3D, but to save the overhead of
    // a lookup we have a special mode for constant values
    bool hasConstantA;
    bool hasConstantB;
    Real constantA;
    Real constantB;

    Real fitScale;

public:
    DistanceGuidedQuatFn(Grid3D* distanceField,  QuatToQuatFn* p, Real a = 300, Real b = 0, Real fitScale = 1):
        distanceField(distanceField), p(p), hasConstantA(true), hasConstantB(true), constantA(a), constantB(b), fitScale(fitScale) {}

    DistanceGuidedQuatFn(Grid3D* distanceField,  QuatToQuatFn* p, FieldFunction3D* a, FieldFunction3D* b, Real fitScale = 1):
        distanceField(distanceField), p(p), a(a), b(b), hasConstantA(false), hasConstantB(false), fitScale(fitScale) {}

    QUATERNION getFieldValue(QUATERNION q) const override {
        // First lookup the radius we're going to project to and save the original

        QUATERNION original = q;
        VEC3F iterateV3(q[0], q[1], q[2]);
        const Real distance = (*distanceField)(iterateV3);

        Real radius = exp(
            hasConstantA ? constantA : a->getFieldValue(iterateV3) *
            (distance - (hasConstantB ? constantB : b->getFieldValue(iterateV3))));

        q = p->getFieldValue(q);

        // If quaternion multiplication fails, revert back to original
        if (q.anyNans()) {
            q = original;
        }

        // Try to normalize q, might produce infs or nans if magnitude
        // is too small or zero if magnitude is too large
        QUATERNION normedIterate = q;
        normedIterate.normalize();

        // Scale radius (accounting for fitScale)
        // Fitscale is a scalar parameter to smoothly (ish) interpolate between the original P and
        // the mapped version, so unless you're doing that it should be always set to 1
        radius = exp((1 - fitScale) * log(q.magnitude()) + (fitScale * log(radius)));

        bool tooSmall = normedIterate.anyNans();
        if (tooSmall) {
            QUATERNION origNorm = original;
            origNorm.normalize();
            if (origNorm.anyNans()) {
                // This should never happen. If it does, your polynomial is
                // probably returing wayyy too large values such that double
                // precision can't hold their inverse. In this case, we need
                // to put the iterate at the specified radius but we don't
                // have any direction info, so we just do the following:

                q = QUATERNION(1,0,0,0) * radius;
            } else {
                q = origNorm * radius;
            }
        } else {
            q = normedIterate;
            q *= radius;
        }

        return q;
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

        QUATERNION normalized = transformed;

        normalized.normalize();

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
