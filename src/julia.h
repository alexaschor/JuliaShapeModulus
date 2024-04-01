#ifndef JULIA_H
#define JULIA_H

#include "SETTINGS.h"
#include "mesh.h"
#include "field.h"
#include "Quaternion/QUATERNION.h"
#include "Quaternion/POLYNOMIAL_4D.h"
#include "PerlinNoise/PerlinNoise.h"

class QuatMap {
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



class SimpleJuliaQuat: public QuatMap {
public:
    QUATERNION c;
    SimpleJuliaQuat(QUATERNION c): c(c) {}

    virtual QUATERNION getFieldValue(QUATERNION q) const override {
        return (q * q) + c;
    }
};

class RationalQuatPoly: public QuatMap {
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
    QuatMap* p;
    int maxIterations;
    Real escape;

public:
    QuaternionJuliaSet(QuatMap* p, int maxIterations = 3, Real escape = 20):
        p(p), maxIterations(maxIterations), escape(escape) {}

    Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION iterate(pos[0], pos[1], pos[2], 0);
        Real magnitude = iterate.magnitude();
        int totalIterations = 0;

        while (magnitude < escape && totalIterations < maxIterations) {
            QUATERNION newIterate = p->getFieldValue(iterate);
            iterate = newIterate;
            magnitude = iterate.magnitude();
            totalIterations++;
        }

        Real out = log(magnitude);
        return out;
    }

};

class DistanceGuidedQuatFn: public QuatMap {
public:
    Grid3D* distanceField;
    QuatMap* p;

    FieldFunction3D* a;
    FieldFunction3D* b;

    // When a or b are constant fields, we could use a
    // ConstantFunction3D, but to save the overhead of
    // a lookup we have a special mode for constant values
    bool hasConstantA;
    bool hasConstantB;
    Real constantA;
    Real constantB;

public:
    DistanceGuidedQuatFn(Grid3D* distanceField,  QuatMap* p, Real a = 300, Real b = 0):
        distanceField(distanceField), p(p), hasConstantA(true), hasConstantB(true), constantA(a), constantB(b) {}

    DistanceGuidedQuatFn(Grid3D* distanceField,  QuatMap* p, FieldFunction3D* a, FieldFunction3D* b):
        distanceField(distanceField), p(p), a(a), b(b), hasConstantA(false), hasConstantB(false) {}

    QUATERNION getFieldValue(QUATERNION q) const override {

        // First lookup the radius we're going to project to and save the original
        QUATERNION original = q;
        VEC3F iterateV3(q[0], q[1], q[2]);
        const Real distance = (*distanceField)(iterateV3);

        Real aValue = (hasConstantA ? constantA : a->getFieldValue(iterateV3));
        Real bValue = (hasConstantB ? constantB : b->getFieldValue(iterateV3)) ;

        Real radius = exp( aValue * (distance - bValue ));

        q = p->getFieldValue(q);

        // If quaternion multiplication fails, revert back to original
        if (q.anyNans()) {
            q = original;
        }

        // Try to normalize q, might produce infs or nans if magnitude
        // is too small or zero if magnitude is too large
        QUATERNION normedIterate = q;
        normedIterate.normalize();

        bool tooSmall = (normedIterate.anyNans() || normedIterate.magnitude() == 0);
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

class R3Map {
public:
    virtual VEC3F getFieldValue(const VEC3F& q) const = 0;

    virtual VEC3F operator()(const VEC3F& q) const {
        return getFieldValue(q);
    }

    virtual void writeCSVPairs(string filename, uint xRes, uint yRes, uint zRes, VEC3F fieldMin, VEC3F fieldMax) {
        ofstream out;
        out.open(filename);
        if (out.is_open() == false)
            return;

        for (uint j = 0; j < xRes; ++j) {
            for (uint k = 0; k < yRes; ++k) {
                for (uint l = 0; l < zRes; ++l) {
                    VEC3F sampleOffset = fieldMax - fieldMin;
                    sampleOffset.x() *= ((float) j / xRes);
                    sampleOffset.y() *= ((float) k / yRes);
                    sampleOffset.z() *= ((float) l / zRes);

                    VEC3F samplePoint = fieldMin + sampleOffset;
                    VEC3F value = getFieldValue(samplePoint);

                    out <<
                        samplePoint.x() << "," <<
                        samplePoint.y() << "," <<
                        samplePoint.z() << "," <<
                        value.x()       << "," <<
                        value.y()       << "," <<
                        value.z()       << "," <<
                        endl;


                }
            }
        }

        printf("Wrote %d x %d x %d field (%d values) to %s\n", xRes, yRes, zRes, (xRes * yRes * zRes), filename.c_str());

        out.close();

    }
};

class R3JuliaSet: public FieldFunction3D {
public:
    R3Map* m;
    int maxIterations;
    Real escape;

public:
    R3JuliaSet(R3Map* m, int maxIterations = 3, Real escape = 20):
        m(m), maxIterations(maxIterations), escape(escape) {}

    Real getFieldValue(const VEC3F& pos) const override {
        VEC3F iterate(pos);
        Real magnitude = iterate.norm();
        int totalIterations = 0;

        while (magnitude < escape && totalIterations < maxIterations) {
            VEC3F newIterate = m->getFieldValue(iterate);
            iterate = newIterate;
            magnitude = iterate.norm();
            totalIterations++;
        }

        Real out = log(magnitude);
        return out;
    }

};

class VersorModulusR3Map: public R3Map {
public:
    R3Map* versor;
    FieldFunction3D* modulus;

    VersorModulusR3Map(R3Map* versor, FieldFunction3D* modulus): versor(versor), modulus(modulus) {};

    VEC3F getFieldValue(const VEC3F& pos) const override {
        return (*versor)(pos) * (*modulus)(pos);
    }
};


class ShapeModulus: public FieldFunction3D {
public:
    Grid3D* distanceField;

    FieldFunction3D* a;
    FieldFunction3D* b;

    // When a or b are constant fields, we could use a
    // ConstantFunction3D, but to save the overhead of
    // a lookup we have a special mode for constant values
    bool hasConstantA;
    bool hasConstantB;
    Real constantA;
    Real constantB;

public:
    ShapeModulus(Grid3D* distanceField, Real a = 300, Real b = 0):
        distanceField(distanceField), hasConstantA(true), hasConstantB(true), constantA(a), constantB(b) {}

    ShapeModulus(Grid3D* distanceField, FieldFunction3D* a, FieldFunction3D* b):
        distanceField(distanceField), a(a), b(b), hasConstantA(false), hasConstantB(false) {}

    Real getFieldValue(const VEC3F& pos) const override {
        Real distance = (*distanceField)(pos);
        Real aValue   = (hasConstantA ? constantA : a->getFieldValue(pos));
        Real bValue   = (hasConstantB ? constantB : b->getFieldValue(pos)) ;
        Real radius   = exp( aValue * (distance - bValue ));

        return radius;
    }

};

class NoiseVersor: public R3Map {
public:
    // siv::PerlinNoise nx{ 1234567u };
    // siv::PerlinNoise ny{ 6543210u };
    // siv::PerlinNoise nz{ 1232101u };

    siv::PerlinNoise nx{ 000u };
    siv::PerlinNoise ny{ 000u };
    siv::PerlinNoise nz{ 000u };

    uint octaves;
    Real scale;

    NoiseVersor(uint octaves, Real scale): octaves(octaves), scale(scale) {
        // nx.reseed(83888u);
        // ny.reseed(39388u);
        // nz.reseed(17474u); // Decent

        // nx.reseed(888u);
        // ny.reseed(388u);
        // nz.reseed(174u); // No nose

        nx.reseed(88u);
        ny.reseed(38u);
        nz.reseed(14u); // Better
    }

    virtual VEC3F getFieldValue(const VEC3F& pos) const override {
        VEC3F p = pos * scale;

        VEC3F v(
            nx.octave3D_01(p.x(), p.y(), p.z(), octaves) * 2 - 1,
            ny.octave3D_01(p.x(), p.y(), p.z(), octaves) * 2 - 1,
            nz.octave3D_01(p.x(), p.y(), p.z(), octaves) * 2 - 1
            );

        return v.normalized();
    }
};

class PortalMap: public R3Map {
public:
    R3Map *map;

    vector<VEC3F>           portalCenters;
    vector<AngleAxis<Real>> portalRotations;

    Real  portalRadius;
    Real  portalScale;
    FieldFunction3D *mask;

    PortalMap(R3Map *map, vector<VEC3F> portalCenters, vector<AngleAxis<Real>> portalRotations, Real portalRadius, Real portalScale, FieldFunction3D *mask = 0): map(map), portalCenters(portalCenters), portalRotations(portalRotations), portalRadius(portalRadius), portalScale(portalScale), mask(mask) {}

    virtual VEC3F getFieldValue(const VEC3F& pos) const override {

        VEC3F closestPortal = portalCenters[0];
        AngleAxis<Real> portalRot = portalRotations[0];

        int i = 0;
        for (auto p : portalCenters) {
            if ((pos - closestPortal).norm() > (pos - p).norm()) {
                closestPortal = p;
                portalRot = portalRotations[i];
            }
            i++;
        }

        Real  dist = (pos - closestPortal).norm();
        VEC3F ang  = (pos - closestPortal).normalized();

        if (dist < portalRadius) {
            if (mask && (*mask)(pos) <= 0) {
                return (*map)(pos);
            }
            VEC3F out = (dist * ang * portalScale);
            out = portalRot * out;
            return out;
        } else {
            return (*map)(pos);
        }

    }
};

// =============== INSPECTION FIELDS =======================

class QuatQuatRotField: public FieldFunction3D {
public:
    QuatMap *func;
    int i;
    QuatQuatRotField(QuatMap* func, int i): func(func), i(i) {}

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
    QuatMap *func;
    QuatQuatMagField(QuatMap* func): func(func) {}

    virtual Real getFieldValue(const VEC3F& pos) const override {
        QUATERNION input(pos[0], pos[1], pos[2], 0);
        QUATERNION transformed = func->getFieldValue(input);
        return transformed.magnitude();
    }

};

#endif
