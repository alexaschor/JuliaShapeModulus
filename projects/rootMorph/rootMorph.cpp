#include "Quaternion/POLYNOMIAL_4D.h"
#include "Quaternion/QUATERNION.h"
#include "SETTINGS.h"
#include "field.h"
#include "mesh.h"
#include "triangle.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// This is a modified version of the Wild Magic DistVector3Triangle3 class.
// This is what is used to generate the distance grid.
//
// The license info is as below:
//
// Wild Magic Source Code
// David Eberle
// http://www.geometrictools.com
// Copyright (c) 1998-2008
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
//
// Version: 4.0.1 (2007/05/06)
//----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////

Real norm2(VEC3F v) {
    return v.squaredNorm();
}

Real norm(VEC3F v) {
    return v.norm();
}

VEC3F cross(VEC3F a, VEC3F b) {
    return a.cross(b);
}

Real pointFaceDistanceSq(Triangle& f, const VEC3F& point) {
    VEC3F v0 = *(f.vertex(0));
    VEC3F v1 = *(f.vertex(1));
    VEC3F v2 = *(f.vertex(2));
    VEC3F kDiff = v0 - point;
    VEC3F kEdge0 = v1 - v0;
    VEC3F kEdge1 = v2 - v0;
    Real fA00 = norm2(kEdge0);
    Real fA01 = kEdge0.dot(kEdge1);
    Real fA11 = norm2(kEdge1);
    Real fB0 = kDiff.dot(kEdge0);
    Real fB1 = kDiff.dot(kEdge1);
    Real fC = norm2(kDiff);
    Real fDet = fabs(fA00*fA11-fA01*fA01);
    Real fS = fA01*fB1-fA11*fB0;
    Real fT = fA01*fB0-fA00*fB1;
    Real fSqrDistance;

    if (fS + fT <= fDet) {
        if (fS < (Real)0.0) {
            if (fT < (Real)0.0)  // region 4
            {
                if (fB0 < (Real)0.0) {
                    fT = (Real)0.0;
                    if (-fB0 >= fA00) {
                        fS = (Real)1.0;
                        fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                    }
                    else {
                        fS = -fB0/fA00;
                        fSqrDistance = fB0*fS+fC;
                    }
                }
                else {
                    fS = (Real)0.0;
                    if (fB1 >= (Real)0.0) {
                        fT = (Real)0.0;
                        fSqrDistance = fC;
                    }
                    else if (-fB1 >= fA11) {
                        fT = (Real)1.0;
                        fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                    }
                    else {
                        fT = -fB1/fA11;
                        fSqrDistance = fB1*fT+fC;
                    }
                }
            }
            else  // region 3
            {
                fS = (Real)0.0;
                if (fB1 >= (Real)0.0) {
                    fT = (Real)0.0;
                    fSqrDistance = fC;
                }
                else if (-fB1 >= fA11) {
                    fT = (Real)1.0;
                    fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                }
                else {
                    fT = -fB1/fA11;
                    fSqrDistance = fB1*fT+fC;
                }
            }
        }
        else if (fT < (Real)0.0)  // region 5
        {
            fT = (Real)0.0;
            if (fB0 >= (Real)0.0) {
                fS = (Real)0.0;
                fSqrDistance = fC;
            }
            else if (-fB0 >= fA00) {
                fS = (Real)1.0;
                fSqrDistance = fA00+((Real)2.0)*fB0+fC;
            }
            else {
                fS = -fB0/fA00;
                fSqrDistance = fB0*fS+fC;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            Real fInvDet = ((Real)1.0)/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
        }
    }
    else {
        Real fTmp0, fTmp1, fNumer, fDenom;

        if (fS < (Real)0.0)  // region 2
        {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if (fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if (fNumer >= fDenom) {
                    fS = (Real)1.0;
                    fT = (Real)0.0;
                    fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                }
                else {
                    fS = fNumer/fDenom;
                    fT = (Real)1.0 - fS;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
                }
            }
            else {
                fS = (Real)0.0;
                if (fTmp1 <= (Real)0.0) {
                    fT = (Real)1.0;
                    fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                }
                else if (fB1 >= (Real)0.0) {
                    fT = (Real)0.0;
                    fSqrDistance = fC;
                }
                else {
                    fT = -fB1/fA11;
                    fSqrDistance = fB1*fT+fC;
                }
            }
        }
        else if (fT < (Real)0.0)  // region 6
        {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if (fTmp1 > fTmp0) {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-((Real)2.0)*fA01+fA11;
                if (fNumer >= fDenom) {
                    fT = (Real)1.0;
                    fS = (Real)0.0;
                    fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                }
                else {
                    fT = fNumer/fDenom;
                    fS = (Real)1.0 - fT;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
                }
            }
            else {
                fT = (Real)0.0;
                if (fTmp1 <= (Real)0.0) {
                    fS = (Real)1.0;
                    fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                }
                else if (fB0 >= (Real)0.0) {
                    fS = (Real)0.0;
                    fSqrDistance = fC;
                }
                else {
                    fS = -fB0/fA00;
                    fSqrDistance = fB0*fS+fC;
                }
            }
        }
        else  // region 1
        {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if (fNumer <= (Real)0.0) {
                fS = (Real)0.0;
                fT = (Real)1.0;
                fSqrDistance = fA11+((Real)2.0)*fB1+fC;
            }
            else {
                fDenom = fA00-2.0f*fA01+fA11;
                if (fNumer >= fDenom) {
                    fS = (Real)1.0;
                    fT = (Real)0.0;
                    fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                }
                else {
                    fS = fNumer/fDenom;
                    fT = (Real)1.0 - fS;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
                }
            }
        }
    }

    // account for numerical round-off error
    if (fSqrDistance < (Real)0.0) {
        fSqrDistance = (Real)0.0;
    }
    return fSqrDistance;
}

VEC3F getBarycentricCoordinates(Triangle& t, VEC3F point) {
    // Project into plane of source triangle
    Real offset = (point.dot(t.normal())) - (t.vertex(0)->dot(t.normal()));
    VEC3F p = point - (offset * t.normal());

    // Make aliases for vertices
    VEC3F a = *t.vertex(0), b = *t.vertex(1), c = *t.vertex(2);

    // Compute barycentric coordinates
    Real area_a = norm(cross((b - p), (c - p))) / 2.0;
    Real area_b = norm(cross((p - a), (c - a))) / 2.0;
    Real area_c = norm(cross((p - a), (b - a))) / 2.0;

    Real area_all = norm(cross((c - a), (b - a))) / 2.0;

    Real alpha = area_a / area_all;
    Real beta = area_b / area_all;
    Real gamma = area_c / area_all;

    return VEC3F(alpha, beta, gamma);
}

int getClosestTriangleIndex(VEC3F point, Mesh& mesh) {
    // Find closest projection to a triangle in the mesh
    size_t closestIdx = 0;
    Real minDistance = std::numeric_limits<Real>::max();
    bool found = false;

    for (size_t i = 0; i < mesh.numFaces(); ++i) {
        Triangle t = mesh.triangle(i);
        Real distance = pointFaceDistanceSq(t, point);

        if (!found || distance < minDistance) {
            VEC3F coords = getBarycentricCoordinates(t, point);
            if (abs((coords[0] + coords[1] + coords[2]) - 1) < 0.1) {
                minDistance = distance;
                closestIdx = i;
                found = true;
            }
        }
    }

    if (!found) {
        cout << "Unable to match point " << point << " to position on mesh (tried " << mesh.numFaces() << " triangles)!" << endl;
        return -1;
    }

    return closestIdx;
}

VEC3F matchTriangleDeformation(VEC3F point, Triangle source, Triangle destination) {
    source.recomputeNormal();
    destination.recomputeNormal();

    // Project into plane of source triangle
    Real offset = (point.dot(source.normal())) - (source.vertex(0)->dot(source.normal()));
    VEC3F p = point - (offset * source.normal());

    // Make aliases for vertices
    VEC3F a = *source.vertex(0), b = *source.vertex(1), c = *source.vertex(2);

    // Compute barycentric coordinates
    Real area_a = norm(cross((b - p), (c - p))) / 2.0;
    Real area_b = norm(cross((p - a), (c - a))) / 2.0;
    Real area_c = norm(cross((p - a), (b - a))) / 2.0;

    Real area_all = norm(cross((c - a), (b - a))) / 2.0;

    Real alpha = area_a / area_all;
    Real beta = area_b / area_all;
    Real gamma = area_c / area_all;


    // Project into destination triangle
    VEC3F p_dest = (alpha * *destination.vertex(0)) + (beta * *destination.vertex(1)) + (gamma * *destination.vertex(2));

    // Offset along destination normal
    return p_dest;

}

VEC3F matchMeshDeformation(VEC3F point, Mesh& m1, Mesh& m2) {

    if (&m1 == &m2) return point;

    int idx = getClosestTriangleIndex(point, m1);

    if (idx == -1) return point; //XXX For now, just return the point as is if it doesn't project to anywhere

    Triangle from = m1.triangle(idx);
    Triangle to   = m2.triangle(idx);

    return matchTriangleDeformation(point, m1.triangle(idx), m2.triangle(idx));
}

void fitMeshToAABB(Mesh& mesh, AABB box) {
    AABB bounds(mesh.vertices[0], mesh.vertices[1]);

    for (VEC3F v : mesh.vertices) {
        bounds.max() = bounds.max().cwiseMax(v);
        bounds.min() = bounds.min().cwiseMin(v);
    }

    VEC3F offset = box.center() - bounds.center();
    VEC3F scale = box.span().cwiseQuotient(bounds.span());

    for (VEC3F &v : mesh.vertices) {
        v += offset;
        v = v.cwiseProduct(scale);
    }

    AABB checkBounds(mesh.vertices[0], mesh.vertices[1]);
    for (VEC3F v : mesh.vertices) {
        checkBounds.max() = checkBounds.max().cwiseMax(v);
        checkBounds.min() = checkBounds.min().cwiseMin(v);
    }

    PRINTV3(checkBounds.min());
    PRINTV3(checkBounds.max());
}

void fitMeshToAABBFromAABB(Mesh& mesh, AABB from, AABB to) {

    VEC3F offset = to.center() - from.center();
    VEC3F scale = to.span().cwiseQuotient(from.span());

    for (VEC3F &v : mesh.vertices) {
        v += offset;
        v = v.cwiseProduct(scale);
    }
}


int main(int argc, char* argv[]) {

    if(argc < 5) {
        cout << "USAGE: " << endl;
        cout << "To generate a sequence of polynomials from a *.poly file, matching the deformation of a mesh sequence" << endl;
        cout << " " << argv[0] << " <*.poly base> <*.obj original> <min x> <min y> <min z> <max x> <max y> <max z> <output path/prefix> <*.obj frame 1> <*.obj frame 2> ... <*.obj frame n>\n";
        exit(-1);
    }

    POLYNOMIAL_4D basePoly{string(argv[1])};

    AABB unitBounds(VEC3F(-0.5, -0.5, -0.5), VEC3F(0.5, 0.5, 0.5));

    AABB origBounds(VEC3F(atof(argv[3]), atof(argv[4]), atof(argv[5])), VEC3F(atof(argv[6]), atof(argv[7]), atof(argv[8])));

    Mesh origMesh;
    origMesh.readOBJ(string(argv[2]));
    fitMeshToAABBFromAABB(origMesh, origBounds, unitBounds);

    for (int i = 10; i < argc; ++i) {
        Mesh frameMesh;
        frameMesh.readOBJ(string(argv[i]));
        fitMeshToAABBFromAABB(frameMesh, origBounds, unitBounds);

        vector<QUATERNION> roots;
        for (QUATERNION oldRoot : basePoly.roots()) {
            VEC3F oldRoot3D(oldRoot[0], oldRoot[1], oldRoot[2]);
            VEC3F newRoot3D = matchMeshDeformation(oldRoot3D, origMesh, frameMesh);

            roots.push_back(QUATERNION(newRoot3D[0], newRoot3D[1], newRoot3D[2], oldRoot[3]));
        }

        POLYNOMIAL_4D framePoly(roots, basePoly.powers());

        char buf[256];
        sprintf(buf, "%s_%04d.poly", argv[9], i-9);

        FILE *out = fopen(buf, "wb");
        framePoly.write(out);

        printf("Wrote polynomial to %s\n", buf);
    }


    return 0;
}
