#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "SETTINGS.h"

using namespace std;
class Triangle {
public:
    Triangle(VEC3F* v0, VEC3F* v1, VEC3F* v2);
    Triangle(Triangle* triangle);
    Triangle();
    void draw();
    void recomputeNormal();
    VEC3F& color(){return _color;}
    void updateColor(const VEC3F &RHS) {_color = RHS;}
    // for equality, sort the indices first and then determine
    // if they are all the same
    bool operator==(const Triangle &RHS) const;
    Triangle& operator=(const Triangle &RHS);

    VEC3F& normal() { return _normal; };
    VEC3F normal() const { return _normal; };
    Real area() const { return _area; };
    Real maxEdgeLength();
    VEC3F*& vertex(int x) { return _vertices[x]; };
    const VEC3F* vertex(int x) const { return _vertices[x]; };
    void setVertex(int x, VEC3F* vertex) { _vertices[x] = vertex; };
    VEC3F centroid();
    vector<VEC3F*>& vertices() { return _vertices; };
    const vector<VEC3F*>& vertices() const { return _vertices; };

    bool positionsEqual(Triangle& RHS);

    bool intersects(Triangle& RHS);

    // intersect with line segment
    bool intersects(VEC3F& start, VEC3F& end);

    // get the bounding box
    void boundingBox(VEC3F& mins, VEC3F& maxs);

private:
    VEC3F _color;

    vector<VEC3F*> _vertices;

    VEC3F _normal;
    Real _area;
};

#endif
