/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _QUATERNION_H
#define _QUATERNION_H

#include "../../src/SETTINGS.h"
#include <iostream>

using namespace std;

class QUATERNION {

public:
  QUATERNION();
  //QUATERNION(const Real& real, const VEC3F& imaginary);
  QUATERNION(Real w, Real x, Real y, Real z);
  QUATERNION(Real w, Real x);
  //QUATERNION(const __m128 v);
  QUATERNION(const VEC3F& vector);
  QUATERNION(const QUATERNION& q);

  // overload operators
  QUATERNION& operator=(const VEC3F& vec);
  //inline QUATERNION& operator=(const QUATERNION q) { _v = q._v; return *this; };
  //inline QUATERNION& operator*=(const Real r) { _v = _mm_mul_ps(_v, _mm_set_ps1(r)); return *this; };
  //inline QUATERNION& operator-=(const QUATERNION& q) { _v = _mm_sub_ps(_v, q._v); return *this; };
  //inline QUATERNION& operator+=(const QUATERNION& q) { _v = _mm_add_ps(_v, q._v); return *this; };
  inline QUATERNION& operator=(const QUATERNION q)   { _w = q.w(); _x = q.x(); _y = q.y(); _z = q.z(); return *this;};
  inline QUATERNION& operator*=(const Real r)        { _w *= r; _x *= r; _y *= r; _z *= r;  return *this; };
  inline QUATERNION& operator*=(const QUATERNION& q) {
    const Real& x = this->_y * q._z - this->_z * q._y + q._w * this->_x + this->_w * q._x;
    const Real& y = this->_z * q._x - this->_x * q._z + q._w * this->_y + this->_w * q._y;
    const Real& z = this->_x * q._y - this->_y * q._x + q._w * this->_z + this->_w * q._z;
    const Real& w = this->_w * q._w - this->_x * q._x - q._y * this->_y - this->_z * q._z;
    _x = x; _y = y; _z = z; _w = w;
    return *this;
  };
  inline QUATERNION& operator-=(const QUATERNION& q) { _w -= q.w(); _x -= q.x(); _y -= q.y(); _z -= q.z(); return *this; };
  inline QUATERNION& operator+=(const QUATERNION& q) { _w += q.w(); _x += q.x(); _y += q.y(); _z += q.z(); return *this; };
  inline Real& operator[](const int x) { return _entries[x]; };
  inline Real operator[](const int x) const { return _entries[x]; };

  // compute the conjugate
  QUATERNION conjugate() const;

  // compute the inverse
  QUATERNION inverse() const;

  // populate a VECTOR of size 4
  VECTOR toVector();

  // normalize
  void normalize();

  // negative the imaginary component, effectively retrieving the transpose
  // of the desired rotation
  void negateIm();

  // accessors
  const Real& x() const { return _x; };
  const Real& y() const { return _y; };
  const Real& z() const { return _z; };
  const Real& w() const { return _w; };

  Real& x() { return _x; };
  Real& y() { return _y; };
  Real& z() { return _z; };
  Real& w() { return _w; };
  //const __m128& v() const { return _v; };

  // 2 norm of components
  inline Real magnitude() const { return sqrt(_w * _w + _x * _x + _y * _y + _z * _z); };

  // in-place equality
  void equals(const QUATERNION& q);

  // check if any components are nan
  inline bool anyNans() { return _x != _x || _y != _y || _z != _z || _w != _w; };

  // do a julia set iteration
  void juliaIteration(const QUATERNION& c);

  void write(FILE* file) const;
  void read(FILE* file);

  inline void multiplyAdd(const QUATERNION& point, const QUATERNION& add)
  {
    const Real x = _x;
    const Real y = _y;
    const Real z = _z;
    const Real w = _w;
    _x = y * point._z - z * point._y + point._w * x + w * point._x + add._x;
    _y = z * point._x - x * point._z + point._w * y + w * point._y + add._y;
    _z = x * point._y - y * point._x + point._w * z + w * point._z + add._z;
    _w = w * point._w - x * point._x - point._y * y - z * point._z + add._w;
  };

  static bool wCompare(const QUATERNION& i, const QUATERNION& j) { return i._w < j._w; };
  static bool xCompare(const QUATERNION& i, const QUATERNION& j) { return i._x < j._x; };
  static bool yCompare(const QUATERNION& i, const QUATERNION& j) { return i._y < j._y; };
  static bool zCompare(const QUATERNION& i, const QUATERNION& j) { return i._z < j._z; };

  // take the exponential
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION exp() const;

  // take the log
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION log() const;

  // take the power
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION pow(const Real& exponent) const;

  // dot against another quaternion
  inline Real dot(const QUATERNION& rhs) const { return _w * rhs._w + _x * rhs._x + _y * rhs._y + _z * rhs._z; };

private:
  // _w and _entries[0] are the real component
  union {
     struct { Real _w,_x,_y,_z; };
     Real _entries[4];
     //__m128 _v;
  };
};

QUATERNION operator/(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator/(const QUATERNION& left, const Real& right);
Real operator^(const QUATERNION& left, const QUATERNION& right);
ostream &operator<<(ostream &out, const QUATERNION& q);

/*
// SSE for this one has been lifted from Eigen!
#ifndef SWIZZLE
#define SWIZZLE(v,p,q,r,s) \
  (_mm_castsi128_ps(_mm_shuffle_epi32( _mm_castps_si128(v), ((s)<<6|(r)<<4|(q)<<2|(p)))))
#endif
inline QUATERNION operator*(const QUATERNION& left, const QUATERNION& right)
{
  const __m128 mask = _mm_castsi128_ps(_mm_setr_epi32(0,0,0,0x80000000));
  const __m128& a = left.v();
  const __m128& b = right.v();
  const __m128 flip1 = _mm_xor_ps(_mm_mul_ps(SWIZZLE(a,1,2,0,2),
                                             SWIZZLE(b,2,0,1,2)), mask);
  const __m128 flip2 = _mm_xor_ps(_mm_mul_ps(SWIZZLE(a,3,3,3,1),
                                             SWIZZLE(b,0,1,2,1)), mask);
  return QUATERNION(_mm_add_ps(_mm_sub_ps(_mm_mul_ps(a,SWIZZLE(b,3,3,3,3)),
                                          _mm_mul_ps(SWIZZLE(a,2,0,1,0), SWIZZLE(b,1,2,0,0))),
                                          _mm_add_ps(flip1,flip2)));
};
*/
inline QUATERNION operator*(const QUATERNION& left, const Real& right) {
  //return QUATERNION(_mm_mul_ps(left.v(), _mm_set_ps1(right)));
  return QUATERNION(left.w() * right, left.x() * right,
                    left.y() * right, left.z() * right);
};
inline QUATERNION operator*(const Real& left, const QUATERNION& right) {
  //return QUATERNION(_mm_mul_ps(right.v(), _mm_set_ps1(left)));
  return QUATERNION(right.w() * left, right.x() * left,
                    right.y() * left, right.z() * left);
};
inline QUATERNION operator-(const QUATERNION& left, const QUATERNION& right) {
  //return QUATERNION(_mm_sub_ps(left.v(), right.v()));
  return QUATERNION(left.w() - right.w(), left.x() - right.x(),
                    left.y() - right.y(), left.z() - right.z());
};
inline QUATERNION operator+(const QUATERNION& left, const QUATERNION& right) {
  //return QUATERNION(_mm_add_ps(left.v(), right.v()));
  return QUATERNION(left.w() + right.w(), left.x() + right.x(),
                    left.y() + right.y(), left.z() + right.z());
};
QUATERNION operator*(const QUATERNION& left, const QUATERNION& right);
#endif
