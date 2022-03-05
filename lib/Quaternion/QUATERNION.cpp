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
#include "QUATERNION.h"

QUATERNION::QUATERNION() :
  _w(0), _x(0), _y(0), _z(0)
{
}

QUATERNION::QUATERNION(Real w, Real x, Real y, Real z) :
  _w(w), _x(x), _y(y), _z(z)
{
}

QUATERNION::QUATERNION(Real w, Real x) :
  _w(w), _x(x), _y(0), _z(0)
{
}

QUATERNION::QUATERNION(const VEC3F& vector) :
  _w(0), _x(vector[0]), _y(vector[1]), _z(vector[2])
{
}

QUATERNION::QUATERNION(const QUATERNION& q) :
  _w(q._w), _x(q._x), _y(q._y), _z(q._z)
{
}

QUATERNION QUATERNION::conjugate() const
{
  QUATERNION final(_w, -_x, -_y, -_z);
  return final;
}

QUATERNION& QUATERNION::operator=(const VEC3F& vec)
{
  this->_x = vec[0];
  this->_y = vec[1];
  this->_z = vec[2];
  this->_w = 0.0;
  return *this;
}

/*
// SSE for this one has been lifted from Eigen!
#define SWIZZLE(v,p,q,r,s) \
  (_mm_castsi128_ps(_mm_shuffle_epi32( _mm_castps_si128(v), ((s)<<6|(r)<<4|(q)<<2|(p)))))
QUATERNION operator*(const QUATERNION& left, const QUATERNION& right)
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
}
*/

QUATERNION operator*(const QUATERNION& left, const QUATERNION& right)
{
  QUATERNION final;
  final.x() = left.y() * right.z() - left.z() * right.y() + right.w() * left.x() + left.w() * right.x();
  final.y() = left.z() * right.x() - left.x() * right.z() + right.w() * left.y() + left.w() * right.y();
  final.z() = left.x() * right.y() - left.y() * right.x() + right.w() * left.z() + left.w() * right.z();
  final.w() = left.w() * right.w() - left.x() * right.x() - right.y() * left.y() - left.z() * right.z();
  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the inverse
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::inverse() const
{
  const Real magnitudeSq = (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] +
                           (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];
  return QUATERNION(conjugate() / magnitudeSq);
}

QUATERNION operator/(const QUATERNION& left, const QUATERNION& right)
{
  // build the reciprocal of the right hand side
  Real magnitudeSq = right[0] * right[0] + right[1] * right[1] +
                     right[2] * right[2] + right[3] * right[3];
  QUATERNION rightInv = right.conjugate() / magnitudeSq;

  return left * rightInv;
}

Real operator^(const QUATERNION& left, const QUATERNION& right)
{
  return left.x() * right.x() + left.y() * right.y() +
         left.z() * right.z() + left.w() * right.w();
}

ostream& operator<<(ostream &out, const QUATERNION& q)
{
  out << "(" << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << ")";
  return out;
}

QUATERNION operator/(const QUATERNION& left, const Real& right)
{
  const Real rightInv = 1.0 / right;
  return left * rightInv;
}

void QUATERNION::normalize()
{
  Real invMagnitude = 1.0 / sqrt(_x * _x + _y * _y + _z * _z + _w * _w);
  _x *= invMagnitude;
  _y *= invMagnitude;
  _z *= invMagnitude;
  _w *= invMagnitude;
}

//////////////////////////////////////////////////////////////////////
// negate the imaginary component -- effetively transpose the rotation
//////////////////////////////////////////////////////////////////////
void QUATERNION::negateIm()
{
  _x *= -1.0;
  _y *= -1.0;
  _z *= -1.0;
}

//////////////////////////////////////////////////////////////////////
// populate a VECTOR of size 4
//////////////////////////////////////////////////////////////////////
VECTOR QUATERNION::toVector()
{
  VECTOR final(4);
  final[0] = _entries[0];
  final[1] = _entries[1];
  final[2] = _entries[2];
  final[3] = _entries[3];
  return final;
}

#if 0
//////////////////////////////////////////////////////////////////////
// 2 norm of the entries
//////////////////////////////////////////////////////////////////////
Real QUATERNION::magnitude() const
{
  /*
  Real final = _entries[0] * _entries[0];
  final += _entries[1] * _entries[1];
  final += _entries[2] * _entries[2];
  final += _entries[3] * _entries[3];

  return sqrt(final);
  */
  return sqrt(_w * _w + _x * _x + _y * _y + _z * _z);
}
#endif

//////////////////////////////////////////////////////////////////////
// in-place equality
//////////////////////////////////////////////////////////////////////
void QUATERNION::equals(const QUATERNION& q)
{
  _entries[0] = q._entries[0];
  _entries[1] = q._entries[1];
  _entries[2] = q._entries[2];
  _entries[3] = q._entries[3];
}

//////////////////////////////////////////////////////////////////////
// do a julia set iteration
//////////////////////////////////////////////////////////////////////
void QUATERNION::juliaIteration(const QUATERNION& c)
{
  const Real copy[] = {_x, _y, _z, _w};

  //const Real copy01 = copy[0] * copy[1];
  //const Real copy12 = copy[1] * copy[2];
  //const Real copy02 = copy[0] * copy[2];
  const Real copy03 = 2.0f * copy[0] * copy[3];
  const Real copy31 = 2.0f * copy[3] * copy[1];
  const Real copy32 = 2.0f * copy[3] * copy[2];

  //_x = copy12 - copy12 + copy03 + copy03 + c[0];
  //_y = copy02 - copy02 + copy31 + copy31 + c[1];
  //_z = copy01 - copy01 + copy32 + copy32 + c[2];
  _x = copy03 + c[0];
  _y = copy31 + c[1];
  _z = copy32 + c[2];
  _w = copy[3] * copy[3] - copy[0] * copy[0] - copy[1] * copy[1] - copy[2] * copy[2] + c[3];
}

void QUATERNION::write(FILE* file) const
{
  if (sizeof(Real) == sizeof(double))
  {
    fwrite((void*)&_entries[0], sizeof(Real), 1, file);
    fwrite((void*)&_entries[1], sizeof(Real), 1, file);
    fwrite((void*)&_entries[2], sizeof(Real), 1, file);
    fwrite((void*)&_entries[3], sizeof(Real), 1, file);
  }
  else
  {
    double entries[] = {_entries[0], _entries[1], _entries[2], _entries[3]};
    fwrite((void*)&entries[0], sizeof(double), 1, file);
    fwrite((void*)&entries[1], sizeof(double), 1, file);
    fwrite((void*)&entries[2], sizeof(double), 1, file);
    fwrite((void*)&entries[3], sizeof(double), 1, file);
  }
}

void QUATERNION::read(FILE* file)
{
  if (sizeof(Real) == sizeof(double))
  {
    fread((void*)&_entries[0], sizeof(Real), 1, file);
    fread((void*)&_entries[1], sizeof(Real), 1, file);
    fread((void*)&_entries[2], sizeof(Real), 1, file);
    fread((void*)&_entries[3], sizeof(Real), 1, file);
  }
  else
  {
    double entries[4];
    fread((void*)&entries[0], sizeof(double), 1, file);
    fread((void*)&entries[1], sizeof(double), 1, file);
    fread((void*)&entries[2], sizeof(double), 1, file);
    fread((void*)&entries[3], sizeof(double), 1, file);

    _entries[0] = entries[0];
    _entries[1] = entries[1];
    _entries[2] = entries[2];
    _entries[3] = entries[3];
  }
}

//////////////////////////////////////////////////////////////////////
// take the exponential
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::exp() const
{
  const Real& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Real magnitude = v.norm();
  const Real exps = std::exp(s);
  const VEC3F vFinal = (v / magnitude) * sin(magnitude);

  return exps * QUATERNION(cos(magnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the log
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::log() const
{
  const Real& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Real qMagnitude = magnitude();
  const Real vMagnitude = v.norm();

  const VEC3F vNormalized = (vMagnitude > 0) ? v / vMagnitude : VEC3F(0,0,0);

  const VEC3F vFinal = vNormalized * acos(s / qMagnitude);

  return QUATERNION(std::log(qMagnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the power
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//
// some care has been taken to optimize this ...
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::pow(const Real& exponent) const
{
  const Real partial = _x * _x + _y * _y + _z * _z;
  const Real qMagnitude = sqrt(partial + _w * _w);
  const Real vMagnitude = sqrt(partial);
  const Real vMagnitudeInv = (vMagnitude > 0) ? 1.0 / vMagnitude : 0;

  const Real scale = exponent * acos(_w / qMagnitude) * vMagnitudeInv;

  const Real magnitude = scale * vMagnitude;
  const Real magnitudeInv = (magnitude > 0) ? 1.0 / magnitude : 0;
  const Real exps = std::exp(exponent * std::log(qMagnitude));

  const Real scale2 = scale * exps * magnitudeInv * sin(magnitude);
  return QUATERNION(exps * cos(magnitude),
                    scale2 * _x,
                    scale2 * _y,
                    scale2 * _z);
}
