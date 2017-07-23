//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_TINYVECTOR_H
#define OHMMS_TINYVECTOR_H

/***************************************************************************
 *
 * \class TinyVector
 * \brief Pooma/AppyTypes/Vecktor is modified to work with PETE.
 *
 * The POOMA Framework
 *
 * This program was prepared by the Regents of the University of
 * California at Los Alamos National Laboratory (the University) under
 * Contract No.  W-7405-ENG-36 with the U.S. Department of Energy (DOE).
 * The University has certain rights in the program pursuant to the
 * contract and the program should not be copied or distributed outside
 * your organization.  All rights in the program are reserved by the DOE
 * and the University.  Neither the U.S.  Government nor the University
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMA for more details
 *
 ***************************************************************************/


#include <iomanip>
#include <algorithm>

namespace qmcplusplus
{
/** Fixed-size array. candidate for array<T,D>
 */
template<class T, unsigned D>
struct TinyVector
{
  typedef T Type_t;
  typedef TinyVector<T,D> This_t;
  enum { Size = D };
  T X[Size];

  // Default Constructor initializes to zero.
  inline TinyVector()
  {
    *this=T(0);
  }

  // Templated TinyVector constructor.
  template<class T1, unsigned D1>
  inline TinyVector(const TinyVector<T1,D1> &rhs)
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    std::copy_n(rhs.begin(),D,X);
  }

  // Constructor from a single T
  inline TinyVector(const T& x00)
  {
    *this=x00;
  }

  // Constructors for fixed dimension
  inline TinyVector(const T& x00, const T& x01)
  {
    X[0] = x00;
    X[1] = x01;
  }
  inline TinyVector(const T& x00, const T& x01, const T& x02)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
  }
  inline TinyVector(const T& x00, const T& x01, const T& x02, const T& x03)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
    X[3] = x03;
  }
  inline TinyVector(const T& x00, const T& x01, const T& x02, const T& x03,
             const T& x10, const T& x11, const T& x12, const T& x13,
             const T& x20, const T& x21, const T& x22, const T& x23,
             const T& x30, const T& x31, const T& x32, const T& x33)
  {
    X[0]  = x00;
    X[1]  = x01;
    X[2]  = x02;
    X[3]  = x03;
    X[4]  = x10;
    X[5]  = x11;
    X[6]  = x12;
    X[7]  = x13;
    X[8]  = x20;
    X[9]  = x21;
    X[10] = x22;
    X[11] = x23;
    X[12] = x30;
    X[13] = x31;
    X[14] = x32;
    X[15] = x33;
  }

  inline TinyVector(const T* restrict base, int offset)
  {
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      X[i]=base[i*offset];
  }

  // Destructor
  ~TinyVector() { }

  inline int size() const
  {
    return D;
  }

  inline int byteSize() const
  {
    return D*sizeof(T);
  }

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i)
  {
    return X[i];
  }
  inline const Type_t& operator[](unsigned int i) const
  {
    return X[i];
  }

  inline Type_t* data()
  {
    return  X;
  }
  inline const Type_t* data() const
  {
    return  X;
  }
  inline Type_t* begin()
  {
    return  X;
  }
  inline const Type_t* begin() const
  {
    return  X;
  }
  inline Type_t* end()
  {
    return  X+D;
  }
  inline const Type_t* end() const
  {
    return  X+D;
  }

  // operator = TinyVector
  template<class T1, unsigned D1>
  inline This_t& operator=(const TinyVector<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    std::copy_n(rhs.begin(),D,X);
    return *this;
  }

  // operator = scalar
  inline This_t& operator=(const T& rhs)
  {
    std::fill_n(X,D,rhs);
    return *this;
  }

  // operator += TinyVector
  template<class T1, unsigned D1>
  inline This_t& operator+=(const TinyVector<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      X[i]+=rhs[i];
    return *this;
  }

  // operator -= TinyVector
  template<class T1, unsigned D1>
  inline This_t& operator-=(const TinyVector<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      X[i]-=rhs[i];
    return *this;
  }

  // operator *= scalar
  template<class T1>
  inline This_t& operator*=(const T1& rhs)
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      X[i]*=rhs;
    return *this;
  }

  // operator /= scalar
  template<class T1>
  inline This_t& operator/=(const T1& rhs)
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      X[i]/=rhs;
    return *this;
  }

  // operator + TinyVector
  template<class T1, unsigned D1>
  inline This_t operator+(const TinyVector<T1,D1>& rhs) const
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    This_t tmp;
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      tmp[i]=X[i]+rhs[i];
    return tmp;
  }

  // operator - TinyVector
  template<class T1, unsigned D1>
  inline This_t operator-(const TinyVector<T1,D1>& rhs) const
  {
    static_assert(D1==D, "Dimession mismatching in TinyVector");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    This_t tmp;
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      tmp[i]=X[i]-rhs[i];
    return tmp;
  }

  // operator * scalar
  template<class T1>
  inline This_t operator*(const T1& rhs) const
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    This_t tmp;
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      tmp[i]=X[i]*rhs;
    return tmp;
  }

  // operator / scalar
  template<class T1>
  inline This_t operator/(const T1& rhs) const
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in TinyVector!");
    This_t tmp;
    #pragma unroll(D)
    for(int i=0; i<D; ++i)
      tmp[i]=X[i]/rhs;
    return tmp;
  }

};

// scalar * TinyVector
template<class T, class T1, unsigned D>
inline TinyVector<T,D> operator*(const T1& lhs, const TinyVector<T,D>& rhs)
{
  static_assert(std::is_convertible<T,T1>::value, "Inconvertible types in TinyVector!");
  TinyVector<T,D> tmp;
  #pragma unroll(D)
  for(int i=0; i<D; ++i)
    tmp[i]=lhs*rhs[i];
  return tmp;
}

//----------------------------------------------------------------------
// dot product
//----------------------------------------------------------------------
template <class T, class T1, unsigned D>
inline T dot(const TinyVector<T,D>& lhs, const TinyVector<T1,D>& rhs)
{
  static_assert(std::is_convertible<T,T1>::value, "Inconvertible types in TinyVector!");
  T tmp=0;
  #pragma unroll(D)
  for(int i=0; i<D; ++i)
    tmp+=lhs[i]*rhs[i];
  return tmp;
}

//----------------------------------------------------------------------
// cross product
//----------------------------------------------------------------------

template<class T, class T1>
inline TinyVector<T,3> cross(const TinyVector<T,3>& a, const TinyVector<T1,3>& b)
{
  static_assert(std::is_convertible<T,T1>::value, "Inconvertible types in TinyVector!");
  TinyVector<T,3> cross;
  cross[0] = a[1]*b[2] - a[2]*b[1];
  cross[1] = a[2]*b[0] - a[0]*b[2];
  cross[2] = a[0]*b[1] - a[1]*b[0];
  return cross;
}

//----------------------------------------------------------------------
// outerProduct product
//----------------------------------------------------------------------

/* TODO
template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
outerProduct(const TinyVector<T1,D> &lhs, const TinyVector<T2,D> &rhs)
{
  return OuterProduct< TinyVector<T1,D> , TinyVector<T2,D> > :: apply(lhs,rhs);
}
*/

//----------------------------------------------------------------------
// I/O
template<class T>
struct printTinyVector {};

// specialized for Vector<TinyVector<T,D> >
template<class T, unsigned D>
struct printTinyVector< TinyVector<T,D>  >
{
  inline static void
  print(std::ostream& os, const TinyVector<T,D>& r)
  {
    for(int d=0; d<D; d++)
      os << std::setw(18) << std::setprecision(10) << r[d];
  }
};

// specialized for Vector<TinyVector<T,2> >
template<class T>
struct printTinyVector<TinyVector<T,2> >
{
  inline static void
  print(std::ostream& os, const TinyVector<T,2>& r)
  {
    os << std::setw(18) << std::setprecision(10) << r[0]
       << std::setw(18) << std::setprecision(10) << r[1];
  }
};

// specialized for Vector<TinyVector<T,3> >
template<class T>
struct printTinyVector<TinyVector<T,3> >
{
  inline static void
  print(std::ostream& os, const TinyVector<T,3>& r)
  {
    os << std::setw(18) << std::setprecision(10) << r[0]
       << std::setw(18) << std::setprecision(10) << r[1]
       << std::setw(18) << std::setprecision(10) << r[2];
  }
};


template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const TinyVector<T,D>& rhs)
{
  printTinyVector<TinyVector<T,D> >::print(out,rhs);
  return out;
}

template<class T, unsigned D>
std::istream& operator>>(std::istream& is, TinyVector<T,D>& rhs)
{
  //printTinyVector<TinyVector<T,D> >::print(out,rhs);
  for(int i=0; i<D; i++)
    is >> rhs[i];
  return is;
}
}

#endif // VEKTOR_H

