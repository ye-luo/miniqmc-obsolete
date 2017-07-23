//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
// 		      Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef	OHMMS_TENSOR_H
#define	OHMMS_TENSOR_H

#include <algorithm>
#include <iostream>
#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"
/***************************************************************************
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



namespace qmcplusplus
{

/** Tensor<T,D>  class for D by D tensor
 *
 * @tparam T datatype
 * @tparm D dimension
 */
template<class T, unsigned D>
class Tensor
{
public:

  typedef T Type_t;
  typedef Tensor<T,D> This_t;
  enum { ElemDim = 2 };
  enum { Size = D*D };

  // Default Constructor
  inline Tensor()
  {
    *this=T(0);
  }

  // Templated Tensor constructor.
  template<class T1, unsigned D1>
  inline Tensor(const Tensor<T1,D1> &rhs)
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    std::copy_n(rhs.begin(),Size,X);
  }

  // constructor from a single T
  inline Tensor(const T& x00)
  {
    *this=x00;
  }

  // constructors for fixed dimension
  inline Tensor(const T& x00, const T& x10, const T& x01, const T& x11)
  {
    X[0] = x00;
    X[1] = x10;
    X[2] = x01;
    X[3] = x11;
  }
  inline Tensor(const T& x00, const T& x10, const T& x20, const T& x01, const T& x11,
         const T& x21, const T& x02, const T& x12, const T& x22)
  {
    X[0] = x00;
    X[1] = x10;
    X[2] = x20;
    X[3] = x01;
    X[4] = x11;
    X[5] = x21;
    X[6] = x02;
    X[7] = x12;
    X[8] = x22;
  }

  // destructor
  ~Tensor() { };

  // Methods

  inline void diagonal(const T& rhs)
  {
    for (int i = 0 ; i < D ; i++ )
      (*this)(i,i) = rhs;
  }

  inline void add2diagonal(T rhs)
  {
    for (int i = 0 ; i < D ; i++ )
      (*this)(i,i) += rhs;
  }

  ///return the size
  inline int len()  const
  {
    return Size;
  }
  ///return the size
  inline int size() const
  {
    return Size;
  }

  /** return the i-th value or assign
   * @param i index [0,D*D)
   */
  inline Type_t &operator[]( unsigned int i )
  {
    return X[i];
  }

  /** return the i-th value
   * @param i index [0,D*D)
   */
  inline const Type_t& operator[]( unsigned int i ) const
  {
    return X[i];
  }

  //TJW: add these 12/16/97 to help with NegReflectAndZeroFace BC:
  // These are the same as operator[] but with () instead:
  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }

  inline const Type_t& operator()(unsigned int i) const
  {
    return X[i];
  }
  //TJW.

  /** return the (i,j) component
   * @param i index [0,D)
   * @param j index [0,D)
   */
  inline Type_t operator()( unsigned int i,  unsigned int j ) const
  {
    return X[i*D+j];
  }

  /** return/assign the (i,j) component
   * @param i index [0,D)
   * @param j index [0,D)
   */
  inline Type_t& operator()( unsigned int i, unsigned int j )
  {
    return X[i*D+j];
  }

  inline TinyVector<T,D> getRow(unsigned int i)
  {
    TinyVector<T,D> res;
    for(int j=0; j<D; j++)
      res[j]=X[i*D+j];
    return res;
  }

  inline TinyVector<T,D> getColumn(unsigned int i)
  {
    TinyVector<T,D> res;
    for(int j=0; j<D; j++)
      res[j]=X[j*D+i];
    return res;
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
    return  X+Size;
  }
  inline const Type_t* end() const
  {
    return  X+Size;
  }

  // assignment operators
  // operator = Tensor
  template<class T1, unsigned D1>
  inline This_t& operator=(const Tensor<T,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    std::copy_n(rhs.begin(),Size,X);
    return *this;
  }

  // operator = scalar
  inline This_t& operator=(const T& rhs)
  {
    std::fill_n(X,Size,rhs);
    return *this;
  }

  // accumulation operators
  // operator += Tensor
  template<class T1, unsigned D1>
  inline This_t& operator+=(const Tensor<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    for(int i=0; i<Size; ++i)
      X[i]+=rhs[i];
    return *this;
  }

  // operator -= Tensor
  template<class T1, unsigned D1>
  inline This_t& operator-=(const Tensor<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    for(int i=0; i<Size; ++i)
      X[i]-=rhs[i];
    return *this;
  }

  // operator *= Tensor
  template<class T1, unsigned D1>
  inline This_t& operator*=(const Tensor<T1,D1>& rhs)
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    for(int i=0; i<D; ++i)
      for(int j=0; j<D; ++j)
      {
        Type_t sum = (*this)(i,0) * rhs(0,j);
        #pragma unroll
        for (int k=1; k<D; ++k)
          sum += (*this)(i,k) * rhs(k,j);
        (*this)(i,j) = sum;
      }
    return *this;
  }

  // operator *= scalar
  template<class T1>
  inline This_t& operator*=(const T1& rhs)
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    for(int i=0; i<Size; ++i)
      X[i]*=rhs;
    return *this;
  }

  // operator /= scalar
  template<class T1>
  inline This_t& operator/=(const T1& rhs)
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    for(int i=0; i<Size; ++i)
      X[i]/=rhs;
    return *this;
  }

  // operator + Tensor
  template<class T1, unsigned D1>
  inline This_t operator+(const Tensor<T1,D1>& rhs) const
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    This_t tmp;
    for(int i=0; i<Size; ++i)
      tmp[i]=X[i]+rhs[i];
    return tmp;
  }

  // operator - Tensor
  template<class T1, unsigned D1>
  inline This_t operator-(const Tensor<T1,D1>& rhs) const
  {
    static_assert(D1==D, "Dimession mismatching in Tensor");
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    This_t tmp;
    for(int i=0; i<Size; ++i)
      tmp[i]=X[i]-rhs[i];
    return tmp;
  }

  // operator * scalar
  template<class T1>
  inline This_t operator*(const T1& rhs) const
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    This_t tmp;
    for(int i=0; i<Size; ++i)
      tmp[i]=X[i]*rhs;
    return tmp;
  }

  // operator / scalar
  template<class T1>
  inline This_t operator/(const T1& rhs) const
  {
    static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
    This_t tmp;
    for(int i=0; i<Size; ++i)
      tmp[i]=X[i]/rhs;
    return tmp;
  }

private:

  // The elements themselves.
  T X[Size];

};

// Tensor * Tensor
template<class T, class T1, unsigned D>
inline Tensor<typename Promote<T,T1>::Type_t,D> operator*(const Tensor<T,D>& lhs, const Tensor<T1,D>& rhs)
{
  static_assert(std::is_convertible<T1,T>::value, "Inconvertible types in Tensor!");
  using Type_t=typename Promote<T,T1>::Type_t;
  Tensor<Type_t,D> tmp;
  for(int i=0; i<D; ++i)
    for(int j=0; j<D; ++j)
    {
      Type_t sum = lhs(i,0) * rhs(0,j);
      #pragma unroll
      for (int k=1; k<D; ++k) 
        sum += lhs(i,k) * rhs(k,j);
      tmp(i,j) = sum;
    }
  return tmp;
}

// scalar * Tensor
template<class T, class T1, unsigned D>
inline Tensor<T,D> operator*(const T1& lhs, const Tensor<T,D>& rhs)
{
  static_assert(std::is_convertible<T,T1>::value, "Inconvertible types in Tensor!");
  Tensor<T,D> tmp;
  for(int i=0; i<D*D; ++i)
    tmp[i]=lhs*rhs[i];
  return tmp;
}

//////////////////////////////////////////////////////////////////////
//
// Free functions
//
//////////////////////////////////////////////////////////////////////

/** trace \f$ result = \sum_k rhs(k,k)\f$
 * @param rhs a tensor
 */
template <class T, unsigned D>
inline T trace(const Tensor<T,D>& rhs)
{
  T result = 0.0;
  for (int i = 0 ; i < D ; i++ )
    result += rhs(i,i);
  return result;
}

/** transpose a tensor
 */
template <class T, unsigned D>
inline Tensor<T,D> transpose(const Tensor<T,D>& rhs)
{
  Tensor<T,D> result; // = Tensor<T,D>::DontInitialize();
  for (int j = 0 ; j < D ; j++ )
    for (int i = 0 ; i < D ; i++ )
      result(i,j) = rhs(j,i);
  return result;
}

/** Tr(a*b), \f$ \sum_i\sum_j a(i,j)*b(j,i) \f$
 */
template <class T1, class T2, unsigned D>
inline T1 trace(const Tensor<T1,D>& a, const Tensor<T2,D>& b)
{
  T1 result = 0.0;
  for (int i = 0 ; i < D ; i++ )
    for(int j=0; j<D; j++)
      result += a(i,j)*b(j,i);
  return result;
}

/** Tr(a^t *b), \f$ \sum_i\sum_j a(i,j)*b(i,j) \f$
 */
template <class T, class T1, unsigned D>
inline T traceAtB(const Tensor<T,D>& a, const Tensor<T1,D>& b)
{
  static_assert(std::is_convertible<T,T1>::value, "Inconvertible types in Tensor!");
  T result(0);
  for (int i = 0 ; i < D*D ; i++ )
    result += a(i)*b(i);
  return result;
}

/** Tensor-Tensor dot product \f$result(i,j)=\sum_k lhs(i,k)*rhs(k,j)\f$
 * @param lhs  a tensor
 * @param rhs  a tensor
 */
template < class T, class T1, unsigned D >
inline Tensor< typename Promote<T,T1>::Type_t, D >
dot(const Tensor<T,D> &lhs, const Tensor<T1,D> &rhs)
{
  return lhs*rhs;
}

/** Vector-Tensor dot product \f$result(i)=\sum_k lhs(k)*rhs(k,i)\f$
 * @param lhs  a vector
 * @param rhs  a tensor
 */
template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D>
dot(const TinyVector<T1,D> &lhs, const Tensor<T2,D> &rhs)
{
  return OTDot< TinyVector<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
}

/** Tensor-Vector dot product \f$result(i)=\sum_k lhs(i,k)*rhs(k)\f$
 * @param lhs  a tensor
 * @param rhs  a vector
 */
template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const Tensor<T1,D> &lhs, const TinyVector<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , TinyVector<T2,D> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// double dot product.
//----------------------------------------------------------------------

// template < class T1, class T2, unsigned D >
// inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
// dotdot(const Tensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
// {
//   return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
//}


//----------------------------------------------------------------------
// Outer product.
//----------------------------------------------------------------------

///** Vector-vector outter product \f$ result(i,j)=v1(i)*v2(j)\f$
// * @param v1 a vector
// * @param v2 a vector
// */
//  template<class T1, class T2, unsigned int D>
//  inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::type,D >
//  outerProduct(const TinyVector<T1,D>& v1, const TinyVector<T2,D>& v2)
//  {
//    typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
////    //#if (defined(POOMA_SGI_CC_721_TYPENAME_BUG) || defined(__MWERKS__))
////    //Tensor<T0,D> ret = Tensor<T0,D>::DontInitialize();
////    //#else
//    Tensor<T0,D> ret = typename Tensor<T0,D>::DontInitialize();
////    //#endif // POOMA_SGI_CC_721_TYPENAME_BUG
//    for (unsigned int i=0; i<D; ++i)
//      for (unsigned int j=0; j<D; ++j)
//        ret(i,j) = v1[i]*v2[j];
//    return ret;
//  }
//
//  template<class T1, unsigned int D>
//  inline Tensor<T1,D >
//  outerProduct(const TinyVector<T1,D>& v1, const TinyVector<T1,D>& v2)
//  {
//    Tensor<T1,D> ret = typename Tensor<T1,D>::DontInitialize();
////    //#endif // POOMA_SGI_CC_721_TYPENAME_BUG
//    for (unsigned int i=0; i<D; ++i)
//      for (unsigned int j=0; j<D; ++j)
//        ret(i,j) = v1[i]*v2[j];
//    return ret;
//  }
//


//----------------------------------------------------------------------
// I/O
template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const Tensor<T,D>& rhs)
{
  if (D >= 1)
  {
    for (int i=0; i<D; i++)
    {
      for (int j=0; j<D-1; j++)
      {
        out << rhs(i,j) << "  ";
      }
      out << rhs(i,D-1) << " ";
      if (i < D - 1)
        out << std::endl;
    }
  }
  else
  {
    out << " " << rhs(0,0) << " ";
  }
  return out;
}

template<class T, unsigned D>
std::istream& operator>>(std::istream& is, Tensor<T,D>& rhs)
{
  for(int i=0; i<D*D; i++)
    is >> rhs[i];
  return is;
}

}

#endif // OHMMS_TENSOR_H

