////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-

/**
 * @file FakeWaveFunction.h
 * @brief Top level wavefunction container
 *
 * Represents a product of wavefunction components (classes based on
 * WaveFunctionComponentBase).
 *
 * Corresponds to QMCWaveFunction/TrialWaveFunction.h in the QMCPACK source.
 */

#ifndef QMCPLUSPLUS_FAKEWAVEFUNCTIONS_H
#define QMCPLUSPLUS_FAKEWAVEFUNCTIONS_H
#include <Configuration.h>
#include <Particle/DistanceTable.h>
#include <QMCWaveFunctions/Jastrow/BsplineFunctor.h>
#include <QMCWaveFunctions/Jastrow/TwoBodyJastrowRef.h>
#include <QMCWaveFunctions/Jastrow/TwoBodyJastrow.h>

namespace qmcplusplus
{
/** A minimal TrialWavefunction
 */
struct FakeWaveFunctionBase
{
  using valT = OHMMS_PRECISION;
  using posT = TinyVector<OHMMS_PRECISION, OHMMS_DIM>;

  valT LogValue;
  DistanceTableData *d_ee;
  DistanceTableData *d_ie;

  inline void setRmax(valT x) { d_ie->setRmax(x); }

  virtual ~FakeWaveFunctionBase() {}
  virtual void evaluateLog(ParticleSet &P) = 0;
  virtual posT evalGrad(ParticleSet &P, int iat) = 0;
  virtual valT ratioGrad(ParticleSet &P, int iat, posT &grad) = 0;
  virtual valT ratio(ParticleSet &P, int iat)      = 0;
  virtual void acceptMove(ParticleSet &P, int iat) = 0;
  virtual void restore(int iat)           = 0;
  virtual void evaluateGL(ParticleSet &P) = 0;
};

struct RefWaveFunction : public FakeWaveFunctionBase
{

  using J2OrbType = TwoBodyJastrowRef<BsplineFunctor<valT>>;
  bool FirstTime;
  J2OrbType *J2;
  PooledData<valT> Buffer;

  RefWaveFunction(ParticleSet &ions, ParticleSet &els);
  ~RefWaveFunction();
  void evaluateLog(ParticleSet &P);
  posT evalGrad(ParticleSet &P, int iat);
  valT ratioGrad(ParticleSet &P, int iat, posT &grad);
  valT ratio(ParticleSet &P, int iat);
  void acceptMove(ParticleSet &P, int iat);
  void restore(int iat);
  void evaluateGL(ParticleSet &P);
};

struct SoAWaveFunction : public FakeWaveFunctionBase
{
  using J2OrbType = TwoBodyJastrow<BsplineFunctor<valT>>;

  bool FirstTime;
  J2OrbType *J2;

  SoAWaveFunction(ParticleSet &ions, ParticleSet &els);
  ~SoAWaveFunction();
  void evaluateLog(ParticleSet &P);
  posT evalGrad(ParticleSet &P, int iat);
  valT ratioGrad(ParticleSet &P, int iat, posT &grad);
  valT ratio(ParticleSet &P, int iat);
  void acceptMove(ParticleSet &P, int iat);
  void restore(int iat);
  void evaluateGL(ParticleSet &P);
};
}
#endif
