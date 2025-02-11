// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#ifndef MC_MOVE_H
#define MC_MOVE_H

#include "IntegratorMSMC.h"

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T,
        class RandomNumberGenerator>
class IntegratorMSMC;

template <class T,
        class RandomNumberGenerator>
class MCMove {
 public:
    MCMove(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC);
  virtual ~MCMove();
  virtual void doTrial() = 0;
  double getStepSize();
  void setStepSize(double sS);

protected:
    IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC;
    double stepSize;
    void adjustStepSize();
    double maxStepSize;
    long numTrials, numAccepted;
    double chiSum;
    double adjustInterval;
    int lastAdjust;
    double adjustStep, minAdjustStep;
    bool verboseAdjust, tunable;
    
};

template <class T,
        class RandomNumberGenerator>
class MCMoveTranslate : public MCMove<T, RandomNumberGenerator> {
public:
      MCMoveTranslate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC);
      ~MCMoveTranslate();
      void doTrial();
 };

template <class T,
        class RandomNumberGenerator>
class MCMoveRotate : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveRotate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC);
    ~MCMoveRotate();
    void doTrial();
};

template <class T,
        class RandomNumberGenerator>
class MCMoveChainVirial : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveChainVirial(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, double sigma);
    ~MCMoveChainVirial();
    void doTrial();
protected:
    double sigma;
};

#endif

