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

#include <vector>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"
#include <cmath>

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T>
class MCMove {
 public:
    MCMove(IntegratorMSMC<T> & integratorMSMC);
  virtual ~MCMove();
  virtual void doTrial() = 0;
  double getStepSize();
  void setStepSize();

protected:
    IntegratorMSMC & integratorMSMC;
    double stepSize;
    void adjustStepSize();
    double maxStepSize;
    long numTrials, numAccepted;
    double chiSum;
    double adjustInterval;
    int lastAdjust;
    double adjustStep, minAdjustStep;
};

template <class T>
MCMove<T>::
MCMove(IntegratorMSMC<T> & integratorMSMC) : integratorMSMC(integratorMSMC)
              {
    numTrials = 0;
    numAccepted = 0;
    chiSum = 0;
    lastAdjust = 0;
    adjustInterval = 100;
    adjustStep = 1.05;
    minAdjustStep = 1;
    verboseAdjust = false;
    tunable = true;
    maxStepSize = 0;
}

template <class T>
MCMove<T>::
  ~MCMove() {
}

template <class T>
class MCMoveTranslate : public MCMove {
public:
      MCMoveTranslate(IntegratorMSMC<T> & integratorMSMC);
      ~MCMoveTranslate();

      void doTrial();

 };

template <class T>
MCMoveTranslate<T>::
MCMoveTranslate(IntegratorMSMC<T> & integratorMSMC) : MCMove(integratorMSMC)
{
    stepSize = cbrt(integratorMSMC.getParticles()[0]->numSpheres())*integratorMSMC.getParticles()[0]->getModel().getSpheres()->at(0).getRadius();
    maxStepSize = 1000;
}

template <class T>
MCMoveTranslate<T>::
~MCMoveTranslate() {
}

template <class T>
void MCMoveTranslate<T>::
doTrial(){
    if (tunable && numTrials >= adjustInterval) {
        adjustStepSize();
    }
    double oldValue = clusterSum.value();
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> step(integratorMSMC.particles->size() - 1);
    for(int j = 0; j < (integratorMSMC.particles->size() - 1); ++j)
    {
        for(int i = 0; i < 3; ++i)
        {
            step[j].set(i, (integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * stepSize);
        }
        integratorMSMC.getParticles()[j+1]->translateBy(step[j]);
    }
    double newValue = clusterSum.value();
    double ratio = newValue / oldValue;

    if((ratio < 1) && (ratio < integratorMSMC.getRandomNumberGenerator()->getRandIn01()))
    {
        for(int j = 0; j < (integratorMSMC.particles->size() - 1); ++j)
        {
            step[j] *= -1;
            integratorMSMC.getParticles()[j+1]->translateBy(step[j]);
        }
    }
}


template <class T>
class MCMoveRotate : public MCMove {
public:
    MCMoveRotate(IntegratorMSMC<T> & integratorMSMC);
    ~MCMoveRotate();

    void doTrial();

};

template <class T>
MCMoveRotate<T>::
MCMoveRotate(IntegratorMSMC<T> & integratorMSMC) : MCMove(integratorMSMC)
{
    stepSize = M_PI/4;
    maxStepSize = M_PI/2;
}

template <class T>
MCMoveRotate<T>::
~MCMoveRotate() {
}

template <class T>
void MCMoveRotate<T>::
doTrial(){
    if (tunable && numTrials >= adjustInterval) {
        adjustStepSize();
    }
    double oldValue = clusterSum.value();
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> axis(integratorMSMC.particles->size() - 1);
    std::vector<double> angle(integratorMSMC.particles->size() - 1);
    for(int j = 0; j < (integratorMSMC.particles->size() - 1); ++j)
    {
        angle[j].set(i, (integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * stepSize);
        integratorMSMC.getRandomUtilities()->setRandomOnSphere(axis[j]);
        integratorMSMC.getParticles()[j+1]->rotateBy(axis[j], angle[j]);
    }
    double newValue = clusterSum.value();
    double ratio = newValue / oldValue;

    if((ratio < 1) && (ratio < integratorMSMC.getRandomNumberGenerator()->getRandIn01()))
    {
        for(int j = 0; j < (integratorMSMC.particles->size() - 1); ++j)
        {
            integratorMSMC.getParticles()[j+1]->rotateBy(axis[j], -angle[j]);
        }
    }
}

template <class T>
void MCMove<T>::
adjustStepSize(){
    double avg = chiSum/numTrials;
    if (avg > 0.5) {
        if (stepSize < maxStepSize) {
            if (lastAdjust < 0) {
                // back and forth
                adjustInterval *= 2;
                adjustStep = sqrt(adjustStep);
            }
            else if (lastAdjust == 5) {
                // sixth consecutive increase; increase adjustment step
                adjustStep *= adjustStep;
                if (adjustStep > 2) {
                    adjustStep = 2;
                }
                lastAdjust = 3;
            }
            stepSize *= adjustStep;
            stepSize = std::min(stepSize, maxStepSize);
            /*if (verboseAdjust) {
                printf("move increasing step size: %f (<chi> = %f)\n", stepSize, avg);
            }*/
            if (lastAdjust < 1) lastAdjust = 1;
            else lastAdjust++;
        }
        /*else if (verboseAdjust) {
            printf("move step size: %f (<chi> = %f\n)", stepSize, avg);
        }*/
    }
    else {
        if (lastAdjust > 0) {
            // back and forth
            adjustInterval *= 2;
            adjustStep = sqrt(adjustStep);
        }
        else if (lastAdjust == -5) {
            // sixth consecutive increase; increase adjustment step
            adjustStep *= adjustStep;
            if (adjustStep > 2) {
                adjustStep = 2;
            }
            lastAdjust = -3;
        }
        stepSize /= adjustStep;
        /*if (verboseAdjust) {
            printf("move decreasing step size: %f (<chi> = %f)\n", stepSize, avg);
        }*/
        if (lastAdjust > -1) lastAdjust = -1;
        else lastAdjust--;
    }
    numTrials = numAccepted = 0;
    chiSum = 0;
}

template <class T>
double MCMove<T>::
getStepSize() {
    return stepSize;
}

template <class T>
void MCMove<T>::
setStepSize(double sS) {
    stepSize = sS;
}

#endif

