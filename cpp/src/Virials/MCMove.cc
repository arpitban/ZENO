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

#include "MCMove.h"

/// Defines particles using assembly of spheres with mutable center and orientation.
///

template <class T,
        class RandomNumberGenerator>
MCMove<T, RandomNumberGenerator>::
MCMove(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC) : integratorMSMC(integratorMSMC)
              {
    stepSize = 0;
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

template <class T,
        class RandomNumberGenerator>
MCMove<T, RandomNumberGenerator>::
  ~MCMove() {
}

template <class T, 
        class RandomNumberGenerator>
MCMoveTranslate<T, RandomNumberGenerator>::
MCMoveTranslate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC) : MCMove<T, RandomNumberGenerator>(integratorMSMC)
{
    MCMove<T, RandomNumberGenerator>::stepSize = cbrt(integratorMSMC.getParticles()[0]->numSpheres())*integratorMSMC.getParticles()[0]->getModel().getSpheres()->at(0).getRadius();
    MCMove<T, RandomNumberGenerator>::maxStepSize = 1000;
}

template <class T,
        class RandomNumberGenerator>
MCMoveTranslate<T, RandomNumberGenerator>::
~MCMoveTranslate() {
}

template <class T,
        class RandomNumberGenerator>
void MCMoveTranslate<T, RandomNumberGenerator>::
doTrial(){
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    double oldValue = MCMove<T, RandomNumberGenerator>::clusterSum.value();
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> step(MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1);
    for(int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1); ++j)
    {
        for(int i = 0; i < 3; ++i)
        {
            MCMove<T, RandomNumberGenerator>::step[j].set(i, (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize);
        }
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j+1]->translateBy(step[j]);
    }
    double newValue = MCMove<T, RandomNumberGenerator>::clusterSum.value();
    double ratio = newValue / oldValue;

    if((ratio < 1) && (ratio < MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01()))
    {
        for(int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1); ++j)
        {
            MCMove<T, RandomNumberGenerator>::step[j] *= -1;
            MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j+1]->translateBy(step[j]);
        }
    }
}

template <class T,
        class RandomNumberGenerator>
MCMoveRotate<T, RandomNumberGenerator>::
MCMoveRotate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC) : MCMove<T, RandomNumberGenerator>(integratorMSMC)
{
    MCMove<T, RandomNumberGenerator>::stepSize = M_PI/4;
    MCMove<T, RandomNumberGenerator>::maxStepSize = M_PI/2;
}

template <class T, class RandomNumberGenerator>
MCMoveRotate<T, RandomNumberGenerator>::
~MCMoveRotate() {
}

template <class T, class RandomNumberGenerator>
void MCMoveRotate<T, RandomNumberGenerator>::
doTrial(){
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    double oldValue = MCMove<T, RandomNumberGenerator>::clusterSum.value();
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> axis(MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1);
    std::vector<double> angle(MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1);
    for(int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1); ++j)
    {
        angle[j] =  (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize;
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomUtilities()->setRandomOnSphere(axis[j]);
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j+1]->rotateBy(axis[j], angle[j]);
    }
    double newValue = MCMove<T, RandomNumberGenerator>::clusterSum.value();
    double ratio = newValue / oldValue;

    if((ratio < 1) && (ratio < MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01()))
    {
        for(int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1); ++j)
        {
            MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j+1]->rotateBy(axis[j], -angle[j]);
        }
    }
}

template <class T, class RandomNumberGenerator>
MCMoveChainVirial<T, RandomNumberGenerator>::
MCMoveChainVirial(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, double sigma): MCMove<T, RandomNumberGenerator>(integratorMSMC), sigma(sigma)
{
}

template <class T, class RandomNumberGenerator>
MCMoveChainVirial<T, RandomNumberGenerator>::
~MCMoveChainVirial() {
}

template <class T, class RandomNumberGenerator>
void MCMoveChainVirial<T, RandomNumberGenerator>::
doTrial() {
    const Vector3<T> rPrev = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[0]->getCenter();
    Vector3<T> sPrev = rPrev;
    for(int j = 1; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.particles->size() - 1); ++j)
    {
        const Vector3<T> r = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j]->getCenter();
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomUtilities()->setRandomInSphere(r);
        Vector3<T> s = r*sigma + sPrev;
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()[j]->setCenter(s);
         sPrev = s;
    }
}

template <class T, class RandomNumberGenerator>
void MCMove<T, RandomNumberGenerator>::
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

template <class T, class RandomNumberGenerator>
double MCMove<T, RandomNumberGenerator>::
getStepSize() {
    return stepSize;
}

template <class T, class RandomNumberGenerator>
void MCMove<T, RandomNumberGenerator>::
setStepSize(double sS) {
    stepSize = sS;
}


