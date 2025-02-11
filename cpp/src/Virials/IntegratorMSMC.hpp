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

#include "IntegratorMSMC.h"


///  Constructs the class to perform a step.
///

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
IntegratorMSMC(Parameters const * parameters,
                int threadNum,
                Timer const * totalTimer,
                RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> *> & models,
                OverlapTester<T> const & overlapTester,
                std::vector<MCMove<T, RandomNumberGenerator>> & mcMoves,
                std::vector<double> & moveProbs,
                ClusterSum<T, RandomNumberGenerator> * clusterSum,
                MeterOverlap<T, RandomNumberGenerator> & meterOverlap) :
              parameters(parameters),
              threadNum(threadNum),
              totalTimer(totalTimer),
              randomNumberGenerator(randomNumberGenerator),
              boundingSpheres(boundingSpheres),
              numParticles(numParticles),
              overlapTester(overlapTester),
              mcMoves(mcMoves),
              moveProbs(moveProbs),
              clusterSum(clusterSum),
              meterOverlap(meterOverlap){
    for(int i = 0; i < numParticles.size(); ++i)
    {
        for (int j =0; j < numParticles[i]; ++j)
        {
            particles.push_back(new Particle<T> (models[i], boundingSpheres[i]));
        }
    }

}

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
  ~IntegratorMSMC() {

}

/// Carries out a step. A step performed by the integrator consists of selecting a MCMove, performing the trial defined by MCMove and collecting data and statistics.
///
template <class T,
  class RandomNumberGenerator>
void
IntegratorMSMC<T,
               RandomNumberGenerator>::
doStep(){
    double random = randomNumberGenerator->getRandIn01();
    double cumProb = 0.0;
    int currentMove = -1;
    for(int i = 0; i < mcMoves.size(); ++i)
    {
        cumProb += moveProbs[i];
        if(random < cumProb)
        {
            currentMove = i;
        }
        mcMoves[i].doTrial();
    }
    meterOverlap.collectData();
}

/// Returns particles.
///
template <class T,
        class RandomNumberGenerator>
std::vector<Particle<T> *>
IntegratorMSMC<T,RandomNumberGenerator>::
getParticles(){
    return particles;
}

/// Returns random number generator.
///
template <class T,
        class RandomNumberGenerator>
RandomNumberGenerator *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomNumberGenerator(){
    return randomNumberGenerator;
}

/// Returns cluster sum.
///
template <class T,
        class RandomNumberGenerator>
ClusterSum<T, RandomNumberGenerator> *
IntegratorMSMC<T,RandomNumberGenerator>::
getClusterSum(){
    return clusterSum;
}

/// Returns random utilities.
///
template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator> *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomUtilities(){
    return & randomUtilities;
}

