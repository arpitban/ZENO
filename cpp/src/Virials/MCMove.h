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
    MCMove(IntegratorMSMC<T> & integratorMsmc);
  virtual ~MCMove();

  virtual void doTrial() = 0;

protected:
    IntegratorMSMC & integratorMsmc;
    double stepSize;
};

template <class T>
MCMove<T>::
MCMove(IntegratorMSMC<T> & integratorMsmc) : integratorMsmc(integratorMsmc)
              {
}

template <class T>
MCMove<T>::
  ~MCMove() {
}

<T>
class MCMoveTranslate : public MCMove {
public:
      MCMoveTranslate(IntegratorMSMC<T> & integratorMsmc);
      ~MCMoveTranslate();

      void doTrial();

 };

template <class T>
MCMoveTranslate<T>::
MCMoveTranslate(IntegratorMSMC<T> & integratorMsmc) : MCMove(integratorMsmc)
{
    stepSize = cbrt(integratorMsmc.getMolecules()[0]->numSpheres())*integratorMsmc.getMolecules()[0]->getModel().getSpheres()->at(0).getRadius();

}

template <class T>
MCMoveTranslate<T>::
~MCMoveTranslate() {
}

<T>
void MCMoveTranslate<T>::
doTrial(){
    double oldValue = clusterSum.value();
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> step(integratorMsmc.molecules->size() - 1);
    for(int j = 0; j < (integratorMsmc.molecules->size() - 1); ++j)
    {
        for(int i = 0; i < 3; ++i)
        {
            step[j].set(i, (integratorMsmc.getRandomNumberGenerator()->getRandIn01() - 0.5) * stepSize);
        }
        integratorMsmc.getMolecules()[j+1]->translateBy(step[j]);
    }
    double newValue = clusterSum.value();
    double ratio = newValue / oldValue;

    if((ratio < 1) && (ratio < integratorMsmc.getRandomNumberGenerator()->getRandIn01()))
    {
        for(int j = 0; j < (integratorMsmc.molecules->size() - 1); ++j)
        {
            step[j] *= -1;
            integratorMsmc.getMolecules()[j+1]->translateBy(step[j]);
        } 
    }
}
/// Defines particles using assembly of spheres with mutable center and orientation.
///

#endif

