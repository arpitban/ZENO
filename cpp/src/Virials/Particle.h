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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cmath>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T>
class Particle {
 public:
    Particle(MixedModel<T> & model,
            Sphere<T> & boundingSphere);
  ~Particle();

  int numSpheres();
  const Vector3<T> getCenter() const;
  void setCenter(Vector3<T> v);
  void translateBy(Vector3<T> step);
  void rotateBy(Vector3<T> axis, T angle);
  const Vector3<T> setFromSpherePosition(int index) const;
  MixedModel<T> * getModel();
  Sphere<T> * getBoundingSphere();

 private:
    MixedModel<T> & model;
    Matrix3x3<T> orientation{1,0,0,0,1,0,0,0,1};
    Vector3<T> center{0,0,0};
    Sphere<T> & boundingSphere;
};

template <class T>
Particle<T>::
Particle(MixedModel<T> & model, Sphere<T> & boundingSphere) : model(model), boundingSphere(boundingSphere)
              {
}

template <class T>
Particle<T>::
  ~Particle() {
}

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T>
int
Particle<T>::
numSpheres(){
    return model.getSpheres() -> size();
}

template <class T>
const Vector3<T>
Particle<T>::
getCenter() const {
    return center;
}

template <class T>
void
Particle<T>::
setCenter(Vector3<T> v) {
    center = v;
}

template <class T>
void
Particle<T>::
translateBy(Vector3<T> step){
    center += step;
}

template <class T>
void
Particle<T>::
rotateBy(Vector3<T> axis, T angle){
    Matrix3x3<T> rotation;
    rotation.setAxisAngle(axis, angle);
    rotation.transform(orientation);
}

template <class T>
const Vector3<T>
Particle<T>::
setFromSpherePosition( int index) const {
    Vector3<T> position = model.getSpheres() -> at(index).getCenter();
    orientation.transform(position);
    position += center;
    return position;
}

template <class T>
MixedModel<T> *
Particle<T>::
getModel(){
    return model;
}

template <class T>
Sphere<T> *
Particle<T>::
getBoundingSphere(){
    return boundingSphere;
}
#endif

