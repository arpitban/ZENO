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

#include <cmath>

#include "RandomUtilities.h"

/// Defines particles using assembly of spheres with mutable center and orientation.
///
template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator>::
RandomUtilities(RandomNumberGenerator * randomNumberGenerator):randomNumberGenerator(randomNumberGenerator){
}

template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator>::
  ~RandomUtilities() {
}

/// Generates random vectors.
///
template <class T,
        class RandomNumberGenerator>
void
RandomUtilities<T, RandomNumberGenerator>::
setRandomOnSphere(Vector3<T> * v) {
    double z1, z2, zsq;
    do  {
        z1 = 2.0 * randomNumberGenerator->getRandIn01() - 1.0;
        z2 = 2.0 * randomNumberGenerator->getRandIn01() - 1.0;
        zsq = z1 * z1 + z2 * z2;
    } while (zsq > 1.0);

    double ranh = 2.0 * sqrt(1.0 - zsq);
    v->setXYZ(z1 * ranh, z2 * ranh, 1.0 - 2.0 * zsq);
}

template <class T,
        class RandomNumberGenerator>
void
RandomUtilities<T, RandomNumberGenerator>::
setRandomInSphere(Vector3<T> * v) {
    double r = cbrt(randomNumberGenerator->getRandIn01());
    double u, w, s;
    do {
        u = 1.0 - 2.0*randomNumberGenerator->getRandIn01();
        w = 1.0 - 2.0*randomNumberGenerator->getRandIn01();
        s = u*u + w*w;
    } while(s > 1);

    double ra = 2 * r * sqrt(1 - s);
    v->setXYZ(ra * u, ra * w, r * (2 * s - 1));

}

