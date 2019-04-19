//
// Created by arpit on 4/19/2019.
//

#ifndef ZENO_DUMMY_H
#define ZENO_DUMMY_H

#include "../Parameters.h"
#include "../Timer.h"
#include "../Geometry/Sphere.h"
#include "../Geometry/MixedModel.h"

template <class T>
class Dummy {
public:
    Dummy();
    ~Dummy();

private:
    Parameters const * parameters;
    int threadNum;
    Timer const * totalTimer;
};


#endif //ZENO_DUMMY_H
