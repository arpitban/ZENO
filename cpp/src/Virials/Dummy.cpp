//
// Created by arpit on 4/19/2019.
//

#include "Dummy.h"

template <class T>
Dummy<T>::Dummy():parameters(nullptr),  threadNum(0),totalTimer(nullptr){

}
template <class T>
Dummy<T>::~Dummy() {}

template class Dummy<double>;