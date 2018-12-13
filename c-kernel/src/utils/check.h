#ifndef CHECK_H
#define CHECK_H

#include <typeinfo>
#include <stdio.h>

template<typename T>
bool type_equal(T val, T type) {
	return typeid(val).name() == typeid(type).name();
}

#endif