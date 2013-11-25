/* 
 * File:   Point.cpp
 * Author: skiloop
 * 
 * Created on 2013年11月22日, 下午4:43
 */

#include <cstdlib>
#include "Point.h"

Point::Point(unsigned i, unsigned j, unsigned k)
: x(i), y(j), z(k) {
}

Point::Point(const Point& orig)
: x(orig.x)
, y(orig.y)
, z(orig.z) {

}

Point::~Point() {
}

bool Point::checkMax(const Point&maxPoint) const {
    return (x < maxPoint.x && y < maxPoint.y && z < maxPoint.z);
}

bool Point::checkMax(unsigned xMax, unsigned yMax, unsigned zMax) const {
    return (x < xMax && y < yMax && z < zMax);
}

void Point::operator=(Point const &obj) {
    if (this != &obj) {
        this->x = obj.x;
        this->y = obj.y;
        this->z = obj.z;
    }
}

void Point::setValue(unsigned i, unsigned j, unsigned k) {
    x = i;
    y = j;
    z = k;
}

