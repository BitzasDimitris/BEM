#ifndef POINT_H
#define POINT_H
#include <math.h>
#include <type_traits>

class Point
{
public:
    Point();
    Point(float x,float y);

    Point operator+(Point p);
    Point& operator+=(Point p);
    Point operator-(Point p);
    Point operator-();
    Point operator*(Point p);
    template <typename T>
    Point operator*(T f){
        return Point(this->x*f,this->y*f);
    }
    Point operator/(float f);
    void operator=(Point p);
    float mag();
    float distance(Point p);

    float x,y;
};

#endif // POINT_H
