#include "point.h"

Point::Point()
{
    x=0;
    y=0;
}

Point::Point(float x//a
             , float y//b
             )
{
    this->x=x;
    //x=a;
    this->y=y;
    //y=b;
}


Point Point::operator+(Point p){
    return Point(this->x+p.x,this->y+p.y);
}

Point& Point::operator+=(Point p){
    this->x+=p.x;
    this->y+=p.y;
    return *this;
}

Point Point::operator-(Point p){
    return Point(this->x-p.x,this->y-p.y);
}

Point Point::operator-(){
    return Point(-this->x,-this->y);
}

Point Point::operator*(Point p){
    return Point(this->x*p.x+this->y*p.y,-this->x*p.y+this->y*p.x);
}

Point Point::operator/(float f){
    return Point(this->x/f,this->y/f);
}

void Point::operator=(Point p){
    this->x=p.x;
    this->y=p.y;
}

float Point::mag(){
    return sqrt(pow(this->x,2)+pow(this->y,2));
}

float Point::distance(Point p){
    return (Point(x,y)-p).mag();
}
