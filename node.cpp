#include "node.h"

node::node()
{
    U=0;
    Q=0;
}


void node::swap(){
    float temp=U;
    U=Q;
    Q=temp;
}
