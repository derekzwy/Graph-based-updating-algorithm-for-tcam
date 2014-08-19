#ifndef __UPDATE_H__
#define __UPDATE_H__
#include <vector>

#include "rulesutils.h"


struct node {
    int index;
    int move;
    int cost;
    bool valid;
    int label;
    struct pc_rule *r;
    std::vector<struct node *> in;
    std::vector<struct node *> out;

    node () : index(-1), move(-1), cost(0),  valid(true), label(-1),  r(nullptr){
        
    }
};





#endif
