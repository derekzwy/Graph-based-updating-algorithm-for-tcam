#ifndef _TREE_CAM_H_
#define _TREE_CAM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <list>
#include <map>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <vector>
#include <arpa/inet.h>
#include <string.h>


#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>



#include "range.h"
#include "rulesutils.h"



struct Node 
{
    static int label;
    int no;
    int dim;
    int dimmask;
    int isLeaf;
    int rulecnt;
    int depth;
    rule_boundary rb;

    std::vector<pc_rule*> pc;
    std::vector<range> segs;
    std::vector<Node *> childnodes;

    Node(): dim(-1), dimmask(0xff), isLeaf(-1), rulecnt(0), depth(-1)
    {
        no = Node::label ++;
        init_boundary(rb);
    }
};

int Node::label = 0;






#endif
