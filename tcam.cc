#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <list>
#include <map>
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
#include <unistd.h>
#include "tcam.h"

using namespace std;

tcam::tcam(vector<pc_rule> &rs) 
{
    pc = rs; 
    valid.resize(rs.size(), true);
}

tcam::tcam() 
{
    
}

int tcam::match(field_type *p) {
    int ret = -1;
    for(size_t i = 0; i < pc.size(); i++) {
        ret = check_rule(&pc[i], p) && valid[i];
        if( ret == 1) {
            ret = i; 
            break;
        }
    }

    return ret;
}

void tcam::write(int i, pc_rule &rule) 
{
    if(i > (int)pc.size()-1) {
        pc.push_back(rule);
        valid.push_back(true);
    }
    pc[i]  = rule;
}

void tcam::del(int i)
{
    if(i > (int)pc.size() -1) 
        return;
    valid[i] = false;
}



