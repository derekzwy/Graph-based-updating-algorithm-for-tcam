#ifndef __TCAM_H__
#define __TCAM_H__


#include<vector>
#include "rulesutils.h"

class tcam {
public:
    std::vector<pc_rule> pc;
    std::vector<bool> valid;
    int match(field_type *p);
    void write(int i, pc_rule &rule);
    void del(int i);
    tcam(vector<pc_rule> &rs);
    tcam();
};




#endif
