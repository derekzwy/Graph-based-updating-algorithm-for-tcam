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

#include "rulesutils.h"
#include "rtrie.h"

void init_boundary(rule_boundary & rb)
{
    rb.field[0].low = 0;
    rb.field[0].high = 0xFFFFFFFF;
    rb.field[1].low = 0;
    rb.field[1].high = 0xFFFFFFFF;
    rb.field[2].low = 0;
    rb.field[2].high = 0xFFFF;
    rb.field[3].low = 0;
    rb.field[3].high = 0xFFFF;
    rb.field[4].low = 0;
    rb.field[4].high = 0xFF;
}

int CheckIPBounds(range fld)
{
    if (fld.low > 0xFFFFFFFF)
    {
        printf("Error: IPRange is buggy!(%llu)\n",fld.low);
        return 1;
    }
    if (fld.high > 0xFFFFFFFF)
    {
        printf("Error: IPRange is buggy!(%llu)\n",fld.high);
        return 1;
    }
    if (fld.low > fld.high)
    {
        printf("Error: IPRange is buggy!(%llu - %llu)\n",fld.low,fld.high);
        return 1;
    }
    return 0;
}

void IP2Range(unsigned ip1,unsigned ip2,unsigned ip3,unsigned ip4,unsigned iplen,pc_rule *rule,int index)
{
    unsigned tmp;
    unsigned Lo,Hi;

    if(iplen == 0){
        Lo = 0;
        Hi = 0xFFFFFFFF;

    }else if(iplen > 0 && iplen <= 8) {
        tmp = ip1 << 24;
        tmp &= (0xffffffff << (32-iplen));

        Lo = tmp;
        Hi = Lo + (1<<(32-iplen)) - 1;
    }else if(iplen > 8 && iplen <= 16){
        tmp =  ip1 << 24; 
        tmp += ip2 << 16;
        tmp &= (0xffffffff << (32-iplen));

        Lo = tmp;
        Hi = Lo + (1<<(32-iplen)) - 1;
    }else if(iplen > 16 && iplen <= 24){
        tmp = ip1 << 24; 
        tmp += ip2 << 16; 
        tmp += ip3 << 8; 
        tmp &= (0xffffffff << (32-iplen));

        Lo = tmp;
        Hi = Lo + (1<<(32-iplen)) - 1;
    }else if(iplen > 24 && iplen <= 32){
        tmp = ip1 << 24; 
        tmp += ip2 << 16; 
        tmp += ip3 << 8; 
        tmp += ip4;
        tmp &= (0xffffffff << (32-iplen));

        Lo = tmp;
        Hi = Lo + (1<<(32-iplen)) - 1;
    }else{
        printf("Error: Src IP length exceeds 32\n");
        exit(1);
    }

    rule->field[index].low  = Lo;
    rule->field[index].high = Hi;

    if (CheckIPBounds(rule->field[index]))
    {
        printf("Error: IP2Range bounds check for %d failed\n",index);
        exit(1);
    }

    //if(rule->field[1].low == 930850816 )
    //    printf("here\n");
    //printf("\t Prefix: %u.%u.%u.%u/%u\n",ip1,ip2,ip3,ip4,iplen);
    //printf("\t Range : %llu : %llu\n",rule->field[index].low,rule->field[index].high);

}

int loadrules(FILE *fp, vector<pc_rule> &classifier) {
    int i = 0;
    int wild = 0;
    unsigned sip1, sip2, sip3, sip4, siplen;
    unsigned dip1, dip2, dip3, dip4, diplen;
    unsigned proto, protomask;
    unsigned junk, junkmask;

    pc_rule rule;

    while(1) {
        wild = 0;
        /*if(fscanf(fp,"@%u.%u.%u.%u/%u\t%u.%u.%u.%u/%u\t%llu : %llu\t%llu : %llu\t%x/%x\t%x/%x\n",
                    &sip1, &sip2, &sip3, &sip4, &siplen, &dip1, &dip2, &dip3, &dip4, &diplen, 
                    &rule.field[2].low, &rule.field[2].high, &rule.field[3].low, &rule.field[3].high,
                    &proto, &protomask, &junk, &junkmask) != 18) break;*/

        if(fscanf(fp,"@%u.%u.%u.%u/%u %u.%u.%u.%u/%u %llu : %llu %llu : %llu %x/%x\n",
                    &sip1, &sip2, &sip3, &sip4, &siplen, &dip1, &dip2, &dip3, &dip4, &diplen, 
                    &rule.field[2].low, &rule.field[2].high, &rule.field[3].low, &rule.field[3].high,
                    &proto, & protomask) != 16) break;
        rule.siplen = siplen;
        rule.diplen = diplen;
        rule.sip[0] = sip1;
        rule.sip[1] = sip2;
        rule.sip[2] = sip3;
        rule.sip[3] = sip4;
        rule.dip[0] = dip1;
        rule.dip[1] = dip2;
        rule.dip[2] = dip3;
        rule.dip[3] = dip4;

        IP2Range(sip1,sip2,sip3,sip4,siplen,&rule,0);
        IP2Range(dip1,dip2,dip3,dip4,diplen,&rule,1);

        if(protomask == 0xFF){
            rule.field[4].low = proto;
            rule.field[4].high = proto;
        }else if(protomask == 0){
            rule.field[4].low = 0;
            rule.field[4].high = 0xFF;
            wild++;
        }else{
            printf("Protocol mask error\n");
            return 0;
        }
        rule.priority = i;
        rule.ruleid = i+1;
        if ((rule.field[0].low == 0) && (rule.field[0].high == 0xffffffff)) {
            wild++;
        }
        if ((rule.field[1].low == 0) && (rule.field[1].high == 0xffffffff)) {
            wild++;
        }
        if ((rule.field[2].low == 0) && (rule.field[2].high == 65535)) {
            wild++;
        }
        if ((rule.field[3].low == 0) && (rule.field[3].high == 65535)) {
            wild++;
        }
        if (wild != 5) {
            classifier.push_back(rule);
            i++;
        }
    }
    return i;
}

range range_in_boundary_1D(range r, range boundary)
{
    range ret;
    if (r.low > boundary.low) {
        ret.low = r.low; 
    }
    else {
        ret.low = boundary.low;
    }

    if (r.high < boundary.high) {
        ret.high = r.high;
    }
    else {
        ret.high = boundary.high;
    }
    return ret;
}


void dump_rule(pc_rule *i, FILE *out)
{
    fprintf(out, "@%d.%d.%d.%d/%d\t%d.%d.%d.%d/%d\t%lld : %lld\t%lld : %lld\t",
            (i)->sip[0], (i)->sip[1], (i)->sip[2], (i)->sip[3], (i)->siplen, \
            (i)->dip[0], (i)->dip[1], (i)->dip[2], (i)->dip[3], (i)->diplen, \
            (i)->field[2].low, (i)->field[2].high, \
            (i)->field[3].low, (i)->field[3].high);

    if((i)->field[4].high == 0xFF && (i)->field[4].low == 0)
        fprintf(out, "0x%x/0x%x\t",0,0);
    else
       fprintf(out, "0x%llx/0x%llx\t", (i)->field[4].low, 0xFFULL);

    fprintf(out, "0x%x/0x%x\n", 0,0);

}

void dump_rules(list<pc_rule*> ruleset, string str)
{
    FILE *fp = fopen(str.c_str(), "w");

    for (list<pc_rule *>::iterator i = ruleset.begin();
            i != ruleset.end();
            i++){
        dump_rule(*i, fp);
    }
    fclose(fp);
} 

bool is_equal(pc_rule rule1, pc_rule rule2, rule_boundary boundary)
{
    int count = 0;
    range r1, r2;
    for (int i = 0;i < MAXDIMENSIONS;i++)
    {
        //if (rule1.field[i].low > boundary.field[i].low) {
        //    r1.low = rule1.field[i].low;
        //} else {
        //    r1.low = boundary.field[i].low;
        //}
        //if (rule1.field[i].high < boundary.field[i].high) {
        //    r1.high = rule1.field[i].high;
        //} else {
        //    r1.high = boundary.field[i].high;
        //}

        r1 = range_in_boundary_1D(rule1.field[i], boundary.field[i]);

        //if (rule2.field[i].low > boundary.field[i].low) {
        //    r2.low = rule2.field[i].low;
        //} else {
        //    r2.low = boundary.field[i].low;
        //}
        //if (rule2.field[i].high < boundary.field[i].high) {
        //    r2.high = rule2.field[i].high;
        //} else {
        //    r2.high = boundary.field[i].high;
        //}

        r2 = range_in_boundary_1D(rule2.field[i], boundary.field[i]);

        if (r1.low <= r2.low && r1.high >= r2.high)

        {
            count++;
        }
    }

    if (count == MAXDIMENSIONS)
        return true;
    else
        return false;
}

bool myequal(pc_rule* first, pc_rule* second)
{
    return (first->priority == second->priority);
}

void init_rnode(rnode *n, rule_boundary &b) 
{
    rg sip(b.field[0].low, b.field[0].high);
    rg dip(b.field[1].low, b.field[1].high);

    n[0].setb(sip);
    n[1].setb(dip);

}


void remove_redund_rt(vector<pc_rule*> &pr, rule_boundary &b, bool large) 
{
    vector<pc_rule*> rulelist;
    //vector<pc_rule*> vpc;
    //l2v(pr, vpc);
    
    range br;
    
    rnode rt[2];
    vector<pc_rule*> set[2];
    //map<int, int> sect;
    vector<pc_rule*> sect;

    init_rnode(rt, b);
    //int count = 0;

    for(auto rule = pr.begin();
            rule != pr.end(); ++ rule) {
        //count ++;
        //cout<<count<<endl;

        //if(large) {
        //    if(!check_range_size((*rule)->field[0], b, 0) 
        //            || !check_range_size((*rule)->field[1], b, 1)) {
        //        rulelist.push_back(*rule);
        //        continue;
        //    }
        //}

        //if((*rule)->priority == 567) {
        //    cout<<"here"<<endl;
        //}

        sect.clear();
        for(int i = 0; i < 2; i++) {
            br = range_in_boundary_1D((*rule)->field[i], b.field[i]); 
            set[i].clear();
            rt_qry_insert(rt+i, rg(br.low, 
                    br.high), set[i], (*rule));
        }

        //if(set[0].size() > 5 * set[1].size()) {
        //    sect = set[1];
        //}
        //else if (set[1].size() > 5 * set[0].size()) {
        //    sect = set[0];

        //}
        //else{ 
        sort(set[0].begin(), set[0].end());
        sort(set[1].begin(), set[1].end());
        set_intersection(set[0].begin(), set[0].end(),
                set[1].begin(), set[1].end(),
                back_inserter(sect));
        //}

        if(sect.empty()) {
            rulelist.push_back(*rule);
        }
        else {
            bool found = false;
            for(auto check_rule = sect.begin();
                    check_rule != sect.end();
                    check_rule++) {
                if(is_equal(**check_rule, **rule, b)) {
                    found = true;
                    //count++;
                    //cout<<(*rule)->priority<<" was hide by " <<*check_rule<<endl;
                    break;
                }
            }

            if(!found) {
                rulelist.push_back(*rule);
            }
        }

    }

    rt_destory(rt[0].left);
    rt_destory(rt[0].right);
    rt_destory(rt[1].left);
    rt_destory(rt[1].right);

    pr = move(rulelist);
    //cout<<"remove "<<count<<endl;
}


void remove_redund(list<pc_rule*> &pr, rule_boundary &rb)
{
    list <pc_rule*> rulelist;
    int i = 0;
    //list <pc_rule*> dellist;
    //list <pc_rule*> cftlist;

    rulelist.clear();

    for (list<pc_rule*>::iterator rule = pr.begin();
            rule != pr.end();++rule)
    {
        int found = 0;
        for (list<pc_rule*>::iterator mule = rulelist.begin();
                mule != rulelist.end();++mule)
        {
            if (is_equal(**mule,**rule, rb) == true)
            {
                found = 1;
                //cout<< (**rule).priority <<" hide by" <<(**mule).priority<<endl;
                i++;
                //dellist.push_back(*rule);
                //cftlist.push_back(*mule);
                break;
            }
        }
        if (found != 1)
        {
            rulelist.push_back(*rule);
        }
    }
    // Now add back the rules 
    pr.clear();
    pr = rulelist;
    //pr.unique(myequal);
    //cout<<"remove "<<i<<endl;
    
    //dump_rules(dellist, "dellist");
}


void load_rule_ptr(vector <pc_rule> &rule_list,list <pc_rule*> &ruleptr_list,int start,int end)
{
    printf("Rule:%d - %d\n",start,end);
    int count = 0;
    for (vector <pc_rule>::iterator i = rule_list.begin();i != rule_list.end();++i) 
    {
        if (count >= start && count <= end)
            ruleptr_list.push_back(&(*i));
        count++;
    }
}


int load_ft(FILE *fpt, field_type *ft) 
{
    int ret;
    static int cnt = 0;
    unsigned int junk, junkmask;
    if(fpt == NULL) 
        return 0;

    ret = fscanf(fpt, "%u\t%u\t%u\t%u\t%u\t%u\t%u\n", &(ft[0]),
            &(ft[1]),
            &(ft[2]),
            &(ft[3]),
            &(ft[4]),
            &junk,
            &junkmask);
    
    if(ret != 7)
        return 0;
    cnt ++;
    return cnt;

}

int check_rule(pc_rule *r, field_type *ft)
{
    if(ft[0] <= r->field[0].high && ft[0] >= r->field[0].low &&
            ft[1] <= r->field[1].high && ft[1] >= r->field[1].low &&
            ft[2] <= r->field[2].high && ft[2] >= r->field[2].low &&
            ft[3] <= r->field[3].high && ft[3] >= r->field[3].low &&
            ft[4] <= r->field[4].high && ft[4] >= r->field[4].low) {
        return 1;
    }
    else {
        return -1;
    }

}


int linear_search(list<pc_rule*> &p_classifier, field_type *ft)
{
    int ret = -1;
    for(list<pc_rule*>::iterator it = p_classifier.begin();
            it != p_classifier.end();
            it++) {
        ret = check_rule(*it, ft);
        if(ret == 1) {
            ret = (*it)->priority;
            break;
        }
    }
    return ret;
}
 
void dump_rule_with_priority(pc_rule *i, FILE *out)
{
    fprintf(out, "@%d.%d.%d.%d/%d\t%d.%d.%d.%d/%d\t%lld : %lld\t%lld : %lld\t",
            (i)->sip[0], (i)->sip[1], (i)->sip[2], (i)->sip[3], (i)->siplen, \
            (i)->dip[0], (i)->dip[1], (i)->dip[2], (i)->dip[3], (i)->diplen, \
            (i)->field[2].low, (i)->field[2].high, \
            (i)->field[3].low, (i)->field[3].high);

    if((i)->field[4].high == 0xFF && (i)->field[4].low == 0)
        fprintf(out, "0x%x/0x%x\t",0,0);
    else
       fprintf(out, "0x%llx/0x%llx\t", (i)->field[4].low, 0xffull);

    fprintf(out, "0x%x/0x%x\t", 0,0);
    fprintf(out, "%d\n", i->priority);

}

void dump_rules_with_priority(list<pc_rule*> ruleset, string file)
{
    FILE *fp = fopen(file.c_str(), "w");

    for (list<pc_rule *>::iterator i = ruleset.begin();
            i != ruleset.end();
            i++){
        dump_rule_with_priority(*i, fp);
    }
    fclose(fp);
} 

#define IPBIN 0.05
#define BIN 0.5

int check_range_size(const range & check, int index)
{
    double field;
    if(index == 0 || index == 1) {
        field = (((double)(check.high - check.low))/0xFFFFFFFF);
        if(field >= IPBIN) {
            return 1;
        }
    } 
    else if (index == 2 || index == 3) {
        field = (((double)(check.high - check.low))/65536);
        if(field >= BIN) {
            return 1;
        }
    }
    else {
        field = (((double)(check.high - check.low))/256);
        if(field >= BIN) {
            return 1;
        }
    }
//usually the Protocol is a fixed value
    return 0;
}


void init_overlap_trees(struct overlap_trees *ots, rule_boundary &rb, vector<pc_rule*> pc) 
{
    for(int i = 0; i < MAXDIMENSIONS; i++) {
        rg r(rb.field[i].low, rb.field[i].high);
        ots->roots[i].setb(r);
    }


    for(auto rule = pc.begin(); rule != pc.end(); rule++) {
        for(int i = 0; i < MAXDIMENSIONS; i++) 
            rt_insert(&(ots->roots[i]), rg((*rule)->field[i].low, (*rule)->field[i].high), *rule);
    }

}

bool overlap(pc_rule *r1, pc_rule *r2) 
{
    bool o = true;
    int i = 0;
    for(;i<MAXDIMENSIONS;i++) {
        o = o && !((r1->field[i].low > r2->field[i].high) || (r1->field[i].high < r2->field[i].low)); 
    }
    return o;
}

void find_overlap_rules_slow(vector<pc_rule*> &pc, pc_rule* rule, int l, int h, vector<int> &set)
{
    for(int i =l ; i < h; i++) {
        if(overlap(pc[i], rule)) {
            set.push_back(pc[i]->priority);
        }
    }

}



void range2prefix(range r, vector<range> &prefixes, uint32_t lo, uint32_t hi)
{
    if(r.high == hi && r.low == lo) {
        prefixes.push_back(r);
        return;
    }
    uint32_t mid = lo + (hi - lo)/2;

    if(r.low <= mid) {
        uint32_t hiend = min(mid, (uint32_t)r.high);
        range2prefix(range(r.low, hiend), prefixes, lo, mid); 
    }

    if(r.high > mid) {
        uint32_t loend = max((uint32_t)r.low, mid+1);
        range2prefix(range(loend, r.high), prefixes, mid+1, hi);
    }
}

void range2prefix(range r, vector<range> &prefixes) 
{
    uint32_t lo = 0;
    uint32_t hi = 0xffffffffUL;
    

    range2prefix(r, prefixes, lo, hi);
}

void newrule_from_ranges(pc_rule *rule, pc_rule *oldrule) 
{ 
    for(int i = 0; i< 4; i++) {
        rule->sip[i] = (rule->field[0].low >>(3-i)*8) & 0xff;
        rule->dip[i] = (rule->field[1].low >>(3-i)*8) & 0xff;
    }
    rule->siplen = 32 - __builtin_popcount(rule->field[0].high - rule->field[1].low);
    rule->diplen = 32 - __builtin_popcount(rule->field[1].high - rule->field[1].low);
    rule->priority = oldrule->priority;
    rule->ruleid = oldrule->ruleid;

}

void extend_rules_dfs(pc_rule *rule, vector<range> &r, vector<pc_rule> &out, int depth)
{

    if (depth == MAXDIMENSIONS)  {
        pc_rule newrule;
        //cout<<"insert ";
        for(int i = 0; i< MAXDIMENSIONS; i++) {
            newrule.field[i] = r[i];
            //cout<<r[i].low<<" "<<r[i].high<<" ";
        }
        newrule_from_ranges(&newrule, rule);
        //cout<<endl;
        out.push_back(newrule);
        return;
    }

    vector<range> curr;
    range2prefix(rule->field[depth], curr); 

    for(size_t i = 0; i< curr.size(); i++) {
        r.push_back(curr[i]);
        extend_rules_dfs(rule, r, out, depth + 1); 
        r.pop_back();
    }

}


void extend_rules(vector<pc_rule*> &in, vector<pc_rule> &out) 
{
    vector<range> r;
    cout<<"Input ruleset: "<<in.size()<<endl;
    for(size_t i = 0; i < in.size(); i++) {
        int previous = out.size();
        extend_rules_dfs(in[i], r, out, 0);
        in[i]->expand = &out[previous];
    }
}

void find_overlap_rules(struct overlap_trees *ots, vector<pc_rule*> &pc, pc_rule * rule, int l, int h, vector<int> &rset) 
{
    int count = 0;
    unordered_map<pc_rule*, int> unionset;

    vector<pc_rule*> fset[MAXDIMENSIONS];
    int v[MAXDIMENSIONS];
    for(int i = 0; i< MAXDIMENSIONS;i ++) {
        v[i] = 0;
        if(!check_range_size(rule->field[i], i)) {
            if(ots->roots[i].rlist.size() < pc.size() - rule->priority) {
                v[i] = 1;
                count++;
            }
        }
    }

    if(count == 0) {
        //fall back to original check 
        find_overlap_rules_slow(pc, rule, l, h, rset);
        return ;
    }

    for(int i = 0; i < MAXDIMENSIONS; i++) {
        if(v[i]) {
            rt_query_or(&(ots->roots[i]), rg(rule->field[i].low, rule->field[i].high), fset[i], l, h); 
            for(auto orule = fset[i].begin(); orule != fset[i].end(); orule++) {
                unionset[*orule] ++;
            }
        }
    }

    for(auto pr = unionset.begin(); pr != unionset.end(); pr++) {
        if(pr->second == count && overlap(pr->first,rule)) {
            rset.push_back(pr->first->priority);
        }
    }
    //cout<<rset.size()<<endl;

}




