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
#include <queue>

#include "update.h"
#include "tcam.h"

using namespace std;

#include "rulesutils.h"
#include "rtrie.h"

FILE *fpr;
FILE *fpt;
int reorder_cnt = 0;

void parseargs(int argc, char *argv[]) {
    int	c;
    while ((c = getopt(argc, argv, "r:t:")) != -1) {
        switch (c) {
            case 'r':
                fpr = fopen(optarg, "r");
                break;
            case 't':
                fpt = fopen(optarg, "r");
                break;
            default:
                break;
        }
    }

    
    if(fpr == NULL){
        printf("can't open ruleset file\n");
        exit(-1);
    }

}


void build_graph(vector<pc_rule*> &pc, vector<struct node *> &G)
{
    G.resize(pc.size());
    for(size_t i = 0; i < pc.size(); i++) {
        pc[i]->n = new node();
        
        G[i] = pc[i]->n;
        //index
        G[i]->index = i;
        G[i]->r = pc[i];
    }

    //G[pc.size()-1].index = pc.size() - 1;

    ssize_t i = (ssize_t)(pc.size()-2);
    ssize_t j;

    for(;i>=0;i--){
        //G[i]->index = i;
        for(j=i+1;j<(ssize_t)pc.size();j++) {
            if(overlap(pc[i], pc[j])) {
                G[i]->out.push_back(G[j]);
                G[j]->in.push_back(G[i]);
            } 
        }
    }
}



void build_graph_fast(vector<pc_rule*> &pc, vector<struct node*> &G)
{
    G.resize(pc.size());
    for(size_t i = 0; i < pc.size(); i++) {
        pc[i]->n = new node();
        
        G[i] = pc[i]->n;
        G[i]->index = i;
        G[i]->r = pc[i];
    }

    struct overlap_trees ots;

    rule_boundary rb;
    init_boundary(rb);

    init_overlap_trees(&ots, rb, pc); 

    int i = pc.size() -2;
    for(; i>=0; i--) {
        vector<int> set;
        //vector<int> set2;
        find_overlap_rules(&ots, pc, pc[i], i+1, pc.size(), set);
        sort(set.begin(), set.end());
        //find_overlap_rules_slow(pc, pc[i], i+1, pc.size(), set);

        //if(set1 != set2) {
        //    cout<<"error"<< i <<endl;
        //}
        for(auto j = set.begin(); j != set.end(); j++) {
            G[i]->out.push_back(G[*j]);
            G[*j]->in.push_back(G[i]);
        }

    }


}

int dfs(struct node *n)
{
    if(n->move == -1) {
        if(n->out.size() != 0)
            n->move = dfs(n->out[0]) + 1;
        else 
            n->move = 1;
    }
    return n->move;
}



void compute_cost(struct node *n) 
{
    int contribute = 0;
    int remove = 0;
    //int node_cnt = 0;
    for(size_t i = 0; i < n->in.size(); i++) {
        auto parent = n->in[i];
        struct node * candidate;
        int gain = 0;

        if(parent->out.size() > 1) {
            candidate = (parent->out[1]);
            gain = candidate->move +1;
        }
        else {
            gain = 1;
        }
        
        if(parent->out[0] == n) {
            contribute += (n->move + 1);
            remove +=  gain;
            //node_cnt ++;
        }
    }
    
    n->cost = contribute - remove;
    //if(node_cnt != 0)
    //    n->cost = (double)(contribute - remove)/node_cnt;
    //else
    //    n->cost = 0;

}

//we have pc of rules, tcam of rules, and Graph index of rules.
//the Graph index of rules and tcam of rules should be the same, 
//because this will reduce the operation of Graph rule movement.
//however, after insert a rule into pc, the index of rules stayed below the 
//inserted rule will all be changed ...
//so we would better have a consistency index system of rules. 
//
//Using the rules on the tcam will be easy
//so remember that the rule indexing in pc does not equal to the postion of rules 
//stored in TCAM. 
//
//


void del_rule(tcam & tmap, vector<pc_rule*> &pc, pc_rule *r)
{
    struct node *n  = r->n;

    for(size_t i = 0; i < n->in.size(); i++) {
        auto parent = n->in[i];
        if(parent->out[0] == n) {
            parent->out.erase(parent->out.begin());
            parent->move = -1;
            dfs(parent);
            //compute_cost(parent);
        }
        else {
            auto tail = remove_if(parent->out.begin(), parent->out.end(), [n](struct node *outn){return outn == n;});
            parent->out.erase(tail, parent->out.end());
        }
    }

    for(size_t i = 0; i < n->out.size(); i++) {
        auto child = n->out[i];
        auto tail = remove_if(child->in.begin(), child->in.end(), [n](struct node *inn) {return inn == n;});
        child->in.erase(tail, child->in.end());
        
    }
    
    tmap.del(n->index);
    
    n->in.clear();
    n->out.clear();
    n->index = -1;
    n->move = -1;
    n->cost = -1;
    n->r = nullptr; 
    n->valid = false;

    //pc.erase(pc.begin() + pos);
    auto tail =  remove_if(pc.begin(), pc.end(), [r](const struct pc_rule *rule){return  r == rule;});
    pc.erase(tail, pc.end());
}

void updfs(struct node *n) 
{
    for(size_t i = 0; i < n->in.size(); i++) {
        struct node *p = n->in[i];
        if(p->out[0] == n) {
            updfs(p);
        }
    }

    dfs(n);
}

void udfs(struct node *n) 
{
    n->move = -1;
    for(size_t i = 0; i < n->in.size(); i++) {
       struct node *p = n->in[i];
       if(p->out[0] == n) {
           udfs(p);
       }
    }

    dfs(n);
    //compute_cost(n);
}

void graph_del(pc_rule *r)
{
    struct node *n  = r->n;

    for(size_t i = 0; i < n->in.size(); i++) {
        auto parent = n->in[i];
        if(parent->out[0] == n) {
            parent->out.erase(parent->out.begin());
            parent->move = -1;
            dfs(parent);
            //compute_cost(parent);
        }
        else {
            auto tail = remove_if(parent->out.begin(), parent->out.end(), [n](struct node *outn){return outn == n;});
            parent->out.erase(tail, parent->out.end());
        }
    }   

    for(size_t i = 0; i < n->out.size(); i++) {
        auto child = n->out[i];
        auto tail = remove_if(child->in.begin(), child->in.end(), [n](struct node *inn) {return inn == n;});
        child->in.erase(tail, child->in.end());        
    }

    n->in.clear();
    n->out.clear();
    n->move = -1;
    n->cost = -1;
    //n->valid = false;
}

bool cmp(const node* n1, const node* n2)
{
    return n1->index < n2->index;
}

void Graph_add(tcam & tmap, vector<pc_rule*> &pc_graph, vector<struct node *> &G , int invaild_pos)
{
    int interval_start;
    if(invaild_pos != -1) {
        interval_start = invaild_pos + 1; 

        for(int i = 0; i < (int)pc_graph.size(); i ++)  
        {
           for(auto parent = pc_graph[i]->n->in.begin(); parent != pc_graph[i]->n->in.end(); parent ++)  
           { 
            if(tmap.valid[(*parent)->index] == true) {
                if((*parent)->out[0] == pc_graph[i]->n) {
                    (*parent)->out.erase((*parent)->out.begin());
                    (*parent)->move = -1;
                    dfs((*parent));
                }
                else {
                    struct node *n= pc_graph[i]->n;
                    auto tail = remove_if((*parent)->out.begin(), (*parent)->out.end(), [n](struct node *outn){return outn == n;});
                    (*parent)->out.erase(tail, (*parent)->out.end());
                }
                auto insert_pos = lower_bound((*parent)->out.begin(), (*parent)->out.end(), pc_graph[i]->n, cmp);
                (*parent)->out.insert(insert_pos, pc_graph[i]->n);
            }
            }
            for(int j = interval_start; j < pc_graph[i]->n->index; j ++) {
                if(overlap(&tmap.pc[j], pc_graph[i]) && tmap.pc[j].ruleid < pc_graph[i]->ruleid && tmap.valid[j] == true) {
                    pc_graph[i]->n->in.push_back(tmap.pc[j].n);
                    auto insert_pos = lower_bound(tmap.pc[j].n->out.begin(), tmap.pc[j].n->out.end(), pc_graph[i]->n, cmp);
                    tmap.pc[j].n->out.insert(insert_pos, pc_graph[i]->n);
                }
           }
           interval_start = pc_graph[i]->n->index + 1;
        }
    }
    else {
        interval_start = pc_graph[0]->n->index + 1;
        
        if(pc_graph.size() > 1) {

        for(int i = 0; i < (int)pc_graph.size(); i ++)  
        {
           for(auto parent = pc_graph[i]->n->in.begin(); parent != pc_graph[i]->n->in.end(); parent ++)  
           { 
            if(tmap.valid[(*parent)->index] == true) {
                if((*parent)->out[0] == pc_graph[i]->n) {
                    (*parent)->out.erase((*parent)->out.begin());
                    (*parent)->move = -1;
                    dfs((*parent));
                    //compute_cost((*parent));
                }
                else {
                    struct node *n= pc_graph[i]->n;
                    auto tail = remove_if((*parent)->out.begin(), (*parent)->out.end(), [n](struct node *outn){return outn == n;});
                    (*parent)->out.erase(tail, (*parent)->out.end());
                }
                auto insert_pos = lower_bound((*parent)->out.begin(), (*parent)->out.end(), pc_graph[i]->n, cmp);
                (*parent)->out.insert(insert_pos, pc_graph[i]->n);
            }
            }
            for(int j = interval_start; j < pc_graph[i]->n->index; j ++) {
                if(overlap(&tmap.pc[j], pc_graph[i]) && tmap.pc[j].ruleid < pc_graph[i]->ruleid && tmap.valid[j] == true) {
                    pc_graph[i]->n->in.push_back(tmap.pc[j].n);
                    auto insert_pos = lower_bound(tmap.pc[j].n->out.begin(), tmap.pc[j].n->out.end(), pc_graph[i]->n, cmp);
                    tmap.pc[j].n->out.insert(insert_pos, pc_graph[i]->n);
                }
           }
           interval_start = pc_graph[i]->n->index + 1;
        }
        }
        
        auto r = pc_graph.begin();

        int pos = (*r)->n->index;
        list<struct node *> update;
    
        for(int i = 0; i < pos; i++) {
        if(overlap(&tmap.pc[i], (*r)) && tmap.valid[i] == true && tmap.pc[i].ruleid < (*r)->ruleid) {
            //if(overlap(&tmap.pc[i], r)) {
                 //cout<<"find overlap in tcam "<<i<<endl;
            if(G[i]->index >= (*r)->n->index){
                cout<<"out cycle "<<G[i]->index<<" >= "<<(*r)->n->index<<endl;
                getchar();
                continue;
            }
            /*for(size_t j = 0; j < G[i]->out.size(); j++) {
                if(pos < (G[i]->out[j])->index ) {
                    G[i]->out.insert(G[i]->out.begin() + j, (*r).n);
                    break;
                }
            }*/
            auto insert_pos = lower_bound(G[i]->out.begin(), G[i]->out.end(), (*r)->n, cmp);
            G[i]->out.insert(insert_pos, (*r)->n);

            if(G[i]->out.size() == 0 || G[i]->out[G[i]->out.size()-1]->index < pos) {
                G[i]->out.push_back((*r)->n);
            }

            (*r)->n->in.push_back(G[i]);

            if(G[i]->out[0] == (*r)->n) {              
                update.push_back(G[i]);
            }
        }
        }

        for(size_t i = pos+1; i < G.size(); i++) {
            if(overlap(&tmap.pc[i], (*r)) && tmap.valid[i] == true && tmap.pc[i].ruleid > (*r)->ruleid) {
                if(G[i]->index <= (*r)->n->index){
                    cout<<"in cycle "<<G[i]->index<<" <= "<<(*r)->n->index<<endl;
                    getchar();
                    continue;
                }
                (*r)->n->out.push_back(G[i]);
                G[i]->in.push_back((*r)->n);
            }
        }

        for(auto n = update.begin(); n != update.end(); n ++) {
            udfs(*n);
        }

        //if no nodes need to update, we should at least update the r.n
        if(update.size() == 0) {
            dfs((*r)->n);
        }
    }
}

bool lazy_check(tcam &tmap, pc_rule *r, int pos, vector<pc_rule*> &pc, vector<struct node*> &G)
{
    //up node
    struct node * un = NULL;
    //down node
    struct node * dn = NULL;
    int upos = 0;
    int dpos = 0;

    for(int i = pos-1; i >=0 ; i--) {
        //cout<<"lazy_check: check overlap before pos"<<endl;
        if(overlap(pc[i], r)) {
            un = pc[i]->n;
            break;
        }
    }

    if(un == NULL) {
        upos = pos;
    }
    else{
        upos = un->index +1;
    }
    
    for(size_t i = pos; i < pc.size(); i++) {
        //cout<<"lazy_check: check overlap after pos"<<endl;
        if(overlap(pc[i], r)){
            dn = pc[i]->n;
            break;
        }
    }

    if(dn == NULL) {
        dpos = pc.size(); 
    }
    else {
        dpos = dn->index;
    }

    for(int i = upos; i < dpos; i++) {
        //cout<<"upos = "<<upos<<" dpos = "<<dpos<<"pc_size = "<<pc.size()<<endl;
        if(tmap.valid[i] == false) {
            cout<<"insert in the empty entry"<<endl;
            //update rule
            r->priority = i;
            //update tcam
            tmap.pc[i] = *r;
            tmap.valid[i] = true;
            //update pc, may repeat
            pc.insert(pc.begin() + pos, r);
            //update graph
            G[i] = r->n;
            G[i]->r = pc[pos];
            //Graph_add(tmap, r, G, i);

            return true;
        }
    }
    return false;
}

int swap_insert(tcam & tmap, vector<struct node*> &G, pc_rule *r, int insert_pos, int up_bound, int invaild_pos) 
{
    struct node *sn = NULL;
    struct node *wn = r->n;
    pc_rule sr; 
    pc_rule wr = *r;

    //vector <pc_rule *> pc_graph;
    vector <pc_rule*> pc_graph;

    int empty_flag = 0;

    while(insert_pos != -1) {
        if(up_bound != -1) {
            for(int i = up_bound; i < insert_pos; i ++) {
                if(tmap.valid[i] == false) {
                    //cout<<"find empty after up_bound"<<endl;
                    //update rule
                    wr.priority = i;
                    wn->index = i;
                    //update Graph
                    G[i] = wn;
                    G[i]->index = i;
                    //update tcam
                    tmap.pc[i] = wr; 
                    tmap.valid[i] = true;      
                    pc_graph.push_back(&(tmap.pc[i]));

                    empty_flag = 1;
                    break;
                }
            }
        }

        if(empty_flag == 1) {
            break;
        }
        else {
            //cout<<"swap insert"<<endl;
            sn = G[insert_pos];
            sr = tmap.pc[insert_pos];

            //update rule
            wr.priority = insert_pos;
            wn->index = insert_pos;
            //update Graph
            G[insert_pos] = wn;
            G[insert_pos]->index = insert_pos;
            //update tcam
            tmap.pc[insert_pos] = wr;       

            //pc_graph.push_back(&tmap.pc[insert_pos]);
            pc_graph.push_back(&(tmap.pc[insert_pos]));

            wn = sn;
            wr = sr;
            up_bound = insert_pos + 1;

            if(sn->out.size() > 0) { //&& tmap.valid[sn->out[0]->index] == true) {
                insert_pos = sn->out[0]->index;
            }
            else {
                break;
            }
        }
    }
    //empty check
    if(empty_flag == 0) {
        for(int i = up_bound; i < tmap.pc.size(); i++) {
            if(tmap.valid[i] == false) {
                //cout<<"find empty before the end"<<endl;
                //update rule
                wr.priority = i;
                wn->index = i;
                //update Graph
                G[i] = wn;
                G[i]->index = i;
                //update tcam
                tmap.pc[i] = wr; 
                tmap.valid[i] = true;      
                pc_graph.push_back(&(tmap.pc[i]));

                empty_flag = 1;
                break;
            }
        }
        if(empty_flag == 0){
            //cout<<"insert in the end"<<endl;
            //update Graph
            G.push_back(wn);
            G[G.size() -1]->index = G.size() - 1;
            //update rule
            wr.priority = G.size() - 1;
            wn->index = G.size() - 1;
            //update tcam
            tmap.pc.push_back(wr);
            tmap.valid.push_back(true);

            //pc_graph.push_back(&tmap.pc[G.size() - 1]);
            pc_graph.push_back(&(tmap.pc[G.size() - 1]));
        }
    }

    Graph_add(tmap, pc_graph, G, invaild_pos);

    return r->n->index;
}

void tcam_insert(tcam & tmap, vector<struct node*> &G, pc_rule *r, int pos, vector<pc_rule*> &pc) 
{
    int insert_pos = -1;
    for(size_t i = pos; i< pc.size(); i++) {
        if(overlap(pc[i], r)) {
            if(insert_pos == -1) {
                insert_pos = pc[i]->n->index;
            }
            else if(pc[i]->n->index < insert_pos) {
                insert_pos = pc[i]->n->index;
            }
        }
    }
    
    int max_overlap_pos = -1;
    for(int i = 0; i < pos; i++){
        if(overlap(pc[i], r) && pc[i]->n->index > insert_pos) {
            //cout<<endl<<"pos "<<pos<<endl;
            //cout<<"reorder old insert pos "<<insert_pos<<" overlap index "<<pc[i]->n->index<<" old ruleid "<<tmap.pc[insert_pos].ruleid<<" new ruleid "<<pc[i]->ruleid<<" r ruleid "<<r->ruleid<<endl;
            if(pc[i]->n->index > max_overlap_pos) {
                max_overlap_pos = pc[i]->n->index;
            }
        }        
    }
    if(max_overlap_pos != -1) {
        //cout<<"insert_pos = "<<insert_pos<<" max_overlap_pos = "<<max_overlap_pos<<endl;
    }

    if(insert_pos == -1) {
        //cout<<"insert_pos "<<insert_pos<<" max_overlap_pos "<<max_overlap_pos<<endl;
        swap_insert(tmap, G, r, insert_pos, max_overlap_pos + 1, -1);
    }
    else if(insert_pos < max_overlap_pos) {
        swap_insert(tmap, G, r, insert_pos, -1, -1);
        while(insert_pos < max_overlap_pos) {
            tmap.valid[insert_pos] = false;
            insert_pos = swap_insert(tmap, G, r, G[insert_pos]->out[0]->index, insert_pos + 1, insert_pos); 
        } 
    }
    else if(insert_pos > max_overlap_pos) {
        //cout<<"insert_pos "<<insert_pos<<" max_overlap_pos "<<max_overlap_pos<<endl;
        swap_insert(tmap, G, r, insert_pos, -1, -1);
    }        
}


void add_rule(tcam & tmap, vector<struct node*> &G, pc_rule *r, int pos, vector<pc_rule*> &pc)
{
    //first add the node in the tcam, perform movement
    //then update the graph.
    
    r->n = new node();
    //first lazy check
    
    //if(!lazy_check(tmap, r, pos, pc, G)) {
        //cout<<"before tcam_insert ruleid "<<r->ruleid<<endl;
        tcam_insert(tmap, G, r, pos, pc); 
        //cout<<"after tcam_insert ruleid "<<r->ruleid<<endl;

    //}

    pc.insert(pc.begin() + pos, r);

    /*cout<<endl;
    int count = 0;
    for(size_t i = 0; i< pc.size(); i++) {
        //cout<< pc[i]->n->index <<" ";
        /*if(pc[i]->ruleid == 816){
            cout<<"count "<<count<<endl;
        }*/
        /*cout<< pc[i]->ruleid <<" ";
        count ++;
    }
    cout<<endl;*/
    //??
    r->n->r = pc[pos];
    
}



void measure(vector<struct node*> &G) 
{
    vector<struct node*> root;

    for(size_t i = 0; i< G.size() ;  i++) {
        if(G[i]->out.size() == 0 && G[i]->valid) {
            G[i]->move = 1;
        }
        if(G[i]->in.size() == 0 && G[i]->valid) {
            root.push_back(G[i]);
        }
    }

    for(size_t i = 0; i < root.size(); i++) {
        dfs(root[i]);
    }
    for(size_t i = 0; i< G.size(); i++) {
        if(G[i]->move == -1 && G[i]->valid)
            updfs(G[i]);
    }
    
    //compute avearage move
    int sum = 0;
    int max_move = 0;
    int valid_cnt = 0;
    for(size_t i = 0; i < G.size(); i++) {
        if(G[i]->valid == true) {
            sum += G[i]->move;
            if(G[i]->move > max_move)
                max_move = G[i]->move;
            valid_cnt ++;
        }
    }
    cout <<"max move "<<max_move<<endl;
    cout <<"average move "<<(double)sum/valid_cnt<<endl;


    //compute cost
    //double max_cost = 0;
    //double min_cost = 0xffffff;
    //for(size_t i = 0; i < G.size(); i++) {
    //    //compute_cost(G[i]);
    //    if(G[i]->cost > max_cost) {
    //        max_cost = G[i]->cost;
    //    }
    //    if(G[i]->cost < min_cost) {
    //        min_cost = G[i]->cost;
    //    }
    //    //cout<<"cost "<<G[i].cost<<endl;
    //}


    //cout<<"cost "<<max_cost<<" "<<min_cost<<endl;
}

void measure_avgmove(vector<struct node*> &G, double &avg_move_r, int &max_move_r) 
{
    int sum =0;
    int max_move = 0;
    int valid_cnt = 0;

    for(size_t i = 0; i< G.size(); i++) {
        if(G[i]->valid == true) {
            sum += G[i]->move;
            valid_cnt ++;
            if(G[i]->move > max_move) 
                max_move = G[i]->move;
        }
    }
    cout<<"sum move "<<sum<<endl;
    cout<<"max move "<<max_move<<endl;

    avg_move_r = (double)sum/valid_cnt;
    max_move_r = max_move;
    cout<<"avg move "<<avg_move_r<<endl;
}

//vector<vector<pc_rule*> > split_ruleset_bound(vector<struct node*> &G, vector<pc_rule*> &pc, int bound, tcam & tmap) 
//{
//    vector<struct node*> root; 
//    for(size_t i = 0; i < G.size(); i++) {
//        if(G[i]->in.size() == 0) {
//            G[i]->type = ROOT;
//            root.push_back(G[i]);
//        }
//    }
//
//    for(size_t i = 0; i< root.size(); i++) {
//
//    }
//
//}

//vector<vector<pc_rule*> >  split_ruleset_cost(vector<struct node*> &G, vector<pc_rule*> &pc, int num_threshold, tcam & tmap)
//{
//    vector<struct node*> SG = G;
//    vector<vector<pc_rule*> > ret;
//
//    vector<pc_rule*> sub;
//    sort(SG.begin(), SG.end(), [](struct node *n1, struct node* n2) { return n1->cost > n2->cost;});
//    
//    int i = 0;
//    for(; i < num_threshold; i++) {
//        struct node * nd = SG[0];
//        SG.erase(SG.begin());
//        sub.push_back(nd->r);
//        del_rule(tmap, pc, nd->r);
//        sort(SG.begin(), SG.end(), [](struct node *n1, struct node* n2) { return n1->cost > n2->cost;});
//    }
//
//    sort(sub.begin(), sub.end(), [](pc_rule *r1, pc_rule *r2) {return r1->priority <  r2->priority;});
//    //double avg_move_r;
//    //int max_move_r;
//    //measure_avgmove(G, avg_move_r, max_move_r);
//
//
//    cout<<"sub ruleset "<<sub.size()<<endl;
//    vector<struct node *> G2;
//    build_graph(sub,G2); 
//    measure(G2);
//    return vector<vector<pc_rule*> >();
//    
//}


vector<pc_rule*> remove_redund_pkg(vector<pc_rule> &pc, vector<pc_rule> &pc_prefix) 
{
    rule_boundary rb;
    init_boundary(rb);
    
    //pc_prefix = pc;
    
    vector<pc_rule*> pcr;
    for_each(pc.begin(), pc.end(), [&pcr](pc_rule &r){pcr.push_back(&r);});
    cout<<"Orignal ruleset "<<pc.size()<<endl;

    remove_redund_rt(pcr, rb, false);
    //return move(pcr);
    extend_rules(pcr, pc_prefix);
    
    cout<<"After prefix explanation: "<< pc_prefix.size()<<endl;

    vector<pc_rule*> pc_pprefix;
    for_each(pc_prefix.begin(), pc_prefix.end(), [&pc_pprefix](pc_rule &r){pc_pprefix.push_back(&r);});


    
    for(int i = 0; i< (int)pcr.size();i ++) {
        pc_pprefix[i]->priority = i;
    }

    return move(pc_pprefix);

}



void print_path(struct node *n) 
{
    struct node *curr = n;
    while(curr->out.size() != 0) {
        cout<<curr->r->priority<< " - ";
        curr = curr->out[0];
    }

    cout<<endl;
}

void print_long_paths(vector<struct node*> &G, int length_path) 
{
    for(size_t i = 0; i < G.size(); i++) {
        if(G[i]->valid && G[i]->move > length_path) {
            print_path(G[i]);
        }
    }
}

//bool label_flags = false;

void print_out_degree(vector<struct node *> &G)
{
    for(size_t i = 0; i< G.size(); i++) {
        if(G[i]->valid) {
            cout<<G[i]->out.size()<<endl;
        }
    }
}

int label_node(struct node *n) 
{
    if(n->in.size() == 0) {
        n->label = 0;
    }

    if(n->label == -1) {
        int zcnt = 0;
        int ocnt = 0;
        for(size_t i = 0; i < n->in.size(); i++) {
            if(label_node(n->in[i]) == 0) {
                zcnt ++;
            }

            if(label_node(n->in[i]) == 1) {
                ocnt ++;
            }
        }
        n->label = (zcnt > ocnt) ? 1 : 0;
        //if (zcnt > ocnt) {
        //    n->label = 1;
        //}
        //else if (zcnt < ocnt) {
        //    n->label = 0;
        //}
        //else {
        //    n->label = (label_flags == true) ? 1: 0;
        //    label_flags = !label_flags;
        //}
    }
    return n->label;
}

void zo_split(vector<struct node *> &G, vector<struct node *> &subG1, vector<struct node *> &subG2)
{
    cout << endl;
    cout << "0 - 1 split" <<endl;
    for(size_t i = 0; i< G.size(); i++) {
        label_node(G[i]);
    }

    vector<pc_rule*> pc_0;
    vector<pc_rule*> pc_1;

    for(size_t i = 0; i < G.size(); i++) {
        if(G[i]->label == 0) {
            pc_0.push_back(G[i]->r);
        }
        if(G[i]->label == 1) {
            pc_1.push_back(G[i]->r);
        }
    }

    vector<struct node *> G0;
    vector<struct node *> G1;

    cout<<"split sets 0 " << pc_0.size()<<endl;
    build_graph(pc_0, G0);
    measure(G0);
    cout<<"split sets 1 " << pc_1.size()<<endl;
    build_graph(pc_1, G1);
    measure(G1);


    subG1 = move(G0);
    subG2 = move(G1);
}

void filtering_largerules(vector<pc_rule*> &in, vector<pc_rule*> &out, vector<pc_rule*> * filtering_out) 
{
    for(auto rule = in.begin(); rule != in.end(); rule++) {
        if(check_range_size((*rule)->field[0], 0) && check_range_size((*rule)->field[1], 1)) {
            if(filtering_out) {
                filtering_out->push_back(*rule);
            }
            continue;
        }
        out.push_back(*rule);
    }

}

void hybrid_arch(vector<pc_rule*> &in)
{
    cout<<endl<<"Hybrid arch" <<endl;
    vector<pc_rule*> filtering;
    vector<pc_rule*> out;
    filtering_largerules(in, out, &filtering);
    cout<<"filtering "<<filtering.size()<<endl;


    vector<struct node*> G;
    vector<struct node*> G2;

    cout<<endl<<"******before******"<<endl;
    build_graph(in, G);
    measure(G);

    cout<<endl<<"*****after******"<<endl;
    build_graph(out, G2);
    measure(G2);

    cout<<endl<<"Split "<<endl;
    vector<struct node *> subG1;
    vector<struct node *> subG2;

    vector<struct node *> subG3;
    vector<struct node *> subG4;

    zo_split(G2, subG1, subG2);
    zo_split(subG2, subG3, subG4);


}

void pure_arch(vector<pc_rule*> &in) 
{
    cout<<endl<<"Pure Arch"<<endl;
    vector<struct node *> G2;
    build_graph(in, G2);

    //cout<<"before 0-1 split, number of set "<<in.size()<<endl;
    //measure(G2);

    vector<struct node *> subG1;
    vector<struct node *> subG2;

    vector<struct node *> subG3;
    vector<struct node *> subG4;

    zo_split(G2, subG1, subG2);
    zo_split(subG2, subG3, subG4);

}

int linear_search_pc(vector <pc_rule*> &pc, unsigned *ft){
    //cout<<"pc priority"<<endl;
    int count = -1;
    for(auto r = pc.begin(); r != pc.end(); r++){
        //cout<<(*r)->priority<<" ";
        count ++;
        if(ft[0] <= (*r)->field[0].high && ft[0] >= (*r)->field[0].low &&
            ft[1] <= (*r)->field[1].high && ft[1] >= (*r)->field[1].low &&
            ft[2] <= (*r)->field[2].high && ft[2] >= (*r)->field[2].low &&
            ft[3] <= (*r)->field[3].high && ft[3] >= (*r)->field[3].low &&
            ft[4] <= (*r)->field[4].high && ft[4] >= (*r)->field[4].low) {
        return (*r)->ruleid;  
    }
    //cout<<endl;
    }
    return -1;
}

int linear_search_tcam(tcam &tmap, unsigned *ft){
    for(size_t i = 0; i < tmap.pc.size(); i++){
        if(tmap.valid[i] == true){
            if(ft[0] <= tmap.pc[i].field[0].high && ft[0] >= tmap.pc[i].field[0].low &&
            ft[1] <= tmap.pc[i].field[1].high && ft[1] >= tmap.pc[i].field[1].low &&
            ft[2] <= tmap.pc[i].field[2].high && ft[2] >= tmap.pc[i].field[2].low &&
            ft[3] <= tmap.pc[i].field[3].high && ft[3] >= tmap.pc[i].field[3].low &&
            ft[4] <= tmap.pc[i].field[4].high && ft[4] >= tmap.pc[i].field[4].low) {
                return tmap.pc[i].ruleid;
            }
        }
    }
    return -1;
}


//void tcamcheck(vector<pc_rule*> &pc, vector<pc_rule> &tcampc)
void tcamcheck(vector<pc_rule*> &pc, tcam &tmap)
{
    field_type header[MAXDIMENSIONS];
    int fid;
    int pc_ret, tcampc_ret;
    int trace_cnt = 0;
    int dismactch_cnt = 0;
    int nf_cnt = 0;
    if(fpt != NULL){
        cout<<"loading trace"<<endl;
        while(fscanf(fpt, "%u %u %d %d %d %d\n", &header[0], &header[1], &header[2], &header[3], &header[4], &fid) == 6){
            trace_cnt ++;
            pc_ret = linear_search_pc(pc, header);
            /*if(pc_ret != fid){
                cout<<"pc_ret != fid"<<endl;
                cout<<"pc_ret = "<<pc_ret<<endl;
                cout<<"fid = "<<fid<<endl;
            }*/
            //tcampc_ret = linear_search_tcam(tcampc, header);
            tcampc_ret = linear_search_tcam(tmap, header);
            if(pc_ret == tcampc_ret && pc_ret == -1){
                nf_cnt++;
                cout<<"rule not found nf count:"<<nf_cnt<<endl;
                cout<<"fid = "<<fid<<endl<<endl;
            }
            if(pc_ret != tcampc_ret){
                dismactch_cnt++;
                cout<<"pc_ret != tcampc_ret dismactch count:"<<dismactch_cnt<<endl;
                cout<<"pc_ret = "<<pc_ret<<endl;
                cout<<"tcampc_ret = "<<tcampc_ret<<endl;
                cout<<"fid = "<<fid<<endl;
                cout<<"trace_cnt = "<<trace_cnt<<endl<<endl;
            }
        }
    }
}


int main(int argc, char *argv[])
{

    vector<pc_rule> pc_orig;
    vector<pc_rule> pc_prefix;
    vector<pc_rule*> pc;

    parseargs(argc, argv);
    loadrules(fpr, pc_orig);

    vector<pc_rule*> pc_r = remove_redund_pkg(pc_orig, pc_prefix);
    cout<<"load "<<pc_r.size()<<" rules"<<endl;

    //pure_arch(pc_r); 
    //hybrid_arch(pc_r);

    tcam tmap;
    vector<struct node*> G;
    for(size_t i = 0; i < (pc_r.size()-1)/2; i++) {
        add_rule(tmap, G, pc_r[2*i+1], i, pc);
        /*if(pc_r[2*i+1]->ruleid == 3976){
            cout<<"first round"<<endl;
            getchar();
        }*/
    }
    for(size_t i = 0; i < pc_r.size()/2; i++) {
        add_rule(tmap, G, pc_r[2*i], 2*i, pc);
        /*if(pc_r[2*i]->ruleid == 3976){
            cout<<"second round"<<endl;
            getchar();
        }*/
    }
    add_rule(tmap, G, pc_r[pc_r.size()-1], pc_r.size()-1, pc);

    cout<<"checking"<<endl;
    tcamcheck(pc, tmap);
    cout<<"pc.size() = "<<pc.size()<<endl;
    cout<<"tmap.pc.size() = "<<tmap.pc.size()<<endl;
    cout<<"check over"<<endl;

    //double avg_move;
    //int max_move;
    //measure_avgmove(G, avg_move, max_move);
    //cout<<endl;
    //cout<<"avg move "<<avg_move<<endl;
    //cout<<"max move "<<max_move<<endl;
    //print_long_paths(G, 10);
    //print_out_degree(G);
    //measure(G);



    //split_ruleset_cost(G, pc, pc.size()/2, tmap);
    //for(size_t i = 0; i< G.size(); i++) {
    //measure_avgmove(G, avg_move, max_move);
    //    if(G[i]->in.size() == G2[i]->in.size() &&
    //            G[i]->out.size() == G2[i]->out.size()
    //            && G[i]->move == G2[i]->move) {

    //        continue;
    //    }i
    //    else {
    //        cout << "error"<< i <<endl;
    //    }
    //}

        //print_out_degree(subG3);
    //print_out_degree(subG2);
    //print_long_paths(subG3, 10);

    return 0;
}

