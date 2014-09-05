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
#include <sys/time.h>

#include "update.h"
#include "tcam.h"

using namespace std;

#include "rulesutils.h"
#include "rtrie.h"

FILE *fpr;
FILE *fpt;
size_t cost;
size_t max_cost;
size_t reorder_cnt;
size_t max_reorder_cnt;

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
        for(j = i+1;j < (ssize_t)pc.size() ; j ++) {
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
    //n->move = -1;
    for(size_t i = 0; i < n->in.size(); i++) {
        struct node *p = n->in[i];
        if(p->out.size() > 0) {
            if(p->out[0] == n) {
                p->move = n->move + 1;
                updfs(p);
            }
        }
    }
}

void Graph_del(pc_rule *r)
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

void Graph_adjust(tcam & tmap, vector<pc_rule> pc_graph, int invaild_pos)
{

    int interval_start;
    int start;
    if(invaild_pos != -1) {
        interval_start = invaild_pos; 
        start = 0;
    }
    else {
        interval_start = pc_graph[0].n->index;
        //interval_start = (*pc_graph.begin())->n->index;
        start = 1;
    }

    for(int i = start; i < (int)pc_graph.size(); i ++)  
    {
       for(auto parent = pc_graph[i].n->in.begin(); parent != pc_graph[i].n->in.end(); parent ++)  
       { 
        if(tmap.valid[(*parent)->index] == true ) {
            int ex_move = -1;
            if((*parent)->out.size() > 0) {
                ex_move = (*parent)->out[0]->move;
            }
            struct node *n= pc_graph[i].n;
            auto tail = remove_if((*parent)->out.begin(), (*parent)->out.end(), [n](struct node *outn){return outn == n;});
            if(tail != (*parent)->out.end()) {
                (*parent)->out.erase(tail, (*parent)->out.end());
                auto insert_pos = lower_bound((*parent)->out.begin(), (*parent)->out.end(), pc_graph[i].n, cmp);
                (*parent)->out.insert(insert_pos, pc_graph[i].n);

                if((*parent)->out[0]->move != ex_move) {
                    (*parent)->move = (*parent)->out[0]->move + 1;
                    updfs((*parent));
                }
            }
            else {
                if((*parent)->index > interval_start && (*parent)->index < pc_graph[i].n->index) {
                    if((*parent)->out.size() > 0) {
                        //ex_move = (*parent)->out[0]->move;
                        auto insert_pos = lower_bound((*parent)->out.begin(), (*parent)->out.end(), pc_graph[i].n, cmp);
                        (*parent)->out.insert(insert_pos, pc_graph[i].n);
                        if(insert_pos == (*parent)->out.begin()) {
                            if(pc_graph[i].n->move != ex_move) {
                                (*parent)->move = pc_graph[i].n->move + 1;
                                updfs((*parent));
                            }
                        }
                    }
                    else {
                        (*parent)->out.insert((*parent)->out.begin(), pc_graph[i].n);
                        (*parent)->move = pc_graph[i].n->move + 1;
                        updfs((*parent));
                    }
                }
            }
        }
    }           
    interval_start = pc_graph[i].n->index;
    }
}

void Graph_add(tcam &tmap, pc_rule *r, vector<pc_rule *> &out_rule, vector<pc_rule *> &in_rule, int insert_pos) 
{
    for(auto out = out_rule.begin(); out != out_rule.end(); out ++) {
        if(tmap.valid[(*out)->n->index] == false)
            getchar();
        (*out)->n->in.push_back(r->n);

        auto pos = lower_bound(r->n->out.begin(), r->n->out.end(), (*out)->n, cmp);
        r->n->out.insert(pos, (*out)->n);
    }

    dfs(r->n);

    for(auto in = in_rule.begin(); in != in_rule.end(); in ++) {
        if(tmap.valid[(*in)->n->index] == false)
            getchar();
        r->n->in.push_back((*in)->n);

        if((*in)->n->index < insert_pos) {
            if((*in)->n->out.size() > 0) {
                int ex_move = (*in)->n->out[0]->move;
                auto pos = lower_bound((*in)->n->out.begin(), (*in)->n->out.end(), r->n, cmp);
                (*in)->n->out.insert(pos, r->n);
                
                if(pos == (*in)->n->out.begin()) {
                    if(r->n->move != ex_move) {
                        (*in)->n->move = r->n->move + 1;
                        updfs((*in)->n);
                    }
                }
            }
            else {
                (*in)->n->out.insert((*in)->n->out.begin(), r->n);
                (*in)->n->move = r->n->move + 1;
                updfs((*in)->n);
            }
        }
    }
}

int swap_insert(tcam & tmap, vector<struct node*> &G, pc_rule *r, int insert_pos, int up_bound, int invaild_pos) 
{
    struct node *sn = NULL;
    struct node *wn = r->n;
    pc_rule sr; 
    pc_rule wr = *r;

    vector <pc_rule> pc_graph;

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
                    //wn->valid = true;
                    G[i] = wn;
                    G[i]->index = i;
                    //G[i]->valid = true;
                    //update tcam
                    tmap.pc[i] = wr; 
                    tmap.valid[i] = true;
                    //tmap.pc[i].n->valid = true; 
                    /*if(G[i]->valid == false) {
                        cout<<"can not change"<<endl;
                        getchar();
                    }*/

                    pc_graph.push_back(tmap.pc[i]);

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
            //G[insert_pos]->valid = true;
            //update tcam
            tmap.pc[insert_pos] = wr;       

            pc_graph.push_back(tmap.pc[insert_pos]);

            wn = sn;
            wr = sr;
            up_bound = insert_pos + 1; 

            if(sn->out.size() > 0) { 
                insert_pos = sn->out[0]->index;
            }
            else {
                break;
            }
        }
    }
    //empty check
    if(empty_flag == 0) {
        for(int i = up_bound; i < (int)tmap.pc.size(); i++) {
            if(tmap.valid[i] == false) {
                //cout<<"find empty before the end"<<endl;
                //update rule
                wr.priority = i;
                wn->index = i;
                //update Graph
                //wn->valid = true;
                G[i] = wn;
                G[i]->index = i;
                //G[i]->valid = true;
                //update tcam
                tmap.pc[i] = wr; 
                tmap.valid[i] = true; 
                //tmap.pc[i].n->valid = true;     
                pc_graph.push_back(tmap.pc[i]);
                /*if(G[i]->valid == false) {
                    cout<<"can not change"<<endl;
                    getchar();
                }*/
                empty_flag = 1;
                break;
            }
        }
        if(empty_flag == 0){
            //cout<<"insert in the end"<<endl;
            //update Graph
            G.push_back(wn);
            G[G.size() -1]->index = G.size() - 1;
            //G[G.size() -1]->valid = true;
            //update rule
            wr.priority = G.size() - 1;
            wn->index = G.size() - 1;
            //update tcam
            tmap.pc.push_back(wr);
            tmap.valid.push_back(true);

            pc_graph.push_back(tmap.pc[G.size() - 1]);
        }
    }

    cost += pc_graph.size() - 1;
    if(pc_graph.size() - 1 > max_cost) {
        max_cost = pc_graph.size() - 1;
    }

    if(pc_graph.size() > 1) {
        Graph_adjust(tmap, pc_graph, invaild_pos);
    }
    return r->n->index;
}

void tcam_insert(tcam & tmap, vector<struct node*> &G, pc_rule *r, int pos, vector<pc_rule*> &pc) 
{
    int insert_pos = -1;
    vector<pc_rule *> out_rule;
    for(size_t i = pos; i< pc.size(); i++) {
        if(overlap(pc[i], r)) {
            out_rule.push_back(pc[i]);
            if(G[pc[i]->n->index]->valid == false) {
                if(tmap.valid[pc[i]->n->index] == false){
                    cout<<"tcam item false"<<endl;
                }
                cout<<"out pc G false"<<endl;
                getchar();
            }
            if(insert_pos == -1) {
                insert_pos = pc[i]->n->index;
            }
            else if(pc[i]->n->index < insert_pos) {
                insert_pos = pc[i]->n->index;
            }
        }
    }
    
    int max_overlap_pos = -1;
    vector<pc_rule *> in_rule;
    for(int i = 0; i < pos; i++){
        if(overlap(pc[i], r)) {
            in_rule.push_back(pc[i]);
            if(G[pc[i]->n->index]->valid == false) {
                if(tmap.valid[pc[i]->n->index] == false){
                    cout<<"tcam item false"<<endl;
                }
                cout<<"in pc G false"<<endl;
                getchar();
            }
            if(pc[i]->n->index > insert_pos) {
                if(pc[i]->n->index > max_overlap_pos) {
                    max_overlap_pos = pc[i]->n->index;
                }
            }
        }        
    }

    if(insert_pos == -1) {
        //cout<<"insert_pos "<<insert_pos<<" max_overlap_pos "<<max_overlap_pos<<endl;
        insert_pos = swap_insert(tmap, G, r, insert_pos, max_overlap_pos + 1, -1);
        Graph_add(tmap, r, out_rule, in_rule, insert_pos);
    }
    else if(insert_pos < max_overlap_pos) {
        swap_insert(tmap, G, r, insert_pos, -1, -1);
        Graph_add(tmap, r, out_rule, in_rule, insert_pos);
        size_t tmp_reorder_cnt = 0; 
        int false_pos;
        while(insert_pos < max_overlap_pos) {
            ///cout<<"reorder insert"<<endl;
            int false_pos = insert_pos;
            /*tmap.valid[insert_pos] = false;
            G[insert_pos]->valid = false;*/
            insert_pos = swap_insert(tmap, G, r, G[insert_pos]->out[0]->index, insert_pos + 1, insert_pos); 
            tmap.valid[false_pos] = false;
            G[false_pos]->valid = false;
            //tmap.pc[false_pos].n->valid = false;
            tmp_reorder_cnt ++;
        } 
        reorder_cnt += tmp_reorder_cnt;
        if(tmp_reorder_cnt > max_reorder_cnt) {
            max_reorder_cnt = tmp_reorder_cnt;
        }
    }
    else if(insert_pos > max_overlap_pos) {
        //cout<<"insert_pos "<<insert_pos<<" max_overlap_pos "<<max_overlap_pos<<endl;
        swap_insert(tmap, G, r, insert_pos, -1, -1);
        Graph_add(tmap, r, out_rule, in_rule, insert_pos);
    } 
}


void add_rule(tcam & tmap, vector<struct node*> &G, pc_rule *r, int pos, vector<pc_rule*> &pc)
{
    //first add the node in the tcam, perform movement
    //then update the graph.
    
    r->n = new node();
    
    tcam_insert(tmap, G, r, pos, pc); 

    pc.insert(pc.begin() + pos, r);

    /*cout<<endl;
    int count = 0;
    for(size_t i = 0; i< pc.size(); i++) {
        cout<< pc[i]->n->index <<" ";
        cout<< pc[i]->ruleid <<" ";
        count ++;
    }
    cout<<endl;*/

    r->n->r = pc[pos];
    
}



double measure(vector<struct node*> &G, int measure_type) 
{
    if(measure_type == 0) {
        vector<struct node*> root;

        for(size_t i = 0; i< G.size() ;  i++) {
            /*if(G[i]->out.size() == 0 && G[i]->valid) {
                G[i]->move = 1;
            }*/
            if(G[i]->in.size() == 0 && G[i]->valid) {
                root.push_back(G[i]);
            }
        }

        for(size_t i = 0; i < root.size(); i++) {
            dfs(root[i]);
        }
        for(size_t i = 0; i < G.size(); i++) {
            if(G[i]->move == -1 && G[i]->valid) {
                dfs(G[i]);
                updfs(G[i]);
            }
        }
    }
    
    //compute avearage move
    int sum = 0;
    int max_move = 0;
    int valid_cnt = 0;
    int invaild_cnt = 0;
    for(size_t i = 0; i < G.size(); i++) {
        if(G[i]->valid == true) {
            if(G[i]->move < 1) {
                cout<<"error: G[i]->move = "<<G[i]->move<<" < 1"<<endl;
                getchar();
            }
            sum += G[i]->move;
            if(G[i]->move > max_move)
                max_move = G[i]->move;
            valid_cnt ++;
        }
        else {
            invaild_cnt ++;
        }
    }
    cout <<"max move "<<max_move<<endl;
    cout <<"average move "<<(double)sum/valid_cnt<<endl;
    cout<<"number of vaild nodes "<<valid_cnt<<endl;
    cout<<"number of invaild nodes "<<invaild_cnt<<endl;
    return (double)sum/valid_cnt;


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

void zo_split(vector<struct node *> &G, vector<struct node *> &subG1, vector<struct node *> &subG2, vector<pc_rule*> &pc_0, vector<pc_rule*> &pc_1)
{
    cout << endl;
    //cout << "0 - 1 split" <<endl;
    for(size_t i = 0; i< G.size(); i++) {
        label_node(G[i]);
    }

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

    //cout<<"split sets 0 " << pc_0.size()<<endl;
    build_graph(pc_0, G0);
    //measure(G0);
    //cout<<"split sets 1 " << pc_1.size()<<endl;
    build_graph(pc_1, G1);
    //measure(G1);


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

/*void hybrid_arch(vector<pc_rule*> &in)
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

}*/

int linear_search_pc(vector <pc_rule*> &pc, unsigned *ft){
    //cout<<"pc priority"<<endl;
    int count = -1;
    for(auto r = pc.begin(); r != pc.end(); r++){
        count ++;
        if(ft[0] <= (*r)->field[0].high && ft[0] >= (*r)->field[0].low &&
            ft[1] <= (*r)->field[1].high && ft[1] >= (*r)->field[1].low &&
            ft[2] <= (*r)->field[2].high && ft[2] >= (*r)->field[2].low &&
            ft[3] <= (*r)->field[3].high && ft[3] >= (*r)->field[3].low &&
            ft[4] <= (*r)->field[4].high && ft[4] >= (*r)->field[4].low) {
                return (*r)->ruleid;  
        }
    }
    return -1;
}

int linear_search_tcam(tcam &tmap, unsigned *ft) {
    for(size_t i = 0; i < tmap.pc.size(); i++) {
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

void tcamcheck(vector<pc_rule*> &pc, tcam &tmap1, tcam &tmap3, tcam &tmap4)
{
    field_type header[MAXDIMENSIONS];
    int junk;
    int fid;
    int pc_ret, tcampc_ret;
    int trace_cnt = 0;
    int dismactch_cnt = 0;  
    int tmp;
    //int nf_cnt = 0;
    if(fpt != NULL){
        cout<<"loading trace"<<endl;
        //while(fscanf(fpt, "%u %u %d %d %d %d\n", &header[0], &header[1], &header[2], &header[3], &header[4], &fid) == 6){
        while(fscanf(fpt, "%u %u %d %d %d %u %d\n", &header[0], &header[1], &header[2], &header[3], &header[4], &junk, &fid) == 7){
            trace_cnt ++;
            pc_ret = linear_search_pc(pc, header);
            /*if(pc_ret != fid){
                cout<<"pc_ret != fid"<<endl;
                cout<<"pc_ret = "<<pc_ret<<endl;
                cout<<"fid = "<<fid<<endl;
            }*/
            tcampc_ret = -1;
            tmp = linear_search_tcam(tmap1, header);
            if(tmp != -1) {
                tcampc_ret = tmp;
            }
            tmp = linear_search_tcam(tmap3, header);
            if(tmp != -1) {
                if(tcampc_ret == -1 || tcampc_ret > tmp) {
                    tcampc_ret = tmp;
                }
            }
            tmp = linear_search_tcam(tmap4, header);
            if(tmp != -1) {
                if(tcampc_ret == -1 || tcampc_ret > tmp) {
                    tcampc_ret = tmp;
                }
            }

            /*if(pc_ret == tcampc_ret && pc_ret == -1){
                nf_cnt++;
                cout<<"rule not found nf count:"<<nf_cnt<<endl;
                cout<<"fid = "<<fid<<endl<<endl;
            }*/

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

bool cmp_r(const pc_rule* p1, const pc_rule* p2)
{
    return p1->ruleid < p2->ruleid;
}

int main(int argc, char *argv[])
{

    vector<pc_rule> pc_orig;
    vector<pc_rule> pc_prefix;
    //vector<pc_rule*> pc;

    parseargs(argc, argv);
    loadrules(fpr, pc_orig);

    vector<pc_rule*> pc_r = remove_redund_pkg(pc_orig, pc_prefix);
    cout<<"load "<<pc_r.size()<<" rules"<<endl;

    //pure_arch(pc_r); 
    //hybrid_arch(pc_r);
    vector <pc_rule*> pc_old;
    vector <pc_rule*> pc_new;
    for(size_t i = 0; i < pc_r.size(); i ++) {
        if(i % 5 != 0) {
            pc_old.push_back(pc_r[i]);
        }
        else {
            pc_new.push_back(pc_r[i]);
        }
    }
    //cout<<pc_old.size()<<" "<<pc_new.size()<<" "<<pc_r.size()<<endl;

    double lower_cst_bnd = 0;
    double min_avg_cost = 0;
    vector<struct node*> G;
    build_graph(pc_old, G);

    vector<pc_rule*> pc_1;
    vector<pc_rule*> pc_2;
    vector<pc_rule*> pc_3;
    vector<pc_rule*> pc_4;
    vector<struct node *> subG1;
    vector<struct node *> subG2; 
    vector<struct node *> subG3;
    vector<struct node *> subG4;
    double cost1 = 0;
    double cost3 = 0;
    double cost4 = 0;

    cout<<endl<<"******before 0-1 split******"<<endl;
    if(measure(G, 0) > lower_cst_bnd) {        
         zo_split(G, subG1, subG2, pc_1, pc_2);
         cout<<"******split set 0******"<<endl;
         cout<<"size of set0 "<<pc_1.size()<<endl;
         cost1 = measure(subG1, 0);
         min_avg_cost = cost1;
         cout<<"******split set 1******"<<endl;
         cout<<"size of set1 "<<pc_2.size()<<endl;
         if(measure(subG2, 0) > lower_cst_bnd) {        
            zo_split(subG2, subG3, subG4, pc_3, pc_4);  
            cout<<"******split set 0******"<<endl;
            cout<<"size of set0 "<<pc_3.size()<<endl;
            cost3 = measure(subG3, 0);
            if(cost3 < min_avg_cost) {
                min_avg_cost = cost3;
            }
            cout<<"******split set 1******"<<endl;
            cout<<"size of set 1 "<<pc_4.size()<<endl; 
            cost4 = measure(subG4, 0);  
            if(cost4 < min_avg_cost) {
                min_avg_cost = cost4;
            }       
         }     
    }

    /*tcam tmap1;
    subG1.clear();
    vector<pc_rule*> pc_1_tmp;
    for(size_t i = 0; i < pc_1.size(); i ++) {
        add_rule(tmap1, subG1 ,pc_1[i], i, pc_1_tmp);
    }
    tcam tmap3;
    subG3.clear();
    vector<pc_rule*> pc_3_tmp;
    for(size_t i = 0; i < pc_3.size(); i ++) {
        add_rule(tmap3, subG3 , pc_3[i], i, pc_3_tmp);
    }
    tcam tmap4;
    subG4.clear();
    vector<pc_rule*> pc_4_tmp;
    for(size_t i = 0; i < pc_4.size(); i ++) {
        add_rule(tmap4, subG4 , pc_4[i], i, pc_4_tmp);
    }*/

    tcam tmap1;
    for(size_t i = 0; i < pc_1.size(); i ++) {
        tmap1.pc.push_back(*(pc_1[i]));
        tmap1.valid.push_back(true);
    }
    tcam tmap3;
    for(size_t i = 0; i < pc_3.size(); i ++) {
        tmap3.pc.push_back(*(pc_3[i]));
        tmap3.valid.push_back(true);
    }
    tcam tmap4;
    for(size_t i = 0; i < pc_4.size(); i ++) {
        tmap4.pc.push_back(*(pc_4[i]));
        tmap4.valid.push_back(true);
    }

    cout<<endl<<"begin to insert new rules"<<endl;
    cost = 0;
    max_cost = 0;
    reorder_cnt = 0;
    max_reorder_cnt = 0;
    int rules_cnt = 0;
    struct timeval start, end;

    size_t check_pos = 0;
    size_t check_interval = 100;
    //size_t check_interval = pc_new.size()/10;

    gettimeofday(&start, NULL);


    while(check_pos < pc_new.size()) {
        int insert_pos;
        if(min_avg_cost == cost4) {
            cout<<"inserting in the tcam4"<<endl;
            for(size_t i = check_pos; i < pc_new.size() && i < check_pos + check_interval; i ++) {
                auto pos = lower_bound(pc_4.begin(), pc_4.end(), pc_new[i], cmp_r);
                insert_pos = pos - pc_4.begin();
                add_rule(tmap4, subG4 , pc_new[i], insert_pos, pc_4);
                rules_cnt ++;
            }
            cout<<"******measure tcam4******"<<endl;
            cout<<"tcam size "<<tmap4.pc.size()<<endl;
            cost4 = measure(subG4, 1);
            min_avg_cost = cost4;
            if(min_avg_cost > cost1) {
                min_avg_cost = cost1;
            }
            if(min_avg_cost > cost3) {
                min_avg_cost = cost3;
            }
            cout<<endl;      
        }
        else if(min_avg_cost == cost3) {
            cout<<"inserting in the tcam3"<<endl;
            for(size_t i = check_pos; i < pc_new.size() && i < check_pos + check_interval; i ++) {
                auto pos = lower_bound(pc_3.begin(), pc_3.end(), pc_new[i], cmp_r);
                insert_pos = pos - pc_3.begin();
                add_rule(tmap3, subG3 , pc_new[i], insert_pos, pc_3);
                rules_cnt ++;
            }
            cout<<"******measure tcam3******"<<endl;
            cout<<"tcam size "<<tmap3.pc.size()<<endl;
            cost3 = measure(subG3, 1);
            min_avg_cost = cost3;
            if(min_avg_cost > cost1) {
                min_avg_cost = cost1;
            }
            if(min_avg_cost > cost4) {
                min_avg_cost = cost4;
            }
            cout<<endl;
        }
        /*else if(min_avg_cost == cost1) {
            for(size_t i = 0; i < pc_new.size(); i ++) {
                auto pos = lower_bound(pc_1_tmp.begin(), pc_1_tmp.end(), pc_new[i], cmp_r);
                insert_pos = pos - pc_1_tmp.begin();
                for(size_t j = 0; j < pc_1_tmp.size(); j ++) {
                    if(pc_new[i]->ruleid < pc_1_tmp[j]->ruleid) {
                        if(insert_pos != j) { 
                            cout<<insert_pos<<" "<<j<<endl;
                            getchar();
                        }
                        break;
                    }
                }
                add_rule(tmap1, subG1 , pc_new[i], insert_pos, pc_1_tmp);
                rules_cnt ++;
            }
        }*/
        else if(min_avg_cost == cost1) {
            cout<<"inserting in the tcam1"<<endl;
            for(size_t i = check_pos; i < pc_new.size() && i < check_pos + check_interval; i ++) {
                auto pos = lower_bound(pc_1.begin(), pc_1.end(), pc_new[i], cmp_r);
                insert_pos = pos - pc_1.begin();
                add_rule(tmap1, subG1 , pc_new[i], insert_pos, pc_1);
                rules_cnt ++;
            }
            cout<<"******measure tcam1******"<<endl;
            cout<<"tcam size "<<tmap1.pc.size()<<endl;
            cost1 = measure(subG1, 1);
            min_avg_cost = cost1;
            if(min_avg_cost > cost3) {
                min_avg_cost = cost3;
            }
            if(min_avg_cost > cost4) {
                min_avg_cost = cost4;
            }
            cout<<endl;
        }
        check_pos += check_interval;
    }


    gettimeofday(&end, NULL);

    double time_use = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout<<"finish inserting "<<rules_cnt<<" rules"<<endl;
    cout<<"use time "<<time_use<<" us"<<endl;
    cout<<"average time "<<time_use/rules_cnt<<" us for inserting one rule"<<endl;

    cout<<endl<<"total number of reorder "<<reorder_cnt<<endl;
    cout<<"average number of reorder "<<(double)reorder_cnt/rules_cnt<<endl;
    cout<<"max number of reorder "<<max_reorder_cnt<<endl;

    cout<<endl<<"total movement cost "<<cost<<endl;
    cout<<"average movement cost "<<(double)cost/rules_cnt<<endl;
    cout<<"max movement cost "<<max_cost<<endl;

    cout<<endl<<"******measure tcam1******"<<endl;
    cout<<"tcam size "<<tmap1.pc.size()<<endl;
    measure(subG1, 1);
    cout<<"******measure tcam3******"<<endl;
    cout<<"tcam size "<<tmap3.pc.size()<<endl;
    measure(subG3, 1);
    cout<<"******measure tcam4******"<<endl;
    cout<<"tcam size "<<tmap4.pc.size()<<endl;
    measure(subG4, 1);

    cout<<endl<<"checking"<<endl;
    //tcamcheck(pc, tmap);
    tcamcheck(pc_r, tmap1, tmap3, tmap4);
    cout<<"pc.size() = "<<pc_r.size()<<endl;
    cout<<"tmap.pc.size() = "<<tmap1.pc.size() + tmap3.pc.size() + tmap4.pc.size()<<endl;
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

