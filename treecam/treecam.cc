#include <queue>
#include <unistd.h>


#include "treecam.h"
#include <assert.h>

using namespace std;

void get_range_map(vector<pc_rule*> &pc, int d, map<range, int> &range_map)
{
    for(auto rule = pc.begin(); rule != pc.end(); rule++) {
        range_map[(*rule)->field[d]] ++;
    }
}

void get_range_map_md(vector<pc_rule*> &pc, vector<map <range, int> > &md_rm)
{
    for(int d = 0; d < MAXDIMENSIONS; d++) {
        map<range, int> range_map;
        get_range_map(pc, d, range_map);
        md_rm.push_back(range_map);
    }
}


int choose_dim(Node *n, vector< map<range, int> > &md_rm)
{
    get_range_map_md(n->pc, md_rm);

    int dim = -1;
    size_t max = 0;
    for(int i = 0; i < MAXDIMENSIONS; i++) {
        if(md_rm[i].size() >= max && (n->dimmask & (1<<i))) {
            max = md_rm[i].size();
            dim = i;
        }
    }

    n->dim = dim;
    return dim;
}

struct range_info {
    unsigned long long high;
    int cnt;

    bool operator < (const range_info &other) const {
        return high > other.high;
    }

    range_info(unsigned long long l, int rule_cnt) : high(l), cnt(rule_cnt) {

    }

};


void split(int dim, Node *n, map<range, int> &range_map, int bin)
{
    int64_t s = 0;
    int64_t e = -1;
    int rules = 0;
    size_t count = 0;
    vector<range> segs;
    vector<range_info> rheap;
    
    for(auto rpair = range_map.begin(); rpair != range_map.end();
            rpair++) {

        range r = rpair->first;
        if(rules + range_map[r] > bin && r.low != (uint64_t)s && rules!=0) {
            e = r.low -1;
            if(e < s) {
                cout<<"error, e is smaller than r"<<endl;
                exit(-1);
            }
            segs.push_back(range(s,e));
            s = e+1;

            while(rheap.size() != 0 && (int64_t)rheap[0].high < s) {
                pop_heap(rheap.begin(), rheap.end());
                rheap.pop_back();
            }
            rules = 0;
            for_each(rheap.begin(), rheap.end(), [&rules](range_info &ri){rules+= ri.cnt;});
            //if(rules > bin){ 
            //    cout<<" The number of rules are large"<<endl;
            //}
        }

        rules += rpair->second;
        
        while(rheap.size() != 0 && rheap[0].high < r.low) {
            pop_heap(rheap.begin(), rheap.end());
            rheap.pop_back();
        }

        rheap.push_back(range_info(r.high, rpair->second));
        push_heap(rheap.begin(), rheap.end());

        count ++;
        if(count  == range_map.size()) {
            e = r.high;
            //searching the heap, may exist some higher items
            for(auto ri = rheap.begin(); ri != rheap.end(); ri ++) {
                if((int64_t)ri->high > e) {
                    e = ri->high;
                }
            }
            segs.push_back(range(s,e));
        }

    }



    n->segs = move(segs);
    for(auto r = n->segs.begin(); r != n->segs.end(); r++) {
        Node *child = new Node();
        child->rb = n->rb;
        child->dimmask &= ~(1<<n->dim);
        child->rb.field[dim] = *r;
        child->depth = n->depth +1;
        for(auto rule = n->pc.begin(); rule != n->pc.end(); rule++) {
            if(range::overlap((*rule)->field[dim], *r)) {
                child->pc.push_back(*rule);
            }
        }
        
        //assert(child->pc.size() != 0);
        remove_redund(child->pc, child->rb);
        n->childnodes.push_back(child);
    }
    n->rulecnt = n->pc.size();
    n->pc.clear();

}


void init_root(Node *root, vector<pc_rule> &pc) 
{
    for_each(pc.begin(), pc.end(), [&root](pc_rule &rule) {root->pc.push_back(&rule);});
    remove_redund(root->pc, root->rb);
    root->depth = 1;
    root->rulecnt = root->pc.size();
}



struct dt_tree_info{
    int max_depth;

    dt_tree_info():max_depth(0)
    {

    }
};

void measure_dt_node(Node *n, dt_tree_info &dti)
{ 
    if(n->depth > dti.max_depth) {
        dti.max_depth = n->depth;
    }
}

void build_tree(Node *root, int bin) 
{
    dt_tree_info dti;

    if(root->pc.size() < (size_t)bin) {
        root->isLeaf = 1;
        return;
    }
    else{
        root->isLeaf = 0;
    }


    list<Node*> q ;
    q.push_back(root);

    while(!q.empty()){
        Node *n = q.front();

        measure_dt_node(n, dti);

        q.pop_front();

        vector< map<range, int > > md_rm;
        //if(n->no ==  3158) {
        //    cout<<"here"<<endl;
        //}
        int dim = choose_dim(n, md_rm);
        split(dim, n, md_rm[dim], bin);

        for(auto child = n->childnodes.begin(); child != n->childnodes.end();
                child ++) {
            if((*child)->pc.size() <= (size_t)bin) {
                (*child)->isLeaf = 1;
                (*child)->rulecnt = (*child)->pc.size();
                continue;
            }
            (*child)->isLeaf = 0;
            q.push_front(*child);
        }
    }


    cout<<"building tree with max depth "<<dti.max_depth<<endl;
    
}

int count_leaf(Node *n, int bin) 
{
    if(n->rulecnt <= bin) {
        return 1;
    }

    int count = 0;

    for(auto child = n->childnodes.begin(); 
            child != n->childnodes.end();
            child++) {
        count += count_leaf(*child, bin);
    }

    return count;
}

void count_coarse_rules(Node *n, int &total,int coarse_bin) 
{
    if(n->rulecnt <= coarse_bin) {
        total += n->rulecnt;
        return;
    }

    for(auto child = n->childnodes.begin(); 
            child != n->childnodes.end();
            child++) {
        count_coarse_rules(*child, total, coarse_bin);
    }

}

struct tree_info
{
    int total_rules;
    int coarse_total_rules;
    int root_rules;

    int total_leaf;

    int max_no_reb;
    int max_local_reb;
    int max_global_reb;

    int max_del_upop;


    int sum_upop;
    int fine_bin;

    vector<vector<Node*> >  upop_del_aux;

    tree_info(int pc_rules, int fine_bin): total_rules(0), coarse_total_rules(0), root_rules(0), 
                        total_leaf(0), max_no_reb(0), max_local_reb(0), max_global_reb(0),max_del_upop(0), sum_upop(0)

    {
        this->fine_bin = fine_bin;
        this->upop_del_aux.resize(pc_rules);
    }
};


int displaced_rules(vector<pc_rule*> &pc, range &bound, int dim) 
{
    size_t size = pc.size();
    int disprules = 0;
    unsigned long long max_s = 0;
    for(int i = 0; i < (int)size; i++) {
        range tmp = range_in_boundary_1D(pc[i]->field[dim], bound); 
        if(tmp.low > max_s) {
            max_s = pc[i]->field[dim].low;
        }
    }

    if(max_s == bound.low) {
        cout<<"[*] max_s == bound.low "<<max_s<<endl;
        return - 1;
    }

    max_s -= 1;

    for(int i = 0; i < (int)size; i++) {
        if(range::within((uint32_t)max_s, pc[i]->field[dim])){
            disprules ++;
        }
    }
    return disprules;
}

int measure_upop_insert(Node *n, Node *p, int index, tree_info &ti) 
{
    int upop = 0;
    if(n->pc.size() < (size_t)ti.fine_bin) {
        upop = n->pc.size()+1;
        if(ti.max_no_reb < upop) {
            ti.max_no_reb = upop;
        }
    }
    
    if(n->pc.size() == (size_t)ti.fine_bin) {
        int dim = p->dim;
        int disprules = displaced_rules(n->pc, p->segs[index], dim);
        if(disprules < ti.fine_bin && disprules != -1) {
            upop = ti.fine_bin * disprules + 1;
            if(ti.max_local_reb < upop) {
                ti.max_local_reb = upop;
            }
        }
        else {
            /*cout<<"[*]Local Rebalancing fails "<<"disprules == "<<disprules<<endl;
            cout<<"[*]splitting the dim "<<dim<<endl;
            for_each(n->pc.begin(), n->pc.end(), [](pc_rule *r)
                    {
                        cout<<r->field[0].low<<" "<<r->field[0].high<<" "
                        <<r->field[1].low<<" "<<r->field[1].high<<" "
                        <<r->field[2].low<<" "<<r->field[2].high<<" "
                        <<r->field[3].low<<" "<<r->field[3].high<<" "
                        <<r->field[4].low<<" "<<r->field[4].high<<" "
                        <<endl;
                    });*/

            build_tree(n, ti.fine_bin -1); 
            int disprules = 0;
            for(size_t i = 1; i < n->childnodes.size(); i++) {
                disprules += n->childnodes[i]->pc.size();
            }
            upop = disprules * ti.fine_bin + 1;
            if(ti.max_local_reb < upop) {
                ti.max_local_reb = upop;
            }
            cout<<"[*]upop "<<upop<<endl;

        }
    }

    ti.sum_upop += upop;

    return upop;
}

void measure_upop_delete(Node *n, tree_info &ti)
{
    for(size_t i = 0; i < n->pc.size(); i++) {
        (ti.upop_del_aux[n->pc[i]->priority]).push_back(n);
    }
}

void compute_upop_del(tree_info &ti)
{
    int max_upop = 0;
    for(size_t i = 0; i < ti.upop_del_aux.size(); i++) {
        int upop = 0;
        for(auto n = ti.upop_del_aux[i].begin();
                n != ti.upop_del_aux[i].end();
                n ++) {
            upop += (*n)->pc.size() -1;
        }
        if(upop > max_upop) {
            max_upop = upop;
        }
    }

    ti.max_del_upop = max_upop;
    cout<<"max del upop "<<max_upop<<endl;
    

}


void count_rules(Node *n, Node *p, int index, tree_info &ti) 
{
    if(n->pc.size()) {
        ti.total_rules += n->pc.size();
        measure_upop_insert(n,p, index,ti);
        measure_upop_delete(n,ti);
        return;
    }

    for(size_t i = 0; i < n->childnodes.size();
            i++) {
        count_rules(n->childnodes[i], n, i, ti);
    }
}

void count_root_rules(Node *n, int &total, int coarse_bin, vector<range> &rv) 
{
    if(n->rulecnt <= coarse_bin) {
        int cnt = 1;
        vector<range> tmp;
        for(auto r = rv.begin(); r != rv.end(); r++) {
            range2prefix(*r, tmp);
            cnt *= tmp.size();
            tmp.clear();
        }
        total += cnt;
        return;
    }

    for(size_t i = 0; i < n->childnodes.size(); i++) {
        rv.push_back(n->segs[i]);
        count_root_rules(n->childnodes[i], total, coarse_bin, rv);
        rv.pop_back();
    }
}


void measure_tree(Node *root, size_t coarse_bin, size_t fine_bin, tree_info &ti)
{
    int coarse_leaf = count_leaf(root, coarse_bin);
    cout<<"coarse tree has "<<coarse_leaf<<" leaf"<<endl;

    int fine_leaf = count_leaf(root, fine_bin);
    cout<<"fine tree has "<<fine_leaf<<" leaf"<<endl;

    count_rules(root, nullptr, -1,  ti);
    count_coarse_rules(root, ti.coarse_total_rules, coarse_bin);
    vector<range> tmp;
    count_root_rules(root, ti.root_rules, coarse_bin, tmp);  

    ti.total_leaf += fine_leaf;


}


Node *root[32];
vector<pc_rule> rs[32];
vector<pc_rule> extend_rs[32];
FILE *fpt;
FILE *fpr;
int ruleset_cnt = 0;


void parseargs(int argc, char *argv[]) {
    int	c;
    while ((c = getopt(argc, argv, "n:t:r:")) != -1) {
        switch (c) {
            case 'n':
                ruleset_cnt = atoi(optarg);
                break;
            case 't':
                fpt = fopen(optarg, "r");
                break;
            case 'r':
                fpr = fopen(optarg, "r");
                break;
            default:
                break;
        }
    }

    
    if(ruleset_cnt == 0){
        printf("can't open ruleset file\n");
        exit(-1);
    }

}

void build_fine_tree(Node *n, int fine_binth)
{

    if(n->isLeaf) {
        build_tree(n, fine_binth); 
        return;
    }

    for(auto child = n->childnodes.begin(); child != n->childnodes.end();
            child ++ ) {
        build_fine_tree(*child, fine_binth);
    }

}


int check_tree(Node *tree, field_type *ft)
{
    Node *n = tree;

    while(!n->pc.size()) {
        bool found = false;
        for(size_t i = 0; i < n->segs.size(); i++) {
            if(range::within(ft[n->dim], n->segs[i])) {
                n = n->childnodes[i];
                found = true;
                break;
            }
        }

        if(!found) {
            return -1;
        }
    }

    return linear_search(n->pc, ft);
}



int check_treecam(Node **tree, int cnt, field_type *ft)
{    
    int final = 0x7fffffff;
    for(int i = 0; i < cnt; i++) {
        int ret = check_tree(tree[i], ft);
        if(ret < final && ret != -1) {
            final = ret;
        } 
    }

    if(final == 0x7fffffff) {
        final = -1;
    }
    return final;
}


 


void check(Node **tree, int cnt, vector<pc_rule> &pc, FILE *fp)
{

    vector<pc_rule*> ppc;
    for_each(pc.begin(), pc.end(), [&ppc](pc_rule &rule){ppc.push_back(&rule);});

    field_type ft[MAXDIMENSIONS];
    int count = 0; 
    while(load_ft(fp, ft)) {
        count ++;
        //if(count == 4704) {
        //    cout<<"here"<<endl;
        //}
        int ret1 = linear_search(ppc, ft);
        int ret2 = check_treecam(tree, cnt, ft);
        if(ret1 != ret2) {
            cout<<"ret1 != ret2 "<<count<<endl;
        }
    }

}



int main(int argc, char *argv[]) 
{
    int rule_cnt = 0;
    int ex_rule_cnt = 0;
    int fine_bin = 8;
    int coarse_bin = 4 * 1024;

    parseargs(argc, argv);

    cout<<"loading "<<ruleset_cnt<<" trees"<<endl;

    for(int i = 0; i < ruleset_cnt; i++) {
        char filename[128];
        sprintf(filename, "f%d", i);
        //cout<<filename<<endl;
        FILE *tmp = fopen(filename, "r");
        if(!tmp) {
            perror("file not found");
            exit(-1);
        }

        //if(i == 2) {
        //    cout<<"ruleset "<<i<<endl;
        //}

        loadrules(tmp, rs[i], 1);
        extend_rules(rs[i], extend_rs[i]); 

        rule_cnt += rs[i].size();
        ex_rule_cnt += extend_rs[i].size();

        root[i] = new Node();
        init_root(root[i], extend_rs[i]);
        build_tree(root[i], coarse_bin);

        cout<<"Building fine tree"<<endl;
        build_fine_tree(root[i], fine_bin);

        //fineroot[i] = new Node();
        //init_root(fineroot[i], extend_rs[i]);
        //build_tree(fineroot[i], 8);
    }

    /*cout<<endl;
    cout<<"original "<<rule_cnt<<" rules"<<endl;
    cout<<"extend "<<ex_rule_cnt<<" rules"<<endl;
    cout<<"extend/original: "<<(double)ex_rule_cnt/rule_cnt<<endl;*/

    tree_info ti(rule_cnt, fine_bin);
    for(int i = 0; i < ruleset_cnt; i++) {
        measure_tree(root[i], coarse_bin, fine_bin, ti);
    }

    cout<<endl;
    cout<<"original "<<rule_cnt<<" rules"<<endl;
    cout<<"extend "<<ex_rule_cnt<<" rules"<<endl;
    cout<<"extend/original: "<<(double)ex_rule_cnt/rule_cnt<<endl;

    cout <<endl<<"In total"<<endl;
    cout <<"We have TreeCAM rules "<<ti.total_rules<<endl;
    cout <<"which is "<<ti.total_rules/(double)rule_cnt<<endl;
    
    cout <<"We have coarse tree "<<ti.coarse_total_rules<<" "<<ti.coarse_total_rules/(double)rule_cnt<<endl;
    cout <<"We have root rules "<<ti.root_rules<<" "<<endl;


    cout <<"max no reblancing "<< ti.max_no_reb << endl;
    cout <<"max local reblancing "<< ti.max_local_reb << endl;
    cout <<"average up op "<<(double)ti.sum_upop/ti.total_leaf<<endl;


    compute_upop_del(ti);

    if(fpt) {
        vector<pc_rule> pc;
        loadrules(fpr, pc, 2);
        check(root, ruleset_cnt, pc, fpt);
    }

    return 0;
}
