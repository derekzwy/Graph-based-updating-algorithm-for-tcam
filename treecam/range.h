#ifndef __RANGE_H__
#define __RANGE_H__

struct range{
    unsigned long long low;
    unsigned long long high;


    range(unsigned long l, 
            unsigned long h): low(l),high(h) {
    }
    range() {
        low = 0;
        high = 0;
    }
    
    bool operator < (const range &r1) const {
        if(low != r1.low) {
            return low < r1.low;
        }
        else {
            return high < r1.high;
        }
    }

    static bool overlap(const range &r1, const range &r2) {
        return (!(r1.high < r2.low || r1.low > r2.high));
    }

    static bool within(uint32_t v, range &r) {
        return v <= r.high && v>= r.low;
    }
    
};

#endif

