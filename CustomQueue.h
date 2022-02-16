#ifndef CustomQueue_h
#define CustomQueue_h

#include <iostream>
#include <vector>

class CustomQueue {
public:
    CustomQueue(long s) {
        size = s;
        xArr.resize(s);
        yArr.resize(s);
    }
    
    void push(long x, long y) {
        xArr[backIndex & size] = x;
        yArr[backIndex & size] = y;
        ++backIndex;
    }
    
    void pop(long &x, long &y) {
        x = xArr[frontIndex & size];
        y = yArr[frontIndex & size];
        ++frontIndex;
    }
    
    void clear() {
        frontIndex = backIndex = 0;
    }
    
    bool empty() {
        return frontIndex == backIndex;
    }
private:
    long size;
    std::vector<long> xArr, yArr;
    uint64_t frontIndex = 0, backIndex = 0;
};

#endif
