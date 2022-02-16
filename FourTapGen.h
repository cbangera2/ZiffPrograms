#ifndef FourTapGen_h
#define FourTapGen_h

class FourTap {
public:
    FourTap(long seed) {
        double a, ee = -1 + 1/2147483648.0;
        
        a = seed/2147483648.0;
        for (;nd <= M; nd++)
        {
            a *= 16807;
            a += ee * (long)(a);
            if (a >= 1) a += ee;
            ra[nd] = (long) (2147483648.0 * a);
        }
        nd = M;
        for(uint32_t i = 0; i<100001; i++) operator()();
        for (uint32_t i=0; i<10; ++i) printf("%10u%15ld\n", i, operator()());
    }
    
    long operator ()() {
        ++nd;
        ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M];
        return ra[nd&M];
    }
    
    constexpr static long min() {
        return 0;
    }
    
    constexpr static long max() {
        return 2147483648.0;
    }

private:
    constexpr static long M = 16383;
    long nd = 0, ra[M+1];
};

#endif
