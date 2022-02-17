#include "stdlib.h"
#include "stdio.h"
#include "math.h"

void randinit(long seed);

#define HEIGHT 16384 //power of two
#define WIDTH HEIGHT
#define RUNSMAX  2000000000
#define PRINTFREQ 1048575  //power of two minus 1
#define OUTFILE "perconebondNP30"
#define PROB 0.02301195 // z=8 0.25036834 
#define H (HEIGHT-1)
#define W (WIDTH-1)
#define R2 18 // square radius
#define S 32767  //for stack (queue) power of two minus 1
#define MAX 131072 //power of two
#define SEED 301
#define M  16383 //for random number generator -- power of two minus 1
#define GetFromStack2(XY) {XY = xylist[gptr & S]; ++gptr;}
#define PutOnStack2(XY)    {xylist[pptr & S]=XY; ++pptr;}
#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long    ra[M+1], nd;            //for random number generator
long    lat2[HEIGHT*WIDTH];     //for lattice
long    xylist[S+1];             //for stack (queue)

FILE *fp1;

int main(void)
{
    long     x, y, xo, yo, xp, yp, dir, prob, gptr, pptr, runs, i, j,k,m,n,dirmin;
    long    delx, dely, bin[64], xmax,z, temp, count[64],count2[64],gs,xl[64],trials,ntrials[64];
    double sum, pro, pr[64], rnd;
    long xy,xyo,xyp,dxy[1024];
    

//xmax is used to make sure the cluster does not grow to the edge of the lattice (at least in the x direction)

    xmax = HEIGHT/2;
    
//calculate the neighbor vectors

    int maxr = sqrt(R2)+1;  //max radius (plus 1)
    k=0;
    for(x=-maxr;x<=maxr;++x)
        for(y=-maxr;y<=maxr;++y)
            if(((x*x+y*y)<=R2)&&(x*x+y*y))
            {
                dxy[k]=x+y*WIDTH;
                printf("x,y %10ld%10ld%15ld\n",x,y,dxy[k]);
                ++k;
            }
    z=k;
    printf("z, R2, maxr %10ld%10d%10d\n",z,R2,maxr);

//calculate the probabilities (pro) and the cumulative probabilities pr[i]

    sum = pro = exp(z*log(1-PROB));
    i=0;
    pr[i]=1-sum;
    printf("%10ld%16.8e%16.8e%16.8e\n",i,pro,1-sum,pr[i]);
    for (i=1;i<=z;++i)
    {
        pro = pro*PROB/(1-PROB)*(z-i+1)*1.0/i;
        sum += pro;
        pr[i] = 1-sum;
        printf("%10ld%16.8e%16.8e%16.8e\n",i,pro,1-sum,pr[i]);
    }
   // if(z<32)pr[z]=0;  ? 

//xmax is used to make sure the cluster does not grow to the edge of the lattice (at least in the x direction)

    xmax = HEIGHT/2;

//initialize arrays 

    for (i = 0; i < 64; ++i)
        bin[i] = count[i] = count2[i] = ntrials[i] = 0;
    
    for(xy=0; xy<HEIGHT*WIDTH; ++xy)
        lat2[xy]=0;

    gs=0;  //gs it the total number of growth sites

//main loop

    for (runs = 1; runs <= RUNSMAX; ++runs)
    {
        xo = yo = HEIGHT/2;
        xyo = (WIDTH/2)+(HEIGHT/2)*WIDTH;
        
        lat2[xyo] = runs;  // here runs indicates site has been visited
 
        {
            gptr = pptr = trials = 0; //stack (queue)
            PutOnStack2(xyo)
            do
            {
                GetFromStack2(xy)  //new growth site
                ++gs; //growthsites

//loop through i times, where i is determined by pr[i] and rnd

                rnd = NewRandomInteger/2147483648.0;
                if(rnd==0) {printf("rnd=0"); rnd=1/2147483648.0;}  // was having problems with rnd=0, which is possible
                i = 0;
                while (rnd < pr[i]) 
                {

//choose dxy[i] by a random permutation of the order of the list of vectors

                    j = i + (int)((z-i)*(NewRandomInteger/2147483648.0));
                    temp = dxy[i];
                    dxy[i] = dxy[j];
                    dxy[j] = temp;
               
                    xyp = xy + dxy[i];
    
                    if ((xyp&W)>xmax) xmax=(xyp&W);
                    if (lat2[xyp] < runs)
                    {
                        PutOnStack2(xyp)
                        lat2[xyp] = runs;
                    }
                    ++i;
                }
                ++count[i]; //count checks if we got the correct number of neighbors
            }
            while ((gptr != pptr) && (pptr<MAX)); //while queue is not empty
            ++bin[(int)(log(pptr)/log(1.9999999))];
        }
        if((runs & PRINTFREQ) == 0)
        {
            double s=runs*HEIGHT*WIDTH;
            printf("%14.8f%12ld%10d%10d%10d%10ld\n", PROB, runs, SEED, HEIGHT, WIDTH, xmax);
            for (i = 0; i < 64; ++i)
                if (bin[i]) printf("%10ld%12ld%18.10e%18.10e%18.10e\n", i, bin[i], bin[i]*1.0/runs,count[i]*1.0/gs, ntrials[i]*1.0/count2[i]);
    
            fp1 = fopen(OUTFILE,"w");
            fprintf(fp1,"%14.8f%12ld%10d%10d%10d%10ld\n", PROB, runs, SEED, HEIGHT, WIDTH, xmax);
            for (i = 0; i < 32; ++i)
            if (bin[i]) fprintf(fp1,"%10ld%12ld%18.10e\n", i, bin[i], bin[i]*1.0/runs);
            fclose(fp1);
      
        } // if runs
   } // for runs
}

void randinit(long seed)   //need to change this to 64 bits (long long)
{
    double a, ee = -1 + 1/2147483648.0;
    long i;
    extern long nd, ra[M+1];
    
    a = seed/2147483648.0;
    for (nd = 0; nd <= M; nd++)
    {
        a *= 16807;
        a += ee * (long)(a);
        if (a >= 1) a += ee;
        ra[nd] = (long) (2147483648.0 * a);
    }
    nd = M;
    for(i = 0; i<100001; i++)
        NewRandomInteger;
    for (i = 0; i < 10;++i) printf("rnd %20ld\n",ra[i]);
}

