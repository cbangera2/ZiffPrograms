#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "FourTapGen.h"
#include "CustomQueue.h"
#include <iostream>

using namespace std;
//g++ percpop21unionwrapRZ.cpp -O3 -std=c++11
//finding nc for finite widths --
//this program is for bond percolation on a square lattice.  Note DIRMAX = 4
//calculates the number of clusters and number of wrapping clusters

#define HEIGHT 1048576      //power of two
#define WIDTH 4 //keep at even width
#define RUNSMAX  100000000
#define PRINTFREQ 200
#define PROB  0.5 //Pc 0.5927460508 square site // 0.3472963553  triangular bond//
#define SEED 146 //change value every time (otherwise same answer for same width)
#define OUTFILE "percpop21triwrap16b"
#define DIRMAX 6 //4=square, 6=triangular lattice
#define H (HEIGHT-1)
#define W (WIDTH-1)
#define S 65535  //for stack (queue) -- don't change
#define M 16383 //for random number generator -- don't change

long    lat[WIDTH][HEIGHT], lat2[WIDTH][HEIGHT];     //for lattice

int dirmax;//new dirmax variable to account for union jack lattice's 2 orientations

int main(void)
{
	/*
	Instructions:
		-put outer point on queue
		-get (x,y,z,w) from queue
		-use Cn to find # of neighbors k
			-Cn is sum of probabilities
			-P_n= combination(z,m)*p^n *(1-p)^(z-n)  //precalculate this at start
		-pseudocode
			set<int> visit;
			queue<int> q;
			for(int i=0; i<k ++i){
				pick a unique vector dx randomly
				from list of 0-59 possible nn vector (int)(60*random#)
				list of previosly visited neigbors
				x'=x+dx;
				check if x' is in set/map to see if visited before
				if(visit.find(num)!=visit.end()){//if visited before
					break;
				}else{
					set.push(num); //add to set
					q.push(num); //add to queue as well
				}
			}
			//after this find# of clusters in each  bin as that is equal to number of occupied sites
	*/
	FourTap gen(SEED);
	long prob = (long)(gen.max() * PROB);

	long dx1[8] = { 1, 1, 0, -1, -1, -1, 0, 1 }, dy1[8] = { 0, 1, 1, 1, 0, -1, -1, -1 }; //dx,dy for unionjack orientation with 8 nearest neighbors
	long dx2[4] = { 1, 0, -1, 0 }, dy2[4] = { 0, 1, 0, -1 }; //dx,dy for unionjack orientation with 4 nearest neighbors
	long nwraps, ymax, ymin;

	long bin[4096], big, nc;
	double  nctot, ncsqrtot;
	FILE* fp1;

	for (int i = 0; i < 4096; ++i) bin[i] = 0;
	nwraps = nctot = 0;

	//clear the lattice

	for (long x = 0; x < WIDTH; x += 1)
		for (long y = 0; y < HEIGHT; y += 1)
			lat2[x][y] = 0;

	for (long runs = 1; runs <= RUNSMAX; ++runs) {
		for (long x = 0; x < WIDTH; x += 1)
			for (long y = 0; y < HEIGHT; y += 1)
				if (gen() < prob) lat[x][y] = 0;
				else lat[x][y] = -1;

		nc = big = 0;
		CustomQueue latticeQueue(S);

		for (long xo = 0; xo < WIDTH; xo += 1)
			for (long yo = 0; yo < HEIGHT; yo += 1)
				if (lat[xo][yo] == 0) {
					//New cluster
					ymax = ymin = yo;
					++nc;
					lat[xo][yo] = nc;
					{
						bool wf = false;
						long count = 1;
						lat2[xo][yo] = xo;
						latticeQueue.clear();
						latticeQueue.push(xo, yo);
						do {
							long x, y;
							latticeQueue.pop(x, y);
							if ((x + y) % 2 == 0) { dirmax = 8; }//if even 8 nearest neighbors
							else { dirmax = 4; } //if odd 4 nearest neighbors
							for (long dir = 0; dir < dirmax; ++dir) {
								long xp, yp;
								if (dirmax == 8) {//using dirmax assign the correct dy and dx matrix
									xp = x + dx1[dir];
									yp = y + dy1[dir];
								}
								else {
									xp = x + dx2[dir];
									yp = y + dy2[dir];
								}
								long xpp = (xp + 16 * WIDTH) % WIDTH;
								long ypp = yp & H;
								if (lat[xpp][ypp] == 0) {
									if (yp > ymax) ymax = yp;
									else if (yp < ymin) ymin = yp;

									count++;

									lat2[xpp][yp & H] = xp;
									lat[xpp][ypp] = nc;
									latticeQueue.push(xp, yp);
								}
								else if (lat[xpp][ypp] == nc) if (xp != lat2[xpp][ypp]) if (!wf) {
									wf = true;
								}//++nwraps;}
							}
						} while (!latticeQueue.empty()); //while queue is not empty
						if (wf) ++bin[(int)(log(ymax - ymin + 1) / log(1.9999999))];
						if (wf) ++nwraps;
						//++bin[ymax-ymin];
					  //  if(wf) ++bin[(int)(log(count)/log(1.9999999))];
						//  ++bin[(int)(log(count)/log(1.9999999))];
					}
				} //if lat < runs

		nctot += nc;
		ncsqrtot += nc * 1.0 * nc;

		if ((runs % PRINTFREQ) == 0) {
			printf("%14.6f%12ld%10d%10d%10d\n", PROB, runs, SEED, HEIGHT, WIDTH);
			printf("wraps %10ld%20.8f\n", nwraps, nwraps * WIDTH * 1.0 / (HEIGHT * runs));
			double error = sqrt(ncsqrtot / runs - (nctot / runs) * (nctot / runs)) / (HEIGHT * WIDTH);
			double sites = runs * 1.0 * HEIGHT * WIDTH;
			printf("clusters %28.16e%28.16e\n", nctot / sites, error / sqrt(runs));
			for (int i = 0; i <= 12; ++i)  printf("%10d%18.10e\n", i, bin[i] * 1.0 / runs);

			fp1 = fopen(OUTFILE, "w");
			fprintf(fp1, "%14.6f%12ld%10d%10d%10d\n", PROB, runs, SEED, HEIGHT, WIDTH);
			fprintf(fp1, "wraps %10ld%20.8f\n", nwraps, nwraps * WIDTH * 1.0 / (HEIGHT * runs));
			fprintf(fp1, "clusters %28.16e%28.16e\n", nctot / sites, error / sqrt(runs));
			for (int i = 0; i <= 12; ++i) fprintf(fp1, "%10d%18.10e\n", i, bin[i] * 1.0 / runs);
			fclose(fp1);
		} // if runs
	} // for runs
}