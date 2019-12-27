// C program to fdnd LCM of two numbers 
#include <stdio.h> 
#include <stdlib.h> 
#include <stdint.h> 



typedef struct {
	uint32_t nbfr;
	int32_t fr[4];
} wspr_fr_set_t;

#define NBFRSET 8

/* wspr freqs
1838100, 3570100, 5366200, 7040100, 10140200, 14097100, 18106100, 21096100, 24926100, 28126100
*/
const wspr_fr_set_t wsprfrsets[NBFRSET]={
{1, { 137500} },
{4, { 1838100, 3570100, 5366200, 7040100 } },
{4, { 3570100, 5366200, 7040100, 10140200 } },
{3, { 7040100, 10140200, 14097100 }  },
{3, { 10140200, 14097100, 18106100 } },
{3, { 14097100, 18106100, 21096100 } },
{3, { 18106100, 21096100, 24926100 } },
{3, { 21096100, 24926100, 28126100 } },
};

// Recursive function to return gcd of a and b 
unsigned int gcd(unsigned int a, unsigned int b) 
{ 
	if (a == 0) 
		return b; 
	return gcd(b % a, a); 
} 


// Driver program to test above function 
int main() 
{ 
  const unsigned int fs=20000000;
  unsigned int fd;
  unsigned int n,s;
  unsigned int min,fm;

  for(s=0;s<NBFRSET;s++) {
  
  min=99999999;fm=0;
  for(fd=wsprfrsets[s].fr[wsprfrsets[s].nbfr-1]+1000000-fs/4;fd<wsprfrsets[s].fr[0]-1000000+fs/4;fd+=100) {
	unsigned int max;
	max=0;
	for(n=0;n<wsprfrsets[s].nbfr;n++) {
		unsigned int v;
		int ff;

		ff=fd-wsprfrsets[s].fr[n]+fs/4;

		if(ff<0) { printf("%d %d\n",fd,ff); exit(1); }
		v=fs/gcd(ff,fs);

		max+=v;
	}
	if (max<min) {
		min=max;
		fm=fd;
	}
  }
   printf("{ %d, %d, { ",wsprfrsets[s].nbfr,fm);
   for(n=0;n<wsprfrsets[s].nbfr;n++) {
		printf("%d, ",wsprfrsets[s].fr[n]);
   }
   printf("} , { ");
   for(n=0;n<wsprfrsets[s].nbfr;n++) {
		int ff,v;

		ff=wsprfrsets[s].fr[n]-fm+fs/4;

		v=fs/gcd(ff,fs);
		printf("%d, ",v);
   }
   printf("}}, \n");
 }
} 

