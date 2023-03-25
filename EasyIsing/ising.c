
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define RAND (((double)rand())/(1.0+RAND_MAX))


char **S;
int L;

void
initSystem (size_t LL, int seed, char *lat)
{
  L = LL;
  srand (seed);
  S = (char **) malloc (L * sizeof (void *));
  int i;
  for (i = 0; i < L; i++)
    S[i] = &(lat[i * L]);
}


int
update (int MCS, double beta, double h)
{
  int c0;
  int c1;
  for (c0 = 0; c0 < MCS; c0++)
    for (c1 = 0; c1 < L * L; c1++)
      {
        int i;
        int j;
        i = (int) floor (RAND * L);
        j = (int) floor (RAND * L);
        double de;
        de = S[i][(j + 1) % L];
        de += S[i][(j - 1 + L) % L];
        de += S[(i + 1) % L][j];
        de += S[(i - 1 + L) % L][j]; 
        de *= 2 * S[i][j];
        de += 2.0*h*S[i][j];
        double p;
        p = exp(-beta*de);
        if (RAND < p)
            S[i][j] *= -1;
      }
  return 0;
}
