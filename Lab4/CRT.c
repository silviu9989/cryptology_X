#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "gmp.h"

#include "xgcd.h"
#include "CRT.h"

/* Given (r0, m0) and (r1, m1), compute n such that
   n mod m0 = r0; n mod m1 = r1.  If no such n exists, then this
   function returns 0. Else returns 1.  The moduli m must all be positive.
*/
int CRT2(mpz_t n, mpz_t r0, mpz_t m0, mpz_t r1, mpz_t m1){
    int status = 1;
/* to be filled in */

	mpz_t g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim;
	mpz_inits(g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim, NULL);
	int auxx=0;
	//gmp_printf("%Zd %Zd %Zd %Zd %Zd\n", n, r0, m0, r1, m1);
	mpz_fdiv_q (aux1, r1, r0);
	mpz_fdiv_q (aux2, m1, m0);
	if(mpz_cmp(aux1, aux2)==0)
	{
		mpz_fdiv_r (aux1, r1, r0);
		mpz_fdiv_r (aux2, m1, m0);
		if(mpz_cmp(aux1, aux2)==0)
		{
			mpz_clears(g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim, NULL);
			return 0;
		}
	}
	XGCD(g, u, v, m0, m1);
	if(mpz_cmp_ui(g, 1)==0)
	{	
		mpz_mul(aux1, u, m0); mpz_mul(aux1, aux1, r1);
		mpz_mul(aux2, v, m1); mpz_mul(aux2, aux2, r0);
		mpz_add(n, aux1, aux2);
		mpz_lcm(mcom, m0, m1);
		mpz_mod(n, n, mcom);
	}
	else
	{
		mpz_lcm(a, m0, m1);
		mpz_tdiv_q(am0, a, m0);
		mpz_tdiv_q(am1, a, m1);
		mpz_sub(aux1, am0, am1);
		mpz_mul(r0_prim, am0, r0);
		mpz_mul(m0_prim, am0, m0);
		mpz_mul(r1_prim, am1, r1);
		mpz_mul(m1_prim, am1, m1);
		mpz_sub(aux2, r0_prim, r1_prim);
		auxx=linear_equation_mod(n, aux1, aux2, m0_prim);
		if(auxx==0) 
		{
			mpz_clears(g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim, NULL);
			return 0;
		}
		mpz_mod(n, n, a);
		if(mpz_cmp_ui(n, 0)==0) 
		{
			mpz_clears(g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim, NULL);
			return 0;
		}
	}
	mpz_clears(g, u, v, aux1, aux2, mcom, a, am0, am1, m0_prim, m1_prim, r0_prim, r1_prim, NULL);
    	return status;
}

/* to be filled in */

/* Given a list S of pairs (r,m), returns an integer n such that n mod
   m = r for each (r,m) in S.  If no such n exists, then this function
   returns 0. Else returns 1.  The moduli m must all be positive.
*/
int CRT(mpz_t n, mpz_t *r, mpz_t *m, int nb_pairs)
{
	int status = 1, i=0;
	mpz_t aux, maux;
	mpz_inits(aux, maux, NULL);
	
    	if(CRT2(aux, r[0], m[0], r[1], m[1])==0) {mpz_clears(aux, maux, NULL);return 0;}
    	else CRT2(n, r[0], m[0], r[1], m[1]);
    	
    	mpz_set(maux, m[0]);
	for(i=1; i<nb_pairs-1; i++)
	{	
		mpz_mul(maux,maux,m[i]);
		if(CRT2(aux, n, maux, r[i+1], m[i+1])==0) {mpz_clears(aux, maux, NULL);return 0;}
		else
		{
			CRT2(n, n, maux, r[i+1], m[i+1]);
		}
	}
    

	mpz_clears(aux, maux, NULL);
   	return status;
   	/* to be filled in */
}
