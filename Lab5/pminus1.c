#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "gmp.h"
#include "utils.h"
#include "pminus1.h"

#define DEBUG 0

static void ReadDifFile(mpz_t p, FILE *file){
    mpz_t q;
    int dp;

    /* we have to read file while read_p <= p */
    mpz_init_set_ui(q, 1);
    while(fscanf(file, "%d", &dp) != EOF){
	mpz_add_ui(q, q, dp << 1);
	if(mpz_cmp(q, p) > 0)
	    break;
    }
    mpz_set(p, q);
    mpz_clear(q);
}

/* Starting from nextprime(p) >= p+1. 
   On unsuccessful exit, p is the smallest prime > bound1.
*/
int PollardPminus1Step1(mpz_t factor, const mpz_t N, long bound1, FILE* ficdp, mpz_t b, mpz_t p)
{
    int status = FACTOR_ERROR;
    long k=0;
    mpz_t q, g, auxp, auxb, aux1, aux2;
    mpz_inits(q, g, auxp, auxb, aux1, aux2, NULL);
    mpz_set(auxp, p);
    if(mpz_cmp_ui(auxp,2)<0) mpz_set_ui(auxp, 2);
    mpz_set_ui(b, 2);
    
    while(mpz_cmp_ui(auxp, bound1) <= 0)
    {
        k=1;
        mpz_pow_ui(aux1, auxp, k);
        mpz_pow_ui(aux2, auxp, k+1);
        while((mpz_cmp_ui(aux2, bound1) < 0))
        {
            k++;
            mpz_pow_ui(aux1, auxp, k);
            mpz_pow_ui(aux2, auxp, k+1);
        }
        mpz_set(q, aux1);
        mpz_powm(b, b, q, N);
        mpz_sub_ui(auxb, b, 1);
        mpz_gcd(g, auxb, N);
        if(mpz_cmp(g, N) == 0)
        {
            status=FACTOR_NOT_FOUND;
        }
        else if(mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, N) != 0)
        {
            mpz_set(factor, g);
            mpz_clears(q, g, auxp, auxb, aux1, aux2, NULL);
            status=FACTOR_FOUND;
            return status;
        }
        mpz_nextprime(auxp, auxp);
    }
    
    mpz_clears(q, g, auxp, auxb, aux1, aux2, NULL);
    return status;
}

int PollardPminus1Step2(mpz_t factor, const mpz_t N, long bound2, FILE* ficdp,
			mpz_t b, mpz_t p){
    mpz_t bm1;
    unsigned long d;
    int dp, status = FACTOR_ERROR;
    int B = (int)log((double)bound2);
    B = B * B;

    mpz_init(bm1);
    ReadDifFile(p, ficdp);
    /* Precomputations */
    mpz_t* precomputations = (mpz_t*)malloc(B * sizeof(mpz_t));
    mpz_t* cursor = precomputations;
    int i;
		
    for(i = 0; i < B; i++, cursor++){
	mpz_init(*cursor);
	mpz_powm_ui(*cursor, b, i, N);
    }
#if DEBUG >= 1
    printf("# Precomputation of phase 2 done.\n");
#endif
    mpz_powm(b, b, p, N);
    while(mpz_cmp_ui(p, bound2) <= 0){
	mpz_sub_ui(bm1, b, 1);
	mpz_gcd(factor, bm1, N);
	if(mpz_cmp_ui(factor, 1) > 0){
	    status = FACTOR_FOUND;
	    break;
	}
	fscanf(ficdp, "%d", &dp);
	d = dp << 1;
	mpz_add_ui(p, p, d);		
	if(d < B){
	    mpz_mul(b, b, precomputations[d]);
	    mpz_mod(b, b, N);
	}
	else{
	    printf("Cramer's rule Failed!\n");
	    printf("WRITE A PAPER!!!\n");
	    return 1;
	}
    }			
    cursor = precomputations;
    for(i = 0; i < B; i++, cursor++){
        mpz_clear(*cursor);
    }
    free(precomputations);
    mpz_clear(bm1);
    if (status != FACTOR_FOUND) {
        status = FACTOR_NOT_FOUND;
    }
    return status;
}

int PollardPminus1(factor_t* res, int *nf, const mpz_t N, long bound1, long bound2, FILE* ficdp)
{
    mpz_t b, p, factor;
    int status=FACTOR_ERROR;
    
    int counter=1, prim_status=0, aux=1;
    mpz_inits(factor, b, p, NULL);
    
    prim_status=PollardPminus1Step1(factor, N, bound1, ficdp, b, p);
    if(prim_status == FACTOR_FOUND)
    {
        aux=mpz_probab_prime_p(factor, 25);
        if(aux >=1) aux=1;
        else aux=FACTOR_IS_COMPOSITE;
        AddFactor(res, factor, counter, aux);
        *nf=1;
        mpz_clears(factor, b, p, NULL);
        status=FACTOR_FOUND;
        return status;
    }
    else 
    {
        status=FACTOR_NOT_FOUND;
        //mpz_clears(factor, b, p, NULL);
    }    
    
    prim_status=PollardPminus1Step2(factor, N, bound2, ficdp, b, p);
    if(prim_status == FACTOR_FOUND)
    {
        aux=mpz_probab_prime_p(factor, 25);
        if(aux >=1) aux=1;
        else aux=FACTOR_IS_COMPOSITE;
        AddFactor(res, factor, counter, aux);
        *nf=1;
        mpz_clears(factor, b, p, NULL);
        status=FACTOR_FOUND;
    }
    else 
    {
        status=FACTOR_NOT_FOUND;
        mpz_clears(factor, b, p, NULL);
    }    
    return status;
}