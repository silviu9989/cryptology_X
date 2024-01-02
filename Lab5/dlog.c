/****************************************************************/
/* dlog.c                                                       */
/* Authors: Alain Couvreur, FMorain                             */
/* alain.couvreur@lix.polytechnique.fr                          */
/* Last modification October 24, 2022                           */
/****************************************************************/

#include <stdio.h>
#include <assert.h>

#include "hash.h"
#include "dlog.h"
#include "gmp.h"
#include "stdlib.h"
#define DEBUG 0

int babySteps(mpz_t result, hash_table H, mpz_t u, mpz_t g, mpz_t p)
{
    int res = 0;
    int *addr;
    addr=(int*)malloc(1*sizeof(int*));
    mpz_t aux, aux2;
    mpz_inits(aux, aux2, NULL);
    while(mpz_cmp(aux, u) < 0)
    {
        mpz_powm(aux2, g, aux, p);
        hash_put_mpz(H, addr, aux2, aux, g, p);
        mpz_add_ui(aux, aux, 1);
    }
    mpz_clears(aux, aux2, NULL);
    free(addr);
    return res+1;
}


int giantSteps(mpz_t result, hash_table H, mpz_t u, mpz_t g, mpz_t p, mpz_t a)
{
    int res = 0;
    int result1;
    mpz_t aux, aux2, c, vz;
    mpz_inits(aux, aux2, c, vz, NULL);
    
    result1 = mpz_invert(aux, g, p);
    if(result1==0) 
    {mpz_clears(aux, aux2, c, vz, NULL); return res;}
    mpz_powm(aux, aux, u, p);
    result1=HASH_NOT_FOUND;
    
    while(result1==HASH_NOT_FOUND && mpz_cmp(c, u) < 0)
    {
        mpz_powm(aux2, aux, c, p);
        mpz_mul(aux2, aux2, a);
        mpz_mod(aux2, aux2, p);
        result1=hash_get_mpz(vz, H, aux2, g, p);
        mpz_add_ui(c, c, 1);
    }
    if(result1==HASH_NOT_FOUND) 
    {mpz_clears(aux, aux2, c, vz, NULL); return res;}
    
    mpz_sub_ui(c, c, 1);
    mpz_mul(aux, u, c);
    mpz_add(vz, vz, aux);
    mpz_set(result, vz);
    mpz_clears(aux, aux2, c, vz, NULL);
    return res+1;
}

int BSGS(mpz_t result, mpz_t a, mpz_t g, mpz_t p)
{
    int res = 0;
    mpz_t aux, aux2, c, vz, u;
    mpz_inits(aux, aux2, c, vz, u, NULL);
    mpz_sub_ui(aux, p, 1);
    mpz_sqrt(u, aux);
    hash_table H;
    H=hash_init(mpz_get_ui(u));

    
    babySteps(aux, H, u, g, p);
    res=giantSteps(result, H, u, g, p, a);
    mpz_clears(aux, aux2, c, vz, u, NULL);
    hash_clear(H);
    return res;
} 