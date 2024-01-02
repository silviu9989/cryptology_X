/****************************************************************/
/* rho.c                                                        */
/* Authors: Alain Couvreur, Maxime Bombar                       */
/* alain.couvreur@lix.polytechnique.fr                          */
/* maxime.bombar@inria.fr                                       */
/* Last modification October 24, 2022                           */
/****************************************************************/

#include <stdio.h>
#include <assert.h>

#include "gmp.h"
#include "utils.h"

#include "rho.h"


/* to be filled in */

int PollardRho_with_long(long *factor, const long N, long nbOfIterations) 
{
    int status = FACTOR_ERROR;
    long x, y, d, i, aux, contor=0;
    x=1; y=1; d=1;
    while(d==1 && contor<nbOfIterations)
    {
        contor++;
        x=(x*x+7)%N;
        y=(y*y+7)%N; y=(y*y+7)%N;
        if(x>y) aux=x-y;
        else aux=y-x;
        for(i=1; i <= aux && i <= N; i++)
        {
            if(aux%i==0 && N%i==0)
                d = i;
        }
    }  
    if(d==N || contor>=nbOfIterations) status=FACTOR_NOT_FOUND;
    else 
    {    
        *factor=d;
        status = FACTOR_FOUND;
    }

    return status;
}

int PollardRhoSteps(mpz_t factor, const mpz_t N, void (*f)(mpz_t, mpz_t, const mpz_t), long nbOfIterations) 
{
    int status = FACTOR_ERROR;
    int contor=0;
    mpz_t x, y, d, aux;
    mpz_inits(x, y, d, aux, NULL);
    mpz_set_ui(x, 1);
    mpz_set_ui(y, 1);
    mpz_set_ui(d, 1);
    
    while(mpz_cmp_ui(d, 1) == 0 && contor<nbOfIterations)
    {
        contor++;
        f(x, x, N);
        f(y, y, N); f(y, y, N);
        if(mpz_cmp(x, y) > 0) 
            mpz_sub(aux, x, y);
        else mpz_sub(aux, y, x);
        mpz_gcd(d, aux, N);
    }
    if(mpz_cmp(d, N) == 0 || contor>=nbOfIterations) 
    {
        status=FACTOR_NOT_FOUND;
        mpz_clears(x, y, d, NULL);
    }
    else
    {
        mpz_set(factor, d);
        mpz_clears(x, y, d, aux, NULL);
        status=FACTOR_FOUND;
    }
    return status;
}

int PollardRho(factor_t *result, int *nf, const mpz_t N, void (*f)(mpz_t, mpz_t, const mpz_t), long nbOfIterations) 
{
    int status = FACTOR_ERROR;
    int counter=1, prim_status=0, aux=1;
    mpz_t factor;
    mpz_inits(factor, NULL);
    
    prim_status=PollardRhoSteps(factor, N, f, nbOfIterations);
    if(prim_status == FACTOR_FOUND)
    {
        aux=mpz_probab_prime_p(factor, 25);
        if(aux >=1) aux=1;
        else aux=FACTOR_IS_COMPOSITE;
        AddFactor(result, factor, counter, aux);
        *nf=1;
        mpz_clears(factor, NULL);
        status=FACTOR_FOUND;
    }
    else 
    {
        status=FACTOR_NOT_FOUND;
        mpz_clears(factor, NULL);
    }
    return status;
}





