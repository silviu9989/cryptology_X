#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gmp.h"
#include "utils.h"
#include "trialdiv.h"

/* OUTPUT: 1 if factorization finished. */

int trialDivision(factor_t* factors, int *nf, mpz_t cof, const mpz_t N, const long bound, uint length, FILE* ficdp)
//NOTE: I do not use the file with prime numbers. Instead, I will use the convenient mpz_nextprime function
{
    int status = FACTOR_NOT_FOUND;
    int counter=0, k=0;
    mpz_t aux, c, aux2;
    mpz_inits(aux, c, aux2, NULL);
    mpz_set_ui(c, 2);
    mpz_set(aux, N);
    
    counter=mpz_probab_prime_p(aux, 20);
    if(counter >= 1)
            {
                mpz_clears(aux, c, aux2, NULL);
                status=FACTOR_NOT_FOUND;
                return status;
            }
            
    counter=0;
    mpz_set(aux, N);
    while(mpz_cmp_ui(c, bound) < 0 && mpz_cmp_ui(aux, 1) != 0)
    {
        mpz_mod(aux2, aux, c);
        while(mpz_cmp_si(aux2, 0) == 0)
        {
            counter++;
            mpz_divexact(aux, aux, c);
            mpz_mod(aux2, aux, c);
        }
        if(counter!=0) 
        {
            AddFactor(factors + k, c, counter, 1);
            k++;
            counter=0;
        }
        if(k == length)
        {
            mpz_clears(aux, c, aux2, NULL);
            status=FACTOR_ERROR;
            return status;
        }
        mpz_nextprime(c, c);
    }
    *nf=k;
    if(mpz_cmp_ui(aux, 1) != 0) 
        {
            mpz_set(cof, aux);
            counter=mpz_probab_prime_p(cof, 20);
            if(counter >=1)
            {
                mpz_clears(aux, c, aux2, NULL);
                status=FACTOR_NOT_FOUND;
                return status;
            }
        }
    else 
        {    
            mpz_set_ui(cof,1);
            mpz_clears(aux, c, aux2, NULL);
            status=FACTOR_FOUND;
            return status;
        }
    return status;
}