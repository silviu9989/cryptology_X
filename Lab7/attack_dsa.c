#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gmp.h"
#include "buffer.h"
#include "sha3.h"

#include "rsa.h"
#include "sign.h"
#include "dsa.h"
#include "attack_dsa.h"

#define DEBUG 1

int linear_equation_mod2(mpz_t x, mpz_t a, mpz_t b, mpz_t m)
{
    int status = 1;
    
	mpz_t aminus1, aux, const0, gcda;
	mpz_inits(aminus1,aux,const0,gcda,NULL);
	
	if(mpz_cmp(a,const0)!=0 && mpz_cmp(b,const0)!=0) 
	{   
	    mpz_gcd(gcda, a, b);
	    mpz_gcd(gcda, gcda, m);
	    mpz_tdiv_q(a,a,gcda);
	    mpz_tdiv_q(b,b,gcda);
	}
	if(mpz_invert(aminus1,a,m)==0) {mpz_clears(aminus1,gcda,aux,const0,NULL); status=0; return status;}
	mpz_invert(aminus1,a,m);
	mpz_mul(x,aminus1,b);
	mpz_mod(x,x,m);
	mpz_clears(aminus1,gcda,aux,const0,NULL);
	return status;
}

void dsa_sign_dummy(buffer_t *msg, mpz_t p, mpz_t q, mpz_t a, mpz_t x, mpz_t r, mpz_t s, mpz_t k)
{
/* to be filled in */
    mpz_t hash, aux, kminus1, const1;
    mpz_inits(hash, aux, kminus1, const1, NULL);
    buffer_t hashaux;
    buffer_init(&hashaux, hash_length(q));
    int aux2;
    mpz_set_ui(const1, 1);

    aux2=linear_equation_mod2(kminus1, k, const1, q);
    if(aux2 == 0) printf("Error: k does not have an inverse ...");
    mpz_powm(r, a, k, p);
    mpz_mod(r, r, q);
    
    buffer_hash(&hashaux, hash_length(q), msg);
    mpz_import(hash, hash_length(q), 1, 1, 1, 0, hashaux.tab);
    mpz_mul(aux, x, r);
    mpz_add(aux, hash, aux);
    mpz_mul(aux, kminus1, aux);
    mpz_mod(s, aux, q);
    
    buffer_clear(&hashaux);
    mpz_clears(hash, aux, kminus1, const1, NULL);
}

void solve_system_modq(mpz_t x, mpz_t r1, mpz_t s1, mpz_t r2, mpz_t s2, mpz_t h1, mpz_t h2, mpz_t q)
{

    /* Solves the system with unknowns k, x:
       s1.k - r1.x = h1
       s2.k - r2.x = h2
       and fills in x
    */
    mpz_t aux1, aux2, aux3, const1;
    mpz_inits(aux1, aux2, aux3, const1, NULL);
    mpz_set_ui(const1, 1);
    int a=0;
    gmp_printf("Linear system :\n%Zd k - %Zd x = %Zd\n%Zd k - %Zd x = %Zd\n\n", s1, r1, h1, s2, r2, h2);
    if(mpz_cmp_ui(s1, 0) == 0 || mpz_cmp_ui(s2, 0) == 0)
    {
        if(mpz_cmp_ui(s1, 0) == 0)
        {
            mpz_mul_si(aux2, r1, -1);
            a=linear_equation_mod2(aux1, aux2, const1, q);
            if(a == 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
            else
            {
                mpz_mul(x, h1, aux1);
                
                mpz_mod(x, x, q);
                gmp_printf("Candidate for secret key obtained from the attack:\nx = %Zd\n\n", x);
                mpz_mul(aux1, x, r2);
                mpz_add(aux1, h2, aux1);
                if(mpz_cmp_ui(s2, 0) == 0 && mpz_cmp_ui(aux1, 0) != 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
                else if(mpz_cmp_ui(s2, 0) != 0)
                {
                    a=linear_equation_mod2(aux2, s2, const1, q);
                    if(a == 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
                }
            }
        }
        else 
        {
            mpz_mul_si(aux2, r2, -1);
            a=linear_equation_mod2(aux1, aux2, const1, q);
            if(a == 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
            else
            {
                mpz_mul(x, h2, aux1);
                mpz_mod(x, x, q);
                gmp_printf("Candidate for secret key obtained from the attack:\nx = %Zd\n\n", x);
                mpz_mul(aux1, x, r1);
                mpz_add(aux1, h1, aux1);
                if(mpz_cmp_ui(s1, 0) == 0 && mpz_cmp_ui(aux1, 0) != 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
                else if(mpz_cmp_ui(s1, 0) != 0)
                {
                    a=linear_equation_mod2(aux2, s1, const1, q);
                    if(a == 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
                }
            }
        }
    }
    else
    {
    a=linear_equation_mod2(aux1, s1, const1, q);
    if(a == 0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
    else
    {
        mpz_mul(aux1, aux1, s2);
        mpz_mul(aux2, aux1, r1);
        mpz_sub(aux2, aux2, r2);
        mpz_mul(aux1, aux1, h1);
        a=linear_equation_mod2(aux3, aux2, const1, q);
        if(a==0) {printf("\nSystem unsolvable\n"); mpz_clears(aux1, aux2, aux3, const1, NULL);}
        else
        {
            mpz_sub(x, h2, aux1);
            mpz_mul(x, x, aux3);
            mpz_mod(x, x, q);
            gmp_printf("Candidate for secret key obtained from the attack:\nx = %Zd\n\n", x);
            mpz_clears(aux1, aux2, aux3, const1, NULL);
        }
    }
    }
/* to be filled in */
}


void dsa_attack(mpz_t x, buffer_t *msg1, buffer_t *msg2, mpz_t p, mpz_t q, mpz_t a, mpz_t r1, mpz_t s1, mpz_t r2, mpz_t s2)
{
/* to be filled in */
    //dsa_sign_dummy(&msg1, p, q, a, x, r1, s1, k);
    //dsa_sign_dummy(&msg2, p, q, a, x, r2, s2, k);
    buffer_t hashaux1;
    buffer_init(&hashaux1, hash_length(q));
    buffer_hash(&hashaux1, hash_length(q), msg1);
    
    buffer_t hashaux2;
    buffer_init(&hashaux2, hash_length(q));
    buffer_hash(&hashaux2, hash_length(q), msg2);
    
    mpz_t hash1, hash2;
    mpz_inits(hash1, hash2, NULL);
    
    mpz_import(hash1, hash_length(q), 1, 1, 1, 0, hashaux1.tab);
    mpz_import(hash2, hash_length(q), 1, 1, 1, 0, hashaux2.tab);
    
    solve_system_modq(x, r1, s1, r2, s2, hash1, hash2, q);
    
    buffer_clear(&hashaux1);
    buffer_clear(&hashaux2);
    mpz_clears(hash1, hash2, NULL);
}

