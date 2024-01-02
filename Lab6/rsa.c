#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"
#include "crt.h"
#include "rsa.h"

#define DEBUG 0
//functia mea
int linear_equation_mod_maftei(mpz_t x, mpz_t a, mpz_t b, mpz_t m){
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
// end functia mea

int is_valid_key(mpz_t p, mpz_t q, mpz_t e, mpz_t d, int nlen, int sec){
    mpz_t lambda, pm1, qm1, pmq, ed, g, bound1, bound2;
    mpz_inits(lambda, pm1, qm1, pmq, ed, g, bound1, bound2, NULL);

    mpz_sub_ui(pm1, p, 1);
    mpz_sub_ui(qm1, q, 1);
    mpz_lcm(lambda, pm1, qm1);

    int status = 1;
    // p, q should be prime
    if(!mpz_probab_prime_p(p, 25) || ! mpz_probab_prime_p(q, 25)){
#if DEBUG
      printf("[is_valid_key] : p, q should be prime.  ");
#endif
      status = 0;
      goto end_valid;
    }

    // e should be odd
    if(mpz_divisible_ui_p(e, 2)){
#if DEBUG
	  printf("[is_valid_key] : e should be odd.  ");
#endif
      status = 0;
      goto end_valid;
    }

    // e should be prime to lambda
    mpz_gcd(g, e, lambda);
    if(mpz_cmp_ui(g, 1) != 0){
#if DEBUG
      printf("[is_valid_key] : e should be prime to lambda. ");
#if DEBUG >= 2
	  gmp_printf("gcd(e, lambda) = %Zd.  ", g);
#endif
#endif
	  status = 0;
      goto end_valid;
    }
    // Bounds on e.
    size_t size_e = mpz_sizeinbase(e, 2);
    if(size_e < 16){
#if DEBUG
      printf("[is_valid_key] : e is too small.   ");
#endif
      status = 0;
      goto end_valid;
    }
    if(size_e > 256){
#if DEBUG
      printf("[is_valid_key] : e is too large.   ");
#endif
      status=0;
      goto end_valid;
    }

    // p, q should be large enough
    mpf_t b_f, tmp_f;
    mpf_inits(b_f, tmp_f, NULL);
    mpf_sqrt_ui(b_f, 2);
    mpf_set_ui(tmp_f, 2);
    mpf_pow_ui(tmp_f, tmp_f, nlen/2 - 1);
    mpf_mul(b_f, b_f, tmp_f);
    mpz_set_f(bound1, b_f);	
    if(mpz_cmp(p, bound1) <= 0 || mpz_cmp(q, bound1) <= 0){
#if DEBUG
	printf("[is_valid_key] : p or q is too small.  ");
#endif
	status = 0;
    }
    mpz_clear(bound1);
    mpf_clears(b_f, tmp_f, NULL);
    if(status == 0)
	  goto end_valid;

    // p, q should not be too close to each other
    mpz_inits(bound2, pmq, NULL);
    mpz_ui_pow_ui(bound2, 2, nlen/2 - sec);
    mpz_sub(pmq, p, q);
    if(mpz_cmpabs(pmq, bound2) <= 0){
#if DEBUG
	  printf("[is_valid_key] : p and q are too close.  ");
#endif
      status = 0;
    }
    mpz_clears(bound2, pmq, NULL);
    if(status == 0)
	  goto end_valid;

    // ed = 1 mod lambda
    mpz_mul(ed, e, d);
    mpz_mod(ed, ed, lambda);
    if(mpz_cmp_ui(ed, 1) != 0){
#if DEBUG
	  printf("[is_valid_key] : ed should be");
	  printf(" congruent to 1 modulo lambda.   ");
#endif
	  status = 0;
      goto end_valid;
    }
 end_valid:
    mpz_clears(lambda, pm1, qm1, ed, g, NULL);

    return status;
}


void RSA_weak_generate_key(mpz_t p, mpz_t q, mpz_t e, mpz_t d, int nlen,
			   gmp_randstate_t state)
{
        mpz_t aux, N, paux, qaux, lambdaN, const1, e1, e2;
	mpz_inits(aux, N, paux, qaux, lambdaN, const1, e1, e2, NULL);
	mpz_set_ui(const1, 1);
	
	int constant=nlen;
	constant=constant/2; constant=constant-1;
	mpz_ui_pow_ui(aux, 2, constant);
	mpz_mul_ui(aux, aux, 3);
	mpz_fdiv_q_ui(aux, aux, 2); //sqrt2 ~ 3/2
        mpz_urandomb(p, state, nlen/2);
	mpz_nextprime(p, p);
	
	while(mpz_cmp(p, aux) <= 0)
	{
	    mpz_urandomb(p, state, nlen/2);
	    mpz_nextprime(p, p);
	}
	
	mpz_urandomb(q, state, nlen/2);
	mpz_nextprime(q, q);
	
	while(mpz_cmp(q, aux) <= 0 || mpz_cmp(p, q) == 0)
	{
	    mpz_urandomb(q, state, nlen/2);
	    mpz_nextprime(q, q);
	}
	mpz_mul(N, p, q);
	mpz_sub_ui(paux, p, 1);
	mpz_sub_ui(qaux, q, 1);
	mpz_lcm(lambdaN, paux, qaux);
	
	mpz_urandomb(e, state, nlen);
	mpz_nextprime(e, e);
	
	mpz_mod(paux, lambdaN, e);
	mpz_ui_pow_ui(e1, 2, 16);
	mpz_ui_pow_ui(e2, 2, 256);
	while(mpz_cmp_ui(paux, 0)==0 || mpz_cmp(e, e1) < 0 || mpz_cmp(e, e2) > 0)
	{
	    mpz_urandomb(e, state, 256);
	    mpz_nextprime(e, e);
	    mpz_mod(paux, lambdaN, e);
	}
	linear_equation_mod_maftei(d, e, const1, lambdaN);
	mpz_clears(aux, N, paux, qaux, lambdaN, const1, e1, e2, NULL);
}


void RSA_generate_key(mpz_t N, mpz_t p, mpz_t q, mpz_t e, mpz_t d,
		      int nlen, int sec, gmp_randstate_t state){
    int done=0;
/* to be filled in */
    do{
	RSA_weak_generate_key(p, q, e, d, nlen, state);
	done=1;
    }while(!is_valid_key(p, q, e, d, nlen, sec) && done);
    mpz_mul(N, p, q);
}



void RSA_encrypt(mpz_t cipher, mpz_t msg, mpz_t N, mpz_t e){
    mpz_powm(cipher, msg, e, N);
}

void RSA_decrypt(mpz_t msg, mpz_t cipher, mpz_t N, mpz_t d){
    mpz_powm(msg, cipher, d, N);
}


/* Use CRT. */
void RSA_decrypt_with_p_q(mpz_t msg, mpz_t cipher, mpz_t N, mpz_t d,
			  mpz_t p, mpz_t q)
{
        mpz_t dP, dQ, qinv, aux, m1, m2, h, hq;
	mpz_inits(dP, dQ, qinv, aux, m1, m2, h, hq, NULL);
	mpz_sub_ui(aux, p, 1);
	mpz_mod(dP, d, aux);
	mpz_sub_ui(aux, q, 1);
	mpz_mod(dQ, d, aux);
	mpz_invert(qinv, q, p);
	
	mpz_powm(m1, cipher, dP, p);
	mpz_powm(m2, cipher, dQ, q);
	
	mpz_sub(h, m1, m2);
	mpz_mul(h, qinv, h);
	mpz_mod(h, h, p);
	mpz_mul(hq, h, q);
	mpz_add(msg, m2, hq);
	mpz_mul(aux, p, q);
	mpz_mod(msg, msg, aux);
	
	mpz_clears(dP, dQ, qinv, aux, m1, m2, h, hq, NULL);
}



void RSA_dummy_generate_key(mpz_t N, mpz_t e, int nlen,
			    gmp_randstate_t state){
	
    mpz_t g, p, q, lambda, pm1, qm1;
    mpz_inits(g, p, q, lambda, pm1, qm1, NULL);

    do{
	mpz_urandomb(p, state, nlen/2);
	mpz_nextprime(p, p);		
	do{
	    mpz_urandomb(q, state, nlen/2);
	    mpz_nextprime(q, q);
	}while(mpz_cmp(p, q) == 0);
		
	mpz_sub_ui(pm1, p, 1);
	mpz_sub_ui(qm1, q, 1);
	mpz_lcm(lambda, pm1, qm1);
	mpz_gcd(g, e, lambda);
	mpz_mul(N, p, q);
	
    }while(mpz_cmp_ui(g, 1) != 0 || mpz_sizeinbase(N, 2) < nlen);
#if DEBUG
    gmp_printf("[RSA_dummy_generate_key] N = %Zd, size : %ld\n",
	       N, mpz_sizeinbase(N, 2));
#endif
	
    mpz_clears(g, p, q, lambda, pm1, qm1, NULL);
}