#include <stdio.h>
#include <assert.h>

#include "gmp.h"

#include "xgcd.h"
#include "CRT.h"
#include "rsa.h"

void generate_probable_prime(mpz_t p, size_t nbits, gmp_randstate_t state){
	mpz_urandomb(p, state, nbits);
	mpz_nextprime(p, p);
}

void RSA_encrypt(mpz_t cipher, mpz_t msg, mpz_t N, mpz_t e)
{
/* to be filled in */
	mpz_powm(cipher, msg, e, N);
}

void RSA_decrypt(mpz_t msg, mpz_t cipher, mpz_t N, mpz_t d)
{
/* to be filled in */
	mpz_powm(msg, cipher, d, N);
}

/* Use CRT. */
void RSA_decrypt_with_p_q(mpz_t msg, mpz_t cipher, mpz_t N, mpz_t d, mpz_t p, mpz_t q)
{
/* to be filled in */
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

/* Generate p and q with nbits. */
void RSA_generate_key(mpz_t N, mpz_t p, mpz_t q, mpz_t e, mpz_t d, size_t nbits, gmp_randstate_t state)
{
	mpz_t const1, paux, qaux, lambdaN;
	mpz_inits(const1, paux, qaux, lambdaN, NULL);
	mpz_set_ui(const1, 1);
	generate_probable_prime(p, nbits, state);
	generate_probable_prime(q, nbits, state);
	while(mpz_cmp(p, q)==0)
	{
		generate_probable_prime(q, nbits, state);
	}
	mpz_mul(N, p, q);
	mpz_sub(paux, p, const1);
	mpz_sub(qaux, q, const1);
	mpz_lcm(lambdaN, paux, qaux);
	generate_probable_prime(e, nbits, state);
	mpz_mod(paux, lambdaN, e);
	while(mpz_cmp_ui(paux, 0)==0)
	{
		generate_probable_prime(e, nbits, state);
		mpz_mod(paux, lambdaN, e);
	}
	linear_equation_mod(d, e, const1, lambdaN);
	
	mpz_clears(const1, paux, qaux, lambdaN, NULL);
/* to be filled in */
}
