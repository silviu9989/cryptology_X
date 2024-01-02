#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gmp.h"
#include "buffer.h"
#include "sha3.h"

#include "rsa.h"
#include "sign.h"
#include "dsa.h"

#define DEBUG 1

void generate_probable_prime(mpz_t p, int psize, gmp_randstate_t state)
{
    do{
	mpz_rrandomb(p, state, psize);
	mpz_nextprime(p, p);
    }while(mpz_sizeinbase(p, 2) < psize);
}

void generate_pq(mpz_t p, mpz_t q, size_t psize, size_t qsize, gmp_randstate_t state)
{
    mpz_t aux;
    mpz_inits(aux, NULL);
    
    int stare=1, n=0;
    size_t size=0;

    generate_probable_prime(q, qsize, state);
    while(stare == 1)
    {
        generate_probable_prime(p, psize, state);
        mpz_fdiv_q(aux, p, q);
        mpz_mul(p, aux, q);
        mpz_add_ui(p, p, 1);
        
        size=mpz_sizeinbase(p, 2);
        n=mpz_probab_prime_p(p, 30);
        if((size == psize) && (n != 0))
        {
            stare=0;
        }
    }
    mpz_clears(aux, NULL);
}
	
void dsa_generate_keys(mpz_t p, mpz_t q, mpz_t a, mpz_t x, mpz_t y, size_t psize, size_t qsize, gmp_randstate_t state)
{
/* to be filled in */
    mpz_t aux, h, auxmod;
    mpz_inits(aux, h, auxmod, NULL);
    
    generate_pq(p, q, psize, qsize, state);
    mpz_sub_ui(aux, p, 1);
    mpz_divexact(aux, aux, q);
    generate_probable_prime(h, qsize, state);
    mpz_add_ui(h, h, 1);
    mpz_powm(a, h, aux, p);
    mpz_mod(auxmod, a, p);
    while(mpz_cmp_ui(auxmod, 1) == 0)
    {
       generate_probable_prime(h, qsize, state);
       mpz_add_ui(h, h, 1);
       mpz_powm(a, h, aux, p); 
    }
    
    mpz_urandomb(x, state, qsize-5);
    mpz_powm(y, a, x, p);
    
    mpz_clears(aux, h, auxmod, NULL);
}


void dsa_generate_key_files(const char* pk_file_name, const char* sk_file_name,
			    size_t psize, size_t qsize,
			    gmp_randstate_t state){
    // 1. INITS
    mpz_t p, q, a, x, y;
    mpz_inits(p, q, a, x, y, NULL);
    FILE* pk = fopen(pk_file_name, "w");
    FILE* sk = fopen(sk_file_name, "w");
	
    // 2. Key generation
    dsa_generate_keys(p, q, a, x, y, psize, qsize, state);

    // 3. Printing files
    fprintf(pk, "#DSA public key (%lu bits, %lu bits):\n", psize, qsize);
    gmp_fprintf(pk, "p = %#Zx\nq = %#Zx\na = %#Zx\ny = %#Zx\n", p, q, a, y);
    fprintf(sk, "#DSA Private Key (%lu bits, %lu bits):\n", psize, qsize);
    gmp_fprintf(sk, "p = %#Zx\nq = %#Zx\na = %#Zx\nx = %#Zx\n", p, q, a, x);
	
    // 4. Cleaning
    mpz_clears(p, q, a, x, y, NULL);
    fclose(pk);
    fclose(sk);
}


void dsa_key_import(const char* key_file_name, mpz_t p, mpz_t q, mpz_t a,
		    mpz_t xy){
    FILE* key = fopen(key_file_name, "r");
	
    // Go to second line, then move from 6 characters to the right
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);

    // Scan the modulus p
    gmp_fscanf(key, "%Zx", p);

    // Same for q
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);
    gmp_fscanf(key, "%Zx", q);

    // Same for a
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);
    gmp_fscanf(key, "%Zx", a);

    // Same for x or y
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);
    gmp_fscanf(key, "%Zx", xy);

    fclose(key);
}

int linear_equation_mod(mpz_t x, mpz_t a, mpz_t b, mpz_t m){
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

void dsa_sign_buffer(buffer_t *msg, mpz_t p, mpz_t q, mpz_t a, mpz_t x, mpz_t r, mpz_t s, gmp_randstate_t state)
{
/* to be filled in */
    mpz_t hash, k, aux, kminus1, const1;
    mpz_inits(hash, k, aux, kminus1, const1, NULL);
    buffer_t hashaux;
    buffer_init(&hashaux, hash_length(q));
    int qsize = mpz_sizeinbase(q, 2);
    int aux2;
    mpz_set_ui(const1, 1);
    
    mpz_urandomb(k, state, qsize-2);
    aux2=linear_equation_mod(kminus1, k, const1, q);
    while(aux2 == 0)
    {
        mpz_urandomb(k, state, qsize-2);
        aux2=linear_equation_mod(kminus1, k, const1, q);
    }
    mpz_powm(r, a, k, p);
    mpz_mod(r, r, q);
    
    buffer_hash(&hashaux, hash_length(q), msg);
    mpz_import(hash, hash_length(q), 1, 1, 1, 0, hashaux.tab);
    mpz_mul(aux, x, r);
    mpz_add(aux, hash, aux);
    mpz_mul(aux, kminus1, aux);
    mpz_mod(s, aux, q);
    
    buffer_clear(&hashaux);
    mpz_clears(hash, k, aux, kminus1, const1, NULL);
}


void dsa_sign(const char* file_name, const char* key_file_name,
	      const char* signature_file_name,
	      gmp_randstate_t state){
    // 1. Initialisation
    mpz_t p, q, a, x, r, s;
    buffer_t msg;
    mpz_inits(r, s, p, q, a, x, NULL);
    buffer_init(&msg, 100);
	
    // 2. Import the message
    buffer_from_file(&msg, file_name);
#if DEBUG
    printf("Length of the message = %lu.\n", msg.length);
#endif

    /* 3. Parse the secret key */
    dsa_key_import(key_file_name, p, q, a, x);
#if DEBUG > 0
    gmp_printf("p = %#Zx\nq = %#Zx\n", p, q);
#endif
	
    /* 4. Sign */
    dsa_sign_buffer(&msg, p, q, a, x, r, s, state);

    /* 5. Write signature in a file */
    FILE* sgn = fopen(signature_file_name, "w");
    gmp_fprintf(sgn, "#DSA signature:\nr = %#Zx\ns = %#Zx\n", r, s);
	
    /* . Cleaning */
    mpz_clears(p, q, a, x, r, s, NULL);
    fclose(sgn);
    buffer_clear(&msg);
}


int dsa_verify_buffer(buffer_t *msg, mpz_t p, mpz_t q, mpz_t a, mpz_t r, mpz_t s, mpz_t y)
{
    int verify = 0;
/* to be filled in */
    
    buffer_t hashaux;
    buffer_init(&hashaux, hash_length(q));
    buffer_hash(&hashaux, hash_length(q), msg);
    mpz_t hash, hash2, w, const1, u1, u2, v1, v2, v;
    mpz_inits(hash, hash2, w, const1, u1, u2, v1, v2, v, NULL);
    mpz_set_ui(const1, 1);
    int aux;
    aux=linear_equation_mod(w, s, const1, q);
    if(aux == 0) {mpz_clears(hash, hash2, w, const1, u1, u2, v1, v2, v, NULL); return verify;}
    
    mpz_import(hash, hash_length(q), 1, 1, 1, 0, hashaux.tab);
    mpz_mul(u1, hash, w);
    mpz_mod(u1, u1, q);
    mpz_mul(u2, r, w);
    mpz_mod(u2, u2, q);
    mpz_powm(v1, a, u1, p);
    mpz_powm(v2, y, u2, p);
    mpz_mul(v, v1, v2);
    mpz_mod(v, v, p);
    mpz_mod(v, v, q);
    
    verify=mpz_cmp(v, r);
    if(verify == 0) 
        verify = 1;
    else verify = 0;
    
    buffer_clear(&hashaux);
    mpz_clears(hash, hash2, w, const1, u1, u2, v1, v2, v, NULL);
    return verify;
}


void dsa_import_signature(mpz_t r, mpz_t s, const char* signature_file_name){
    FILE* sgn = fopen(signature_file_name, "r");
    while(fgetc(sgn) != '\n');
    fseek(sgn, 6, SEEK_CUR);
    gmp_fscanf(sgn, "%Zx", r);
	
    while(fgetc(sgn) != '\n');
    fseek(sgn, 6, SEEK_CUR);
    gmp_fscanf(sgn, "%Zx", s);
	
    fclose(sgn);
}


int dsa_verify(const char* file_name, const char* key_file_name,
	       const char* signature_file_name){
    // 1. INIT
    mpz_t p, q, a, y, r, s;
    buffer_t msg;
    mpz_inits(p, q, a, y, r, s, NULL);
    buffer_init(&msg, 100);
	
    // 2. Imports the message
    buffer_from_file(&msg, file_name);
	
    // 3. Parse the public key 
    dsa_key_import(key_file_name, p, q, a, y);

#if DEBUG > 0
    gmp_printf("\n\np = %#Zx\nq = %#Zx\n\n", p, q);
#endif
	
    // 4. Parse the signature 
    dsa_import_signature(r, s, signature_file_name);
    int verify = dsa_verify_buffer(&msg, p, q, a, r, s, y);
	
    // 5. Cleaning and return
    mpz_clears(p, q, a, y, r, s, NULL);
    buffer_clear(&msg);
    return verify;
}