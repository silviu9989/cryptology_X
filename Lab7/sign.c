#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gmp.h"

#include "buffer.h"
#include "sha3.h"
#include "rsa.h"
#include "sign.h"


#define DEBUG 0


void RSA_generate_key_files(const char *pk_file_name,
			    const char *sk_file_name,
			    size_t nbits, int sec, gmp_randstate_t state){

    mpz_t N, p, q, e, d;
    mpz_inits(N, p, q, e, d, NULL);
    FILE *f1;
    f1=fopen(pk_file_name, "w");
    FILE *f2;
    f2=fopen(sk_file_name, "w");
    RSA_generate_key(N, p, q, e, d, nbits, sec, state);
    
    gmp_fprintf(f1,"#RSA Public key (%i bits):\nN = %#Zx\ne = %#Zx", nbits, N, e);
    gmp_fprintf(f2,"#RSA Secret key (%i bits):\nN = %#Zx\nd = %#Zx", nbits, N, d);
    
    mpz_clears(N, p, q, e, d, NULL);
    fclose(f1);
    fclose(f2);
}


void RSA_key_import(mpz_t N, mpz_t ed, const char *key_file_name){
    FILE *key = fopen(key_file_name, "r");
    /* Go to second line, then move from 6 characters to the right */
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);

    /* Scan the modulus N */
    gmp_fscanf(key, "%Zx", N);

    /* Same for e or d*/
    while(fgetc(key) != '\n');
    fseek(key, 6, SEEK_CUR);
    gmp_fscanf(key, "%Zx", ed);

    fclose(key);
}


int hash_length(mpz_t N){
    int bit_size_N = mpz_sizeinbase(N, 2);
    return (bit_size_N % BYTE_SIZE == 0) ?
	bit_size_N / BYTE_SIZE - 1 : (bit_size_N / BYTE_SIZE);
}


void RSA_sign_buffer(mpz_t sgn, buffer_t *msg, mpz_t N, mpz_t d)
{
    mpz_t hash;
    mpz_inits(hash, NULL);

    buffer_t hashaux;
    buffer_init(&hashaux, hash_length(N));
    
    buffer_hash(&hashaux, hash_length(N), msg);
    mpz_import(hash, hash_length(N), 1, 1, 1, 0, hashaux.tab);

    mpz_powm(sgn, hash, d, N);
    buffer_clear(&hashaux);
    mpz_clears(hash, NULL);
}


int RSA_verify_signature(mpz_t sgn, buffer_t *msg, mpz_t N, mpz_t e)
{
    int verify = 0;
    
    buffer_t hashaux;
    buffer_init(&hashaux, hash_length(N));
    buffer_hash(&hashaux, hash_length(N), msg);
    
    mpz_t hash, hash2;
    mpz_inits(hash, hash2, NULL);
    
    mpz_powm(hash, sgn, e, N);
    
    mpz_import(hash2, hash_length(N), 1, 1, 1, 0, hashaux.tab);
    
    verify=mpz_cmp(hash, hash2);
    if(verify == 0) 
        verify = 1;
    else verify = 0;
    
    buffer_clear(&hashaux);
    mpz_clears(hash, hash2, NULL);
    return verify;
}


void RSA_signature_import(mpz_t S, const char* signature_file_name){
    FILE *sgn = fopen(signature_file_name, "r");
    while(fgetc(sgn) != '\n');
    fseek(sgn, 6, SEEK_CUR);
    gmp_fscanf(sgn, "%Zx", S);
    fclose(sgn);
}


void RSA_sign(const char* file_name, const char* key_file_name,
	      const char* signature_file_name){
    // 1. Initialisation
    buffer_t msg;
    mpz_t N, d, signature;
    mpz_inits(N, d, signature, NULL);
    buffer_init(&msg, 100);

    // 2. Import the message in a buffer
    buffer_from_file(&msg, file_name);
	
    // 3. Parse the secret key
    RSA_key_import(N, d, key_file_name);

    // 4. Sign the buffer
    RSA_sign_buffer(signature, &msg, N, d);

    // 5. Exports the signature in a file
    FILE* sgn = fopen(signature_file_name, "w");
    gmp_fprintf(sgn, "#RSA signature\nS = %#Zx\n", signature);
	
    // 6. Close and free
    fclose(sgn);
    mpz_clears(N, d, signature, NULL);
    buffer_clear(&msg);
}


int RSA_verify(const char* file_name, const char* key_file_name,
	       const char* signature_file_name){
    // 1. Initialisation
    buffer_t msg;
    mpz_t N, e, S;	
    buffer_init(&msg, 100);
    mpz_inits(N, e, S, NULL);	

	
    // 2. Import the message into a buffer
    buffer_from_file(&msg, file_name);

    // 3. Import the public key
    RSA_key_import(N, e, key_file_name);

    // 4. Parse the signature
    RSA_signature_import(S, signature_file_name);
	
    // 5. Verify
    int verify = RSA_verify_signature(S, &msg, N, e);
	
    // 6. Close, free and return
    mpz_clears(S, N, e, NULL);
    buffer_clear(&msg);
    return verify;
}