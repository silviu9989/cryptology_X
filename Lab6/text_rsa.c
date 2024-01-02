#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "buffer.h"
#include "rsa.h"
#include "text_rsa.h"

#define DEBUG 0


void lengths(int *block_length, int *cipher_length, int *last_block_size, buffer_t *msg, mpz_t N)
{
    int aux=mpz_sizeinbase(N, 2);
    *block_length=aux/8-1;
    
    int aux2=*block_length;
    int aux3=msg->length;
    float aux4=(float)aux3/aux2;
    int aux5=aux3/aux2;
    if(aux4==aux5)
    {
        *cipher_length=aux5;
        *last_block_size=*block_length;
    }
    else 
    {
        *cipher_length=aux5+1;
        aux5=aux5*aux2;
        aux5=aux3-aux5;
        *last_block_size=aux5;
    }
 
}


void RSA_text_encrypt(mpz_t *cipher, int block_length,
		      int cipher_length, int last_block_size,
		      buffer_t *msg, mpz_t N, mpz_t e)
{
    // cipher is a table of mpz_t of length cipher_length.
    // Memory allocation and initialisation of the cells is
    // already done.

    // block_length denotes the size of blocks of uchar's
    // which will partition the message.
    // last_block_size denotes the size of the last block. It may
    // be 0.
    
    int i=0;
    mpz_t aux;
    mpz_inits(aux, NULL);
    uchar* msgs;
    msgs=(uchar*)malloc(msg->length*sizeof(char));
    msgs=string_from_buffer(msg);
    
    for(i=0; i<cipher_length-1; i++)
    {
        mpz_import(cipher[i], block_length*sizeof(uchar), 1, 1, 0, 0, &msgs[i*block_length]);
    }
    
    i=cipher_length-1;
    if(last_block_size == block_length || last_block_size == 0) 
        mpz_import(cipher[i], block_length*sizeof(uchar), 1, 1, 0, 0, &msgs[i*block_length]);
    else
    {
        mpz_import(cipher[i], last_block_size*sizeof(uchar), 1, 1, 0, 0, &msgs[i*block_length]);
    }
    
    for(i=0; i<cipher_length; i++)
        mpz_powm(cipher[i], cipher[i], e, N);
    free(msgs);

}



void RSA_text_decrypt(buffer_t *decrypted, mpz_t *cipher,
		      int cipher_length, int block_length,
		      int last_block_size,
		      mpz_t N, mpz_t d)
{

    // buffer decrypted is supposed to be initialised.
    buffer_reset(decrypted);
     
    int i=0, length=0;
    
    if(last_block_size != 0 || last_block_size != block_length)
        length=block_length*(cipher_length-1)+last_block_size;
    else if(last_block_size == 0)
        length=block_length*(cipher_length-1);
    else length=block_length*cipher_length;
    
    mpz_t* msg;
    msg=(mpz_t*)malloc(cipher_length*sizeof(mpz_t));
    
    for(i=0;i<cipher_length;i++)
    {
        mpz_init(msg[i]);
    }
    for(i=0;i<cipher_length;i++)
    {
        mpz_powm(msg[i], cipher[i], d, N);
    }
    
    uchar* msgs;
    msgs=(uchar*)malloc(length*sizeof(char));
    
    for(i=0; i<cipher_length; i++)
    {
        mpz_export(&msgs[i*block_length], NULL, 1, 1, 0, 0, msg[i]);
    }

    buffer_from_string(decrypted, msgs, length);
        
    free(msg);
}