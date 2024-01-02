#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"



void fermat(int a0, int pmin, int pmax, int composites){
/* to be filled in */

	mpz_t a0m, pminm, pmaxm, const2, const1, const0, im, aux, im_aux, r;
	mpz_init(a0m);mpz_init(r);mpz_init(pminm);mpz_init(pmaxm);mpz_init(const2);mpz_init(const1);mpz_init(const0);mpz_init(im);mpz_init(aux);mpz_init(im_aux);
	
	mpz_init_set_str(const2, "2", 10);
	mpz_init_set_str(const0, "0", 10);
	mpz_init_set_str(const1, "1", 10);
	mpz_set_ui(a0m, a0);
	mpz_set_ui(pminm,pmin);
	mpz_set_ui(pmaxm,pmax);
	mpz_set(im, pminm);
	mpz_mod(r,im,const2);
	if(mpz_cmp(r,const0)==0) mpz_add(im,im,const1);
	if(composites==0)
	{
		while(mpz_cmp(im,pmaxm)<0)
		{
			mpz_sub(im_aux,im,const1);
			mpz_powm(aux, a0m, im_aux, im);
			if(mpz_cmp(aux,const1)==0)
				gmp_printf("%Zd %Zd\n",a0m, im);
			mpz_add(im,im,const2);
		}
	}
	if(composites==1)
	{
		while(mpz_cmp(im,pmaxm)<0)
		{
			mpz_sub(im_aux,im,const1);
			mpz_powm(aux, a0m, im_aux, im);
			if(mpz_cmp(aux,const1)==0 && mpz_probab_prime_p(im,50)==0)
				gmp_printf("%Zd %Zd\n",a0m, im);
			mpz_add(im,im,const2);
		}
	}
	mpz_clears(a0m, pminm, pmaxm, const2, const1, const0, im, aux, im_aux, r, NULL);
}


void Usage(char *cmd){
    fprintf(stderr, "Usage: %s a pmin pmax [0|1]\n", cmd);
}


int main (int argc, char *argv[]){
    if(argc < 4){
	Usage(argv[0]);
	return 0;
    }
    int composites = 0;
    if(argc == 5)
	composites = atoi(argv[4]);
    fermat(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), composites);
    return 0;
}
