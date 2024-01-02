#include <stdio.h>
#include <assert.h>

#include "gmp.h"

#include "xgcd.h"

#define DEBUG 0


/* to be filled in */

/* compute g, u and v s.t. a*u+b*v = g = gcd(a, b) */
void XGCD_long(long *g, long *u, long *v, long a, long b){
/* to be filled in */

	long u1,v1,u2,v2,r1,r2,q,raux,uaux,vaux;
	u1=1; v1=0; u2=0; v2=1; r1=a; r2=b;
	while(r2!=0)
	{
		q=(int)r1/r2;
		raux=r1; uaux=u1; vaux=v1;
		r1=r2; u1=u2; v1=v2;
		r2=raux-q*r2;
		u2=uaux-q*u2; 
		v2=vaux-q*v2;
	}
	
	*g=r1;
	*u=u1;
	*v=v1;	
}

/* compute g, u and v s.t. a*u+b*v = g = gcd(a, b) */
void XGCD(mpz_t g, mpz_t u, mpz_t v, mpz_t a, mpz_t b)
{
	mpz_t u1,v1,u2,v2,r1,r2,q,raux,uaux,vaux,const1,const0,aux;
	mpz_inits(u1,v1,u2,v2,r1,r2,q,raux,uaux,vaux,const1,const0,aux,NULL);
	mpz_add_ui(const1,const0,1);
	mpz_add_ui(u1,const0,1);
	mpz_add_ui(v2,const0,1);
	mpz_set(r1, a);
	mpz_set(r2, b);
	
	while(mpz_cmp(r2,const0)!=0)
	{
		mpz_tdiv_q(q,r1,r2);
		mpz_set(raux, r1); mpz_set(uaux, u1); mpz_set(vaux, v1);
		mpz_set(r1, r2); mpz_set(u1, u2); mpz_set(v1, v2);
		mpz_mul(aux,q,r2); mpz_sub(r2,raux,aux);
		mpz_mul(aux,q,u2); mpz_sub(u2,uaux,aux);
		mpz_mul(aux,q,v2); mpz_sub(v2,vaux,aux);
	}
	
	mpz_set(g, r1);
	mpz_set(u, u1);
	mpz_set(v, v1);
	mpz_clears(u1,v1,u2,v2,r1,r2,q,raux,uaux,vaux,const1,const0,aux,NULL);
/* to be filled in */
}


/* Solve a*x=b mod m if possible. 
   
 */
int linear_equation_mod(mpz_t x, mpz_t a, mpz_t b, mpz_t m){
    int status = 1;
/* to be filled in */

    
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
 
