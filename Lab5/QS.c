/****************************************************************/
/* QS.c                                                         */
/* Author : F. Morain                                           */
/* morainr@lix.polytechnique.fr                                 */
/* Last modification October 24, 2017                           */
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "gmp.h"
#include "hash.h"
#include "utils.h"
#include "QS.h"

#define DEBUG 0

int trial_div(char *tabex, mpz_t cof, const mpz_t Px, int* B, int cardB){
    int status = FACTOR_NOT_FOUND;
    int i=0, counter=0;
    mpz_t aux;
    mpz_inits(aux, NULL);
    mpz_set(cof, Px);
    if(mpz_cmp_si(cof, 0)<0)
        {
            tabex[0]=1;
            mpz_mul_si(cof, cof, -1);
        }
    else tabex[0]=0;
        
    for(i=1;i<cardB;i++)
    {
        tabex[i]=0;
        mpz_mod_ui(aux, cof, B[i]);
        while(mpz_cmp_si(aux, 0) == 0)
        {
            counter++;
            mpz_divexact_ui(cof, cof, B[i]);
            mpz_mod_ui(aux, cof, B[i]);
        }
        if(counter!=0) 
        {
            tabex[i]=counter;
            counter=0;
        }
    }
    if(mpz_cmp_ui(cof, 1)!=0)
    {mpz_clears(aux, NULL);status=0;}
    else {status = 1;mpz_clears(aux, NULL);}
    return status;
}

void StoreRelation(relation_t *rel, mpz_t kN, mpz_t g, 
		   int cardB, int x, char *tabex){
    char *tmp = (char *)malloc(cardB * sizeof(char));

    memcpy(tmp, tabex, cardB);
    mpz_init_set_si(rel->y, x);
    mpz_add(rel->y, rel->y, g);
    mpz_mod(rel->y, rel->y, kN);
    rel->tabex = tmp;
}

void AddRelation(relation_t *tabrels, mpz_t kN, mpz_t g, 
		 int cardB, int i, int x, char *tabex){
    StoreRelation(tabrels+i, kN, g, cardB, x, tabex);
#if 1
    printf("%d:", i);
    for(i = 0; i < cardB; i++)
	printf(" %d", tabex[i]);
    printf("\n");
#endif
}

/* OUTPUT: the actual number of relations found <= nrelsmax */
int FindRelationsUsingTrialDivision(relation_t *tabrels, mpz_t kN, mpz_t g,
				    int *B, int cardB, int M, int nrelsmax){
    char *tabex = (char *)malloc(cardB * sizeof(char));
    int x, nrels = 0;
    mpz_t Px, cof;

    mpz_inits(Px, cof, NULL);
    for(x = -M; x <= M; x++){
	mpz_set_si(Px, x);
	mpz_add(Px, Px, g);
	mpz_mul(Px, Px, Px);
	mpz_sub(Px, Px, kN);
	if(trial_div(tabex, cof, Px, B, cardB) != 0){
	    gmp_printf("x=%d Px=%Zd\n", x, Px);
	    AddRelation(tabrels, kN, g, cardB, nrels, x, tabex);
	    nrels++;
	    if(nrels == nrelsmax)
		break;
	}
    }
    free(tabex);
    mpz_clears(Px, cof, NULL);
    return nrels;
}

/* to be filled in */

/* Sieving over [-M, M].
   OUTPUT: the actual number of relations found <= nrelsmax */
int FindRelationsUsingSieving(relation_t *tabrels, mpz_t kN, mpz_t g,
			      int *B, int cardB, int lpB, int M, int nrelsmax){
    int nrels = 0;
/* to be filled in */
    return nrels;
}

/* OUTPUT: the actual number of relations found <= nrelsmax */
int FindRelations(relation_t *tabrels, mpz_t kN, mpz_t g, int *B, 
		  int cardB, int M, int nrelsmax){
#if 1
    return FindRelationsUsingTrialDivision(tabrels, kN, g, B, cardB, M, nrelsmax);
#else
    int lpB = (mpz_cmp_ui(kN, 2000) <= 0 ? 0 : 100 * B[cardB-1]);
    return FindRelationsUsingSieving(tabrels,kN,g,B,cardB,lpB,M,nrelsmax);
#endif
}

void PrintMatrix(char **mat, int nrows, int ncols){
    int i, j;

    for(i = 0; i < nrows; i++){
	for(j = 0; j < ncols; j++)
	    printf(" %d", mat[i][j]);
	printf("\n");
    }
}

void PrintMatrices(char **mat, char **C, int nrows, int ncols){
    int i=0, j=0;
    for(i=0;i<nrows;i++)
    {
        for(j=0;j<ncols;j++)
        {
            printf("%i ", mat[i][j]);
        }
        for(j=0;j<=i;j++)
            printf("%i ", C[i][j]);
        printf("\n");
    }
}

/* mat[i1] += mat[i2]; C[i1] += C[i2]. */
void AddRows(char **mat, char **C, int ncols, int i1, int i2, int j){
    
    int i=0;
    for(i=0;i<ncols;i++)
        mat[i1][i]=(mat[i1][i]+mat[i2][i])%2;
    for(i=0; i<=i1 && i<=i2; i++)
        C[i1][i]=(C[i1][i]+C[i2][i])%2;
}

void Gauss(char **mat, char **C, int nrows, int ncols){
    int i=0, j=0, k=0, counter=0, pivot=0, flag=0;
    for(k=0;k<ncols;k++)
    {
        flag=0;
        for(i=k;i<nrows-1;i++)
        {
            if(mat[i][k]==1)
            {
                for(j=i+1;j<nrows;j++)
                {
                    if(mat[j][k]==1)
                    {
                        if(flag == 0)
                        {
                            pivot=i;
                            flag=1;
                        }
                        AddRows(mat, C, ncols, j, i, j);
                    }
                }
            }
        }
        if(flag == 0) {pivot=7;}
      printf("pivot[%i]=%i\n",counter, pivot);
      counter++;
      PrintMatrices(mat, C, nrows, ncols);  
    }
}

char **MatrixFromRelations(relation_t *tabrels, int nrows, int ncols){
    char **mat = (char **)malloc(nrows * sizeof(char *));
    int i,j;
    for(i=0;i<nrows;i++)
    {
        mat[i]=(char*)malloc(ncols*sizeof(char));
        for(j=0;j<ncols;j++)
        {
            mat[i][j]=tabrels[i].tabex[j]%2;
        }
    }
    return mat;
}

char **BuildCompanionMatrix(int nrows){
    char **C = (char **)malloc(nrows * sizeof(char *));
    int i=0;
    for(i=0;i<nrows;i++)
        {
        C[i]=(char*)calloc((i+1), sizeof(char));
        C[i][i]=1;
        }
    return C;
}

int FinishFactorization(factor_t *tabf, int *nf, mpz_t N, mpz_t kN, mpz_t g, 
			relation_t *tabrels, char **mat, char **C,
			int nrelsmax, int *B, int cardB){
    int status = FACTOR_NOT_FOUND;
    int i=0,j=0,flag=0,k=0;
    int aux=0;
    mpz_t X, Y, auxX, auxX2, auxY, comdiv;
    mpz_inits(X, Y, auxX, auxX2, auxY, comdiv, NULL);
    mpz_set_ui(Y, 1);
    mpz_set_ui(comdiv, 1);
    mpz_set_ui(X, 1);
    mpz_set_ui(auxX, 1);
    mpz_set_ui(auxX2, 1);
    char *vector;
    vector=(char*)calloc(cardB,sizeof(char));
    for(i=0;i<nrelsmax;i++)
    {
        flag=0;
        for(j=0;j<cardB;j++)
        {
            if(mat[i][j]==1) 
                flag=1;
        }
        if(flag == 0)
        {
            printf("dep: ");
            for(j=0;j<=i;j++)
            {
                if(C[i][j]==1) 
                {
                    status=1;
                    for(k=0;k<cardB;k++)
                    {
                        vector[k]=vector[k]+tabrels[j].tabex[k];
                    }
                    printf("%i ", j);
                }
            }
            printf("-> ");
            mpz_set_ui(Y, 1);
            mpz_set_ui(X, 1);
            mpz_set_ui(auxX, 1);
            mpz_set_ui(auxX2, 1);
            for(k=0;k<cardB;k++) 
            {
                printf("%i ", vector[k]);
                aux=pow(B[k],vector[k]);
                mpz_mul_ui(Y, Y, aux);
            }
            mpz_set(auxY, Y);
            mpz_sqrt(auxY, auxY);
            mpz_mod(auxY, auxY, kN);
            while(1)
            {
                mpz_mul(auxX2, auxX, auxX);
                mpz_sub(auxX2, auxX2, Y);
                mpz_mod(auxX2, auxX2, kN);
                if(mpz_cmp_ui(auxX2,0)==0 && mpz_cmp(auxX,auxY)!=0)
                {break;}
                mpz_add_ui(auxX, auxX, 1);
            }
            mpz_set(X, auxX);
            mpz_sqrt(Y, Y);
            mpz_mod(Y, Y, kN);
            mpz_sub(auxX, X, Y);
            mpz_gcd(comdiv, auxX, N);
            if(mpz_cmp_ui(comdiv,1)!=0) gmp_printf("-> X=%Zd, Y=%Zd -> %Zd\n", X, Y, comdiv);
            else
            {
            mpz_add(auxX, X, Y);
            mpz_gcd(comdiv, auxX, N);
            if(mpz_cmp_ui(comdiv,1)!=0) gmp_printf("-> X=%Zd, Y=%Zd -> %Zd\n", X, Y, comdiv);
            }
            mpz_set_ui(Y, 1);
        }
        for(k=0;k<cardB;k++)
        {
            vector[k]=0;
        }
    }
    mpz_clears(X, Y, auxX, auxX2, auxY, comdiv, NULL);
    return status;
}

int QS_aux(factor_t *tabf, int *nf, mpz_t N, mpz_t kN, mpz_t g, int *B,
	   int cardB, int M, int phase){
    int nrelsmax = cardB+2, nrels, i, status = FACTOR_NOT_FOUND;
    relation_t *tabrels = (relation_t *)malloc(nrelsmax * sizeof(relation_t));
    char **mat, **C;

    nrels = FindRelations(tabrels, kN, g, B, cardB, M, nrelsmax);
    if(nrels < cardB){
	printf("Not enough relations: %d // %d\n", nrels, nrelsmax);
	return -1;
    }
    if(phase == 1)
	return 0;
    mat = MatrixFromRelations(tabrels, nrels, cardB);
#if 1
    PrintMatrix(mat, nrels, cardB);
#endif
    if(phase == 2)
	return 0;
    C = BuildCompanionMatrix(nrels);
#if 1
    PrintMatrices(mat, C, nrels, cardB);
#endif
    if(phase == 3)
	return 0;
    Gauss(mat, C, nrels, cardB);
    if(phase == 4)
	return 0;
    status = FinishFactorization(tabf, nf, N, kN, g, tabrels, mat, C, nrels, B, cardB);
    free(tabrels);
    for(i = 0; i < nrels; i++){
	free(mat[i]);
	free(C[i]);
    }
    free(mat);
    free(C);
    return status;
}

/* Source: Silverman87. */
int FindMultiplier(mpz_t N){
    int kopt = -1;
/* to be filled in */
    return kopt;
}

/* OUTPUT: NULL if some problem occurred, a factor base otherwise of size
   cardB, starting {-1, 2, ...}. */
int *BuildFactorBase(mpz_t kN, int k, int cardB, FILE *file){
    int *B = NULL;
/* to be filled in */
    return B;
}

/* This is from Silverman87, but for MPQS. */
void ChooseParameters(int *cardB, int *M, mpz_t N){
    size_t dd = mpz_sizeinbase(N, 10);
    int thresh[] = {10, 24, 30, 36, 42, 48, 54, 60, 66, 0};
    int tcardB[] = {50, 100, 200, 400, 900, 1200, 2000, 3000, 4500, 0};
    int tM[] = {1000, 5000, 25000, 25000, 50000, 100000, 250000, 350000, 500000, 0};
    int i;

    *cardB = -1; *M = -1;
    for(i = 0; thresh[i] != 0; i++){
	if(dd <= thresh[i]){
	    *cardB = tcardB[i];
	    *M = tM[i];
	    break;
	}
    }
}

int QS(factor_t *tabf, int *nf, mpz_t N, int k, int cardB, int M, FILE *file,
       int phase){
    int *B, status = FACTOR_NOT_FOUND;
    mpz_t kN, g;

    mpz_inits(kN, g, NULL);
    if(k == 0){
	k = FindMultiplier(N);
	printf("Best multiplier: %d\n", k);
    }
    mpz_mul_ui(kN, N, k);
    /* g = trunc(sqrt(k*N)) */
    mpz_sqrt(g, kN);
    if(cardB == 0 || M == 0){
	int cardB0 = cardB, M0 = M;
	
	ChooseParameters(&cardB, &M, kN);
	if(cardB0 != 0)
	    cardB = cardB0;
	if(M0 != 0)
	    M = M0;
    }
    B = BuildFactorBase(kN, k, cardB, file);
    if(B == NULL)
	return FACTOR_ERROR;
    printf("cardB=%d, M=%d\n", cardB, M);
    status = QS_aux(tabf, nf, N, kN, g, B, cardB, M, phase);
    
    mpz_clears(kN, g, NULL);
    free(B);
    return status;
}