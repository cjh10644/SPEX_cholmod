//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Tcov/tcov_test.c: test coverage for SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/*
 * test coverage
 */

#define SPEX_FREE_ALL                            \
{                                                \
    SPEX_FREE(option);                           \
    SPEX_matrix_free(&L);                        \
    SPEX_matrix_free(&U);                        \
    SPEX_matrix_free(&A);                        \
    SPEX_vector_free(&vk);                       \
    spex_delete_mpz_array(&d, n);                \
    spex_delete_mpz_array(&sd, n);               \
    spex_delete_mpq_array(&S, 3*n);              \
    SPEX_FREE(P);                                \
    SPEX_FREE(P_inv);                            \
    SPEX_FREE(Q);                                \
    SPEX_FREE(Q_inv);                            \
    SPEX_MPZ_CLEAR(tmpz);                        \
    SPEX_MPQ_CLEAR(tmpq);                        \
    SPEX_finalize() ;                            \
}


#include "tcov_malloc_test.h"

#define TEST_CHECK(method)                       \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_OK)                         \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        SPEX_FREE_ALL;printf("n=%ld\n",n);                           \
        return 0 ;/*continue; */                               \
    }                                            \
}

#define TEST_CHECK_FAILURE(method)               \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_INCORRECT_INPUT && info != SPEX_SINGULAR) \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        SPEX_FREE_ALL ;                          \
        continue ;                               \
    }                                            \
    else                                         \
    {                                            \
        printf("Expected failure at line %d\n", __LINE__);\
    }                                            \
}

#define MAX_MALLOC_COUNT 1000

int64_t Aold[16] = {3, 5, 6, 7, 11, 2, -7, 10, 8, 3, -2, -2, 7, 5, 1, -6};//CSC
int64_t Lold[10] = {3, 5, 6, 7, -49, -87, -47, -17, 527, 884};// CSC
int64_t Uold[10] = {3, 11, 8, 7, -49, -31, -20, -17, 57,884};//CSR
int64_t Ak_new[4]= {1,4,7,11}; // new column
#include <assert.h>

int main( int argc, char* argv[])
{
    SPEX_info info;
    int64_t n = 4, n_max, i, j, p;//TODO debug n=90 10
    GOTCHA;

    if (argc >= 2)
    {
        n = atoi(argv[1]);
    }
    n_max = n+1;
    if (argc >= 3)
    {
        n_max = atoi(argv[2]);
    }

    for (; n < n_max;n++)
    {
    for (int64_t loop = 0; loop < 100; loop++)
    {
        printf("\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n+++++++++++++++++++++++++++++++++case %ld+++++++++++++++++++++++++++++++++\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n",loop);
        //------------------------------------------------------------------
        // Initialize SPEX CHOLMOD process
        //------------------------------------------------------------------

        SPEX_initialize_expert (tcov_malloc, tcov_calloc,
            tcov_realloc, tcov_free) ;

        info = SPEX_initialize ( ) ;
        assert (info == SPEX_PANIC) ;
    GOTCHA;

        //------------------------------------------------------------------
        // Allocate memory
        //------------------------------------------------------------------

        SPEX_options* option = SPEX_create_default_options();
        if (!option) return 0;//{continue;}

        SPEX_matrix *L = NULL, *U = NULL, *A = NULL;
        SPEX_vector *vk = NULL;
        mpz_t *d = NULL, *sd = NULL;
        mpq_t *S = NULL;
        int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
        mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
        mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
    GOTCHA;

        TEST_CHECK(SPEX_matrix_alloc(&L, n, n, true));
        printf("l0=[\n");
        for (j = 0; j < n; j++)
        {
            TEST_CHECK(SPEX_vector_realloc(L->v[j], n-j));
            TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[0], 1)); // diagnal entry
            L->v[j]->i[0] = j;
            p = 1;
            for (i = 0; i < j; i++){printf("0 ");}
                    printf("1 ");
            for (i = j+1; i < n; i++)
            {
                if (rand() > RAND_MAX/2)
                {
                    printf("1 ");
                    TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[p], 1));
                    L->v[j]->i[p] = i;
                    p++;
                }
                else
                {
                    printf("0 ");
                }
            }
            printf("%s;%%   --> %ld(%ld) \n",j==n-1?"]'":"",j,p);
            L->v[j]->nz = p;
        }
        printf("\nu0=[\n");
        TEST_CHECK(SPEX_matrix_alloc(&U, n, n, true));
        for (j = 0; j < n; j++)
        {
            TEST_CHECK(SPEX_vector_realloc(U->v[j], n-j));
            TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[0], 1)); // diagnal entry
            U->v[j]->i[0] = j;
            p = 1;
            for (i = 0; i < j; i++){printf("0 ");}
                    printf("1 ");
            for (i = j+1; i < n; i++)
            {
                if (rand() > RAND_MAX/2)
                {
                    printf("1 ");
                    TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[p], 1));
                    U->v[j]->i[p] = i;
                    p++;
                }
                else
                {
                    printf("0 ");
                }
            }
            printf("%s;%%   --> %ld(%ld) \n",j==n-1?"]":"",j,p);
            U->v[j]->nz = p;
        }
        TEST_CHECK(SPEX_matrix_alloc(&A, n, n, true));
        for (i = 0; i < n; i++)
        {
            TEST_CHECK(SPEX_vector_realloc(A->v[i], n));
            for (p = 0; p < n; p++)
            {
                TEST_CHECK(SPEX_mpz_set_si(A->v[i]->x[p], 0));
                A->v[i]->i[p] = p;
            }
            A->v[i]->nz = n;
        }
        for (int64_t iter = 0; iter < n; iter++)
        {
            for (int64_t Up = 0; Up < U->v[iter]->nz; Up++)
            {
                j = U->v[iter]->i[Up];
                for (int64_t Lp = 0; Lp < L->v[iter]->nz; Lp++)
                {
                    i = L->v[iter]->i[Lp];
                    TEST_CHECK(SPEX_mpz_add(A->v[j]->x[i],
                                            A->v[j]->x[i], L->v[0]->x[0]));
                }
            }
        }

        int64_t k = rand()%n;
        printf("A0=l0*u0;\nk=%ld;%%1-based, for matlab code \nvk=[",k+1);
        TEST_CHECK(SPEX_vector_alloc(&vk, n, true));
        p = 0;
        for (i = 0; i < n; i++)
        {
            if (i == k || rand() > RAND_MAX/2)
            {
                TEST_CHECK(SPEX_mpz_set_si(vk->x[p], 1));
                vk->i[p] = i;
                p++;
                printf("1 ");
            }
            else
            {
                printf("0 ");
            }
        }
        printf("]';\n\n\n");
        vk->nz = p;

        d  = spex_create_mpz_array(n);
        sd = spex_create_mpz_array(n);
        S  = spex_create_mpq_array(3*n);
        P     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
        Q     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
        P_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
        Q_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
        if (!d || !sd || !S || !P || !Q || !P_inv || !Q_inv)
        {
            SPEX_FREE_ALL;
            return 0;
        }
        for (i = 0; i < n; i++)
        {
            TEST_CHECK(SPEX_mpz_set_ui(d[i],  1));
            TEST_CHECK(SPEX_mpz_set_ui(sd[i], 1));
            TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, i), 1, 1));
            TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
            TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));
            P[i] = i;
            Q[i] = i;
            P_inv[i] = i;
            Q_inv[i] = i;
        }
        GOTCHA;
        info=SPEX_LUU(A, L, U, d, sd, S, P, P_inv, Q, Q_inv, &vk, k, NULL);
        if(info!=SPEX_SINGULAR && info!=SPEX_OK){TEST_CHECK(info);}

        /*TEST_CHECK(SPEX_mpz_init(tmpz));
        TEST_CHECK(SPEX_mpq_init(tmpq));
        for (i = 0; i < n; i++)
        {
            TEST_CHECK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 1, i), SPEX_2D(S, 3, i)));

            for (p = 0; p < L->v[i]->nz; p++)
            {
                j = L->v[i]->i[p];
                TEST_CHECK(SPEX_mpz_divexact(tmpz, L->v[i]->x[p],
                    SPEX_MPQ_DEN(tmpq)));
                TEST_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                TEST_CHECK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
            }
            TEST_CHECK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
        }
        for (i = 0; i < n; i++)
        {
            TEST_CHECK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 2, i), SPEX_2D(S, 3, i)));

            for (p = 0; p < U->v[i]->nz; p++)
            {
                j = U->v[i]->i[p];
                TEST_CHECK(SPEX_mpz_divexact(tmpz, U->v[i]->x[p],
                    SPEX_MPQ_DEN(tmpq)));
                TEST_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                TEST_CHECK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
            }
            TEST_CHECK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
        }
        printf("P=[%ld %ld %ld %ld]\n",P[0], P[1], P[2], P[3]);
        printf("Q=[%ld %ld %ld %ld]\n",Q[0], Q[1], Q[2], Q[3]);*/
        SPEX_FREE_ALL;
    }
    }

    return 0;
}

