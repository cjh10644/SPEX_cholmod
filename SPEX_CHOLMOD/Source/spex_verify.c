//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_verify.c: verify if A=LD^(-1)U
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to verify if A=L(P,:)D^(-1)U(:,Q)

#define SPEX_FREE_ALL                 \
    SPEX_matrix_free(&b);             \
    SPEX_matrix_free(&x);             \
    SPEX_matrix_free(&b2);            \
    SPEX_MPQ_CLEAR(x_scale);

#include "spex_internal.h"

SPEX_info spex_verify
(
    bool *correct,         // indicate if the verification is passed
    const SPEX_matrix *L,  // lower triangular matrix
    const SPEX_matrix *U,  // upper triangular matrix
    const SPEX_matrix *A,  // Input matrix
    int64_t *h,            // history vector
    const mpz_t *sd,       // array of scaled pivots
    mpz_t *d,              // array of unscaled pivots
    mpq_t *S,              // the pending scale factor matrix
    const int64_t *P,      // row permutation
    const int64_t *P_inv,  // inverse of row permutation
    const int64_t *Q,      // column permutation
    const int64_t *Q_inv,  // inverse of column permutation
    const int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Ucp,    // col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,    // the value of k-th entry is found as 
                           // U->v[Uci[k]]->x[Ucx[k]]
    const SPEX_options *option// command options
)
{
    SPEX_info info;
    int64_t tmp, i, n = L->n;
    int sgn;
    mpq_t x_scale; SPEX_MPQ_SET_NULL(x_scale);
    SPEX_matrix *b = NULL; // the dense right-hand-side matrix to be generated
    SPEX_matrix *x = NULL; // the dense solution matrix to be generated
    SPEX_matrix *b2 = NULL; // the dense matrix to store the result of A*x

    SPEX_CHECK(SPEX_mpq_init(x_scale));
    SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));
    SPEX_CHECK(SPEX_matrix_alloc(&b , 1, n, false));
    SPEX_CHECK(SPEX_matrix_alloc(&b2, 1, n, false));

    // -------------------------------------------------------------------------
    // generate random right-hand-size vector
    // -------------------------------------------------------------------------
    // initialize random number generator
    int seed = 10;
    srand(seed);
    for (i = 0; i < n; i++)
    {
        tmp = i+1;//rand(); //TODO
        SPEX_CHECK(SPEX_mpz_set_si(b->v[0]->x[i], tmp));
    }

    // -------------------------------------------------------------------------
    // solve LD^(-1)Ux = b for x
    // -------------------------------------------------------------------------
    SPEX_CHECK(SPEX_solve(&x, b, h, A, L, U, S, sd, d, Ldiag, Ucp, Ucx,
        P, P_inv, Q, Q_inv, true, option));

    // -------------------------------------------------------------------------
    // compute b2 = A*x
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->v[0]->x[i]));
        if (sgn == 0) { continue;}

        for (int64_t p = 0; p < A->v[i]->nz; p++)
        {
            int64_t j = A->v[i]->i[p];
            // b2[j] += x[i]*A(j,i)
            SPEX_CHECK(SPEX_mpz_addmul(b2->v[0]->x[j],
                                       x->v[0]->x[i], A->v[i]->x[p]));
        }
    }
    // update b2->scale = x->scale*A->scale
    SPEX_CHECK(SPEX_mpq_mul(b2->scale, x->scale, A->scale));

    // -------------------------------------------------------------------------
    // check if b2 == b
    // -------------------------------------------------------------------------
    *correct = true;
    // set b2->scale = b2->scale/b->scale since we only want to compare the
    // integer values b2*b->scale and b*b->scale. b*b->scale are the values
    // stored in b->v->x. It can be shown that the resulted b2->scale = sd[n-1]
    SPEX_CHECK(SPEX_mpq_div(b2->scale, b2->scale, b->scale));
#ifdef SPEX_DEBUG
    SPEX_CHECK(SPEX_mpq_cmp_z(&sgn, b2->scale, sd[n-1]));
    ASSERT(sgn == 0);
#endif
    for (i = 0; i < n; i++)
    {
        // set b2[i] = b2[i]/sd[n-1] = b2[i]/b2->scale
        // This division will be exact since b are integers
        SPEX_CHECK(SPEX_mpz_divexact(b2->v[0]->x[i],
                                     b2->v[0]->x[i], SPEX_MPQ_NUM(b2->scale)));

        SPEX_CHECK(SPEX_mpz_cmp(&sgn, b2->v[0]->x[i], b->v[0]->x[i]));
        if (sgn != 0)
        {
            *correct = false;
            break;
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
