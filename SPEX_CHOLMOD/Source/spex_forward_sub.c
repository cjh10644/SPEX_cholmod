//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_forward_sub.c: sparse forward substitution (x = (LD)\v)
//------------------------------------------------------------------------------
    
// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------
    
// Purpose: This function is to perform sparse forward substitution, which is
// essentially the same as sparse REF triangular solve for LDx=v, but with v
// as a dense vector. This function assumes v is in the same row permutation as
// L. This function take v as input using x_input, and the solution is stored
// in x_output. In case when v has multiple column, simply iterate this function
// and solve each column at each iteration.

#include "spex_internal.h"

SPEX_info spex_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    int64_t *h,         // history vector for x
    const SPEX_matrix *L,// matrix L
    const SPEX_matrix *U,// matrix U
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Ucp, // col pointers for col-wise nnz pattern of U
    const int64_t *Ucx, // the value of k-th entry is found as 
                        // U->v[Uci[k]]->x[Ucx[k]]
    mpq_t *S,           // the pending scale factor matrix
    const mpz_t *sd,    // array of scaled pivots
    mpz_t *d,           // array of unscaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t *Q    // column permutation
)
{
    SPEX_info info;
    int sgn;
    if (!x || !h || !L || !U || !Ldiag || !Ucp || !Ucx ||
        !S || !P || !P_inv || !Q || !sd)
    {
        return SPEX_INCORRECT_INPUT;
    }

    int64_t i, n = L->n;
    for (i = 0; i < n; i++)
    {
        h[i] = -1; 
    }

    for (i = 0; i < n; i++)
    {
        // skip if x(P[i]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x[P[i]]));
        if (sgn == 0)       { continue; }

        // perform i-th IPGE update for x
        SPEX_CHECK(spex_ipge(x, h, NULL, L->v[i], P, P_inv, sd,
            d, U->v[i]->x[Ucx[Ucp[Q[i]+1]-1]],
            SPEX_2D(S, 1, i), SPEX_2D(S, 3, i), SPEX_2D(S, 2, i), i, Ldiag[i]));
    }

    return SPEX_OK; 
}
