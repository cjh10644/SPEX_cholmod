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
    mpq_t x_scale,      // pending scale for x, initially 1
    int64_t *h,         // history vector for x
    const SPEX_matrix *L,// matrix L
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const mpq_t *S,     // the pending scale factor matrix
    const mpz_t *sd,    // array of scaled pivots
    const int64_t *P,   // row permutation
    const int64_t *P_inv// inverse of row permutation
)
{
    SPEX_info info;
    int sgn;
    if (!x || !h || !L || !Ldiag || !S || !P || !P_inv || !sd)
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
        SPEX_CHECK(spex_ipge(x, x_scale, h, NULL, L->v[i], P,
            P_inv, sd, SPEX_2D(S, 1, i), SPEX_2D(S, 3, i), Ldiag[i], i));
    } 
    // apply x_scale to x[P[n]] and set x_scale to 1
    SPEX_CHECK(SPEX_mpz_divexact(x->x[P[n-1]],
                            x->x[P[n-1]], SPEX_MPQ_DEN(x_scale)));
    SPEX_CHECK(SPEX_mpz_mul(x->x[P[n-1]],
                            x->x[P[n-1]], SPEX_MPQ_NUM(x_scale)));
    SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));

    return SPEX_OK; 
}
