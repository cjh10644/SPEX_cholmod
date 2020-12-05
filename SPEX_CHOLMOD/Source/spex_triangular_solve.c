//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_triangular_solve.c: perform REF triangular solve up to
// specified iteration, additional history update for certain entries should
// be done after calling this function.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform REF triangular solve for LDx=v up to
// specified IPGE iteration when both L and v are sparse. Additional history
// update should be done for certain entries based on the history vector. This
// function is used internally when computing the inserted column for U. For
// forward solving LDUx=b where b is mostly considered and stored as dense
// vector, we will use spex_forward_sub.

#include "spex_internal.h"

SPEX_info spex_triangular_solve // perform REF triangular solve for LDx=v
(
    spex_scattered_vector *sv_x,// the scattered version of solution for LDx=v,
                        // using the first k-1 columns of L
    mpq_t x_scale,      // pending scale for x
    int64_t *h,         // history vector for x
    int64_t *last_update,// the number of finished IPGE iterations, which is
                        // also the number of columns in L used last time
    int64_t *i_2ndlast, // i_2ndlast is the index of the found last nnz entry
                        // of x[P] less than n, this could be NULL if not needed
    const int64_t k,    // compute x up to k-th IPGE iteration, that is, using
                        // the first k-1 columns of L
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
    if (!sv_x || !h || !last_update || !L || !S || !P || !P_inv || !sd)
    {
        return SPEX_INCORRECT_INPUT;
    }

    if (*last_update < k-1)
    {
        // TODO iterate sorted nnz pattern?
        for (int64_t j = *last_update+1; j < k; j++)
        {
            // skip if x(P[j]) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[P[j]]));
            if (sgn == 0)       { continue; }

            // perform j-th IPGE update for x
            SPEX_CHECK(spex_ipge(sv_x, x_scale, h, i_2ndlast, L->v[j], P,
                P_inv, sd, SPEX_2D(S, 1, j), SPEX_2D(S, 3, j), Ldiag[j], j));
        }
        *last_update = k-1;
    }

    // check if x[P[k]] is zero
    SPEX_CHECK(SPEX_mpz_sgn(&sgn, sv_x->x[P[k]]));
    if (sgn == 0)
    {
        *last_update = k;
    }

    return SPEX_OK;
}
