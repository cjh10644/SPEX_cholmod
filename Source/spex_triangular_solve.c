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

#define SPEX_FREE_ALL               \
    SPEX_MPQ_CLEAR(Lvj_scale);

#include "spex_internal.h"

SPEX_info spex_triangular_solve // perform REF triangular solve for LDx=v
(
    mpz_t *x,           // the scattered version of solution for LDx=v, whose
                        // entries have been initialized to be same as v
    mpq_t x_scale,      // pending scale for x
    int64_t *h,         // history vector for x
    SPEX_vector *v,     // compressed form of sparse vector v
    int64_t *last_update,// the number of  finished IPGE iterations
    int64_t *i_2ndlast, // i_2ndlast is the index of the found the 2nd last nnz 
                        // entry of x[P], this could be NULL if not needed
    const SPEX_matrix *L,// matrix L
    const mpq_t *S,     // the pending scale factor matrix
    const int64_t k,    // compute x up to k-th IPGE iteration
    const int64_t *P,   // row permutation
    const int64_t *P_inv,// inverse of row permutation
    const mpz_t *sd     // array of scaled pivots
)
{
    SPEX_info info;
    if (!x || !v || !h || !last_update || !L || !S || !P || !P_inv || !sd)
    {
        return SPEX_INCORRECT_INPUT;
    }

    int sgn;
    mpq_t Lvj_scale; SPEX_MPQ_SET_NULL(Lvj_scale);

    if (*last_update < k-1)
    {
        SPEX_CHECK(SPEX_mpq_init(Lvj_scale));
        // TODO iterate sorted nnz pattern?
        for (int64_t j = *last_update+1; j < k; j++)
        {
            // skip if x(P[j]) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[P[j]]));
            if (sgn == 0)       { continue; }

            // compute the pending scale for L->v[j]
            SPEX_CHECK(SPEX_mpq_mul(Lvj_scale,
                SPEX_2D(S, 1, j), SPEX_2D(S, 3, j)));
            // perform j-th IPGE update for x
            SPEX_CHECK(spex_ipge(x, v, x_scale, h, NULL, i_2ndlast, L->v[j], P,
                P_inv, sd, Lvj_scale, Ldiag[j], j));
        }
        *last_update = k-1;
    }

    // check if x[P[k]] is zero
    SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[P[k]]));
    if (sgn == 0)
    {
        *last_update = k;
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
