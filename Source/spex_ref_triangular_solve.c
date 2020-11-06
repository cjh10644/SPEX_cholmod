//------------------------------------------------------------------------------
//SPEX_CHOLMOD/spex_triangular_solve.c: perform REF triangular solve up to
//specified iteration.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform REF triangular solve for LDx=v up to
// specified IPGE iteration when both L and v are sparse. This function is used
// when computing the inserted column for U. For forward solving LDUx=b where b
// is mostly considered and stored as dense vector, we will use
// spex_forward_sub. This function returns when error occurs, or when the
// specified entry x[target_index] (if target_index > -1) turns out to be
// exactly cancelled, or when all entries in x are computed
// up to specified IPGE iteration.

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);

#include "spex_internal.h"

SPEX_info spex_triangular_solve
(
    mpz_t *x,
    mpq_t Sx,
    bool *target_IS_zero,
    SPEX_matrix *L,
    SPEX_vector *v,
    int64_t k,
    int64_t target_index,// the target entry will be x[P[target_index]]
    int64_t *h,
    int64_t *P,
    mpz_t *d,
    mpz_t *sd,
    int64_t *last_update
)
{
    SPEX_info info;

    if (*last_update < -1)
    {
        for (p = 0; p < v->nz; p++)
        {
            i = v->i[p];
            SPEX_CHECK(SPEX_mpz_swap(x[i], v->x[p]));
        }
        *last_update = -1;
    }

    if (*last_update < k-1)
    {
        for (j = *last_update+1; j < k; j++)// TODO iterate sorted nnz pattern?
        {
            // skip if x(P[j]) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[P[j]]));
            if (sgn == 0)       { continue; }

            // perform j-th IPGE update for x
            SPEX_CHECK(spex_ipge());
            *last_update = j;
        }

        // check if the target entry is zero
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[P[target_index]]));
        if (sgn == 0)
        {
            *target_IS_zero = true;
        }
    }
    SPEX_FREE_WORK;
    return SPEX_OK;
}
