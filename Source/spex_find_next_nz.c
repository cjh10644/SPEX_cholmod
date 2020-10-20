//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_find_next_nz.c: find the next index of next nz in given
// vector
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to find the row index of fisrt off diagonal of
//          L(P,k) or the column index of the first off diagonal of U(k,Q)


#include "spex_internal.h"

SPEX_info spex_find_next_nz
(
    int64_t *next,
    SPEX_matrix *A,
    SPEX_vector *Ak_dense,
    int64_t *perm_inv
    int64_t k
)
{
    SPEX_info info;
    inext = n;
    for (p = 0; p < A->v[k]->nz; p++)
    {
        i = A->v[k]->i[p];
        if (k == perm_inv[i]) { continue; }
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Ak_dense[i]));
        if (sgn != 0 && perm_inv[i] < inext)
        {
            inext = i;
        }
    }
    return SPEX_OK;
}
