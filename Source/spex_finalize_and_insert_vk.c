//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_finalize_and_insert_vk.c: perform history update for
// entries that would be in L and insert entries that would in U to
// corresponding row of U.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform history update for entries that
// would be in L and insert entries that would in U to corresponding row of U.
// Entries that would appear in U are removed from the nnz pattern, while
// entries that would be in L are moved to SPEX_vector vk and remained in nnz
// pattern. The k-th entry (i.e., the new k-th pivot) is kept in both vk and
// d[k] but not inserted to U.

#include "spex_internal.h"

SPEX_info spex_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                      // LDx=v using the first k-1 columns of L
    int64_t *h,       // history vector for vk_dense
    SPEX_matrix *U,   // matrix U
    SPEX_matrix *L,   // matrix L
    mpq_t *S,         // array of size 3*n that stores pending scales
    mpz_t *d,         // array of unscaled pivots
    const int64_t *Ldiag,// L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const mpz_t *sd,  // array of scaled pivots
    const int64_t *Q, // the column permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t k,  // the column index in L that vk_dense will be inserted
    const int64_t diag,// the index of entry in vk_dense that will be diagonal
    const mpq_t one   // a const mpq number, just to avoid constantly alloc
)
{
    SPEX_info info;
    int64_t i, p = 0, vk_nz = vk_dense->nz, n = U->n;

    // move entries to U
    while(p < vk_nz)
    {
        i = vk_dense->i[p];
        if (P_inv[i] < k)
        {
            SPEX_CHECK(spex_insert_new_entry(vk_dense->x[i], U->v[i],
                SPEX_2D(S, 2, i), L->v[i], SPEX_2D(S, 1, i), SPEX_2D(S, 3, i),
                d[i], Q[k], Ldiag[i], one));
            vk_nz--;
            vk_dense->i[p] = vk_dense->i[vk_nz];
        }
        else
        {
            p++;
        }
    }

    // check if L->v[k] needs more space for all remaining entries
    if (vk_nz > L->v[k]->nzmax)
    {
        SPEX_CHECK(spex_expand_vector(L->v[k], vk_nz));
    }

    // move the remaining entries to L->v[k]
    for (p = 0; p < vk_nz; p++)
    {
        i = vk_dense->i[p];
        h[i] = SPEX_FLIP(h[i]);
        if (h[i] < k-1)
        {
            SPEX_CHECK(SPEX_mpz_mul(vk_dense->x[i], vk_dense->x[i], sd[k-1]));
            if (h[i] > -1)
            {
                SPEX_CHECK(SPEX_mpz_divexact(vk_dense->x[i],
                                             vk_dense->x[i], sd[h[i]]));
            }
        }
        if (P_inv[i] == diag)
        {
            SPEX_CHECK(SPEX_mpz_set(d[k], vk_dense->x[i]));
        }
        SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[p], vk_dense->x[i]));
        L->v[k]->i[p] = i;
    }
    L->v[k]->nz = vk_nz;

    return SPEX_OK;
}
