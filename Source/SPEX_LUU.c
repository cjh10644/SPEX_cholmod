//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_LUU.c: perform LU update for column replacement
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform LU update for column replacement

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);         \
    SPEX_FREE(h);

#include "spex_internal.h"

SPEX_info SPEX_LUU
(
    SPEX_matrix *L,
    SPEX_matrix *U,
    mpz_t *d,
    mpz_t *sd,
    int64_t *P,
    int64_t *Q,
    mpz_t *vk,
    int64_t *k
)
{
    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t ks;

    // get the row-wise nnz pattern for L and column-wise nnz pattern for U
    SPEX_CHECK(spex_get_nnz_pattern());

    k = Q_inv(k);
    // build Lk_dense_col and Uk_dense_row
    SPEX_CHECK(spex_get_dense_Ak(Lk_dense_col, L, inext));
    SPEX_CHECK(spex_get_dense_Ak());

    // push column k to position n-1
    while (k < n-1)
    {
        if (inext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[P[inext]]));
            if (sgn == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&inext, L, Lk_dense_col, P_inv,k));
            }
        }
        if (jnext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_col[Q[jnext]]));
            if (sgn == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&jnext, U, Uk_dense_col,
                Q_inv,k));
            }
        }

        if (inext == n)
        {
            // if L(:,k) has zero off-diagonal, then only perform dppu
            ks = n-1;
            while (ks > k+1)
            {
                if (L_row_offdiag[ks] <= k && U_col_offdiag[ks] <= k)
                {
                    break;
                }
                ks--;
            }
            if (jnext > ks)
            {
                SPEX_CHECK(spex_dppu1());
            }
            else
            {
                SPEX_CHECK(spex_dppu2());
            }
        }
        else
        {
            // remaining entries in current row of U are 0s, check inserted col
            if (jnext == n)
            {
                if (L(k,n) == 0) //TODO
                {
                    SPEX_FREE_WORK;
                    return SPEX_SINGULAR;
                }
                else
                {
                    ks = n;
                    SPEX_CHECK(spex_cppu());
                }
            }
            else if (U_col_offdiag[jnext] == k || inext == k+1)
            {
                ks = jnext;
                SPEX_CHECK(spex_cppu());
            }
            else
            {
                ks = (inext < jnext) ? inext: jnext-1;
                while (ks > k+1)
                {
                    if (L_row_offdiag[ks] <= k && U_col_offdiag[ks] <= k)
                    {
                        break;
                    }
                    ks--;
                }
                SPEX_CHECK(spex_dppu1());
                if (inext == n)
            }
        }

        k = ks;
    }
    // move the entry from Lk_dense_col

    // update row/column permutation
    Q[n] = Q[n+1];
    
    SPEX_FREE_WORK;
    return SPEX_OK;
}
