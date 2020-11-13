//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_LUU.c: perform LU update for column replacement
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform LU update for column replacement

#define SPEX_FREE_ALL                \
    SPEX_FREE(&h);                   \
    SPEX_FREE(&h_for_vk);            \
    SPEX_FREE(&Ldiag);               \
    SPEX_FREE(&Lr_offdiag);          \
    SPEX_FREE(&Uci);                 \
    SPEX_FREE(&Ucp);                 \
    SPEX_FREE(&Ucx);                 \
    SPEX_FREE(&Lk_dense_col);        \
    SPEX_FREE(&Uk_dense_row);        \
    SPEX_FREE(&vk_dense);            \
    SPEX_MPZ_CLEAR();            \
    SPEX_MPQ_CLEAR(vk_scale);             \
    SPEX_MPQ_CLEAR();   \
    SPEX_MPQ_CLEAR();         \

#include "spex_internal.h"

SPEX_info SPEX_LUU
(
    SPEX_matrix *L,         // stored in compressed-column form
    SPEX_matrix *U,         // stored in comptessed-row form
    mpz_t *d,               // an array of size n that stores the unscaled pivot
    mpz_t *sd,              // an array of size n that stores the scaled pivot
    mpq_t *S,               // an array of size 3*n that stores pending scales
    int64_t *P,             // row permutation
    int64_t *P_inv,         // inverse of row permutation
    int64_t *Q,             // column permutation
    int64_t *Q_inv,         // inverse of column permutation
    SPEX_vector *vk,        // the inserted column
    bool keep_vk,           // indicate if the vector vk will be keep unchanged
    int64_t k,              // the column index that vk will be inserted
    const SPEX_option *option// command parameters
)
{
    // initialize workspace
    SPEX_info info;
    if (!spex_initialized()) {return SPEX_PANIC;}

    if (L->n != U->n)
    {
        return SPEX_INCORRECT_INPUT;
    }
    int sgn;
    int64_t ks, p, i, j, inext, jnext, n = L->n;
    int64_t *h = NULL, *h_for_vk = NULL, *Ldiag = NULL, *Lr_offdiag = NULL,
        *Uci = NULL, *Ucp = NULL, *Ucx = NULL, *map = NULL;
    spex_scattered_vector *Lk_dense_col = NULL, *Uk_dense_row = NULL,
        *vk_dense = NULL;
    mpq_t vk_scale; SLIP_MPQ_SET_NULL(vk_scale);
    SPEX_CHECK(SPEX_mpq_init(vk_scale));

    h        = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    h_for_vk = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    map      = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    if (!h || !h_for_vk || !map)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // get the row-wise nnz pattern for L and column-wise nnz pattern for U
    SPEX_CHECK(spex_get_nnz_pattern(&Ldiag, &Lr_offdiag, &Uci, &Ucp, &Ucx,
        L, U, option));

    k = Q_inv(k);
    // build Lk_dense_col and Uk_dense_row, remove explicit 0 if k == 0
    SPEX_CHECK(spex_get_scattered_v(Lk_dense_col, L->v[k], (k == 0), false));
    SPEX_CHECK(spex_get_scattered_v(Uk_dense_row, U->v[k], (k == 0), false));

    // initialize environment for the inserted column
    SPEX_CHECK(spex_get_scattered_v(vk_dense, vk, true, keep_vk));
    int64_t last_update = -1;
    int64_t vk_2ndlastnz = -1;
    SPEX_CHECK(SPEX_mpq_set_ui(vk_scale, 1, 1));
    for (p = 0; p < vk_dense->nz; p++)
    {
        i = vk_dense->i[p];
        if (P_inv[i] > vk_2ndlastnz && P_inv[i] != n-1)
        {
            vk_2ndlastnz = P_inv[i];
        }
        h_for_vk[i] = -1;
    }

    // initialize certain variables required by the loop
    int64_t maximum, last_max_ks = k;
    int64_t use_col_n = 0; // 0: unknown; 1: use; -1: don't use
    SPEX_CHECK(spex_find_next_nz(&inext, Lk_dense_col, P_inv, k));
    SPEX_CHECK(spex_find_next_nz(&jnext, Uk_dense_col, Q_inv, k));

    // push column k to position n-1
    while (k < n-1)
    {
        // get the k-th IPGE update of inserted column
        SPEX_CHECK(spex_ref_triangular_solve(vk_dense, vk_scale, h_for_vk, vk,
            &last_update, &vk_2ndlastnz, L, S, k, P, sd));

        if (inext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[P[inext]]));
            if (sgn == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&inext, Lk_dense_col, P_inv, k));
            }
        }
        if (jnext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_col[Q[jnext]]));
            if (sgn == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&jnext, Uk_dense_col, Q_inv, k));
            }
        }
        // report singular if 
        // - remaining entries in current row of U are 0s and the current row
        //   of vk is also 0 (if last_update == k, vk_dense[P[k]] is zero).
        // - OR all entries below (k-1)-th row in vk are zeros
        // - OR all off-diagonal entries below (k-1)-th row in vk
        //   and column n-1 of U are zeros
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk_dense[P[n-1]]));
        if ((jnext == n && last_update == k) ||
            (vk_2ndlastnz <= k &&
             (sgn == 0 || Ucp[Q[n-1]+1]-Ucp[Q[n-1]] == 1 ||
              Uci[Ucp[Q[n-1]+1]-2] < k                     )))
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }

        //----------------------------------------------------------------------
        // if L(:,k) has zero off-diagonal, then only perform dppu, which will
        // maintain the sparsity of L(:,k). Use dppu1 if possible.
        // When arriving the last iteration, always use the inserted column
        // if possible, since we can perform less IPGE iterations for it.
        //----------------------------------------------------------------------
        if (inext == n)
        {
            // build the map to find ks in the first loop
            if (n-2 > last_max_ks)
            {
                for (j = k+1; j <= n-2; j++)
                {
                    if (Ucp[Q[j]+1] - Ucp[Q[j]] == 1)
                    {
                        // no off-diag in column Q[j] of U
                        maximum = SPEX_MAX(Lr_offdiag[P[j]], -1);
                    }
                    else
                    {
                        maximum = SPEX_MAX(Lr_offdiag[P[j]],Uci[Ucp[Q[j]+1]-2]);
                    }
                    for (i = k; i < j;)
                    {
                        if (maximum <= i)
                        {
                            map[i] = j;
                            break;
                        }
                        i = map[i];
                    }
                }
                last_max_ks = n-2;
            }

            // use the inserted column only when its last entry is nnz and
            // using it instead of column n-1 can make a bigger jump.
            if (use_col_n == 0 && Lr_offdiag[P[n-1]] < k)
            {
                if (sgn == 0)// vk_dense[P[n-1]] == 0
                {
                    // the inserted column cannot be used
                    use_col_n = -1;
                }
                else
                {
                    if (U_col_offdiag[Q[n-1]] > vk_2ndlastnz)
                    {
                        use_col_n = 1;
                    }
                    else
                    {
                        use_col_n = -1;
                    }
                }
                maximum = SPEX_MAX(L_row_offdiag[P[n-1]],
                        use_col_n == -1 ? U_col_offdiag[Q[n-1]] : vk_2ndlastnz);
                for (i = k; i < n-1;)
                {
                    if (maximum <= i)
                    {
                        map[i] = n-1;
                        break;
                    }
                    i = map[i];
                }

                last_max_ks = n-1;
            }

            // get ks from the map
            ks = map[k];
            if (ks = n-1 && use_col_n == 1)
            {
                // apply history update and pending scaling factor to vk[P[k]]
                real_h = SPEX_FLIP(h_for_vk[P[k]]);
                if (real_h < k-1)
                {
                    SPEX_CHECK(SPEX_mpz_mul(vk[P[k]], vk[P[k]], sd[k-1]));
                    if (real_h > -1)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(vk[P[k]],
                                                     vk[P[k]], sd[real_h]));
                    }
                }
                SPEX_CHECK(SPEX_mpz_divexact(vk[P[k]],
                                        vk[P[k]], SPEX_MPQ_DEN(vk_scale)));
                SPEX_CHECK(SPEX_mpz_mul(vk[P[k]],
                                        vk[P[k]], SPEX_MPQ_NUM(vk_scale)));

                ks = n;
            }
            if (jnext > ks || (ks == n && jnext == n))
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
            // if jnext == n, swapping columns k and n will be more efficient,
            // since there is no need to backtrack column n (we just perform
            // k-th IPGE iteration for column n), and column n-1 can be updated
            // by scaling after CPPU.
            // if jnext == n-1 and U(k,n) != 0, swapping columns k and n will
            // be more efficient, since no matter column n or n-1 is swapped,
            // column n needs n IPGE iterations, while if column n-1 is
            // swapped, additional backtracking for column n-1 needs to be
            // performed.
            // if jnext == n-1 and U(k,n) == 0, use column n only when
            // U_col_offdiag[n] < k && L_row_offdiag[n-1] <= k && inext >= n-1.
            if (use_col_n == 0 && jnext >= n-1)
            {
                if (last_update == k) // vk[P[k]] == 0
                {
                    if (sgn != 0 /*vk[P[n-1]] != 0*/ && vk_2ndlastnz <= k &&
                        L_row_offdiag[P[n-1]] <= k && inext >= n-1)
                    {
                        // swap columns n and n-1
                        ks = n;
                        SPEX_CHECK(spex_dppu1());
                        break;
                    }
                    else
                    {
                        // the inserted column cannot be used
                        use_col_n = -1;
                    }
                }
                else
                {
                    use_col_n = 1;
                    ks = n;
                    SPEX_CHECK(SPEX_cppu());
                }
            }
            // this implicitly includes the case of jnext == k+1
            if (inext == k+1 ||
                (use_col_n <= 0 && U_col_offdiag[Q[jnext]] == k))
            {
                ks = jnext;
                SPEX_CHECK(spex_cppu());
            }
            else
            {
                ks = SPEX_MIN(n-2, (inext < jnext) ? inext: jnext-1);
                // build the map to find ks if current map is out of date.
                // all the swaps (i.e., pivot updates) except the last one
                // using this map will not change the nnz patter of current
                // frame, since only scaling will be involved.
                if (ks > last_max_ks)
                {
                    SPEX_ASSERT(k+1 >= last_max_ks);
                    for (j = k+1; j <= ks; j++)
                    {
                        maximum = SPEX_MAX(L_row_offdiag[P[j]],
                                           U_col_offdiag[Q[j]]);
                        for (i = k; i < j;)
                        {
                            if (maximum <= i)
                            {
                                map[i] = j;
                                break;
                            }
                            i = map[i];
                        }
                    }
                    last_max_ks = ks;
                }

                // get ks from the map
                ks = map[k];

                SPEX_CHECK(spex_dppu1());
                if (inext == n)
            }
        }

        // update k
        if(ks != n)
        {
            k = ks;
        }
        else
        {
            break;
        }
    }
    // swap inserted column vk with column k
    // update d[k]=vk[k], sd[k]=U(k,k)=vk[k]*v_scale
    // S(:,k)=[v_scale;1;1]

    // move the entry from Uk_dense_row except entry in col Q[k]
    // set U(n-1,n-1)=L(n-1,n-1)= Uk_dense_row[Q[n-1]] when using_col_n;
    // delete Lk_dense_col


    // update row/column permutation
    Q[n] = Q[n+1];
    
    SPEX_FREE_ALL;
    return SPEX_OK;
}
