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
    SPEX_vector *vk,
    int64_t *k
)
{
    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t ks;

    h_for_vk = SPEX_calloc(h, n*sizeof(int64_t));
    if (h_for_vk == NULL)
    {
        SPEX_FREE_WORK;
        return SPEX_OUT_OF_MEMORY;
    }

    // get the row-wise nnz pattern for L and column-wise nnz pattern for U
    SPEX_CHECK(spex_get_nnz_pattern());

    k = Q_inv(k);
    // build Lk_dense_col and Uk_dense_row
    SPEX_CHECK(spex_get_dense_Ak(Lk_dense_col, L, inext));
    SPEX_CHECK(spex_get_dense_Ak());

    // initialize history vector for the inserted column
    vk_2ndlastnz = -1;
    for (p = 0; p < vk->nz; p++)
    {
        i = vk->i[p];
        if (i > vk_2ndlastnz && i != n-1) { vk_2ndlastnz = i;}
        h_for_vk[i] = -1;
    }
    // update the nnz pattern up to k-th IPGE
    for (j = 0; j < k; j++)// TODO iterate sorted nnz pattern?
    {
        if (h_for_vk[P[j]] == -1) // if this entry is nnz
        {
            for (p = 0; p < L->v[j]->nz; p++)
            {
                i = L->v[j]->i[p];
                if (h_for_vk[i] != -1) // if this entry isn't in the nnz pattern
                {
                    if (i > vk_2ndlastnz && i != n-1) { vk_2ndlastnz = i;}
                    h_for_vk[i] = -1;
                    vk->i[vk->nz] = i;
                    vk->nz ++;
                }
            }
        }
    }

    // initialize certain variables required by the loop
    last_max_ks = k;
    use_col_n = 0; // 0: unknown; 1: use; -1: don't use

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
        // report singular if remaining entries in current row of U are 0s and
        // the current row of the inserted col is also 0
        if (jnext == n && h_for_vk[k] > -1)
        {
            SPEX_FREE_WORK;
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
                last_max_ks = n-2;
            }

            // use the inserted column only when its last entry is nnz and
            // using it instead of column n-1 can make a bigger jump.
            if (use_col_n == 0 && L_row_offdiag[P[n-1]] < k)
            {
                if (U(n-1, n) == 0)
                {
                    use_col_n = -1;
                }
                else
                {
                    // get the k-th IPGE update of inserted column
                    SPEX_CHECK(spex_ref_triangular_solve());
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk[P[n-1]]));
                    if ((target_IS_zero &&
                         (jnext == n || (vk_2ndlastnz <= k && sgn == 0))) ||
                        (U_col_offdia[Q[n-1]] < k && vk_2ndlastnz < k)      ) 
                    {
                        // report singularity if U(k,n) is exactly cancelled
                        // or all entries below (k-1)-th row in vk are zeros
                        // or all off-diagonal entries below (k-1)-th row in vk
                        // and column n-1 of U are zeros
                        SPEX_FREE_WORK;
                        return SPEX_SINGULAR;
                    }
                    else if (sgn == 0)
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
                // TODO handle column n in dppu and cppu
                // get the k-th IPGE update of inserted column
                SPEX_CHECK(spex_ref_triangular_solve());
                ks = n;
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
                // prepare the inserted column to be swapped with k-th
                // column, i.e., perform (k-1)-th IPGE update for the
                // inserted column.
                SPEX_CHECK(spex_ref_triangular_solve());
                if (target_IS_zero)
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk[P[n-1]]));
                    if (jnext == n || (vk_2ndlastnz <= k && sgn == 0))
                    {
                        // report singularity if U(k,n) is exactly cancelled
                        // or all entries below (k-1)-th row in vk are zeros
                        SPEX_FREE_WORK;
                        return SPEX_SINGULAR;
                    }
                    else
                    {
                        if (sgn != 0 && vk_2ndlastnz < k &&
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

        // update the history vector for the inserted column by adding fillin
        for (j = k; j < ks; j++)// TODO iterate sorted nnz pattern?
        {
            if (h_for_vk[P[j]] == -1) // if this entry is nnz
            {
                for (p = 0; p < L->v[j]->nz; p++)
                {
                    i = L->v[j]->i[p];
                    if (h_for_vk[i] > -1) // add only if not in the nnz pattern
                    {
                        if (i > vk_2ndlastnz && i != n-1) { vk_2ndlastnz = i;}
                        h_for_vk[i] = -1;
                        vk->i[vk->nz] = i;
                        vk->nz ++;
                    }
                }
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
    
    SPEX_FREE_WORK;
    return SPEX_OK;
}
