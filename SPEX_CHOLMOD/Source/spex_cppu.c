//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_cppu.c: perform column permutation pivot update
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform column permutation pivot update
// when the submatrix (formed by rows and columns k to ks) has the following
// pattern
//       x 0 0 0 x       <- row k
//       . x . . .
//       . . x . .
//       . . . x .
//       . . . . x       <- row ks
//
//       ^       ^
//       |       |
//     col k   col ks
//       
// This function will swap columns k and ks in L and U. Noted that the columns
// of U are permuted implicitly via the permutation matrix based on Q.

#define SPEX_FREE_ALL                \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(tmpq);            \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_internal.h"

SPEX_info spex_cppu
(
    SPEX_matrix *L,  // matrix L
    SPEX_matrix *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    const mpq_t vk_scale,// scale factor for newly inserted column vk, which
                     // should be in col k of L in the last iteration when used.
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    const int64_t *P,// row permutation
    const int64_t *P_inv,// inverse of row permutation
    int64_t *Ldiag,  // L(P(k),k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                     // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k, // current column index 0 <= k < n
    const int64_t ks // index of the diagonal to be swapped with, [0,n)
)
{
    printf("using cppu\n");
    // initialize workspace
    SPEX_info info;
    int sgn, r;
    int64_t pk, ck, pks, cks, pi, ci, i, j, n = U->n;
    // the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
    int64_t Ucp_k_ks;
    *inext = n;
    *jnext = n;

    mpq_t pending_scale, tmpq, one;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(tmpq);
    mpz_t Uiks, tmpz; SPEX_MPZ_SET_NULL(Uiks); SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(tmpq));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Uiks));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    if (ks == n)
    {
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping column k with
        // column n-1. Then we set sd[k] to d[k]=vk[P[k]] with pending scaling
        // factor applied
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Q[k]], sd[k]));
        SPEX_CHECK(SPEX_mpz_divexact(sd[k],  d[k], SPEX_MPQ_DEN(vk_scale)));
        SPEX_CHECK(SPEX_mpz_mul     (sd[k], sd[k], SPEX_MPQ_NUM(vk_scale)));

        // get the scale for entries between frames k and n-1
        // pending_scale = sd(k)/Uk_dense_row[Q[k]]
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        // if the inserted column is used, we won't need to perform
        // backtracking. Instead, we just need to perform RwSOP. If U(k,
        // Q[n-1]) == 0, RwSOP is simplified as pure scaling for all frame from
        // k+1:n-1.
        //
        // Otherwise (which won't happen due to heuristic), we need to first
        // compute the (n-1)-th IPGE iteration for the column k, and use the
        // result to perform RwSOP on column n-1 of U.
#if 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[Q[n-1]]));
        if (sgn == 0)
        {
#endif
            // just need to perform RwSOP by scaling all frame k+1:n-1
            for (j = k+1; j < n; j++)
            {
                // S(3,k+1:n-1) = S(3,k+1:n-1)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                        SPEX_2D(S, 3, j), pending_scale));
                // sd(k+1:n-1) = sd(k+1:n-1)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }

            // move data from Uk_dense_row, but there is no entry that needs
            // to move
            U->v[k]->nz = 0;
#if 0
        }
        else
        {
            // S(2,k) = S(2,k)*S(3,k)
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, k),
                                    SPEX_2D(S, 2, k), SPEX_2D(S, 3, k)));
            // U(k,Q[n-1]) = U(k, Q[n-1])*S(2,k)
            SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row->x[Q[n-1]],
                      Uk_dense_row->x[Q[n-1]], SPEX_MPQ_DEN(SPEX_2D(S, 2, k))));
            SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row->x[Q[n-1]],
                      Uk_dense_row->x[Q[n-1]], SPEX_MPQ_NUM(SPEX_2D(S, 2, k))));

            // perform 1 IPGE iteration on Lk_dense_col using vk (which has
            // been inserted as L->v[k]) and update the history vector. Then
            // use L to perform the remaining IPGE update till (n-1)-th
            // iteration. Finally, use the result to update column n-1 of U.
            //
            // It should be noted that, since the resulted column in the
            // (n-1)-th IPGE iteration is computed using column k of L instead
            // of the inserted column, its sign should be flipped when applying
            // RwSOP.
            // initialize history vector
            for (pk = 0; pk < Lk_dense_col->nz; pk++)
            {
                ck = Lk_dense_col->i[pk];
                h[ck] = SPEX_FLIP(k-1);
            }
            for (pk = 0; pk < L->v[k]->nz; pk++)
            {
                ck = L->v[k]->i[pk];
                if (ck == P[k])     { continue; }// skip updating L(P(k),k)
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[k]->x[pk]));
                if (sgn == 0)       { continue;   }

                // L(ck,k) = (L(ck, k)*vk(P[k])-L(P[k],k)*vk(ck))/sd[k-1]
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn != 0)
                {
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                            Lk_dense_col->x[ck], d[k]));
                }
                else if (h[ck] >= -1) //this entry wasn't in the nnz pattern
                {
                    // insert new entry in the nonzero pattern
                    Lk_dense_col->i[Lk_dense_col->nz] = ck;
                    Lk_dense_col->nz++;
                }
                SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[P[k]], L->v[k]->x[pk]));
                if (k > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck], sd[k-1]));
                }
                h[ck] = SPEX_FLIP(k);
            }
            // S(1,k) = vk_scale, the existing pending scale for column k of L
            // can be ignored when performing RwSOP
            SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 1, k), vk_scale));

            // perform IPGE and RwSOP
            for (i = k+1; i < n; i++)
            {
                // S(2,i) = S(2, i)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, i),
                                        SPEX_2D(S, 3, i), pending_scale));

                // if Lk_dense_col[P[i]] == 0, just perform RwSOP
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[P[i]]));
                if (sgn == 0) 
                {
                    SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                            SPEX_MPQ_DEN(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                            SPEX_MPQ_NUM(pending_scale)));
                    continue;
                }

                // perform RwSOP for row i with flipped-sign entries in
                // Lk_dense_col. All entries in row i of U must be SCALEUP such
                // that S(:,i)=[S(1,i)*S(3,i);1;1]

                // set S(:,i) = [S(1,i)*S(3,i); S(2,i)*S(3,i); 1]
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, i),
                                        SPEX_2D(S, 1, i), SPEX_2D(S, 3,i)));
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, i),
                                        SPEX_2D(S, 2, i), SPEX_2D(S, 3,i)));
                SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));

                if (i != n-1)
                {
                    // perform i-th IPGE update for Lk_dense_col
                    SPEX_CHECK(spex_ipge(Lk_dense_col, SPEX_2D(S, 1, k), h,
                        NULL, L->v[i], P, P_inv, (const mpz_t*) sd,
                        SPEX_2D(S, 1, i), one, Ldiag[i], i));
                    GOTCHA;
                }
                else
                {
                    // finish history update for last entry of Lk_dense_col
                    ck = P[n-1];
                    int64_t real_h = SPEX_FLIP(h[ck]);
                    if (real_h < n-1) // require history update
                    {
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck], sd[n-1]));
                        if (real_h > -1)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck], sd[real_h]));
                        }
                    }
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                        Lk_dense_col->x[ck], SPEX_MPQ_DEN(SPEX_2D(S, 1, k))));
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                        Lk_dense_col->x[ck], SPEX_MPQ_NUM(SPEX_2D(S, 1, k))));
                    GOTCHA;
                }

                int64_t p = -1; // the pointer to U(i,Q(n-1))
                // iterate all nnz in row i of U
                for (pi = 0; pi < U->v[i]->nz; pi++)
                {
                    ci = U->v[i]->i[pi];
                    if (Q_inv[ci] < n-1)
                    {
                        // apply S(2,i) to U(i,Q(i:n-2))
                        SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi],
                               U->v[i]->x[pi], SPEX_MPQ_DEN(SPEX_2D(S, 2, i))));
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi],
                               U->v[i]->x[pi], SPEX_MPQ_NUM(SPEX_2D(S, 2, i))));
                        // set sd[i] = U(i, Q(i))
                        if (ci == Q[i])
                        {
                            SPEX_CHECK(SPEX_mpz_set(sd[i], U->v[i]->x[pi]));
                        }
                    }
                    else
                    {
                        p = pi;
                    }
                }
                // perform RwSOP to U(i,Q(n-1)), POSSIBLE FILLIN
                // sign is changed here due to column swap
                // U(i,ci)= U(i,ci)*S(2,i) +
                //        Lk_dense_col(P[i])*U(k,ci)/Lk_dense_col[P[k]]
                ci = Q[n-1];
                if (p > -1)
                {
                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_NUM(SPEX_2D(S, 2, i))));
                    SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_DEN(SPEX_2D(S, 2, i))));
                    SPEX_CHECK(SPEX_mpz_mul(tmpz,
                               Lk_dense_col->x[P[i]], Uk_dense_row->x[ci]));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                               Lk_dense_col->x[P[k]]));
                    SPEX_CHECK(SPEX_mpz_add(U->v[i]->x[p],U->v[i]->x[p], tmpz));
                }
                else // U(i,Q(n-1)) was not in the nnz pattern
                {
                    p = U->v[i]->nz;
                    // reallocate the nonzero pattern if needed
                    if (p == U->v[i]->nzmax)
                    {
                        SPEX_CHECK(SPEX_vector_realloc(U->v[i],
                            2*(U->v[i]->nzmax)));
                    }
                    // insert new entry in the nonzero pattern
                    U->v[i]->i[p] = ci;
                    U->v[i]->nz++;

                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                                  Lk_dense_col->x[P[i]], Uk_dense_row->x[ci]));
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p],
                                  U->v[i]->x[p], Lk_dense_col->x[P[k]]));
                }
                SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
            }
            // L(P(n-1),n-1) = d(n-1) = sd(n-1) = U(n-1, Q(n-1))
            SPEX_CHECK(SPEX_mpz_set(L->v[n-1]->x[0], U->v[n-1]->x[0]));
            SPEX_CHECK(SPEX_mpz_set(   d[n-1]      , U->v[n-1]->x[0]));
            SPEX_CHECK(SPEX_mpz_set(  sd[n-1]      , U->v[n-1]->x[0]));
            // S(:,n-1) = ones
            SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 1, n-1), one));
            SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 2, n-1), one));
            SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 3, n-1), one));

            // move data from Uk_dense_row, there is only one entry that needs
            // to move, which is U(k,Q[n-1]) and has no pending scale
            SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[0], Uk_dense_row->x[Q[n-1]]));
            U->v[k]->i[0] = Q[n-1];
            U->v[k]->nz = 1;
        }
#endif

        SPEX_FREE_ALL;
        return SPEX_OK;
    }
    //-------------------------------------------------------------------------
    // Backtracking column ks of L and U, the backtracking result will be moved
    // to column k of L.
    //-------------------------------------------------------------------------
    // When U(k+1:ks-1,Q(ks)) are all zero(s), there is no need to make a copy
    // of Lk_dense_col. Instead, we can direct operate on Lk_dense_col.
    // Otherwise, Lk_dense_col will firstly be moved to L->v[k]->x in a
    // compressed-column form. Then Lk_dense_col will make a copy of ks-th
    // column of scaled L and scaled U(k:ks-1,Q(ks)). Backtracking will be
    // performed on Lk_dense_col for each nonzero in U(k:ks-1, Q(ks)).
    //-------------------------------------------------------------------------
    // initialized to be 2nd last entry in the Q[ks]-th col of U
    Ucp_k_ks = Ucp[Q[ks]+1]-2;
    if (Uci[Ucp_k_ks] == k) // 2nd last nnz in U(:,Q(ks)) is at row k
    {
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        // backtracking jumbled sparse column ks of L using 'dense' column k of
        // L and store the result in Lk_dense_col. This will introduce new
        // entry to Lk_dense_col.
        // Assuming explicit zeros in L(:, ks) resulted from exact
        // cancellation in IPGE update were not removed (SLIP LU keeps those
        // zeros in output L and U), nonzero pattern of L(P(k+1:n+1),k) should
        // be a subset of U(:,ks). Therefore, the backtracking will need to 
        // simply iterate all nonzero in the L(:,ks), and the final Lk_dense_col
        // will have mostly the same nnz pattern as L(:,ks), except L(P(k),k).
        //
        // REMARK:
        // when the IPGE update results in exact cancellation and the resulted
        // zero is removed from L or U, new entry will be inserted to
        // Lk_dense_col during backtracking. In this case, we need to
        // additionally iterate across all nonzeros in Lk_dense_col to find if
        // any row index of nonzero is untouched, then perform the following
        // L(i,k) = L(i,k)*U(k, Q[ks])/L(P(k),k)
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // Uiks = U(k,Q(ks))*S(2,k)*S(3,k)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, k), SPEX_2D(S, 2, k),
                                SPEX_2D(S, 3, k)));
        SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row->x[Q[ks]],
                                SPEX_MPQ_DEN(SPEX_2D(S, 2, k))));
        SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                SPEX_MPQ_NUM(SPEX_2D(S, 2, k))));
        // update d[k] = U(k, Q(ks))
        SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row->x[Q[ks]]));

        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k-1]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
        }
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks-1]));
        // remove common factor in mpq_den and mpq_num
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_2D(S,1,ks)));
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_2D(S,3,ks)));
        // perform backtracking for each nonzero in col ks of L and store
        // results in Lk_dense_col
        // NOTE: this will cause fillin in the k(th) column of L
        // make sure L->v[k] has enough space
        if (L->v[k]->nzmax < L->v[ks]->nzmax+1)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], L->v[ks]->nzmax+1));
        }
        pk = 0;
        int64_t Lk_nz = Lk_dense_col->nz;
        for (pks = 0; pks < L->v[ks]->nz; pks++)
        {
            // row index in column ks of L
            cks = L->v[ks]->i[pks];

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // L(cks,k) = L(cks,k)*Uiks/L(P[k],k)+L(cks,ks)*pending_scale
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            if (sgn != 0)
            {
                // L(cks,k) = floor(L(cks,k)*Uiks/L(P[k],k))
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                           Lk_dense_col->x[cks], Uiks));
                SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[cks],
                           Lk_dense_col->x[cks], Lk_dense_col->x[P[k]]));

                // tmpz = ceil(L(cks,ks)*pending_scale)
                SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[ks]->x[pks],
                           SPEX_MPQ_NUM(pending_scale)));
                SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                           SPEX_MPQ_DEN(pending_scale)));
                // L(cks,k) = L(cks,k)+tmpz
                SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[cks],
                           Lk_dense_col->x[cks], tmpz));
                Lk_nz--;
            }
            else  // faster 
            {
                // L(cks,k) = L(cks,ks)*pending_scale
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                           L->v[ks]->x[pks], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                           Lk_dense_col->x[cks], SPEX_MPQ_NUM(pending_scale)));
            }

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // move nonzero entry from Lk_dense_col to L->v[k]->x
            // NOTE: explicit zero due to exact cancellation in backtracking
            //       is removed
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            if (sgn != 0)
            {
                L->v[k]->i[pk] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[cks]));
                pk++;
            }
        }

        // continue the backtracking process in case explicit zeros are removed
        // in column ks, which means certain entries in column k remain
        // untouched.
        if (Lk_nz > 1) // there must be 1 nonzero left untouched, L(P[k],k)
        {
            if (L->v[k]->nzmax < L->v[k]->nz+Lk_nz)
            {
                SPEX_CHECK(SPEX_vector_realloc(L->v[k], L->v[k]->nz+Lk_nz));
            }
            for (pk = 0; Lk_nz > 1 && pk < Lk_dense_col->nz; pk++)
            {
                // row index in scattered form of column k of L
                ck = Lk_dense_col->i[pk];

                // skip for L(P[k],k)
                if (ck == P[k]) {continue;}
                // remove entries L(P[k+1:ks-1],k), since these should not be
                // appeared in the backtracking result
                if (ck < P[ks] && ck != P[k])
                {
                    SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[ck], 0));
                    Lk_nz--;
                    continue;
                }

                // update entries that is not in column ks
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn == 0) {continue;}

                // L(ck,k) = L(ck,k)*Uiks/L(P[k],k)
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[ck], Uiks));
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck],
                                             Lk_dense_col->x[P[k]]));
                L->v[k]->i[pk] = ck;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[ck]));
                pk++;
                Lk_nz--;
            }
        }
        // update L(P(k),k) and reset Lk_dense_col[P[k]]
        L->v[k]->i[pk] = P[k];
        SPEX_CHECK(SPEX_mpz_set(sd[k], Uiks));
        SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Uiks));
        // no need to put L(P[k],k) to U(k,Q(ks)) since k-th col will be deleted
        SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[P[k]], 0));
        Ldiag[k] = pk;
        pk++;
        // update number of nnz
        L->v[k]->nz = pk;

        // get the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
        // Ucp_k_ks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
    }
    else  // U(ks+1:k-1,Q(ks)) contains nnz(s)
    {
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // construct column k of L based on Lk_dense_col
        // explicit zeros are not removed/skipped, since they should be rarely
        // found here and also this vector will be updated right after the
        // backtracking process.
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        if (Lk_dense_col->nz > L->v[k]->nzmax)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], Lk_dense_col->nz));
        }
        for (pk = 0; pk < Lk_dense_col->nz; pk++)
        {
            ck = Lk_dense_col->i[pk];
            if (ck == P[k])
            {
                Ldiag[k] = pk;  // update Ldiag[k]
            }
            // swap the entries in the Lk_dense_col and L->v[k]->x
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[ck]));
            L->v[k]->i[pk] = ck;
        }
        L->v[k]->nz = Lk_dense_col->nz;
        Lk_dense_col->nz = 0;

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // initialize history vector h and copy L->v[ks]->x to Lk_dense_col.
        // Explicit zero(s) are kept if exist
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        for (pks = 0; pks < L->v[ks]->nz; pks++) 
        { 
            cks = L->v[ks]->i[pks]; 
            // copy value from L->v[ks]->x to Lk_dense_col
            SPEX_CHECK(SPEX_mpz_set(Lk_dense_col->x[cks], L->v[ks]->x[pks]));
            Lk_dense_col->i[pks] = cks;
            // Assuming that zero resulted from exact cancellation is kept in
            // the nnz pattern, then performing backtracking will not introduce
            // new fillin. That is, for any i such that U(i,Q[ks])!=0, the nnz
            // pattern of L(i+1:n, i) is a subset of the nnz pattern of
            // Lk_dense_col(i+1:n). Therefore, the final nnz pattern can be
            // found by L->v[ks]->i and column-wise nnz pattern of U(:,ks).
            // We can just initialize the history vector as:
            //
            // h[cks] = ks-1;
            //
            // However, when explicit zero(s) are always eleminated, the
            // following initialization should be used instead:
            //
            // h[cks] = SPEX_FLIP(ks-1); 
            //
            // With such initialization, entry with h > -1 is clearly not in
            // nnz pattern and any entry in the nnz pattern with h = -1 must be
            // nonzero. In all, any explicit zero with h >= -1 must not be in
            // the nnz pattern.  In this way, we can determine if a zero entry
            // in Lk_dense_col is in the nnz pattern.
            h[cks] = SPEX_FLIP(ks-1);
        }
        Lk_dense_col->nz = L->v[ks]->nz;

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // backtrack column ks of L and U for each nnz in U(k:ks-1,Q(ks))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        pks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
        i = Uci[pks];         // row index
        while (i >= k)
        {
            // skip if U(i, Q[ks]) turns out to be explicit zero
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx[pks]]));
            if (sgn == 0)   { continue; }

            // Uiks = U(i,Q(ks))*S(2,k)*S(3,k)
            if (i == k)
            {
                // store this pointer for later use
                Ucp_k_ks = pks;

                // update S(2,k)
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, k), SPEX_2D(S, 2, k),
                                        SPEX_2D(S, 3, k)));
                // update d[k]
                SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row->x[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row->x[Q[ks]],
                                        SPEX_MPQ_DEN(SPEX_2D(S, 2, k))));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                        SPEX_MPQ_NUM(SPEX_2D(S, 2, k))));
            }
            else
            {
                SPEX_CHECK(SPEX_mpq_mul(pending_scale, SPEX_2D(S, 2, i),
                                        SPEX_2D(S, 3, i)));
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, U->v[i]->x[Ucx[pks]],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                        SPEX_MPQ_NUM(pending_scale)));
            }

            // r = (S(1,i)*S(3,i) == 1)
            SPEX_CHECK(SPEX_mpq_mul(pending_scale, SPEX_2D(S, 1, i),
                                    SPEX_2D(S, 3, i)));
            SPEX_CHECK(SPEX_mpq_equal(&r, pending_scale, one));

            // pending_scale for Lk_dense_col, this will be used when
            // 1. there needs history update, OR
            // 2. S(1,i)*S(3,i) != 1, OR
            // 3. the first loop
            // pending_scale = S(1,ks)*S(3,ks)*sd[i-1] for first loop, OR
            // pending_scale = 1
            if (i > 0)
            {
                SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[i-1]));
            }
            else
            {
                SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
            }
            // for the 1st loop, additional scale S(1,ks)*S(3,ks) needs to apply
            if (pks == Ucp[Q[ks]+1]-2)
            {
                SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale,
                                        SPEX_2D(S, 1, ks)));
                SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale,
                                        SPEX_2D(S, 3, ks)));
            }
            for (pi = 0; pi < L->v[i]->nz; pi++)
            {
                // exclude L(P[i], i)
                if (pi == Ldiag[i]) { continue; }

                // row index of entry in column i of L
                ci = U->v[i]->i[pi];

                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ci]));
                if (sgn != 0)
                {
                    h[ci] = SPEX_FLIP(h[ci]);
                    if (i < h[ci] || r == 0 || pks == Ucp[Q[ks]+1]-2)
                    {
                        // Lk_dense_col = Lk_dense_col*pending_scale/sd(h[ci])
                        //                + L(:,i)*Uiks/L(P(i),i)
                        // use L(P(i),i) instead of sd[i] to avoid scaling
                        // for L(:,i)
                        SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[i]->x[pi], Uiks));
                        SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                                   L->v[i]->x[Ldiag[i]]));

                        SPEX_CHECK(SPEX_mpq_set_ui(tmpq, 1, 1));
                        SPEX_CHECK(SPEX_mpq_set_den(tmpq, sd[h[ci]]));
                        SPEX_CHECK(SPEX_mpq_mul(tmpq, tmpq, pending_scale));

                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], SPEX_MPQ_NUM(tmpq)));
                        SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], SPEX_MPQ_DEN(tmpq)));

                        SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], tmpz));

                        // update h[ci]
                        h[ci] = SPEX_FLIP(i-1);
                        continue;
                    }
                    else if (i > 0)// more efficient
                    {
                        // Lk_dense_col = (Lk_dense_col*sd(i-1)
                        //                                 + L(:,i)*Uiks)/sd(i)
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ci],
                                                Lk_dense_col->x[ci], sd[i-1]));
                    }
                }
                else if (h[ci] >= -1)
                {
                    // this is a fill-in
                    Lk_dense_col->i[Lk_dense_col->nz] = ci;
                    Lk_dense_col->nz++;
                }
                SPEX_CHECK(SPEX_mpz_addmul(Lk_dense_col->x[ci],
                                           L->v[i]->x[pi], Uiks));
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ci],
                                             Lk_dense_col->x[pi], sd[i]));

                // update h[ci]
                h[ci] = SPEX_FLIP(i-1);
            }

            // move Uiks (scaled U(i,Q(ks)) to Lk_dense_col
            SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col->x[P[i]], Uiks));
            Lk_dense_col->i[Lk_dense_col->nz] = P[i];
            Lk_dense_col->nz++;
            // update corresponding entry in the history vector
            h[P[i]] = SPEX_FLIP(i-1);

            // get next nnz
            pks --;
            if (pks < 0) { break;}
            i = Uci[pks];         // row index
        }

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // 1. Iterate across all nnz in Lk_dense_col, perform history update if
        //    needed, then move all nonzero entry from Lk_dense_col to
        //    L->v[k]->x
        // NOTE: explicit zero due to exact cancellation in backtracking
        //       is removed.
        // 2. Swap values from L->v[ks]->x and Lk_dense_col
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // reallocate the nonzero pattern if needed
        if (L->v[k]->nzmax < Lk_dense_col->nz)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], Lk_dense_col->nz));
        }
        pk = 0;
        for (pks = 0; pks < Lk_dense_col->nz; pks++) 
        {
            cks = Lk_dense_col->i[pks]; 
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            h[cks] = SPEX_FLIP(h[cks]);
            if (sgn != 0)
            {
                // check if need to perform history update
                if (h[cks] != k-1)
                {
                    if (k > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                                   Lk_dense_col->x[cks], sd[k-1]));
                    }
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                               Lk_dense_col->x[cks], sd[h[cks]]));
                }
                if (cks == P[k]) // find the diagnal of L(P[k],k)
                {
                    SPEX_CHECK(SPEX_mpz_set(sd[k], Lk_dense_col->x[cks]));
                    Ldiag[k] = pk;
                }
                L->v[k]->i[pk] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[cks]));
                pk++;
            }
        }
        // update number of nnz in column k of L
        L->v[k]->nz = pk;
    }
    // update S(1,k) and S(2,k) as 1, since all entry in L(:,k) are scaled
    SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 1, k), one));
    SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 3, k), one));

    //-------------------------------------------------------------------------
    // swap values from L->v[ks]->x and Lk_dense_col
    //-------------------------------------------------------------------------
    pk = 0;
    for (pks = 0; pks < L->v[ks]->nz; pks++)
    {
        cks = L->v[ks]->i[pks];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[ks]->x[pks]));
        if (sgn != 0 && P_inv[cks] < *inext && P_inv[cks] > ks)
        {
            *inext = P_inv[cks];
        }
        SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col->x[cks], L->v[ks]->x[pks]));
        Lk_dense_col->i[pk] = cks;
        pk++;
    }
    Lk_dense_col->nz = pk;

    //-------------------------------------------------------------------------
    // copy entries to U(k,Q(k+1:ks-1))
    //-------------------------------------------------------------------------
    if (Uk_dense_row->nz-1 > U->v[k]->nzmax)
    {
        SPEX_CHECK(SPEX_vector_realloc(U->v[k], Uk_dense_row->nz-1));
    }
    int64_t Uk_nz = 0;
    pk = 0;
    while (pk < Uk_dense_row->nz)
    {
        j = Uk_dense_row->i[pk];
        // remove explicit zeros
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[j]));
        if (sgn == 0)
        {
            Uk_dense_row->nz--;
            Uk_dense_row->i[pk] = Uk_dense_row->i[Uk_dense_row->nz];
            continue;
        }

        if (Q_inv[j] == k || Q_inv[j] >= ks)
        {
            pk++;
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[Uk_nz], Uk_dense_row->x[j]));
            U->v[k]->i[Uk_nz] = j;
            Uk_nz++;
            Uk_dense_row->nz--;
            Uk_dense_row->i[pk] = Uk_dense_row->i[Uk_dense_row->nz];
        }
    }

    //-------------------------------------------------------------------------
    // RwSOP. This could cause fill-ins!
    //-------------------------------------------------------------------------
    // pending_scale = U(k, Q(ks))/ U(k, Q(k))
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, Uk_dense_row->x[Q[ks]]));
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    int64_t count = 0;
    // iterate for all nnz in U(k+1:ks,Q(ks))
    for (pks = Ucp_k_ks+1; pks < Ucp[Q[ks]+1]; pks++)
    {
        i = Uci[pks];         // row index

        // skip scaling for frames between iterations
        // REMARK: Uci[pks-1] must be at least equal to k
        for (int64_t i1 = Uci[pks-1]+1; i1 < i; i1++)
        {
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, i1),
                                    SPEX_2D(S, 3, i1), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i1], sd[i1],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i1], sd[i1],
                                    SPEX_MPQ_NUM(pending_scale)));
        }

        // skip scaling if U(i, Q(ks)) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx[pks]]));
        if (sgn == 0)
        {
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, i),
                                    SPEX_2D(S, 3, i), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                    SPEX_MPQ_NUM(pending_scale)));
            continue;
        }

        for (pk = 0; pk < Uk_dense_row->nz; pk++)
        {
            j = Uk_dense_row->i[pk];
            h[j] = -2;
        }
        for (pi = 0; pi < U->v[i]->nz; pi++)
        {
            ci = U->v[i]->i[pi];
            if (Q_inv[ci] == ks)
            {
                // handle U(i, Q[ks]) after iteration
                continue;
            }
            else if (Q_inv[ci] > ks)
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[ci]));
            }
            else
            {
                sgn = 0;
            }

            if (sgn == 0)
            {
                // apply pending_scale to U(i,Q(i:ks-1))
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                        SPEX_MPQ_NUM(pending_scale)));
                // set sd[i] = U(i, Q(i))
                if (ci == Q[i])
                {
                    SPEX_CHECK(SPEX_mpz_set(sd[i], U->v[i]->x[pi]));
                }
            }
            else
            {
                // perform RwSOP to U(i,Q(ks+1:n+1))
                // U(i,ci)= (U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci))/U(k,Q(k))
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row->x[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_submul(U->v[i]->x[pi], U->v[i]->x[Ucx[pks]],
                                        Uk_dense_row->x[ci]));
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row->x[Q[k]]));
                count ++;
                h[ci] = -1;
            }
        }

        // CPPU shortcuts
        pi = Ucx[pks];  // pointer for U(i,Q(ks))
        int64_t Ui_nz = U->v[i]->nz;
        if (i != ks)
        {
            // Mathematically, we should flip sign for U(i, Q(ks)). However,
            // since column ks of U will become column k after permutation,
            // which will be deleted when finished, we will instead delete
            // U(i,Q(ks)).
            if (count < Uk_dense_row->nz-2)
            {
                // store U(i,Q(ks)) if there is fillin to be added
                SPEX_CHECK(SPEX_mpz_neg(Uiks, U->v[i]->x[pi]));
            }
            Ui_nz--;
            U->v[i]->i[pi] = U->v[i]->i[Ui_nz];
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[pi], U->v[i]->x[Ui_nz]));
            // skip scaling for L(:,i) by setting S(1,i) = S(1,i)*pending_scale
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, i),
                                    SPEX_2D(S, 1, i), pending_scale));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_neg(U->v[i]->x[pi], U->v[i]->x[pi]));
            SPEX_CHECK(SPEX_mpq_neg(SPEX_2D(S, 1, i), SPEX_2D(S, 1, i)));
            SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
        }

        // finish RwSOP by checking if there is fillin that should be added
        if (count < Uk_dense_row->nz-2)
        {
            // allocate additional space if needed
            int64_t num_of_fillin = Uk_dense_row->nz-2-count;
            if (U->v[i]->nzmax < Ui_nz+num_of_fillin)
            {
                SPEX_CHECK(SPEX_vector_realloc(U->v[i], Ui_nz+num_of_fillin));
            }
            // add FILLIN
            for (pk = 0; count < Uk_dense_row->nz-2 && pk < Uk_dense_row->nz;)
            {
                ci = Uk_dense_row->i[pk];
                if (h[ci] == -2)
                {
                    U->v[i]->i[Ui_nz] = ci;
                    // U(i,ci)= -U(i,Q(ks))*U(k,ci)/U(k,Q(k))
                    //        = Uiks       *U(k,ci)/U(k,Q(k))
                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[Ui_nz], Uiks,
                        Uk_dense_row->x[ci]));
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[Ui_nz],
                        U->v[i]->x[Ui_nz], Uk_dense_row->x[Q[k]]));
                    Ui_nz++;
                    h[ci] = -1;
                    count++;
                }
                pk++;
            }
        }
        U->v[i]->nz = Ui_nz;
    }

    //-------------------------------------------------------------------------
    // finish copying the remaining entries in U(k,:)
    // and copy U(ks,:) to Uk_dense_row
    //-------------------------------------------------------------------------
    for (pk = 0; pk < Uk_dense_row->nz; pk++)
    {
        j = Uk_dense_row->i[pk];
        if (j == Q[k])
        {
            // no need to copy U(k,Q[k])
            SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[j], 0));
            continue;
        }

        SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[Uk_nz], Uk_dense_row->x[j]));
        U->v[k]->i[Uk_nz] = j;
        Uk_nz++;
    }
    U->v[k]->nz = Uk_nz;

    // no need to do this when ks>=n-1
    if (ks < n-1)
    {
        // construct a scattered vector for U->v[ks]
        for (pks = 0; pks < U->v[ks]->nz; pks++)
        {
            j = U->v[ks]->i[pks];
            if (j == Q[ks])
            {
                // this entry should be considered as the IPGE update of the
                // entry in the Q[k]-th column
                j = Q[k];
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[ks]->x[pks]));
                if (sgn != 0 && Q_inv[j] < *jnext)
                {
                    *jnext = Q_inv[j];
                }
            }
            Uk_dense_row->i[pks] = j;
            SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[j], U->v[ks]->x[pks]));
        }
        Uk_dense_row->nz = U->v[ks]->nz;
    }

    //-------------------------------------------------------------------------
    // update column permutation
    //-------------------------------------------------------------------------
    int64_t tmpi;
    tmpi = Q[k];      Q[k] = Q[ks];          Q[ks] = tmpi;
    Q_inv[Q[k]] = k;  Q_inv[Q[ks]] = ks;

    //-------------------------------------------------------------------------
    // flip sign for columns and rows ks+1 to n and update corresponding sd
    //-------------------------------------------------------------------------
    for (i = ks+1; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpq_neg(SPEX_2D(S, 3, i), SPEX_2D(S, 3, i)));
        SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
