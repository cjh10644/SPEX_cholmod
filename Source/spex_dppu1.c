//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_dppu1.c: perform diagonal permutation pivot update
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform diagonal permutation pivot update
// when the submatrix (formed by rows and columns k to ks) has the following
// pattern
//       x 0 0 0 0
//       0 x . . 0
//       0 . x . 0
//       0 . . x 0
//       . 0 0 0 x
// This function will swap rows and columns k and ks in L and U. Noted that the
// rows of L and columns of U are permuted implicitly via the permutation
// matrices based on P and Q.

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Lksk);            \
    SPEX_MPQ_CLEAR(tmp_mpq);         \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);         \
    SPEX_FREE(h);

#include "spex_internal.h"

SPEX_info spex_dppu1
(
    SPEX_matrix *L,
    SPEX_matrix *U,
    SPEX_matrix *S
    mpz_t *d
    mpz_t *sd,
    mpz_t *Lk_dense_col,
    mpz_t *Uk_dense_row,
    int64_t LksxPks,         // L->v[ks]->x[LksxPks] = L(P(ks),ks)
    int64_t *P,
    int64_t *Q,
    int64_t k,             // current column index 0 <= k < n
    int64_t ks,            // index of the diagonal to be swapped with, [0,n)
    bool Lksk_IS_zero      // simply check if L_row_offdiag(ks)<k
)
{

    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t pk, ck, pks, cks, tmp_int;
    spex_vector *v;
    bool *h;  // used to mark whether a entry in Lk_dense_col is updated

    mpq_t pending_scale, tmp_mpq;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(tmp_mpq);
    mpz_t Lksk, tmp_mpz; SPEX_MPZ_SET_NULL(Lksk); SPEX_MPZ_SET_NULL(tmp_mpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(tmp_mpq));
    SPEX_CHECK(SPEX_mpz_init(Lksk));
    SPEX_CHECK(SPEX_mpz_init(tmp_mpz));

    h = (bool*) SPEX_calloc(n, sizeof(bool));
    if (h == NULL)
    {
        SPEX_FREE_WORK;
        return SPEX_OUT_OF_MEMORY;
    }
    // -------------------------------------------------------------------------
    // handle the special case when swapping with the inserted column. Since it
    // is only in the k-th IPGE iteration, there is no need to perform
    // backtracking for the inserted column. Therefore, only need to perform
    // backtracking for U(n-1, n-1).
    // NOTE: U(k,ks) must be nnz, otherwise it would cause singularity, since
    //       U(k,k+1:n) will be all zeros
    // -------------------------------------------------------------------------
    if (ks == n)
    {
        // backtrack U(n-1,n-1)
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpz_mul(U->v[n-1]->x[0], sd[n-1], sd[k-1]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_set(U->v[n-1]->x[0], sd[n-1]));
        }
        if (Lksk_IS_zero)
        {
            SPEX_CHECK(SPEX_mpz_divexact(U->v[n-1]->x[0],
                                         U->v[n-1]->x[0], sd[n-2]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[n-1]->x[0]
                                       U->v[n-1]->x[0], sd[n-2]));
            // tmp_mpq = S(1,k)*S(3,k)
            SPEX_CHECK(SPEX_mpq_mul(tmp_mpq, SPEX_2D(S, 1, k),
                                    SPEX_2D(S, 3, k) ));
            // assign L(P(n-1),k) = L(P(n-1),k)*tmp_mpq, which should be integer
            SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[P[n-1]],
                                   Lk_dense_col[P[n-1]],SPEX_MPQ_DEN(tmp_mpq)));
            SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[P[n-1]],
                                   Lk_dense_col[P[n-1]],SPEX_MPQ_NUM(tmp_mpq)));
            // tmp_mpz = ceil(U(k,Q(n-1))*L(P(n-1),k)/U(k,Q(k))
            SPEX_CHECK(SPEX_mpz_mul(tmp_mpz, Uk_dense_row[Q[n-1]],
                                    Lk_dense_col[P[n-1]]));
            SPEX_CHECK(SPEX_mpz_cdiv_q(tmp_mpz,
                                       tmp_mpz, Uk_dense_row[Q[k]]));
            // U(n-1,n-1) = U(n-1,n-1)+tmp_mpz
            SPEX_CHECK(SPEX_mpz_add(U->v[n-1]->x[0], U->v[n-1]->x[0], tmp_mpz));
        }

        // ---------------------------------------------------------------------
        // scale entries in frames k+1:n-2
        // ---------------------------------------------------------------------
        if (n > k+2) // n-1 > k+1
        {
            // get the scale for entries between frames k and n-1 % O(1) time
            // pending_scale = vk(P[n-1])/sd(k);
            // TODO make sure vk[P[n-1]] is scaled
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, vk[P[n-1]]));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

            for (j = k+1; j < n-1; j++)
            {
                // S(3,k+1:n-2) = S(3,k+1:n-2)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                        SPEX_2D(S, 3, j), pending_scale));
                // sd(k+1:n-2) = sd(k+1:n-2)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }
        }
        // ---------------------------------------------------------------------
        // set d, sd and S for frame n-1
        // sd[k], d[k] and S(:,k) will be updated when vk is inserted
        // ---------------------------------------------------------------------
        // since U(k,n-1) will be U(n-1,n-1) after permutation, we will keep
        // the scaling factor for frame n-1 same as that for row k of U, while
        // multiply with the scaling factor due to the IPGE update.
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, n-1),
                                SPEX_2D(S, 3, k), SPEX_2D(S, 2, k)));
        // get the scale for IPGE update: pending_scale = sd(n-2)/sd(k-1);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[n-2]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, n-1),
                                SPEX_2D(S, 3, n-1), pending_scale));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, n-1), 1, 1));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, n-1), 1, 1));

        // U(k,n-1) will be U(n-1,n-1)
        // TODO set U(n-1,n-1)=L(n-1,n-1)=d[n-1]
        SPEX_CHECK(SPEX_mpz_set(d[n-1],  Uk_dense_row[Q[n-1]]));
        SPEX_CHECK(SPEX_mpz_divexact(sd[n-1],
                                 d[n-1], SPEX_MPQ_DEN(SPEX_2D(S, 3, n-1))));
        SPEX_CHECK(SPEX_mpz_mul(sd[n-1],
                                sd[n-1], SPEX_MPQ_NUM(SPEX_2D(S, 3, n-1))));

        // ---------------------------------------------------------------------
        // Columns of k and n will be swapped after calling this function, we
        // only need to swap rows of k and n-1
        // ---------------------------------------------------------------------
        // update the nnz pattern of column k of U, Uk_x and Uk_i
        Uk_i[0] = ks; // k-th row will become row ks after swap

        // swap rows k and n-1 of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[n-1];    U->v[n-1] = v;

        // update row permutation to swap rows of L implicitly
        tmp_int = P[k];    P[k] = P[n-1];          P[n-1] = tmp_int;

        SPEX_FREE_WORK;
        return SPEX_OK;
    }

    // -------------------------------------------------------------------------
    // perform backtracking
    // -------------------------------------------------------------------------
    // find the scale for backtracking: pending_scale = sd(k-1)/sd(ks-1)
    if (k == 0)
    {
        SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
    }
    else
    {
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k-1]));
    }
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks-1]));
    // remove common factor in mpq_den and mpq_num
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    if (Lksk_IS_zero) // L(P(ks),k) == 0
    {
        // S(3, ks) = pending_scale*S(3,ks)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                                SPEX_2D(S, 3, ks), pending_scale));

        // sd(ks) = sd(ks)*pending_scale
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks],
                                     sd[ks], SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));
    }
    else
    {
        // tmp_mpq = S(1,k)*S(3,k)
        SPEX_CHECK(SPEX_mpq_mul(tmp_mpq, SPEX_2D(S, 1, k), SPEX_2D(S, 3, k) ));
        // assign Lksk = L(P(ks),k)*tmp_mpq, which should be integer
        SPEX_CHECK(SPEX_mpz_divexact(Lksk, Lk_dense_col[P[ks]],
                                     SPEX_MPQ_DEN(tmp_mpq)));
        SPEX_CHECK(SPEX_mpz_mul(Lksk, Lksk, SPEX_MPQ_NUM(tmp_mpq)));

        // update the scale for col ks of L due to backtracking
        // S(1,ks) = tmp_mpq*pending_scale = S(1,ks)*S(3,ks)*sd(k-1)/sd(ks-1)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, ks), tmp_mpq, pending_scale));

        // pending_scale = pending_scale*S(2,ks)*S(3,ks)
        //               = S(2,ks)*S(3,ks)*sd(k-1)/sd(ks-1)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                pending_scale, SPEX_2D(S, 2, ks)));
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                pending_scale, SPEX_2D(S, 3, ks)));

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        // backtracking jumbled sparse row ks of U using scattered row k of U.
        // Assuming explicit zeros in U(ks, :) resulted from exact cancellation
        // in IPGE update were not removed (SLIP LU keeps those zeros in output
        // L and U), nonzero pattern of U(k,Q(k+1:n+1)) should be a subset of
        // U(ks,:).
        //
        // REMARK:
        // when the IPGE update results in exact cancellation and the resulted
        // zero is removed from L or U, new entry should be inserted to row ks
        // after backtracking. In this case, we need to iterate across all
        // nonzeros in Uk_dense_row to find if any column index of nonzero is
        // untouched, then a new nonzero should be added.
        // U(ks,cks) = U(k,cks)*Lksk/U(k,Q(k))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        pks = 0;
        while (pks < U->v[ks]->nz)
        {
            // column index in row ks of U
            cks = U->v[ks]->i[pks];

            // U(ks,cks) = U(ks,cks)*pending_scale+U(k,cks)*Lksk/U(k,Q(k))
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row[cks]));
            if (sgn != 0)  // U(k,cks) != 0
            {
                // tmp_mpz = ceil(U(k,cks)*Lksk/U(k,Q(k))
                SPEX_CHECK(SPEX_mpz_mul(tmp_mpz, Uk_dense_row[cks], Lksk));
                SPEX_CHECK(SPEX_mpz_cdiv_q(tmp_mpz,
                                        tmp_mpz, Uk_dense_row[Q[k]]));
                // U(ks,cks) = floor(U(ks,cks)*pending_scale)
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(pending_scale)));
                SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_DEN(pending_scale)));
                // U(ks,cks) = U(ks,cks)+tmp_mpz
                SPEX_CHECK(SPEX_mpz_add(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        tmp_mpz));

                // remove this entry if it becomes zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[ks]->x[pks]));
                if (sgn == 0)
                {
                    U->v[ks]->nz --;
                    U->v[ks]->i[pks] = U->v[ks]->i[U->v[ks]->nz];
                    continue;
                }
            }
            else           // U(k,cks) == 0
            {
                // U(ks,cks) = U(ks,cks)*pending_scale
                SPEX_CHECK(SPEX_mpz_divexact(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(pending_scale)));

                // update sd[ks] = U(ks, Q(ks))
                if (cks == Q[ks])
                {
                    SPEX_CHECK(SPEX_mpz_set(sd[ks], U->v[ks]->x[pks]));
                }
            }
            pks ++;
        }

        // d[ks] = L(P(ks), ks)
        SPEX_CHECK(SPEX_mpz_set(d[ks], L->v[ks]->x[Ldiag[ks]]));
        // reset S(2,ks) and S(3,ks)
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, ks), 1, 1));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, ks), 1, 1));

        // Mathematically, we should insert new entry at U(ks, Q[k]) and swap
        // its value with Lksk. However, since the value of this entry will not
        // be used beyond this point, and column Q[k] of U will be deleted when
        // finished, we will skipped adding U(ks, Q[k]) here.
        // On the other hand, set Lksk = L(P(ks), k), which has same pending
        // scale factor as the rest of entries in column k of L, so that the
        // skipped scaling for col k of L can still be skipped when performing
        // IPGE update
        SPEX_CHECK(SPEX_mpq_set(Lksk, Lk_dense_col[P[ks]]));
    }
    
    // ------------------------------------------------------------------------
    // swap rows and columns of k and ks
    // ------------------------------------------------------------------------
    // swap columns k and ks of L        % O(1) time
    v = L->v[k];       L->v[k] = L->v[ks];    L->v[ks] = v;
    // swap rows k and ks of U           % O(1) time
    v = U->v[k];       U->v[k] = U->v[ks];    U->v[ks] = v;

    // update row/column permutation to swap rows of L and cols of U implicitly
    tmp_int = Q[k];    Q[k] = Q[ks];          Q[ks] = tmp_int;
    tmp_int = P[k];    P[k] = P[ks];          P[ks] = tmp_int;

    // swap entries in d, sd and S
    SPEX_CHECK(SPEX_mpz_swap(sd[k], sd[ks]));
    SPEX_CHECK(SPEX_mpz_swap(d[k],  d[ks]));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 1, ks), SPEX_2D(S, 1, k)));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 2, ks), SPEX_2D(S, 2, k)));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 3, ks), SPEX_2D(S, 3, k)));

    // ------------------------------------------------------------------------
    // scale entries in frames k+1:ks-1
    // ------------------------------------------------------------------------
    if (ks > k+1)
    {
        // get the scale for entries between frames k and ks % O(1) time
        // pending_scale = sd(k)/sd (ks);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        for (j = k+1; j < ks; j++)
        {
            // S(3,k+1:ks-1) = S(3,k+1:ks-1)*pending_scale;
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                    SPEX_2D(S, 3, j), pending_scale));
            // sd(k+1:ks-1) = sd(k+1:ks-1)*pending_scale;
            SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                    sd[j], SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                    sd[j], SPEX_MPQ_NUM(pending_scale)));
        }
    }

    // ------------------------------------------------------------------------
    // perform IPGE for frame ks, whose entries are in Lk_dense_col and
    // Uk_dense_row.
    // NOTE: in the last iteration, i.e., swapping column k with last column,
    // the IPGE update column k of L is useless, since the updated column will
    // be deleted and replaced. Therefore, its IPGE update in the last
    // iteration can be treated same as row k of U for better efficiency.
    // ------------------------------------------------------------------------
    if (Lksk_IS_zero || ks == n-1)
    {
        // get the scale for IPGE update: pending_scale = sd(ks-1)/sd (k-1);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }

        // S(3, ks) = S(3, ks)*pending_scale
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                                SPEX_2D(S, 3, ks), pending_scale));
        // sd(ks) = sd(ks)*pending_scale;
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks], sd[ks],
                                     SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));
    }
    else
    {
        // skip scaling for U for 1 IPGE iteration
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        // S(2,ks) = S(2,ks)*pending_scale
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, ks),
                                SPEX_2D(S, 2, ks), pending_scale));

        // initialize history vector
        for (pks = 0; pks < L->v[ks]->nz; pks++)
        {
            cks = L->v[ks]->i[pks];
            // formally, we should set h[cks] = SPEX_FLIP(k-1), so we will know
            // the entries in L(:,ks) are in (k-1)-th IPGE iteration. However,
            // since we need to perform only one IPGE iteration, we just need
            // to know whether the corresponding entry is updated. Therefore,
            // the initialization for history vector is set as
            h[cks] = -2; // only entry in the nnz patter has h < -1
        }
        SPEX_ASSERT(pks == L->v[ks]->nz);

        // NOTE: this will cause fillin in the ks(th) column of L
        //       There is no subset relation between nnz pattern in L(:,ks)
        //       and L(:, k). Both could have explicit zero(s).
        //       L(:,k) can be jumbled.
        for (pk = 0; pk < L->v[k]->nz; pk++)
        {
            // row index in column k of L
            ck = L->v[k]->i[pk];
            // exclude L(P[k],k)
            if (ck == P[k])
            {
                continue;
            }

            // L(ck,ks) = (L(ck,ks)*d(k)-L(ck,k)*Lksk)/sd(k-1);
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[ck]));
            if (sgn != 0) // L(ck, ks) != 0
            {
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[ck],
                                        Lk_dense_col[ck], d[k]));
            }
            else if (h[ck] >= -1) // this entry was not in nnz pattern
            {
                // TODO: check if this will 0?
                // check if this will be the first off diagonal entry in L(P,ks)
                if (P_inv[ck] < inext)
                {
                    // inext is the row index of the found first off-diagonal
                    // entry in L(P,ks)
                    inext = P_inv[ck];
                }

                // reallocate the nonzero pattern if needed
                if (pks == L->v[ks]->max_nnz)
                {
                    SPEX_CHECK(spex_expand(&(L->v[ks])));
                }
                // insert new entry in the nonzero pattern
                L->v[ks]->i[pks] = ck;
                pks++;
            }
            SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col[ck],
                                       L->v[k]->x[pk], Lksk));
            if (k > 0)
            {
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[ck],
                                             Lk_dense_col[ck], sd[k-1]));
            }
            
            // update h[ck] to mark Lk_dense_col[ck] need no further update
            h[ck] = -1;
        }
        // update the number of nnz with pks-1 since one entry will be deleted
        L->v[ks]->nz = pks-1;

        for (pks = 0; pks < L->v[ks]->nz; pks ++)
        {
            // row index in column ks of L
            cks = L->v[ks]->i[pks];

            if (h[cks] < -1) //only need to update entries that were not updated
            {
                // L(P(k), ks) should be removed from nnz pattern
                if (cks == P[k])
                {
                    // move the row index of last nonzero to current position
                    cks = L->v[ks]->i[L->v[ks]->nz];
                    L->v[ks]->i[pks] = cks;
                    if (h[cks] >= -1) { continue; }
                }
                // L(ck,ks) = (L(ck,ks)*d(k))/sd(k-1);
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[cks],
                                        Lk_dense_col[cks], d[k]));
                if (k > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[cks],
                                                 Lk_dense_col[cks], sd[k-1]));
                }
                // reset history vector to any value >= -1
                h[cks] = -1;
            }
        }

        // S(1,ks) = S(1,ks)*S(1,k)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, ks),
                                SPEX_2D(S, 1, ks), SPEX_2D(S, 1, k)));
        
        // skip the rest of IPGE iterations
        // pending_scale = sd(ks-1)/sd(k);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        // S(3, ks) = S(3, ks)*pending_scale;
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                                SPEX_2D(S, 3, ks), pending_scale));
        // sd(ks) = L(P(ks),ks)*S(1,ks)*S(3,ks);
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                SPEX_2D(S, 3, ks), SPEX_2D(S, 1, ks)));
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks], sd[ks],
                                SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], Lk_dense_col[P[ks]],
                                SPEX_MPQ_NUM(pending_scale)));
    }
    SPEX_FREE_WORK;
    return SPEX_OK;
}
