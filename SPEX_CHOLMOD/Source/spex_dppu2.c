//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_dppu2.c: perform diagonal permutation pivot update
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform diagonal permutation pivot update
// when the submatrix (formed by columns k to ks) has the following
// pattern
//       . . . . .  (row 1)
//       . . . . .  (rows 2 to k-1)
//       x . . . .  (row k)
//       0 x . . 0  (row k+1)
//       0 . x . 0  ( .... )
//       0 . . x 0  (row ks-1)
//       0 0 0 0 x  (row ks)
//       0 . . . .  (row ks+1 to n-1)
//       0 . . . .  (row n)
// This function will swap rows and columns k and ks in L and U. Noted that the
// rows of L and columns of U are permuted implicitly via the permutation
// matrices based on P and Q.

#define SPEX_FREE_ALL                \
    SPEX_MPZ_CLEAR(Lksk);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);

#include "spex_internal.h"

SPEX_info spex_dppu2
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
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                       // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
)
{
    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t tmpi, j, n = U->n, tmp_ks = SPEX_MIN(ks, n-1);
    SPEX_vector *v;

    mpq_t pending_scale, one;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(one);
    mpz_t Lksk; SPEX_MPZ_SET_NULL(Lksk);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Lksk));

    if (ks == n)
    {
        //----------------------------------------------------------------------
        // backtrack U(n-1,n-1)
        //----------------------------------------------------------------------
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpz_mul(U->v[n-1]->x[0], sd[n-1], sd[k-1]));
            SPEX_CHECK(SPEX_mpz_divexact(U->v[n-1]->x[0],
                                         U->v[n-1]->x[0], sd[n-2]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_divexact(U->v[n-1]->x[0], sd[n-1], sd[n-2]));
        }

        //----------------------------------------------------------------------
        // update entries in frames between k and n-1
        //----------------------------------------------------------------------        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping columns and rows of
        // k and n-1. Then we set sd[k] to d[k] with pending scaling factor
        // applied
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Q[k]], sd[k]));
        SPEX_CHECK(SPEX_mpz_divexact(sd[k],  d[k], SPEX_MPQ_DEN(vk_scale)));
        SPEX_CHECK(SPEX_mpz_mul     (sd[k], sd[k], SPEX_MPQ_NUM(vk_scale)));

        if (n > k+2) // n-1 > k+1
        {
            // pending_scale = sd(k)/Uk_dense_row[Q[k]]
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

            // scale entries in frames k+1:n-2
            for (j = k+1; j < n-1; j++)
            {
                // S(3,j) = S(3,j)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                        SPEX_2D(S, 3, j), pending_scale));
                // sd(j) = sd(j)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }
        }

        //----------------------------------------------------------------------
        // swap rows and columns k and n-1 of L and U
        //----------------------------------------------------------------------
        // swap rows k and n-1 of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[n-1];    U->v[n-1] = v;

        // update row permutation to swap rows of L implicitly
        tmpi = P[k];    P[k] = P[n-1];          P[n-1] = tmpi;
        P_inv[P[k]] = k;   P_inv[tmpi] = n-1;

        // U(k,Q(k)) and S(:,k) will be updated after calling this function
        // S(:,n-1) = [1;1;S(3,k)]; since S(:,k) and S(:,n-1) not swapped
        SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 3, n-1), SPEX_2D(S, 3, k))); 
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, n-1), 1, 1));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, n-1), 1, 1));
    }
    else
    {
        //----------------------------------------------------------------------
        // perform backtracking for frame ks
        //----------------------------------------------------------------------
        // find the scale for backtracking
        if (k == 0)
        {
            // pending_scale = 1/sd(ks-1)
            SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
        }
        else
        {
            // pending_scale = sd(k-1)/sd(ks-1)
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k-1]));
        }
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[tmp_ks-1]));
        // remove common factor in mpq_den and mpq_num
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        // S(3,ks) = pending_scale*S(3,ks)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, tmp_ks),
                                SPEX_2D(S, 3, tmp_ks), pending_scale));
        // sd(ks) = sd(ks)*pending_scale
        SPEX_CHECK(SPEX_mpz_divexact(sd[tmp_ks],
                                sd[tmp_ks], SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[tmp_ks],
                                sd[tmp_ks], SPEX_MPQ_NUM(pending_scale)));

        //----------------------------------------------------------------------
        // swap rows and columns k and ks of L and U
        //----------------------------------------------------------------------
        // swap entries in d, sd and S
        SPEX_CHECK(SPEX_mpz_swap(sd[k], sd[tmp_ks]));
        SPEX_CHECK(SPEX_mpz_swap(d[k],  d[tmp_ks]));
        SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 1, tmp_ks), SPEX_2D(S, 1, k)));
        SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 2, tmp_ks), SPEX_2D(S, 2, k)));
        SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 3, tmp_ks), SPEX_2D(S, 3, k)));

        // swap columns k and ks of L        % O(1) time
        v = L->v[k];       L->v[k] = L->v[tmp_ks];    L->v[tmp_ks] = v;
        // swap rows k and ks of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[tmp_ks];    U->v[tmp_ks] = v;

        // update row/col permutation to swap rows of L and cols of U implicitly
        tmpi = Q[k];       Q[k] = Q[tmp_ks];          Q[tmp_ks] = tmpi;
        Q_inv[Q[k]] = k;   Q_inv[tmpi] = tmp_ks;
        tmpi = P[k];       P[k] = P[tmp_ks];          P[tmp_ks] = tmpi;
        P_inv[P[k]] = k;   P_inv[tmpi] = tmp_ks;

        // update Ldiag[k] = Ldiag[ks]
        Ldiag[k] = Ldiag[tmp_ks];

        //----------------------------------------------------------------------
        // update entries in frames between k and ks
        //----------------------------------------------------------------------
        if (tmp_ks > k+1)
        {
            // get the scale for entries between frames k and ks % O(1) time 
            // pending_scale = sd(k)/sd (ks); 
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k])); 
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[tmp_ks])); 
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
            // scale entries in frames k+1:ks-1
            for (j = k+1; j < tmp_ks; j++)
            {
                // S(3,j) = S(3,j)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                        SPEX_2D(S, 3, j), pending_scale));
                // sd(j) = sd(j)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }
        }
    }

    //--------------------------------------------------------------------------
    // perform IPGE for row ks, skip IPGE for column since it is all zero
    //--------------------------------------------------------------------------
    int64_t pks, cks, last_nz_b4_ks = k-1;
    // initialize history vector h
    for (pks = 0; pks < Uk_dense_row->nz; pks++)
    {
        cks = Uk_dense_row->i[pks];
        // Lk_dense_col or Uk_dense_row are initialized with no explicit zero
        // for column/row 0 (may contain explicit zero for column/row j>0).  And
        // entries in h are maintained to be >= -1. Therefore, with such
        // initialization, entry with h > -1 is clearly not in nnz pattern and
        // any entry in the nnz pattern with h = -1 must be nonzero.  In all,
        // any explicit zero with h >= -1 must not be in the nnz pattern.
        //
        // REMARK:
        // This is useful only for IPGE update, or when there is possibility of
        // fillin.
        h[cks] = SPEX_FLIP(k-1);
    }
    for (j = k; j < tmp_ks; j++)
    {
        if (ks == n && j == k)
        {
            // find the pointer to the entry L(P[ks],k), there are only two
            // entries in L->v[k]
            pks = (P[L->v[k]->i[0]] == tmp_ks) ? 0 : 1;
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[k]->x[pks]));
            if (sgn == 0) { continue; }

            // only need to perform IPGE for U(ks, Q(ks)) since there is only
            // one off-diagonal nnz in row k of U, which is U(k,Q(ks))
            // U(ks,Q(ks)) = (U(ks, Q(ks))*d[k]-
            //                                     vk(P[ks])*U(k,Q(ks)))/sd[k-1]
            SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row->x[Q[tmp_ks]],
                                    Uk_dense_row->x[Q[tmp_ks]], d[k]));
            SPEX_CHECK(SPEX_mpz_addmul(Uk_dense_row->x[Q[tmp_ks]],
                                    U->v[k]->x[0], L->v[k]->x[pks]));
            if (k > 0)
            {
                SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row->x[Q[tmp_ks]],
                                    Uk_dense_row->x[Q[tmp_ks]], sd[k-1]));
            }
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, tmp_ks),
                                    SPEX_2D(S, 3, tmp_ks), vk_scale));

            // update history vector
            h[Q[tmp_ks]] = SPEX_FLIP(k);
            continue;
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[Q[j]]));
        }
        // skip if U(ks, Q[j]) == 0
        if (sgn == 0) { continue; }

        // perform j-th IPGE update for U(ks,:)
        SPEX_CHECK(spex_ipge(Uk_dense_row, SPEX_2D(S, 3, tmp_ks), h, NULL,
            U->v[j], Q, Q_inv, (const mpz_t*)sd, SPEX_2D(S, 2, j),
            SPEX_2D(S, 3, j), Ucx[Ucp[Q[j]+1]-1], j));
        // update last_nz_b4_ks
        last_nz_b4_ks = j;

        // insert new entry L(P(ks), j) to L and swap its value with U(ks, Q(j))
        SPEX_CHECK(spex_insert_new_entry(Uk_dense_row->x[Q[j]], L->v[j],
            SPEX_2D(S, 1, j), U->v[j], SPEX_2D(S, 2, j), SPEX_2D(S, 3, j), d[j],
            P[tmp_ks], Ucx[Ucp[Q[j]+1]-1], one));
    }

    // There must be at least one nonzero in U(ks, Q(k:ks-1)) for this case.
    // Otherwise, this will be handled by spex_dppu1(...). Therefore,
    // last_nz_b4_ks must have updated at least once and be greater than k-1
    ASSERT(last_nz_b4_ks > k-1);

    if (ks == n)
    {
        h[Q[k]]  = SPEX_FLIP(h[Q[k]]);
        h[Q[tmp_ks]] = SPEX_FLIP(h[Q[tmp_ks]]);
        last_nz_b4_ks = h[Q[tmp_ks]];
    }
    else
    {
        // perform history update up to (last_nz_b4_ks)-th IPGE iteration
        // and remove zero from row ks of U
        pks = 0;
        *jnext = n;
        while (pks < Uk_dense_row->nz)
        {
            // column index in row ks of U
            cks = Uk_dense_row->i[pks];
            h[cks] = SPEX_FLIP(h[cks]);
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[cks]));
            if (sgn == 0) {continue;}

            if (h[cks] < last_nz_b4_ks) // require history update
            {
                // U(ks, Q(k:ks-1)) were set to 0 but remained in nnz pattern
                if (Q_inv[cks] < tmp_ks)
                {
                    // move the column index of last nonzero to current position
                    Uk_dense_row->nz--;
                    Uk_dense_row->i[pks] = Uk_dense_row->i[Uk_dense_row->nz];
                    continue;
                }
                else if (Q_inv[cks] > tmp_ks && Q_inv[cks] < *jnext)
                {
                    // update the index of next off-diagonal nnz entry
                    *jnext = Q_inv[cks];
                }

                // U(ks,cks) = (U(ks,cks)*sd(last_nz_b4_ks))/sd(h[cks]);
                SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row->x[cks],
                                      Uk_dense_row->x[cks], sd[last_nz_b4_ks]));
                if (h[cks] > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row->x[cks],
                                     Uk_dense_row->x[cks], sd[h[cks]]));
                }
            }
            pks++;
        }
    }

    // update S(3,ks)= S(3,ks)*sd(ks-1)/sd(last_nz_b4_ks)
    if (last_nz_b4_ks != tmp_ks-1)
    {
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[tmp_ks-1])); 
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[last_nz_b4_ks])); 
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, tmp_ks),
                                SPEX_2D(S, 3, tmp_ks), pending_scale));
    }

    // d(ks)       = U(ks,Q(ks));
    SPEX_CHECK(SPEX_mpz_set(d[tmp_ks], Uk_dense_row->x[Q[tmp_ks]]));
    // sd(ks)      = U(ks,Q(ks))*S(3,ks)
    SPEX_CHECK(SPEX_mpz_divexact(sd[tmp_ks],
                                  d[tmp_ks],SPEX_MPQ_DEN(SPEX_2D(S,3,tmp_ks))));
    SPEX_CHECK(SPEX_mpz_mul     (sd[tmp_ks],
                                 sd[tmp_ks],SPEX_MPQ_NUM(SPEX_2D(S,3,tmp_ks))));

    // no need to update L(P(ks), ks), which will be updated in the
    // last iteration
    //SPEX_CHECK(SPEX_mpz_set(Lk_dense_col[P[tmp_ks]],Uk_dense_row[Q[tmp_ks]]));
    if (ks == n)
    {
        // move data from Uk_dense_row, there is only one entry that needs
        // to move, which is U(k,Q[n-1])
        // set U(n-1,n-1)=L(n-1,n-1)=Uk_dense_row[Q[n-1]]
        SPEX_CHECK(SPEX_mpz_swap(U->v[n-1]->x[0], Uk_dense_row->x[Q[n-1]]));
        SPEX_CHECK(SPEX_mpz_set (L->v[n-1]->x[0], U->v[n-1]->x[0]     ));
        U->v[n-1]->i[0] = Q[n-1];
        U->v[n-1]->nz = 1;
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}