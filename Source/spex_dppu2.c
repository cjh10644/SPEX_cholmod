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

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Lksk);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);         \
    SPEX_FREE(h);

#include "spex_internal.h"

SPEX_info spex_dppu2
(
    SPEX_matrix *L,
    SPEX_matrix *U,
    SPEX_matrix *S
    mpz_t *d
    mpz_t *sd,
    mpz_t *Lk_dense_col,
    mpz_t *Uk_dense_row,
    int64_t *Udiag,         // U->v[k]->x[Udiag[k]] = U(k,Q(k))
    int64_t LksxPks,         // L->v[ks]->x[LksxPks] = L(P(ks),ks)
    int64_t *P,
    int64_t *Q,
    int64_t k,             // current column index 0 <= k < n
    int64_t ks             // index of the diagonal to be swapped with, [0,n)
)
{
    // initialize workspace
    SPEX_info info;
    int64_t real_h;

    mpq_t pending_scale, one;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(one);
    mpz_t Lksk, tmp_mpz; SPEX_MPZ_SET_NULL(Lksk); SPEX_MPZ_SET_NULL(tmp_mpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Lksk));
    SPEX_CHECK(SPEX_mpz_init(tmp_mpz));

    h = (int64_t*) SPEX_malloc(n, sizeof(int64_t));
    if (h == NULL)
    {
        SPEX_FREE_WORK;
        return SPEX_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // perform backtracking for frame ks
    //--------------------------------------------------------------------------
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
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks-1]));
    // remove common factor in mpq_den and mpq_num
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    // S(3,ks) = pending_scale*S(3,ks)
    SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                            SPEX_2D(S, 3, ks), pending_scale));
    // sd(ks) = sd(ks)*pending_scale
    SPEX_CHECK(SPEX_mpz_divexact(sd[ks], sd[ks], SPEX_MPQ_DEN(pending_scale)));
    SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));

    //--------------------------------------------------------------------------
    // swap rows and columns k and ks of L and U
    //--------------------------------------------------------------------------
    spex_vector *v;
    // swap columns k and ks of L        % O(1) time
    v = L->v[k];       L->v[k] = L->v[ks];    L->v[ks] = v;
    // swap rows k and ks of U           % O(1) time
    v = U->v[k];       U->v[k] = U->v[ks];    U->v[ks] = v;

    int64_t tmp_int;
    // update row/column permutation to swap rows of L and cols of U implicitly
    tmp_int = Q[k];    Q[k] = Q[ks];          Q[ks] = tmp_int;
    tmp_int = P[k];    P[k] = P[ks];          P[ks] = tmp_int;

    // swap entries in d, sd and S
    SPEX_CHECK(SPEX_mpz_swap(sd[k], sd[ks]));
    SPEX_CHECK(SPEX_mpz_swap(d[k],  d[ks])); 
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 1, ks), SPEX_2D(S, 1, k))); 
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 2, ks), SPEX_2D(S, 2, k))); 
    SPEX_CHECK(SPEX_mpq_swap(SPEX_2D(S, 3, ks), SPEX_2D(S, 3, k)));

    // TODO:? update Udiag and Ldiag
    tmp_int = Udiag[k]; Udiag[k] = Udiag[ks]; Udiag[ks] = tmp_int;
    //tmp_int = Ldiag[k]; Ldiag[k] = Ldiag[ks]; Ldiag[ks] = tmp_int;

    //--------------------------------------------------------------------------
    // update entries in frames between k and ks
    //--------------------------------------------------------------------------
    if (ks > k+1)
    {
        // get the scale for entries between frames k and ks % O(1) time 
        // pending_scale = sd(k)/sd (ks); 
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k])); 
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks])); 
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        // scale entries in frames k+1:ks-1
        for (int64_t j = k+1; j < ks; j++)
        {
            // S(3,j) = S(3,j)*pending_scale;
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, j),
                                    SPEX_2D(S, 3, j), pending_scale));
            // sd(j) = sd(j)*pending_scale;
            SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                         sd[j], SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[j], sd[j], SPEX_MPQ_NUM(pending_scale)));
        }
    }

    //--------------------------------------------------------------------------
    // perform IPGE for row ks, skip IPGE for column since it is all zero
    //--------------------------------------------------------------------------
    // S(:,ks)     = [1;1;S(2,ks)*S(3,ks)];
    SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                            SPEX_2D(S, 3, ks), SPEX_2D(S, 2, ks))); 
    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, ks), 1, 1)); 
    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, ks), 1, 1)); 
    int64_t pks, cks, pj, last_nz_b4_ks = k-1;
    // initialize history vector h
    for (pks = 0; pks < U->v[ks]->nz; pks++)
    {
        cks = U->v[ks]->i[pks];
        // Lk_dense_col or Uk_dense_col are initialized with no explicit zero
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
    for (int64_t j = k; j < ks; j++)
    {
        // skip if U(ks, Q[j]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row[Q[j]]));
        if (sgn == 0) { continue; }

        // perform j-th IPGE update for U(ks,:)
        SPEX_CHECK(spex_ipge());
        // update last_nz_b4_ks
        last_nz_b4_ks = j;

        // SCALEUP: finish all skipped scaling for column j of L
        // mpq_equal is said to be faster than mpq_cmq
        SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_2D(S, 3, j), one));
        if (r == 0) // S(3,j) != 1
        {
            // S(:,j) = [S(1,j)*S(3,j); S(2,j)*S(3,j); 1]
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, j),
                                    SPEX_2D(S, 1, j), SPEX_2D(S, 3, j)));
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, j),
                                    SPEX_2D(S, 2, j), SPEX_2D(S, 3, j))); 
            SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, j), 1, 1));
        }

        SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_2D(S, 1, j), one));
        if (r == 0) // S(1,j) != 1
        {
            for (pj = 0; pj < L->v[j]->nz; pj++)
            {
                if (pj == Ldiag[j])
                {
                    // L(P(j), j) = sd[j]
                    SPEX_CHECK(SPEX_mpz_set(L->v[j]->x[pj], sd[j]));
                }
                else
                {
                    // L(cj, j) = L(cj,j)*S(1,j)
                    SPEX_CHECK(SPEX_mpz_divexact(L->v[j]->x[pj],
                                            L->v[j]->x[pj],
                                            SPEX_MPQ_DEN(SPEX_2D(S, 1, j))));
                    SPEX_CHECK(SPEX_mpz_mul(L->v[j]->x[pj], L->v[j]->x[pj],
                                            SPEX_MPQ_NUM(SPEX_2D(S, 1, j))));
                }
            }
            SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, j), 1, 1)); 
        }
        // d[j] = U(j,Q(j))
        SPEX_CHECK(SPEX_mpz_set(d[j], U->v[j]->x[Udiag[j]]));

        // insert new entry L(P(ks), j) to L and swap its value with U(ks, Q(j))
        pj = L->v[j]->nz;
        // reallocate the nonzero pattern if needed
        if (pj == L->v[j]->max_nnz)
        {
            SPEX_CHECK(spex_expand(&(L->v[j])));
        }
        SPEX_CHECK(SPEX_mpz_init(L->v[j]->x[pj]));
        // L(P(ks), j) = U(ks, Q(j))
        SPEX_CHECK(SPEX_mpz_swap(L->v[j]->x[pj], Uk_dense_row[Q[j]]));
        L->v[j]->i[pj] = P[ks];
        L->v[j]->nz ++;
    }

    // There must be at least one nonzero in U(ks, Q(k:ks-1)) for this case.
    // Otherwise, this will be handled by spex_dppu1(...). Therefore,
    // last_nz_b4_ks must have updated at least once and be greater than k-1
    SPEX_ASSERT(last_nz_b4_ks > k-1);

    // perform history update up to (last_nz_b4_ks)-th IPGE iteration
    // and remove zero from row ks of U
    pks = 0;
    while (pks < U->v[ks]->nz)
    {
        // column index in row ks of U
        cks = U->v[ks]->i[pks];
        h[cks] = SPEX_FLIP(h[cks]);

        if (h[cks] < last_nz_b4_ks) // require history update
        {
            // U(ks, Q(k:ks-1)) were set to 0 but not removed from nnz pattern
            if (Q_inv[cks] < ks)
            {
                // move the column index of last nonzero to current position
                cks = U->v[ks]->i[U->v[ks]->nz];
                U->v[ks]->i[pks] = cks;
                U->v[ks]->nz--;
                continue;
            }

            // U(ks,cks) = (U(ks,cks)*sd(last_nz_b4_ks))/sd(h[cks]);
            SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row[cks],
                                    Uk_dense_row[cks], sd[last_nz_b4_ks]));
            if (h[cks] > 0)
            {
                SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row[cks],
                                             Uk_dense_row[cks], sd[h[cks]]));
            }
        }
        pks++;
    }

    // update S(3,ks)= S(3,ks)*sd(ks-1)/sd(last_nz_b4_ks)
    if (last_nz_b4_ks != ks-1)
    {
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1])); 
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[last_nz_b4_ks])); 
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 3, ks),
                                SPEX_2D(S, 3, ks), pending_scale));
    }

    // L(P(ks),ks) = U(ks,Q(ks));
    SPEX_CHECK(SPEX_mpz_set(Lk_dense_col[P[ks]], Uk_dense_row[Q[ks]]));
    // d(ks)       = U(ks,Q(ks));
    SPEX_CHECK(SPEX_mpz_set(       d[ks]       , Uk_dense_row[Q[ks]]));
    // sd(ks)      = U(ks,Q(ks))*S(3,ks)
    SPEX_CHECK(SPEX_mpz_divexact(sd[ks], Uk_dense_row[Q[ks]], 
                            SPEX_MPQ_DEN(SPEX_2D(S, 3, ks))));
    SPEX_CHECK(SPEX_mpz_mul(sd[ks],  sd[ks],
                            SPEX_MPQ_NUM(SPEX_2D(S, 3, ks))));

    SPEX_FREE_WORK;
    return SPEX_OK;
}
