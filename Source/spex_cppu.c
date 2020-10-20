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
//       x 0 0 0 x
//       . x . . .
//       . . x . .
//       . . . x .
//       . . . . x
// This function will swap columns k and ks in L and U. Noted that the columns
// of U are permuted implicitly via the permutation matrix based on Q.

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);         \
    SPEX_FREE(h);

#include "spex_internal.h"

SPEX_info spex_cppu
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
    int64_t ks             // index of the diagonal to be swapped with, [0,n)
)
{
    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t pk, ck, pks, cks;
    // the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
    int64_t Ucp_k_ks;
    int64_t *h;  // history vector

    mpq_t pending_scale, tmp_mpq;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(tmp_mpq);
    mpz_t Uiks, tmp_mpz; SPEX_MPZ_SET_NULL(Uiks); SPEX_MPZ_SET_NULL(tmp_mpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Uiks));
    SPEX_CHECK(SPEX_mpz_init(tmp_mpz));

    h = (int64_t*) SPEX_malloc(n, sizeof(int64_t));
    if (h == NULL)
    {
        SPEX_FREE_WORK;
        return SPEX_OUT_OF_MEMORY;
    }

    //-------------------------------------------------------------------------
    // Backtracking column ks of L and U, the backtracking result will be moved
    // to column k of L and entries in column ks of L are moved to
    // Lk_dense_col.
    //-------------------------------------------------------------------------
    // when ks = k+1, there is no need to make a copy of Lk_dense_col. Instead,
    // we can direct operate on Lk_dense_col.
    // when ks > k+1, Lk_dense_col will firstly be moved to L->v[k]->x in a
    // compressed-column form. Then Lk_dense_col will make a copy of ks-th
    // column of scaled L and scaled U(k:ks-1,Q(ks)). Backtracking will be
    // performed on Lk_dense_col for each nonzero in U(k:ks-1, Q(ks)).
    //-------------------------------------------------------------------------
    if (ks == k+1)
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
        SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row[Q[ks]],
                                SPEX_MPQ_DEN(SPEX_2D(S, 2, k))));
        SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                SPEX_MPQ_NUM(SPEX_2D(S, 2, k))));
        // update d[k]
        SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row[Q[ks]]));

        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpz_set_z(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
            // remove common factor in mpq_den and mpq_num
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        else
        {
            SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
        }
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_2D(S,1,ks)));
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_2D(S,3,ks)));
        // perform backtracking for each nonzero in col ks of L and store
        // results in Lk_dense_col
        // NOTE: this will cause fillin in the k(th) column of L
        pk = 0;
        for (pks = 0; pks < L->v[ks]->nz; pks++)
        {
            // row index in column ks of L
            cks = L->v[ks]->i[pks];

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // L(cks,k) = L(cks,k)*Uiks/L(P[k],k)+L(cks,ks)*pending_scale
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[cks]));
            if (sgn != 0)
            {
                // L(cks,k) = floor(L(cks,k)*Uiks/L(P[k],k))
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[cks],
                                        Lk_dense_col[cks], Uiks));
                SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col[cks],
                                      Lk_dense_col[cks], L->v[k]->x[Ldiag[k]]));

                // tmp_mpz = ceil(L(cks,ks)*pending_scale)
                SPEX_CHECK(SPEX_mpz_mul(tmp_mpz, L->v[ks]->x[pks],
                           SPEX_MPQ_NUM(pending_scale)));
                SPEX_CHECK(SPEX_mpz_cdiv_q(tmp_mpz, tmp_mpz,
                           SPEX_MPQ_DEN(pending_scale)));
                // L(cks,k) = L(cks,k)+tmp_mpz
                SPEX_CHECK(SPEX_mpz_add(Lk_dense_col[cks], Lk_dense_col[cks],
                                        tmp_mpz));
            }
            else  // faster 
            {
                // L(cks,k) = L(cks,ks)*pending_scale
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[cks],
                           L->v[ks]->x[pks], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[cks], Lk_dense_col[cks],
                           SPEX_MPQ_NUM(pending_scale)));
            }

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // move nonzero entry from Lk_dense_col to L->v[k]->x
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[cks]));
            if (sgn != 0)
            {
                // reallocate the nonzero pattern if needed
                if (pk == L->v[k]->max_nnz)
                {
                    SPEX_CHECK(spex_expand(&(L->v[k])));
                }
                L->v[k]->i[pk] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col[cks]));
                pk++;
            }
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // swap values from L->v[ks]->x and Lk_dense_col
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col[cks], L->v[ks]->x[pks]));
        }
        // update L(P(k),k) and reset Lk_dense_col[P[k]]
        // reallocate the nonzero pattern if needed
        if (pk == L->v[k]->max_nnz)
        {
            SPEX_CHECK(spex_expand(&(L->v[k])));
        }
        L->v[k]->i[pk] = P[k];
        SPEX_CHECK(SPEX_mpz_set(sd[k], Uiks));
        SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Uiks));
        SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col[P[k]], 0));
        pk++;
        // update number of nnz
        L->v[k]->nz = pk;

        // get the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
        Ucp_k_ks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
    }
    else  // ks > k+1
    {
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // construct column k of L based on Lk_dense_col
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        for (pk = 0; pk < L->v[k]->nz; pk++)
        {
            ck = L->v[k]->i[pk];
            if (ck == P[k])
            {
                Ldiag[k] = pk;  // update Ldiag[k]
            }
            // swap the entries in the Lk_dense_col and L->v[k]->x
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col[ck]));
        }

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // initialize history vector h and copy L->v[ks]->x to Lk_dense_col.
        // Explicit zero(s) are kept if exist
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        for (pks = 0; pks < L->v[ks]->nz; pks++) 
        { 
            cks = L->v[ks]->i[pks]; 
            // copy value from L->v[ks]->x to Lk_dense_col
            SPEX_CHECK(SPEX_mpz_set(Lk_dense_col[cks], L->v[ks]->x[pks]));
            // Assuming that zero resulted from exact cancellation is kept in
            // the nnz pattern, then performing backtracking will not introduce
            // new fillin. That is, for any i such that U(i,Q[ks])!=0, the nnz
            // pattern of L(i+1:n, i) is a subset of the nnz pattern of
            // Lk_dense_col(i+1:n). Therefore, the final nnz pattern can be
            // found by L->v[ks]->i and column-wise nnz pattern of U(:,ks)
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
            h[cks] = ks-1;
        }

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // backtrack column ks of L and U for each nnz in U(k:ks-1,Q(ks))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        pks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
        i = Uci[pks];         // row index
        while (i >= k)
        {
            // skip if U(i, Q[ks]) turns out to be explicit zero
            // no need to initialize history vector for this entry
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx[pks]]));
            if (sgn == 0)   { continue; }

            // Uiks = U(i,Q(ks))*S(2,k)*S(3,k)
            if (i == k)
            {
                // update S(2,k)
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, k), SPEX_2D(S, 2, k),
                                        SPEX_2D(S, 3, k)));
                // update d[k]
                SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row[Q[ks]],
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
                SPEX_CHECK(SPEX_mpz_set_z(pending_scale, sd[i-1]));
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

                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[ci]));
                if (sgn != 0)
                {
                    if (i < real_h || r == 0 || pks == Ucp[Q[ks]+1]-2)
                    {
                        // Lk_dense_col = Lk_dense_col*pending_scale/sd(h[ci])
                        //                + L(:,i)*Uiks/L(P(i),i)
                        // use L(P(i),i) instead of sd[i] to avoid scaling
                        // for L(:,i)
                        SPEX_CHECK(SPEX_mpz_mul(tmp_mpz, L->v[i]->x[pi], Uiks));
                        SPEX_CHECK(SPEX_mpz_cdiv_q(tmp_mpz, tmp_mpz,
                                   L->v[i]->x[Ldiag[i]]));

                        SPEX_CHECK(SPEX_mpq_set_ui(tmp_mpq, 1, 1));
                        SPEX_CHECK(SPEX_mpq_set_den(tmp_mpq, sd[h[ci]]));
                        SPEX_CHECK(SPEX_mpq_mul(tmp_mpq,tmp_mpq,pending_scale));

                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[ci],
                                   Lk_dense_col[ci], SPEX_MPQ_NUM(tmp_mpq)));
                        SPEX_CHECK(SPEX_mpz_fdiv(Lk_dense_col[ci],
                                   Lk_dense_col[ci], SPEX_MPQ_DEN(tmp_mpq)));

                        SPEX_CHECK(SPEX_mpz_add(Lk_dense_col[ci],
                                   Lk_dense_col[ci], tmp_mpz));

                        // update h[ci]
                        h[ci] = i-1;
                        continue;
                    }
                    else if (i > 0)
                    {
                        // Lk_dense_col = (Lk_dense_col*sd(i-1)
                        //                                 + L(:,i)*Uiks)/sd(i)
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[ci],
                                                Lk_dense_col[ci], sd[i-1]));
                    }
                }
                SPEX_CHECK(SPEX_mpz_addmul(Lk_dense_col[ci],
                                           L->v[i]->x[pi], Uiks));
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[ci],
                                             Lk_dense_col[pi], sd[i]));

                // update h[ci]
                h[ci] = i-1;
            }

            // move Uiks (scaled U(i,Q(ks)) to Lk_dense_col
            SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col[P[i]], Uiks));
            // update corresponding entry in the history vector
            h[P[i]] = i-1;

            // get next nnz
            pks --;
            if (pks < 0) { break;}
            i = Uci[pks];         // row index
        }

        //
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // perform history update if needed
        // move nonzero entry from Lk_dense_col to L->v[k]->x
        // and swap values from L->v[ks]->x and Lk_dense_col
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        pk = 0;
        for (pks = 0; pks < L->v[ks]->nz; pks++) 
        {
            cks = L->v[ks]->i[pks]; 
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[cks]));
            if (sgn != 0)
            {
                // check if need to perform history update
                if (h[cks] != k-1)
                {
                    if (k > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[cks],
                                   Lk_dense_col[cks], sd[k-1]));
                    }
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[cks],
                               Lk_dense_col[cks], sd[h[cks]]));
                }
                // reallocate the nonzero pattern if needed
                if (pk == L->v[k]->max_nnz)
                {
                    SPEX_CHECK(spex_expand(&(L->v[k])));
                }
                L->v[k]->i[pk] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col[cks]));
                pk++;
            }
            SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col[cks], L->v[ks]->x[pks]));
        }
        for (pks = Ucp[Q[ks]+1]-2; pks >= Ucp[Q[ks]]; pks--)
        {
            i = Uci[pks];         // row index
            ci = P[i];
            if (i == k)
            {
                // reallocate the nonzero pattern if needed
                if (pk == L->v[k]->max_nnz)
                {
                    SPEX_CHECK(spex_expand(&(L->v[k])));
                }
                L->v[k]->i[pk] = ci;
                SPEX_CHECK(SPEX_mpz_set(sd[k], Lk_dense_col[ci]));
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col[ci]));
                pk++;
                Ucp_k_ks = pks;
                break;  // this is the last nonzero
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col[ci]));
                if (sgn != 0)
                {
                    // check if need to perform history update
                    if (h[ci] != k-1)
                    {
                        if (k > 0)
                        {
                            SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col[ci],
                                       Lk_dense_col[ci], sd[k-1]));
                        }
                        SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col[ci],
                                   Lk_dense_col[ci], sd[h[ci]]));
                    }
                    // reallocate the nonzero pattern if needed
                    if (pk == L->v[k]->max_nnz)
                    {
                        SPEX_CHECK(spex_expand(&(L->v[k])));
                    }
                    L->v[k]->i[pk] = ci;
                    SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col[ci]));
                    pk++;
                }
            }
        }
        // update number of nnz
        L->v[k]->nz = pk;
    }
    // update S(1,k) and S(2,k) as 1, since all entry in L(:,k) are scaled
    SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 1, k), one));
    SPEX_CHECK(SPEX_mpq_set(SPEX_2D(S, 3, k), one));

    //-------------------------------------------------------------------------
    // RwSOP
    // NOTE: If explicit zero(s) from exact cancelllation of IPGE update are
    // kept, no fillin will be added.
    //-------------------------------------------------------------------------
    // pending_scale = U(k, Q(ks))/ U(k, Q(k))
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, Uk_dense_row[Q[ks]]));
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row[Q[k]]));
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    // iterate for all nnz in U(k+1:ks,Q(ks))
    for (pks = Ucp_k_ks+1; pks < Ucp[Q[ks]+1]; pks++)
    {
        i = Uci[pks];         // row index

        // skip scaling for frames between iterations
        // REMARK: Uci[pks-1] must be at least equal to k
        for (i1 = Uci[pks-1]+1; i1 < i; i1++)
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

        for (pi = 0; pi < U->v[i]->nz; pi++)
        {
            ci = U->v[i]->i[pi];
            // handle U(i, Q[ks]) after iteration
            if (ci == Q[ks]) { continue;}
            else if (Q_inv[ci] < ks)
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
                // perform RwSOP to U(i,Q(ks+1:n+1)), NO FILLIN
                // U(i,ci)= (U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci))/U(k,Q(k))
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_submul(U->v[i]->x[pi], U->v[i]->x[Ucx[pks]],
                                        Uk_dense_row[ci]));
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row[Q[k]]));
            }

        }
        pi = Ucx[pks];  // pointer for U(i,Q(ks))
        // flip sign for U(i, Q(ks))
        SPEX_CHECK(SPEX_mpz_neg(U->v[i]->x[pi], U->v[i]->x[pi]));
        // update the column index of U(i, Q(ks)) to Q(k)
        U->v[i]->i[pi] = Q[k];

        // skip scaling for L(:,i) by setting S(1,i) = S(1,i)*pending_scale
        if (i != ks)
        {
            SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, i),
                                    SPEX_2D(S, 1, i), pending_scale));
        }
        else
        {
            SPEX_CHECK(SPEX_mpq_neg(SPEX_2D(S, 1, i), SPEX_2D(S, 1, i)));
            SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
        }
    }

    //-------------------------------------------------------------------------
    // update column permutation
    //-------------------------------------------------------------------------
    int64_t tmp_int;
    tmp_int = Q[k];    Q[k] = Q[ks];          Q[ks] = tmp_int;

    //-------------------------------------------------------------------------
    // flip sign for columns and rows ks+1 to n and update corresponding sd
    //-------------------------------------------------------------------------
    for (i = ks+1; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpq_neg(SPEX_2D(S, 3, i), SPEX_2D(S, 3, i)));
        SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
    }

    SPEX_FREE_WORK;
    return SPEX_OK;
}
