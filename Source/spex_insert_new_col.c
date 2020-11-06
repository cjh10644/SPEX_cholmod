//------------------------------------------------------------------------------
//SPEX_CHOLMOD/spex_insert_new_col.c: finish history update if needed and
//insert the given vector to column k of L and U.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to finish any remaining history update for the
// given vector x and insert x to column k of L and U, while maintaining the
// scale factors for L and U up to date after the insertion.

#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);

#include "spex_internal.h"

SPEX_info spex_insert_new_col
(
    mpz_t *x,
    mpq_t Sx,
    SPEX_matrix *L,
    SPEX_vector *v,
    int64_t k,
    int64_t *h,
    int64_t *P,
    mpz_t *d,
    mpz_t *sd
)
{
    p = 0;
    while (p < v->nz)
    {
        i = v->i[p];
        if (P[i] < =k) // to be appeared in U which won't need history update
        {
            // Since a nonzero UNSCALED entry will be inserted to the row i of
            // U, we need to SCALEUP: scale up row i of U if the scaling factor
            // is not 1.
            // d[i] = L(P[i],i) = d[i]*S(2,i)
            SPEX_CHECK(SPEX_mpz_set(d[i], L->v[i]->x[Ldiag[i]]));
            // mpq_equal is said to be faster than mpq_cmq
            SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_2D(S, 3, i), one));
            if (r == 0) // S(3,i) != 1
            {
#if 0
                // find the gcd of inserted entry and numerator of scaling
                // factor for row i of U and column i of L
                SPEX_CHECK(SPEX_mpz_gcd(gcd, x[P[i]],
                                        SPEX_MPQ_NUM(SPEX_2D(S, 3, i))));
                SPEX_CHECK(SPEX_mpz_divexact(x[P[i]], x[P[i]], gcd));
                SPEX_CHECK(SPEX_mpz_divexact(SPEX_MPQ_NUM(SPEX_2D(S, 3, i)),
                                             SPEX_MPQ_NUM(SPEX_2D(S, 3, i)),
                                             gcd));
#endif
                // S(:,i) = [S(1,i)*S(3,i); S(2,i)*S(3,i); 1]
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 1, i),
                                        SPEX_2D(S, 1, i), SPEX_2D(S, 3, i)));
                SPEX_CHECK(SPEX_mpq_mul(SPEX_2D(S, 2, i),
                                        SPEX_2D(S, 2, i), SPEX_2D(S, 3, i)));
#if 0
                SPEX_CHECK(SPEX_mpq_set_z(SPEX_2D(S, 3, i), gcd));
#else
                SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));
#endif
            }

            SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_2D(S, 2, i), one));
            if (r == 0) //S(2, i) != 1
            {
#if 0
                // find the gcd of inserted entry and numerator of scaling
                // factor for row i of U
                SPEX_CHECK(SPEX_mpz_gcd(gcd, x[P[i]],
                                        SPEX_MPQ_NUM(SPEX_2D(S, 2, i))));
                SPEX_CHECK(SPEX_mpz_divexact(x[P[i]], x[P[i]], gcd));
                SPEX_CHECK(SPEX_mpz_divexact(SPEX_MPQ_NUM(SPEX_2D(S, 2, i)),
                                             SPEX_MPQ_NUM(SPEX_2D(S, 2, i)),
                                             gcd));
#endif
                for (p = 0; p < U->v[i]->nz; p++)
                {
                    // Since U(i, cp) will be integer after scale, we can
                // perform division first to make it small, and this
                // division will preserve integer propety
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p], U->v[i]->x[p],
                                           SPEX_MPQ_DEN(SPEX_2D(S, 2, i))));
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p], U->v[i]->x[p],
                                           SPEX_MPQ_NUM(SPEX_2D(S, 2, i))));
                }
#if 0
                SPEX_CHECK(SPEX_mpq_set_z(SPEX_2D(S, 2, i), gcd));
#else
                SPEX_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
#endif
            }
            // d[i] = L(P(i),i)
            SPEX_CHECK(SPEX_mpz_set(d[i], L->v[i]->x[Ldiag[i]]));

            // append x[P[i]] to row i of U
            pi = U->v[i]->nz;
            // reallocate the nonzero pattern if needed
            if (pi == U->v[i]->max_nnz)
            {
                SPEX_CHECK(spex_expand(&(U->v[i])));
            }
            U->v[i]->nz ++;
            U->v[i]->i[pi] = k;
            SPEX_CHECK(SPEX_mpz_init(U->v[i]->x[pi]));
            if (P[i] == k) // diagonal entry
            {
                // TODO history update?
                SPEX_CHECK(SPEX_mpz_set(U->v[i]->x[pi], x[P[i]]));
                SPEX_CHECK(SPEX_mpz_swap(v->x[p], x[i]));
                p++;
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[pi], x[P[i]]));

                // move the column index of last nonzero to current position
                v->i[p] = v->i[v->nz];
                v->nz --;
            }
        }
        else // to be appeared in L, which will need history update
        {
            h[i] = SPEX_FLIP(h[i]);
            if (h[i] < k)
            {
                // x[i] = x[i]*sd[k]/sd[h[i]]
                SPEX_CHECK(SPEX_mpz_mul(x[i], x[i], sd[k]));
                if (h[i] > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(x[i], x[i], sd[h[i]]));
                }
            }
            SPEX_CHECK(SPEX_mpz_swap(v->x[p], x[i]));
            p++;
        }
    }
    SPEX_vector *tmp_v;
    tmp_v = L->v[n]; /*TODO*/ L->v[n] = v; v = tmp_v;

    SPEX_FREE_WORK;
    return SPEX_OK;
}
