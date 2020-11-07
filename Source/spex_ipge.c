//------------------------------------------------------------------------------
//SPEX_CHOLMOD/spex_ipge.c: perform one iteration of IPGE with effort to skip
//                          possible scaling process.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform one iteration of IPGE when the
// involving vectors have pending scale factors that are wished not to be
// applied. In addition, this function should be used in successive IPGE
// update, where history update could be involved for certain entry. Therefore,
// a history vector is required as a input/output. In case of a single IPGE
// iteration with no need for history update, this function should not be used.
// This function is called by the following three functions:
// spex_dppu2: successive IPGE update for row k of U after swapping with row
//             ks of U.
// spex_cppu: compute the (n-1)-th IPGE update of column k after vk is inserted.
// spex_triangular_solve: the REF triangular solve for LDx=v when L and v are
//             sparse.
// spex_forward_sub: the forward substitution when solving LDUx=b, which is
//             essensially REF triangular solve for LDx=b when b is dense.
//
// Algorithm explanation:
// Consider when performing the j-th IPGE update for x[i], and the j-th pivot
// is sd[j] and the j-th vector (row/column) is v. When all entries in v have
// no pending scale factor, v(perm[j])=sd[j], and we can perform IPGE for x[i]
// as
// x[i] = (x[i]*sd[j]-v[i]*x[perm[j]])/sd[j-1].
// In order to skip scaling for v when its pending scaling factor v_scale != 1,
// we can compute x[i] as
// x[i] = (x[i]*v[perm[j]]-v[i]*x[perm[j]])/sd[j-1],
// since v(perm[j])*v_scale=sd[j]. When the update finished, the pending scale
// for x(perm(j+1:n)) needs to multiply with v_scale.
//
// In addition, in case of history update is needed before the IPGE update for
// x[i] and/or x[perm[j]], the equation becomes
// x[i] = x[i]*v[perm[j]]/sd[h[i]]- v(i)*x[perm[j]])/sd[h[perm[j]]].
//
// When the IPGE update finished, all entries x[perm[1:j]] will be scaled such
// that the pending scaling factors become 1, while x[perm[j+1:n]] has common
// factor x_scale=x_scale*v_scale.


#define SPEX_FREE_WORK               \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPQ_CLEAR(tmp_mpz);

#include "spex_internal.h"

SPEX_info spex_ipge
(
    mpz_t *x,
    SPEX_vector *x_sparse,
    SPEX_vector *v,
    int64_t *h,
    int64_t *perm,
    mpz_t *sd,
    mpq_t x_scale,
    mpq_t v_scale,
    int64_t diag_j,
    int64_t j
)
{
    SPEX_info info;
    int64_t p, i, real_hj, real_hi, x_nz;
    int sgn;
    if (x_sparse != NULL)
    {
        x_nz = x_sparse->nz;
    }
    else
    {
        return SPEX_INCORRECT_INPUT;
    }

    // pending_scale = x[perm[j]]/sd[h[perm[j]]]
    real_hj = SPEX_FLIP(h[perm[j]]);
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, x[perm[j]]));
    if (real_hj > -1)
    {
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[real_hj]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
    }
    // NOTE: this could cause fillin in x 
    for (p = 0; p < v->nz; p++)
    {
        // exclude v(perm[j])
        if (p == diag_j) // same as (i == perm[j])
        {
            continue;
        }
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
        if (sgn == 0)    // v[i] == 0
        {
            continue;
        }
        // column/row index in v
        i = v->i[p];
        real_hi = SPEX_FLIP(h[i]);

        // x[i] = floor(x[i]*v[perm[j]]/sd[h[i]])
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[i]));
        if (sgn != 0)    // x[i] != 0
        {
            // x[i] = x[i]*v[perm[j]]
            SPEX_CHECK(SPEX_mpz_mul(x[i], x[i], v->x[diag_j]));
            if (real_hi != real_hj && real_hi > -1)
            {
                SPEX_CHECK(SPEX_mpz_fdiv_q(x[i], x[i], sd[real_hi]));
            }
        }
        else if (h[i] >= -1) // this entry was not in nnz pattern
        {
            // TODO remove this and think about how to handle it efficiently
            // in spex_dppu2
            // TODO: check if this will 0?
            // check if this will be the first off diagonal entry in U(ks,Q)
            if (Q_inv[i] < jnext)
            {
                // jnext is the column index of the found first off-diagonal
                // entry in U(ks,Q)
                jnext = Q_inv[i];
            }

            // reallocate the nonzero pattern if needed
            if (x_nz == x_sparse->max_nnz)
            {
                SPEX_CHECK(spex_expand(&x_sparse));
            }
            // insert new entry in the nonzero pattern
            x_sparse->i[x_nz] = i;
            x_nz++;
        }
        if (real_hi != real_hj)
        {
            // tmp_mpz = floor(v(i)*pending_scale)
            SPEX_CHECK(SPEX_mpz_mul(tmp_mpz, v->x[p],
                                    SPEX_MPQ_NUM(pending_scale)));
            // TODO change or fdiv and cdiv to divexact?
            SPEX_CHECK(SPEX_mpz_fdiv_q(tmp_mpz, tmp_mpz,
                                    SPEX_MPQ_DEN(pending_scale)));
            // x[i] = x[i]- tmp_mpz
            SPEX_CHECK(SPEX_mpz_sub(x[i], x[i], tmp_mpz));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_addmul(x[i], v->x[p], x[perm[j]]));
            SPEX_CHECK(SPEX_divexact(x[i], x[i], sd[real_hi]));
        }

        // update h[i] and last_nz_b4_ks
        h[i] = SPEX_FLIP(j);
    }

    // update # of nnz in x
    x_sparse->nz = x_nz;

    // x(perm(j)) = x[perm[j]]*x_scale*history_update
    if (j-1 > real_hj) // require history update
    {
        SPEX_CHECK(SPEX_mpz_mul(x[perm[j]], x[perm[j]], sd[j-1]));
        if (real_hj > -1)
        {
            SPEX_CHECK(SPEX_mpz_divexact(x[perm[j]], x[perm[j]], sd[real_hj]));
        }
    }
    SPEX_CHECK(SPEX_mpz_divexact(x[perm[j]], x[perm[j]],SPEX_MPQ_DEN(x_scale)));
    SPEX_CHECK(SPEX_mpz_mul(x[perm[j]], x[perm[j]], SPEX_MPQ_NUM(x_scale)));

    // update scaling for x, since the scaling for v is skipped
    // x_scale=x_scale*v_scale.
    SPEX_CHECK(SPEX_mpq_mul(x_scale, x_scale, v_scale));

    return SPEX_OK;
}
