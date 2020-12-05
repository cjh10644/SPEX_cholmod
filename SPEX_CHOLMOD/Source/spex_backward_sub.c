//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_backward_sub: sparse REF backward substitution (x = U\b)
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF backward substitution, solving
 * the system Ux = b. x is internally multiplied by the determinant of A to
 * maintain integral. x_scale is set to 1/det(A) in the output. Once x_scale is
 * applied to x, the solution vector will not be in integer domain.
 *
 * U is a sparse mpz matrix, which is stored by row, and x is a dense mpz
 * vector.
 *
 * The input argument x contains b on input, and it is overwritten on output
 * by the solution x.
 */

#define SPEX_FREE_ALL                \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_internal.h"

SPEX_info spex_backward_sub  // performs sparse REF backward substitution
(
    mpz_t *x,               // right hand side vector
    mpq_t x_scale,          // pending scale for x, applying this x_scale to
                            // x will result in a non-integer value
    const SPEX_matrix *U,   // input upper triangular matrix
    const mpq_t *S,         // the pending scale factor matrix
    const mpz_t *sd,        // array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv    // inverse of column permutation
)
{
    SPEX_info info ;
    int sgn;
    mpq_t pending_scale; SPEX_MPQ_SET_NULL(pending_scale);// TODO make input
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    // Start at x[n-1], since x[n] will remain the same
    for (int64_t i = U->n-2; i >= 0; i--)
    {
        // tmpz = 0
        SPEX_CHECK(SPEX_mpz_set_ui(tmpz, 0));

        for (int64_t p = 0; p < U->v[i]->nz; p++)
        {
            int64_t j = U->v[i]->i[p];// the real col index is Q_inv[j]

            // skip if diagonal or corresponding entry in x is zero
            if (Q_inv[j] == i) {continue;}
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, x[P[Q_inv[j]]]));
            if (sgn == 0) {continue;}

            // tmpz -= U(i,j)*x[j]
            SPEX_CHECK(SPEX_mpz_submul(tmpz, U->v[i]->x[p], x[P[Q_inv[j]]]));
        }
        // pending_scale = S(2, i)*S(3, i)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                SPEX_2D(S, 2, i), SPEX_2D(S, 3, i)));
        // tmpz = tmpz*pending_scale
        SPEX_CHECK(SPEX_mpz_divexact(tmpz, tmpz, SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(pending_scale)));

        // x[i] = x[i]*sd[n-1]+tmpz
        SPEX_CHECK(SPEX_mpz_mul(x[i], x[i], sd[n-1]));
        SPEX_CHECK(SPEX_mpz_add(x[i], x[i], tmpz));

        // x[i] = x[i]/sd[i]
        SPEX_CHECK(SPEX_mpz_divexact(x[i], x[i], sd[i]));
    }

    // set x_scale = 1/sd[n]
    SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));
    SPEX_CHECK(SPEX_mpz_set(SPEX_MPQ_DEN(x_scale), sd[n-1]));

    SPEX_FREE_ALL;
    return (SPEX_OK) ;
}