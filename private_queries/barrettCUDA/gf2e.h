// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of version 2 of the GNU General Public License as
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// There is a copy of the GNU General Public License in the COPYING file
// packaged with this plugin; if you cannot find it, write to the Free
// Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA

#ifndef __GF2E_H__
#define __GF2E_H__

#include <stdint.h>

typedef uint8_t GF28_Element;

typedef uint16_t GF216_Element;

/* Use the same representation as AES */
static const GF28_Element GF28_Generator = 0x1b;

extern const GF28_Element GF28_mult_table[256][256];

extern const GF28_Element GF28_inv_table[256];

extern const GF216_Element GF216_exp_table[131070];

extern const GF216_Element GF216_log_table[65536];

//Multiply two GF(2^E) elements
//
template <typename E>
inline E multiply_GF2E(E x, E y);
// specialization for GF(2^8)
template <>
inline GF28_Element multiply_GF2E<GF28_Element>(GF28_Element x, GF28_Element y) {
    return GF28_mult_table[x][y];
}
// specialization for GF(2^16)
template <>
inline GF216_Element multiply_GF2E<GF216_Element>(GF216_Element x, GF216_Element y) {
    if(x == 0 || y == 0) return 0;

    GF216_Element log_x = GF216_log_table[x];
    GF216_Element log_y = GF216_log_table[y];

    return GF216_exp_table[log_x+log_y];
}

//Compute the inverse of a GF(2^E) element
//
template <typename E>
inline E inverse_GF2E(E x);
// specialization for GF(2^8)
template <>
inline GF28_Element inverse_GF2E<GF28_Element>(GF28_Element x) {
    return GF28_inv_table[x];
}
// specialization for GF(2^16)
template <>
inline GF216_Element inverse_GF2E<GF216_Element>(GF216_Element x) {
    if (x == 0) return 0;
    GF216_Element log_x = GF216_log_table[x];
    return GF216_exp_table[65535-log_x];
}

// Evaluate a given polynomial over GF(2^E) at the given index
template <typename E>
E evalpoly_GF2E(E *coeffs, unsigned short degree,
	E index)
{
    E res = 0;

    int i = degree + 1;

    while (i > 0) {
        --i;
        res = multiply_GF2E<E>(res, index);
        res ^= coeffs[i];
    }

    return res;
}

// return f(alpha), where f is the polynomial of degree numpoints-1 such
// that f(indices[i]) = values[i] for i=0..(numpoints-1)
template <typename E>
E interpolate_GF2E(const E *indices, const E *values, 
        unsigned short numpoints, const E alpha) {
    E res;

    E *lagrange = new E[numpoints];
    unsigned short i,j;

    for (i=0;i<numpoints;++i) {
        E numer = 1, denom = 1;
        for (j=0;j<numpoints;++j) {
            if (j==i) continue;
            E numerdiff = indices[j] ^ alpha;
            E denomdiff = indices[j] ^ indices[i];
            numer = multiply_GF2E<E>(numer, numerdiff);
            denom = multiply_GF2E<E>(denom, denomdiff);
        }
        lagrange[i] = multiply_GF2E<E>(numer, inverse_GF2E<E>(denom));
    }

    res = 0;
    for (i=0;i<numpoints;++i) {
        res ^= multiply_GF2E<E>(lagrange[i], values[i]);
    }

    delete[] lagrange;
    return res;
}

// Multiply the polynomial in1 of degree deg1 (a pointer to
// deg1+1 coefficients) by the polynomial in2 of degree deg2
// (a pointer to deg2+1 coefficients) and put the output
// into result (a pointer to deg1+deg2+1 coefficients).
// Note that you cannot multiply "in place"; that is, the
// memory regions pointed to by result and in1 (and result
// and in2) cannot overlap.  in1 and in2 may overlap, however.
// This function requires that deg1 and deg2 each be >= 0 (so you
// can't multiply by the zero polynomial with no coefficients; you
// can of course multiply by the constant (degree 0) polynomial whose
// only coefficient is 0, however).
template <typename E>
inline void multpoly_GF2E(E *result, const E *in1, unsigned short deg1,
                                     const E *in2, unsigned short deg2)
{
    unsigned short resdeg = deg1 + deg2;

    for (unsigned short d=0; d<=resdeg; ++d) {
        result[d] = 0;

        /* lower = max(0,d-deg2) */
        unsigned short lower = 0;
        if (d > deg2) lower = d-deg2;

        /* upper = min(d,deg1) */
        unsigned short upper = d;
        if (deg1 < d) upper = deg1;

        for (unsigned short i=lower; i<=upper; ++i) {
            result[d] ^= multiply_GF2E<E>(in1[i],in2[d-i]);
        }
    }
}

// Create the polynomial whose roots are the ones in the given array
// The output pointer must point to pre-allocated memory of size
// (numroots+1)*sizeof(E)
template <typename E>
inline void polyfromroots_GF2E(E *result, const E *roots,
	unsigned short numroots)
{
    result[numroots] = 1;
    for (unsigned short r=0; r<numroots; ++r) {
	result[numroots-r-1] = 0;
	for (unsigned short i=0; i<=r; ++i) {
	    result[numroots-r+i-1] ^=
		multiply_GF2E<E>(result[numroots-r+i],roots[r]);
	}
    }
}

#endif
