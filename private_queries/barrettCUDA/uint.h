// This file is part of BarrettCUDA v0.1.
// 
// BarrettCUDA is a fast(ish) implementation of finite field sparse
// matrix-vector multiplication (SpMV) for Nvidia GPU devices, written
// in CUDA C++. BarrettCUDA supports SpMV for matrices expressed in
// the 'compressed column storage' (CCS) sparse matrix representation
// over (i) the field of integers modulo an arbitrary multi-precision
// prime, or (ii) either of the binary fields GF(2^8) or GF(2^16).
// 
// Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz.
// 
// BarrettCUDA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// BarrettCUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with BarrettCUDA. If not, see <http://www.gnu.org/licenses/>.

#ifndef __UINT_H__
#define __UINT_H__

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

#define BITS_PER_LIMB	(sizeof(uint) * 8)
#define BITS_IN(limbs)	((limbs) * BITS_PER_LIMB)

typedef unsigned int uint;
typedef uint  uint32;
typedef uint2 uint64;
typedef uint3 uint96;
typedef uint4 uint128;

template <typename T>
struct uintXp
{
    T lo;
    uint hi;
};

template <typename T>
static inline NTL::ZZ to_ZZ(const T & n)
{
    //std::cout<<"to_ZZ inside uint.h"<<"\n";
    //std::cout<<sizeof(T)<<"\n";
    return NTL::ZZFromBytes((unsigned char *)&n, sizeof(T));
}

template <typename T>
static inline T & to_uint(const NTL::ZZ & n, T & ret)
{
    NTL::BytesFromZZ((unsigned char *)&ret, n, sizeof(T));
    return ret;
}

template <typename T>
static inline T & to_uint(const NTL::ZZ_p & n, T & ret)
{
    NTL::BytesFromZZ((unsigned char *)&ret, rep(n), sizeof(T));
    return ret;
}

template <typename T>
static inline NTL::ZZ_p to_ZZ_p(const T & n)
{
    return NTL::to_ZZ_p(to_ZZ(n));
}

template <typename T>
static inline void print(const T & n)
{
    std::cout << to_ZZ<T>(n);
}

template <typename T>
static inline void print_limbs(const T & n, const uint limbs)
{
    if (!limbs) return;
    const uint * _n = (uint *)&n;
    std::cout << _n[0];
    for (int i = 1; i < limbs; ++i) std::cout << "." << _n[i];
	std::cout<<"\n";
}

template <typename T>
static inline void print_limbs(const NTL::ZZ & n, const uint limbs)
{
    T ret;
    print_limbs(to_uint(n, ret), limbs);
}

template <typename T>
static inline void print_limbs(const NTL::ZZ_p & n, const uint limbs)
{
    T ret;
    print_limbs(to_uint(n, ret), limbs);
}

template <typename T>
__device__ static inline void _print_limbs(const T & n, const uint limbs)
{
    if (!limbs) return;
    const uint * _n = (uint *)&n;
    printf("%u", _n[0]);
    for (int i = 1; i < limbs; ++i) printf(".%u", _n[i]);
	printf("\n");
}

#endif
