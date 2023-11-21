## This file is part of BarrettCUDA v0.1.
## 
## BarrettCUDA is a fast(ish) implementation of finite field sparse
## matrix-vector multiplication (SpMV) for Nvidia GPU devices, written
## in CUDA C++. BarrettCUDA supports SpMV for matrices expressed in
## the 'compressed column storage' (CCS) sparse matrix representation
## over (i) the field of integers modulo an arbitrary multi-precision
## prime, or (ii) either of the binary fields GF(2^8) or GF(2^16).
## 
## Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz.
## 
## BarrettCUDA is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License,
## or (at your option) any later version.
## 
## BarrettCUDA is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with BarrettCUDA. If not, see <http://www.gnu.org/licenses/>.

LIMBS=4

NVCC=nvcc
CXXFLAGS=-Wall -g -O2 -pedantic -std=c++11
NVCCFLAGS=-std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -Xcompiler -fopenmp -lpthread
GCCFLAGS=-std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -fopenmp -lpthread
LDLIBS=-lntl -lgmp -L$(HOME)/sw/lib -I$(HOME)/sw/include

barrett: VspM_barrett_GPU.cu barrett.h uint.h uintX.h gf2e.h gf2e.cc
	$(NVCC) $(NVCCFLAGS) -o $@ VspM_barrett_GPU.cu gf2e.cc $(LDLIBS)
	
genuint: genuint.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

run-genuint: genuint
	./genuint $(LIMBS) > uintX.h

.PHONY: clean run-genuint

clean:
	-rm -f barrett genuint
