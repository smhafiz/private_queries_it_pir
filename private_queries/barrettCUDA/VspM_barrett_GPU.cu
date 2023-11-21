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

#include <atomic>
#include <chrono>
#include <ratio>
#include <fstream>

#include "barrett.h"
#include "uintX.h"

#define DEBUG true

#define THREADS_PER_BLOCK(n) (n >= 512 ? 512 : n)
//#define THREADS_PER_BLOCK(n) (n >= 32 ? 32 : n)
#define NUM_BLOCKS(n) ((n + THREADS_PER_BLOCK(n) - 1) / THREADS_PER_BLOCK(n))

NTL_CLIENT

// specialization for uintX
template <typename T> struct _SpMV_specializer<T,0>
{
    static __device__ void device_SpMV(T * response, const T * query,
	const uint nvals, const T * vals, const uint ncols, const uint * cols,
	const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	uintXp<T> res_lo = { 0 };
	T res_hi = { 0 };
	uint overflow = { 0 };

	// do the SpMV
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    mad(res_lo.lo, res_hi, overflow, vals[j], query[rows[j]]);
	}

	T * subtrahends = (T *)d_subtrahends;
	uintXp<T> * mu = (uintXp<T> *)d_mu;
	T * modulus = (T *)d_modulus;

	// do the Barrett reduction
	normalize(res_lo.lo, res_hi, subtrahends[overflow], (overflow ? -1: 0));
	uintXp<T> q = get_q(res_lo.lo, res_hi, *mu);
	uintXp<T> r2 = get_r2(q, *modulus);
	res_lo.hi = sub(res_lo.lo, res_hi, r2);
	if (res_lo.hi) sub_modulus(res_lo, *modulus);
	if (res_lo.hi) sub_modulus(res_lo, *modulus);

	// write final result to global memory
	response[i] = res_lo.lo;
    }
};

// specialization for GF28_Element
template <> struct _SpMV_specializer<GF28_Element,0>
{
    static __device__ void device_SpMV(GF28_Element * response,
	const GF28_Element * query, const uint nvals, const GF28_Element * vals,
	const uint ncols, const uint * cols, const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	GF28_Element res = 0;
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    res ^= d_GF28_mult_table[vals[j]][query[rows[j]]];
	}
	response[i] = res;
    }
};

// specialization for GF216_Element
template <> struct _SpMV_specializer<GF216_Element,0>
{
    static __device__ void device_SpMV(GF216_Element * response,
	const GF216_Element * query, const uint nvals,
	const GF216_Element * vals, const uint ncols, const uint * cols,
	const uint * rows)
    {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= ncols) return;

	GF216_Element res = 0;
	for (int j = cols[i]; j < cols[i+1]; ++j)
	{
	    GF216_Element log_x = d_GF216_log_table[vals[j]];
	    GF216_Element log_y = d_GF216_log_table[query[rows[j]]];
	    res ^= d_GF216_exp_table[log_x+log_y];
	}
	response[i] = res;
    }
};

template <typename T>
void SpMV_ntl(NTL::vec_ZZ_p & response, const T * query, 
    const SparseMatrix<T> & matrix)
{
    for (int i = 0; i < matrix.ncols; i++)
    {
	response[i] = NTL::to_ZZ_p(0);
	for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; ++j)
	{
	    response[i] += to_ZZ_p(matrix.l_vals[j])
			 * to_ZZ_p(query[matrix.l_rows[j]]);
	}
    }
    //created by Warren
    //std::cout<<"SpMV_ntl output: \n";
    //std::cout<<response<<"\n";
}

#ifdef DEBUG
    template <typename T>
    void SpMV_ntl_barrett(NTL::vec_ZZ_p & response, const T * query,
	const SparseMatrix<T> & matrix, struct BarrettParams<T> & barrett)
    {
	NTL::vec_ZZ response_ZZ(INIT_SIZE, matrix.ncols);
	for (int i = 0; i < matrix.ncols; i++)
	{
	    response_ZZ[i] = NTL::to_ZZ(0);

	    for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; ++j)
	    {
		response_ZZ[i] += to_ZZ(matrix.l_vals[j]) * to_ZZ(query[matrix.l_rows[j]]);
	    }
	    uint overflow = (uint)NTL::trunc_long(response_ZZ[i] >> 2*BITS_IN(LIMBS_PER_UINTX), BITS_IN(sizeof(uint)));
	    response_ZZ[i] -= barrett.l_subtrahends[overflow];
	    NTL::ZZ q1 = response_ZZ[i] >> BITS_IN(LIMBS_PER_UINTX-1);
	    NTL::ZZ q2 = q1 * barrett.l_mu;
	    NTL::ZZ q3 = q2 >> BITS_IN(LIMBS_PER_UINTX+1);
	    NTL::ZZ r1 = response_ZZ[i] % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    NTL::ZZ r2 = q3 * barrett.l_modulus % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    NTL::ZZ r = (r1 - r2) % NTL::power2_ZZ(BITS_IN(LIMBS_PER_UINTX+1));
	    response[i] = NTL::to_ZZ_p(r);
	}
    //created by Warren
    //std::cout<<"query in SpMV_ntl_barrett: \n";
    //for(int i = 0; i <4; i++){std::cout<<to_ZZ(query[i])<<"\n";}
    //std::cout<< "Response in SpMV_ntl_barrett: \n";
    //std::cout<< response <<"\n";
    }
#endif // DEBUG

template <typename T>
void SpMV(T * l_response, const T * l_query,
	T * d_response, T * d_query, const cudaStream_t & stream,
	const SparseMatrix<T> & matrix)
{
    gpuErrchk(cudaMemcpyAsync(d_query, l_query, matrix.nrows * sizeof(T),
	cudaMemcpyHostToDevice, stream));

    const dim3 Dg(NUM_BLOCKS(matrix.ncols), 1, 1);
    const dim3 Db(THREADS_PER_BLOCK(matrix.ncols), 1, 1);
    const size_t Ns = 0;

    SpMV_kernel<T> <<< Dg, Db, Ns, stream >>> (d_response, d_query,
	matrix.nvals, matrix.d_vals, matrix.ncols, matrix.d_cols, matrix.d_rows);
    gpuErrchk(cudaPeekAtLastError());

    gpuErrchk(cudaMemcpyAsync(l_response, d_response, matrix.ncols * sizeof(T),
	cudaMemcpyDeviceToHost, stream));
}

int main(int argc, char ** argv)
{
    int nstreams = 4;

    if (argc < 4)
    {
	std::cout << "Usage: " << argv[0] << " VALUES ROWS COLS\n\n";
	return 1;
    }

    time_t t0 = time(0);
    t0 = 1481878728;
    //NTL::SetSeed(to_ZZ(t0));
    //std::cout << "seed: " << t0 << "\n";

    struct SparseMatrix<uintX> matrix = { 0 };
    NTL::ZZ modulus;

    uint max_overflow;
    initMatrix(argv[1], argv[2], argv[3], modulus, matrix, max_overflow);
    NTL::ZZ_p::init(modulus);
    //std::cout<<"Line 201"<<std::endl;
    struct BarrettParams<uintX> barrett;
    initBarrett<uintX>(modulus, barrett, max_overflow);

    uintX * l_query, * d_query;
    gpuErrchk(cudaMallocHost((void**)&l_query,
	nstreams * matrix.nrows * sizeof(uintX)));
    gpuErrchk(cudaMalloc((void**)&d_query,
	nstreams * matrix.nrows * sizeof(uintX)));

    uintX * l_response, * d_response;
    gpuErrchk(cudaMallocHost((void**)&l_response,
	nstreams * matrix.ncols * sizeof(uintX)));
    gpuErrchk(cudaMalloc((void**)&d_response,
	nstreams * matrix.ncols * sizeof(uintX)));

    cudaStream_t * streams = new cudaStream_t[nstreams];
    for (int i = 0; i < nstreams; ++i) cudaStreamCreate(&streams[i]);

    NTL::vec_vec_ZZ_p responses(INIT_SIZE, nstreams,
	NTL::vec_ZZ_p(INIT_SIZE, matrix.ncols));

    for (int i = 0; i < nstreams * matrix.nrows; i++)
    {
    NTL::ZZ_p temp = NTL::random_ZZ_p();
    //NTL::ZZ_p temp = to_ZZ_p(0); // inserted by Warren
    //std::cout<<temp<<"\n";
	to_uint<uintX>(temp, l_query[i]);
    }
    //to_uint<uintX>(to_ZZ_p(1), l_query[3]); // by Warren
     
    /*
    //input by Warren
    std::cout<<modulus<<"\n";
    std::cout<<"nrows: " << matrix.nrows<<"\n";
    std::cout<<"uintX size: " << sizeof(uintX)<<"\n";
    std::cout<<"l_query size: " << sizeof(l_query)<<"\n";
    std::cout<<"LIMBS_PER_UINTX: "<<LIMBS_PER_UINTX<<"\n";
    //std::cout << "Matrix: " << matrix. << "\n";
    std::cout << "l_query w/o p: \n";
    
    for(int i = 0;i< matrix.nrows;i++){
        //print_limbs<uintX>(l_query[i], LIMBS_PER_UINTX);
        std::cout<<to_ZZ(l_query[i])<<"\n";
    } 
    std::cout << "l_query in p: \n";
    
    for(int i = 0;i< matrix.nrows;i++){
        //print_limbs<uintX>(l_query[i], LIMBS_PER_UINTX);
        std::cout<<to_ZZ_p(l_query[i])<<"\n";
    }
    
    //std::cout<<to_string(l_query).name()<<"\n";
    //std::cout<<strcat(l_query)<<"\n";
    std::cout<<"\n";
    */
    //end Warren print

    std::atomic<int> cnt = ATOMIC_VAR_INIT(0);
    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds onesec{1000000000};
    //std::cout<<"Line 231"<<std::endl;
    while (std::chrono::duration_cast<std::chrono::duration<int,std::nano>>(std::chrono::high_resolution_clock::now() - start) < onesec)
    {
    // by Warren: commented line below to solve CPU problem uncomment to use GPU mode
	//#pragma omp parallel
	for (int i = 0; i < nstreams; i++)
	{
	    uintX * __l_response = l_response + i * matrix.ncols;
	    uintX * __d_response = d_response + i * matrix.ncols;
	    uintX * __l_query = l_query + i * matrix.nrows;
	    uintX * __d_query = d_query + i * matrix.nrows;		

	    SpMV<uintX>(__l_response, __l_query, __d_response,
		__d_query, streams[i], matrix);
 
	    //SpMV_ntl(responses[i], __l_query, matrix);
        
	    //SpMV_ntl_barrett(responses[i], __l_query, matrix, barrett);
       
        
	    std::atomic_fetch_add(&cnt, 1);
        
//	    for (int j = 0; j < matrix.nrows; j++)
//	    {
//		to_uint<uintX>(NTL::rep(NTL::random_ZZ_p()), __l_query[j]);
//	    }
	}
    }
    // Output
    //std::cout <<"Completed " << cnt << " SpMV per second\n";
    std::cout << cnt<<"\n"; // output for bash file data

    /*
    //Input by Warren
    std::cout << "l_response w/0 p: \n";
    for(int i = 0;i< 8;i++){
        //print_limbs<uintX>(l_response[i], LIMBS_PER_UINTX);
        std::cout<<to_ZZ_p(l_response[i])<<"\n";
    }
    */
    // end Warren Print

    // cleanup
    for (int i = 0; i < nstreams; ++i) gpuErrchk(cudaStreamDestroy(streams[i]));
    delete [] streams;
    responses.kill();

    gpuErrchk(cudaFreeHost(l_query));
    gpuErrchk(cudaFree(d_query));
    gpuErrchk(cudaFreeHost(l_response));
    gpuErrchk(cudaFree(d_response));

    freeBarrett<uintX>(barrett);
    freeMatrix<uintX>(matrix);

    return 0;
}


template <typename T>
void initMatrix(const char * valfile, const char * rowfile,
	const char * colfile, NTL::ZZ & modulus,
	struct SparseMatrix<T> & matrix, uint & max_overflow)
{
    std::ifstream valstream(valfile, std::ifstream::in);
    if (!valstream) { cerr << "Error: opening VALS files\n"; exit(-1); }
    std::ifstream rowstream(rowfile, std::ifstream::in);
    if (!rowstream) { cerr << "Error: opening ROWS files\n"; exit(-1); }
    std::ifstream colstream(colfile, std::ifstream::in);
    if (!colstream) { cerr << "Error: opening COLS files\n"; exit(-1); }
    NTL::ZZ tmp_zz;

    valstream >> tmp_zz;

    modulus = NTL::trunc_ZZ(tmp_zz,sizeof(T)*8);

    rowstream >> matrix.nrows;
    rowstream >> matrix.nvals;
    matrix.l_rows = (uint *)malloc(matrix.nvals * sizeof(uint));
    gpuErrchk(cudaMalloc((void**)&matrix.d_rows, matrix.nvals * sizeof(uint)));
    matrix.l_vals = (T *)malloc(matrix.nvals * sizeof(T));
    gpuErrchk(cudaMalloc((void**)&matrix.d_vals, matrix.nvals * sizeof(T)));

    colstream >> matrix.ncols;
    matrix.l_cols = (uint *)malloc((matrix.ncols+1) * sizeof(uint));
    gpuErrchk(cudaMalloc((void**)&matrix.d_cols,
	(matrix.ncols+1) * sizeof(uint)));

    NTL::ZZ_pPush p(modulus);
    for (int i = 0; i < matrix.nvals; i++)
    {
	NTL::ZZ_p tmp;
	valstream >> tmp;
	to_uint<T>(NTL::rep(tmp), matrix.l_vals[i]);
	rowstream >> matrix.l_rows[i];
    }
    valstream.close();
    rowstream.close();
    gpuErrchk(cudaMemcpy(matrix.d_vals, matrix.l_vals, matrix.nvals * sizeof(T),
	cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(matrix.d_rows, matrix.l_rows,
	matrix.nvals * sizeof(uint), cudaMemcpyHostToDevice));

    for (int i = 0; i < matrix.ncols+1; i++)
    {
	colstream >> matrix.l_cols[i];
    }
    colstream.close();

    NTL::ZZ max_col = NTL::to_ZZ(0);
    for (int i = 0; i < matrix.ncols; i++)
    {
	NTL::ZZ this_col = NTL::to_ZZ(0);
	for (int j = matrix.l_cols[i]; j < matrix.l_cols[i+1]; j++)
	{
	    this_col += to_ZZ(matrix.l_vals[j]);
	}
	max_col = (max_col > this_col) ? max_col : this_col;
    }
    max_col *= (modulus-1);
    max_col >>= (2*BITS_PER_LIMB*LIMBS_PER_UINTX);
    max_overflow = (uint)trunc_long(max_col, 32);

    gpuErrchk(cudaMemcpy(matrix.d_cols, matrix.l_cols,
	(matrix.ncols+1) * sizeof(uint), cudaMemcpyHostToDevice));
}

template <typename T>
void freeMatrix(struct SparseMatrix<T> & matrix)
{
    free(matrix.l_vals);
    free(matrix.l_rows);
    free(matrix.l_cols);
    gpuErrchk(cudaFree(matrix.d_vals));
    gpuErrchk(cudaFree(matrix.d_cols));
    gpuErrchk(cudaFree(matrix.d_rows));
}

template <typename T>
void initBarrett(const NTL::ZZ & modulus_zz, BarrettParams<T> & barrett,
	const uint max_overflow)
{
    barrett.l_modulus = modulus_zz;
    barrett.l_mu = NTL::power2_ZZ(2*BITS_PER_LIMB*LIMBS_PER_UINTX) / modulus_zz;
    T modulus;
    to_uint<T>(modulus_zz, modulus);
    uintXp<T> mu;
    to_uint<uintXp<T>>(barrett.l_mu, mu);

    barrett.l_subtrahends.SetLength(max_overflow+1);
    T * subtrahends = (T *)malloc((max_overflow+1) * sizeof(T));
    for (int i = 0; i <= max_overflow; ++i)
    {
	barrett.l_subtrahends[i] = ((NTL::to_ZZ(i)
	    << (2*BITS_PER_LIMB*LIMBS_PER_UINTX)) / modulus_zz) * modulus_zz;
	NTL::BytesFromZZ((unsigned char *)&subtrahends[i],
	    barrett.l_subtrahends[i], LIMBS_PER_UINTX * sizeof(uint));
    }

    gpuErrchk(cudaMalloc((void**)&barrett.d_modulus, sizeof(T)));
    gpuErrchk(cudaMemcpy(barrett.d_modulus, &modulus, sizeof(T),
	cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_modulus, &barrett.d_modulus, sizeof(T *)));

    gpuErrchk(cudaMalloc((void**)&barrett.d_mu, sizeof(uintXp<T>)));
    gpuErrchk(cudaMemcpy(barrett.d_mu, &mu, sizeof(uintXp<T>),
	cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_mu, &barrett.d_mu, sizeof(uintXp<T> *)));

    gpuErrchk(cudaMalloc((void**)&barrett.d_subtrahends,
	(max_overflow+1) * sizeof(T)));
    gpuErrchk(cudaMemcpy(barrett.d_subtrahends, subtrahends,
	(max_overflow+1) * sizeof(T), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(d_subtrahends, &barrett.d_subtrahends,
	sizeof(T *)));

    free(subtrahends);
}

template<typename T>
void freeBarrett(struct BarrettParams<T> & barrett)
{
    barrett.l_subtrahends.kill();
    barrett.l_modulus.kill();
    barrett.l_mu.kill();
    gpuErrchk(cudaFree(barrett.d_modulus));
    gpuErrchk(cudaFree(barrett.d_mu));
    gpuErrchk(cudaFree(barrett.d_subtrahends));
}
