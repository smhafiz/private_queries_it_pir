// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>
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

#include <vec_ZZ_p.h>
#include <string.h>
#include <stdio.h>
#include <set>
#include "itserver.h"
#include "xor.h"
#include <ctime>
#include <sys/time.h>
//#include "chrono_io"

#ifdef SPIR_SUPPORT
/*
#include "../PolyCommit/PolyCommitCommon.h"
#include "pspir_crypt.h"
*/
#include "spirserver.h"
#endif

// for regular itpir comment all nonzero_columns below
// dbsize_t nonzero_columns = 1004128;   // uncomment for twitter, all, all attr. and essen. attr. (row 1)
// dbsize_t nonzero_columns = 6297;      // uncomment for twitter, filt., all attr. and essen. attr. (row 2)
// dbsize_t nonzero_columns = 58975;     // uncomment for mimic3, all, all attr. and essen. attr. (row 3)
// dbsize_t nonzero_columns = 124026;    // uncomment for mimic3, filt., all attr. and essen. attr. (row 4)
// dbsize_t nonzero_columns = 150345;    // uncomment for yelp, all, all attr. and essen. attr. (row 5)
// dbsize_t nonzero_columns = 601383;    // uncomment for yelp, aug., all attr. and essen. attr. (row 6)

////// ZZ_p Server //////

PercyServer_ZZ_p::PercyServer_ZZ_p (DataStore * datastore, 
	const PercyServerParams * params, PercyStats * stats)
:
    PercyServer(datastore, params, stats),
    params(static_cast<const ZZ_pParams*>(params->percy_params()))
{}

PercyServer_ZZ_p::~PercyServer_ZZ_p ()
{}

// TODO: Fix SPIR to match new function call.
bool PercyServer_ZZ_p::handle_request_impl (std::vector<unsigned char*> requests,
	std::vector<unsigned char*> responses)
{
    //std::cerr<<"itserver.cc line 51 11-08-2022\n";
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }

    // Set the appropriate modulus
    unsigned int modulus_bytes;
    if (params->hybrid()) {
        params->mod_modulussq();
        modulus_bytes = params->modulussq_bytes();
    } else {
        params->mod_modulus();
        modulus_bytes = params->modulus_bytes();
    }

    // Read some values from the params
    bool hybrid_protection = params->hybrid();
    bool tau_independent = params->tau();
    dbsize_t words_per_block = params->words_per_block();
    dbsize_t num_blocks = params->num_blocks();
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t virtual_block_size = params->virtual_block_size();
    dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
    dbsize_t bytes_per_word = (tau_independent ? params->modulus_bytes() : 
	    params->word_size() / 8);

#ifdef USE_W160_OPT
    // Use a faster method if w160 is selected and there is only one query
    bool use_w160_method = !tau_independent && bytes_per_word == 20 
	&& modulus_bytes == 21 && num_queries == 1 && num_blocks < (1ULL<<32)
	&& virtual_block_size == 1;
#endif

#ifdef SPIR_SUPPORT
    PolyCommitParams * pcparamsp = params->get_pcparamsp();

    //Save ZZ_p context
    ZZ_pContext savectx;
    if (params->spir()) {
        savectx.save();
        ZZ_p::init(pcparamsp->get_order());
    }

    // Read in the servers' shared secret key from a file. This too should probably 
    // be read in from somewhere else. Later we combine this with a hash to seed 
    // NTL's PRNG.
    unsigned char prng_seed[24];
    ifstream infile2("server_key.txt");
    if(!infile2.is_open()) {
        std::cerr << "Error: Cannot open server key file." << std::endl;
        exit(1);
    }
    infile2.read((char*)prng_seed, 16);
    infile2.close();

    SPIRServerQuery * spir_queries = new SPIRServerQuery[num_queries];
#endif

    //Read in the input
    vec_ZZ_p * inputvector = new vec_ZZ_p[num_queries]; //Used for non-w160, and w160 with SPIR
    for(nqueries_t q = 0; q < num_queries; ++q) {
	unsigned char * query_request = requests[q];
        //std::cerr << "Receiving query " << (q+1) << " of " << num_queries << " from client...";
#ifdef SPIR_SUPPORT
        if (params->spir()) {
	    bool ret = spir_queries[q].read_spir_input(params, is);
	    if (!ret) {
		delete[] spir_queries;
		return false;
	    }
        }
#endif
#ifdef USE_W160_OPT
        //If the w160 method is used and SPIR is 
        //being used, copy request data into the inputvector as well
        if (!use_w160_method || (use_w160_method && params->spir())) {
#endif
	    inputvector[q].SetLength(num_virtual_blocks);
	    for(dbsize_t i = 0; i < num_virtual_blocks; ++i) {
		ZZ inputz;
		ZZFromBytes(inputz, query_request + i*modulus_bytes, modulus_bytes);
		inputvector[q][i] = to_ZZ_p(inputz);
	    }
#ifdef USE_W160_OPT
        }
#endif
    }

#ifdef SPIR_SUPPORT
    if (params->spir()) {
        for(unsigned int q = 0; q < num_queries; ++q) {
	    bool ret = spir_queries[q].query_verification(params, inputvector[q]);
	    if (!ret) {
		delete[] spir_queries;
		return false;
	    }
        }

        for(unsigned int q = 0; q < num_queries; ++q) {
	    spir_queries[q].init_randomization(params, prng_seed);
        }

        savectx.restore();
    }
#endif

#ifdef USE_W160_OPT
    // Compute the output vector and send it back to the client
    if (use_w160_method) {
        // Get a pointer to the database data
        const unsigned char *dbdata = datastore->get_data();

        // Initialize the output to 0.  Each word of the output will
        // temporarily be stored as 5 128-bit values, and we'll do the
        // modular reduction right at the end.  The 5 128-bit values
        // represent the number
        // a[0] + a[1]<<56 + a[2]<<112 + a[3]<<168 + a[4]<<224.
        __uint128_t out128s[5*words_per_block];
        memset(out128s, 0, sizeof(out128s));

        // Process each block of the database
        unsigned char *multdata = request[0];
        for (dbsize_t j = 0; j < num_blocks;
                multdata += modulus_bytes, ++j) {
            __uint128_t multiplier[3];
            multiplier[0] = (*(uint64_t*)multdata) & 0x00ffffffffffffffULL;
            multiplier[1] = (*(uint64_t*)(multdata+7)) & 0x00ffffffffffffffULL;
            multiplier[2] = (*(uint64_t*)(multdata+14)) & 0x00ffffffffffffffULL;

            __uint128_t *outval = out128s;

            // Process each word in the block
            for (dbsize_t c = 0; c < words_per_block-1;
                    dbdata += bytes_per_word, outval += 5, ++c) {
                uint64_t dbval[3];
                dbval[0] = (*(uint64_t*)dbdata) & 0x00ffffffffffffffULL;
                dbval[1] = (*(uint64_t*)(dbdata+7)) & 0x00ffffffffffffffULL;
                dbval[2] = (*(uint64_t*)(dbdata+14)) & 0x0000ffffffffffffULL;

                // outval += multiplier * dbval
                outval[0] += multiplier[0] * dbval[0];
                outval[1] += multiplier[0] * dbval[1]
                    + multiplier[1] * dbval[0];
                outval[2] += multiplier[0] * dbval[2]
                    + multiplier[1] * dbval[1]
                    + multiplier[2] * dbval[0];
                outval[3] += multiplier[1] * dbval[2]
                    + multiplier[2] * dbval[1];
                outval[4] += multiplier[2] * dbval[2];

                // Do the carries
                outval[1] += (outval[0]>>56);
                outval[0] &= 0x00ffffffffffffffULL;
                outval[2] += (outval[1]>>56);
                outval[1] &= 0x00ffffffffffffffULL;
                outval[3] += (outval[2]>>56);
                outval[2] &= 0x00ffffffffffffffULL;
                outval[4] += (outval[3]>>56);
                outval[3] &= 0x00ffffffffffffffULL;
            }
            // Last block: we may not be able to read past the end of
            // the database
            {
                uint64_t dbval[3];
                dbval[0] = (*(uint64_t*)dbdata) & 0x00ffffffffffffffULL;
                dbval[1] = (*(uint64_t*)(dbdata+7)) & 0x00ffffffffffffffULL;
                dbval[2] = 0;
                memmove(dbval+2, dbdata+14, 6);

                // outval += multiplier * dbval
                outval[0] += multiplier[0] * dbval[0];
                outval[1] += multiplier[0] * dbval[1]
                    + multiplier[1] * dbval[0];
                outval[2] += multiplier[0] * dbval[2]
                    + multiplier[1] * dbval[1]
                    + multiplier[2] * dbval[0];
                outval[3] += multiplier[1] * dbval[2]
                    + multiplier[2] * dbval[1];
                outval[4] += multiplier[2] * dbval[2];

                // Do the carries
                outval[1] += (outval[0]>>56);
                outval[0] &= 0x00ffffffffffffffULL;
                outval[2] += (outval[1]>>56);
                outval[1] &= 0x00ffffffffffffffULL;
                outval[3] += (outval[2]>>56);
                outval[2] &= 0x00ffffffffffffffULL;
                outval[4] += (outval[3]>>56);
                outval[3] &= 0x00ffffffffffffffULL;

                dbdata += bytes_per_word;
            }
        }

        // Now convert the 5 128-bit values in each output word into a
        // single ZZ, reduce it mod the modulus, and output the
        // resulting 21-byte value
        __uint128_t *outval = out128s;
        //std::cerr << "Sending response to client...";
	vec_ZZ_p response;
	response.SetLength(words_per_block);
        for (dbsize_t c = 0; c < words_per_block; outval += 5, ++c) {
            unsigned char valbuf[4*7+16];
            memmove(valbuf, (unsigned char*)(outval), 7);
            memmove(valbuf+7, (unsigned char*)(outval+1), 7);
            memmove(valbuf+14, (unsigned char*)(outval+2), 7);
            memmove(valbuf+21, (unsigned char*)(outval+3), 7);
            memmove(valbuf+28, (unsigned char*)(outval+4), 16);
            response[c] = to_ZZ_p(ZZFromBytes(valbuf,sizeof(valbuf)));
	}

#ifdef SPIR_SUPPORT
	if (params->spir()) {
	    spir_queries[0].randomize_response(response);
	}
#endif

	for (dbsize_t c = 0; c < words_per_block; ++c) {
            BytesFromZZ(response + c * modulus_bytes, rep(response[c]), 
		    modulus_bytes);
        }
    } else {
#endif
	
	vec_ZZ_p * responses_vec = new vec_ZZ_p[num_queries];
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    responses_vec[q].SetLength(words_per_virtual_block);
	}

	const unsigned char *data = datastore->get_data();

	dbsize_t last_vblock_words = words_per_block *
		((num_blocks - 1) % virtual_block_size + 1);
        //std::cerr<<"itserver.cc line 285 11-08-2022\n";
        //timer start
        //note: chose most accurate timer
        compute_all(data, responses_vec, hybrid_protection, num_virtual_blocks,
                words_per_virtual_block, bytes_per_word, last_vblock_words, 
		num_queries, inputvector);
        //timer end
        //print time in ms
#ifdef SPIR_SUPPORT
	if (params->spir()) {
	    for (nqueries_t q = 0; q < num_queries; ++q) {
		spir_queries[q].randomize_response(responses_vec[q]);
	    }
	}
#endif

	for (nqueries_t q = 0; q < num_queries; ++q) {
	    unsigned char * response_location = responses[q];
	    for (dbsize_t c = 0; c < words_per_virtual_block; ++c) {
		BytesFromZZ(response_location, rep(responses_vec[q][c]), 
			modulus_bytes);
		response_location += modulus_bytes;
	    }
	}

	delete[] responses_vec;
#ifdef USE_W160_OPT
    }
#endif

#ifdef SPIR_SUPPORT
    delete[] spir_queries;
#endif
    delete[] inputvector;

    return true;
}

void PercyServer_ZZ_p::compute_all(const unsigned char * data, 
	vec_ZZ_p * responses, bool hybrid_protection, dbsize_t num_blocks, 
	dbsize_t words_per_block, dbsize_t bytes_per_word, 
	dbsize_t last_block_words, nqueries_t num_queries, 
	const vec_ZZ_p * inputvector)
{
    //std::cerr<<"itserver.cc line 326 11-08-2022\n";
    if(strassen_max_depth == PercyServer::STRASSEN_OPTIMAL) {
	strassen_max_depth = optimal_strassen_depth(num_queries,
		num_blocks, words_per_block, bytes_per_word);
    }
    
    if( hybrid_protection || strassen_max_depth == 0) {
        timeval t1, t2;
        double elapsedTime;
        gettimeofday(&t1, NULL); // required for all cases
        for (dbsize_t c = 0, skip_c = 0; c < words_per_block; ++c) { //words_per_block = # of columns
            
            ZZ_p * reswrds = new ZZ_p[num_queries];
            compute_one(data, reswrds, hybrid_protection, 
            num_blocks - (c < last_block_words ? 0 : 1),
                    words_per_block, bytes_per_word, num_queries, 
            inputvector, c);
            for (unsigned int q = 0; q < num_queries; ++q) {
                responses[q][c] = reswrds[q];
            }
            delete[] reswrds;
        }
        gettimeofday(&t2, NULL); // required for all cases
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
        //std::cerr<<"time: " << elapsedTime <<" miliseconds\n";
        std::cerr<<elapsedTime<<"\n";
    } else {
#ifdef STRASSEN_COUNT_OPERATIONS
	    nr_multiplications = 0;
	    nr_additions = 0;
#endif
	    compute_mult_strassen(data, responses, num_blocks, words_per_block, 
		    bytes_per_word, last_block_words, num_queries, inputvector);
#ifdef STRASSEN_COUNT_OPERATIONS
	    std::cerr << "Multiplications: " << nr_multiplications <<
		", additions: " << nr_additions << endl;
#endif
    }
}

nqueries_t PercyServer_ZZ_p::optimal_strassen_depth(
	nqueries_t num_queries, dbsize_t num_blocks,
	dbsize_t words_per_block, dbsize_t bytes_per_word)
{
    // These numbers have been determined emperically
    // the algorithm will reduce until the dimensions
    // are no bigger than the given maximums
    nqueries_t min_num_queries = 8;
    nqueries_t max_num_queries = 3;
    nqueries_t max_num_blocks = 4;
    nqueries_t max_words_per_block = 4;

    // Don't do Strassen for super small problems
    if( num_blocks < 32 || words_per_block < 32 ) {
	return 0;
    }

    if( num_blocks <= 512 ) {
	min_num_queries = 6;
	if (bytes_per_word <= 4) {
	    max_num_blocks = 8;
	}
    } else {
	// num_blocks > 512
	if( bytes_per_word <= 4 ) {
	    max_num_queries = 2;
	    max_num_blocks = 64;
	}
    }

    if( num_queries < min_num_queries ) {
	return 0;
    }

    nqueries_t depth = 0;
    while(num_queries > max_num_queries &&
	    num_blocks > max_num_blocks &&
	    words_per_block > max_words_per_block) {
	num_queries /= 2;
	num_blocks /= 2;
	words_per_block /= 2;

	depth++;
    }

    return depth;
}

void PercyServer_ZZ_p::compute_mult_strassen(const unsigned char * data, 
	vec_ZZ_p * responses, dbsize_t num_blocks, dbsize_t words_per_block,
	dbsize_t bytes_per_word, dbsize_t last_block_words, 
	nqueries_t num_queries, const vec_ZZ_p * inputvector)
{
	vec_ZZ_p * data_matrix = new vec_ZZ_p[words_per_block];
	for (dbsize_t c = 0; c < words_per_block; ++c) {
	    data_matrix[c].SetLength(num_blocks);
	}

        ZZ tmp_wrd;
        for(dbsize_t r = 0; r < num_blocks-1; ++r) {
            for (dbsize_t c = 0; c < words_per_block; ++c) {
                ZZFromBytes(tmp_wrd, data + (r * words_per_block + c) * bytes_per_word,
                        bytes_per_word);
                data_matrix[c][r] = to_ZZ_p(tmp_wrd);
            }
        }
	// Last block
	dbsize_t r = num_blocks-1;
	for (dbsize_t c = 0; c < last_block_words; ++c) {
	    ZZFromBytes(tmp_wrd, data + (r * words_per_block + c) * bytes_per_word,
		    bytes_per_word);
	    data_matrix[c][r] = to_ZZ_p(tmp_wrd);
	}
	for (dbsize_t c = last_block_words; c < words_per_block; ++c) {
	    data_matrix[c][r] = to_ZZ_p(ZZ::zero());
	}

        strassen_split(inputvector, num_queries, data_matrix, words_per_block, responses, 0);
}

// matrix_a[0] is the first row
// matrix_b_transp[0] is the first column!!!
void PercyServer_ZZ_p::matrix_mult_naive(const vec_ZZ_p * matrix_a,
	dbsize_t num_rows_a, const vec_ZZ_p * matrix_b_transp,
	dbsize_t num_cols_b, vec_ZZ_p * matrix_res,
	dbsize_t res_row_offset, dbsize_t res_col_offset,
	dbsize_t inner_offset, dbsize_t inner_dim)
{
    ZZ_p tmp;

    if(inner_dim == 0) {
	inner_dim = matrix_a[0].length();
    }

    for(dbsize_t r = 0; r < num_rows_a; ++r) {
        for(dbsize_t c = 0; c < num_cols_b; ++c) {
            for(dbsize_t i = inner_offset; i < inner_dim; ++i) {
                mul(tmp, matrix_a[r][i], matrix_b_transp[c][i]);
                matrix_res[r + res_row_offset][c + res_col_offset] +=
		    tmp;
            }
        }
    }

#ifdef STRASSEN_COUNT_OPERATIONS
    nr_multiplications += num_rows_a * num_cols_b * inner_dim;
    nr_additions += num_rows_a * num_cols_b * inner_dim;
#endif
}

// matrix_a[0] is the first row
// matrix_b_transp[0] is the first column!!!
void PercyServer_ZZ_p::strassen_split(const vec_ZZ_p * matrix_a,
	dbsize_t num_rows_a, const vec_ZZ_p * matrix_b_transp,
	dbsize_t num_cols_b,
	vec_ZZ_p * matrix_res, nqueries_t depth)
{
    dbsize_t q, r, s, qrem, rrem, srem;
    dbsize_t inner_dimension = matrix_a[0].length();

    // Stop Strassen processes when any dimension reaches this size
    const dbsize_t CUTOFF_DIM = 2;

    if(depth > strassen_level_reached) {
	strassen_level_reached = depth;
    }

    if( num_rows_a <= CUTOFF_DIM || num_cols_b <= CUTOFF_DIM || inner_dimension <= CUTOFF_DIM || depth >= strassen_max_depth ) {
	matrix_mult_naive(matrix_a, num_rows_a, matrix_b_transp, num_cols_b, matrix_res);
	return;
    }

    qrem = num_rows_a % 2;
    rrem = inner_dimension % 2;
    srem = num_cols_b % 2;

    q = num_rows_a - qrem;
    r = inner_dimension - rrem;
    s = num_cols_b - srem;

    // A21 * B11
    if(qrem > 0) {
	matrix_mult_naive( matrix_a + q, qrem, matrix_b_transp, s,
		matrix_res, q, 0, 0, r);
    }

    // A11 * B12
    if(srem > 0) {
	matrix_mult_naive( matrix_a, q, matrix_b_transp + s, srem,
		matrix_res, 0, s, 0, r);
    }

    // A21 * B12
    if(qrem > 0 && srem > 0) {
	matrix_mult_naive( matrix_a + q, qrem, matrix_b_transp + s, srem,
		matrix_res, q, s, 0, r);
    }

    if(rrem > 0) {
	// A12 * B21
	matrix_mult_naive(matrix_a, q, matrix_b_transp, s,
		matrix_res, 0, 0, r);

	// A22 * B21
	if(qrem > 0) {
	    matrix_mult_naive( matrix_a + q, qrem, matrix_b_transp, s,
		    matrix_res, q, 0, r);
	}

	// A12 * B22
	if(srem > 0) {
	    matrix_mult_naive(matrix_a, q, matrix_b_transp + s, srem,
		    matrix_res, 0, s, r);
	}

	// A22 * B22
	if(qrem > 0 && srem > 0) {
	    matrix_mult_naive( matrix_a + q, qrem, matrix_b_transp + s,
		    srem, matrix_res, q, s, r);
	}
    }

    strassen_mult(matrix_a, num_rows_a, matrix_b_transp, num_cols_b,
	    matrix_res, depth);
}

void PercyServer_ZZ_p::strassen_mult(const vec_ZZ_p * matrix_a,
	dbsize_t num_rows_a, const vec_ZZ_p * matrix_b_transp,
	dbsize_t num_cols_b, vec_ZZ_p * matrix_res, nqueries_t depth)
{
    dbsize_t inner_dimension = matrix_a[0].length();

    vec_ZZ_p * matrix_left = new vec_ZZ_p[num_rows_a / 2];
    for (dbsize_t r = 0; r < num_rows_a/2; ++r) {
        matrix_left[r].SetLength(inner_dimension/2);
    }

    vec_ZZ_p * matrix_right = new vec_ZZ_p[num_cols_b / 2];
    for (dbsize_t c = 0; c < num_cols_b/2; ++c) {
        matrix_right[c].SetLength(inner_dimension/2);
    }

    vec_ZZ_p * matrix_subresult = new vec_ZZ_p[num_rows_a / 2];
    for (dbsize_t r = 0; r < num_rows_a/2; ++r) {
	matrix_subresult[r].SetLength(num_cols_b/2);
    }

#ifdef STRASSEN_COUNT_OPERATIONS
    nr_additions += (num_rows_a / 2) * (inner_dimension / 2) * 5;
    nr_additions += (num_cols_b / 2) * (inner_dimension / 2) * 5;
    nr_additions += (num_rows_a / 2) * (num_cols_b / 2) * 10;
#endif

    const ZZ_p &zero = ZZ_p::zero();

    for(int i = 1; i <= 7; ++i) {
        // Fill up left matrix
        for(dbsize_t r1 = 0; r1 < num_rows_a / 2; ++r1) {
            dbsize_t r2 = r1 + num_rows_a / 2;
            for(dbsize_t c1 = 0; c1 < inner_dimension/2; ++c1) {
                dbsize_t c2 = c1 + inner_dimension/2;
                switch(i) {
                case 1: add(matrix_left[r1][c1], matrix_a[r1][c1], matrix_a[r2][c2]); break;
                case 2: add(matrix_left[r1][c1], matrix_a[r2][c1], matrix_a[r2][c2]); break;
                case 3: matrix_left[r1][c1] = matrix_a[r1][c1]                      ; break;
                case 4: matrix_left[r1][c1] = matrix_a[r2][c2]                      ; break;
                case 5: add(matrix_left[r1][c1], matrix_a[r1][c1], matrix_a[r1][c2]); break;
                case 6: sub(matrix_left[r1][c1], matrix_a[r2][c1], matrix_a[r1][c1]); break;
                case 7: sub(matrix_left[r1][c1], matrix_a[r1][c2], matrix_a[r2][c2]); break;
                }
            }
        }
        // Fill up right matrix
        for(dbsize_t c1 = 0; c1 < num_cols_b/2; ++c1) {
            dbsize_t c2 = c1 + num_cols_b/2;
            for(dbsize_t r1 = 0; r1 < inner_dimension/2; ++r1) {
                dbsize_t r2 = r1 + inner_dimension/2;
                switch(i) {
                case 1: add(matrix_right[c1][r1], matrix_b_transp[c1][r1], matrix_b_transp[c2][r2]); break;
                case 2: matrix_right[c1][r1] = matrix_b_transp[c1][r1]                             ; break;
                case 3: sub(matrix_right[c1][r1], matrix_b_transp[c2][r1], matrix_b_transp[c2][r2]); break;
                case 4: sub(matrix_right[c1][r1], matrix_b_transp[c1][r2], matrix_b_transp[c1][r1]); break;
                case 5: matrix_right[c1][r1] = matrix_b_transp[c2][r2]                             ; break;
                case 6: add(matrix_right[c1][r1], matrix_b_transp[c1][r1], matrix_b_transp[c2][r1]); break;
                case 7: add(matrix_right[c1][r1], matrix_b_transp[c1][r2], matrix_b_transp[c2][r2]); break;
                }
            }
        }

	// Calculate subresult
	strassen_split(matrix_left, num_rows_a/2, matrix_right,
		num_cols_b/2, matrix_subresult, depth + 1);

	for(dbsize_t r1 = 0; r1 < num_rows_a/2; ++r1) {
	    dbsize_t r2 = r1 + num_rows_a/2;
	    for(dbsize_t c1 = 0; c1 < num_cols_b/2; ++c1) {
		dbsize_t c2 = c1 + num_cols_b/2;
		switch(i) {
		case 1:
		    add(matrix_res[r1][c1], matrix_res[r1][c1], matrix_subresult[r1][c1]);
		    add(matrix_res[r2][c2], matrix_res[r2][c2], matrix_subresult[r1][c1]);
		    break;
		case 2:
		    add(matrix_res[r2][c1], matrix_res[r2][c1], matrix_subresult[r1][c1]);
		    sub(matrix_res[r2][c2], matrix_res[r2][c2], matrix_subresult[r1][c1]);
		    break;
		case 3:
		    add(matrix_res[r1][c2], matrix_res[r1][c2], matrix_subresult[r1][c1]);
		    add(matrix_res[r2][c2], matrix_res[r2][c2], matrix_subresult[r1][c1]);
		    break;
		case 4:
		    add(matrix_res[r1][c1], matrix_res[r1][c1], matrix_subresult[r1][c1]);
		    add(matrix_res[r2][c1], matrix_res[r2][c1], matrix_subresult[r1][c1]);
		    break;
		case 5:
		    sub(matrix_res[r1][c1], matrix_res[r1][c1], matrix_subresult[r1][c1]);
		    add(matrix_res[r1][c2], matrix_res[r1][c2], matrix_subresult[r1][c1]);
		    break;
		case 6:
		    add(matrix_res[r2][c2], matrix_res[r2][c2], matrix_subresult[r1][c1]);
		    break;
		case 7:
		    add(matrix_res[r1][c1], matrix_res[r1][c1], matrix_subresult[r1][c1]);
		    break;
		}

		matrix_subresult[r1][c1] = zero;
	    }
	}
    }

    delete[] matrix_left;
    delete[] matrix_right;
    delete[] matrix_subresult;
}

void PercyServer_ZZ_p::compute_one(const unsigned char * data, ZZ_p * value, 
	bool hybrid_protection, dbsize_t num_blocks, dbsize_t words_per_block,
	dbsize_t bytes_per_word, nqueries_t num_queries, 
	const vec_ZZ_p * inputvector, dbsize_t c)
{

    for(nqueries_t q = 0; q < num_queries; ++q) {
        value[q] = hybrid_protection ? 1 : 0;
    }

    ZZ_p prod;
    ZZ wrd_raw;
    ZZ_p wrd;
    
    for (dbsize_t j = 0, skip_c=0; j < num_blocks; ++j) {
        // The cth word of the jth block

        // if (j >= nonzero_columns){continue;} // uncomment for all attr. and essen. attr.

	ZZFromBytes(wrd_raw, data + (j * words_per_block + c) * bytes_per_word,
		bytes_per_word);

        if (hybrid_protection) {
            for(nqueries_t q = 0; q < num_queries; ++q) {
                value[q] = value[q] * power(inputvector[q][j], wrd_raw);
            }
        } else {
	    wrd = to_ZZ_p(wrd_raw);
            for(nqueries_t q = 0; q < num_queries; ++q) {
                //value[q] = value[q] + inputvector[q][j] * to_ZZ_p(wrd);
		mul(prod, inputvector[q][j], wrd);
		add(value[q], value[q], prod);
            }
        }
    }

    if (byzantine) {
        // Produce a *consistent* incorrect value for maximal client
        // confusion.
        for(nqueries_t q = 0; q < num_queries; ++q) {
            value[q] += 1;
        }
    }
}

void PercyServer_ZZ_p::combine_results (unsigned char * result, 
	std::vector<unsigned char*> worker_results)
{
    nservers_t num_workers = worker_results.size();
    dbsize_t word_bytes = params->modulus_bytes();
    dbsize_t num_words = params->words_per_block() * params->virtual_block_size();
    vec_ZZ_p result_ZZ_p;
    result_ZZ_p.SetLength(num_words);
    for (nservers_t i = 0; i < num_workers; ++i) {
	unsigned char * wresult = worker_results[i];
	for (dbsize_t w = 0; w < num_words; ++w) {
	    result_ZZ_p[w] += to_ZZ_p(
		    ZZFromBytes(wresult + w * word_bytes, word_bytes));
	}
    }
    for (dbsize_t w = 0; w < num_words; ++w) {
	BytesFromZZ(result + w * word_bytes, rep(result_ZZ_p[w]), word_bytes);
    }
}


////// Chor Server //////

PercyServer_Chor::PercyServer_Chor (DataStore * datastore, 
	const PercyServerParams * params, PercyStats * stats)
:
    PercyServer(datastore, params, stats),
    params(static_cast<const ChorParams*>(params->percy_params()))
{}

PercyServer_Chor::~PercyServer_Chor ()
{}

bool PercyServer_Chor::handle_request_impl (std::vector<unsigned char*> requests,
	std::vector<unsigned char*> responses)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }

    // Read some values from the params
    dbsize_t block_size = params->block_size();
    dbsize_t num_bytes = (params->num_virtual_blocks() - 1) / 8 + 1;
    dbsize_t virtual_block_bytes = block_size * params->virtual_block_size();
    nqueries_t q;


    /*
    std::cerr << "Server " << ": ";
    printBS_debug(inputvector, num_blocks/8);
    std::cerr << std::endl;
    */

    const unsigned char *data = datastore->get_data();

    // Compute the output vector and send it back to the client
#ifdef VERBOSE_CHOR
    struct timeval ts, te;
    gettimeofday(&ts, NULL);
#endif
    for (q = 0; q < num_queries; q++) {
	unsigned char * query_request = requests[q];
	unsigned char * query_response = responses[q];
	memset(query_response, '\0', virtual_block_bytes);
	for (unsigned int i = 0; i < num_bytes; i++) {
	    unsigned char query_byte = query_request[i];
	    if (query_byte & 1)
		XOR_equal(query_response, data + 8*i*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 2)
		XOR_equal(query_response, data + (8*i+1)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 4)
		XOR_equal(query_response, data + (8*i+2)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 8)
		XOR_equal(query_response, data + (8*i+3)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 16)
		XOR_equal(query_response, data + (8*i+4)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 32)
		XOR_equal(query_response, data + (8*i+5)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 64)
		XOR_equal(query_response, data + (8*i+6)*virtual_block_bytes, 
			virtual_block_bytes);
	    if (query_byte & 128)
		XOR_equal(query_response, data + (8*i+7)*virtual_block_bytes, 
			virtual_block_bytes);
	}
	if (byzantine) {
	    for (dbsize_t i = 0; i < virtual_block_bytes; ++i) {
		query_response[i]++;
	    }
	}
    }

#ifdef VERBOSE_CHOR
    gettimeofday(&te, NULL);
    int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
    fprintf(stderr, "Chor: %d.%03d msec computation\n", td/1000, td%1000);
#endif

    /*
    std::cerr << "Server: ";
    printBS_debug(outputvector+(num_queries-1)*block_size, block_size);
    std::cerr << std::endl;
    */

    return true;
}

void PercyServer_Chor::combine_results (unsigned char * result, 
	std::vector<unsigned char*> worker_results)
{
    dbsize_t block_size = params->block_size() * params->virtual_block_size();
    for (nservers_t i = 0; i < worker_results.size(); ++i) {
	XOR_equal(result, worker_results[i], block_size);
    }
}

