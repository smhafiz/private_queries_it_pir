// This file is part of Private expressive queries v0.1.
// 
// 
// Copyright (C) 2016, Ryan Henry and Syed Mahbub Hafiz.
// 
// Private expressive queries is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Private expressive queries is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Private expressive queries. If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <chrono>

#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_vec_ZZ_p.h>

NTL_CLIENT

int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();
    //std::chrono::nanoseconds onesec{1000000000};

    NTL::ZZ modulus = NTL::RandomPrime_ZZ(atoi(argv[1])*8);//to support GF(2^8)
    
    //for(int d=0;d<argc;d++)cout << argv[d] << endl;
    NTL::ZZ_p::init(modulus);
    //int modulus = 65537;
    int u = argc - 3;

    std::ifstream xcoords(argv[2], std::ifstream::in);
    std::ifstream * irow = new std::ifstream[u];
    std::ifstream * icol = new std::ifstream[u];

    // determine how many servers, and what evals they will hold
    int num_servers;
    xcoords >> num_servers;
    NTL::vec_ZZ_p server_xcoord(INIT_SIZE, num_servers);
    for (int i = 0; i < num_servers; ++i)
    {
    	uint16_t tmp;
	xcoords >> tmp;
	server_xcoord[i] = to_ZZ_p(tmp);
    }

    // open the files containing the sparse permutation-like matrices
    for (int i = 1; i <= u; ++i)
    {
	char * infile = new char[strlen(argv[i+2]) + 4];
	sprintf(infile, "%s.row", argv[i+2]);
	irow[i - 1].open(infile, std::ifstream::in);
	sprintf(infile, "%s.col", argv[i+2]);
	icol[i - 1].open(infile, std::ifstream::in);
	delete [] infile;
    }

    // determine the dimensions (and ensure that all files agree!)
    int p, r;
    irow[0] >> p;
    icol[0] >> r;
    for (int i = 1; i < u; ++i)
    {
	int _p,_r;
	irow[i] >> _p;
	icol[i] >> _r;
	if (p != _p || r != _r)
	{
		std::cout << "size mismatch!\n";
		exit(1);
	}
    }

    // open the output files
    std::ofstream orow;
    std::ofstream ocol;
    std::ofstream opoly;
    std::ofstream * oval = new std::ofstream[num_servers];
    orow.open("out.row", std::ofstream::out);
    ocol.open("out.col", std::ofstream::out);
    opoly.open("out.polys", std::ofstream::out);
    char * outfile = new char[9];
    for (int i = 0; i < num_servers; ++i)
    {
	sprintf(outfile, "val.%lu", trunc_long(rep(server_xcoord[i]), 16));
	oval[i].open(outfile, std::ofstream::out);
	oval[i] << modulus << " ";
    }
    delete [] outfile;

    orow << p << " ";
    long opos = orow.tellp();
    orow << "           ";
    ocol << r << " 0 ";
    opoly << modulus << " ";

    // init NTL and precompute Lagrange polynomials
    ZZ_p::init(to_ZZ(modulus));

    ZZ_pX upoly(INIT_SIZE, u);
    vec_ZZ_pX lagrange(INIT_SIZE, u, upoly);
    const ZZ_pX X(INIT_MONO, 1);
    ZZ_pX tmp(INIT_SIZE, u);
    ZZ_p w(INIT_ALLOC);

    vec_ZZ_p interp_xcoords(INIT_SIZE, u);
    ZZ_pX zeros(INIT_SIZE, u + 1);
    for (int i = 0; i < u; ++i) interp_xcoords[i] = to_ZZ_p(i);
    NTL::BuildFromRoots(zeros, interp_xcoords);
    interp_xcoords.kill();

    for (int i = 0; i < u; i++)
    {
	NTL::set(w);
	for (int j = 0; j < u; ++j)
	{
	    if (i == j) continue;
	    w *= (i - j);
	}
	NTL::div(tmp, zeros, (X - i));
	NTL::mul(lagrange[i], tmp, NTL::inv(w));
    }

    vec_vec_ZZ_p precomp(INIT_SIZE, num_servers, vec_ZZ_p(INIT_SIZE, u));
    for (int i = 0; i < num_servers; i++)
    {
    	for (int j = 0; j < u; j++)
    	{
    	    precomp[i][j] = eval(lagrange[j], server_xcoord[i]);
    	}
    }

    // read the input matrices and create the index in CCS format, on the fly
    int * prev_col = new int[u];
    memset(prev_col, 0, sizeof(int)*u);

    int num_vals = 0;
    vec_vec_ZZ_p buffer(INIT_SIZE, num_servers, vec_ZZ_p(INIT_SIZE, p));

    //vec_ZZ_pX polybuf(INIT_SIZE, p, upoly);
    int next_col = 0;
    bool * nz = new bool[p];

    for (int i = 0; i < p; i++) nz[i] = false;
    for (int j = 0; j < r; ++j)
    {
    
	for (int i = 0; i < u; ++i)
	{
        
	    int to_read = prev_col[i];
	    icol[i] >> prev_col[i];
	    to_read = prev_col[i] - to_read;
        
	    while (to_read--)
	    {
		int which_row;
		irow[i] >> which_row;
		//NTL::add(polybuf[which_row], polybuf[which_row], lagrange[i]);

		for (int k = 0; k < num_servers; k++)
		{
            //problem line
            //std::cout <<"line:157, k: " << k << " | i: " << i << " | which_row: " << which_row << " | p: " << p << "\n";
            NTL::add(buffer[k][which_row], buffer[k][which_row], precomp[k][i]);

		}
		nz[which_row] = true;
	    }
        
	}
	for (int i = 0; i < p; ++i)
	{
	    if (nz[i])
	    {
		next_col++;
		orow << i << " ";
		//opoly << polybuf[i] << " ";
		for (int k = 0; k < num_servers; ++k)
		{
		    oval[k] << buffer[k][i] << " ";
		    NTL::clear(buffer[k][i]);
		}
                num_vals++;
	    }
	    nz[i] = false;
	    //NTL::clear(polybuf[i]);
	}
    
	ocol << next_col << " ";
    }
    
    orow.seekp(opos);
    
    orow << num_vals;
    
    // cleanup
    delete [] prev_col;
    for (int i = 0; i < u; ++i)
    {
	irow[i].close();
	icol[i].close();
    }
    for (int i = 0; i < num_servers; ++i)
    {
	oval[i].close();
    }
    delete [] irow;
    delete [] icol;
    delete [] oval;
    orow.close();
    ocol.close();
    opoly.close();
    lagrange.kill();
    buffer.kill();
    delete [] nz;

    //auto time_elapsed = ;
    //std::cout << "Time requires to batch CCS files: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << " milliseconds\n";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "\n";
    return 0;
}
