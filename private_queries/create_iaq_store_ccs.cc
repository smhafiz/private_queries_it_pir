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

#include <iostream>
#include <stdlib.h>
#include <set>
#include <fstream>

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <random>
#include <chrono>
#include <ratio>
#include <thread>
#include <mutex>

#include <cstring>


static const int NumberOfTrials = 1;


void printMatrix(int **a,int r, int c);
void printVector(unsigned *a,int c);
void printVector2(unsigned long long int *a,int c);

int main(int argc, char ** argv)
{

	//std::ofstream myfile;
	std::ofstream ccsRow, ccsCol;
	//myfile.open ("ccs-out.txt");
	uint64_t seed = time(NULL);
	//myfile << "Seed: " << seed ;
//	std::cout << "Seed: " << seed ;
	srand(seed);
	char * infile = new char[strlen(argv[1]) + 4];
	sprintf(infile, "%s.row", argv[1]);
	ccsRow.open (infile);
	char * infile2 = new char[strlen(argv[1]) + 4];
	sprintf(infile2, "%s.col", argv[1]);
	ccsCol.open (infile2);
	
	long p,r,u;
	p = 1L << (atoi(argv[2]));
	r = 1L << (atoi(argv[3]));
	u = atoi(argv[4]);
	//std::cout << "r: " << r << " p: " << p << std::endl ;
	//const long max_u =2, r = 1L << 4;
	//for (long p = 2; p <= r; p <<= 1)
	//{
		//for (long u = 1; u <= max_u; u <<=1)
		//{
			//myfile << "\n\np: " << p << "\tr: " << r << "\tu: " << u ;
			ccsRow << p;
			ccsCol << r;
			for(int trials=0;trials<NumberOfTrials;trials++){

				//top:


				std::chrono::high_resolution_clock::time_point pruTotalTimeStart = std::chrono::high_resolution_clock::now();
				std::vector<std::set<long>> cols;
				std::vector<std::set<long>*> cols2;
				cols.resize(r);
				for (auto it = begin(cols); it != end(cols); ++it) cols2.push_back(&(*it));

				for (long i = 1; i <= p; ++i)
				{
					for (long j = 1; j <= u; )
					{
						long c = rand() % cols2.size();
						if (cols2[c]->size() < u && cols2[c]->insert(i).second)
						{
							j++;
						}
						else
						{
							long a = rand() % r;
							if (cols[a].size() > 0 && cols[a].find(i) == end(cols[a]))
							{
								auto elt = begin(cols[a]);
								std::advance(elt, rand() % cols[a].size());
								long tmp = *elt;
								if (cols2[c]->find(tmp) == end(*(cols2[c])))
								{
									cols[a].erase(elt);
									cols[a].insert(i);
									cols2[c]->insert(tmp);
									j++;
								}
							}
						}
						if (cols2[c]->size() == u) cols2.erase(begin(cols2) + c);
					}
				}

				int numberOfNonZeroElements = p*u;
				int lengthOfColumnPtr = r+1;
			
				int *rowIndex = (int*)malloc(sizeof(int)*numberOfNonZeroElements);
				int *columnPtr = (int*)malloc(sizeof(int)*(lengthOfColumnPtr));

				int t=0;
				//int sum=0;
				columnPtr[0] = 0;
				int lengthOfCPReduced = 0;
				ccsCol << " " << 0;
				for (int i = 0; i < r; ++i)
				{
					for (auto it = begin(cols[i]); it != end(cols[i]); ++it)
					{
						rowIndex[t++] = (*it)-1;
						ccsRow << " " << rowIndex[t-1];
					}
					//if (cols[i].size())
					//{
						columnPtr[lengthOfCPReduced+1]=columnPtr[lengthOfCPReduced]+cols[i].size();
						ccsCol << " " << columnPtr[lengthOfCPReduced+1];
						lengthOfCPReduced++;
					//}
					//sum+=cols[i].size();
				}

				//std::cout << "\n\nCol (" << cols.size() <<"): ";
				std::chrono::high_resolution_clock::time_point pruMatrixCreationEnd = std::chrono::high_resolution_clock::now();
				std::chrono::duration<int,std::milli> pruMatrixCreationDuration = std::chrono::duration_cast<std::chrono::duration<int,std::milli>>(pruMatrixCreationEnd-pruTotalTimeStart);
//				std::cout<<"\nTrial: "<<trials+1<< "\tMatrix Creation: "<< (double) pruMatrixCreationDuration.count()/1000 <<"s";
				//myfile <<"\nTrial: "<<trials+1<< "\tMatrix Creation: "<< (double) pruMatrixCreationDuration.count()/1000 <<"s";
				/*
				 * CUDA started
				 *
				 **/
				std::chrono::high_resolution_clock::time_point pruMatrixCopyStart = std::chrono::high_resolution_clock::now();
				int lengthOfResultVectorReduced = lengthOfCPReduced-1;
			}

		//}
	//}
	//myfile.close();
	ccsCol.close();
	ccsRow.close();
}


void printMatrix(int **a,int r, int c) {
	int i=0,j=0;
	for(;i<r;i++){
		for(j=0;j<c;j++){
			printf("%d ",a[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(unsigned *a,int c) {
	int j=0;
	for(j=0;j<c;j++){
		printf("%d ",a[j]);
	}
	printf("\n");
}

void printVector2(unsigned long long int *a,int c) {
	int j=0;
	for(j=0;j<c;j++){
		printf("%lld ",a[j]);
	}
	printf("\n");
}


