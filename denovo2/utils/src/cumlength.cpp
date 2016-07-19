/*
 * SOLiD Sequence Assembilng Project
 * Copyright 2010 Life Technologies
 *
 * Dima Brinza
 * Staff Scientist, Bioinformatics
 * Life Technologies
 * 850 Lincoln Center Dr. M/S 408-2
 * Foster City, CA 94404
 *
 * office:  800 Building, room 8200-08
 * phone: (650)-554-2052
 * e-mail: dima@appliedbiosystems.com 
 *
 */

/*
   Outputs contig lengths in dicreasing order
  
Build:

g++ -g -O3 cumlength.cpp -o cumlength

Usage:

cumlength <contigs_file> <min_length> > <cum.len.txt>

Input:

 contigs_file  -  fasta file with assembled contigs/scaffolds.
 min_length    -  minimal length of contig to be included in analysis.
 cum.len.txt   -  output file name with sorted in decreasing order contigs sizes.

Output:

 cum.len.txt   -  file with sorted in decreasing order contigs sizes.

*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <algorithm>

using namespace std;
//---------------------------------------------------------------------------------
void GenerateCummulativeContigsLength(string FileName, int minLength);
//---------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  GenerateCummulativeContigsLength(argv[1], atoi(argv[2]));
  return 0;
}
//---------------------------------------------------------------------------------
void GenerateCummulativeContigsLength(string FileName, int minLength)
{
	ifstream infile;
	infile.open(FileName.c_str());

	string line;
    if (!infile) {
        cerr << "Problem with input file " << FileName << endl;
        exit(1);
    }

  string s ="";
  vector<int> contig_length;
  while(getline(infile, line))
      {
       if(line[0]=='>')
       {
        if(s.length()>minLength)contig_length.push_back(s.length());
        s="";
       }
       else
       s+=line;
     }
     
  if(s.length()>minLength)contig_length.push_back(s.length());
  infile.close();     
  
  sort(contig_length.begin(),contig_length.end());
  
  int sum = 0;
  for(int i=contig_length.size()-1;i>=0;i--)
  {
  sum += contig_length[i];
  cout << sum << "\t" << contig_length[i] << endl;
  }
 infile.close();
}
//---------------------------------------------------------------------------------
