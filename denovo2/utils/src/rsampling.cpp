/*
 * SOLiD Reads De Novo Assembly Project
 * Copyright (c) 2010 Life Technologies
 *
 * Dumitru Brinza
 * Staff Scientist, Bioinformatics
 * Life Technologies
 * 850 Lincoln Center Dr. M/S 408-2
 * Foster City, CA 94404
 *
 * office:  800 Building, room 8200-08
 * phone: (650)-554-2052
 * e-mail: Dumitru.Brinza@lifetech.com
 *
 */


#define HEADMESSAGE "                                                        \n\
                                                                             \n\
SOLiD reads sampling v.0.2                                                   \n\
Copyright (2010) by Life Technologies                                        \n\
***************************************************                           "
#define HELPMESSAGE "                                                        \n\
Usage:                                                                       \n\
                                                                             \n\
 - for fragment library data run:                                            \n\
                                                                             \n\
 ./rsampling <f3_csfasta> <f3_qual> <refLength> [-options]                   \n\
                                                                             \n\
 - for mate-paired library data run:                                         \n\
                                                                             \n\
 ./rsampling <f3_csfasta> <f3_qual> <refLength> -r3 r3_csfasta -r3qv r3_qual [-options]\n\
                                                                             \n\
Input:                                                                       \n\
                                                                             \n\
                                                                             \n\
 f3_csfasta - csfasta/fasta file with reads.                                 \n\
              in case of mate-paired data this is the file with F3 reads.    \n\
              for fragment data the title of each read is irrelevant.        \n\
              header of the file may contain comments and descriptions.      \n\
                                                                             \n\
 f3_qual    - filename with quality values (if available). notice, that order\n\
              of reads in csfasta file should be the same as in quality value file.\n\
              if file is not available then input \"none\".                  \n\
 refLength  - expected length of sequenced DNA region, e.g., 4600000 for     \n\
              E.Coli 4.6Mb genome.                                           \n\
                                                                             \n\
                                                                             \n\
Options required for mate-paired data:                                       \n\
                                                                             \n\
 -r3 r3_csfasta     csfasta/fasta file with R3 reads.                        \n\
 -r3qv r3_qual      file with quality values for R3 reads (if available). notice, that order \n\
                    of reads in csfasta file should be the same as in quality values file.   \n\
                    do not include this option if quality file is not available or if        \n\
                    f3_qual is \"none\".                                                     \n\
                                                                                             \n\
Output:                                                                                      \n\
                                                                                             \n\
 ./subreads.csfasta   - csfasta file with subsampled F3 reads.                               \n\
 ./subreads.qual      - qual file of subsampled F3 reads.                                    \n\
 ./subreads2.csfasta  - csfasta file with subsampled R3 reads.                               \n\
 ./subreads2.qual     - qual file of subsampled R3 reads.                                    \n\
 ./.param             - file with one text line with information on coverage and read length.\n\
                                                                                             \n\
Options:                                                                                     \n\
                                                                                             \n\
 -outdir dir        outputs files into \"dir\" directory (default \".\").                    \n\
 -maxcov c          indicates the coverage formed by sub-sampled reads (default c=300).      \n\
                    if \"c\" is higher than existing coverage then whole dataset is          \n\
                    considered and value of \"c\" is changed to reflect actual coverage.     \n\
-nosampling         computes coverage and read_length, and excludes sub-sampling.            \n\
-log                include if execution should not be outputted.                            \n"

//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <cstring>

bool LOG = true;

using namespace std;
//---------------------------------------------------------------------------------
void ExtractSubsetOfReads(string csfasta_file, string quality_file, string csfasta_r3, string quality_r3, 
                          string outdir, int reduction_factor);
                          
double get_number_of_bases(string csfasta_file, int & read_length);
//---------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  
  if(argc<4) {   cout << HEADMESSAGE << endl;
                 cout << HELPMESSAGE << endl; exit(1);};
  
// Set default values
  int acnt = 1;
  string csfasta_file = argv[acnt++];
  string quality_file = argv[acnt++];
  unsigned long refLength = atol(argv[acnt++]);
  string outdir = ".";
  bool nosampling = false;
  string csfasta_r3 = "";
  string quality_r3 = ""; 
  unsigned long maxcov = 300;

// Read input options
  while(acnt<argc)
  {
   if(strcmp(argv[acnt], "-outdir")==0) outdir = argv[++acnt];
   else
   if(strcmp(argv[acnt], "-r3")==0) csfasta_r3 = argv[++acnt];
   else
   if(strcmp(argv[acnt], "-r3qv")==0) quality_r3 = argv[++acnt];
   else
   if(strcmp(argv[acnt], "-maxcov")==0) maxcov = atol(argv[++acnt]);
   else
   if(strcmp(argv[acnt], "-nosampling")==0) nosampling = true;   
   else
   if(strcmp(argv[acnt], "-log")==0) LOG = false;
   else {
    cout << HELPMESSAGE << endl; exit(1);}
    ++acnt;
  }

 if(LOG){ cout << HEADMESSAGE << endl; }

// Check feasibility bounds
  if(refLength<100 || refLength>10000000000){
   cout << "ERROR: Expected reference Length= " << refLength << " is out of bounds [100,1000000000]" << endl; 
   exit(1);
   }
  
  if(LOG) cout << "Expected reference Length= " << refLength << endl;   
  
  double reduction_factor = 0;
  int read_length;
  if(!nosampling) {
    reduction_factor = get_number_of_bases(csfasta_file, read_length)/refLength/maxcov;
    if(csfasta_r3!="")
     {
      int read2_length;
      reduction_factor += get_number_of_bases(csfasta_r3, read2_length)/refLength/maxcov;
     }
  }   

  int modselect = (int)round(reduction_factor);
  
  if(LOG) cout << "reduction factor " << modselect << endl;
  
  string ofparam = outdir + "/.param";
	ofstream param;
	param.open(ofparam.c_str());
	if(!param) { cout << " ERROR: Can't create param file " << endl; exit(1);}  

  if(modselect<2) { if(LOG) cout << "resulting coverage " << (int)(maxcov*reduction_factor) << endl; param << 0 << " " << (int)(maxcov*reduction_factor) << " " << read_length; 
  param.close();
  return 0;}
  
  ExtractSubsetOfReads(csfasta_file, quality_file, csfasta_r3, quality_r3, outdir, modselect );
	
  if(modselect>1) { 	if(LOG) cout << "resulting coverage " << maxcov << endl; param << 1 << " " << maxcov << " " << read_length;}
  
  param.close();		
	return 0;

}
//---------------------------------------------------------------------------------
void ExtractSubsetOfReads(string csfasta_file, string quality_file, string csfasta_r3, string quality_r3, string outdir, int reduction_factor)
{

	ifstream infile1, infile2, infile3, infile4;
	infile1.open(csfasta_file.c_str());
	if(quality_file!="none") infile2.open(quality_file.c_str());
	if(csfasta_r3!="") infile3.open(csfasta_r3.c_str());
  if(quality_r3!="") infile4.open(quality_r3.c_str());

	string line1,line2,line3,line4;

  if(!infile1) { cout << "ERROR: " << csfasta_file << " does not exist" << endl; exit(1);}
  if(!infile3) { cout << "ERROR: " << csfasta_r3 << " does not exist" << endl; exit(1);}
  if(quality_file!="none" && !infile2) { cout << "ERROR: " << quality_file << " does not exist" << endl; exit(1);}
  if(quality_r3!="" && !infile4) { cout << "ERROR: " << quality_r3 << " does not exist" << endl; exit(1);}  
  if(csfasta_r3!="" && quality_file!="none" && quality_r3==""){ cout << "ERROR: Please specify quality file for R3" << endl; exit(1);}
  
  
  string title1 = "", title2 ="";
  int count = 0;
  
  string ofreads = outdir + "/subreads.csfasta";
	string ofqual =  outdir + "/subreads.qual";
	string ofreads2 = outdir + "/subreads2.csfasta";
	string ofqual2 = outdir + "/subreads2.qual";

	ofstream freads,fqual,freads2,fqual2;
	freads.open(ofreads.c_str());
	if(!freads) { cout << " ERROR: Can't create " << ofreads << endl; exit(1);}

  if(quality_file!="none") {
                         fqual.open(ofqual.c_str());
                         if(!fqual)   { cout << " ERROR: Can't create " << ofqual << endl; exit(1);}
 }
  if(csfasta_r3!="") {
                         freads2.open(ofreads2.c_str()); 
                       	 if(!freads2) { cout << " ERROR: Can't create " << ofreads2 << endl; exit(1);}
 }
  if(quality_r3!="") {
                         fqual2.open(ofqual2.c_str()); 
                         if(!fqual2)  { cout << " ERROR: Can't create " << ofqual2 << endl; exit(1);}
 }
	 
   while(getline(infile1, line1) && (quality_file=="none" || getline(infile2, line2))  && (csfasta_r3=="" || getline(infile3, line3))  && (quality_r3=="" || getline(infile4, line4)))
      {
       if(line1[0]=='#') while(getline(infile1,line1) && (line1.length()==0 || line1[0]=='#'));
       if(quality_file!="none" && line2[0]=='#') while(getline(infile2,line2) && (line2.length()==0 || line2[0]=='#'));

       if(csfasta_r3!="" && line3[0]=='#') while(getline(infile3,line3) && (line3.length()==0 || line3[0]=='#'));
       if(quality_r3!="" && line4[0]=='#') while(getline(infile4,line4) && (line4.length()==0 || line4[0]=='#'));       
       
       if(line1[0]=='>')
       {
        if(quality_file!="none" && line1!=line2) { cout << "ERROR: "<< csfasta_file << " is not synchronizaed with "<< quality_file << endl; exit(1); }
        if(csfasta_r3!="" && quality_r3!="" && line3!=line4) { cout << "ERROR: "<< csfasta_r3 << " is not synchronizaed with "<< quality_r3 << endl; exit(1); }
        title1 = line1;
        if(csfasta_r3!="")
         {
          title2 = line3;
          string::size_type loc, loc2;
          
          if((loc = title1.find("_",0))==string::npos || (loc = title1.find("_",loc+1))==string::npos || (loc = title1.find("_",loc+1))==string::npos) 
           { cout << "ERROR: "<< csfasta_file << " unsupported format for read title >X_Y_Z_ instead: " << title1 << endl; exit(1);  }
           
          if((loc2 = title2.find("_",0))==string::npos || (loc2 = title2.find("_",loc2+1))==string::npos || (loc2 = title2.find("_",loc2+1))==string::npos) 
           { cout << "ERROR: "<< csfasta_r3 << " unsupported format for read title >X_Y_Z_ instead: " << title2 << endl; exit(1);  }
          
           string t1 = title1.substr(1,loc-1);
           string t2 = title2.substr(1,loc2-1);
           
           if(t1<t2) cout << "WARN: " << title1 << " has missing mate." << endl;
           if(t1>t2) cout << "WARN: " << title2 << " has missing mate." << endl;
           
           bool eof1=false, eof2=false;       
           while(t1!=t2 && (!eof1 || !eof2))
           {
              eof1=true;eof2=true;
              while(t1<t2 && getline(infile1, line1) && (quality_file=="none" || getline(infile2, line2))){
               eof1 = false;
               if(line1[0]=='>')
               {
               if(quality_file!="none" && line1!=line2) { cout << "ERROR: "<< csfasta_file << " is not synchronizaed with "<< quality_file << endl; exit(1); }
               title1 = line1;
               if((loc = title1.find("_",0))==string::npos || (loc = title1.find("_",loc+1))==string::npos || (loc = title1.find("_",loc+1))==string::npos) 
               { cout << "ERROR: "<< csfasta_file << " unsupported format for read title >X_Y_Z_ instead: " << title1 << endl; exit(1);  }
               t1 = title1.substr(1,loc-1);
               if(t1<t2) cout << "WARN: " << title1 << " has missing mate." << endl;
               }
              }
              
               while(t1>t2 && getline(infile3, line3) && (quality_r3=="" || getline(infile4, line4))){
               eof2 = false;
               if(line3[0]=='>')
               {
               if(quality_r3!="none" && line3!=line4) { cout << "ERROR: "<< csfasta_r3 << " is not synchronizaed with "<< quality_r3 << endl; exit(1); }
               title2 = line3;
               if((loc2 = title2.find("_",0))==string::npos || (loc2 = title2.find("_",loc2+1))==string::npos || (loc2 = title2.find("_",loc2+1))==string::npos) 
                { cout << "ERROR: "<< csfasta_r3 << " unsupported format for read title >X_Y_Z_ instead: " << title2 << endl; exit(1);  }
               t2 = title2.substr(1,loc2-1);
               if(t1>t2) cout << "WARN: " << title2 << " has missing mate." << endl;
               }
              }
           }
          
          }
        continue;
       }
       else
       {
        if(reduction_factor<2 || rand()%reduction_factor==0){
        freads << title1 << endl << line1 << endl; 
        if(csfasta_r3!="") freads2 << title2 << endl << line3 << endl;
        if(quality_file!="none") fqual << title1 << endl << line2 << endl;
        if(quality_r3!="") fqual2 << title2 << endl << line4 << endl;
        }
       }
       }
 freads.close(); infile1.close();  
 if(csfasta_r3!="") {freads2.close(); infile3.close();}
 if(quality_file!="none") {fqual.close(); infile2.close();}
 if(quality_r3!="") {fqual2.close(); infile4.close();}
}
//---------------------------------------------------------------------------------
double get_number_of_bases(string csfasta_file, int & read_length)
{
   int count = 0;
  
  ifstream infile;
	infile.open(csfasta_file.c_str());
	string line;
		
    if (!infile) {
        cerr << "Unable to read " <<  csfasta_file << endl;
        exit(1);
    }
 
 unsigned long begin = infile.tellg();
  
 string sread = "", title = "";
 unsigned long cumseq = 0;
 unsigned long cumread = 0;
 
 while(getline(infile, line) && count < 10)
	{
	 if(line.length()>0)
	 {
	  if(line[0]=='>')
	  {
	  if(title.length()>1) 
	  {
     cumseq +=title.length();
     cumread += sread.length();
	   count++;
	  }
	  sread="";
	  title=line;
	  } 
     else
 	   if(line[0]!='#')sread+=line;
	}
 }
 
   infile.seekg (0, ios::end);
 unsigned long end = infile.tellg();
 end-=begin;
 
 infile.close();
 cumseq+=cumread;
 
 read_length = cumread/count-1;
 
 
 if(cumread<20) return 0;
 else
 return 1.0*(end)/(((double)cumseq)/count)*((double)read_length);
 
 
}
//-------------------------------------------------------------------------


