/*
 * SOLiD Sequence Assembilng Project
 * Copyright 2009 Life Technologies
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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>

/*  Extracts a number of good (no ".") reads from 4 files
 *  Notice that there no correspondence between F3/R3
 *  
 */

using namespace std;

struct TagPosition
{
 TagPosition(string t, int p){ tag = t; pos = p;}
 string tag;
 int pos;
};

//---------------------------------------------------------------------------------
void ExtractSubsetOfReads(string FileNameFW, string FileNameRW, string FileNameQVFW, string FileNameQVRW,
                           string outFW, string outRW, string outQVFW, string outQVRW, long unsigned number);

//---------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
 ExtractSubsetOfReads(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], atoi(argv[9]));
	return 0;
}
//---------------------------------------------------------------------------------
void ExtractSubsetOfReads(string FileNameFW, string FileNameRW, string FileNameQVFW, string FileNameQVRW,
                           string outFW, string outRW, string outQVFW, string outQVRW, long unsigned number)
{
	ifstream infile1, infile2, infile3, infile4;
	infile1.open(FileNameFW.c_str());
	infile2.open(FileNameRW.c_str());
	infile3.open(FileNameQVFW.c_str());
	infile4.open(FileNameQVRW.c_str());
	
	string line1,line2,line3,line4;
    if (!infile1||!infile2||!infile3||!infile4) {
        cerr << "Problem with one if the input files" << endl;
        exit(1);
    }

  ofstream ostr1, ostr2, ostr3, ostr4;
  ostr1.open(outFW.c_str(),ofstream::app);
  ostr2.open(outRW.c_str(),ofstream::app);
  ostr3.open(outQVFW.c_str(),ofstream::app);
  ostr4.open(outQVRW.c_str(),ofstream::app);
  
  string title1 = "", title2 ="";
  long unsigned count = 0;
  int pos1 = 0, pos2 = 0;
  int a1,a2,a3,b1=0,b2,b3;
  string::size_type pos;
  
  while(count<number && getline(infile1, line1) && getline(infile3, line3))
      {
       if(line1[0]=='#') continue;
       if(line1[0]=='>')
       {
        if(line1!=line3) { cout << "Error in files synchronization" << endl; exit(0); }
        title1 = line1.substr(1,line1.length()-3);
        continue;
       }
       
       if((pos=line1.find("."))==string::npos || (pos=line1.find(".",pos+1))==string::npos  || (pos=line1.find(".",pos+1))==string::npos)
       {
        pos1 = title1.find("_",1);
        pos2 = title1.find("_",pos1+1);
        a1 = atoi(title1.c_str());
        a2 = atoi(title1.substr(pos1+1).c_str());
        a3 = atoi(title1.substr(pos2+1).c_str());
        
        if(b1>0)
        {
        pos1 = title2.find("_",1);
        pos2 = title2.find("_",pos1+1);
        b1 = atoi(title2.c_str());
        b2 = atoi(title2.substr(pos1+1).c_str());
        b3 = atoi(title2.substr(pos2+1).c_str());
        }
        
        while((b1<a1||(b1==a1&&b2<a2)||(b1==a1&&b2==a2&&b3<a3))&&getline(infile2, line2) && getline(infile4, line4))
        {
         if(line2[0]=='#') continue;
          if(line2[0]=='>')
          {
           if(line2!=line4) { cout << "Error in files synchronization" << endl; exit(0); }
           title2 = line2.substr(1,line2.length()-6);
           //cout << title1 << " " << title2 << endl;
          
        pos1 = title2.find("_",1);
        pos2 = title2.find("_",pos1+1);
        b1 = atoi(title2.c_str());
        b2 = atoi(title2.substr(pos1+1).c_str());
        b3 = atoi(title2.substr(pos2+1).c_str());
        continue;
        }
        
        }
        if(title2==title1  && getline(infile2, line2) && getline(infile4, line4) && rand()%20==0)
        {       
        {
        count++;
        ostr1 << ">" << title1 << "F3" << endl << line1 << endl;
       	ostr2 << ">" << title2 << "F5-BC" << endl << line2 << endl;
       	ostr3 << ">" << title1 << "F3"  << endl << line3 << endl;
       	ostr4 << ">" << title2 << "F5-BC"  << endl << line4 << endl;
       	}
       	
       	if(count%(number/100/20)==0) cout << 100.0*count/number/20 << endl; 
       }
     }
     }
 ostr1.close();    
 ostr2.close();
 ostr3.close();
 ostr4.close();
 infile1.close();
 infile2.close();
 infile3.close();
 infile4.close();     
}
//---------------------------------------------------------------------------------

