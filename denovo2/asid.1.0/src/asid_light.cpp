/*
 * SOLiD Reads De Novo Assembly Project
 * Copyright (c) 2010-2011 Life Technologies
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
Assembly Assistant for SOLiD(TM) (ASiD) v.0.2                               \n\
Tool for filling gaps between contigs in scafolds                            \n\
and base-space conversion of final assembly                                  \n\
Copyright (2010-2011) by Life Technologies                                   \n\
***************************************************                           "
#define HELPMESSAGE "                                                        \n\
Usage:                                                                       \n\
                                                                             \n\
 - for collecting reads from gap regions run:                                \n\
                                                                             \n\
 ./asid_light -collect <global_graph> <de_scaffolds> <de_reads> <ins_length> \n\
              <tmp_dir>                                                      \n\
                                                                             \n\
 - for joining global contigs using local contigs and (i) generating final   \n\
   double encoded sequence or (ii) generate ma compatible file containing    \n\
   reads from global and local contigs run:                                  \n\
                                                                             \n\
 ./asid_light -merge <global_graph> <de_scaffolds> <color_reads> <reads_idx> \n\
  <tmp_dir> <reads_ma> <asid_scaff_de> <graph2ma | fixed2de | fixed2ma>      \n\
  <min_cnt_len> <lib_type>                                                   \n\
                                                                             \n\
 - for converting assembly (alignment of reads in ma format) into base-space:\n\
                                                                             \n\
 ./asid_light -convert <reads_ma> <line_size> > <nt_contigs>                 \n\
                                                                             \n\
- for joining global contigs using local contigs and generating final        \n\
  base-space sequence                                                        \n\
                                                                             \n\
 ./asid_light -combine <nt_contigs> <join_nt_contigs> <join_nt_scaffolds>    \n\
  <line_size> <min_cnt_length>                                               \n\
                                                                             \n\
Input:                                                                       \n\
                                                                             \n\
 global_graph  - LastGraph file produced by Velvet.                           \n\
 de_scaffolds  - de sequence of scaffolds produced by Velvet (contigs.fa).    \n\
 de_reads      - ordered list of reads created by Velvet (Sequences)          \n\
 ins_length    - library insert size for paired-end or mate-pair data         \n\
 tmp_dir       - location of temporary files containing de reads from gap     \n\
                 regions and localy assembled graphs.                         \n\
 color_reads   - ordered list of reads in color space corresponding to de list\n\
                 used by Velvet (Sequencies).                                 \n\
 reads_ma      - (Output/Input) alignment of reads in ma like format.         \n\
 min_cnt_len   - minimum contig length produced by assembler.                 \n\
 lib_type      - library type: fragment | mates | paired.                     \n\
 line_size     - characters per line in the output file.                      \n\
 graph2ma      - convert global graph into aligned reads in ma format.        \n\
 fixed2ma      - convert global graph & local graphs into one ma file.        \n\
 fixed2de      - merge sequences of global and local graphs and output de contigs.\n\
                                                                              \n\
Output:                                                                       \n\
                                                                              \n\
 asid_scaff_de - de sequence of scaffolds treated with ASiD.                  \n\
 reads_idx     - name for internaly created index of reads file.              \n\
 nt_contigs    - contigs in base-space.                                       \n\
 join_nt_contigs,                                                             \n\
 join_nt_scaffolds - sequence of joined contigs and scaffolds in base-space   \n\
                                                                              \n\
 Example:                                                                     \n\
                                                                              \n\
 ./asid_light -collect velvet/LastGraph velvet/contigs.fa velvet/Sequences postprocessor/gap_reads\n\
                                                                              \n\
 ./asid_light -merge velvet/LastGraph velvet/contigs.fa preprocessor/colorspace_input.csfasta\n\
   postprocessor/colorspace_input.idx postprocessor/gap_reads/ postprocessor/color_reads.ma\n\
   asid_scaffolds.de fixed2ma 100                                              \n\
                                                                               \n\
 ./asid_light -convert postprocessor/color_reads.ma 70 > postprocessor/asid_ntcontigs.tmp\n\
                                                                               \n\
 ./asid_light -combine postprocessor/asid_ntcontigs.tmp nt_contigs.fa nt_scaffolds.fa 70 100\n\
                                                                             \n"


#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "util.h"
#include "zutil.h"


using namespace std;
//---------------------------------------------------------------------------------
struct TRead
{
 TRead(long unsigned i, long int o, int s){id = i; offset = o; startcoord = s;}
 long unsigned id;
 long int offset;
 int startcoord;
   bool operator<(const TRead& data) const
	 {
	  return (offset < data.offset);
	 }
};

struct TInsert
{
 int start;
 int end;
 int shift;
 int insert_start;
 int insert_end;
 TInsert(int s, int e, int sh, int is, int ie){ start = s; end = e; shift = sh; insert_start = is; insert_end = ie; }
};

struct TGap
{
 int start,end;
 TGap(int s, int e){ start = s; end = e; };
 vector<string> frag_set;
 
};


struct TNode
{
 string fw_seq,rv_seq;
 vector<TRead> fw_reads, rv_reads;
 vector<TGap> gaps;
};

struct TArc
{
 int start_node, end_node, multiplicity;
};

struct TGapReads
{
  int read_id;
  int gap_id;
  int node_id;
  char status;
  
  // 0 - fw read on left anchor
  // 1 - fw read on right anchor
  // 2 - fw read which is mate of rw read from right anchor
  // 3 - rw read which is mate of fw read from left anchor
  // 4 - rw read on left anchor
  // 5 - rw read on right anchor
  
  TGapReads(int rid, int gid, int nid, char c ){read_id = rid; gap_id = gid; node_id = nid; status = c;}
  
  bool operator<(const TGapReads& data) const
	 {
	  return (read_id < data.read_id);
	 }
};

struct TGraph
{
 int kmer_sz;
 vector<TNode> nodes;
 vector<TGapReads> greads;
 int load_from_file(string filename, int fid = 0,int min_cnt_len = 25);
 void load_gaps_coordinates(string filename);
 int set_gaps_list(int node_id, string & seq);
 void extract_gap_reads(string filename);
 void load_reads_from_gaps(unsigned ins_size);
 void output_gaps_reads(string temp_folder);
 void fill_gaps_between_contigs(string temp_folder, string seqfilename, string indxfilename, string outmapfile, 
                                 string mode,int significant_contig_len = 100, char lib_type = 'f');
 void output_graph_contigs(string temp_folder, int min_len);
};

//---------------------------------------------------------------------------------
void CombineContigsIntoScaffolds(string FileName, string out_contigs, string out_scaffolds, int minLength = 70, int significant_contig_len = 100);
void convert_aligned_reads_into_ntcontigs(string readsfile, int outSeqLen, string base_reads = "", int colorCost = 1, int BaseCost = 1);
bool covers_gap(string & contig, string & leftSeq, string & rightSeq, int & left_match, 
                 string::size_type & leftpos, int & right_match, string::size_type & rightpos);
void add_reads_to_gaps(TGap & gp, int read_id, char read_status, string & seq);
void generate_indexfile(string seqfilename, string indxfilename);
void ReadGenome(string &seq, string filename);
inline string reverse(string s);
string itos(int i);	// convert int to string

bool LOG = true;

//---------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  
  if(argc<2) { cout << HEADMESSAGE << endl; cout << HELPMESSAGE << endl; exit(1);};
  if(strcmp(argv[1],"-convert")==0) { } else {cout << HEADMESSAGE << endl;};
   

  if(strcmp(argv[1],"-collect")==0)
   {
    TGraph g;
    unsigned ins_size = atoi(argv[5]);
    g.load_from_file(argv[2]);
    g.load_gaps_coordinates(argv[3]);
    cout << "Loading reads/gaps attribute." << endl;
    g.load_reads_from_gaps(ins_size);
    cout << "Extract reads sequence." << endl;
    g.extract_gap_reads(argv[4]);
    cout << "Output gaps sequence in right files." << endl;  
    g.output_gaps_reads(argv[6]);
    return 0;
   }   
   
 // system("$denovo/velvet_0.7.55/");
 
 if(strcmp(argv[1],"-merge")==0)
 {
  TGraph g;
  int min_cnt_len = atoi(argv[10]);
  cout << "minimum contig length = " << min_cnt_len << endl;
  g.load_from_file(argv[2], 0, min_cnt_len);
  g.load_gaps_coordinates(argv[3]);
  cout << "Building index for sequence file .." << endl;
  generate_indexfile(argv[4],argv[5]);
  char lib_type = 'f';
  if(strcmp(argv[11],"paired")==0) lib_type = 'p';
  else
  if(strcmp(argv[11],"mates")==0) lib_type = 'm';
  g.fill_gaps_between_contigs(argv[6],argv[4],argv[5],argv[7],argv[9], min_cnt_len,lib_type);
  if(strcmp(argv[9],"fixed2de")==0)g.output_graph_contigs(argv[8],atoi(argv[10]));
  return 0;
 }
 
 if(strcmp(argv[1],"-combine")==0)
 {
  CombineContigsIntoScaffolds(argv[2], argv[3], argv[4] ,atoi(argv[5]), atoi(argv[6]));
  return 0;
 }

 if(strcmp(argv[1],"-convert")==0)
 {
  convert_aligned_reads_into_ntcontigs(argv[2], atoi(argv[3]));
} 
	return 0;
}
//---------------------------------------------------------------------------------
void de2color(string & sequence)
{
 for(int i=0;i<sequence.length();i++)
 {
   switch(toupper(sequence[i]))
   {
    case 'A': sequence[i] = '0'; break;
    case 'C': sequence[i] = '1'; break;
    case 'G': sequence[i] = '2'; break;
    case 'T': sequence[i] = '3'; break;
   }
 }
}
//---------------------------------------------------------------------------------
void color2de(string & sequence)
{
 for(int i=0;i<sequence.length();i++)
 {
   switch(sequence[i])
   {
    case '0': sequence[i] = 'A'; break;
    case '1': sequence[i] = 'C'; break;
    case '2': sequence[i] = 'G'; break;
    case '3': sequence[i] = 'T'; break;
   }
 }
}
//---------------------------------------------------------------------------------
void load_reads_id(string filename, vector<long int> & reads_id)
{

 	ifstream infile;
	infile.open(filename.c_str());

	string line;
    if (!infile) {
        cerr << "Can't read input file: " << filename << endl;
        exit(1);
    }
  
  while(getline(infile, line))
  {
  if(line.length()>1 && line[0]=='>')
  { 
   if(line[1]=='F') reads_id.push_back(atoi(&line[0]+2));
   else reads_id.push_back(-(atoi(&line[0]+2)));
  // cout << line << ">" << atoi(&line[0]+2) << endl;
  }
  }
  infile.close();
}
long unsigned MAX_READ_ID = 0;
//---------------------------------------------------------------------------------
void generate_indexfile(string seqfilename, string indxfilename)
{ 
  ifstream seqfile;
	seqfile.open(seqfilename.c_str());

	string line;
    if (!seqfile) {
        cerr << "Can't read sequence input file: " << seqfilename << endl;
        exit(1);
    }
 
	FILE * wFile;
  wFile = fopen ( indxfilename.c_str() , "wb" );
  if (wFile==NULL)
	  	  {
	  	   cout << "ERROR: Can not write index file. \n" << indxfilename << endl;
	  	   exit(1);
	  	  }
 
  streampos p = seqfile.tellg();
  while(getline(seqfile, line))
  {
   if(line.length()>1 && line[0]=='>')  {fwrite(&p,sizeof(streampos),1,wFile);MAX_READ_ID++;}
   p = seqfile.tellg();
  }

  fclose(wFile);
  seqfile.close();
  
}
//---------------------------------------------------------------------------------
inline void output_read_alignment(unsigned id, long int pzF3, long int pzR3, ifstream & idxFile, ifstream & seqfile, ofstream & ostr, char lib_type = 'f', bool addlen = false)
{
    string line,line2;
    streampos p,q;    
    q = (streampos)id*sizeof(streampos);
    if(id>MAX_READ_ID || id<0) {cout << "Internal ERROR: indexing requested id > than MAX_READ_ID=" << MAX_READ_ID << endl; exit(1);}
    idxFile.seekg(q);
	  idxFile.read((char *)(&p),sizeof(streampos));  	  
	  seqfile.seekg(p);
	  getline(seqfile,line);getline(seqfile,line2);
	  if(line.find("F3",0)!=string::npos)
	   {
	    if( (lib_type == 'm' && !addlen) || (lib_type == 'p' && addlen) ) pzF3 -= (long(line2.length()-2));
	    if(pzF3>0 || pzF3 <= -(long(line2.length()-2))) ostr << line << "," << pzF3 << endl << line2 << endl;
	   }
	  else
	   {
	   if(addlen)pzR3 -= (long(line2.length()-2));
	    if(pzR3>0 || pzR3 <= -(long(line2.length()-2)))  ostr << line << "," << pzR3 << endl << line2 << endl;
	   }
}
//---------------------------------------------------------------------------------
void TGraph::fill_gaps_between_contigs(string temp_folder, string seqfilename, string indxfilename, string outmapfile, string mode, int significant_contig_len, char lib_type)
{

  // mode
  // graph2ma (default)
  // fixed2ma
  // fixed2de
  
  bool fix_gaps = true;
  int total_gaps = 0;
  int number_of_fixed_gaps = 0;
  int number_of_switches = 0;
  
  //opening seq file which we will travers to extract reads
  
  ifstream seqfile;
	seqfile.open(seqfilename.c_str());

    if (!seqfile) {
        cerr << "ERROR: Can't read sequence input file: " << seqfilename << endl;
        exit(1);
    }
   
  ifstream idxFile; 
  idxFile.open( indxfilename.c_str() , ios::in | ios::binary);  
    if (!idxFile)
	  	  {
	  	   cerr << "ERROR: Can't open index file. \n" << indxfilename << endl;
	  	   exit(1);
	  	  }
	
	ofstream ostr;
	if(mode != "fixed2de")
	{
	  ostr.open(outmapfile.c_str());  
	  
	  if (!ostr) {
        cerr << "ERROR: Can't write output file: " << outmapfile << endl;
        exit(1);
    }
   }
   
  for(int node_id = 1; node_id < nodes.size(); node_id++)
    {
    if(LOG && node_id%1000==0) cout << " processed " << node_id << " nodes." << endl;
    if(nodes[node_id].fw_seq.length()+kmer_sz-1>=significant_contig_len)
     {
       vector<TNode> gap_graph;
       vector<int> gp_kmer_sz;
       int insert_shift = 0;

      for(int gap_id = 0; gap_id < nodes[node_id].gaps.size(); gap_id++) 
      {
       total_gaps++;
   
       if(mode == "fixed2ma" || mode == "fixed2de")
       {
	       bool fixed = false;
	       int file_id = 1;
         
         string::size_type leftpos, rightpos;
         int left_match, right_match;
         
         int left_start = nodes[node_id].gaps[gap_id].start - significant_contig_len + insert_shift;
         string leftSeq = nodes[node_id].fw_seq.substr(left_start<0?0:left_start,left_start<0?(significant_contig_len+left_start):significant_contig_len);
         int right_start = nodes[node_id].gaps[gap_id].end + insert_shift;
         string rightSeq = nodes[node_id].fw_seq.substr(right_start, significant_contig_len);
         if(leftSeq.length()<significant_contig_len || rightSeq.length()<significant_contig_len) fixed = false;
         else
	       while(!fixed)
	       {

	       //read mini-assembly results for this gap
	       string graph_file = temp_folder + "/reads_" + itos(node_id) + "_" + itos(gap_id) + ".de." + itos(file_id) + ".graph";
	  
	       TGraph gp;
	       bool lg = LOG; LOG = false;
	       if(gp.load_from_file(graph_file,file_id)!=0) { LOG = lg; break;}
         LOG = lg;
         
         // iterating through mini-assembled contigs
          for(int gp_node_id = 1; gp_node_id < gp.nodes.size(); gp_node_id++)
           if(gp.nodes[gp_node_id].fw_seq.length()>30)
           {
            if(covers_gap(gp.nodes[gp_node_id].fw_seq, leftSeq, rightSeq, left_match, leftpos, right_match, rightpos))
            {
             gp_kmer_sz.push_back(gp.kmer_sz);
             gap_graph.push_back(gp.nodes[gp_node_id]);
             
            
              if(mode == "fixed2de")
               {
                nodes[node_id].fw_seq = nodes[node_id].fw_seq.substr(0, left_start<0?left_match:(left_start+left_match)) +
                                 gp.nodes[gp_node_id].fw_seq.substr(leftpos,rightpos-leftpos) +
                                 nodes[node_id].fw_seq.substr(right_start + right_match);
               insert_shift += left_start<0?(rightpos-leftpos - (right_start + right_match-left_match)):(rightpos-leftpos - (right_start + right_match-left_start-left_match));
               }
               
             fixed = true;   
             number_of_fixed_gaps++;
             break;
            } 
           }
           file_id+=1;
          }
      
         if(!fixed)
          {
           gp_kmer_sz.push_back(0);
           gap_graph.push_back(TNode()); 
           
           if(mode == "fixed2de")
            {
             nodes[node_id].fw_seq = nodes[node_id].fw_seq.substr(0, left_start+significant_contig_len) +
                                 "NN" +
                                 nodes[node_id].fw_seq.substr(right_start);
             insert_shift += 2 - (right_start-left_start-significant_contig_len);
            }        
          }
        }
       }
    
        if(mode == "fixed2de") continue;
           
        //output aligned reads in corresponding contigs
        
        sort(nodes[node_id].fw_reads.begin(), nodes[node_id].fw_reads.end());
        sort(nodes[node_id].rv_reads.begin(), nodes[node_id].rv_reads.end());
        
        int contig_len = nodes[node_id].fw_seq.length();
        nodes[node_id].gaps.push_back(TGap(contig_len,contig_len));
        gp_kmer_sz.push_back(0);
        
        int gap_id = 0;
        int contig_start = 0;
        ostr << "-end" << endl << "=NODE_" << node_id << "-" << gap_id << "." << (nodes[node_id].gaps[gap_id].end - contig_start) << endl;
        
        int rv_rid = nodes[node_id].rv_reads.size() - 1;
        int fw_rid = 0;
        unsigned shift = 50; // shifts all graph alignmnets by 100, designed to include reads with start pos out of reference
        
        while(fw_rid < nodes[node_id].fw_reads.size() || rv_rid >=0)
        {
         while(gap_id<nodes[node_id].gaps.size() && 
         (fw_rid >=nodes[node_id].fw_reads.size() || nodes[node_id].fw_reads[fw_rid].offset > nodes[node_id].gaps[gap_id].start) 
         && (rv_rid < 0 || contig_len - nodes[node_id].rv_reads[rv_rid].offset > nodes[node_id].gaps[gap_id].start)) 
         {
         
         if(mode == "fixed2ma" && gp_kmer_sz[gap_id]>0)
         {
          
          // output additional reads from gapped region        
         
           string reads_file = temp_folder + "/reads_" + itos(node_id) + "_" + itos(gap_id) + ".de";
                   
           vector<long int> reads_id;
           reads_id.push_back(0);
           load_reads_id(reads_file,reads_id);
           
           //this could be removed
           //sort(gap_graph[gap_id].fw_reads.begin(), gap_graph[gap_id].fw_reads.end());
           //sort(gap_graph[gap_id].rv_reads.begin(), gap_graph[gap_id].rv_reads.end());
           
           long unsigned gap_contig_len = gap_graph[gap_id].fw_seq.length();          
           ostr << "-end" << endl << "=NODE_" << node_id << "-" << gap_id << "+." << gap_contig_len << endl;
           
           for(int idx=0;idx<gap_graph[gap_id].fw_reads.size(); idx++)
            {
               if(reads_id[gap_graph[gap_id].fw_reads[idx].id]>0)
                {
                  //regular reads
                  long int offset =  gap_graph[gap_id].fw_reads[idx].offset;
                  int startcoord = gap_graph[gap_id].fw_reads[idx].startcoord;
                  if(offset + startcoord >= 0) output_read_alignment(reads_id[gap_graph[gap_id].fw_reads[idx].id]-1,
                                 (lib_type=='p')?(offset - startcoord):(-offset + startcoord), offset - startcoord, idxFile,seqfile, ostr, lib_type);          
                  }
                   else
                  { 
                  //flipped reads
                  long int offset = gap_graph[gap_id].fw_reads[idx].offset;
                  int startcoord = gap_graph[gap_id].fw_reads[idx].startcoord;
                  if(offset + startcoord >= 0) output_read_alignment(-reads_id[gap_graph[gap_id].fw_reads[idx].id]-1,
                             (lib_type=='p')?(-offset + startcoord):offset, -offset + startcoord , idxFile,seqfile, ostr, lib_type, 1);
                  }
                 }
          }
         
         contig_start = nodes[node_id].gaps[gap_id].end;
         gap_id++;
         if(gap_id<nodes[node_id].gaps.size()) ostr << "-end" << endl << "=NODE_" << node_id << "-" << gap_id << "." << (nodes[node_id].gaps[gap_id].end - nodes[node_id].gaps[gap_id-1].start) << endl;
         }

         if(gap_id>=nodes[node_id].gaps.size()) break;
       
         if(fw_rid < nodes[node_id].fw_reads.size() && nodes[node_id].fw_reads[fw_rid].offset<=nodes[node_id].gaps[gap_id].start) 
          {
          long int pz = nodes[node_id].fw_reads[fw_rid].offset - nodes[node_id].fw_reads[fw_rid].startcoord - contig_start + (gap_id>0?(kmer_sz-1):0);
          //cout << "place 1 " << pz << " " << contig_len << " " << nodes[node_id].rv_reads[rv_rid].offset << " " <<  nodes[node_id].rv_reads[rv_rid].startcoord << " " << contig_start << " " << (gap_id>0?(kmer_sz-1):0) << endl;
          output_read_alignment(nodes[node_id].fw_reads[fw_rid].id-1,(lib_type!='m')?(shift + pz):-pz,(lib_type!='m')?(shift + pz):pz, idxFile,seqfile, ostr, lib_type);          
          fw_rid++;
          }
          
          if(rv_rid >=0 && contig_len - nodes[node_id].rv_reads[rv_rid].offset<=nodes[node_id].gaps[gap_id].start) 
          {
          long int pz = contig_len-nodes[node_id].rv_reads[rv_rid].offset + nodes[node_id].rv_reads[rv_rid].startcoord - contig_start + (gap_id>0?(kmer_sz-1):0);
          //cout << "place 4 " << pz << " " << contig_len << " " << nodes[node_id].rv_reads[rv_rid].offset << " " <<  nodes[node_id].rv_reads[rv_rid].startcoord << " " << contig_start << " " << (gap_id>0?(kmer_sz-1):0) << endl;
          output_read_alignment(nodes[node_id].rv_reads[rv_rid].id-1,(lib_type!='m')?-(shift + pz):pz, (lib_type!='m')?-(shift + pz):-pz, idxFile,seqfile, ostr, lib_type);
          rv_rid--;          
          } 
        
        //char a;
        //cin >> a;
        
        }
  }
 }
 
 seqfile.close();
 idxFile.close(); 
 ostr.close(); 

 if(LOG) cout << "Total gaps: " << total_gaps << endl;
 if(LOG) cout << "Number of fixed gaps: " << number_of_fixed_gaps << endl;
 
}
//---------------------------------------------------------------------------------
bool covers_gap(string & contig, string & leftSeq, string & rightSeq, int & left_match, string::size_type & leftpos, int & right_match, string::size_type & rightpos)
{
 for(left_match=0;left_match<(long(leftSeq.length()))-25; left_match+=10)
  {
   if((leftpos=contig.find(leftSeq.substr(left_match,25),0))!=string::npos)
    {
     for(right_match=0;right_match<(long(rightSeq.length()))-25; right_match+=10)
      {
       if((rightpos=contig.find(rightSeq.substr(right_match,25),0))!=string::npos) break;
      }
    }
   if(leftpos!=string::npos && rightpos!=string::npos) break;
  }
 if(leftpos!=string::npos && rightpos!=string::npos) return true;
 return false;
} 
//---------------------------------------------------------------------------------
void TGraph::output_gaps_reads(string temp_folder)
{
  DIR *dir = opendir(temp_folder.c_str());
	if (dir == NULL) mkdir(temp_folder.c_str(), 0777);
	
  for(int node_id = 1; node_id < nodes.size(); node_id++)
  {
   if(node_id%1000==0) cout << " processed " << node_id << " nodes." << endl;
   for(int gap_id = 0; gap_id < nodes[node_id].gaps.size(); gap_id++) 
   {
	
	  string outFile = temp_folder + "/reads_" + itos(node_id) + "_" + itos(gap_id) + ".de";
		ofstream ostr;
	  ostr.open(outFile.c_str());  
	  if(ostr==NULL)
	  	  {
	  	   cout << "ERROR: Can not write output file. \n" << outFile << endl;
	  	   exit(1);
	  	  }
	  for(int i=0;i<nodes[node_id].gaps[gap_id].frag_set.size();i++) ostr << nodes[node_id].gaps[gap_id].frag_set[i] << endl;
	  ostr.close();
	 }
	} 
}	
//---------------------------------------------------------------------------------
void TGraph::load_reads_from_gaps(unsigned ins_size)
{
  int anchor_size = 100;
  int inclusion_size = unsigned(1.5*ins_size);
  if(inclusion_size<anchor_size)inclusion_size=anchor_size;
  int left_gap_end, left_gap_start, left_gap_start2, right_gap_start, right_gap_end, right_gap_end2;
  
  for(int node_id = 1; node_id < nodes.size(); node_id++)
  {
   if(node_id%1000==0) cout << " processed " << node_id << " nodes." << endl;
   int contig_len = nodes[node_id].fw_seq.length();
   for(int gap_id = 0; gap_id < nodes[node_id].gaps.size(); gap_id++) 
   {
   
    left_gap_end    = nodes[node_id].gaps[gap_id].start;  
    left_gap_start  = left_gap_end - anchor_size;
    left_gap_start2 = left_gap_end - inclusion_size;
    right_gap_start = nodes[node_id].gaps[gap_id].end;
    right_gap_end   = right_gap_start + anchor_size;
    right_gap_end2  = right_gap_start + inclusion_size;
    
    // -- analysizng reads mapped to fw sequence
    
    for(int idx = 0; idx < nodes[node_id].fw_reads.size(); idx++)
    {
    if(nodes[node_id].fw_reads[idx].offset <= left_gap_end)
     {
      if(nodes[node_id].fw_reads[idx].offset >= left_gap_start)
       greads.push_back(TGapReads(nodes[node_id].fw_reads[idx].id, gap_id, node_id, 0)); // 0 - fw read on left anchor
    
      if(nodes[node_id].fw_reads[idx].offset >= left_gap_start2)
       {
        if(nodes[node_id].fw_reads[idx].id % 2 == 1) // F3 read
           greads.push_back(TGapReads(nodes[node_id].fw_reads[idx].id + 1, gap_id, node_id, 3)); // 3 - rw read which is mate of fw read from left anchor 
           else  // R3 read
           greads.push_back(TGapReads(nodes[node_id].fw_reads[idx].id - 1, gap_id, node_id, 3)); // 3 - rw read which is mate of fw read from left anchor 
       }    
     }  
       
    if(nodes[node_id].fw_reads[idx].offset >= right_gap_start && nodes[node_id].fw_reads[idx].offset <= right_gap_end) 
       greads.push_back(TGapReads(nodes[node_id].fw_reads[idx].id, gap_id, node_id, 1)); // 1 - fw read on right anchor
     
    }
    
    // -- analyzing reads mapped to rw sequence
    
    
    for(int idx = 0; idx < nodes[node_id].rv_reads.size(); idx++)
    {
    if(nodes[node_id].rv_reads[idx].offset >= contig_len - left_gap_end && nodes[node_id].rv_reads[idx].offset <= contig_len - left_gap_start)
       greads.push_back(TGapReads(nodes[node_id].rv_reads[idx].id, gap_id, node_id, 4));  // 4 - rw read on left anchor

    if(nodes[node_id].rv_reads[idx].offset <= contig_len - right_gap_start)
      {
        if(nodes[node_id].rv_reads[idx].offset >= contig_len - right_gap_end)
        greads.push_back(TGapReads(nodes[node_id].rv_reads[idx].id, gap_id, node_id, 5)); // 5 - rw read on right anchor
        
        if(nodes[node_id].rv_reads[idx].offset >= contig_len - right_gap_end2)
        {
         if(nodes[node_id].rv_reads[idx].id % 2 == 1) // F3 read
         greads.push_back(TGapReads(nodes[node_id].rv_reads[idx].id+1, gap_id, node_id, 2)); // 2 - fw read which is mate of rw read from right anchor
         else // R3 read
         greads.push_back(TGapReads(nodes[node_id].rv_reads[idx].id-1, gap_id, node_id, 2)); // 2 - fw read which is mate of rw read from right anchor
        }
        
      }
    }
    // ----
          
   }
  }
  sort(greads.begin(), greads.end());
}  
//---------------------------------------------------------------------------------
void TGraph::extract_gap_reads(string filename)
{

 	ifstream infile;
	infile.open(filename.c_str());

	string line;
    if (!infile) {
        cerr << "Can't read input file: " << filename << endl;
        exit(1);
    }
  
  string::size_type pos;
  string seq = "";
  string title = "";
  int read_id = 0,idx = 0;
  cout << "number of reads to be collected: " << greads.size() << endl;
  
  while(getline(infile, line) && idx < greads.size())
  {
    if(line.length()>0)
     { 
      if(line[0]=='>')
      {
        if(title!="")
        {

         if((pos = title.find("\t",0))!=string::npos) read_id = atoi(&title[0]+pos+1);
          else
          { 
           cerr << "ERROR: sequence id is missing in " << title << endl;
           exit(1);
          }

          while(read_id == greads[idx].read_id)
          {
           add_reads_to_gaps(nodes[greads[idx].node_id].gaps[greads[idx].gap_id], read_id, greads[idx].status, seq);
           idx++;
          }
          
        }
        title = line; seq = ""; }
      else
      seq += line;
     }
  }         
 
 infile.close();
 // repeated the same as above for the last read
         
        if(title!="")
        {

         if((pos = title.find("\t",0))!=string::npos) read_id = atoi(&title[0]+pos+1);
          else
          { 
           cerr << "ERROR: sequence id is missing in " << title << endl;
           exit(1);
          }
          
          while(read_id == greads[idx].read_id)
          {
           add_reads_to_gaps(nodes[greads[idx].node_id].gaps[greads[idx].gap_id], read_id, greads[idx].status, seq);
           idx++;
          }
        }

   cout << "number of recorded reads: " << idx << endl;
}
//---------------------------------------------------------------------------------
void add_reads_to_gaps(TGap & gp, int read_id, char read_status, string & seq)
{
   // 0 - fw read on left anchor
   // 1 - fw read on right anchor
   // 2 - fw read which is mate of rw read from right anchor

   if(read_status<3) {gp.frag_set.push_back(">F"+itos(read_id)); gp.frag_set.push_back(seq);}              
   else
   {gp.frag_set.push_back(">R"+itos(read_id)); gp.frag_set.push_back(reverse(seq));}   
 
   // 3 - rw read which is mate of fw read from left anchor
   // 4 - rw read on left anchor
   // 5 - rw read on right anchor
}
//---------------------------------------------------------------------------------

int TGraph::load_from_file(string filename, int fid, int min_cnt_len)
{
 	ifstream infile;
	infile.open(filename.c_str());

	string line;
    if (!infile) {
        if(fid<=1) 
        {
        cerr << "Can't read input file: " << filename << endl;
        }
        if(fid==0) exit(1);
        return 1;
    }

  int num_nodes = 0;

  string::size_type pos,loc;
  if(getline(infile, line))
  {
    if(line.length()>0){
     num_nodes = atoi(line.c_str());
     kmer_sz = atoi(line.substr(line.find("\t",line.find("\t",0)+1)+1).c_str());
    } 
  }
  
  if(LOG) cout << fid << " loading " << num_nodes << " nodes of the graph " << endl;
  if(LOG) cout << "k-mer size " << kmer_sz << endl;  
  
  if(num_nodes>0) 
  {
   nodes.resize(num_nodes+1);
  } 

  int num_loaded_nodes = 0, node_id;
    
  while(getline(infile, line) && (line.length()<2 || line.substr(0,2)!="NR"))
  {
   if(line.length()>0 && line.substr(0,4)=="NODE")
   {
    node_id = atoi(&line[0]+5);
    getline(infile, nodes[node_id].fw_seq);
    getline(infile, nodes[node_id].rv_seq);
    if(++num_loaded_nodes % 10000 == 0) cout << num_loaded_nodes << endl;
   }
   
  }
  
  if(LOG) cout << "loaded " << num_loaded_nodes << " nodes." << endl;
  
  if(LOG) cout << "loading read ids into memory" << endl;

  unsigned long long num_loaded_reads = 0;

  if(line.length()>=2 && line.substr(0,2)=="NR") 
  {
   node_id = atoi(line.c_str()+3);

   while(getline(infile, line))
   {
    if(line.length()>0)
     {
      if(line.substr(0,2)=="NR") node_id = atoi(&line[0]+3);
      else
      {
       if((pos = line.find("\t",0))!=string::npos && (loc = line.find("\t",pos+1))!=string::npos)
        {
         if(node_id < 0) { 
         if(nodes[-node_id].fw_seq.length()+kmer_sz>=min_cnt_len) 
         {
         nodes[-node_id].rv_reads.push_back(TRead(atoi(line.c_str()), atoi(line.substr(pos+1).c_str()), atoi(line.substr(loc+1).c_str()) ));
         num_loaded_reads++;
         }
         }
         else
         { 
          if(nodes[node_id].fw_seq.length()+kmer_sz>=min_cnt_len)  
          {
          nodes[node_id].fw_reads.push_back(TRead(atoi(line.c_str()), atoi(line.substr(pos+1).c_str()), atoi(line.substr(loc+1).c_str()) ));
          num_loaded_reads++;
          }
          }
         
        } 
      }
     }
   }    
  } 

 if(LOG) cout << "loaded " << num_loaded_reads << " read ids." << endl;

infile.close();
return 0;
}

//---------------------------------------------------------------------------------
void TGraph::load_gaps_coordinates(string filename) {
	
 	ifstream infile;
	infile.open(filename.c_str());

	string line;
    if (!infile) {
        cerr << "Can't read input file: " << filename << endl;
        exit(1);
    }

  int node_id = -1;
  string scaf_seq = "";
  int num_gaps = 0;
  
  while(getline(infile, line))
  {
   if(line.length()>0)
   {
    if(line.substr(0,5)==">NODE")
    {
     if(scaf_seq!="" && node_id>0) num_gaps+=set_gaps_list(node_id,scaf_seq);
     node_id = atoi(&line[0]+6);   
     scaf_seq = "";
    } else scaf_seq+=line;
   }
  }
  if(scaf_seq!="" && node_id>0) num_gaps+=set_gaps_list(node_id,scaf_seq);  
  infile.close();
   
   if(LOG) cout << "Total number of gaps " << num_gaps << endl; 
 } 

//---------------------------------------------------------------------------------
int TGraph::set_gaps_list(int node_id, string & seq) {
    nodes[node_id].fw_seq = seq.substr(0,kmer_sz-1) + nodes[node_id].fw_seq;
    int num_gaps = 0;
	  string::size_type pos = 0;
	  while((pos=seq.find( "N", pos ))!=string::npos)
	  {
	   int start = pos;
	   int end = pos;
	   while(seq[++end]=='N');
	   pos = end;
	   nodes[node_id].gaps.push_back(TGap(start, end));
	   num_gaps++;
    }
    return num_gaps;
}  	
//---------------------------------------------------------------------------------
string itos(int i)	// convert int to string
{
		stringstream s;
		s << i;
		return s.str();
}
//---------------------------------------------------------------------------------
inline string reverse(string s){
	string x=s;
	for(int i=s.length()-1;i>=0;i--)
		{
		 x[s.length()-1-i] = s[i];
		}
	return x;
	}
//---------------------------------------------------------------------------------
void ReadGenome(string &seq, string filename)
{
	ifstream infile;
	infile.open(filename.c_str());
	string line;
    if (!infile) {
        cerr << "Unable to read " <<  filename  << endl;
        exit(1);
    }
    getline(infile, line);
    while(getline(infile, line)) seq+=line;
    infile.close();
} 
//---------------------------------------------------------------------------------
void TGraph::output_graph_contigs(string outFile, int min_len)
{
		ofstream ostr;
	  ostr.open(outFile.c_str());  
	  if(ostr==NULL)
	  	  {
	  	   cerr << "ERROR: Can not write output file. \n" << outFile << endl;
	  	   exit(1);
	  	  }
	  	  
	   for(int node_id = 1; node_id < nodes.size(); node_id++)
     { 	  
	  	if(nodes[node_id].fw_seq.length()>=min_len)
	  	{
	  	 ostr << ">"<<node_id << "_length_"<<nodes[node_id].fw_seq.length() << endl;
	  	 for(int i=0;i<nodes[node_id].fw_seq.length();i+=70)
	  	 {
	  	  ostr << nodes[node_id].fw_seq.substr(i,70) << endl;
	  	 }
	  	}
	  }
	  ostr.close();
}	
//---------------------------------------------------------------------------------
inline void outputs(string s, int block_sz, ofstream & ostr, string title)
{
 ostr << title << endl;
 int i=0;
 while(i<s.length()){ostr << s.substr(i,block_sz) << endl; i+=block_sz;}
}

//---------------------------------------------------------------------------------
inline void output(string & s, int block_sz, ofstream & ostrc, ofstream & ostrs, string & title)
{
 	string::size_type prev = 0, dotpos=0;
 	int idx = 0;
	while(prev<s.length() && (dotpos = s.find(".", prev))!=string::npos)
	{
	s[dotpos] = 'N';
	 outputs(s.substr(prev, dotpos-prev-1), block_sz, ostrc, title + "-" + itos(idx));
	 prev=dotpos+1;
	 idx++;
	}
	if(prev<s.length()) outputs(s.substr(prev), block_sz, ostrc, title + "-" + itos(idx));
	
	 outputs(s, block_sz, ostrs, title);
}
//---------------------------------------------------------------------------------
void load_seq(ifstream & infile, string & seq, string & next_title)
{
 next_title = "";
 string line;
 while(getline(infile, line))
	{
	 if(line.length()>0)
	 {
	  if(line[0]=='>')
	  {
	   next_title = line;
	   return;
	  }
	  else
 	  if(line[0]!='#')seq+=line;
	 }
 }
}
//---------------------------------------------------------------------------------
void CombineContigsIntoScaffolds(string FileName,  string out_contigs, string out_scaffolds, int minLength, int significant_contig_len)
{
  int total_scaffolds = 0;
  int total_gaps = 0;
  int fixed_gaps = 0;


	ifstream infile;
	infile.open(FileName.c_str());

	string line;
    if (!infile) {
        cerr << "Problem with input file " << FileName << endl;
        exit(1);
    }

  ofstream ostrc, ostrs;
	ostrc.open(out_contigs.c_str());  
	if(ostrc==NULL)
	 {
	  cerr << "ERROR: Can not write output file. \n" << out_contigs << endl;
	  exit(1);
	 }

	ostrs.open(out_scaffolds.c_str());  
	if(ostrs==NULL)
	 {
	  cerr << "ERROR: Can not write output file. \n" << out_scaffolds << endl;
	  exit(1);
	 }

 string title = "";
 while(getline(infile, line)&&(line.length()<1 || line[0]!='>') );
 if(line.length()>1 && line[0]=='>') title = line;
 string next_title = "", seq;
 do
 {
   load_seq(infile, seq, next_title);
   
   if(seq.length()>0) 
	   {
	    if(next_title == "" || atoi(title.substr(title.find("_",0)+1).c_str())!=atoi(next_title.substr(next_title.find("_",0)+1).c_str()))
	    {
	
	     total_scaffolds++;
      
	     output(seq,minLength, ostrc, ostrs, title);
	     seq = "";
	     title=next_title; 
	    }  
	    else
	    {
	       total_gaps++;
	       
	     while(next_title!="" && next_title.find("+",0)!=string::npos)
	       {
	        string gap_seq = "";
	        string gap_title = next_title;
	        string next_contig = "";
	        string next_contig_title = "";
	        
	        load_seq(infile, gap_seq, next_contig_title);
	        load_seq(infile, next_contig, next_title);
	        
	       int left_start = seq.length()-significant_contig_len;
         string leftSeq = seq.substr(left_start<0?0:left_start,left_start<0?(significant_contig_len+left_start):significant_contig_len);
         int right_start = 0;
         string rightSeq = next_contig.substr(right_start, significant_contig_len);
         
         string::size_type leftpos, rightpos;
         int left_match, right_match;
         
         if(covers_gap(gap_seq, leftSeq, rightSeq, left_match, leftpos, right_match, rightpos))
          {
           seq = seq.substr(0, left_start<0?left_match:(left_start+left_match)) +
                            gap_seq.substr(leftpos,rightpos-leftpos) +
                            next_contig.substr(right_start + right_match);
           fixed_gaps++;
          } 
           else
           seq+="."+next_contig;
           
           total_gaps++;
          }
          
       if(next_title == "") {output(seq,minLength,ostrc,ostrs,title); }
       else 
       if(atoi(title.substr(title.find("_",0)+1).c_str())==atoi(next_title.substr(next_title.find("_",0)+1).c_str())){
        total_gaps++;
        seq+="."; 
        }
       else
       {
       total_scaffolds++;
       output(seq,minLength,ostrc,ostrs,title);
	     seq = "";
	     title=next_title; 
	     }   
	    }
	 }
 }
 while(next_title!="");

 ostrs.close(); 
 ostrc.close(); 
 infile.close();
 
  if(LOG)
 {
  cout << "# scaffolds  : " << total_scaffolds << endl;
  cout << "# gaps       : " << total_gaps << endl;
  cout << "# fixed gaps : " << fixed_gaps << endl;
  }

}
//---------------------------------------------------------------------------------
// Zheng's code for base translation
//---------------------------------------------------------------------------------
class Color2Seq {
  public:
    Color2Seq(unsigned cpl){ChrPerLine = cpl;};
    ~Color2Seq();
    int ChrPerLine;
    int  inputReads(FILE *file, int cost);
    int  inputBaseReads(FILE *file, int cost);
    void outputSeq();
    void init(int alen);
  protected:
    int assemlen;
    int **Co, **Ba;
    char *seq, *cseq, **color, *notN;
    int **pa;
    char color2seq[128][128];
    char code[128];
    int first, last;
    int numR;
}; 

void Color2Seq::init(int alen)
{
    assemlen = alen;
    first = alen+1;
    last = 0;	
    Co = new int*[3*assemlen];
    Ba = Co+assemlen;
    pa = Ba+assemlen;
    seq = new char[assemlen*2+2];
    cseq = seq+assemlen+1;
    Co[0] = new int[4*3*assemlen];
    memset(Co[0], 0, sizeof(int)*4*3*assemlen);
    Ba[0] = Co[0]+4*assemlen;
    pa[0] = Ba[0]+4*assemlen; 
    color = new char*[assemlen];
    color[0] = new char[5*assemlen];
    notN = color[0]+4*assemlen;
    memset(color[0], 0, 5*assemlen);
    int i;
    for (i = 1; i < assemlen; i++) {
	Co[i] = Co[i-1]+4; 
  	Ba[i] = Ba[i-1]+4;
	pa[i] = pa[i-1]+4;
	color[i] = color[i-1]+4;
    }

    color2seq['A']['0'] = 'A';
    color2seq['A']['1'] = 'C';
    color2seq['A']['2'] = 'G';
    color2seq['A']['3'] = 'T';
    color2seq['C']['0'] = 'C';
    color2seq['C']['1'] = 'A';
    color2seq['C']['2'] = 'T';
    color2seq['C']['3'] = 'G';
    color2seq['G']['0'] = 'G';
    color2seq['G']['1'] = 'T';
    color2seq['G']['2'] = 'A';
    color2seq['G']['3'] = 'C';
    color2seq['T']['0'] = 'T';
    color2seq['T']['1'] = 'G';
    color2seq['T']['2'] = 'C';
    color2seq['T']['3'] = 'A';
    init_scode(code);
}

Color2Seq::~Color2Seq()
{
    delete [] Co[0];
    delete [] seq;
    delete [] Co;
    delete [] color[0];
    delete color;
}

void Color2Seq::outputSeq()
{
    int i;
    int CinLine=0;
    if (numR == 0) return;
    if (last >= assemlen) fatal("contig longer than maximum\n");
    pa[last+1][0] = pa[last+1][1] =pa[last+1][2]=pa[last+1][3]=0;
    char c[] = "ACGT";
    for (i = last; i >= first; i--) 
      {
	     int j;
	     for (j = 0; j < 4; j++) 
	      {
	       char B = c[j];
	       char N;
	       int k, min=500000000, index = 0;
	      
	        for (k = 0; k < 4; k++) 
	         {
		N = color2seq[B]['0'+k];
		int x = code[N];
	        	int sc = pa[i+1][x]+Co[i][k];
	        	if (k==0 || sc < min){ min = sc; index =x;} 
	         }  
	      pa[i][j] = min+Ba[i][j];
	      color[i][j] = index;
	     }
	    if (pa[i][0] > 10000000) 
	    {
	     pa[i][1] -= pa[i][0];
	     pa[i][2] -= pa[i][0];
                    pa[i][3] -= pa[i][0];
    	     pa[i][0] = 0;
	    }	    
    }
    int min = pa[first][0], w = 0;
    for (i = 1; i < 4; i++) {
	if (min > pa[first][i]) { min = pa[first][i]; w = i;}
    }
    printf("%c", c[w]); CinLine = 1;
    
    for (i = first+1; i <= last; i++) {
	  w = color[i-1][w];
        if (CinLine >= ChrPerLine) { printf("\n"); CinLine=0;}
        if(!notN[i]) { printf("%c",'N'); while(++i<=last && !notN[i]);i--;}
        else printf("%c", c[w]);
	  CinLine++;
    }
    printf("\n");
}

int Color2Seq::inputBaseReads(FILE *fp, int cost)
{
    if (fp == NULL) return 0;
    char line[10000], line1[1000];
    int y = 0;
    //numR = 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] != '=') continue;
        //printf(">%s", line+1);
        y = 1;
        break;
    }
    if (y == 0) return 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '#') continue;
        if (line[0]== '-') {
            break;
        }
        if (line[0] != '>') fatal("format of input wrong\n");
        if (!fgets(line1, sizeof line1, fp)) break;
        char *c = strchr(line, ',');
        if (!c) continue;
        int pos = atoi(c+1);
        int p = abs(pos);
        int len = strcspn(line1, "\n\r");
        if (pos >= 0) {
            if (pos < first) first = pos;
            if (pos+len-1>last) last = pos+len-1;
        } else {
            if (-pos > last) last = -pos;
            if (-pos-len+1 < first) first = -pos-len+1;
        }
	if (first < 0 || last >= assemlen) 
                fatal("contig exceed maximum length or reads match to negative position\n");

	int i;
	for (i = 0; i < len; i++) {
	    if (pos >= 0) { 
		p = pos+i;
	    } else {
		p = -pos-i;
	    }
	    Ba[p][0]+=cost; Ba[p][1]+=cost; Ba[p][2]+=cost; Ba[p][3]+=cost;
	    int x = code[line1[i]];
	    if (pos < 0) x = 3-x;
	    Ba[p][x]-=cost;
	}
	numR++;
    }
    return 1;
}


int Color2Seq::inputReads(FILE *fp, int cost)
{
    char line[10000], line1[1000];
    int y = 0;
    numR = 0;
    int node_id = 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] != '=') continue;
	printf(">%s", line+1);
	
	string sline(line);
	string::size_type dotpos, id_pos;
	if((dotpos = sline.find("."))!=string::npos) init(atoi(line + dotpos+1)+200);
	else init(5000000); 
	y = 1;
	break;
    }
    if (y == 0) return 0;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') continue;
	if (line[0]== '-') {
	    break;
	}
        if (line[0] != '>') fatal("format of input wrong\n");
        if (!fgets(line1, sizeof line1, fp)) break;
	char *c = strchr(line, ',');
	if (!c) continue;
	int pos = atoi(c+1);
	int p = abs(pos);
               int len = strcspn(line1, "\n\r");
 
       if (pos >= 0) {
	    if (pos < first) first = pos;
	    if (pos+len-2>last) last = pos+len-2; 
	} else {
	    if (-pos > last) last = -pos;
	    if (-pos-len+2 < first) first = -pos-len+2;
	}
	Ba[p][0]+=cost; Ba[p][1]+=cost; Ba[p][2]+=cost; Ba[p][3]+=cost;
	int x = code[color2seq[line1[0]][line1[1]]];
	if (pos < 0) x = 3-x;
	Ba[p][x]-=cost; 
	notN[p] = 1;
	int i;
	for (i = 2; i < len; i++) {
	    if (pos >= 0) { 
		    p = pos +i-2;
	                   notN[p+1] = 1;
	                   } 
	               else 
	                  {
                                   p = -pos-i+1;
		     notN[p] = 1;
	                  }
	     if (p < assemlen  && p >= 0) 
	      {
	     	if (line1[i] > '3' || line1[i] < '0') continue;
		    Co[p][0]+=cost; Co[p][1]+=cost; Co[p][2]+=cost; Co[p][3]+=cost;
		    Co[p][line1[i]-'0']-=cost;
	    } else {
		cerr << "contig exceed maximum length or reads match to negative position\n";
	    }
	   }
	numR++;
 }
 
return 1;
}

//----------------------------------------------------------------------
void convert_aligned_reads_into_ntcontigs(string readsfile, int outSeqLen, string base_reads, int colorCost, int BaseCost)
{
  // The code was inherited from denovoadp (V2.3 Oct 21,2009) readsfile max_assem_length [B=base_reads][C=##][M=##]
  // Here C is the penalty for color, M is the penalty for base in a base read, defaults for both 1

    FILE *fp = ckopen(readsfile.c_str(), "r"); 
    FILE *fpp = (base_reads != ""?fpp =ckopen(base_reads.c_str(), "r"):NULL);

    do { 
       	Color2Seq *g = new Color2Seq(outSeqLen);
      	int n;
    	  if ((n=g->inputReads(fp, colorCost))==0) break; 
	      g->inputBaseReads(fpp, BaseCost);
    	  g->outputSeq();
	      delete g;
    } while (true);
}
//-----------------------------------------------------------------------


