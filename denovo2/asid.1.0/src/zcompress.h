#ifndef __zcompress_h__
#define __zcompress_h__ 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <math.h>
#include <zlib.h>
#include <stdarg.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using std::string;
using std::vector;
using std::ofstream;
using std::ostream;
using std::ifstream;

//#include "boost/fdstream.hpp"

  /* gzFile is the Zlib equivalent of FILE from stdio */
  //gzFile file;
  /* Opens up the file with Zlib */
  //file = gzopen("penguin.gz", "r");

class genFile {
 public:
  genFile();
  bool open (const char *fname,const char *direction);
  bool open (const string fname,const char *direction) { return open(fname.c_str(),direction); };
  bool openW (const char *fname,const char *direction);
  bool openR (const char *fname,const char *direction);

  bool getIsNull() const { return isNull; };
  bool getIsZ() const { return isZ; };
  bool getUsePipe() const { return usePipe; };
  std::string getFilename() const { return filename; };

  void setUsePipe() { usePipe=true; };
  void setUsePipe(const int setting) { usePipe=((setting==0)?false:true); };
  void* gets(char *line, int lineSize);
  z_off_t tello ();

  void gFrewind ();
  int seek (long offset, int origin);

  int printf(const char *format,...); // not implemented using zlib, use pipe for gz

  bool seekToPanelId (int panelId);  //return 0/false if unsucessful.

  void close();
  FILE *getFILE() { return ncFILE; };
  FILE *get_ncFILE() { return ncFILE; };
  gzFile &get_gzFILE() { return gzFILE; };

 private:
  std::string filename;
  bool isNull;
  bool isZ;
  bool usePipe;
  gzFile gzFILE;
  FILE *ncFILE;
  //static const
};

bool remove_gz_ext (const char *fname, char *no_gz_fname);

template <class T>
void printVector(const std::vector<T> vect, string seperator, ostream &out) {
  for(typename vector< T >::const_iterator i = vect.begin(); i!= vect.end(); ++i) {
    if (i != vect.end() - 1) {
      out << *i << seperator;
    } else {
      out << *i;
    }
  }
}

template <class T>
void printVector(const std::vector<T> vect, string seperator,
		 typename std::vector<T>::const_iterator begin,
		 typename std::vector<T>::const_iterator end) {
  for(typename vector< T >::const_iterator i = begin; i!= end; ++i) {
    if (i != end - 1) {
      std::cout << *i << seperator;
    } else {
      std::cout << *i;
    }
  }
}

//Some utility functions for vectors
template <class T1, class T2>
void printVectorPair(const std::vector<T1> vect1,const std::vector<T2> vect2, ostream &out) {
  typename vector< T1 >::const_iterator i1 = vect1.begin();
  typename vector< T2 >::const_iterator i2 = vect2.begin();
  for(; i1!= vect1.end() && i2!= vect2.end();) {
    out << *i1 << " " << *i2 << std::endl;
    i1++; i2++;
  }
}

template <class T1, class T2, class T3>
void printVectorTrio(const std::vector<T1> vect1,const std::vector<T2> vect2,const std::vector<T3> vect3, ostream &out) {
  typename vector< T1 >::const_iterator i1 = vect1.begin();
  typename vector< T2 >::const_iterator i2 = vect2.begin();
  typename vector< T3 >::const_iterator i3 = vect3.begin();
  for(; i1!= vect1.end() && i2!= vect2.end() && i3!= vect3.end();) {
    out << *i1 << " " << *i2 << " " << *i3 << std::endl;
    i1++; i2++; i3++;
  }
}

template <class T1, class T2, class T3, class T4>
void printVectorQuad(const std::vector<T1> vect1,const std::vector<T2> vect2,const std::vector<T3> vect3, std::vector<T4> vect4, ostream &out) {
  typename vector< T1 >::const_iterator i1 = vect1.begin();
  typename vector< T2 >::const_iterator i2 = vect2.begin();
  typename vector< T3 >::const_iterator i3 = vect3.begin();
  typename vector< T4 >::const_iterator i4 = vect4.begin();
  for(; i1!= vect1.end() && i2!= vect2.end() && i3!= vect3.end() && i4!=vect4.end();) {
    out << *i1 << " " << *i2 << " " << *i3 << " " << *i4 << std::endl;
    i1++; i2++; i3++; i4++;
  }
}

template <class T>
T sum(const vector<T> &vect) {
  typename vector< T >::const_iterator itr=vect.begin();
  T sumRes = 0;
  for(;itr!=vect.end();itr++) {
    sumRes += (*itr);
  }
  return sumRes;
}

template <class T>
T sum(const vector<T> &vect, T nan_value) {
  typename vector< T >::const_iterator itr=vect.begin();
  T sumRes = 0;
  for(;itr!=vect.end();itr++) {
    if(*itr!=nan_value)
      sumRes += (*itr);
  }
  return sumRes;
}

class PanelIndexing {
 public:
  PanelIndexing();
  ~PanelIndexing();

  //Set parameters
  void setMinMaxPanelId(int min, int max) {
    min_panel_id = min; max_panel_id = max;
    if (min>max) {
      std::cerr << "Can't set min greater than max for panel id.\n";
      exit(-1);
    }
  };
  void setMaxNumberPanels(int maxNum) {
    max_number_of_panels = maxNum; };

  //Maker Index functions
  void determinePanelFilePositions(genFile &csFasta_File);
  void makeIndexFile();
  void readIndexFile(const string &in_index_fn);
  void ensureIndexFile(genFile &csFasta_File);

  //Modifiers (call before grouping)
  void reduceToSubset(const vector<int> &panel_ids_to_keep);
  void reduceToSubset(int subsetSize);

  vector<int> reduceToRandomSubset(int subsetSize);
  void intersection(PanelIndexing &panelIndexing2);

  //Grouping functions
  void determinePanelStartsEnds(const int numSplits);
  void determinePanelStartsEnds(const int numSplits, ostream &out);
  void keepNonPairedPanels();

  //output functions
  void printAllPanelsInGroup( ostream &out, bool printFilePos = false ) const;
  void printStartEndGrouping( ostream &out, bool printFilePos = false ) const;
  void streamOutReadsFile(genFile &csfasta_file, int minPanelNum, int maxPanelNum) const;

  //query functions
  int get_vect_pos_panel(int panel_id) const;
  int get_vect_pos_greater_equal_panel_id(int panel_id) const;
  int get_last_panel_id() const { return last_panel_id; };

  size_t get_panel_file_position(int pos) const
  {   return panel_file_positions.at(pos); };
  size_t get_panel_file_position_by_panel(int panel_id) const
  {   return panel_file_positions.at(get_vect_pos_panel(panel_id)); };
  size_t get_panel_end_file_position_by_panel(int panel_id) const
  {   return panel_end_file_positions.at(get_vect_pos_panel(panel_id)); };
  int get_panel_id(int pos) const {   return panel_ids.at(pos); };

  int get_num_beads(int pos) const { return num_beads_in_panel.at(pos); };
  int get_num_beads_panel(int panel_id) const {
    return get_num_beads(get_vect_pos_panel(panel_id));
  }
  int getNumbOfPanels() const { return panel_ids.size(); };

 private:
  std::vector<int> panel_ids;
  std::vector<size_t> panel_file_positions;
  std::vector<size_t> panel_end_file_positions;
  // to be precise, each end_file_position is the start position of the next group
  // but this is to keep consistency with other functions that define range as
  // (inclusive start) to (exclusive end).  The last position here is the size of the
  // file.

  std::vector<int> num_beads_in_panel;
  std::vector<int> panel_group_starts;
  std::vector<int> panel_group_ends;
  int last_panel_id;

  string index_fn;
  int min_panel_id;
  int max_panel_id;
  int max_number_of_panels;
};

std::vector<string> split_string (const string &inputString, const string &delim);

std::vector<int> split_int (const string &inputString, const string &delim);

// * zcompress library *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/zcompress.cxx

// * remduphitZ program *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE src/remduphitsZ.cxx -o ./remduphitsZ zutil.o zcompress.o -lz

// * fasta-io library *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/fasta-ioZ.cxx

// * mapZ library *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/mapZ.cxx
// * mapZ program *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE src/mapmain.cxx mapZ.o -o ./mapZ zcompress.o util.o zutil.o fasta-ioZ.o -lz

// * mapreads program *
//g++ -g -O2 -march=nocona -pipe -ffast-math  -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE zcompress.o src/mapreadsZ.cxx zutil.o util.o fasta-ioZ.o -o ./mapreadsZ -lz


#endif
