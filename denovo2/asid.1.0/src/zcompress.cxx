#include "zcompress.h"

genFile::genFile ()
{
  isNull=true;
  usePipe=false;
  ncFILE = NULL;
  isZ = false;
}

bool genFile::open (const char *fname,const char *direction) {
  filename.assign(fname);
  if (direction[0]=='r') {return openR(fname,direction);}
  if (direction[0]=='w') {return openW(fname,direction);}
  return false;
}

bool genFile::openW (const char *fname,const char *direction) {
  std::string fn(fname);
  isZ = false;
  isNull = false;
  if (strcmp(fname, "stdout")==0) { ncFILE = stdout; return true;}
  char baseFn[5000];
  isZ = remove_gz_ext(fname,baseFn);  //see if it ends with .gz

  if(isZ) {
    if(usePipe) {
      std::string cmd = "gzip -c -1 > " + fn;
      ncFILE = popen(cmd.c_str(), "w");
      //fprintf(stderr,"detected as gz. (pipe mode)\n");
    } else {
      gzFILE = gzopen(fname, direction);
      if(gzFILE == NULL) {
	isNull = true;
	return false;
      }
    }
  } else {
    ncFILE = fopen(fname, direction);
    if(ncFILE == NULL) {
      isNull = true;
      return false;
    }
  }
  return true;
}

bool genFile::openR (const char *fname,const char *direction)
{
  ncFILE = fopen(fname, direction);
  if (strcmp(fname, "stdin") ==0) {
	isNull = false; ncFILE = stdin;
	return true;
  }
  std::string fn(fname);  //Not doing "using std" b/c std functions may be
                          //custom defined as something else
  if(ncFILE == NULL) {
    int lng = fn.length();
    fn += ".gz";  // see if fn.gz exists
    ncFILE = fopen(fn.c_str(), "r");
    if(ncFILE == NULL && lng>=4) {
      fn = fname;
      lng = fn.length();
      std::string fn_ext = fn.substr(lng-3,3);
      if(fn_ext==".gz") {
	fn = fn.substr(0,lng-3); //see if fn without .gz exists
	ncFILE = fopen(fn.c_str(), "r");
      }
    }
  }
  if( ncFILE == NULL ) {
    fprintf (stderr, "Can't open file '%s'\n", fn.c_str());
    isNull = true;
    return false;
  }

  unsigned char magic[2];
  static unsigned char gzmagic[] = { '\037', 0213 };
  fread(magic, 1, 2, ncFILE);
  if(memcmp(magic, gzmagic, 2) == 0) {
    // File is gzip compressed.
    fclose (ncFILE);
    // close standard stream and open gz stream

    //fprintf(stderr,"fname: %s, real-fn(if different): %s, fmode: %s, was ",
	    //fname,fn.c_str(),direction);

    if(usePipe) {
      std::string cmd = "gzip -cd " + fn;
      ncFILE = popen(cmd.c_str(), "r");
      //fprintf(stderr,"detected as gz. (pipe mode)\n");
      isNull=(ncFILE==NULL);
    } else {
      gzFILE = gzopen(fn.c_str(), "rb");
      //fprintf(stderr,"detected as gz (gzopen mode)\n");
      isNull=(gzFILE==NULL);
    }
    isZ=true;
    //fprintf(stderr,"m_file.isNull: %d\n", isNull);
    return !isNull;
  } else {
    //fprintf(stderr,"detected as text.\n");
    rewind(ncFILE);
    isZ=false;
    isNull=false;
    return true;
  }
}

void* genFile::gets(char *line, int lineSize) {
  if (isNull) return NULL;
  if (isZ && !usePipe) {
    return gzgets(gzFILE,line,lineSize);
  } else {
    return fgets(line,lineSize,ncFILE);
  }
}

z_off_t genFile::tello () {
  if (isZ && !usePipe)
    return gztell(gzFILE);
  else
    return ftello(ncFILE);
}

void genFile::gFrewind () {
  if (isZ && !usePipe) {
    gzrewind(gzFILE);
  } else {
    rewind(ncFILE);
  }
}

int genFile::seek (long offset, int origin) {
  if (isZ && !usePipe) {
    return gzseek(gzFILE,offset,origin);
  } else {
    return fseek(ncFILE,offset,origin);
  }
}

int genFile::printf(const char *format,...) { //this really doesn't work
  va_list argList;                   //as planned, b/c no vgzprintf function
  va_start(argList,format);
  if (isZ && !usePipe) {
    printf("Cannot write to gz file yet.\n");
    return 0; //vgzprintf(gzFILE,format,argList);
  } else {
    return vfprintf(ncFILE,format,argList);
  }
}

bool genFile::seekToPanelId (int panelId) { //, char* lineBUFFER) {

  if (panelId==0) {
	  return true;  //no seeking required
  }
  size_t file_position = 0; //will seek to this

  string index_fn = filename+".idx";
  ifstream index(index_fn.c_str());
  if(isZ || !index.is_open()) {
    //If it's gz, or if there is no index, then do a linear search for the bead id

    //Note, could do autogen of index here by using PanelIndexing::ensureIndexFile,
    //but that could cause performance issues at best (and wrong results at worst)
    //in a parallel environment.

    char line[3000000];  //*sing* We program like it's 1999. *end* (or rather 1983)
    bool found = false;

    gFrewind();

    while ( gets(line, sizeof line ) ) {
      if (line[0] != '#') break;  //get past header line if there
      file_position += strlen(line);
    }

    seek(file_position,0);  //go back to before line was read

    int curPanelId;
    while ( gets(line, sizeof line ) ) {
      if (line[0] == '>') {
	std::string bead_id(line);
	string::size_type loc = bead_id.find( '_', 0);
	bead_id = bead_id.substr(1,loc-1);
	curPanelId = atoi(bead_id.c_str());
	if (curPanelId==panelId) {
	  found = true;
	  break;
	}
      }
      file_position += strlen(line);
    }
    if(found)
      seek(file_position,0);
    return found;
  } else {
    std::cerr << "In seekToPanelId, using index.\n";
    PanelIndexing panelIndex;
    panelIndex.readIndexFile(index_fn);
    int vectLoc = //panelIndex.get_vect_pos_greater_equal_panel_id(panelId);
      panelIndex.get_vect_pos_panel(panelId);
    file_position = panelIndex.get_panel_file_position(vectLoc);
    seek (file_position,0);
  }
  index.close();

  return true;
}

void genFile::close() {
  if(isNull!=true) {
    isNull=true;
    if (isZ && !usePipe ) {
      gzclose(gzFILE);
    } else {
	if (ncFILE == stdout || ncFILE == stdin) return;
	fflush(ncFILE);
      fclose(ncFILE);
    }
  }
}

bool remove_gz_ext (const char *fname, char *no_gz_fname) {
  bool isGzExt = false;
  std::string fn = fname;
  int lng = fn.length();
  std::string fn_ext = fn.substr(lng-3,3);
  if(fn_ext==".gz") {
    fn = fn.substr(0,lng-3); //see if fn without .gz exists
    isGzExt = true;
  }
  strcpy(no_gz_fname,fn.c_str());
  return isGzExt;
}

PanelIndexing::PanelIndexing() {
  min_panel_id = -1;
  max_panel_id = -1;
  max_number_of_panels = -1;
}

PanelIndexing::~PanelIndexing() {

}

void PanelIndexing::streamOutReadsFile(genFile &csfasta_file, int minPanelNum, int maxPanelNum) const {
  if(minPanelNum>0) {
    std::vector<int>::const_iterator panel_idsItr=panel_ids.begin();
    for(;panel_idsItr!=panel_ids.end();panel_idsItr++) {
      if(*panel_idsItr >= minPanelNum ) {
	minPanelNum = *panel_idsItr;
	break;
      }
    }
    if(panel_idsItr==panel_ids.end()) {
      std::cerr << "Panel id " << minPanelNum << " is past the end of the file.\n";
      exit(1);
    }
   csfasta_file.seekToPanelId(minPanelNum);
  }

  char line[1000000];
  while ( csfasta_file.gets(line, sizeof line ) ) {
    if (line[0] == '>') {
      int curPanelId;
      std::string bead_id(line);
      string::size_type loc = bead_id.find( '_', 0);
      bead_id = bead_id.substr(1,loc-1);
      curPanelId = atoi(bead_id.c_str());
      if(maxPanelNum>=0 && curPanelId>maxPanelNum) break;
    }
    std::cout << line;
  }
}


void PanelIndexing::determinePanelFilePositions(genFile &csfasta_file) {

  size_t file_position = 0;
  char line[1000000];
  csfasta_file.gFrewind();

  //size_t last_file_position = file_position;
  while ( csfasta_file.gets(line, sizeof line ) ) {
    if (line[0] != '#') break;  //get past header line if there
    file_position += strlen(line);
  }

  csfasta_file.seek(file_position,0);

  int lastPanelId = -1;
  int curNumbBeads = -2;

  while ( csfasta_file.gets(line, sizeof line ) ) {
    if (line[0] == '>') {
      int curPanelId;
      curNumbBeads++;
      std::string bead_id(line);
      string::size_type loc = bead_id.find( '_', 0);
      bead_id = bead_id.substr(1,loc-1);
      curPanelId = atoi(bead_id.c_str());
      if (curPanelId!=lastPanelId) {
	//std::cerr << curPanelId << " " << lastPanelId << " "
	//	  << file_position << " "
	//	  << line << "\n";
	lastPanelId=curPanelId;
	panel_ids.push_back(curPanelId);
	if (curNumbBeads!=-1) { //this is calc. 1 off phase from rest
	  num_beads_in_panel.push_back(curNumbBeads);
	  panel_end_file_positions.push_back(file_position);
	}
	panel_file_positions.push_back(file_position);
	curNumbBeads = 0;
      }
    }
    file_position += strlen(line);
  }

  num_beads_in_panel.push_back(curNumbBeads); //last group's number of beads
  panel_end_file_positions.push_back(file_position); // and end file position

  last_panel_id = lastPanelId;

  index_fn = csfasta_file.getFilename() + ".idx";
  csfasta_file.gFrewind();
}

void PanelIndexing::makeIndexFile() {
  std::cerr << "Making index " << index_fn << "\n";

  string index_part_fn = index_fn+".part";
  ofstream index(index_part_fn.c_str());
  if(!index.is_open()) {
    std::cerr << "Cannot open file " << index_fn << "\n";
    exit(1);
  }
  index << "# Index created as " << index_fn << "\n";
  index << "##PanelId Numb_beads_in_panel FilePositionNumber FilePositionEnd\n";
  printVectorQuad(panel_ids,num_beads_in_panel,
		  panel_file_positions,
		  panel_end_file_positions,
		  index);
  if(index.fail()) {
    std::cerr << "Writing failed to file " << index_fn << "\n";
    exit(1);
  }
  index.close();
  int didCompleteTag = rename(index_part_fn.c_str(),index_fn.c_str());
  if( didCompleteTag != 0 ) {
    std::cerr << "Cannot rename file " << index_part_fn << " to " << index_fn << ".\n";
    exit(1);
  }
}

void PanelIndexing::readIndexFile(const string &in_index_fn) {
  panel_ids.clear();
  panel_file_positions.clear();
  panel_end_file_positions.clear();
  index_fn = in_index_fn;
  ifstream index(index_fn.c_str());
  if(!index.is_open()) {
    std::cerr << "Cannot open file " << index_fn << "\n";
    exit(1);
  }

  string line;
  while (! index.eof() )   {
    getline (index,line);
    if (line.size()>0) {
      if (line.at(0) != '#') {  //ignore comment lines
	//good fast split function could be used here instead
	size_t parse_loc = line.find_first_of(" ");
	size_t parse_loc2 = line.substr(parse_loc+1,line.size()).find_first_of(" ")
	  + parse_loc + 1;
	size_t parse_loc3 = line.find_last_of(" ");
	if(parse_loc==parse_loc2 || parse_loc2==parse_loc3 || parse_loc==parse_loc3) {
	  std::cerr << in_index_fn << " is corrupted.  Try deleting it\n";
	  exit(1);
	}
	int panel_id = atoi(line.substr(0,parse_loc).c_str());
	//int num_beads_in_panel = atol(line.substr(parse_loc+1,parse_loc2).c_str());
	size_t file_position = atol(line.substr(parse_loc2+1,parse_loc3).c_str());
	size_t end_file_position = atol(line.substr(parse_loc3+1,line.size()).c_str());

	//std::cerr << parse_loc << "Regular line: " << line << "\n";
	//std::cerr << "'" << parse_loc << "' '" << panel_id << "' '" << position << "'\n";

	bool doInsert = true;
	if(panel_id < min_panel_id) doInsert = false;
	if(max_panel_id > 0 && max_panel_id < panel_id) doInsert = false;
	if (doInsert) {

	  panel_ids.push_back(panel_id);
	  panel_file_positions.push_back(file_position);
	  panel_end_file_positions.push_back(end_file_position);

	  if(max_number_of_panels==(int)panel_ids.size()) {
	    index.close();
	    return;
	  }
	}
      }
    }
  }
  index.close();

  //printVectorPair(panel_ids,panel_file_positions,std::cerr);
}

void PanelIndexing::ensureIndexFile(genFile &csFasta_File) {
  //If there's an index file, then it will read it, if not, it will make it.
  string index_fn_in = csFasta_File.getFilename() + ".idx";
  ifstream index(index_fn_in.c_str());
  if(!index.is_open()) { //no file, so create it
    std::cerr << "Creating index.\n";
    determinePanelFilePositions(csFasta_File);
    makeIndexFile();
    readIndexFile(index_fn_in);
  } else { //just read it
    std::cerr << "Reading index.\n";
    readIndexFile(index_fn_in);
  }
  //printVectorPair(panel_ids,panel_file_positions,std::cerr);
}

void PanelIndexing::reduceToSubset(const vector<int> &panel_ids_to_keep) {
  vector<size_t> new_panel_file_positions, new_panel_end_file_positions;
  vector<int>::const_iterator panel_ids_to_keepItr = panel_ids_to_keep.begin();
  for(;panel_ids_to_keepItr!=panel_ids_to_keep.end();panel_ids_to_keepItr++) {
    int vect_pos = get_vect_pos_panel(*panel_ids_to_keepItr);
    new_panel_file_positions.push_back(panel_file_positions.at(vect_pos));
    new_panel_end_file_positions.push_back(panel_end_file_positions.at(vect_pos));
  }
  panel_file_positions = new_panel_file_positions;
  panel_end_file_positions = new_panel_end_file_positions;
  panel_ids = panel_ids_to_keep;
}

void PanelIndexing::reduceToSubset(int subsetSize) {
  if(subsetSize>(int)panel_ids.size()) {
    std::cerr << "Note, subsetSize (" << subsetSize <<
      ") is larger than total size, using all panels (within min/max panel id if specified).\n";
    return;
  }
  panel_file_positions.erase(panel_file_positions.begin()+subsetSize,panel_file_positions.end());
  panel_end_file_positions.erase(panel_end_file_positions.begin()+subsetSize,panel_end_file_positions.end());
  panel_ids.erase(panel_ids.begin()+subsetSize,panel_ids.end());
}

vector<int> PanelIndexing::reduceToRandomSubset(int subsetSize) {
  int numToErase = panel_ids.size() - subsetSize;
  if (numToErase<0) {
    std::cerr << "In making random subset of panel ids, there's only " << panel_ids.size() << " panels present, so can't use a subset size of " << subsetSize << "\n";
    exit(-1);
  }
  vector<int> rand_panel_ids = panel_ids;
  random_shuffle ( rand_panel_ids.begin(), rand_panel_ids.end() );
  rand_panel_ids.erase(rand_panel_ids.end()-numToErase,rand_panel_ids.end());
  sort(rand_panel_ids.begin(),rand_panel_ids.end());

  reduceToSubset(rand_panel_ids);
  return rand_panel_ids;
}

void PanelIndexing::intersection(PanelIndexing &panelIndexing2) {
  vector<int> panel_ids_common;
  for(int i=0; i<(int)panel_ids.size(); i++) panel_ids_common.push_back(-1);
  vector<int>::iterator commonItr;
  commonItr = set_intersection (panel_ids.begin(),panel_ids.end(),
				panelIndexing2.panel_ids.begin(),
				panelIndexing2.panel_ids.end(),
				panel_ids_common.begin());
  panel_ids_common.erase( commonItr,panel_ids_common.end() );
  reduceToSubset(panel_ids_common);
  panelIndexing2.reduceToSubset(panel_ids_common);
}

void PanelIndexing::determinePanelStartsEnds( const int numSplits ) {
  //Make groups of panels that have as even as possible number of panels each

  panel_group_starts.clear();
  panel_group_ends.clear();

  int numDiffPanels = panel_ids.size();
  if (numSplits>numDiffPanels) {
    std::cerr << "There are not enough panels (" << numDiffPanels << ") to split "
	      << numSplits << " ways in " << index_fn << "\n";
    exit(1);
  }

  //0-0, 1-1 , 2-2,  3-3
  //0-1, 2-3,  4-5,  6-7

  int remainder = numDiffPanels % numSplits;
  int minNumPanelsPerGroup = (numDiffPanels - remainder) / numSplits;

  //grouping evenly by number of panels, could also do it weighted for # of beads
  int groupNum;
  int remainder_adder = 0;
  for(groupNum=0;groupNum<numSplits; groupNum++) {
    int startPanel = panel_ids.at(groupNum * minNumPanelsPerGroup + remainder_adder);
    if (remainder_adder< remainder) remainder_adder++;
    int endPanel = panel_ids.at((groupNum+1) * minNumPanelsPerGroup + remainder_adder - 1);
    panel_group_starts.push_back(startPanel);
    panel_group_ends.push_back(endPanel);
  }

  std::cerr << "# Panels, remainder, min#/group: " << numDiffPanels << " "
	    << remainder << " " << minNumPanelsPerGroup << "\n";
  return;
}

void PanelIndexing::keepNonPairedPanels () {
	vector<int>::iterator panelStartItr = panel_group_starts.begin();
    vector<int>::iterator panelEndsItr = panel_group_ends.begin();
    for(;panelStartItr!=panel_group_starts.end();) {
    	if (panelStartItr == panel_group_starts.begin()) {
    		*panelStartItr = 0;
    	}
    	if (panelEndsItr != panel_group_ends.end()-1) {
    		*panelEndsItr = (*(panelStartItr+1)) -1;
    	} else {
    		*panelEndsItr = 99999999;
    	}
        panelStartItr++;
        panelEndsItr++;
    }
}

void PanelIndexing::determinePanelStartsEnds( const int numSplits, ostream &out) {
  determinePanelStartsEnds(numSplits);
  printVectorPair(panel_group_starts,panel_group_ends,out);
}

void PanelIndexing::printStartEndGrouping( ostream &out, bool printFilePos ) const {
  if(printFilePos) {
    vector<int>::const_iterator panelStartItr = panel_group_starts.begin();
    vector<int>::const_iterator panelEndsItr = panel_group_ends.begin();
    for(;panelStartItr!=panel_group_starts.end();) {
      out << *panelStartItr << " " << *panelEndsItr << " "
	  << get_panel_file_position_by_panel(*panelStartItr) << " "
	  << get_panel_end_file_position_by_panel(*panelStartItr) << " "
	  << "\n";
      panelStartItr++;
      panelEndsItr++;
    }
  } else {
    printVectorPair(panel_group_starts,panel_group_ends,out);
  }
}

void PanelIndexing::printAllPanelsInGroup( ostream &out, bool printFilePos ) const {
  std::vector<int>::const_iterator panelStartItr = panel_group_starts.begin();
  std::vector<int>::const_iterator panelEndsItr = panel_group_ends.begin();
  for(;panelStartItr!=panel_group_starts.end();) {
    std::vector<int>::const_iterator startItr =
      panel_ids.begin() + get_vect_pos_panel(*panelStartItr);
    std::vector<int>::const_iterator endItr =
      panel_ids.begin() + get_vect_pos_panel(*panelEndsItr);
    printVector(panel_ids," ", startItr, endItr+1);
    if (printFilePos) {
      out << " " << get_panel_file_position_by_panel(*panelStartItr);
      out << " " << get_panel_end_file_position_by_panel(*panelStartItr);
    }
    out << "\n";
    panelStartItr++;
    panelEndsItr++;
  }
}

int PanelIndexing::get_vect_pos_panel(int panel_id) const {
  //int num_beads = -1;
  std::vector<int>::const_iterator panelLocIter = find(panel_ids.begin(),
						       panel_ids.end(),
						       panel_id);
  if (panelLocIter == panel_ids.end()) {
    std::cerr << "Cannot find panel number " << panel_id
	      << ".\nExiting.\n";
    exit(1);
  }
  std::vector<int>::const_iterator counter = panel_ids.begin();

  int panelLoc = 0;
  for(;counter!=panelLocIter;counter++) panelLoc++;
  return panelLoc;
  //returns where in the vector (i.e. the 17th one) a panel number is at
}

bool comparison_gte ( int i1, int i2 ) {
  return ( i1 >= i2 );
}

int PanelIndexing::get_vect_pos_greater_equal_panel_id ( int search_panel_id ) const {

  int search[] = { search_panel_id };

  vector<int>::const_iterator panel_idsItr;
  panel_idsItr = find_first_of(panel_ids.begin(),panel_ids.end(),search,search+1,
			       comparison_gte );

  ;//find_if(panel_ids.begin(),panel_ids.end(),greater);

  std::vector<int>::const_iterator counter = panel_ids.begin();

  int panelLoc = 0;
  for(;counter!=panel_idsItr;counter++) panelLoc++;
  return panelLoc;
}

std::vector<string> split_string (const string &inputString, const string &delim) {
  vector<string> parts;
  size_t curLoc = 0;
  while( true ) {
    size_t delimLoc = inputString.find( delim[0], curLoc );
    //cout << delimLoc << "\n";
    if( delimLoc == inputString.npos ) {
      parts.push_back( inputString.substr( curLoc, inputString.size() ));
      break;
    }
    parts.push_back( inputString.substr( curLoc, delimLoc - curLoc )) ;
    curLoc = delimLoc + 1;
  }
  //cout << "input string: " << inputString << " delim: " << delim << "\n";
  return parts;
}

std::vector<int> split_int (const string &inputString, const string &delim) {
  vector<string> parts = split_string(inputString,delim);

  vector<int> int_parts;
  int_parts.reserve(parts.size());

  vector<string>::const_iterator partsItr = parts.begin();
  for(;partsItr!=parts.end();partsItr++) {
    int_parts.push_back(atoi(partsItr->c_str()));
  }

  return int_parts;
}
