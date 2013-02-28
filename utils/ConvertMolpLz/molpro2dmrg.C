#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <vector>

using namespace std;
using namespace boost;

enum Int_t { AA, AB, BA, BB };

void explain()
{
  cout << "reads in molpro FCIDUMP and returns dmrg integral file (in ab/ab/ab ... order) <<  usage: exe dumpfile newoneints newtwoints  " << endl;
}

void ReadMeaningfulLine(ifstream& input, string& msg, int msgsize)
{
  msg.resize(0);
  bool readmore = true;
  while (readmore && !input.eof()) {
    char msgctr[msgsize];
    input.getline(msgctr, msgsize+1);

    msg=string(msgctr);
    if(msg.size() == msgsize) {
      cerr << "in the process of reading line begining with "<<endl<<msg<<endl;
      cerr<< "this line is too long"<<endl;
      abort();
    }
    int pos = msg.find("!");
    msg = msg.substr(0, pos);
    trim(msg);
    if (msg.size() != 0)
      readmore = false;
  }
}

int main(int argc, char** argv)
{
  cout << argc<<endl;
  
  if (argc != 4) { explain(); exit(-1); }
  //  int norbs = atoi(argv[1]);
  //int spatnorbs = norbs / 2;
  ifstream intsfile(argv[1]);
  ofstream oneintsf(argv[2]);
  ofstream twointsf(argv[3]);

  string line;
  int msgsize = 5000;

  ReadMeaningfulLine(intsfile, line, msgsize);
  vector<string> tok;
  boost::split(tok, line, is_any_of(" \t=,"), token_compress_on);
  int norbs = atoi(tok[2].c_str());
  cout <<"# orbitals = "<< norbs<<endl;

  //now read the next three dummy lines
  ReadMeaningfulLine(intsfile, line, msgsize);
  boost::split(tok, line, is_any_of(" \t=,"), token_compress_on);
  while (!iequals(tok[0], "&END")) {
    ReadMeaningfulLine(intsfile, line, msgsize);
    boost::split(tok, line, is_any_of(" \t=,"), token_compress_on);
  }

  int i, j, k, l;
  double val;
  oneintsf.precision(20);
  twointsf.precision(20);
  oneintsf<<norbs<<endl;
  twointsf<<norbs<<endl;
  bool switched = false;

  while (intsfile >> val >> i >> j >> k >> l)
    /* note i, j, k, l (i.e. lowercase) are the MOLPRO indices which start from 1 */
    {
      if ((i == 0) && (j == 0) && (k == 0) && (l == 0)) break;
	
      /* DMRG indices (starting from 0) */
      int I = i - 1;
      int J = j - 1;
      int K = k - 1;
      int L = l - 1;

      /* AA one ints */
      if ((k == 0) && (l == 0)) // one e ints
	{
	  oneintsf<<I<<" "<< J<<" "<<val<<endl;
	}
      /* AA two ints */
      else
	{
	  twointsf<<I<<" "<< K<<" "<<J<<" "<<L<<" "<<val<<endl;
	}
    }
  twointsf.close();
  oneintsf.close();
}

