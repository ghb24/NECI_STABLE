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
  if (argc != 4) { explain(); exit(-1); }
  //  int norbs = atoi(argv[1]);
  //int spatnorbs = norbs / 2;
  ofstream intsfile(argv[3]);
  ifstream oneintsf(argv[1]);
  ifstream twointsf(argv[2]);

  string line;
  int msgsize = 5000;


  int i, j, k, l;
  double val;
  intsfile.precision(20);
  bool switched = false;

  int norbs;
  twointsf >> norbs;
  oneintsf >>norbs;
  while (twointsf >> i >> j >> k >> l >> val)
    /* note i, j, k, l (i.e. lowercase) are the MOLPRO indices which start from 1 */
    {
	
      /* DMRG indices (starting from 0) */
      int I = i + 1;
      int J = j + 1;
      int K = k + 1;
      int L = l + 1;

      intsfile <<val<<"  "<< I<<"  "<<K<<"  "<<J<<"  "<<L<<"  "<<endl;
    }
  while(oneintsf >> i>> j>> val)
    {
      int I = i+1;
      int J = j+1;
      intsfile<<val<<"  "<<I<<" "<< J<<" 0  0"<<endl;
	
    }
  intsfile.close();
}

