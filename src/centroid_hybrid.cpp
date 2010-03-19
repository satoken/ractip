#include <iostream>
#include <string>
#include <list>
#include <boost/multi_array.hpp>
#include "fa.h"

extern "C" {
#include "pf_duplex.h"
};

typedef unsigned int uint;

class Centroid
{
public:
  Centroid() {}
  ~Centroid() {}

  template < class PMatrix >
  void execute(const PMatrix& pm, uint n1, uint n2,
               float gamma, std::string& r1, std::string& r2) const
  {
    enum {M=0, IX=1, IY=2};

    float th=1.0/(gamma+1.0);
    boost::multi_array<float, 2> dp(boost::extents[n1+1][n2+1]);
    boost::multi_array<int, 2> bt(boost::extents[n1+1][n2+1]);

    // initialize
    dp[1][n2]=pm[1][n2];
    bt[1][n2]=pm[1][n2]>th ? M : IX;
    // j=n2
    for (uint i=2; i!=n1+1; ++i)
    {
      dp[i][n2]=dp[i-1][n2];
      bt[i][n2]=IX;
      if (pm[i][n2]>th && dp[i][n2]<pm[i][n2])
      {
        dp[i][n2]=pm[i][n2];
        bt[i][n2]=M;
      }
    }
    // i=1
    for (uint j=n2-1; j!=0; --j)
    {
      dp[1][j]=dp[1][j+1];
      bt[1][j]=IY;
      if (pm[1][j]>th && dp[1][j]<pm[1][j])
      {
        dp[1][j]=pm[1][j];
        bt[1][j]=M;
      }
    }
    // others
    for (uint i=2; i!=n1+1; ++i)
    {
      for (uint j=n2-1; j!=0; --j)
      {
        dp[i][j]=dp[i-1][j];
        bt[i][j]=IX;
        if (dp[i][j]<dp[i][j+1])
        {
          dp[i][j]=dp[i][j+1];
          bt[i][j]=IY;
        }
        if (pm[i][j]>th && dp[i][j]<dp[i-1][j+1]+pm[i][j])
        {
          dp[i][j]=dp[i-1][j+1]+pm[i][j];
          bt[i][j]=M;
        }
      }
    }

    // trace back
    r1.resize(n1);
    r2.resize(n2);
    std::fill(r1.begin(), r1.end(), '.');
    std::fill(r2.begin(), r2.end(), '.');

    uint i=n1, j=1;
    while (i>0 && j<=n2)
    {
      switch (bt[i][j])
      {
        case M:
          r1[i-1]='[';
          r2[j-1]=']';
          i--; j++;
          break;
        case IX:
          i--;
          break;
        case IY:
          j++;
          break;
      }
    }
  }
};

int
main(int argc, char* argv[])
{
  float gamma=1e10;
  std::list<Fasta> fa1, fa2;
  Fasta::load(fa1, argv[1]);
  Fasta::load(fa2, argv[2]);
  std::string s1=fa1.front().seq();
  std::string s2=fa2.front().seq();
  uint n1=s1.size();
  uint n2=s2.size();

  // calculate posterior probabilities and copy them
  pf_duplex(s1.c_str(), s2.c_str());
  Centroid cent;
  std::string r1, r2;
  cent.execute(pr_duplex, n1, n2, gamma, r1, r2);
  free_pf_duplex();

  std::cout << ">" << fa1.front().name() << std::endl
            << s1 << std::endl << r1 << std::endl
            << ">" << fa2.front().name() << std::endl
            << s2 << std::endl << r2 << std::endl;

  return 0;
}
