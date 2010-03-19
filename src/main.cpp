#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <boost/multi_array.hpp>

extern "C" {
#include "pf_duplex.h"
};

typedef unsigned int uint;

int
main(int argc, char* argv[])
{
  enum {M=0, IX=1, IY=2};
  
  //float gamma=1.0;
  //float th=1/(gamma+1.0);
  float th=0.0;
  std::string s1="acguacgu";
  std::string s2="ugcaugca";

  uint n1 = s1.size();
  uint n2 = s2.size();

  // calculate posterior probabilities and copy them
  pf_duplex(s1.c_str(), s2.c_str());
  boost::multi_array<float, 2> ptable(boost::extents[n1+1][n2+1]);
  for (uint i=1; i!=n1+1; ++i)
    for (uint j=1; j!=n2+1; ++j)
      ptable[i][j] = pr_duplex[i][j];
  free_pf_duplex();

  boost::multi_array<float, 2> dp(boost::extents[n1+1][n2+1]);
  boost::multi_array<int, 2> bt(boost::extents[n1+1][n2+1]);

  // initialize
  dp[1][n2]=ptable[1][n2];
  bt[1][n2]=ptable[1][n2]>th ? M : IX;
  // j=n2
  for (uint i=2; i!=n1+1; ++i)
  {
    dp[i][n2]=dp[i-1][n2];
    bt[i][n2]=IX;
    if (ptable[i][n2]>th && dp[i][n2]<ptable[i][n2])
    {
      dp[i][n2]=ptable[i][n2];
      bt[i][n2]=M;
    }
  }
  // i=1
  for (uint j=n2-1; j!=0; --j)
  {
    dp[1][j]=dp[1][j+1];
    bt[1][j]=IY;
    if (ptable[1][j]>th && dp[1][j]<ptable[1][j])
    {
      dp[1][j]=ptable[1][j];
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
      if (ptable[i][j]>th && dp[i][j]<dp[i-1][j+1]+ptable[i][j])
      {
        dp[i][j]=dp[i-1][j+1]+ptable[i][j];
        bt[i][j]=M;
      }
    }
  }

  // trace back
  std::string r1(s1.size(), '.');
  std::string r2(s2.size(), '.');
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

  std::cout << s1 << std::endl << r1 << std::endl
            << s2 << std::endl << r2 << std::endl;

  return 0;
}
