/*
 * $Id$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of RactIP.
 *
 * RactIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RactIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RactIP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "cmdline.h"
#include <unistd.h>
#include <cstdlib>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include "fa.h"
#include "ip.h"

#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/DuplexEngine.hpp"
#include "contrafold/Defaults.ipp"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/part_func_up.h>
#include <ViennaRNA/part_func_co.h>
#include <ViennaRNA/utils.h>
#ifndef HAVE_VIENNA20
#include <ViennaRNA/PS_dot.h>
  extern int eos_debug;
#endif
#include "pf_duplex.h"
  extern void read_parameter_file(const char fname[]);
};
};

namespace uShuffle {
extern "C" {
#include "ushuffle.h"
};
};

extern "C" {
#include "boltzmann_param.h"
};

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
#ifdef HAVE_VIENNA18
typedef Vienna::plist pair_info;
#else
typedef Vienna::pair_info pair_info;
#endif

class RactIP
{
public:
  RactIP()
    : alpha_(0.5),
      beta_(0.0),
      th_hy_(0.2),
      th_ss_(0.5),
      th_ac_(0.0),
      max_w_(0),
      min_w_(0),
      enable_zscore_(0),
      num_shuffling_(0),
      seed_(0),
      in_pk_(true),
      use_contrafold_(true),
      use_pf_duplex_(true),
      stacking_constraints_(true),
      show_energy_(false),
      allow_concat_(true),
      run_with_modena_(false),
      n_th_(1),
      rip_file_(),
      param_file_(),
      fa1_(),
      fa2_()
  {
  }
  
  ~RactIP()
  {
  }

  RactIP& parse_options(int& argc, char**& argv);
  int run();
  float solve(const std::string& s1, const std::string& s2,
              std::string& r1, std::string& r2);

  static void calculate_energy(const std::string s1, const std::string& s2,
                               const std::string r1, const std::string& r2,
                               float& e1, float& e2, float& e3);

private:
  void contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const;
  void contraduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const;
  void rnafold(const std::string& seq, VF& bp, VI& offset) const;
  void rnafold(const std::string& seq, VF& bp, VI& offset, VVF& up, uint max_w) const;
  void rnaduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const;
  void load_from_rip(const char* filename,
                     const std::string& s1, const std::string& s2,
                     VF& bp1, VI& offset1, VF& bp2, VI& offset2, VVF& hp) const;

private:
  // options
  float alpha_;                // weight for the hybridization score
  float beta_;                 // weight for unpaired bases
  float th_hy_;                // threshold for the hybridization probability
  float th_ss_;                // threshold for the base-pairing probability
  float th_ac_;                // threshold for the accessible probability
  int max_w_;                  // maximum length of accessible regions
  int min_w_;                  // mimimum length of accessible regions
  int enable_zscore_;          // flag for calculating z-score
  int num_shuffling_;          // the number of shuffling for calculating z-score
  uint seed_;                  // seed for random()
  bool in_pk_;                 // allow internal pseudoknots or not
  bool use_contrafold_;        // use CONTRAfold model or not
  bool use_pf_duplex_;
  bool stacking_constraints_;
  bool show_energy_;
  bool allow_concat_;
  bool run_with_modena_;
  int n_th_;                   // the number of threads
  std::string rip_file_;
  std::string param_file_;
  std::string fa1_;
  std::string fa2_;
};

void
RactIP::
contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const
{
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  VF w = GetDefaultComplementaryValues<float>();
  bp.resize((seq.size()+1)*(seq.size()+2)/2);
  std::fill(bp.begin(), bp.end(), 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(0, bp, offset);

  const uint L=seq.size();
  up.resize(L, VF(1, 1.0));
  for (uint i=0; i!=L; ++i)
  {
    for (uint j=0; j<i; ++j)
      up[i][0] -= bp[offset[j+1]+(i+1)];
    for (uint j=i+1; j<L; ++j)
      up[i][0] -= bp[offset[i+1]+(j+1)];
    up[i][0] = std::max(0.0f, up[i][0]);
  }
}    

void
RactIP::
contraduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const
{
  hp.resize(seq1.size()+1, VF(seq2.size()+1));
  SStruct ss1("unknown", seq1), ss2("unknown", seq2);
  ParameterManager<float> pm;
  DuplexEngine<float> en(false);
  VF w = GetDefaultComplementaryValues<float>();
  VF ip;
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss1, ss2);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(th_hy_, ip);
  for (uint i=0; i!=seq1.size(); ++i)
    for (uint j=0; j!=seq2.size(); ++j)
      hp[i+1][j+1] = ip[en.GetOffset(i+1)+(j+1)];
}

void
RactIP::
rnafold(const std::string& seq, VF& bp, VI& offset) const
{
  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];
  Vienna::free_pf_arrays();
}

void
RactIP::
rnafold(const std::string& seq, VF& bp, VI& offset, VVF& up, uint max_w) const
{
  uint L=seq.size();
  bp.resize((L+1)*(L+2)/2);
  offset.resize(L+1);
  for (uint i=0; i<=L; ++i)
    offset[i] = i*((L+1)+(L+1)-i-1)/2;
#if 0
  std::string str(seq.size()+1, '.');
  float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &str[0]);
  float sfact = 1.07;
  float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
#else
  Vienna::pf_scale = -1;
#endif
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
      bp[offset[i+1]+(j+1)] = pr[iindx[i+1]-(j+1)];

  up.resize(L, VF(max_w));
  Vienna::pu_contrib* pu = Vienna::pf_unstru(const_cast<char*>(seq.c_str()), max_w);
  assert((int)L==pu->length);
  for (uint i=0; i!=L; ++i)
    for (uint j=0; j!=max_w; ++j)
      up[i][j]=pu->H[i+1][j] + pu->I[i+1][j] + pu->M[i+1][j] + pu->E[i+1][j];
#ifdef HAVE_VIENNA20
  Vienna::free_pu_contrib_struct(pu);
#else
  Vienna::free_pu_contrib(pu);
#endif
  Vienna::free_pf_arrays();
}

void
RactIP::
rnaduplex(const std::string& s1, const std::string& s2, VVF& hp) const
{
  if (use_pf_duplex_)
  {
    Vienna::pf_scale = -1;
    hp.resize(s1.size()+1, VF(s2.size()+1));
    Vienna::pf_duplex(s1.c_str(), s2.c_str());
    for (uint i=0; i!=s1.size(); ++i)
      for (uint j=0; j!=s2.size(); ++j)
        hp[i+1][j+1] = Vienna::pr_duplex[i+1][j+1];
    Vienna::free_pf_duplex();
  }
  else
  {
    hp.clear();
    hp.resize(s1.size()+1, VF(s2.size()+1, 0.0));
    std::string s=s1+s2;
    std::string c(s.size(), 'e');
    Vienna::pf_scale = -1;
    Vienna::cut_point = s1.size()+1;
    Vienna::co_pf_fold(const_cast<char*>(s.c_str()), const_cast<char*>(c.c_str()));
    pair_info* pi = NULL;
#ifdef HAVE_VIENNA20
    Vienna::assign_plist_from_pr(&pi, Vienna::export_co_bppm(), s.size(), th_hy_);
#else
    pi = Vienna::get_plist((pair_info*)malloc(sizeof(*pi)*s.size()), s.size(), th_hy_);
#endif
    for (uint k=0; pi[k].i!=0; ++k)
      if (pi[k].i<Vienna::cut_point && pi[k].j>=Vienna::cut_point && pi[k].p>th_hy_)
        hp[pi[k].i][pi[k].j-Vienna::cut_point+1]=pi[k].p;
    if (pi) free(pi);
    Vienna::free_co_pf_arrays();
    Vienna::cut_point = -1;
  }
}

void
RactIP::
load_from_rip(const char* filename,
              const std::string& s1, const std::string& s2,
              VF& bp1, VI& offset1, VF& bp2, VI& offset2, VVF& hp) const
{
  enum { NONE, TABLE_R, TABLE_S, TABLE_I };
  uint st=NONE;

  uint L1=s1.size();
  bp1.resize((L1+1)*(L1+2)/2);
  offset1.resize(L1+1);
  for (uint i=0; i<=L1; ++i)
    offset1[i] = i*((L1+1)+(L1+1)-i-1)/2;

  uint L2=s2.size();
  bp2.resize((L2+1)*(L2+2)/2);
  offset2.resize(L2+1);
  for (uint i=0; i<=L2; ++i)
    offset2[i] = i*((L2+1)+(L2+1)-i-1)/2;

  hp.resize(L1+1, VF(L2+1));

  std::ifstream is(filename);
  std::string l;
  while (std::getline(is, l))
  {
    if (strncmp(l.c_str(), "Table R:", 8)==0) st=TABLE_R;
    else if (strncmp(l.c_str(), "Table S:", 8)==0) st=TABLE_S;
    else if (strncmp(l.c_str(), "Table I:", 8)==0) st=TABLE_I;
    else if (st!=NONE && l[0]>='0' && l[0]<='9')
    {
      int i,j;
      double p;
      std::istringstream s(l.c_str());
      s >> i >> j >> p;
      switch (st)
      {
        case TABLE_R:
          bp1[offset1[i]+j] = p;
          break;
        case TABLE_S:
          bp2[offset2[L2-j+1]+L2-i+1] = p;
          break;
        case TABLE_I:
          hp[i][L2-j+1] = p;
          break;
        default:
          break;
      }
    }
    else st=NONE;
  }
}

float
RactIP::
solve(const std::string& s1, const std::string& s2, std::string& r1, std::string& r2)
{
  IP ip(IP::MAX, n_th_);
  VF bp1, bp2;
  VI offset1, offset2;
  VVF hp;
  VVF up1, up2;
  bool enable_accessibility = min_w_>1 && max_w_>=min_w_;

  // calculate posterior probability matrices
  if (!rip_file_.empty())
  {
    load_from_rip(rip_file_.c_str(), s1, s2, bp1, offset1, bp2, offset2, hp);
  }
  else if (use_contrafold_)
  {
    contrafold(s1, bp1, offset1, up1);
    contrafold(s2, bp2, offset2, up2);
    //contraduplex(s1, s2, hp);
    rnaduplex(s1, s2, hp);
  }
  else
  {
#if 0
    if (enable_accessibility)
    {
      rnafold(s1, bp1, offset1, up1, max_w_);
      rnafold(s2, bp2, offset2, up2, max_w_);
    }
    else
    {
      rnafold(s1, bp1, offset1);
      rnafold(s2, bp2, offset2);
    }
#else
    rnafold(s1, bp1, offset1, up1, std::max(1, max_w_));
    rnafold(s2, bp2, offset2, up2, std::max(1, max_w_));
#endif
    rnaduplex(s1, s2, hp);
  }
  
  // make objective variables with their weights
  VVI x(s1.size(), VI(s1.size(), -1));
  VVI xx(s1.size());
  for (uint j=1; j!=s1.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      const float& p=bp1[offset1[i+1]+(j+1)];
      if (p>th_ss_)
      {
        x[i][j] = ip.make_variable(p);
        xx[i].push_back(j);
      }
    }
  }

  VVI y(s2.size(), VI(s2.size(), -1));
  VVI yy(s2.size());
  for (uint j=1; j!=s2.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      const float& p=bp2[offset2[i+1]+(j+1)];
      if (p>th_ss_)
      {
        y[i][j] = ip.make_variable(p);
        yy[i].push_back(j);
      }
    }
  }

  VVI z(s1.size(), VI(s2.size(), -1));
  std::vector< std::vector<int> > zz(s1.size());
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      const float& p=hp[i+1][j+1];
      if (p>th_hy_ && (min_w_==1 && up1[i][0]>th_ac_ && up2[j][0]>th_ac_ || min_w_!=1))
      {
        z[i][j] = ip.make_variable(alpha_*p+beta_*(up1[i][0]+up2[j][0]));
        zz[i].push_back(j);
      }
    }
  }

  VI v, w;
  std::vector< std::pair<uint,uint> > vv, ww;
  if (enable_accessibility)
  {
    for (uint i=0; i!=up1.size(); ++i)
      for (uint j=min_w_-1; j<up1[i].size(); ++j)
        if (up1[i][j]>th_ac_)
        {
          v.push_back(ip.make_variable(0.0));
          vv.push_back(std::make_pair(i,i+j));
        }

    for (uint i=0; i!=up2.size(); ++i)
      for (uint j=min_w_-1; j<up2[i].size(); ++j)
        if (up2[i][j]>th_ac_)
        {
          w.push_back(ip.make_variable(0.0));
          ww.push_back(std::make_pair(i,i+j));
        }
  }

  ip.update();

  // constraint 1: each a_i is paired with at most one base
  for (uint i=0; i!=s1.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    if (enable_accessibility)
    {
      int row_ac = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j!=v.size(); ++j)
        if (vv[j].first<=i && i<=vv[j].second)
        {
          ip.add_constraint(row, v[j], 1);
          ip.add_constraint(row_ac, v[j], 1);
        }
      for (uint j=0; j!=s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row_ac, z[i][j], -1);
    }
    else
    {
      for (uint j=0; j!=s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }

    for (uint j=0; j<i; ++j)
      if (x[j][i]>=0)
        ip.add_constraint(row, x[j][i], 1);
    for (uint j=i+1; j<s1.size(); ++j)
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  if (enable_accessibility && !allow_concat_)
  {
    for (uint i=0; i!=v.size(); ++i)
      for (uint j=i+1; j!=v.size(); ++j)
        if (vv[i].second+1==vv[j].first || vv[j].second+1==vv[i].first)
        {
          int row = ip.make_constraint(IP::UP, 0, 1);
          ip.add_constraint(row, v[i], 1);
          ip.add_constraint(row, v[j], 1);
        }
  }

  // constraint 2: each b_i is paired with at most one
  for (uint i=0; i!=s2.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    if (enable_accessibility)
    {
      int row_ac = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j!=w.size(); ++j)
        if (ww[j].first<=i && i<=ww[j].second)
        {
          ip.add_constraint(row, w[j], 1);
          ip.add_constraint(row_ac, w[j], 1);
        }
      for (uint j=0; j!=s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row_ac, z[j][i], -1);
    }
    else
    {
      for (uint j=0; j!=s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row, z[j][i], 1);
    }

    for (uint j=0; j<i; ++j)
      if (y[j][i]>=0)
        ip.add_constraint(row, y[j][i], 1);
    for (uint j=i+1; j<s2.size(); ++j)
      if (y[i][j]>=0)
        ip.add_constraint(row, y[i][j], 1);
  }

  if (enable_accessibility && !allow_concat_)
  {
    for (uint i=0; i!=w.size(); ++i)
      for (uint j=i+1; j!=w.size(); ++j)
        if (ww[i].second+1==ww[j].first || ww[j].second+1==ww[i].first)
        {
          int row = ip.make_constraint(IP::UP, 0, 1);
          ip.add_constraint(row, w[i], 1);
          ip.add_constraint(row, w[j], 1);
        }
  }

  // constraint 3: disallow external pseudoknots
  for (uint i=0; i<zz.size(); ++i)
  {
    for (uint k=i+1; k<zz.size(); ++k)
    {
      for (uint p=0; p<zz[i].size(); ++p)
      { 
        uint j=zz[i][p];
        for (uint q=0; q<zz[k].size(); ++q)
        {
          uint l=zz[k][q];
          if (j<l)
          {
            int row = ip.make_constraint(IP::UP, 0, 1);
            ip.add_constraint(row, z[i][j], 1);
            ip.add_constraint(row, z[k][l], 1);
          }
        }
      }
    }
  }

  if (in_pk_)
  {
    // constraint 4: disallow internal pseudoknots in a
    for (uint i=0; i<xx.size(); ++i)
    {
      for (uint p=0; p<xx[i].size(); ++p)
      {
        uint j=xx[i][p];
        for (uint k=i+1; k<j; ++k)
        {
          for (uint q=0; q<xx[k].size(); ++q)
          {
            uint l=xx[k][q];
            if (j<l)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, x[i][j], 1);
              ip.add_constraint(row, x[k][l], 1);
            }
          }
        }
      }
    }

    // constraint 5: disallow internal pseudoknots in b
    for (uint i=0; i<yy.size(); ++i)
    {
      for (uint p=0; p<yy[i].size(); ++p)
      {
        uint j=yy[i][p];
        for (uint k=i+1; k<j; ++k)
        {
          for (uint q=0; q<yy[k].size(); ++q)
          {
            uint l=yy[k][q];
            if (j<l)
            {
              int row = ip.make_constraint(IP::UP, 0, 1);
              ip.add_constraint(row, y[i][j], 1);
              ip.add_constraint(row, y[k][l], 1);
            }
          }
        }
      }
    }
  }

  if (stacking_constraints_)
  {
    // upstream of s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<i; ++j)
        if (x[j][i]>=0)
          ip.add_constraint(row, x[j][i], -1);
      if (i>0)
        for (uint j=0; j<i-1; ++j)
          if (x[j][i-1]>=0)
            ip.add_constraint(row, x[j][i-1], 1);
      if (i+1<s1.size())
        for (uint j=0; j<i+1; ++j)
          if (x[j][i+1]>=0)
            ip.add_constraint(row, x[j][i+1], 1);
    }

    // downstream of s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=i+1; j<s1.size(); ++j)
        if (x[i][j]>=0)
          ip.add_constraint(row, x[i][j], -1);
      if (i>0)
        for (uint j=i; j<s1.size(); ++j)
          if (x[i-1][j]>=0)
            ip.add_constraint(row, x[i-1][j], 1);
      if (i+1<s1.size())
        for (uint j=i+2; j<s1.size(); ++j)
          if (x[i+1][j]>=0)
            ip.add_constraint(row, x[i+1][j], 1);
    }

    // upstream of s2
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<i; ++j)
        if (y[j][i]>=0)
          ip.add_constraint(row, y[j][i], -1);
      if (i>0)
        for (uint j=0; j<i-1; ++j)
          if (y[j][i-1]>=0)
            ip.add_constraint(row, y[j][i-1], 1);
      if (i+1<s2.size())
        for (uint j=0; j<i+1; ++j)
          if (y[j][i+1]>=0)
            ip.add_constraint(row, y[j][i+1], 1);
    }

    // downstream of s1
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=i+1; j<s2.size(); ++j)
        if (y[i][j]>=0)
          ip.add_constraint(row, y[i][j], -1);
      if (i>0)
        for (uint j=i; j<s2.size(); ++j)
          if (y[i-1][j]>=0)
            ip.add_constraint(row, y[i-1][j], 1);
      if (i+1<s2.size())
        for (uint j=i+2; j<s2.size(); ++j)
          if (y[i+1][j]>=0)
            ip.add_constraint(row, y[i+1][j], 1);
    }

    // for s2
    for (uint i=0; i<s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<s1.size(); ++j)
        if (z[j][i]>=0)
          ip.add_constraint(row, z[j][i], -1);
      if (i>0)
        for (uint j=0; j<s1.size(); ++j)
          if (z[j][i-1]>=0)
            ip.add_constraint(row, z[j][i-1], 1);
      if (i+1<s2.size())
        for (uint j=0; j<s1.size(); ++j)
          if (z[j][i+1]>=0)
            ip.add_constraint(row, z[j][i+1], 1);
    }

    // for s1
    for (uint i=0; i<s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<s2.size(); ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], -1);
      if (i>0)
        for (uint j=0; j<s2.size(); ++j)
          if (z[i-1][j]>=0)
            ip.add_constraint(row, z[i-1][j], 1);
      if (i+1<s1.size())
        for (uint j=0; j<s2.size(); ++j)
          if (z[i+1][j]>=0)
            ip.add_constraint(row, z[i+1][j], 1);
    }
  }

  // execute optimization
  float ea = ip.solve();

  // build the resultant structure
  r1.resize(s1.size());
  r2.resize(s2.size());
  std::fill(r1.begin(), r1.end(), '.');
  std::fill(r2.begin(), r2.end(), '.');
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      if (z[i][j]>=0 && ip.get_value(z[i][j])>0.5)
      {
        r1[i]='['; r2[j]=']';
      }
    }
  }
  if (in_pk_)
  {
    for (uint i=0; i<s1.size(); ++i)
    {
      for (uint j=i+1; j<s1.size(); ++j)
      {
        if (x[i][j]>=0 && ip.get_value(x[i][j])>0.5)
        {
          assert(r1[i]=='.'); assert(r1[j]=='.');
          r1[i]='('; r1[j]=')';
        }
      }
    }
    for (uint i=0; i<s2.size(); ++i)
    {
      for (uint j=i+1; j<s2.size(); ++j)
      {
        if (y[i][j]>=0 && ip.get_value(y[i][j])>0.5)
        {
          assert(r2[i]=='.'); assert(r2[j]=='.');
          r2[i]='('; r2[j]=')';
        }
      }
    }
  }

#if 0
  std::cout << "v: ";
  for (uint i=0; i!=v.size(); ++i)
  {
    if (ip.get_value(v[i])>0.5)
      std::cout << "(" << vv[i].first << "," << vv[i].second << ","
                << up1[vv[i].first][vv[i].second-vv[i].first]<< "), ";
  }
  std::cout << std::endl;

  std::cout << "w: ";
  for (uint i=0; i!=w.size(); ++i)
  {
    if (ip.get_value(w[i])>0.5)
      std::cout << "(" << ww[i].first << "," << ww[i].second << ","
                << up2[ww[i].first][ww[i].second-ww[i].first]<< "), ";
  }
  std::cout << std::endl;
#endif

  return ea;
}

RactIP&
RactIP::
parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);

  alpha_ = args_info.alpha_arg;
  beta_ = args_info.beta_arg;
  th_ss_ = args_info.fold_th_arg;
  th_hy_ = args_info.hybridize_th_arg;
  th_ac_ = args_info.acc_th_arg;
  max_w_ = args_info.max_w_arg;
  min_w_ = args_info.min_w_arg;
  enable_zscore_ = args_info.zscore_arg;
  num_shuffling_ = args_info.num_shuffling_arg;
  seed_ = args_info.seed_arg;
  in_pk_ = args_info.no_pk_flag==0;
  use_contrafold_ = args_info.mccaskill_flag==0;
  //use_pf_duplex_ = args_info.pf_duplex_flag;
  stacking_constraints_ = args_info.allow_isolated_flag==0;
  //allow_concat_ = args_info.allow_concat_flag;
  //run_with_modena_ = args_info.modena_flag;
  n_th_ = 1; // args_info.n_th_arg;
  if (args_info.rip_given) rip_file_ = args_info.rip_arg;
  show_energy_ = args_info.show_energy_flag==1;
  if (args_info.param_file_given) param_file_ = args_info.param_file_arg;

  if (args_info.inputs_num==0 ||
      (min_w_!=0 && max_w_!=0 && (min_w_>max_w_ || use_contrafold_)))
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }
  if (args_info.inputs_num>=1)
    fa1_ = args_info.inputs[0];
  if (args_info.inputs_num>=2)
    fa2_ = args_info.inputs[1];

  cmdline_parser_free(&args_info);
  return *this;
}

// static
void
RactIP::
calculate_energy(const std::string s1, const std::string& s2,
                 const std::string r1, const std::string& r2,
                 float& e1, float& e2, float& e3)
{
#ifdef HAVE_VIENNA20
  e1=Vienna::energy_of_structure(s1.c_str(), r1.c_str(), -1);
  e2=Vienna::energy_of_structure(s2.c_str(), r2.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  e1=Vienna::energy_of_struct(s1.c_str(), r1.c_str());
  e2=Vienna::energy_of_struct(s2.c_str(), r2.c_str());
#endif
  std::string ss(s1+"NNN"+s2);

  std::string r1_temp(r1);
  for (std::string::iterator x=r1_temp.begin(); x!=r1_temp.end(); ++x)
  {
    switch (*x)
    {
      case '(': case ')': *x='.'; break;
      case '[': *x='('; break;
      default: break;
    }
  }

  std::string r2_temp(r2);
  for (std::string::iterator x=r2_temp.begin(); x!=r2_temp.end(); ++x)
  {
    switch (*x)
    {
      case '(': case ')': *x='.'; break;
      case ']': *x=')'; break;
      default: break;
    }
  }

  std::string rr(r1_temp+"..."+r2_temp);
  //std::cout << ss << std::endl << rr << std::endl;
#ifdef HAVE_VIENNA20
  e3=Vienna::energy_of_structure(ss.c_str(), rr.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  e3=Vienna::energy_of_struct(ss.c_str(), rr.c_str());
#endif
}

int
RactIP::
run()
{
  // set the energy parameters
  copy_boltzmann_parameters();
  if (!param_file_.empty())
    Vienna::read_parameter_file(param_file_.c_str());
  
  // read sequences
  Fasta fa1, fa2;
  if (!fa1_.empty() && !fa2_.empty())
  {
    std::list<Fasta> l1, l2;
    if (Fasta::load(l1, fa1_.c_str())==0)
      throw (fa1_+": Format error").c_str();
    if (Fasta::load(l2, fa2_.c_str())==0)
      throw (fa2_+": Format error").c_str();
    fa1=l1.front();
    fa2=l2.front();
  }
  else if (!fa1_.empty())
  {
    std::list<Fasta> l1;
    if (Fasta::load(l1, fa1_.c_str())<2)
      throw (fa1_+": Format error").c_str();
    std::list<Fasta>::const_iterator x=l1.begin();
    fa1=*(x++);
    fa2=*(x++);
  }
  else { throw "unreachable"; }

  // predict the interation
  std::string r1, r2;
  float ea = solve(fa1.seq(), fa2.seq(), r1, r2);

  // display the result
  std::cout << ">" << fa1.name() << std::endl
            << fa1.seq() << std::endl << r1 << std::endl
            << ">" << fa2.name() << std::endl
            << fa2.seq() << std::endl << r2 << std::endl;

  // show energy of the joint structure
  if (show_energy_ || enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
  {
    float e1, e2, e3;
    calculate_energy(fa1.seq(), fa2.seq(), r1, r2, e1, e2, e3);
    if (show_energy_)
      std::cout << "(E: S1=" << e1 << ", "
                << "S2=" << e2 << ", "
                << "H=" << e3 << ", "
                << "JS=" << e1+e2+e3 << ")" << std::endl;

    if (enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
    {
      float sum=0.0, sum2=0.0;
      std::string s1(fa1.seq());
      std::string s2(fa2.seq());
      if (seed_==0)
      {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	seed_ = (unsigned long) tv.tv_sec;
      }
      srandom(seed_);
      uShuffle::set_randfunc((uShuffle::randfunc_t) random);
      for (uint i=0; i!=num_shuffling_; ++i)
      {
        if (enable_zscore_==1 || enable_zscore_==12)
          uShuffle::shuffle(fa1.seq().c_str(), &s1[0], fa1.seq().size(), 2);
        if (enable_zscore_==2 || enable_zscore_==12)
          uShuffle::shuffle(fa2.seq().c_str(), &s2[0], fa2.seq().size(), 2);
        solve(s1, s2, r1, r2);
        float ee1, ee2, ee3;
        calculate_energy(s1, s2, r1, r2, ee1, ee2, ee3);
        float ee=ee1+ee2+ee3;
#if 0
        std::cout << s1 << std::endl << r1 << std::endl
                  << s2 << std::endl << r2 << std::endl;
        std::cout << ee1 << " + " << ee2 << " + " << ee3 << " = " << ee << std::endl;
#endif
        sum += ee;
        sum2 += ee*ee;
      }
      float m=sum/num_shuffling_;
      float v=sum2/num_shuffling_-m*m;
      if (v<0.0) v=0.0;
#if 0
      std::cout << m << " " << v << std::endl;
#endif
      std::cout << "z-score: " << (e1+e2+e3-m)/sqrt(v) << std::endl;
    }
  }

  return 0;
}

int
main(int argc, char* argv[])
{
  try
  {
    RactIP ractip;
    return ractip.parse_options(argc, argv).run();
  }
  catch (const char* msg)
  {
    std::cout << msg << std::endl;
  }
  catch (std::logic_error err)
  {
    std::cout << err.what() << std::endl;
  }
  return 1;
}
