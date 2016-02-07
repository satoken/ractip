/*
 * $Id$
 * 
 * Copyright (C) 2010-2013 Kengo Sato
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
#include <stack>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstring>
#include "fa.h"
#include "ip.h"

#ifdef USE_CONTRAFOLD
#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/DuplexEngine.hpp"
#include "contrafold/Defaults.ipp"
#endif

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

#ifndef FLT_OR_DBL
typedef Vienna::FLT_OR_DBL FLT_OR_DBL;
#endif

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
      acc_max_(false),
      acc_max_ss_(false),
      acc_num_(0),
      max_w_(0),
      min_w_(0),
      enable_zscore_(0),
      num_shuffling_(0),
      seed_(0),
      in_pk_(true),
      use_contrafold_(false),
      use_pf_duplex_(false),
      use_constraint_(false),
      force_constraint_(false),
      stacking_constraints_(true),
      show_energy_(false),
      run_with_modena_(false),
      n_th_(1),
      rip_file_(),
      param_file_(),
      use_bl_param_(true),
      fa1_(),
      fa2_()
  {
  }
  
  ~RactIP()
  {
  }

  RactIP& parse_options(int& argc, char**& argv);
  int run();
  float solve(const Fasta& fa1, const Fasta& fa2,
              std::string& r1, std::string& r2,
              float* e1=NULL, float* e2=NULL, float* e3=NULL);
  float solve_ss(const std::string& s, const VF& bp, const VI& offset,
                 const std::vector<bool>& u, std::string& r, float* e=NULL);
  float solve_ss(const std::string& s, const VF& bp, const VI& offset,
                 std::string& r, float* e=NULL);

  static float energy_of_duplex(const std::string s1, const std::string& s2,
                                const std::string r1, const std::string& r2);

private:
#ifdef USE_CONTRAFOLD
  void contrafold(const std::string& seq, VF& bp, VI& offset, VVF& up) const;
  void contraduplex(const std::string& seq1, const std::string& seq2, VVF& hp) const;
#endif
  void rnafold(const Fasta& fa, VF& bp, VI& offset) const;
  void rnafold(const Fasta& fa, VF& bp, VI& offset, VVF& up, uint max_w) const;
  void rnaduplex(const Fasta& fa1, const Fasta& fa2, VVF& hp) const;
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
  bool acc_max_;               // optimize for accessibility instead of internal secondary structures
  bool acc_max_ss_;            // additional prediction of interanal secondary structures
  int acc_num_;                // the number of accessible regions
  int max_w_;                  // maximum length of accessible regions
  int min_w_;                  // mimimum length of accessible regions
  int enable_zscore_;          // flag for calculating z-score
  uint num_shuffling_;          // the number of shuffling for calculating z-score
  uint seed_;                  // seed for random()
  bool in_pk_;                 // allow internal pseudoknots or not
  bool use_contrafold_;        // use CONTRAfold model or not
  bool use_pf_duplex_;
  bool use_constraint_;
  bool force_constraint_;
  bool stacking_constraints_;
  bool show_energy_;
  bool run_with_modena_;
  int n_th_;                   // the number of threads
  std::string rip_file_;
  std::string param_file_;
  bool use_bl_param_;
  std::string fa1_;
  std::string fa2_;

  VF bp1_;
  VF bp2_;
  VI offset1_;
  VI offset2_;
  VVF hp_;
  VVF up1_;
  VVF up2_;
};

#ifdef USE_CONTRAFOLD
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
#endif

void
RactIP::
rnafold(const Fasta& fa, VF& bp, VI& offset) const
{
  const std::string& seq=fa.seq();
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
  if (use_constraint_) 
  {
    int fc_back=Vienna::fold_constrained;
    Vienna::fold_constrained=1;
    const std::string& str=fa.str();
    std::string c(seq.size(), '.');
    for (uint i=0; i!=str.size() && i!=seq.size(); ++i)
    {
      switch (str[i])
      {
        case '[':
        case ']':
        case 'e':
          c[i]='x'; break;
        default: 
          c[i]=str[i]; break;
      }
    }
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), const_cast<char*>(&c[0]));
    Vienna::fold_constrained=fc_back;
  }
  else
  {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  }
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
rnafold(const Fasta& fa, VF& bp, VI& offset, VVF& up, uint max_w) const
{
  const std::string& seq=fa.seq();
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
  if (use_constraint_) 
  {
    int fc_bak=Vienna::fold_constrained;
    Vienna::fold_constrained=1;
    const std::string& str=fa.str();
    std::string c(seq.size(), '.');
    for (uint i=0; i!=str.size() && i!=seq.size(); ++i)
    {
      switch (str[i])
      {
        case '[':
        case ']':
        case 'e':
          c[i]='x'; break;
        default: 
          c[i]=str[i]; break;
      }
    }
#if 0
    std::cout << c << std::endl;
#endif
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), const_cast<char*>(&c[0]));
    Vienna::fold_constrained=fc_bak;
  }
  else
  {
    Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  }
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

  // accessibility
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
rnaduplex(const Fasta& fa1, const Fasta& fa2, VVF& hp) const
{
  const std::string& s1=fa1.seq();
  const std::string& s2=fa2.seq();
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
    int fc_bak=Vienna::fold_constrained;
    Vienna::fold_constrained=1;
    hp.clear();
    hp.resize(s1.size()+1, VF(s2.size()+1, 0.0));
    std::string s=s1+s2;
    std::string c(s.size(), '.');
    if (use_constraint_)
    {
      const std::string str1=fa1.str();
      const std::string str2=fa2.str();
      for (uint i=0; i!=str1.size() && i!=s1.size(); ++i)
      {
        switch (str1[i])
        {
          case '[': c[i]='('; break;
          case '(':
          case ')':
          case 'l':
          case 'x': c[i]='x'; break;
          default: break;
        }
      }
      for (uint i=0; i!=str2.size() && i!=s2.size(); ++i)
      {
        switch (str2[i])
        {
          case ']': c[s1.size()+i]=')'; break;
          case '(':
          case ')':
          case 'l':
          case 'x': c[s1.size()+i]='x'; break;
          default: break;
        }
      }
#if 0
      std::cout << c.substr(0, s1.size()) << std::endl
                << c.substr(s1.size()) << std::endl;
#endif
    }

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
    Vienna::fold_constrained=fc_bak;
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
solve(const Fasta& fa1, const Fasta& fa2, std::string& r1, std::string& r2,
      float* e1 /* =NULL */, float* e2 /* =NULL */, float* e3 /* =NULL */)
{
  const std::string& s1=fa1.seq();
  const std::string& s2=fa2.seq();
  const std::string& str1=fa1.str();
  const std::string& str2=fa2.str();
  IP ip(IP::MAX, n_th_);
  bool enable_accessibility = min_w_>1 && max_w_>=min_w_;
  bool enable_structure_s1 = !acc_max_;
  bool enable_structure_s2 = !acc_max_;

  // calculate posterior probability matrices
  if (!rip_file_.empty())
  {
    load_from_rip(rip_file_.c_str(), s1, s2, bp1_, offset1_, bp2_, offset2_, hp_);
  }
#ifdef USE_CONTRAFOLD
  else if (use_contrafold_)
  {
    contrafold(s1, bp1_, offset1_, up1_);
    contrafold(s2, bp2_, offset2_, up2_);
    //contraduplex(s1_, s2_, hp_);
    rnaduplex(fa1, fa2, hp_);
  }
#endif
  else
  {
    rnafold(fa1, bp1_, offset1_, up1_, std::max(1, max_w_));
    rnafold(fa2, bp2_, offset2_, up2_, std::max(1, max_w_));
    rnaduplex(fa1, fa2, hp_);
  }
  
  // make objective variables with their weights
  VVI x(s1.size(), VI(s1.size(), -1)); // internal base-pairs in s1
  VVI xx(s1.size());
  VI x_un(s1.size(), -1);       // internally-unpaired bases in s1
  if (enable_structure_s1)
  {
    for (uint j=1; j!=s1.size(); ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp1_[offset1_[i+1]+(j+1)];
        if (p>th_ss_)
        {
          x[i][j] = x[j][i] = ip.make_variable(p-th_ss_);
          xx[i].push_back(j);
        }
      }
    }
    for (uint i=0; i!=s1.size(); ++i)
      x_un[i] = ip.make_variable(0.0);
  }

  VVI y(s2.size(), VI(s2.size(), -1)); // internal base-pairs in s2
  VVI yy(s2.size());
  VI y_un(s2.size(), -1);       // internally-unpaired bases in s2
  if (enable_structure_s2)
  {
    for (uint j=1; j!=s2.size(); ++j)
    {
      for (uint i=j-1; i!=-1u; --i)
      {
        const float& p=bp2_[offset2_[i+1]+(j+1)];
        if (p>th_ss_)
        {
          y[i][j] = y[j][i] = ip.make_variable(p-th_ss_);
          yy[i].push_back(j);
        }
      }
    }
    for (uint i=0; i!=s2.size(); ++i)
      y_un[i] = ip.make_variable(0.0);
  }

  VVI z(s1.size(), VI(s2.size(), -1)); // external base-pairs
  VVI zz(s1.size());
  VI z_un1(s1.size(), -1);      // externally-unpaired bases in s1
  VI z_un2(s2.size(), -1);      // externally-unpaired bases in s2
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      const float& p=hp_[i+1][j+1];
      if (p>th_hy_)
      {
        z[i][j] = ip.make_variable(alpha_*(p-th_hy_));
        zz[i].push_back(j);
      }
    }
  }
  for (uint i=0; i!=s1.size(); ++i)
    z_un1[i] = ip.make_variable(0.0);
  for (uint i=0; i!=s2.size(); ++i)
    z_un2[i] = ip.make_variable(0.0);

  VI v;                         // accessible regions in s1
  std::vector< std::pair<uint,uint> > vv;
  VI v_st(s1.size(), -1);       // # of accessible regions that start at i in s1
  VI v_en(s1.size(), -1);       // # of accessible regions that end at i in s1
  if (enable_accessibility)
  {
    for (uint i=0; i!=up1_.size(); ++i)
      for (uint j=min_w_-1; j<up1_[i].size(); ++j)
        if (up1_[i][j]>th_ac_)
        {
          v.push_back(ip.make_variable(beta_*(up1_[i][j]-th_ac_)));
          vv.push_back(std::make_pair(i,i+j));
        }
  }
  for (uint i=0; i!=s1.size(); ++i)
  {
    v_st[i] = ip.make_variable(0.0);
    v_en[i] = ip.make_variable(0.0);
  }

  VI w;                         // accessible regions in s2
  std::vector< std::pair<uint,uint> > ww;
  VI w_st(s2.size(), -1);       // # of accessible regions that start at i in s2
  VI w_en(s2.size(), -1);       // # of accessible regions that end at i in s2
  if (enable_accessibility)
  {
    for (uint i=0; i!=up2_.size(); ++i)
      for (uint j=min_w_-1; j<up2_[i].size(); ++j)
        if (up2_[i][j]>th_ac_)
        {
          w.push_back(ip.make_variable(beta_*(up2_[i][j]-th_ac_)));
          ww.push_back(std::make_pair(i,i+j));
        }
  }
  for (uint i=0; i!=s2.size(); ++i)
  {
    w_st[i] = ip.make_variable(0.0);
    w_en[i] = ip.make_variable(0.0);
  }

  if (force_constraint_)
  {
    std::stack<int> st_x, st_y, st_z;
    for (uint i=0; i!=str1.size(); ++i)
    {
      switch (str1[i]) 
      {
        case '(': 
          st_x.push(i); 
          break;
        case ')':
          if (x[st_x.top()][i]<0)
          {
            const int j=st_x.top();
            const float& p=bp1_[offset1_[j+1]+(i+1)];
            x[i][j] = x[j][i] = ip.make_variable(p-th_ss_);
            xx[j].push_back(i);
          }
          st_x.pop();
          break;
        case '[':
          st_z.push(i);
          break;
        default:
          break;
      }
    }
    for (uint i=0; i!=str2.size(); ++i)
    {
      switch (str2[i]) 
      {
        case '(': 
          st_y.push(i); 
          break;
        case ')': 
          if (y[st_y.top()][i]<0)
          {
            const int j=st_y.top();
            const float& p=bp2_[offset2_[j+1]+(i+1)];
            y[i][j] = y[j][i] = ip.make_variable(p-th_ss_);
            yy[j].push_back(i);
          }
          st_y.pop();
          break;
        case ']':
          if (z[st_z.top()][i]<0)
          {
            const int j=st_z.top();
            const float& p=hp_[j+1][i+1];
            z[j][i] = ip.make_variable(alpha_*(p-th_hy_));
            zz[j].push_back(i);
          }
          st_z.pop(); 
          break;
        default:
          break;
      }
    }
  }

  ip.update();

  // constraints for helper variables
  if (enable_structure_s1)
  {
    for (uint i=0; i!=s1.size(); ++i)
    {
      // sum_j x[i][j] + x_un[i] = 1
      int row = ip.make_constraint(IP::FX, 1, 1);
      ip.add_constraint(row, x_un[i], 1);
      for (uint j=0; j!=s1.size(); ++j)
        if (x[i][j]>=0)
          ip.add_constraint(row, x[i][j], 1);
    }
  }

  for (uint i=0; i!=s1.size(); ++i)
  {
    // sum_j z[i][j] + z_un1[i] = 1
    int row = ip.make_constraint(IP::FX, 1, 1);
    ip.add_constraint(row, z_un1[i], 1);
    for (uint j=0; j!=s2.size(); ++j)
      if (z[i][j]>=0)
        ip.add_constraint(row, z[i][j], 1);
  }

  if (enable_structure_s2)
  {
    for (uint i=0; i!=s2.size(); ++i)
    {
      // sum_j y[i][j] + y_un[i] = 1
      int row = ip.make_constraint(IP::FX, 1, 1);
      ip.add_constraint(row, y_un[i], 1);
      for (uint j=0; j!=s2.size(); ++j)
        if (y[i][j]>=0)
          ip.add_constraint(row, y[i][j], 1);
    }
  }

  for (uint i=0; i!=s2.size(); ++i)
  {
    // sum_j z[j][i] + z_un1[i] = 1
    int row = ip.make_constraint(IP::FX, 1, 1);
    ip.add_constraint(row, z_un2[i], 1);
    for (uint j=0; j!=s1.size(); ++j)
      if (z[j][i]>=0)
        ip.add_constraint(row, z[j][i], 1);
  }

  if (enable_accessibility)
  {
    // sum_q v[p][q] = v_st[p]
    // sum_p v[p][q] = v_en[q]
    VI row_v_st(s1.size(), -1);
    VI row_v_en(s1.size(), -1);
    for (uint i=0; i!=s1.size(); ++i)
    {
      row_v_st[i] = ip.make_constraint(IP::FX, 0, 0);
      ip.add_constraint(row_v_st[i], v_st[i], -1);
      row_v_en[i] = ip.make_constraint(IP::FX, 0, 0);
      ip.add_constraint(row_v_en[i], v_en[i], -1);
    }
    for (uint i=0; i!=v.size(); ++i)
    {
      ip.add_constraint(row_v_st[vv[i].first], v[i], 1);
      ip.add_constraint(row_v_en[vv[i].second], v[i], 1);
    }

    // sum_q w[p][q] = w_st[p]
    // sum_p w[p][q] = w_en[q]
    VI row_w_st(s2.size(), -1);
    VI row_w_en(s2.size(), -1);
    for (uint i=0; i!=s2.size(); ++i)
    {
      row_w_st[i] = ip.make_constraint(IP::FX, 0, 0);
      ip.add_constraint(row_w_st[i], w_st[i], -1);
      row_w_en[i] = ip.make_constraint(IP::FX, 0, 0);
      ip.add_constraint(row_w_en[i], w_en[i], -1);
    }
    for (uint i=0; i!=w.size(); ++i)
    {
      ip.add_constraint(row_w_st[ww[i].first], w[i], 1);
      ip.add_constraint(row_w_en[ww[i].second], w[i], 1);
    }
  }


  if (!enable_accessibility)
  {
    if (enable_structure_s1)
    {
      // each base is paired with at most one base
      // x_un[i]+z_un1[i]>=0
      //  -> sum_j x[i][j] + sum_k z[i][k] <= 1
      for (uint i=0; i!=s1.size(); ++i)
      {
        int row = ip.make_constraint(IP::LO, 1, 0);
        ip.add_constraint(row, x_un[i], 1);
        ip.add_constraint(row, z_un1[i], 1);
      }
    }

    if (enable_structure_s2)
    {
      // each base is paired with at most one base
      // y_un[i]+y_un1[i]>=0
      //  -> sum_j y[i][j] + sum_k z[k][i] <= 1
      for (uint i=0; i!=s2.size(); ++i)
      {
        int row = ip.make_constraint(IP::LO, 1, 0);
        ip.add_constraint(row, y_un[i], 1);
        ip.add_constraint(row, z_un2[i], 1);
      }
    }
  }
  else
  {
    // each internal base pair x[i][j] is not accessible.
    // -x_un[i]+sum_{p<=i<=q} v[p][q] <=0
    //  -> sum_j x[i][j] + sum_{p<=i<=q} v[p][q] <=1
    if (enable_structure_s1)
    {
      VI row(s1.size());
      for (uint i=0; i!=s1.size(); ++i)
      {
        row[i] = ip.make_constraint(IP::UP, 0, 0);
        ip.add_constraint(row[i], x_un[i], -1);
      }
      for (uint j=0; j!=v.size(); ++j)
        for (uint i=vv[j].first; i<=vv[j].second; ++i)
          ip.add_constraint(row[i], v[j], 1);
    }

    // each external base pair z[i][j] is accessible in s1.
    // z_un1[i] +  sum_{p<=i<=q} v[p][q] >=1
    //  -> sum_{p<=i<=q} v[p][q] >= sum_k z[i][k]
    {
      VI row(s1.size(), 0);
      for (uint i=0; i!=s1.size(); ++i)
      {
        row[i] = ip.make_constraint(IP::LO, 1, 0);
        ip.add_constraint(row[i], z_un1[i], 1);
      }
      for (uint j=0; j!=v.size(); ++j)
        for (uint i=vv[j].first; i<=vv[j].second; ++i)
          ip.add_constraint(row[i], v[j], 1);
    }

    // each internal base pair x[i][j] is not accessible.
    // -y_un[i]+sum_{p<=i<=q} w[p][q] <=0
    //  -> sum_j y[i][j] + sum_{p<=i<=q} w[p][q] <=1
    if (enable_structure_s2)
    {
      VI row(s2.size());
      for (uint i=0; i!=s2.size(); ++i)
      {
        row[i] = ip.make_constraint(IP::UP, 0, 0);
        ip.add_constraint(row[i], y_un[i], -1);
      }
      for (uint j=0; j!=w.size(); ++j)
        for (uint i=ww[j].first; i<=ww[j].second; ++i)
          ip.add_constraint(row[i], w[j], 1);
    }

    // each external base pair z[i][j] is accessible in s2.
    // z_un2[i] +  sum_{p<=i<=q} w[p][q] >=1
    //  -> sum_{p<=i<=q} w[p][q] >= sum_k z[k][i]
    {
      VI row(s2.size());
      for (uint i=0; i!=s2.size(); ++i)
      {
        row[i] = ip.make_constraint(IP::LO, 1, 0);
        ip.add_constraint(row[i], z_un2[i], 1);
      }
      for (uint j=0; j!=w.size(); ++j)
        for (uint i=ww[j].first; i<=ww[j].second; ++i)
          ip.add_constraint(row[i], w[j], 1);
    }

    // each position is included in at most one region
    // sum_{p<=i<=q} v[p][q] <= 1
    {
      VI row(s1.size(), -1);
      for (uint i=0; i!=s1.size(); ++i)
        row[i] = ip.make_constraint(IP::UP, 0, 1);
      for (uint j=0; j!=v.size(); ++j)
        for (uint i=vv[j].first; i<=vv[j].second; ++i)
          ip.add_constraint(row[i], v[j], 1);
    }
    
    // no pairs of accessible resgions adjoin each other
    // v_en[i-1] + v_st[i] <= 1
    //  -> v[i][j] + v[j+1][k] <= 1
    for (uint i=1; i!=s1.size(); ++i)
    {
      int row = ip.make_constraint(IP::UP, 0, 1);
      ip.add_constraint(row, v_en[i-1], 1);
      ip.add_constraint(row, v_st[i], 1);
    }

    // each position is included in at most one region
    // sum_{p<=i<=q} w[p][q] <= 1
    {
      VI row(s2.size(), -1);
      for (uint i=0; i!=s2.size(); ++i)
        row[i] = ip.make_constraint(IP::UP, 0, 1);
      for (uint j=0; j!=w.size(); ++j)
        for (uint i=ww[j].first; i<=ww[j].second; ++i)
          ip.add_constraint(row[i], w[j], 1);
    }

    // no pairs of accessible resgions adjoin each other
    // w_en[i-1] + w_st[i] <= 1
    //  -> w[i][j] + w[j+1][k] <= 1
    for (uint i=1; i!=s2.size(); ++i)
    {
      int row = ip.make_constraint(IP::UP, 0, 1);
      ip.add_constraint(row, w_en[i-1], 1);
      ip.add_constraint(row, w_st[i], 1);
    }

    if (beta_>0.0)
    {
      // each accessible region v[p][q] contains at least one z[i][j]
      // sum_{p<=i<=q} z_un1[i] + v[p][q] <= p-q+1
      //  -> sum_{p<=i<=q} sum_k z[i][k] >= v[p][q]
      for (uint j=0; j!=v.size(); ++j)
      {
        int row = ip.make_constraint(IP::UP, 0, vv[j].second-vv[j].first+1);
        ip.add_constraint(row, v[j], 1);
        for (uint i=vv[j].first; i<=vv[j].second; ++i)
          ip.add_constraint(row, z_un1[i], 1);
      }

      // each accessible region w[p][q] contains at least one z[i][j]
      // sum_{p<=i<=q} z_un2[i] + w[p][q] <= p-q+1
      //  -> sum_{p<=i<=q} sum_k z[k][i] >= w[p][q]
      for (uint j=0; j!=w.size(); ++j)
      {
        int row = ip.make_constraint(IP::UP, 0, ww[j].second-ww[j].first+1);
        ip.add_constraint(row, w[j], 1);
        for (uint i=ww[j].first; i<=ww[j].second; ++i)
          ip.add_constraint(row, z_un2[i], 1);
      }
    }

#if 0
    // sum_{p,q} v[p][q] = sum_{p,q} w[p][q]
    {
      int row = ip.make_constraint(IP::FX, 0, 0);
      for (uint i=0; i!=v.size(); ++i)
        ip.add_constraint(row, v[i], 1);
      for (uint i=0; i!=w.size(); ++i)
        ip.add_constraint(row, w[i], -1);
    }
#endif
    if (acc_num_>0)
    {
      // sum_{p,q} v[p][q] <= acc_num_
      int row = ip.make_constraint(IP::UP, 0, acc_num_);
      for (uint i=0; i!=v.size(); ++i)
        ip.add_constraint(row, v[i], 1);

      // sum_{p,q} w[p][q] <= acc_num_
      row = ip.make_constraint(IP::UP, 0, acc_num_);
      for (uint i=0; i!=w.size(); ++i)
        ip.add_constraint(row, w[i], 1);
    }
  }

  if (enable_accessibility && acc_num_>0)
  {
    int row = ip.make_constraint(IP::UP, 0, acc_num_);
    for (uint i=0; i!=v.size(); ++i)
      ip.add_constraint(row, v[i], 1);

    row = ip.make_constraint(IP::UP, 0, acc_num_);
    for (uint i=0; i!=w.size(); ++i)
      ip.add_constraint(row, w[i], 1);
  }

  // disallow external pseudoknots
  for (uint i=0; i<zz.size(); ++i)
    for (uint k=i+1; k<zz.size(); ++k)
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

  if (in_pk_)
  {
    // disallow internal pseudoknots in a
    if (enable_structure_s1)
    {
      for (uint i=0; i<xx.size(); ++i)
        for (uint p=0; p<xx[i].size(); ++p)
        {
          uint j=xx[i][p];
          for (uint k=i+1; k<j; ++k)
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

    // disallow internal pseudoknots in b
    if (enable_structure_s2)
    {
      for (uint i=0; i<yy.size(); ++i)
        for (uint p=0; p<yy[i].size(); ++p)
        {
          uint j=yy[i][p];
          for (uint k=i+1; k<j; ++k)
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

  if (stacking_constraints_)
  {
    if (enable_structure_s1)
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
    }

    if (enable_structure_s2)
    {
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

  if (force_constraint_)
  {
    std::stack<int> st_x, st_y, st_z;
    for (uint i=0; i!=str1.size(); ++i)
    {
      switch (str1[i]) 
      {
        case '(': 
          st_x.push(i); 
          break;
        case ')':
          if (x[st_x.top()][i]>=0)
          {
            int row = ip.make_constraint(IP::FX, 1, 1);
            ip.add_constraint(row, x[st_x.top()][i], 1);
          }
          st_x.pop();
          break;
        case '[':
          st_z.push(i);
          break;
        default:
          break;
      }
    }
    for (uint i=0; i!=str2.size(); ++i)
    {
      switch (str2[i]) 
      {
        case '(': 
          st_y.push(i); 
          break;
        case ')': 
          if (y[st_y.top()][i]>=0)
          {
            int row = ip.make_constraint(IP::FX, 1, 1);
            ip.add_constraint(row, y[st_y.top()][i], 1);
          }
          st_y.pop();
          break;
        case ']':
          if (z[st_z.top()][i]>=0)
          {
            int row = ip.make_constraint(IP::FX, 1, 1);
            ip.add_constraint(row, z[st_z.top()][i], 1);
          }
          st_z.pop(); 
          break;
        default:
          break;
      }
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
    for (uint j=0; j!=s2.size(); ++j)
      if (z[i][j]>=0 && ip.get_value(z[i][j])>0.5)
      {
        r1[i]='['; r2[j]=']';
      }
  if (e3!=NULL)
    *e3 = energy_of_duplex(s1, s2, r1, r2);

  if (enable_structure_s1)
  {
    if (in_pk_)
      for (uint i=0; i<s1.size(); ++i)
        for (uint j=i+1; j<s1.size(); ++j)
          if (x[i][j]>=0 && ip.get_value(x[i][j])>0.5)
          {
            assert(r1[i]=='.'); assert(r1[j]=='.');
            r1[i]='('; r1[j]=')';
          }
    if (e1!=NULL)
    {
#ifdef HAVE_VIENNA20
      *e1=Vienna::energy_of_structure(s1.c_str(), r1.c_str(), -1);
#else
      Vienna::eos_debug = -1;
      *e1=Vienna::energy_of_struct(s1.c_str(), r1.c_str());
#endif
    }      
  }
  else 
  {
    if (acc_max_ss_)
    {
      std::vector<bool> u1(s1.size(), true);
      for (uint j=0; j!=v.size(); ++j)
        if (ip.get_value(v[j])>0.5)
          for (uint i=vv[j].first; i<=vv[j].second; ++i)
            u1[i]=false;
      ea += solve_ss(s1, bp1_, offset1_, u1, r1);
    }
    if (e1!=NULL)
    {
#ifdef GASCONST
      const double kT = (Vienna::temperature+K0)*GASCONST/1000.;
#else
      const double kT = (Vienna::temperature+273.15)*1.98717/1000.;
#endif
      *e1=0.0;
      for (uint j=0; j!=v.size(); ++j)
        if (ip.get_value(v[j])>0.5)
          *e1+=-log(up1_[vv[j].first][vv[j].second-vv[j].first])*kT;
    }
  }

  if (enable_structure_s2)
  {
    if (in_pk_)
      for (uint i=0; i<s2.size(); ++i)
        for (uint j=i+1; j<s2.size(); ++j)
          if (y[i][j]>=0 && ip.get_value(y[i][j])>0.5)
          {
            assert(r2[i]=='.'); assert(r2[j]=='.');
            r2[i]='('; r2[j]=')';
          }
    if (e2!=NULL)
    {
#ifdef HAVE_VIENNA20
      *e2=Vienna::energy_of_structure(s2.c_str(), r2.c_str(), -1);
#else
      Vienna::eos_debug = -1;
      *e2=Vienna::energy_of_struct(s2.c_str(), r2.c_str());
#endif
    }      
  }
  else 
  {
    if (acc_max_ss_)
    {
      std::vector<bool> u2(s2.size(), true);
      for (uint j=0; j!=w.size(); ++j)
        if (ip.get_value(w[j])>0.5)
          for (uint i=ww[j].first; i<=ww[j].second; ++i)
            u2[i]=false;
      ea += solve_ss(s2, bp2_, offset2_, u2, r2);
    }
    if (e2!=NULL)
    {
#ifdef GASCONST
      const double kT = (Vienna::temperature+K0)*GASCONST/1000.;
#else
      const double kT = (Vienna::temperature+273.15)*1.98717/1000.;
#endif
      *e2=0.0;
      for (uint j=0; j!=w.size(); ++j)
        if (ip.get_value(w[j])>0.5)
          *e2+=-log(up2_[ww[j].first][ww[j].second-ww[j].first])*kT;
    }
  }

#if 0
  std::cout << v.size() << " "  << w.size() << std::endl;
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

float
RactIP::
solve_ss(const std::string& s, const VF& bp, const VI& offset,
         std::string& r, float* e /*=NULL*/)
{
  std::vector<bool> u(s.size(), true);
  r.resize(s.size());
  std::fill(r.begin(), r.end(), '.');
  return solve_ss(s, bp, offset, u, r, e);
}

float
RactIP::
solve_ss(const std::string& s, const VF& bp, const VI& offset,
         const std::vector<bool>& u, std::string& r, float* e /*=NULL*/)
{
  IP ip(IP::MAX, n_th_);

  // make objective variables with their weights
  VVI x(s.size(), VI(s.size(), -1));
  VVI xx(s.size());
  for (uint j=1; j!=s.size(); ++j)
  {
    if (!u[j]) continue;
    for (uint i=j-1; i!=-1u; --i)
    {
      if (!u[i]) continue;
      const float& p=bp[offset[i+1]+(j+1)];
      if (p>th_ss_)
      {
        x[i][j] = x[j][i] = ip.make_variable(p-th_ss_);
        xx[i].push_back(j);
      }
    }
  }

  ip.update();

  // each a_i is paired with at most one base
  for (uint i=0; i!=s.size(); ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j=0; j<s.size(); ++j)
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  if (stacking_constraints_)
  {
    // upstream of s1
    for (uint i=0; i<s.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=0; j<i; ++j)
        if (x[j][i]>=0)
          ip.add_constraint(row, x[j][i], -1);
      if (i>0)
        for (uint j=0; j<i-1; ++j)
          if (x[j][i-1]>=0)
            ip.add_constraint(row, x[j][i-1], 1);
      if (i+1<s.size())
        for (uint j=0; j<i+1; ++j)
          if (x[j][i+1]>=0)
            ip.add_constraint(row, x[j][i+1], 1);
    }

    // downstream of s1
    for (uint i=0; i<s.size(); ++i)
    {
      int row = ip.make_constraint(IP::LO, 0, 0);
      for (uint j=i+1; j<s.size(); ++j)
        if (x[i][j]>=0)
          ip.add_constraint(row, x[i][j], -1);
      if (i>0)
        for (uint j=i; j<s.size(); ++j)
          if (x[i-1][j]>=0)
            ip.add_constraint(row, x[i-1][j], 1);
      if (i+1<s.size())
        for (uint j=i+2; j<s.size(); ++j)
          if (x[i+1][j]>=0)
            ip.add_constraint(row, x[i+1][j], 1);
    }
  }

  // execute optimization
  float ea = ip.solve();
  
  for (uint i=0; i<s.size(); ++i)
  {
    for (uint j=i+1; j<s.size(); ++j)
    {
      if (x[i][j]>=0 && ip.get_value(x[i][j])>0.5)
      {
        assert(r[i]=='.'); assert(r[j]=='.');
        r[i]='('; r[j]=')';
      }
    }
  }

  if (e!=NULL)
  {
#ifdef HAVE_VIENNA20
      *e=Vienna::energy_of_structure(s.c_str(), r.c_str(), -1);
#else
      Vienna::eos_debug = -1;
      *e=Vienna::energy_of_struct(s.c_str(), r.c_str());
#endif
  }

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
  acc_max_ = args_info.acc_max_flag==1;
  acc_max_ss_ = args_info.acc_max_ss_flag==1;
  acc_num_ = args_info.acc_num_arg;
  max_w_ = args_info.max_w_arg;
  min_w_ = args_info.min_w_arg;
  enable_zscore_ = args_info.zscore_arg;
  num_shuffling_ = args_info.num_shuffling_arg;
  seed_ = args_info.seed_arg;
  in_pk_ = args_info.no_pk_flag==0;
  use_contrafold_ = args_info.contrafold_flag==1;
  use_pf_duplex_ = args_info.duplex_flag==1;
  use_constraint_ = args_info.use_constraint_flag==1;
  force_constraint_ = args_info.force_constraint_flag==1;
  stacking_constraints_ = args_info.allow_isolated_flag==0;
  //run_with_modena_ = args_info.modena_flag;
  n_th_ = 1; // args_info.n_th_arg;
  if (args_info.rip_given) rip_file_ = args_info.rip_arg;
  show_energy_ = args_info.show_energy_flag==1;
  if (args_info.param_file_given) param_file_ = args_info.param_file_arg;
  use_bl_param_ = args_info.no_bl_flag==0;

#ifndef USE_CONTRAFOLD
  if (use_contrafold_)
  {
    std::cout << "ERROR: This binary is not configured for using CONTRAfold." 
              << std::endl << std::endl;
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }
#endif
  
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
float
RactIP::
energy_of_duplex(const std::string s1, const std::string& s2,
                 const std::string r1, const std::string& r2)
{
  float e;
  std::string ss(s1+s2);
  std::string rr(r1+r2);
  for (std::string::iterator x=rr.begin(); x!=rr.end(); ++x)
  {
    switch (*x)
    {
      case '(': case ')': *x='.'; break;
      case '[': *x='('; break;
      case ']': *x=')'; break;
      default: break;
    }
  }
  //std::cout << ss << std::endl << rr << std::endl;

  Vienna::cut_point = s1.size()+1;
#ifdef HAVE_VIENNA20
  e=Vienna::energy_of_structure(ss.c_str(), rr.c_str(), -1);
#else
  Vienna::eos_debug = -1;
  e=Vienna::energy_of_struct(ss.c_str(), rr.c_str());
#endif
  Vienna::cut_point = -1;

  return e;
}

int
RactIP::
run()
{
  // set the energy parameters
  if (use_bl_param_)
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
  std::string r1, r2, r1s, r2s;
  float /*ea, */ e1, e2, e3, e1s, e2s;
  if (show_energy_ || enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
  {
    /*ea = */ solve(fa1, fa2, r1, r2, &e1, &e2, &e3);
    solve_ss(fa1.seq(), bp1_, offset1_, r1s, &e1s);
    solve_ss(fa2.seq(), bp2_, offset2_, r2s, &e2s);
  }
  else
    /*ea = */ solve(fa1, fa2, r1, r2);

  // display the result
  std::cout << ">" << fa1.name() << std::endl
            << fa1.seq() << std::endl << r1 << std::endl
            << ">" << fa2.name() << std::endl
            << fa2.seq() << std::endl << r2 << std::endl;

  // show energy of the joint structure
  if (show_energy_ || enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
  {
    if (show_energy_)
      std::cout << "(E: JS= " << e1+e2+e3
                << " = " << e1
                << (e2>=0.0 ? "+" : "") << e2 
                << (e3>=0.0 ? "+" : "") << e3 << ", "
                << "S1+S2= " << e1s+e2s << " = " << e1s
                << (e2s>=0.0 ? "+" : "") << e2s
                << ")" << std::endl;

    if (enable_zscore_==1 || enable_zscore_==2 || enable_zscore_==12)
    {
      float sum=0.0, sum2=0.0;
      float sums=0.0, sum2s=0.0;
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
        float ee1, ee2, ee3, ee1s, ee2s;
        solve(Fasta(fa1.name(), s1), Fasta(fa2.name(), s2), r1, r2, &ee1, &ee2, &ee3);
        solve_ss(s1, bp1_, offset1_, r1s, &ee1s);
        solve_ss(s2, bp2_, offset2_, r2s, &ee2s);
        float ee=ee1+ee2+ee3;
        float ees=ee-ee1s-ee2s;
#if 0
        std::cout << s1 << std::endl << r1 << std::endl
                  << s2 << std::endl << r2 << std::endl;
        std::cout << ee1 << " + " << ee2 << " + " << ee3 << " = " << ee << std::endl;
#endif
        sum += ee; sum2 += ee*ee;
        sums += ees; sum2s += ees*ees;
      }
      float m=sum/num_shuffling_;
      float v=sum2/num_shuffling_-m*m;
      if (v<0.0) v=0.0;
      float ms=sums/num_shuffling_;
      float vs=sum2s/num_shuffling_-ms*ms;
      if (vs<0.0) vs=0.0;
#if 0
      std::cout << m << " " << v << std::endl;
#endif
      std::cout << "z-score: "
                << (e1+e2+e3-m)/sqrt(v) << ", "
                << (e1+e2+e3-e1s-e2s-ms)/sqrt(vs) << std::endl;
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
