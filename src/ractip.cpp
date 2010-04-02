#include "config.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <boost/multi_array.hpp>
#include "fa.h"
#include "gurobi_c++.h"

#include "contrafold/SStruct.hpp"
#include "contrafold/InferenceEngine.hpp"
#include "contrafold/DuplexEngine.hpp"
#include "contrafold/Defaults.ipp"

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include "part_func.h"
#include "pf_duplex.h"
};
};

typedef unsigned int uint;

class RNAIP
{
public:
  RNAIP(float th_hy, float th_ss, float alpha, bool in_pk,
        bool use_contrafold, bool use_stacked_bp, bool stacking_constraints, int n_th, const char* rip_file)
    : env_(NULL),
      model_(NULL),
      th_hy_(th_hy),
      th_ss_(th_ss),
      alpha_(alpha),
      in_pk_(in_pk),
      use_contrafold_(use_contrafold),
      use_stacked_bp_(use_stacked_bp),
      stacking_constraints_(stacking_constraints),
      n_th_(n_th),
      rip_file_(rip_file)
  {
    env_ = new GRBEnv;
    env_->set(GRB_IntParam_Threads, n_th_); // # of threads
    model_ = new GRBModel(*env_);
  }

  ~RNAIP()
  {
    delete model_;
    delete env_;
  }

  void solve(const std::string& s1, const std::string& s2,
             std::string& r1, std::string& r2);

private:
  void contrafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<bool, 2>& e);
  void contrafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
                  boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee);
  void contraduplex(const std::string& seq1, const std::string& seq2, boost::multi_array<GRBVar, 2>& v,
                    boost::multi_array<bool, 2>& e);
  void contraduplex(const std::string& seq1, const std::string& seq2,
                    boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
                    boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee);
  void rnafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<bool, 2>& e);
  void rnafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
               boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee);
  void rnaduplex(const std::string& seq1, const std::string& seq2, boost::multi_array<GRBVar, 2>& v,
                 boost::multi_array<bool, 2>& e);
  void rnaduplex(const std::string& seq1, const std::string& seq2,
                 boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
                 boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee);
  void load_from_rip(const char* filename);

private:
  GRBEnv* env_;
  GRBModel* model_;
  
  // options
  float th_hy_;                // threshold for the hybridization probability
  float th_ss_;                // threshold for the base-pairing probability
  float alpha_;                // weight for the hybridization score
  bool in_pk_;                 // allow internal pseudoknots or not
  bool use_contrafold_;        // use CONTRAfold model or not
  bool use_stacked_bp_;
  bool stacking_constraints_;
  int n_th_;                   // the number of threads
  const char* rip_file_;
  
  // binary variables
  boost::multi_array<GRBVar, 2> x_;
  boost::multi_array<GRBVar, 2> y_;
  boost::multi_array<GRBVar, 2> z_;
  boost::multi_array<GRBVar, 2> xx_;
  boost::multi_array<GRBVar, 2> yy_;
  boost::multi_array<GRBVar, 2> zz_;
  std::vector<GRBVar> x_up_;
  std::vector<GRBVar> x_down_;
  std::vector<GRBVar> y_up_;
  std::vector<GRBVar> y_down_;
  std::vector<GRBVar> z_up_;
  std::vector<GRBVar> z_down_;

  boost::multi_array<bool, 2> ex_;
  boost::multi_array<bool, 2> ey_;
  boost::multi_array<bool, 2> ez_;
  boost::multi_array<bool, 2> exx_;
  boost::multi_array<bool, 2> eyy_;
  boost::multi_array<bool, 2> ezz_;
};

void
RNAIP::
contrafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<bool, 2>& e) 
{
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  std::vector<float> bp((seq.size()+1)*(seq.size()+2)/2, 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(th_ss_, bp);

  for (uint j=1; j!=seq.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      float p=bp[en.GetOffset(i+1)+(j+1)];
      if (p>th_ss_)
      {
        v[i][j] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
        e[i][j] = true;
      }
    }
  }
}    

void
RNAIP::
contrafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
           boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee) 
{
  SStruct ss("unknown", seq);
  ParameterManager<float> pm;
  InferenceEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  std::vector<float> bp((seq.size()+1)*(seq.size()+2)/2, 0.0);
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior2(th_ss_, bp);

  for (uint j=1; j!=seq.size(); ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      if (i+1<j-1)
      {
        float p=bp[en.GetOffset(i+1)+(j+1)];
        if (p>th_ss_)
        {
          v[i][j] = model_->addVar(0, 1, 0, GRB_BINARY);
          vv[i][j] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
          e[i][j] = ee[i][j] = true;
          if (!e[i+1][j-1])
          {
            v[i+1][j-1] = model_->addVar(0, 1, 0, GRB_BINARY);
            e[i+1][j-1] = true;
          }
        }
      }
    }
  }
}    

void
RNAIP::
contraduplex(const std::string& seq1, const std::string& seq2, boost::multi_array<GRBVar, 2>& v,
             boost::multi_array<bool, 2>& e)
{
  SStruct ss1("unknown", seq1), ss2("unknown", seq2);
  ParameterManager<float> pm;
  DuplexEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  std::vector<float> ip;
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss1, ss2);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior(th_hy_, ip);
  for (uint i=0; i!=seq1.size(); ++i)
  {
    for (uint j=0; j!=seq2.size(); ++j)
    {
      float p=ip[en.GetOffset(i+1)+(j+1)];
      if (p>th_hy_)
      {
        v[i][j] = model_->addVar(0.0, 1.0, -p*alpha_, GRB_BINARY);
        e[i][j] = true;
      }
    }
  }
}

void
RNAIP::
contraduplex(const std::string& seq1, const std::string& seq2,
             boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
             boost::multi_array<bool,2>& e, boost::multi_array<bool,2>& ee)
{
  SStruct ss1("unknown", seq1), ss2("unknown", seq2);
  ParameterManager<float> pm;
  DuplexEngine<float> en(false);
  std::vector<float> w = GetDefaultComplementaryValues<float>();
  std::vector<float> ip;
  en.RegisterParameters(pm);
  en.LoadValues(w);
  en.LoadSequence(ss1, ss2);
  en.ComputeInside();
  en.ComputeOutside();
  en.ComputePosterior();
  en.GetPosterior2(th_hy_, ip);
  for (uint i=1; i!=seq1.size(); ++i)
  {
    for (uint j=seq2.size()-1; j!=-1u; --j)
    {
      float p=ip[en.GetOffset(i+1)+(j+1)];
      if (p>th_hy_)
      {
        v[i][j] = model_->addVar(0.0, 1.0, 0, GRB_BINARY);
        vv[i][j] = model_->addVar(0.0, 1.0, -p*alpha_, GRB_BINARY);
        e[i][j] = ee[i][j] = true;
        if (!e[i-1][j+1])
        {
          v[i-1][j+1] = model_->addVar(0.0, 1.0, 0, GRB_BINARY);
          e[i-1][j+1] = true;
        }
      }
    }
  }
}


void
RNAIP::
rnafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<bool, 2>& e)
{
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  for (uint i=0; i!=seq.size()-1; ++i)
  {
    for (uint j=i+1; j!=seq.size(); ++j)
    {
      float p=Vienna::pr[Vienna::iindx[i+1]-(j+1)];
      if (p>th_ss_)
      {
        v[i][j] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
        e[i][j] = true;
      }
    }
  }
  Vienna::free_pf_arrays();
}

void
RNAIP::
rnafold(const std::string& seq, boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
        boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee)
{
  Vienna::init_pf_fold(seq.size());
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);
  for (uint i=0; i!=seq.size()-1; ++i)
  {
    for (uint j=i+1; j!=seq.size(); ++j)
    {
      if (i+1<j-1)
      {
        float p=Vienna::pr2[Vienna::iindx[i+1]-(j+1)];
        if (p>th_ss_)
        {
          v[i][j] = model_->addVar(0, 1, 0, GRB_BINARY);
          vv[i][j] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
          e[i][j] = ee[i][j] = true;
          if (!e[i+1][j-1])
          {
            v[i+1][j-1] = model_->addVar(0, 1, 0, GRB_BINARY);
            e[i+1][j-1] = true;
          }
        }
      }
    }
  }
  Vienna::free_pf_arrays();
}

void
RNAIP::
rnaduplex(const std::string& s1, const std::string& s2, boost::multi_array<GRBVar, 2>& v,
          boost::multi_array<bool, 2>& e)
{
  Vienna::pf_duplex(s1.c_str(), s2.c_str());
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      float p=Vienna::pr_duplex[i+1][j+1];
      if (p>th_hy_)
      {
        v[i][j] = model_->addVar(0.0, 1.0, -p*alpha_, GRB_BINARY);
        e[i][j] = true;
      }
    }
  }
  Vienna::free_pf_duplex();
}

void
RNAIP::
rnaduplex(const std::string& s1, const std::string& s2,
          boost::multi_array<GRBVar, 2>& v, boost::multi_array<GRBVar, 2>& vv,
          boost::multi_array<bool, 2>& e, boost::multi_array<bool, 2>& ee)
{
  Vienna::pf_duplex(s1.c_str(), s2.c_str());
  for (uint i=1; i!=s1.size(); ++i)
  {
    for (uint j=s2.size()-1; j!=-1u; --j)
    {
      float p=Vienna::pr_duplex2[i+1][j+1];
      if (p>th_hy_)
      {
        v[i][j] = model_->addVar(0.0, 1.0, 0, GRB_BINARY);
        vv[i][j] = model_->addVar(0.0, 1.0, -p*alpha_, GRB_BINARY);
        e[i][j] = ee[i][j] = true;
        if (!e[i-1][j+1])
        {
          v[i-1][j+1] = model_->addVar(0.0, 1.0, 0, GRB_BINARY);
          e[i-1][j+1] = true;
        }
      }
    }
  }
  Vienna::free_pf_duplex();
}

void
RNAIP::
load_from_rip(const char* filename)
{
  enum { NONE, TABLE_R, TABLE_S, TABLE_I };
  uint st=NONE;
  uint y_len=y_.size();
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
          if (p>th_ss_)
          {
            x_[i-1][j-1] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
            ex_[i-1][j-1] = true;
          }
          break;
        case TABLE_S:
          if (p>th_ss_)
          {
            y_[y_len-j][y_len-i] = model_->addVar(0, 1, -p*(1-alpha_), GRB_BINARY);
            ey_[y_len-j][y_len-i] = true;
          }
          break;
        case TABLE_I:
          if (p>th_hy_)
          {
            z_[i-1][y_len-j] = model_->addVar(0.0, 1.0, -p*alpha_, GRB_BINARY);
            ez_[i-1][y_len-j] = true;
          }
          break;
        default:
          break;
      }
    }
    else st=NONE;
  }
}

void
RNAIP::
solve(const std::string& s1, const std::string& s2, std::string& r1, std::string& r2)
{
  x_.resize(boost::extents[s1.size()][s1.size()]);
  xx_.resize(boost::extents[s1.size()][s1.size()]);
  ex_.resize(boost::extents[s1.size()][s1.size()]);
  exx_.resize(boost::extents[s1.size()][s1.size()]);
  y_.resize(boost::extents[s2.size()][s2.size()]);
  yy_.resize(boost::extents[s2.size()][s2.size()]);
  ey_.resize(boost::extents[s2.size()][s2.size()]);
  eyy_.resize(boost::extents[s2.size()][s2.size()]);
  z_.resize(boost::extents[s1.size()][s2.size()]);
  zz_.resize(boost::extents[s1.size()][s2.size()]);
  ez_.resize(boost::extents[s1.size()][s2.size()]);
  ezz_.resize(boost::extents[s1.size()][s2.size()]);

  if (rip_file_)
  {
    load_from_rip(rip_file_);
  }
  else
  {
    if (use_contrafold_)
    {
      if (use_stacked_bp_)
      {
        contrafold(s1, x_, xx_, ex_, exx_);
        contrafold(s2, y_, yy_, ey_, eyy_);
        //contraduplex(s1, s2, z_, zz_, ez_, ezz_);
        rnaduplex(s1, s2, z_, zz_, ez_, ezz_);
      }
      else
      {
        contrafold(s1, x_, ex_);
        contrafold(s2, y_, ey_);
        //contraduplex(s1, s2, z_, ez_);
        rnaduplex(s1, s2, z_, ez_);
      }
    }
    else
    {
      Vienna::pf_scale = -1;
      if (use_stacked_bp_)
      {
        rnafold(s1, x_, xx_, ex_, exx_);
        rnafold(s2, y_, yy_, ey_, eyy_);
        rnaduplex(s1, s2, z_, zz_, ez_, ezz_);
      }
      else
      {
        rnafold(s1, x_, ex_);
        rnafold(s2, y_, ey_);
        rnaduplex(s1, s2, z_, ez_);
      }
    }
  }

  if (stacking_constraints_)
  {
    x_up_.resize(s1.size());
    x_down_.resize(s1.size());
    for (uint i=0; i!=s1.size(); ++i)
    {
      x_up_[i] = model_->addVar(0, 1, 0, GRB_BINARY);
      x_down_[i] = model_->addVar(0, 1, 0, GRB_BINARY);
    }

    y_up_.resize(s2.size());
    y_down_.resize(s2.size());
    for (uint i=0; i!=s2.size(); ++i)
    {
      y_up_[i] = model_->addVar(0, 1, 0, GRB_BINARY);
      y_down_[i] = model_->addVar(0, 1, 0, GRB_BINARY);
    }
  
    z_up_.resize(s2.size());
    z_down_.resize(s1.size());
    for (uint j=0; j!=s2.size(); ++j)
      z_up_[j] = model_->addVar(0, 1, 0, GRB_BINARY);
    for (uint i=0; i!=s1.size(); ++i)
      z_down_[i] = model_->addVar(0, 1, 0, GRB_BINARY);
  }

  model_->update();

  // constraint 1: each a_i is paired with at most one base
  for (uint i=0; i!=s1.size(); ++i)
  {
    GRBLinExpr c1;
    for (uint j=0; j!=s2.size(); ++j) if (ez_[i][j]) c1 += z_[i][j];
    for (uint j=0; j<i; ++j) if (ex_[j][i]) c1 += x_[j][i];
    for (uint j=i+1; j<s1.size(); ++j) if (ex_[i][j]) c1 += x_[i][j];
    model_->addConstr(c1 <= 1);
  }

  // constraint 2: each b_i is paired with at moat one
  for (uint i=0; i!=s2.size(); ++i)
  {
    GRBLinExpr c2;
    for (uint j=0; j!=s1.size(); ++j) if (ez_[j][i]) c2 += z_[j][i];
    for (uint j=0; j<i; ++j) if (ey_[j][i]) c2 += y_[j][i];
    for (uint j=i+1; j<s2.size(); ++j) if (ey_[i][j]) c2 += y_[i][j];
    model_->addConstr(c2 <= 1);
  }

  // constraint 3: disallow external pseudoknots
  for (uint i=0; i<s1.size(); ++i)
    for (uint k=i+1; k<s1.size(); ++k)
      for (uint j=0; j<s2.size(); ++j)
        if (ez_[i][j])
          for (uint l=j+1; l<s2.size(); ++l)
            if (ez_[k][l])
              model_->addConstr(z_[i][j]+z_[k][l] <= 1);

  if (in_pk_)
  {
    // constraint 4: disallow internal pseudoknots in a
    for (uint i=0; i<s1.size(); ++i)
      for (uint k=i+1; k<s1.size(); ++k)
        for (uint j=k+1; j<s1.size(); ++j)
          if (ex_[i][j])
            for (uint l=j+1; l<s1.size(); ++l)
              if (ex_[k][l])
                model_->addConstr(x_[i][j]+x_[k][l] <= 1);

    // constraint 5: disallow internal pseudoknots in b
    for (uint i=0; i<s2.size(); ++i)
      for (uint k=i+1; k<s2.size(); ++k)
        for (uint j=k+1; j<s2.size(); ++j)
          if (ey_[i][j])
            for (uint l=j+1; l<s2.size(); ++l)
              if (ey_[k][l])
                model_->addConstr(y_[i][j]+y_[k][l] <= 1);
  }

  if (stacking_constraints_)
  {
    for (uint i=0; i<s1.size(); ++i)
    {
      GRBLinExpr c_up1;
      for (uint j=0; j<i; ++j) if (ex_[j][i]) c_up1 += x_[j][i];
      model_->addConstr(c_up1-x_up_[i] == 0);
      GRBLinExpr c_down1;
      for (uint j=i+1; j<s1.size(); ++j) if (ex_[i][j]) c_down1 += x_[i][j];
      model_->addConstr(c_down1-x_down_[i] == 0);

      GRBLinExpr c_up2 = 1-x_up_[i];
      if (i>0) c_up2 += x_up_[i-1];
      if (i+1<s1.size()) c_up2 += x_up_[i+1];
      model_->addConstr(c_up2 >= 1);
      GRBLinExpr c_down2 = 1-x_down_[i];
      if (i>0) c_down2 += x_down_[i-1];
      if (i+1<s1.size()) c_down2 += x_down_[i+1];
      model_->addConstr(c_down2 >= 1);
    }

    for (uint i=0; i<s2.size(); ++i)
    {
      GRBLinExpr c_up1;
      for (uint j=0; j<i; ++j) if (ey_[j][i]) c_up1 += y_[j][i];
      model_->addConstr(c_up1-y_up_[i] == 0);
      GRBLinExpr c_down1;
      for (uint j=i+1; j<s2.size(); ++j) if (ey_[i][j]) c_down1 += y_[i][j];
      model_->addConstr(c_down1-y_down_[i] == 0);

      GRBLinExpr c_up2 = 1-y_up_[i];
      if (i>0) c_up2 += y_up_[i-1];
      if (i+1<s2.size()) c_up2 += y_up_[i+1];
      model_->addConstr(c_up2 >= 1);
      GRBLinExpr c_down2 = 1-y_down_[i];
      if (i>0) c_down2 += y_down_[i-1];
      if (i+1<s2.size()) c_down2 += y_down_[i+1];
      model_->addConstr(c_down2 >= 1);
    }

    for (uint i=0; i<s2.size(); ++i)
    {
      GRBLinExpr c_up1;
      for (uint j=0; j<s1.size(); ++j) if (ez_[j][i]) c_up1 += z_[j][i];
      model_->addConstr(c_up1-z_up_[i] == 0);

      GRBLinExpr c_up2 = 1-z_up_[i];
      if (i>0) c_up2 += z_up_[i-1];
      if (i+1<s2.size()) c_up2 += z_up_[i+1];
      model_->addConstr(c_up2 >= 1);
    }

    for (uint i=0; i<s1.size(); ++i)
    {
      GRBLinExpr c_down1;
      for (uint j=0; j<s2.size(); ++j) if (ez_[i][j]) c_down1 += z_[i][j];
      model_->addConstr(c_down1-z_down_[i] == 0);

      GRBLinExpr c_down2 = 1-z_down_[i];
      if (i>0) c_down2 += z_down_[i-1];
      if (i+1<s1.size()) c_down2 += z_down_[i+1];
      model_->addConstr(c_down2 >= 1);
    }
  }

  // stacked base-pairing probabilities
  if (use_stacked_bp_)
  {
    for (uint j=1; j!=s1.size(); ++j)
      for (uint i=j-1; i!=-1u; --i)
        if (exx_[i][j])
          model_->addConstr(2*xx_[i][j] <= x_[i][j]+x_[i+1][j-1]);
    for (uint j=1; j!=s2.size(); ++j)
      for (uint i=j-1; i!=-1u; --i)
        if (eyy_[i][j])
          model_->addConstr(2*yy_[i][j] <= y_[i][j]+y_[i+1][j-1]);
    for (uint i=1; i<s1.size(); ++i)
      for (uint j=0; j<s2.size()-1; ++j)
        if (ezz_[i][j])
          model_->addConstr(2*zz_[i][j] <= z_[i][j]+z_[i-1][j+1]);
  }

  // execute optimization
  model_->optimize();

  // build the resultant structure
  r1.resize(s1.size());
  r2.resize(s2.size());
  std::fill(r1.begin(), r1.end(), '.');
  std::fill(r2.begin(), r2.end(), '.');
  for (uint i=0; i!=s1.size(); ++i)
  {
    for (uint j=0; j!=s2.size(); ++j)
    {
      if (use_stacked_bp_)
      {
        if (ezz_[i][j] && zz_[i][j].get(GRB_DoubleAttr_X)>0.5)
        {
          r1[i]='['; r2[j]=']';
          r1[i-1]='['; r2[j+1]=']';
        }
      }
      else
      {
        if (ez_[i][j] && z_[i][j].get(GRB_DoubleAttr_X)>0.5)
        {
          r1[i]='['; r2[j]=']';
        }
      }
    }
  }
  if (in_pk_)
  {
    for (uint i=0; i<s1.size(); ++i)
    {
      for (uint j=i+1; j<s1.size(); ++j)
      {
        if (use_stacked_bp_)
        {
          if (exx_[i][j] && xx_[i][j].get(GRB_DoubleAttr_X)>0.5)
          {
            assert(r1[i]=='.' || r1[i]=='('); assert(r1[j]=='.' || r1[j]==')');
            assert(r1[i+1]=='.'); assert(r1[j-1]=='.');
            r1[i]='('; r1[j]=')';
            r1[i+1]='('; r1[j-1]=')';
          }
        }
        else
        {
          if (ex_[i][j] && x_[i][j].get(GRB_DoubleAttr_X)>0.5)
          {
            assert(r1[i]=='.'); assert(r1[j]=='.');
            r1[i]='('; r1[j]=')';
          }
        }
      }
    }
    for (uint i=0; i<s2.size(); ++i)
    {
      for (uint j=i+1; j<s2.size(); ++j)
      {
        if (use_stacked_bp_)
        {
          if (eyy_[i][j] && yy_[i][j].get(GRB_DoubleAttr_X)>0.5)
          {
            assert(r2[i]=='.' || r2[i]=='('); assert(r2[j]=='.' || r2[j]==')');
            assert(r2[i+1]=='.'); assert(r2[j-1]=='.');
            r2[i]='('; r2[j]=')';
            r2[i+1]='('; r2[j-1]=')';
          }
        }
        else
        {
          if (ey_[i][j] && y_[i][j].get(GRB_DoubleAttr_X)>0.5)
          {
            assert(r2[i]=='.'); assert(r2[j]=='.');
            r2[i]='('; r2[j]=')';
          }
        }
      }
    }
  }
}

void
usage(const char* progname)
{
  std::cout << progname << ": [options] fasta1 fasta2" << std::endl
            << " -h:        show this message" << std::endl
            << " -p:        use the constraints for interenal pseudoknots" << std::endl
            << " -a alpha:  weight for hybridation probabilities" << std::endl
            << " -t th_bp:  threshold of base-pairing probabilities" << std::endl
            << " -u th_hy:  threshold of hybridazation probabilities" << std::endl
            << " -c:        use CONTRAfold model (default: McCaskill model)" << std::endl
            << " -i:        allow isolated base-pairs" << std::endl
            << " -s:        use stacked base-pairing probabilities instead of standard base-pairing probabilities" << std::endl
            << " -n n_th:   specify the number of threads" << std::endl;
}

int
main(int argc, char* argv[])
{
  char* progname=argv[0];
  // parse options
  char ch;
  float alpha=1.0;
  float th_bp=0.0;
  float th_hy=0.0;
  bool in_pk=false;
  bool isolated_bp=false;
  bool stacked_bp=false;
  bool use_contrafold=false;
  int n_th=1;
  const char* rip_file=NULL;
  while ((ch=getopt(argc, argv, "a:t:u:pcisn:r:h"))!=-1)
  {
    switch (ch)
    {
      case 'c':
        use_contrafold=true;
        break;
      case 'a':
        alpha=atof(optarg);
        break;
      case 't':
        th_bp=atof(optarg);
        break;
      case 'u':
        th_hy=atof(optarg);
        break;
      case 'p':
        in_pk=true;
        break;
      case 's':
        stacked_bp=true;
        break;
      case 'i':
        isolated_bp=true;
        break;
      case 'n':
        n_th=atoi(optarg);
        break;
      case 'r':
        rip_file=optarg;
        break;
      case 'h': case '?': default:
        usage(progname);
        return 1;
        break;
    }
  }
  argc -= optind;
  argv += optind;

  // read sequences
  Fasta fa1, fa2;
  if (argc==2)
  {
    std::list<Fasta> l1, l2;
    Fasta::load(l1, argv[0]);
    Fasta::load(l2, argv[1]);
    fa1=l1.front();
    fa2=l2.front();
  }
  else if (argc==1)
  {
    std::list<Fasta> l1;
    Fasta::load(l1, argv[0]);
    if (l1.size()<2) { usage(progname); return 1; }
    std::list<Fasta>::const_iterator x=l1.begin();
    fa1=*(x++);
    fa2=*(x++);
  }
  else { usage(progname); return 1;}

  // predict the interation
  std::string r1, r2;
  //try {
    RNAIP rnaip(th_hy, th_bp, alpha, in_pk,
                use_contrafold, stacked_bp, !isolated_bp, n_th, rip_file);
    rnaip.solve(fa1.seq(), fa2.seq(), r1, r2);
#if 0
  } catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    throw;
  }
#endif

  // display the result
  std::cout << ">" << fa1.name() << std::endl
            << fa1.seq() << std::endl << r1 << std::endl
            << ">" << fa2.name() << std::endl
            << fa2.seq() << std::endl << r2 << std::endl;

  return 0;
}
