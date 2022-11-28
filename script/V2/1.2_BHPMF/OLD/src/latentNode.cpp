//
//  latentNode.cpp
//  
//
//  Created by Farideh Fazayeli on 7/15/13.
//
//

#include "latentNode.h"
#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

NodeLatent::NodeLatent(int num_child, int num_feat, int level) :
  num_parent_(0),
  num_child_(num_child),
  num_feat_(num_feat),
  level_(level),
  last_child_(-1),
  num_observ_(0) {    
  if (num_child_) {
    //        prior_ = new double[num_feat_];
    children_ = new NodeLatent*[num_child_];
    prior_ = std::unique_ptr<double[]>(new double[num_feat_]);
  }
    
  //    latent_factor_ = new double[num_feat_];
  latent_factor_ = std::unique_ptr<double[]>(new double[num_feat_]);
  //fill_n(latent_factor_, num_feat_, 0);
  for(int i = 0; i < num_feat_; i++) {
    latent_factor_[i] = 2 * unif_rand() -1;
  }

  old_latent_factor_ = std::unique_ptr<double[]>(new double[num_feat_]);
  gradient_ = std::unique_ptr<double[]>(new double[num_feat_]);
  //    old_latent_factor_ = new double[num_feat_];
}

NodeLatent::NodeLatent(int num_parent, int num_child, int par_idx, NodeLatent ***tree, int num_feat, int level, int par_level) :
  num_parent_(num_parent),
  num_child_(num_child),
  num_feat_(num_feat),
  level_(level),
  last_child_(-1),
  num_observ_(0) {
  //if (par) 
  //        parents_ = new NodeLatent*[par];

  if (num_child)
    children_ = new NodeLatent*[num_child];
     
  for (int ind = 0; ind < num_parent; ind++){
    //         parents_[ind] = tree[parLevel_][parIdx[ind]-1];
    parents_ = tree[par_level][par_idx-1];
    parents_->children_[++parents_->last_child_] = this;
  }
    
  if (num_parent_ || num_child_) { 
	  prior_ = std::unique_ptr<double[]>(new double[num_feat_]);
  }
  //      prior_ = new double[num_feat_];

  latent_factor_ = std::unique_ptr<double[]>(new double[num_feat_]);
  //    latent_factor_ = new double[num_feat_];
  //fill_n(latent_factor_, num_feat_, 0); 
  for(int i = 0; i < num_feat_; i++)
    latent_factor_[i] = 2 * unif_rand() -1;
		
  //    old_latent_factor_ = new double[num_feat_];  
  old_latent_factor_ = std::unique_ptr<double[]>(new double[num_feat_]);
  gradient_ = std::unique_ptr<double[]>(new double[num_feat_]);
}



NodeLatent::~NodeLatent(){
    
  //    SafeDelete(prior);
  //  SafeDelete(latent_factor_);
  //  SafeDelete(old_latent_factor_);

  //  delete[] old_latent_factor_;
  // delete[] latent_factor_;
  // delete[] prior_;
    
    //    for (int ind = 0; ind < num_parent_; ind++)
      //        SafeDelete(parents_[ind]); 
    //  delete parents_[ind];
    
  /*
    for(int ind = 0; ind < num_child_; ind++)
      delete children_[ind];
  */
  SafeDelete(children_);
    //    SafeDelete(parents_);
    
}

void NodeLatent::set_observ(float observ, int ind){
    observ_.push_back(observ);
    observ_idx_.push_back(ind);
    num_observ_++;
}


void NodeLatent::UpdatePrior(int mode) {
    const int inc_one = 1;
    const double cc = 1;

	std::unique_ptr<double[]> prior_parent(new double[num_feat_]);
	std::unique_ptr<double[]> prior_child(new double[num_feat_]);

    fill_n(prior_parent.get(), num_feat_, 0);
    fill_n(prior_child.get(), num_feat_, 0);  
    fill_n(old_latent_factor_.get(), num_feat_, 0);

    for (int ind = 0; ind < num_parent_; ind++) {
		F77_NAME(daxpy)(&num_feat_, &cc, parents_->get_latent_factor(), &inc_one, prior_parent.get(), &inc_one);
    }    
    for (int ind = 0; ind < num_child_; ind++) {
      F77_NAME(daxpy)(&num_feat_, &cc, children_[ind]->get_latent_factor(), &inc_one, prior_child.get(), &inc_one);
    }
    
    if (mode == 1) {
      for (int ii = 0; ii < num_feat_; ii++) {
		  latent_factor_[ii] = prior_parent[ii]/num_parent_ + rnorm(0, 0.01);    // Plant feature vectors
      }
    } else {
      for (int ii = 0; ii < num_feat_; ii++) {
		  latent_factor_[ii] = prior_child[ii]/num_child_ + rnorm(0, 0.01); // Plant feature vectors        
      }
    }
    
    fill_n(prior_.get(), num_feat_, 0);
    F77_NAME(daxpy)(&num_feat_, &cc, prior_parent.get(), &inc_one, prior_.get(), &inc_one);
    F77_NAME(daxpy)(&num_feat_, &cc, prior_child.get(), &inc_one, prior_.get(), &inc_one);
}

void NodeLatent::UpdateGradient(int freq, double err, double *vv) {
  int inc_one = 1;
  if (freq == 0) {
    fill_n(gradient_.get(), num_feat_, 0);
  }
  F77_NAME(daxpy)(&num_feat_, &err, vv, &inc_one, gradient_.get(), &inc_one);  // dE_du[i] = \sum_j^N_i err_ij * V_j
}

void NodeLatent::UpdateLatentFactor(double coef, double momentum, double epsilon) {
    //can be improved by writting a simple for loop to do all this operations!!! ???????????????????????????????
    const int inc_one = 1;
    double cc = -1 * coef;
        
    F77_NAME(daxpy)(&num_feat_, &cc, prior_.get(), &inc_one, gradient_.get(), &inc_one);    // deV[j] -= coef * ( sum(v_pars) + sum(v_chs) )  ;
    
    coef *= (num_parent_+num_child_);  //numPrior = (num_parent_ + numCh)  ;  coef = N_j * lambda * (num_parent_ + numCh)_v
    F77_NAME(daxpy)(&num_feat_, &coef, latent_factor_.get(), &inc_one, gradient_.get(), &inc_one);  // deV[j] += coef * v[j]
    
    F77_NAME(dscal)(&num_feat_, &momentum, old_latent_factor_.get(), &inc_one);
    F77_NAME(daxpy)(&num_feat_, &epsilon, gradient_.get(), &inc_one, old_latent_factor_.get(), &inc_one);  // old = old * momentum + epsilon * deV[j]

    cc = -1;
    F77_NAME(daxpy)(&num_feat_, &cc, old_latent_factor_.get(), &inc_one, latent_factor_.get(), &inc_one);  // latent_factor_ -= old;
}

void InverseMat(double* mat, int dim)
{
	std::unique_ptr<int[]> ipiv(new int [dim+1]);
    int mat_size = dim * dim;
	std::unique_ptr<double[]> tmp(new double [mat_size]);
    int info;

    F77_NAME(dgetrf)(&dim, &dim, mat, &dim, ipiv.get(), &info);
    F77_NAME(dgetri)(&dim, mat, &dim, ipiv.get(), tmp.get(), &mat_size, &info);
}


void NodeLatent::GibbsUpdate(NodeLatent*** v_tree, double sig_inv, double sig_u_inv, int level){

  /************************************
	BLAS, Lapack Parameter
  *************************************/
  const char *lower = "L";
  const char *ntrans = "N";
  const char *ytrans = "T";
  // int num_feat_L = num_feat_;
  int info;

  const double zero = 0;
  const double one = 1;
  const int inc_one = 1;

  std::unique_ptr<double[]> prior_parent(new double[num_feat_]);
  std::unique_ptr<double[]> prior_child(new double[num_feat_]);
    
  fill_n(prior_parent.get(), num_feat_, 0);
  fill_n(prior_child.get(), num_feat_, 0);
  fill_n(old_latent_factor_.get(), num_feat_, 0);
    
  for (int ind = 0; ind < num_parent_; ind++)
    F77_NAME(daxpy)(&num_feat_, &one, parents_->get_latent_factor(), &inc_one, prior_parent.get(), &inc_one);
    
  //if(level_ < 4){
  for (int ind = 0; ind < num_child_; ind++)
    F77_NAME(daxpy)(&num_feat_, &one, children_[ind]->get_latent_factor(), &inc_one, prior_child.get(), &inc_one);
    
  fill_n(prior_.get(), num_feat_, 0);
  F77_NAME(daxpy)(&num_feat_, &one, prior_parent.get(), &inc_one, prior_.get(), &inc_one);
  F77_NAME(daxpy)(&num_feat_, &one, prior_child.get(), &inc_one, prior_.get(), &inc_one);
        
  /************************************
	Setup variables
  *************************************/
  double numPrior = (double) (num_child_ + num_parent_);  

  std::unique_ptr<double[]> mu(new double[num_feat_]);
  std::unique_ptr<double[]> mu_u(new double[num_feat_]);
  std::unique_ptr<double[]> cov(new double[num_feat_ * num_feat_]);
  std::unique_ptr<double[]> latent_trans(new double[num_observ_ * num_feat_]);

  // double *mu = new double[num_feat_];
  // double *mu_u = new double[num_feat_];
  // double *cov = new double [num_feat_ * num_feat_];
  // double *latentT = new double[num_observ_*num_feat_];

  fill_n(mu.get(), num_feat_, 0);
    
  for (int ind = 0; ind < num_observ_; ind++) {
    double beta = observ_[ind] * sig_inv;
    double *tmp = v_tree[level][observ_idx_[ind]]->get_latent_factor(); //?????????????????
        
    copy(tmp, tmp+num_feat_, latent_trans.get()+ind*num_feat_);
    //latentT[ind*num_feat_] = tmp[0];            
        
    F77_NAME(daxpy)(&num_feat_, &beta, tmp, &inc_one, mu.get(), &inc_one);     //mu += v_j * r_{ij} / sig
  }
    
  //mu += prior / sigU ; prior = sum(u_ch) + u_pr
  F77_NAME(daxpy)(&num_feat_, &sig_u_inv, prior_.get(), &inc_one, mu.get(), &inc_one);
    

  //cov = Identity / sigU
  for(int ii = 0; ii < num_feat_*num_feat_; ii++) 
    cov[ii] = 0.0;
  for(int ii = 0; ii < num_feat_; ii++)
    cov[ii*num_feat_+ii] = sig_u_inv;
	
  //cov = V_i' * V_i + (|ch| + 1)/sigUInv I;    V_i' = latentT, V_i sub matrix of V for non zero values of i
  F77_NAME(dgemm)(ntrans, ytrans, &num_feat_, &num_feat_, &num_observ_, &sig_inv, latent_trans.get(), &num_feat_, latent_trans.get(), &num_feat_, &numPrior, cov.get(), &num_feat_);
    
  //cov = inverse cov
  //dpotrf_(lower, &num_feat_L, cov, &num_feat_L, &info); if(info != 0) {cerr << "C++ error: Cholesky failed: " << info << endl;  exit (1);}
  //dpotri_(lower, &num_feat_L, cov, &num_feat_L, &info); if(info != 0) {cerr << "C++ error: Cholesky inverse failed." << endl; exit(1);}
    
  InverseMat(cov.get(), num_feat_);

  F77_NAME(dgemv)(ntrans, &num_feat_, &num_feat_, &one, cov.get(), &num_feat_, mu.get(), &inc_one, &zero, mu_u.get(), &inc_one);

  F77_NAME(dpotrf)(lower, &num_feat_, cov.get(), &num_feat_, &info); 
  if(info != 0) {
    cerr << "C++ error: Cholesky failed. " << endl;
    exit (1);
  }
	
  mvrnorm(latent_factor_.get(), mu_u.get(), cov.get(), num_feat_, false);	
}
