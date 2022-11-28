//
//  Latentnode.h
//  
//
//  Created by Farideh Fazayeli on 7/15/13.
//
//

#ifndef ____latentNode__
#define ____latentNode__

#include <iostream>
#include <memory>
#include "utillity.h"
using namespace std;
class NodeLatent{
    
private:
	std::unique_ptr<double[]> latent_factor_;  //a vector of size num_feat_ storing latent factor for this node
	std::unique_ptr<double[]> old_latent_factor_;   //store the old value for gradient descent update
	std::unique_ptr<double[]> prior_;      // a vector of size num_feat_ storing latent factor for parents and children (add together)
	std::unique_ptr<double[]> gradient_; // a vecor of size num_feat_ storing the gradient for the gradient descent update
  //    int numPrior;       // Number of parents + children
    
  int num_parent_;         //number of parents
  int num_child_;       //number of children
  int last_child_;
  int num_observ_;
  vector<float> observ_;
  vector<int> observ_idx_;

  int num_feat_;
  int level_;
 
  NodeLatent* parents_;
  NodeLatent** children_;
    
  NodeLatent() { 
    num_parent_ = 0; 
    num_child_ = 0; 
  }
    
 public:
  ~NodeLatent();
  NodeLatent(int child, int feat, int lev);
    
  NodeLatent(int par, int child, int parIdx, NodeLatent ***tree, int feat, int lev, int parLevel);    //getting num parent and num children and set numPrior
    
  void UpdateGradient(int freq, double err, double *vv);
  void UpdateLatentFactor(double coef, double momentum, double epsilon);
    
  void UpdatePrior(int mode);
    
  inline double *get_latent_factor() { 
    return latent_factor_.get(); 
  }
    
  inline void PrintVars(int i) {
    cout << "par: " << parents_ << endl;
    PrintArr(parents_->get_latent_factor(), num_feat_);
  }
    
  inline void PrintAdd(){
    cout << "add: " << this << endl;
  }
    
  inline void PrintPrior(){
    PrintArr(prior_.get(), num_feat_);
  }
    
  void set_observ(float obs, int ind);
    
  void GibbsUpdate(NodeLatent*** vTree, double sigInv, double sigUInv, int level);
};

//in preprocessing, we initialize the following of each node, based on phyinfo
//numPrior: num_parent + num_child
//link to parents and children.


#endif /* defined(____latentNode__) */

