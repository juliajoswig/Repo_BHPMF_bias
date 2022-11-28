//
//  HPMF.h
//  
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#ifndef ____HPMF__
#define ____HPMF__

#include <iostream>
#include <string>
#include <algorithm>    // std::move_backward
#include <fenv.h>
//#include <array>        // std::array
//#include <random>       // std::default_random_engine
//#include <chrono>       // std::chrono::system_clock

#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <numeric>
//#include <Rinternals.h>

#include "utillity.h"
#include "latentNode.h"
#include "omp.h"

#include <vector>

class HPMF{
private:
    int numPlants; //number of plants
    int numTraits; //number of traits
   
    //latent matrices
    NodeLatent ***uTree;
    NodeLatent ***vTree;

	// hiearchy vector
	matVecInt hierarchy_;
	
	vector<int> num_parents_;
	vector<int> num_childern_;
//    vector<int> num_nodes_per_level_;

	vector<int> nodePerLevel;

//	bool countMaxV;
//	bool countMaxU;

//    gsl_rng *rGen;
    
    double rmseVal; //rmse value on validation data
    
    int lenTrainData;
    int lenValData;
    int lenTestData;
    unsigned char mode;
	int gaps_;
	int burn_;
	int num_samples_;
    //    int *nodePerLevel;
    
    int numFeat;  //size of latent factor $k$
    int taxonomyLevels;
    int numTaxonomyLevels;
	int num_effective_samples_;
    int tunedLevel;
	int predict_level_;
	int num_hierarchy_level_;
	int used_hierarchy_level_;
    int cvInd;
    int datasetId;
    const char* dataPath;
    int iterations;
    int opt;
    int saveFileFlag;
    bool outWholeFlag;
    SEXP env;
    matVecFloat all_preds_mean_;
    matVecFloat all_preds_std_;
	matVecFloat test_data_;
	vector<float> test_preds_mean_;
    vector<float> test_preds_std_;
	int num_row_test_;
	vector<double> num_nodes_per_level_;

	const char* meanFilePath_;
	const char* stdFilePath_;    
    int epsilon;
    double lambdaV;
    double lambdaU;
    double momentum;
//	double test_preds[20][49210];
  
 //   char * levelNames[5] = {"phylo", "family", "genus", "species", "tree"};
//    const string levelNames[5] = {"phylo", "family", "genus", "species", "tree"};
    
    void gibbsTrain(int level);
    void hpmfTrain(int level, matVecFloat trainData, matVecFloat valData);
    void intializeTree(NodeLatent*** tree, int numNodes);
	void LoadData(int level, matVecFloat* trainData, matVecFloat* valData);
	void Predict(int ind, int level, double* err, double *uu, double *vv, matVecFloat data);
    void setObservation();
    void updateMats(double ***uMat, double ***vMat, double ***pred, int iter);

public:
    
    ~HPMF();
    HPMF(const char *dataPath, int numTraits, int numFeat, int taxonomyLevels, int tunedLevel, int used_hierarchy_level, int cvInd, int datasetId, int iterations, int epsilon, double momentum, double lambdaV, double lambdaU, int opt, int saveFileFlag, int outWholeFlag, int gaps, int burn, int num_samples, const char* meanFilePath, const char* stdFilePath, double* num_nodes_per_level_, SEXP env);
    void hpmfTopdown();
    
    void hpmfDowntop();

    double runHPMF(bool saveFlag, int &numPlants);
    double runHPMF(int  numGibbs);
    double runHPMF(double *allRMSE, double *allRMSE2, int &numPlants);
    double findTestErr();
    double findTestErr(FILE * outfile);
	double findTestErr(vector<double> &sum_test_preds, int iter);
    void FindMeanStdAllFields (int iter);
	void FindMeanStdTestData(int iter);
    void setHyperParam(int epsilon, double lambdaV, double lambdaU, double momentum, int numFeat, int cvInd);
	void intialize_hierarchy();
};


inline int randWrapper(const int n) { return floor(unif_rand()*n); }

#endif /* defined(____HPMF__) */
