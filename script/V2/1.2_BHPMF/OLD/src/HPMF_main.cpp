//
//  HPMF_main.cpp
//  
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#include <iostream>

#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <time.h>
#include <Rdefines.h>

#include "HPMF.h"
#include <math.h> 
using namespace std;

int XDemo(int numSamples, double *allRMSE, double *allRMSE2,  int datasetId, double &rmse, SEXP args, const char* meanFilePath, const char* stdFilePath, SEXP numNodesPerLevelR) {
    
    int testLen = -1, numTraits = 13, numFeat = 15, taxonomyLevels = 4, tunedLevel, cvInd= 0, iterations = 1, gaps, burn;
    int epsilon = 50; bool cvFlag;
    double momentum = 0, lambdaU = 0.01, lambdaV = 0.01;
    int opt = 2;  //Gibbs Sampling or SGD
    int tune = 0;
    int saveFlag = 0;
	int used_num_hierarchy_level;
    int outWholeFlag = 0; 
    int numFolds = 2;

	double *numNodesPerLevel =  REAL(numNodesPerLevelR);

    burn = INTEGER(getListElement(args, "Burn"))[0];
    gaps = INTEGER(getListElement(args, "Gaps"))[0];
	saveFlag = INTEGER(getListElement(args, "SaveFileFlag"))[0];
	outWholeFlag = INTEGER(getListElement(args, "OutWholeFlag"))[0];
    numTraits = INTEGER(getListElement(args, "NumTraits"))[0];
	numFeat = INTEGER(getListElement(args, "NumFeats"))[0];
    taxonomyLevels = INTEGER(getListElement(args, "NumHierarchyLevel"))[0];
    used_num_hierarchy_level = INTEGER(getListElement(args, "UsedNumHierarchyLevel"))[0];
	tunedLevel = INTEGER(getListElement(args, "PredictLevel"))[0];
    opt = INTEGER(getListElement(args, "Opt"))[0];
    SEXP env = getListElement(args, "Env");
	SEXP dirEXP = getListElement(args, "InputDir");
    const char *input_dir = CHAR(STRING_ELT(dirEXP, 0));
	cout << "name: " << input_dir << endl;


/*
	int nRow, nCol = 3;
    nRow = INTEGER(getListElement(args, "nrowsTrainData"))[0];

    mat = (double *) R_alloc(nRow * nCol, sizeof(double));
    mat = REAL(VECTOR_ELT(args, 0));

    //copy the array into a matrix                                                                                                                                                                                                         
    matVecFloat training_data; // = (float**) R_alloc(nRow, sizeof(float*));
	training_data.resize(nRow);
    for (int ii = 0; ii < nRow; ii++) {
        training_data[ii].resize(nCol); //(float*) R_alloc(nCol,sizeof(float));
        for (int jj = 0; jj < nCol; jj++) {
            training_data[ii][jj] = (float) mat[ii*nCol+jj];
		}
    }
*/  
//    char *dataPath = (char *) R_alloc(100, sizeof(char));
//    dataPath = "../data";

//    tunedLevel = taxonomyLevels + 1;
    int numPlants;
    GetRNGstate();

	cout << meanFilePath << endl;
    HPMF *obj = new HPMF(input_dir, numTraits, numFeat, taxonomyLevels, tunedLevel, used_num_hierarchy_level, cvInd, datasetId, iterations, epsilon, momentum, lambdaV, lambdaU, opt, saveFlag, outWholeFlag, gaps, burn, numSamples, meanFilePath, stdFilePath, numNodesPerLevel, env);
       
    //Hyperparameter tuning for HPMF with SGD
    if (tune) {
      int epsilonSet[] = {50, 40, 30, 20, 10};
      double momentumSet[] = {0.7, 0.6, 0.5 };
      double lambdaVSet[] = {0.01, 0.001};
      double numFeatSet[] = {10, 15, 30, 50};
       
      double minTestErr = 10000;  
        for (int epsId = 0; epsId < 5; epsId++) {
            for (int momId = 0 ; momId < 3; momId++) {
                for (int lambId = 0; lambId < 2; lambId++) {
                    for (int featId = 0; featId < 4; featId++) {
                        double avgTestErr = 0;
                        for (cvInd = 1; cvInd <= numFolds; cvInd++) {
                            
			  obj->setHyperParam(epsilonSet[epsId], lambdaVSet[lambId], lambdaVSet[lambId], momentumSet[momId], numFeatSet[featId], cvInd);
			  double testErr = obj->runHPMF(false, numPlants);
                            avgTestErr += testErr;
                        }
                        
                        if (avgTestErr < minTestErr) {
                            minTestErr = avgTestErr;
                            epsilon = epsilonSet[epsId];
                            momentum = momentumSet[momId];
                            lambdaU = lambdaVSet[lambId];
                            lambdaV = lambdaVSet[lambId];
                            numFeat = numFeatSet[featId];
                        }
                    }
                }
            }
        }
    }
    
    //find the best latent size for Gibbs sampling
    if(tune == 1){
    if (opt == 2) {
        
    int ss = (numTraits - 5)/5;
	double numFeatSet[ss+4];
   
	
	int ind = 0;
	for(int ii = 5; ii < numTraits+6; ii+=5) numFeatSet[ind++] =  ii;

	//        double numFeatSet[] = {10, 15, 30, 50};
        
        double minTestErr = 10000;
        double testErr;
        
        for (int featId = 0; featId < ind; featId++) {
            double avgTestErr = 0;
            for (int cvInd = 1; cvInd <= numFolds; cvInd++) {
                obj->setHyperParam(epsilon, lambdaV, lambdaU, momentum, numFeatSet[featId], cvInd);
                int numSam = floor(numSamples/3);
                testErr = obj->runHPMF(numSam);
                avgTestErr += testErr;
            }
            
            if (avgTestErr < minTestErr) {
                minTestErr = avgTestErr;
                numFeat = numFeatSet[featId];
            }


        }
    }
    }        

 
    cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << numFeat << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
    obj->setHyperParam(epsilon, lambdaV, lambdaU, momentum, numFeat, 0);
	cout << "set the hyper param: " << endl;
    testLen = 0;
    if (opt == 1)
      obj->runHPMF(true, numPlants);
    else{
		rmse = obj->runHPMF(allRMSE, allRMSE2, numPlants);
//		cout << "done!" << endl;
//		cout << "rmse: " << err << endl;
    }

//	cout << "rmse in Demo: " << rmse << endl;

    PutRNGstate();

//   cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
//    cout << numPlants << endl;
//    cout << "++++++++++++++++++++++++++++++++++++++++" << endl;

    //    return opt;
    return numPlants;
}




extern "C" {

  	SEXP DemoHPMF(SEXP args, SEXP meanFilePathSexp, SEXP stdFilePathSexp, SEXP numNodesPerLevelR){
 
  	    int numSamples = INTEGER(getListElement(args, "NumSamples"))[0];
		int datasetId = INTEGER(getListElement(args, "DatasetId"))[0];
		SEXP env = getListElement(args, "Env");

		const char *meanFilePath = CHAR(STRING_ELT(meanFilePathSexp, 0));
		const char *stdFilePath = CHAR(STRING_ELT(stdFilePathSexp, 0));

		double rmse;
    	double *allRMSE = (double *) R_alloc(numSamples, sizeof(double));
    	double *allRMSE2 = (double *) R_alloc(numSamples, sizeof(double));
	//        int opt = XDemo(numSamples, allRMSE, datasetId, env);
        int numPlants = XDemo(numSamples, allRMSE, allRMSE2,  datasetId, rmse, args, meanFilePath, stdFilePath, numNodesPerLevelR);
		cout << "RMSE: " << rmse << endl;
	        
        //        SEXP result, resultNames;
        //
        //        //make return object
        //        SEXP testPred_r;
        //        int nProtect = 0;
        //        int nResultListObjs = 1;
        //
        //        PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
        //        PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
        //        PROTECT(testPred_r = allocMatrix(REALSXP, numSamples, testLen)); nProtect++;
        //
        //        // double *rans = REAL(testPred_r);
        //        //rans = testPred;
        //
        //        // rans = REAL(testPred_r);
        //        for(int ii = 0; ii < numSamples; ii++) {
        //            for(int jj = 0; jj < testLen; jj++)
        //                REAL(testPred_r)[ii + numSamples*jj] = testPred[ii][jj];
        //        }
        //
        //
        //        //samples
        //        SET_VECTOR_ELT(result, 0, testPred_r);
        //        SET_VECTOR_ELT(resultNames, 0, mkChar("pred"));
        //
        
        //        namesgets(result, resultNames);
        //
        //        //unprotect
        //        UNPROTECT(nProtect);
        //
        //
        //
        //        R_FlushConsole();
        //        R_ProcessEvents();
        //
        //
        //
        //
        //
        //	   return result;
        
        SEXP result, resultNames, pred_r, opt_r, rmse_r, numPlants_r, all_rmse_r, all_rmse_r2;
        
        // Allocating storage space:
        //        cout << numSamples*testLen << endl;
        /*PROTECT(pred_r = NEW_NUMERIC(testLen*numSamples));

	//        double *testPred_r  = new double[numSamples*testLen];
	double *testPred_r  = (double *) R_alloc(numSamples*testLen, sizeof(double));

        testPred_r = NUMERIC_POINTER(pred_r);
        
        for(int ii = 0; ii < numSamples; ii++) {
            //          cout << ii << endl;
            
            for(int jj = 0; jj < testLen; jj++){
                testPred_r[ii + numSamples*jj] = testPred[ii][jj];
            }
	    }
	*/

	int nProtected = 0;
    PROTECT(all_rmse_r = NEW_NUMERIC(numSamples));
    ++nProtected;

    double *allRMSE_r  =  (double *) R_alloc(numSamples, sizeof(double));
    allRMSE_r = NUMERIC_POINTER(all_rmse_r);
    for(int ii = 0; ii < numSamples; ii++) {
		allRMSE_r[ii] = allRMSE[ii];
    }      

    PROTECT(all_rmse_r2 = NEW_NUMERIC(numSamples));
    ++nProtected;

    double *allRMSE_r2  =  (double *) R_alloc(numSamples, sizeof(double));
    allRMSE_r2 = NUMERIC_POINTER(all_rmse_r2);
    for(int ii = 0; ii < numSamples; ii++) {
		allRMSE_r2[ii] = allRMSE2[ii];
    }      
       
        // PROTECT(opt_r = NEW_INTEGER(1));
	// ++nProtected;        
	// int *optMeth_r = (int *) R_alloc(1, sizeof(int));
        // optMeth_r = INTEGER_POINTER(opt_r);        
        // optMeth_r[0] = opt;

    PROTECT(numPlants_r = NEW_INTEGER(1));
	++nProtected;        
	int *numP_r = (int *) R_alloc(1, sizeof(int));
    numP_r = INTEGER_POINTER(numPlants_r);        
    numP_r[0] = numPlants;

    PROTECT(rmse_r = NEW_NUMERIC(1));
	++nProtected;        
	double *RMSE_r = (double *) R_alloc(1, sizeof(double));
    RMSE_r = NUMERIC_POINTER(rmse_r);        
    RMSE_r[0] = rmse;


	//	cout << "done init" << endl;

	//        char *names[2] = {"allRMSE", "opt"};
	char *names[4] = {"allRMSE", "numPlants", "RMSE", "allRMSE2"};
        
    PROTECT(resultNames = allocVector(STRSXP, 4));
	++nProtected;        

    for(int i = 0; i < 4; i++)
		SET_STRING_ELT(resultNames,i,mkChar(names[i]));

    //	SET_STRING_ELT(list_names,i,mkChar(names[i]));
        
    PROTECT(result = allocVector(VECSXP, 4));
	++nProtected;

    // attaching myint vector to list:
	//        SET_VECTOR_ELT(result, 0, pred_r);
    SET_VECTOR_ELT(result, 0, all_rmse_r);
    SET_VECTOR_ELT(result, 1, numPlants_r);
    SET_VECTOR_ELT(result, 2, rmse_r);
    SET_VECTOR_ELT(result, 3, all_rmse_r2);
    setAttrib(result, R_NamesSymbol, resultNames);
    UNPROTECT(nProtected);
        
    return result;
         
  }

  void writeFilledGap(SEXP args, SEXP inName, SEXP outName){

    int numSamples = INTEGER(getListElement(args, "numSamples"))[0];
    int numTraits = INTEGER(getListElement(args, "numTraits"))[0];
    int numPlants = INTEGER(getListElement(args, "numPlants"))[0];    
    int gap = INTEGER(getListElement(args, "gap"))[0];    
    int burn = INTEGER(getListElement(args, "burn"))[0];    
    //    char *name = STRING_ELT(getListElement(args, "fileName"),0)[0];
    const char *inFilePath = CHAR(STRING_ELT(inName, 0));
    const char *outFilePath = CHAR(STRING_ELT(outName, 0));

    cout << "inputFile: " << inFilePath << endl;
    cout << "outpurFile: " << outFilePath << endl;
    //    char *inFilePath = "../data/output/testPred_dataset2_opt2_NotBin.txt";
    //    char *inFilePath = "../data/output/testPred_dataset2_opt2_NotBinTab.txt";

    //    char *inFilePath = "../data/output/testWrite.txt";
    //    char *outFilePath = "../data/output/GapFilled_opt2.txt";

    //    dlmwriteVec(inFilePath, outFilePath, numPlants, numTraits, gap, burn, numSamples);

    dlmwriteVec(inFilePath, outFilePath, numPlants, numTraits, gap, burn, numSamples);
    cout << "finish writing!" << endl;
  }
    
}
