//
//  HPMF.cpp
//  
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#include "HPMF.h"

bool countMaxV = true;
bool countMaxU = true;

HPMF::HPMF(const char *dataPath, int numTraits, int numFeat, int taxonomyLevels, int tunedLevel, int used_hierarchy_level, int cvInd, int datasetId, int iterations, int epsion, double momentum, double lambdaV, double lambdaU, int opt, int saveFileFlag, int outWholeFlag, int gaps, int burn, int num_samples, const char* meanFilePath, const char* stdFilePath, double* num_nodes_per_level, SEXP env) {
    if(taxonomyLevels > tunedLevel) {
        cerr << "Error: wrong number of parameters: taxonomyLevel should be smaller than tunedLevel" << endl;
        exit(1);
    }
    cout << "Initializing HPMF/BHPMF paramters...";

    /***************************** 
     Hyper-parameters
    *****************************/
    this->epsilon = epsilon;
    this->lambdaV = lambdaV;
    this->lambdaU = lambdaU;
    this->momentum = momentum;   
    this->numFeat = numFeat;                   //size of latent factor $k$
	this->gaps_ = gaps;
	this->burn_ = burn;
	this->num_samples_ = num_samples;
	this->num_effective_samples_ = floor((num_samples_ - burn_) / gaps_) + 1;

    /*****************************
     Other parameters      
    *****************************/
    this->numTraits = numTraits;
    this->cvInd = cvInd;                     //used for prediction
    this->datasetId = datasetId;             //the datasetId if using CV to provide average RMSE 
    this->dataPath = dataPath;              //data path 
    this->iterations = iterations;          //number of iterations for HPMF (number of top-down, and down-top learning)
    this->opt = opt;                        //opt = 1 : Gradient Descent, opt = 2: Gibbs sampling
    this->env = env;
    this->saveFileFlag = saveFileFlag;    
    this->outWholeFlag = outWholeFlag;
	cout << meanFilePath << endl;
    this->meanFilePath_ = meanFilePath;
    this->stdFilePath_ = stdFilePath;
	this->num_nodes_per_level_.assign(num_nodes_per_level, num_nodes_per_level+taxonomyLevels+1);

    /*****************************
    Taxonomy information: Phylogenetic Group, Famlily, Genus, Species
    taxonomyLevels = 4 : use all of taxonomy information (Phylogenetic Group + Famlily + Genus + Species)
    taxonomyLevels = 3 : use 3 levels of taxonomy information (Phylogenetic Group + Famlily + Genus )
    taxonomyLevels = 2 : use 2 levels of taxonomy information (Phylogenetic Group + Famlily)
    taxonomyLevels = 1 : use 1 level of taxonomy information (Phylogenetic Group)
    tunedLevel = 5 : prediction at the Plant level
    tunedLevel = 4 : prediction at the Species level
    tunedLevel = 3 : prediction at the Genus level
    tunedLevel = 2 : prediction at the Family level
    tunedLevel = 1 : prediction at the Phylogenetic Group level
    tunedLevel > taxonomyLevels
    *****************************/
    this->taxonomyLevels = taxonomyLevels;   //number of taxonomy level   
    this->tunedLevel = tunedLevel;           //prediction level (eg. if looking at prediction at the plant level = taxonomyLevels+1)
	this->predict_level_ = tunedLevel;
	this->used_hierarchy_level_ = taxonomyLevels;
	this->num_hierarchy_level_ = taxonomyLevels;

//	this->num_parents = num_parents;
//	this->num_children = num_children;
//	this->num_nodes_per_level_ = num_nodes_per_level_;

	/*****************************
     Initializing the hierarchy
	 *****************************/ 
	intialize_hierarchy();

    /*****************************
     Initializing Latent matrices
	*****************************/
    nodePerLevel.resize(taxonomyLevels+1);
    uTree = new NodeLatent**[taxonomyLevels+1];
    intializeTree(uTree, 0);
	cout << "initializing Vtree: " << endl;

    vTree = new NodeLatent**[taxonomyLevels+1];
    intializeTree(vTree, numTraits);
    
    if (opt == 2) {  //gibbs sampling
		cout << "\n Initialized BHPMF using Gibbs Sampling." << endl;
    } else{
		cout << "\n Initialized HPMF using Block Gradient Descent Optimization" << endl;
    }

    /*****************************
    Initializing all_pred                                                                                                                                                      
	*****************************/
    all_preds_mean_.resize(num_nodes_per_level_[taxonomyLevels]);
    all_preds_std_.resize(num_nodes_per_level_[taxonomyLevels]);
    for (int ii = 0; ii < num_nodes_per_level_[taxonomyLevels]; ii++) {
        all_preds_mean_[ii].resize(numTraits);
        all_preds_std_[ii].resize(numTraits);
        for (int jj = 0; jj < numTraits; jj++) {
            all_preds_mean_[ii][jj] = 0;
            all_preds_std_[ii][jj] = 0;
        }
    }


    /************************************
     Initialize test_pred
	*************************************/
    char testVar[100];
    if(cvInd) {
        sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytest%d.txt", dataPath, datasetId, cvInd, predict_level_); //level+1);
    } else {
        sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, predict_level_); //level+1);
    }
    // Rprintf("testVar: %s : %d\n", testVar, numFeat);
    dlmreadVec(testVar, test_data_, num_row_test_, 3);
	test_preds_std_.resize(num_row_test_);
	test_preds_mean_.resize(num_row_test_);

}

/*****************************          
 Free Heap variables       
******************************/
HPMF::~HPMF() {
	for (int ii = taxonomyLevels; ii > 0; ii--) {
		for (int jj = 0; jj < numTraits; jj++) {
			SafeDelete(vTree[ii][jj]);	
		}
		for(int jj = 0; jj < num_nodes_per_level_[ii]; jj++) {
			SafeDelete(uTree[ii][jj]);
		}
		SafeDelete(uTree[ii]);
		SafeDelete(vTree[ii]);
    }
	SafeDelete(uTree);
	SafeDelete(vTree);
}

/*****************************
Initialize the hierarchy_
Update the hierarchy info based
on number of taxnomy to be used
******************************/
void HPMF::intialize_hierarchy() {

/*	num_nodes_per_level_.resize(num_hierarchy_level_+1);
	num_nodes_per_level_[0] = 6;
	num_nodes_per_level_[1] = 358;
	num_nodes_per_level_[2] = 3793;
	num_nodes_per_level_[3] = 14320;
	num_nodes_per_level_[4] = 78300;
*/
	int nRow = accumulate(num_nodes_per_level_.begin(), num_nodes_per_level_.end(), 0);

	// read num_parents info
	char file_path[25];
    sprintf (file_path, "%s/num_parents.txt", dataPath);
	dlmreadVec(file_path, num_parents_);
	// cout << endl << "num_nodes at num_parents: " << num_parents_.size() << endl;

	// read num_children info
    sprintf (file_path, "%s/num_children.txt", dataPath);
	dlmreadVec(file_path, num_childern_);
	// cout << "num_nodes at num_children: " << num_childern_.size() << endl;

	// read the hierarchy info
    sprintf (file_path, "%s/hierarchy_info.txt", dataPath);
	dlmreadVec(file_path, hierarchy_, nRow, 1);// num_parents_);

	// cout << hierarchy_[0][0] << "\t" << hierarchy_[2][0] << "\t" << hierarchy_[4][0] << endl;
	// cout << "nRow: " << nRow << endl;

	int last_ind = 0;

	if ((used_hierarchy_level_ == predict_level_-1) ||
		(used_hierarchy_level_ == num_hierarchy_level_)) return;

	cout << "changed hierachy info: " << endl;
	// vector<int> new_num_childern;
	// new_num_childern.resize(num_nodes_per_level_[used_hierarchy_level_]);

	int parent_level_start_ind = accumulate(
		num_nodes_per_level_.begin(),
		num_nodes_per_level_.begin()+used_hierarchy_level_, 0);

	int used_num_hierarchy_level_start_ind =
		parent_level_start_ind - num_nodes_per_level_[used_hierarchy_level_-1];

	cout << parent_level_start_ind << "\t" << used_num_hierarchy_level_start_ind << "\t" << num_nodes_per_level_[used_hierarchy_level_-1] << endl;
	vector<bool> flag_if_first_seen(nRow, true);

	// update parents of nodet at lower levels to used_hierarchy_level_
	for (int level = used_hierarchy_level_+1; level < predict_level_; level++) {
		int start_ind = parent_level_start_ind+num_nodes_per_level_[level-1];
		int last_ind = start_ind +num_nodes_per_level_[level];

		cout << level << "\t" << start_ind << "\t" << last_ind << endl;

		for (int node_ind = start_ind; node_ind < last_ind; node_ind++) {
			vector<int> new_parents_list;
			// int parent;
			for (int parent_ind = 0; parent_ind < num_parents_[node_ind]; parent_ind++) {
				int parent = hierarchy_[node_ind][parent_ind];				
				vector<int> grant_parent = hierarchy_[parent_level_start_ind+parent-1];
				new_parents_list.insert(new_parents_list.end(), grant_parent.begin(), grant_parent.end()); //point to parents of your parents
//				cout << node_ind-start_ind << "\t" << parent << "\t" << grant_parent[0] << "\t" << parent_level_start_ind+parent-1 << endl;
			}
			hierarchy_[node_ind] = new_parents_list;
			num_parents_[node_ind] = new_parents_list.size();
//			cout << new_parents_list[0] << "\t" << new_parents_list.size() << "\t" << hierarchy_[node_ind][0] << "\t" << num_parents_[node_ind] << endl;

		}
		parent_level_start_ind = start_ind;
	}

	parent_level_start_ind = accumulate(num_nodes_per_level_.begin(), num_nodes_per_level_.begin()+predict_level_-2, 0);
//	used_num_hierarchy_level_start_ind = parent_level_start_ind - num_nodes_per_level_[used_hierarchy_level_-1];

	// update num_childern at used_num_hierarchy_level to the num_children at parent of predict_level
//	cout << parent_level_start_ind << "\t" << parent_level_start_ind+num_nodes_per_level_[predict_level_-2] << endl;
	for (int node_ind = parent_level_start_ind; node_ind < parent_level_start_ind+num_nodes_per_level_[predict_level_-2]; node_ind++) {
		for (int parent_ind = 0; parent_ind < num_parents_[node_ind]; parent_ind++) {
			int parent = hierarchy_[node_ind][parent_ind];
			int grant_parent_ind = used_num_hierarchy_level_start_ind + parent - 1;
			if (flag_if_first_seen[grant_parent_ind]) {
				num_childern_[grant_parent_ind] = num_childern_[node_ind];
				flag_if_first_seen[grant_parent_ind] = false;
			}
			else {
				num_childern_[grant_parent_ind] += num_childern_[node_ind];
			}
//			cout << node_ind-parent_level_start_ind << "\t" << parent << "\t" << grant_parent_ind << "\t" << num_childern_[grant_parent_ind] << endl;

		}
	}

//	num_childern_[used_hierarchy_level_] = new_num_childern;
}

// void HPMF:intializeHierarchy() {

// 	// read the hierarchy info
// 	char file_path[25];
//     sprintf (file_path, "%s/hierarchy_info.txt", dataPath);
// 	dlmreadVec(file_path, hierarchy_, num_nodes_per_level_, num_parents_);

// 	vector<int> new_num_childern;
// 	new_num_childern.resize(num_nodes_per_level_[used_num_hiearchy_level_];
// 	for (int level = used_num_hierarchy_level_-2; level >= predict_level_; level--) {
// 		for (int node_ind = 0; node_ind < num_nodes_per_level_[level]; node_ind++) {
// 			set<int> new_parents_list;
// 			for (int parent_ind = 0; parent_ind < num_parents_[level][node_ind]; parent_ind++) {
// 				int parent = hierarchy_[level][node_ind][parent_ind];
// 				vector<int> tmp = hierarchy_[level+1][parent];
// 				new_parents_list.insert(new_parents_list.end(), tmp.begin(), tmp.end());				
// 			}
// //			std::vector<int>::iterator it;
// //			it = std::unique (new_parents_list.begin(), new_parents_list.end());   // 10 20 30 ?  ?  ?  ?
// //			new_parents_list.resize(std::distance(new_parents_list.begin(),it)); // 10 20 30			
// 			hierarchy_[level][node_ind] = new_parents_list;
// 			num_parents_[level][node_ind] = new_parents_list.size();

// 			// update num_child at used_num_hierarchy_level
// 			if (level != predict_level_-1) {
// 				continue;
// 			}
			
// 			for (int parent_ind = 0; parent_ind < num_parents_[level][node_ind]; parent_ind++) {
// 				int parent = hierarchy_[level][node_ind][parent_ind];
// 				new_num_childern[parent] += num_childern_[level][node_ind];
// 			}
// 		}
// 	}

// 	num_childern_[used_num_hiearchy_level_] = new_num_childern;
// }


/*****************************
 Intializing Latent Matrix represented 
 as tree to capture the parent and child info as well 
*****************************/
void HPMF::intializeTree(NodeLatent*** tree, int numNodes) {

	int last_level_ind = 0;
    for (int level = 0; level <=  taxonomyLevels; level++) {

//        Rprintf("\nlevel: %i \n", level);      
        int parLevel, lev;
        if (level == taxonomyLevels){
			for (int ii = level; ii < predict_level_-1; ii++) {
				if (!numNodes)
					last_level_ind += num_nodes_per_level_[ii];
			} 
        lev = tunedLevel-1;
            parLevel = taxonomyLevels-1;
        } else {
            lev = level;
            parLevel = level - 1;
        }
        
        int nRow, nCol = 3;
		matVecInt phyInfo;
        /*****************************
         Reading hierarchy informaiton  
		*****************************/
        if (!numNodes) {
//			nRow = num_nodes_per_level_[level];
            char filePath[25];
            sprintf (filePath, "%s/phyinfo%d_%d", dataPath, lev+1, parLevel+1);
//			dlmreadVec(filePath, phyInfo, nRow, nCol);  // each row is a list of: first number is number of childern, second number is number of
			// parents, rest of numbers are indexes of parents    
			nRow = num_nodes_per_level_[lev];
			nodePerLevel[level] = nRow;
        } else {
            nRow = numNodes;
			phyInfo.resize(nRow);
            for (int trait = 0; trait < numNodes; trait++) {
				phyInfo[trait].resize(3);
                phyInfo[trait][0] = (level == taxonomyLevels) ?  0 : 1;
                phyInfo[trait][1] = 1;
                phyInfo[trait][2] = trait+1;
            }
        }

//		cout << "last_level_ind: " << last_level_ind <<"\t" << nRow << endl;
	   
		// for (int node = 0; node < nRow; node++) {
		// 	num_childern_.push_back(phyInfo[node][0]);
		// 	num_parents_.push_back(phyInfo[node][1]);
		// }

        tree[level] = (NodeLatent **) R_alloc(nRow, sizeof(NodeLatent*));
        for (int node = 0; node < nRow; node++) {
			int num_children;
			int num_parents;
			vector<int> parents;
			if (numNodes) {  // treeV
				// cout << level << "\t" << node << endl;
				num_children = (level == taxonomyLevels) ?  0 : 1;
				num_parents = 1;
				parents.push_back(node+1);
			} else { //treeU
				num_children = (level == taxonomyLevels) ?  0 : num_childern_[node+last_level_ind];
				num_parents = num_parents_[node+last_level_ind];
				parents = hierarchy_[node+last_level_ind];
			}

//			cout << node << "\t" << num_children << "\t" << num_parents << "\t" << parents[0] << endl;

			
            if (!level){  // if it is at the root of tree with no parent
//				tree[level][node] = new NodeLatent(phyInfo[node][0], numFeat, level);
				tree[level][node] = new NodeLatent(num_children, numFeat, level);
			} else {
//                int numChild = (!numNodes && level == taxonomyLevels) ?  0 : phyInfo[node][0];
//                int numChild = (!numNodes && level == taxonomyLevels) ?  0 : num_childern_[node+last_level_ind];
//				if(level == 4) {
//					if(numChild != 0)
//						cout << "node: " << numChild << "\t" << numNodes << endl; 
//				}
//				tree[level][node] = new NodeLatent(phyInfo[node][1], numChild, phyInfo[node][2], tree, numFeat, level, parLevel); 
  				tree[level][node] = new NodeLatent(num_parents, num_children, parents[0], tree, numFeat, level, parLevel);
          }
        }
		if (!numNodes) last_level_ind += num_nodes_per_level_[level];
    }
}

void HPMF::setObservation() { // ????????????? If it is for gibbs training, can we use validation data here
    for (int level = 0; level <= taxonomyLevels; level++) {
        int lev;
        if (level == taxonomyLevels) {
            lev = tunedLevel-1;
		} else {
            lev = level;
        }
        // cout << "level "<< lev+1 <<" * trait \n";
        
        /************************************ 
         reading the training data
         check for cross validation or learning
        *************************************/ 
        char trainVar[100], valVar[100];
        
        if (cvInd > 0) {
			sprintf (trainVar, "%s/fold%d/Tunning/cv%d/Ytrain%d.txt", dataPath, datasetId, cvInd, lev+1);  
        } else if (cvInd == 0) {
			sprintf (trainVar, "%s/fold%d/Ytrain%d.txt", dataPath, datasetId, lev+1);
        }
		cout << "trainVar: " << trainVar << endl;

        int nRow, nCol = 3;
		matVecFloat trainData;
        dlmreadVec(trainVar, trainData, nRow, nCol);

        lenTrainData = nRow;

        for (int ind = 0; ind < lenTrainData; ind++) {   
            int ii = (int) trainData[ind][0]-1;
            int jj = (int) trainData[ind][1]-1;
            float obs = trainData[ind] [2];
            
            uTree[level][ii]->set_observ(obs, jj);
            vTree[level][jj]->set_observ(obs, ii);
        }
    }
}


// void HPMF::setObservation() { // ????????????? If it is for gibbs training, can we use validation data here

//         char trainVar[100], numObsVar[100];
        
// //        if (cvInd > 0) {
// //			sprintf (trainVar, "%s/fold%d/Tunning/cv%d/Ytrain%d.txt", dataPath, datasetId, cvInd, lev+1);  
// //        } else if (cvInd == 0) {
// 			sprintf (trainVar, "%s/YtrainAll.txt", dataPath);
// 			//       }
// 		cout << "trainVar: " << trainVar << endl;

//         int nRow, nCol = 3, lenTrainData;
// 		matVecFloat trainData;
//         dlmreadVec(trainVar, trainData, nRow, nCol);

// //		sprintf(numObsVar, "%s/lenTrainDatas.txt", dataPath);

// 		vector<int> num_train_data_per_level;
// 		num_train_data_per_level.resize(5);
// 		num_train_data_per_level[0] = 63;
// 		num_train_data_per_level[1] = 2308;
// 		num_train_data_per_level[2] = 13776;
// 		num_train_data_per_level[3] = 32987; 
// 		num_train_data_per_level[4] = 120366;
// //		dlmreadVec(numObsVar, lenTrainDatas);

// 	int last_level_ind = 0;
//     for (int level = 0; level <= taxonomyLevels; level++) {
//         int lev;
//         if (level == taxonomyLevels) {
//             for (int ii = level; ii < predict_level_-1; ii++) {
//                 last_level_ind += num_train_data_per_level[ii];
//             }
//             lev = tunedLevel-1;
// 		} else {
//             lev = level;
//         }
//         cout << "level "<< lev+1 <<" * trait \n";
        
//         /************************************ 
//          reading the training data
//          check for cross validation or learning
//         *************************************/ 

//         lenTrainData = num_train_data_per_level[lev];

// 		cout << lenTrainData << endl;
//         for (int ind = 0; ind < lenTrainData; ind++) {   
//             int ii = (int) trainData[last_level_ind+ind][0]-1;
//             int jj = (int) trainData[last_level_ind+ind][1]-1;
//             float obs = trainData[last_level_ind+ind] [2];
            
//             uTree[level][ii]->set_observ(obs, jj);
//             vTree[level][jj]->set_observ(obs, ii);
//         }

// 		last_level_ind += lenTrainData;
//     }
// }

/*****************************                   
 Setting Hyper-parameters
*****************************/

void HPMF::setHyperParam(int epsilon, double lambdaV, double lambdaU, double momentum, int numFeat, int cvInd) {
	this->epsilon = epsilon;
	this->lambdaV = lambdaV;
	this->lambdaU = lambdaU;
	this->momentum = momentum;
	this->numFeat = numFeat;                  //size of latent factor $k$                                      
	this->cvInd = cvInd;

	if(opt == 2) {
		setObservation();
	}
}

double HPMF::runHPMF(bool saveFlag, int &noPlants) {	
	noPlants =  nodePerLevel[predict_level_-1];
	if(opt == 1) {  //grdient descent    
//    if (taxonomyLevels == 0) {
//        
//        int level = tunedLevel-1;
//        
//        Rprintf("level: %i \n", level+1);
//    //    cout << "level " << level+1 << " * trait" << endl;
//        
//        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        //reading the training data
//        //check for cross validation or learning
//        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        char trainVar[100], testVar[100], valVar[100];
//        
//        if (cvInd > 0){
//            
//            sprintf (trainVar, "%s/cv%d/Ytrain%d", dataPath, cvInd, level+1);
//            sprintf (testVar, "%s/cv%d/Ytest%d", dataPath, cvInd, level+1);
//            sprintf (valVar, "%s/cv%d/Yv%d", dataPath, cvInd, level+1);
//            
//        }
//        else if (cvInd == 0){
//            
//            sprintf (trainVar, "%s/Ytrain%d", dataPath, level+1);
//            sprintf (testVar, "%s/Ytest%d", dataPath, level+1);
//            sprintf (valVar, "%s/Yv%d", dataPath, level+1);
//            
//        }
//        
//        int nRow, nCol = 3;
//        float **trainData = dlmread2(trainVar, nRow, nCol, '\t', 1);
//        lenTrainData = nRow;
//        
//        float **valData = dlmread2(valVar, nRow, nCol, '\t', 1);
//        lenValData = nRow;
//        
//        
//        numPlants = nodePerLevel[level];
//        
//        
//        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        //train hpmf
//        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        hpmfTrain(level, trainData, valData);
//        
//        
//        for (int ii = 0; ii < lenTrainData ; ii++)
//            safeDelete(trainData[ii]);
//        safeDelete(trainData);


//        
//        for (int ii = 0; ii < lenValData; ii++) {
//            safeDelete(valData[ii]);
//        }
//        safeDelete(valData);
//
//    }
//    
//    else{
//    
        for (int ii = 0; ii < iterations; ii++) {
            Rprintf("***************** \n" );
            Rprintf("iter: %i \n", ii+1);
            
            hpmfTopdown();
            hpmfDowntop();
			Rprintf("Test Err: %f \n", findTestErr()); 
        }
//    }
	}

	double testErr;
	if(saveFlag) {
		FILE * outFile;                         
		char filePath[100];
		sprintf (filePath, "%s/output/testPred_dataset%d_opt%d.txt", dataPath, datasetId, opt);
		cout << "output is written in " << filePath << endl;
		outFile = fopen(filePath, "w"); 
		testErr =  findTestErr(outFile);
		fclose(outFile);
	} else {
		testErr =findTestErr();
	}
	return testErr;
}

double HPMF::runHPMF(double *allRMSE, double *allRMSE2, int &noPlants) {
	cout << "predict_level_: " << predict_level_ << endl;
	noPlants =  num_nodes_per_level_[predict_level_-1];
//	cout << noPlants << endl;
//	FILE * outFile;	
//	if (saveFileFlag) {
		// char filePath[100];
		// sprintf (filePath, "%s/output/testPred_dataset%d_opt%d_feat_%d.txt", dataPath, datasetId, opt, numFeat);
		// cout << "output is written in " << filePath << endl;
		// outFile = fopen(filePath, "w");
//	}


	int testLen;
	double RMSE;


	char testVar[20];
    if(cvInd) {
        sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytest%d.txt", dataPath, datasetId, cvInd, predict_level_); //level+1);                  
    } else {
        sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, predict_level_); //level+1);                                      
    }
    // Rprintf("testVar: %s : %d\n", testVar, numFeat);

    int len_test_data, nCol = 3;
    matVecFloat test_data;
    dlmreadVec(testVar, test_data, len_test_data, nCol);

	vector<double> sum_test_preds(len_test_data, 0);
	if(opt == 2) {
		int count = 1;
		numPlants = nodePerLevel[tunedLevel-1];
		cout << "Test Err: " << findTestErr() << endl;
		for(int iter = 0; iter < num_samples_; iter++) {
			double t1 = omp_get_wtime();
			Rprintf("iter:  %i \n", iter); 
			Rprintf("------TopDown----- \n");
			for (int indLevel = 1; indLevel <= taxonomyLevels ; indLevel++ ) {  // for each level from down to top
		   		int level;
				if ( indLevel == taxonomyLevels || taxonomyLevels == 1) {  // at the lowest level
					level = tunedLevel-1;
				} else {
					level = indLevel;
				}
				
				level = indLevel;
				Rprintf("level %i * trait \n", level+1);
				numPlants = nodePerLevel[level];
				// cout << level << "\t" << numPlants << endl;
				gibbsTrain(level);
			}
		
			Rprintf("------DownTop----- \n");
			for (int level = taxonomyLevels-1; level >= 0 ; level-- ) {   // for each level from down to top
				Rprintf("level %i * trait \n", level+1);
				numPlants = nodePerLevel[level];
				gibbsTrain(level);
			}

            if (iter >= burn_ & (iter-burn_) % gaps_ == 0) {
				if (outWholeFlag) { 
					cout << "fill gaps\t" << count << endl;
					FindMeanStdAllFields(count);
				}
				cout << "predict tests iteration\t" << count << endl;
				FindMeanStdTestData(count);
					// findTestErr(outFile);
				
				count++;
            }
            double err = findTestErr();
            Rprintf("Test Err: %f \n", err);
            allRMSE[iter] = err;
            double t2 = omp_get_wtime();
            Rprintf("it tooks %f seconds \n",  t2 - t1);
		}

        // calculate the average 
		if (outWholeFlag) {
		    FILE* outMeanFile;
			FILE* outStdFile;
			cout << meanFilePath_ << "\t" << stdFilePath_ << endl;
			outMeanFile = fopen(meanFilePath_, "w"); 
			outStdFile = fopen(stdFilePath_, "w");
			for (int ii = 0; ii < nodePerLevel[taxonomyLevels]; ii++) {
				for (int jj = 0; jj < numTraits; jj++) {
					all_preds_std_[ii][jj] /= (count-1);
					fprintf(outMeanFile,"%5.4f\t", all_preds_mean_[ii][jj]);
					fprintf(outStdFile,"%5.4f\t", sqrt(all_preds_std_[ii][jj]));
				}
				fprintf(outMeanFile,"\n");
				fprintf(outStdFile,"\n");
			}
			fclose(outMeanFile);
			fclose(outStdFile);
//		    double testErr =findTestErr();
//			Rprintf("Test Err: %f \n", testErr);
//			return RMSE;
		} // else {
			FILE* outMeanFileTest;
            FILE* outStdFileTest;
            if (saveFileFlag) {
				outMeanFileTest = fopen(meanFilePath_, "w");
				outStdFileTest = fopen(stdFilePath_, "w");
			}
			double sum_res = 0;
			for (int ind = 0; ind < num_row_test_;  ind++) {
				double res = test_preds_mean_[ind] - test_data_[ind][2];
				sum_res += pow(res, 2);
				test_preds_std_[ind] /= (count-1);
                if (saveFileFlag) {
					fprintf(outMeanFileTest,"%5.4f\t", test_preds_mean_[ind]);
					fprintf(outStdFileTest,"%5.4f\t", sqrt(test_preds_std_[ind]));
				}
            }
			if (saveFileFlag) {
				fprintf(outMeanFileTest,"\n");
				fprintf(outStdFileTest,"\n");
				fclose(outMeanFileTest);
				fclose(outStdFileTest);
			}
			double RMSE = sqrt(sum_res/num_row_test_);
		    return RMSE;
			// }
	}
	return -1;
//	cout << sum_test_preds[0] << "\t" << sum_test_preds[1] << "\t" << sum_test_preds[1000] << endl;
//	cout << "RMSE in C: " << RMSE << endl; 
//	return RMSE;
}

double HPMF::runHPMF(int numGibbs) { 
	int testLen;
	if(opt == 2) {
		numPlants = nodePerLevel[tunedLevel-1];
		cout << "Test Err: " << findTestErr() << endl;
		for(int iter = 0; iter < numGibbs; iter++) {
			double t1 = omp_get_wtime();
			printf("iter:  %d \n", iter); 
			printf("------TopDown----- \n");
			for (int indLevel = 1; indLevel <= taxonomyLevels ; indLevel++ ) {  // for each level from down to top
				int level;
				if( indLevel == taxonomyLevels || taxonomyLevels == 1) {  // at the lowest level
					level = tunedLevel-1;
				} else {
					level = indLevel;
				}
				printf("level %i * trait \n", level+1);
				numPlants = nodePerLevel[level];
				gibbsTrain(level);
			}
		
			printf("------DownTop----- \n");
			for (int level = taxonomyLevels-1; level >= 0 ; level-- ) {  // for each level from down to top
				printf("level %i * trait \n", level+1);
				numPlants = nodePerLevel[level];
				gibbsTrain(level);
			}
            
			double err =  findTestErr();
			cout << "Test Err: " << err << endl;
            double t2 = omp_get_wtime();
            cout << "it tooks " << t2 - t1 << "seconds" << endl;
		}	
	}
	double testErr =findTestErr();
    Rprintf("Test Err: %f \n", testErr);
    return testErr;
}

void HPMF::LoadData(int level, matVecFloat* trainData, matVecFloat* valData) {
	Rprintf("level %i * trait \n", level+1);
	/************************************ 
    reading the training data
    check for cross validation or learning
	*************************************/ 
	char trainVar[100], testVar[100], valVar[100];
        
	if (cvInd > 0) {   
		sprintf (trainVar, "%s/fold%d/Tunning/cv%d/Ytrain%d.txt", dataPath, datasetId, cvInd, level+1);
		sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytestn%d.txt", dataPath, datasetId, cvInd, level+1);
		sprintf (valVar, "%s/fold%d/Tunning/cv%d/Yv%d.txt", dataPath, datasetId, cvInd, level+1);
	} else if (cvInd == 0) {
		sprintf (trainVar, "%s/fold%d/Ytrain%d.txt", dataPath, datasetId, level+1);
		sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, level+1);
		sprintf (valVar, "%s/fold%d/Yv%d.txt", dataPath, datasetId, level+1);
	}

	cout << epsilon << "\t" <<  momentum << "\t" << lambdaV << "\t" << numFeat << endl;
	cout << "trainVar: " << trainVar << endl;
         
	int nRow, nCol = 3;
	dlmreadVec(trainVar, *trainData, nRow, nCol);
	lenTrainData = nRow;
                
	dlmreadVec(valVar, *valData, nRow, nCol);

	lenValData = nRow;
        
	numPlants = nodePerLevel[level];
	cout << "numPlants: " << numPlants << endl;
}

void HPMF::hpmfDowntop() {
    mode = 2; // down-top mode
    Rprintf("------DownTop----- \n");
    
    for (int level = taxonomyLevels-1; level >= 0 ; level-- ) {  // for each level from down to top
		matVecFloat trainData;
		matVecFloat valData;		
		LoadData(level, &trainData, &valData);  // load the data
		hpmfTrain(level, trainData, valData);  // train HPMF
    }
}

void HPMF::hpmfTopdown() {
    mode = 1;  // top-down mode
    Rprintf("------TopDown-----\n");
    
    for (int indLevel = 1; indLevel <= taxonomyLevels ; indLevel++ ) {  // for each level from down to top
        int level;
        if( indLevel == taxonomyLevels || taxonomyLevels == 1) {  // at the lowest level
            level = tunedLevel-1;
		} else {
            level = indLevel;
        }
		matVecFloat trainData;
		matVecFloat valData;
		LoadData(level, &trainData, &valData);  // load the data        
        hpmfTrain(level, trainData, valData);  // train HPMF
    }
}

void HPMF::Predict(int ind, int level, double* err, double *uu, double *vv, matVecFloat data) {
	const int inc_one = 1;
	int ii  = (int) data[ind][0]-1;
	int jj  = (int) data[ind][1]-1;
	float obs = data[ind][2];
	cout << ii << "\t" << jj << "\t" << obs << endl;
	uu = uTree[level][ii]->get_latent_factor();  // uMat[ii]->u;  //assuming           
	vv = vTree[level][jj]->get_latent_factor();  // vMat[jj]->v;  //assuming
	for (int i = 0; i < numFeat; i++)
		cout << uu[i] << "\t";
	cout << endl;
 
	double pred = F77_NAME(ddot)(&numFeat, uu, &inc_one, vv, &inc_one);  // pred = uu .* vv (inner product)    
	*err = pred - obs;
	cout << pred << "\t" << *err << endl;
}

void HPMF::hpmfTrain(int level, matVecFloat trainData, matVecFloat valData) {
    const int inc_one = 1;
	cout << "hpmf Train: " << endl;
//	lenTrainData = 10;
    int batch_size = 10000; // number training triplets per batch
    int num_batches = ceil( (double) lenTrainData / batch_size);
    int max_epoch = 50;
    double  old_err = -1;    
    std::unique_ptr<int[]> freq_u(new int[numPlants]);
    std::unique_ptr<int[]> freq_v(new int[numTraits]);
    vector<int> shuffle_idx(lenTrainData);
    for (int ii = 0; ii < lenTrainData; ii++) {
		shuffle_idx[ii] = ii;
    }
	double Inf = pow(2,63)-1;
	// update the prior for all the nodes.
    for (int node = 0; node < numPlants; node++) {
		uTree[level][node]->UpdatePrior(mode);
    }
    for (int node = 0; node < numTraits; node++) {
		vTree[level][node]->UpdatePrior(mode);
    }
    
//	feenableexcept(FE_INVALID | FE_OVERFLOW);
    for (int epoch = 0; epoch < max_epoch; epoch++) {
		random_shuffle(shuffle_idx.begin(), shuffle_idx.begin()+lenTrainData);  //??????
		int cN;
		for (int batch = 0; batch < num_batches; batch++) {
			int last_ind;
			/************************************ 
            determin the entries to be passed at the current batch
			*************************************/ 
			if ((batch+1) * batch_size <= lenTrainData) {
				cN = batch_size;
				last_ind = (batch+1) * batch_size;
			} else {
				last_ind = lenTrainData;
				cN = lenTrainData - (batch) * batch_size;
			}
			fill_n(freq_u.get(), numPlants, 0);
			fill_n(freq_v.get(), numTraits, 0);

			/************************************ 
            Update Gradient for nodes in current batch
			*************************************/
			for (int ind = batch * batch_size; ind < last_ind; ind++) {   
				int obs_ind = shuffle_idx[ind];
				double err;
				double *uu;
				// 	cout << ind << endl;
//    std::unique_ptr<double[]> uu(new double[numFeat]);
//				cout << uu[0] << "\t" << uu[1] << endl;
				double *vv;
				int ii  = (int) trainData[obs_ind][0]-1;
				int jj  = (int) trainData[obs_ind][1]-1;
//				Predict(obs_ind, level, &err, uu.release(), vv, trainData);
//				cout << uu[0] << "\t" << uu[1] << endl;

				float obs = trainData[obs_ind][2];

				uu = uTree[level][ii]->get_latent_factor();  // uMat[ii]->u;  //assuming           
				vv = vTree[level][jj]->get_latent_factor();  // vMat[jj]->v;  //assuming

				double pred = F77_NAME(ddot)(&numFeat, uu, &inc_one, vv, &inc_one);  // pred = uu .* vv (inner product)  
				err = pred - obs;
				uTree[level][ii]->UpdateGradient(freq_u[ii], err, vv);  // dU[ii] += err*v[jj]
				vTree[level][jj]->UpdateGradient(freq_v[jj], err, uu);  // dV[jj] += err*u[ii] 
 
				 
				double maxU = -Inf;
				double minU = Inf;
				for (int iu = 0; iu < numFeat ; iu ++) {
					if (uu[iu] > maxU)  maxU = uu[iu];
					if (uu[iu] < minU)  minU = uu[iu];
				}
				double maxV = -Inf;
				double minV = Inf;
                for (int iu = 0; iu < numFeat ; iu ++) {
					if (vv[iu] > maxV)  maxV = vv[iu];
                    if (vv[iu] < minV)  minV = vv[iu];
                }
/*
                if ((maxU > 6 || minU < -6) && (maxV < 6 ||minV > -6)){
                    cout << "U becomes bigger!" << endl;
                    for (int iu = 0; iu < numFeat ; iu ++)
                        cout << uu[iu] << "\t";
                    cout << endl;
                    for (int iu = 0; iu < numFeat ; iu ++)
                        cout << vv[iu] << "\t";
                    cout << endl;
                    exit(1);
                }

				if ((maxV > 6 || minV < -6) && (maxU < 6 ||minU > -6)){
                    cout << "V becomes bigger!" << endl;
					for (int iu = 0; iu < numFeat ; iu ++)
						cout << uu[iu] << "\t";
                    cout << endl;
                    for (int iu = 0; iu < numFeat ; iu ++)
						cout << vv[iu] << "\t";
                    cout << endl;
                    exit(1);
                }
*/
				freq_u[ii]++;
				freq_v[jj]++;
			}

			/************************************ 
            Update Latent Factors
			*************************************/ 
			for (int jj = 0; jj < numTraits; jj++) {  // loop in hash map in vv 
				if (!freq_v[jj]) {
					continue;  // skip traits with no observation in this batch
				}  
				double coef = lambdaV * freq_v[jj];  // vMat[jj]->prior = \lambda * ( sum(v_pars) + sum(v_chs))
				vTree[level][jj]->UpdateLatentFactor(coef, momentum, (double) epsilon/(batch_size));
			}
            
			for (int ii = 0; ii < numPlants; ii++) {  // loop in hash map in uu  
				if (!freq_u[ii]) {
					continue;  // skip plants with no observation in this batch
				}
				double coef = lambdaU * freq_u[ii];  // uMat[ii]->prior = \lambda * ( sum(u_pars) + sum(u_chs))
				uTree[level][ii]->UpdateLatentFactor(coef, momentum, (double) epsilon/(batch_size));
			}            
		}
                
		/************************************ 
	    Find validation error, early stop if over
		*************************************/ 
		double sumRes = 0;
		lenValData = 1;
		for (int ind = 0; ind < lenValData; ind++) {
			double res, *uu, *vv;
				int ii  = (int) valData[ind][0]-1;
				int jj  = (int) valData[ind][1]-1;
//				cout << "ii: " << ii << "\t" << jj << endl;
				// Predict(obs_ind, level, &err, uu, vv, trainData);
	float obs = valData[ind][2];

	uu = uTree[level][ii]->get_latent_factor();  // uMat[ii]->u;  //assuming           
	vv = vTree[level][jj]->get_latent_factor();  // vMat[jj]->v;  //assuming
 
	double pred = F77_NAME(ddot)(&numFeat, uu, &inc_one, vv, &inc_one);  // pred = uu .* vv (inner product)    
/*	for (int iu = 0; iu < numFeat ; iu ++)
		cout << uu[iu] << "\t";
	cout << endl;
    for (int iu = 0; iu< numFeat ; iu ++)
		cout <<vv[iu] << "\t";
	cout << endl;
*/// cout << obs << "\t" << pred << endl; 
   res = pred - obs;
   
			// Predict(ind, level, &res, uu, vv, valData);           
			sumRes += pow(res, 2);
		}
		sumRes /= lenValData;
		double val_err = sqrt(sumRes);
        
		if (epoch > 5 && val_err > old_err) {
			break;
		}        
		old_err = val_err;
    }
    
    Rprintf("valErr = %f\n", old_err);
}




void HPMF::gibbsTrain(int level) {
	 
	/************************************
    Initial Parameter
	*************************************/
    double sigInv = 1/0.01;
    double sigUInv = 1/0.01;
    double sigVInv = 1/0.01;
       
    //for each plant
    //~ cout << "palnts \n";
    #pragma omp parallel for
    for (int node = 0; node < numPlants; node++){
        uTree[level][node]->GibbsUpdate(vTree, sigInv, sigUInv, level);
    }
	
    //for each trait
    //~ cout << "traits: \n";
    #pragma omp parallel for
    for (int node = 0; node < numTraits; node++){
        vTree[level][node]->GibbsUpdate(uTree, sigInv, sigVInv, level);
    }
}

double HPMF::findTestErr() {
    int level = used_hierarchy_level_; // tunedLevel-1;
    char testVar[100];
    const int incOne = 1;

    if(cvInd)
		sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytest%d.txt", dataPath, datasetId, cvInd, predict_level_); //level+1);
    else
		sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, predict_level_); //level+1);

    Rprintf("testVar: %s\n", testVar);
    
    int nRow, nCol = 3;
    matVecFloat testData;
    dlmreadVec(testVar, testData, nRow, nCol);
    lenTestData = nRow;
    
    /************************************
    Find validation error, early stop if over
    *************************************/
    double sumRes = 0;
    for (int ind = 0; ind < lenTestData; ind++) {
		double res, *uu, *vv;
		// Predict(ind, level, &res, uu, vv, testData);      
						int ii  = (int) testData[ind][0]-1;
				int jj  = (int) testData[ind][1]-1;
				// cout << "ii: " << ii << "\t" << jj << endl;
				// Predict(obs_ind, level, &err, uu, vv, trainData);
	float obs = testData[ind][2];

	uu = uTree[level][ii]->get_latent_factor();  // uMat[ii]->u;  //assuming           
	vv = vTree[level][jj]->get_latent_factor();  // vMat[jj]->v;  //assuming
 
	double pred = F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);  // pred = uu .* vv (inner product)    
   res = pred - obs;
//   cout << ind << "\t" << pred << "\t" << level << endl;
        sumRes += pow(res, 2);
    }
   
    sumRes /= lenTestData;
    double testErr = sqrt(sumRes);   
    return testErr;
}

double HPMF::findTestErr(vector<double> &sum_test_preds, int iter) {
    int level = used_hierarchy_level_; //tunedLevel-1;
    char testVar[100];
    const int incOne = 1;

    if(cvInd) {
        sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytest%d.txt", dataPath, datasetId, cvInd, predict_level_); //level+1)
    } else {
        sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, predict_level_); //level+1);                   
    }
    Rprintf("testVar: %s : %d %d\n", testVar, numFeat, num_effective_samples_);

    int nRow, nCol = 3;
    matVecFloat testData;
    dlmreadVec(testVar, testData, nRow, nCol);
    lenTestData = nRow;

    /************************************                                                                                                                
     Find validation error, early stop if ove
	*************************************/
    double sum_res_sq_per_iter = 0;

    for (int ind = 0; ind < lenTestData;  ind++) {
        int ii  = (int) testData[ind][0]-1;
        int jj  = (int) testData[ind][1]-1;
        float obs = testData[ind][2];

        double *uu = uTree[level][ii]->get_latent_factor();           //uMat[ii]->u;  //assuming                            
        double *vv = vTree[level][jj]->get_latent_factor();           //vMat[jj]->v;  //assuming                            

        double pred =  F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);      //pred = uu .* vv (inner product)
		
		double res = pred - obs;
        sum_res_sq_per_iter += pow(res, 2);

		sum_test_preds[ind] += pred;
		
	}

	if (iter == num_samples_-1) {
		double sum_res_sq_all = 0;
		for (int ind = 0; ind < lenTestData;  ind++) {
			double avg_pred = sum_test_preds[ind] / num_effective_samples_;
			double res = testData[ind][2] - avg_pred;
			sum_res_sq_all += pow(res, 2);
		}
		return sqrt(sum_res_sq_all/lenTestData);
	}

	return sqrt(sum_res_sq_per_iter/lenTestData);
}


//double HPMF::findTestErr(int iter, int &testLen, FILE *outFile){
double HPMF::findTestErr(FILE *outFile) {
    int level = used_hierarchy_level_; //tunedLevel-1;
    char testVar[100];
    const int incOne = 1;
    
    if(cvInd) {
		sprintf (testVar, "%s/fold%d/Tunning/cv%d/Ytest%d.txt", dataPath, datasetId, cvInd, predict_level_); //level+1);  
	} else {
		sprintf (testVar, "%s/fold%d/Ytest%d.txt", dataPath, datasetId, predict_level_); //level+1);  
	}
//    Rprintf("testVar: %s : %d\n", testVar, numFeat);
    
    int nRow, nCol = 3;
    matVecFloat testData;
    dlmreadVec(testVar, testData, nRow, nCol);
    lenTestData = nRow;
    
    //    testPred[iter] = (double *) R_alloc(lenTestData, sizeof(double));

    /************************************
     Find validation error, early stop if over
    *************************************/
    double sumRes = 0;

	double test_preds[lenTestData];

    for (int ind = 0; ind < lenTestData;  ind++) { 
        int ii  = (int) testData[ind][0]-1;
        int jj  = (int) testData[ind][1]-1;
        float obs = testData[ind][2];
        
        double *uu = uTree[level][ii]->get_latent_factor();           //uMat[ii]->u;  //assuming
        double *vv = vTree[level][jj]->get_latent_factor();           //vMat[jj]->v;  //assuming
        
        double pred =  F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);      //pred = uu .* vv (inner product)
		//  testPred[iter][ind] = pred;
		test_preds[ind] = pred;

		if(outWholeFlag){
			fprintf(outFile,"%5.4f\t",pred);
		}


//		fwrite(&pred, sizeof(double), 1, outFile);
        double res = pred - obs;
        
        sumRes += pow(res, 2);
    }    

//	if(outWholeFlag) {
// 	FILE * fp;
//     char filePath[100];
//     sprintf (filePath, "%s/output/testPred_dataset%d_opt%d_feat_%d.bin", dataPath, datasetId, opt, numFeat);
//     cout << "output is written in " << filePath << endl;
//     fp = fopen(filePath, "ab");
// 	fwrite(test_preds, sizeof(double), lenTestData, fp);
// //	}
// 	fclose(fp);

//		for (int ind = 0; ind < 5; ind++) cout << test_preds[ind] << "\t";
//		cout << endl;

/*    if(outWholeFlag){

		 for (int ii = 0; ii < nodePerLevel[taxonomyLevels]; ii++) {
		 for(int jj = 0; jj < numTraits; jj++){

		 //	float obs = testData[ind][2];

		 double *uu = uTree[level][ii]->get_latent_factor();           //uMat[ii]->u;  //assuming                                                                                                                  
		 double *vv = vTree[level][jj]->get_latent_factor();           //vMat[jj]->v;  //assuming                                                                                                                  
		 double pred =  F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);      //pred = uu .* vv (inner product)                                                                                         
		 //  testPred[iter][ind] = pred;                                                                                                                                                                   
		 fprintf(outFile,"%5.4f\t",pred);

		 //double res = pred - obs;

		 //sumRes += pow(res, 2);
		 }
		 }
	}
*/    
    
    fprintf(outFile,"\n");
    sumRes /= lenTestData;
    double testErr = sqrt(sumRes);
    //testLen = lenTestData;   
    
    return testErr;
}


void HPMF::FindMeanStdTestData(int iter) {
	int level = predict_level_-1;
    const int incOne = 1;
    /************************************
     Find validation error, early stop if over
	*************************************/
    double sumRes = 0;
    #pragma omp parallel for
    for (int ind = 0; ind < num_row_test_;  ind++) {
        int ii  = (int) test_data_[ind][0]-1;
        int jj  = (int) test_data_[ind][1]-1;
        float obs = test_data_[ind][2];

        double *uu = uTree[level][ii]->get_latent_factor();           // uMat[ii]->u;
        double *vv = vTree[level][jj]->get_latent_factor();           // vMat[jj]->v;
        double pred =  F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);      // pred = uu.*vv (inner product) 
        if (iter == 1) {
			test_preds_mean_[ind] = pred;
            test_preds_std_[ind] = 0;
         } else {
            float prev_mean = test_preds_mean_[ind];
            test_preds_mean_[ind] += (pred - prev_mean)/iter;
            test_preds_std_[ind] += (pred - prev_mean) * (pred - test_preds_mean_[ind]);
		}
    }
}

void HPMF::FindMeanStdAllFields(int iter) {
    // if(!outWholeFlag) return;

    int level = predict_level_-1;
    const int incOne = 1;

    #pragma omp parallel for
    for (int ii = 0; ii < nodePerLevel[taxonomyLevels]; ii++) {
        for(int jj = 0; jj < numTraits; jj++) {
            double *uu = uTree[level][ii]->get_latent_factor();           // uMat[ii]->u;
            double *vv = vTree[level][jj]->get_latent_factor();           // vMat[jj]->v;

            double pred =  F77_NAME(ddot)(&numFeat, uu, &incOne, vv, &incOne);      // pred = uu.*vv (inner product)
            if (iter == 1) {
                all_preds_mean_[ii][jj] = pred;
                all_preds_std_[ii][jj] = 0;
            } else {
                float prev_mean = all_preds_mean_[ii][jj];
                all_preds_mean_[ii][jj] += (pred - prev_mean)/iter;
                all_preds_std_[ii][jj] += (pred - prev_mean) * (pred - all_preds_mean_[ii][jj]);
            }
        }
    }
}
