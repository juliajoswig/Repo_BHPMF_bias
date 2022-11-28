//
//  utillity.cpp
//  
//
//  Created by Farideh Fazayeli on 7/18/13.
//
//

#include "utillity.h"

void mvrnorm(double *des, double *mu, double *cholCov, int dim, bool upper){

	int i;
	int inc = 1;
	double one = 1.0;
	double zero = 0.0;
	//make some std norm draws
	for(i = 0; i < dim; i++)
		des[i] = rnorm(0.0, 1.0);
    
	//mult this vector by the lower triangle of the cholCov
	if(upper)
		F77_NAME(dtrmv)("U", "T", "N", &dim, cholCov, &dim, des, &inc);
	else
		F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);

	//add the mean to the result
	F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);


}

double getcputime(){
    struct timeval tim;
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    tim=ru.ru_utime;
    double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
    tim=ru.ru_stime;
    t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
    return t;
}

/*
  Print an array, assumes array dimension is known
*/
void PrintArr(double *myArr, int nRow){

    for (int ii = 0; ii < nRow; ii++) {
        
        Rprintf("%5.3f ", myArr[ii]);
        //cout << myArr[ii] << " ";
    }
    Rprintf("\n");
    //cout << endl;
}

void PrintArr(int *myArr, int nRow){
    
    for (int ii = 0; ii < nRow; ii++) {
        Rprintf("%5.3f ", myArr[ii]);
        //cout << myArr[ii] << " ";
    }
    Rprintf("\n");
    //cout << endl;
}

/*
  Print a matrix stored in an array, assumes matrix dimension is known
*/
void printMat(float *myArr, int nRow, int nCol){
    
    for (int ii = 0; ii < nRow; ii++) {
        for (int jj = 0; jj < nCol; jj++) {
            Rprintf("%5.3f ", myArr[ii*nCol+jj]);
            //cout << myArr[ii*nCol+jj] << " ";
        }
        Rprintf("\n");
        //cout << endl;
    }
}

/*
  Print a matrix from a given matrix, assumes matrix dimension is known
*/
void printMat(float **myMat, int nRow, int nCol){
    
    for (int ii = 0; ii < nRow; ii++) {
        for (int jj = 0; jj < nCol; jj++) {
            Rprintf("%5.3f ", myMat[ii][jj]);
            //cout << myMat[ii][jj] << " ";
        }
        Rprintf("\n");
        //cout << endl;
    }
}


void printMat(int **myMat, int nRow, int nCol){
    
    for (int ii = 0; ii < nRow; ii++) {
        for (int jj = 0; jj < nCol; jj++) {
            Rprintf("%5.3f ", myMat[ii][jj]);
            //cout << myMat[ii][jj] << " ";
        }
        Rprintf("\n");
        //cout << endl;
    }
}



SEXP getListElement (SEXP list, char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

float **dlmreadRFloat(SEXP env, char *fileName, int &nRow, int &nCol, char *delim)
{
	int nProtected = 0;
	SEXP argsList, Rfcall, listRout;
	double *mat;

	if (!isEnvironment(env))
		error("'env' should be an environment");

	PROTECT(argsList = allocVector(STRSXP, 2));
	++nProtected;
	SET_STRING_ELT(argsList, 0, mkChar(fileName));
	SET_STRING_ELT(argsList, 1, mkChar(delim));


	PROTECT(Rfcall = lang2(install("readFileNumeric"), R_NilValue));
	++nProtected;

	SETCADR(Rfcall, argsList);
	listRout = eval(Rfcall, env);

	nRow = INTEGER(getListElement(listRout, "nrows"))[0];
	nCol = INTEGER(getListElement(listRout, "ncols"))[0];

	//mat = new double[nRow * nCol];
	mat = (double *) R_alloc(nRow * nCol, sizeof(double)); 
	mat = REAL(VECTOR_ELT(listRout, 0));

	//copy the array into a matrix                                                                                                                                    
	//  float **outMat = new float*[nRow];
	float **outMat = (float**) R_alloc(nRow, sizeof(float*));
	for (int ii = 0; ii < nRow; ii++) {
		//    outMat[ii] = new float[nCol];
		outMat[ii] = (float*) R_alloc(nCol,sizeof(float));
		for (int jj = 0; jj < nCol; jj++) {
			outMat[ii][jj] = (float) mat[ii*nCol+jj];
		}
	}

	//      safeDelete(mat);                                                                                                                                                  
	UNPROTECT(nProtected);
	Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
	return outMat;
}

int **dlmreadRInt(SEXP env, char *fileName, int &nRow, int &nCol, char *delim)
{
	int nProtected = 0;
	SEXP argsList, Rfcall, listRout;

	if (!isEnvironment(env))
		error("'env' should be an environment");

	PROTECT(argsList = allocVector(STRSXP, 2));
	++nProtected;
	SET_STRING_ELT(argsList, 0, mkChar(fileName));
	SET_STRING_ELT(argsList, 1, mkChar(delim));

	PROTECT(Rfcall = lang2(install("readFileInteger"), R_NilValue));
	++nProtected;

	SETCADR(Rfcall, argsList);
	listRout = eval(Rfcall, env);

	nRow = INTEGER(getListElement(listRout, "nrows"))[0];
	nCol = INTEGER(getListElement(listRout, "ncols"))[0];

	int *mat;
	//  mat = new int[nRow * nCol];
	mat = (int*) R_alloc(nRow * nCol, sizeof(int));
	mat = INTEGER(VECTOR_ELT(listRout, 0));

	UNPROTECT(nProtected);
	//copy the array into a matrix                                                                                                                                    
	//  int **outMat = new int*[nRow];
	int **outMat = (int**) R_alloc(nRow, sizeof(int*));
	for (int ii = 0; ii < nRow; ii++) {
		//outMat[ii] = new int[nCol];
		outMat[ii] = (int*) R_alloc(nCol, sizeof(int));
		for (int jj = 0; jj < nCol; jj++) {
			outMat[ii][jj] = mat[ii*nCol+jj];
		}
	}

	//safeDelete(mat);                                                                                                                                                        
	Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
	return outMat;

}

/*                                                                                                                                                                                
																																												  void dlmreadVec(char *filePath, matVecFloat &myMat, int &nRow, int nCol)                                                                                                         
                                                                                                                                                                                  
																																												  Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.                                                                         
*/
void dlmreadVec(char *filePath, matVecFloat &myMat, int &nRow, int nCol){

	vector<float>  myArr;
	FILE * pFile;
	pFile = fopen (filePath, "r");


	if (!pFile)
    {
		cout << "Cannot open file "<< filePath << endl;
		//cerr << "Failed to open " << pFile << endl;                                                                                                                         
		//            exit(1);  //abort program                                                                                                                                           
    }

	cout << "ncol: " << nCol << endl;

	nRow = 0;
	float num;

	while(fscanf(pFile, "%f", &num) != EOF){
		myArr.push_back(num);
	}

	fclose(pFile);
	nRow = myArr.size()/nCol;

	myMat.resize(nRow);
	for(int ii = 0; ii < nRow; ii++){
		myMat[ii].resize(nCol);
		for (int jj = 0; jj < nCol; jj++)
			myMat[ii][jj] = myArr[ii*nCol+jj];
	}
	Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
}


/*                                                                                                                                                                          
																																											void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol)                                                                                                           
                                                                                                                                                                                  
																																											Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.                                                                         
*/
void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol){

	vector<int> myArr;
	FILE * pFile;
	pFile = fopen (filePath, "r");


	if (!pFile)
    {
		cout << "Cannot open file "<< filePath << endl;
		//cerr << "Failed to open " << pFile << endl;                                                                                                                         
		//            exit(1);  //abort program                                                                                                                                           
    }

	cout <<"nCol: " <<  nCol << endl;

	nRow = 0;
	int num;

	while(fscanf(pFile, "%d", &num) != EOF){
		myArr.push_back(num);
	}

	fclose(pFile);

	nRow = myArr.size()/nCol;

	myMat.resize(nRow);
	for(int ii = 0; ii < nRow; ii++){
		myMat[ii].resize(nCol);
		for (int jj = 0; jj < nCol; jj++)
			myMat[ii][jj] = myArr[ii*nCol+jj];
	}

	Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
}



/*                                                                                                                                                                          
																																											void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol)                                                                                                           
                                                                                                                                                                                  
																																											Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.                                                                         
*/
void dlmreadVec(char *filePath, vector<int> &myArr){

	FILE * pFile;
	pFile = fopen (filePath, "r");


	if (!pFile)
    {
		cout << "Cannot open file "<< filePath << endl;
		//cerr << "Failed to open " << pFile << endl;                                                                                                                         
		//            exit(1);  //abort program                                                                                                                                           
    }

	int num;

	while(fscanf(pFile, "%d", &num) != EOF){
		myArr.push_back(num);
	}

	fclose(pFile);

}



/*
	void dlmreadVec(char *filePath, matVecInt &myMat, int nRow, vector<int> nCols)
	Read a text file into a matrix, a very simple code with assumption number of rows and number of columns of each row of the matrix is known.
*/
void dlmreadVec(char *filePath, matVecInt &myMat, int nRow, vector<int> nCols){

	vector<int> myArr;
	FILE * pFile;
	pFile = fopen (filePath, "r");


	if (!pFile)
    {
		cout << "Cannot open file "<< filePath << endl;
		//cerr << "Failed to open " << pFile << endl;                                                                                                                         
		//            exit(1);  //abort program                                                                                                                                           
    }

//	cout <<"nCol: " <<  nCol << endl;

	nRow = 0;
	int num;

	while(fscanf(pFile, "%d", &num) != EOF){
		myArr.push_back(num);
	}

	fclose(pFile);

	myMat.resize(nRow);
	for(int ii = 0; ii < nRow; ii++){
		myMat[ii].resize(nCols[ii]);
		for (int jj = 0; jj < nCols[ii]; jj++)
			myMat[ii][jj] = myArr[ii*nCols[ii]+jj];
	}

//	Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
}





/*																						  
																						  void dlmwriteVec(char *inFilePath, char *outFilePath, int nRow,  int nCol, int gap, int burn, int nSams)                                                                                       
                                                                                                                                                                                  
																						  write the Gibbs sampler output into a text file, a very simple code with assumption number of rows of the matrix is known.                                                                         
*/
void dlmwriteVec(const char *inFilePath, const char *outFilePath, int nRow, int nCol, int gap, int burn, int nSams){

	matVecFloat myMat;

	//  vector<float>  myArr;
	FILE * ipFile;
	ipFile = fopen (inFilePath, "r");

	if (!ipFile)
	{
		cout << "Cannot open file "<< inFilePath << endl;
		//cerr << "Failed to open " << pFile << endl;                                                                                                                         
		//            exit(1);  //abort program
	}

	float num;

	myMat.resize(nRow);
	for(int ii = 0; ii < nRow; ii++)
		myMat[ii].resize(nCol);

	cout << nRow << "\t" << nCol << "\t" << nSams << endl;

	int count = 1;
	//  while(fscanf(pFile, "%f", &num) != EOF){
	for(int sam = 0; sam < nSams; sam++){
		for(int ii = 0; ii < nRow; ii++){
			//  if(!sam && !ii)   myMat[ii].resize(nCol);
			for(int jj = 0; jj < nCol; jj++){
				fscanf(ipFile, "%f", &num);
				//	cout << sam << "\t" << ii << "\t" << jj << "\t" << num << endl; //<< "\t" << myMat[ii][jj]<< endl;
				if(sam == gap){
					myMat[ii][jj] = num;
				}
				else if (sam > gap & (sam-gap) % burn == 0){
					//	  cout << sam << "\t" << ii << "\t" << jj << "\t" << num << "\t" << myMat[ii][jj]<< endl;
					myMat[ii][jj] += num;
					if (!ii & !jj) count++;
				}
			}
		}
	}

	cout << "effective number of samples: " << count << endl;

	fclose(ipFile);
	//  nRow = myArr.size()/nCol;

	FILE * opFile;
	opFile = fopen (outFilePath, "w");

	if (!opFile)
	{
		cout << "Cannot open file "<< outFilePath << endl;                                                                                                                   
		//            exit(1);  //abort program
	}

	for(int ii = 0; ii < nRow; ii++){
		for (int jj = 0; jj < nCol; jj++){
			//  cout  << ii << "\t" << jj << "\t" << myMat[ii][jj]/count << endl;
			fprintf(opFile,"%5.4f \t", myMat[ii][jj]/count);
		}
		fprintf(opFile,"\n");
	}

	fclose(opFile);

}


/*
  double dlmread2(char *filePath, int nRow, int nCol, char delim)
 
  Read a text file into a matrix, assumes matrix dimension is known.
  Stringstream is used for splitting the string.
*/
float **dlmread2(char *filePath, int &nRow, int &nCol, char delim){
    
    //read the whole file into a buffer
    const int BUFFSIZE = 1000;
    
    
    FILE *pfile = fopen(filePath, "r");
	if (!pfile) {
		Rprintf("Cannot open file %s \n", filePath);
		exit( EXIT_FAILURE);
	}
	fseek(pfile, 0, SEEK_END); // read file
	int fileLen = ftell(pfile);
	rewind(pfile);
    char *buffer = new char[fileLen];

	fread((void *) &buffer[0], sizeof(char), fileLen, pfile);
	fclose(pfile);
    
    
//    ifstream pFile;
//    
//    if (!pFile.is_open()) //check is file has been opened
//    {
//        pFile.open (filePath, ios::in | ios::out);
//        
//        if (!pFile)
//        {
//            Rprintf("Cannot open file %s \n", filePath);
//            //cerr << "Failed to open " << pFile << endl;
//            exit(EXIT_FAILURE);  //abort program
//        }
//    }
//    
//    int fileLen;
//    char * buffer;
//    string lines;
//    
//    pFile.seekg(0, ios::end);
//    fileLen = pFile.tellg();
//    pFile.seekg (0, ios::beg);
//    
//    buffer = new char [fileLen];
//    
//    pFile.read (buffer, fileLen);
    
    //store the buffer into an matrix
    stringstream mat, ss;
    mat << buffer;
    char *line = new char[100];
    
    float **myMat = new float*[nRow];
    int cols, row;
    
    for (row = 0; mat.getline( line,  BUFFSIZE ) && row < nRow; row++ ) {
        // copy the entire buffered line into the stringstream
        ss << line;
        
        myMat[row] = new float[nCol];
        
        int tmpCol = 1;
        for (int col = 0; ss.getline( line, 20, delim ) && col < nCol; col++ ){
            myMat[row][col] = atof(line);
            tmpCol++;
        }
        
        if (row == 0)
            cols = tmpCol;
        if (tmpCol != cols){
            Rprintf( "Wrong Number of Filds in line %i \n", nRow );
            //cerr << "Wrong Number of Filds in line " << nRow << endl;
            exit(1);
        }
        
        ss << "";
        ss.clear();
    }
    
    nRow = row;
    nCol = cols;
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
    //cout << "There are exactly " << nRow << " rows and " << nCol << " columns in the matrix" << endl;
    return myMat;
}










/*
  double dlmread2(char *filePath, char delim)
 
  Read a text file into a matrix. First, matrix dimension is determined.
  Return the matrix, update dimension in nRow, nCol
  Stringstream is used for splitting the string.
*/
float **dlmread2(char *filePath, int &nRow, int &nCol, char delim, int flag){
    
    //read the whole file into a buffer
    const int BUFFSIZE = 1000;
    
    
    FILE *pfile = fopen(filePath, "r");
	if (!pfile) {
		Rprintf("Cannot open file %s \n", filePath);
		exit( EXIT_FAILURE);
	}
	fseek(pfile, 0, SEEK_END); // read file
	int fileLen = ftell(pfile);
	rewind(pfile);
    char *buffer = new char[fileLen];
    
	fread((void *) &buffer[0], sizeof(char), fileLen, pfile);
	fclose(pfile);
    
//    ifstream pFile;
//    
//    if (!pFile.is_open()) //check is file has been opened
//    {
//        pFile.open (filePath, ios::in | ios::out);
//        
//        if (!pFile)
//        {
//            Rprintf("Cannot open file %s \n", filePath);
//            //cerr << "Failed to open " << filePath << endl;
//            exit(EXIT_FAILURE);  //abort program
//        }
//    }
//    
//    int fileLen;
//    char * buffer;
//    string lines;
//    
//    pFile.seekg(0, ios::end);
//    fileLen = pFile.tellg();
//    pFile.seekg (0, ios::beg);
//    
//    buffer = new char [fileLen+10];
//    pFile.read (buffer, fileLen);
    
    
    //find the matrix dimension
    int ind = 0;
    nRow = 0; nCol = 1;
    
    for (int ind = 0; ind < fileLen; ind++) {
        int tmpCol = 1;
        for ( ; (ind < fileLen - 1) && (buffer[ind] != '\r') && (buffer[ind] != '\n'); ind++){
            if (buffer[ind] == delim)   tmpCol++;
        }
        
        if ((ind < fileLen - 1) && (buffer[ind] == '\r') || (buffer[ind] == '\n')){
			nRow++; //blank lines
            
            if (nRow == 1)  nCol = tmpCol;
            
            if (tmpCol != nCol){
				Rprintf( "Wrong Number of Filds in line %i \n", nRow );
                //cerr << "Wrong Number of Filds in line " << nRow << endl;
                exit(1);
            }
        }
    }
    
    
    //store the buffer into an matrix
    stringstream mat, ss;
    mat << buffer;
    SafeDelete(buffer);
//    delete[] buffer; buffer = NULL;
    
    char *line = new char[100];
    
    float **myMat = new float*[nRow];
    
    for (int row = 0; mat.getline( line,  BUFFSIZE ) && row < nRow; row++ ) {
        // copy the entire buffered line into the stringstream
        ss << line;
        
        myMat[row] = new float[nCol];
        
        for (int col = 0; ss.getline( line, 20, delim ) && col < nCol; col++ )
            myMat[row][col] = atof(line);
        
        ss << "";
        ss.clear();
    }
    
    //delete[] line; line = NULL;
    SafeDelete(line);
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
    //cout << "There are exactly " << nRow << " rows and " << nCol << " columns in the matrix" << endl;
    
    return myMat;
    
}




/*
  double dlmread2(char *filePath, char delim)
 
  Read a text file into a matrix. First, matrix dimension is determined.
  Return the matrix, update dimension in nRow, nCol
  Stringstream is used for splitting the string.
*/
int **readFile(char *filePath, int &nRow, char delim){
    
    //read the whole file into a buffer
    const int BUFFSIZE = 1000;
    
    FILE *pfile = fopen(filePath, "r");
	if (!pfile) {
		Rprintf("Cannot open file %s \n", filePath);
		exit( EXIT_FAILURE);
	}
	fseek(pfile, 0, SEEK_END); // read file
	int fileLen = ftell(pfile);
	rewind(pfile);
    char *buffer = new char[fileLen+10];
    
	fread((void *) &buffer[0], sizeof(char), fileLen, pfile);
	fclose(pfile);
	
//	cout << buffer << endl;

    
//    ifstream pFile;
//    
//    if (!pFile.is_open()) //check is file has been opened
//    {
//        pFile.open (filePath, ios::in | ios::out);
//        
//        if (!pFile)
//        {
//            Rprintf("Cannot open file %s \n", filePath);
//            //cerr << "Failed to open " << filePath << endl;
//            exit(EXIT_FAILURE);  //abort program
//        }
//    }
//    
//    int fileLen;
//    char * buffer;
//    string lines;
//    
//    pFile.seekg(0, ios::end);
//    fileLen = pFile.tellg();
//    pFile.seekg (0, ios::beg);
//    
//    buffer = new char [fileLen+10];
//    pFile.read (buffer, fileLen);
    
    //find the matrix dimension
    int ind = 0;
//    nRow = 0;
    
    stringstream mat, ss;
    mat << buffer;
    
    
    char *line = new char[100];
//    size_t len = 0;
//    
//    FILE *fp = fopen(filePath, "r");
//    if (fp == NULL)
//        exit(EXIT_FAILURE);
//    
//    while(!feof(fp)) {
//        fgets(line, sizeof line, fp);
//        ////getline(&line, &len, fp)
//        if(strcmp (line, "") != 0 || strcmp (line, "\n") != 0 || strcmp (line, "\r
//                                                                         ") != 0){
//            cout << "ll" << line << "ll" << endl;
//            nRow++;
//        }
//    }
//    cout << "nRow: " << nRow << endl;
//    
//    fclose(fp);
    
    
/*    for (nRow = 0; mat.getline( line,  BUFFSIZE); ){
	  if(strcmp (line, "") == 0 || strcmp (line, " ") == 0)
	  continue;
	  cout << "ll" << line << "ll" << endl;
	  nRow++;
	  }
*/
/*
  mat.getline( line,  BUFFSIZE);
  ss << line;
  cout << line << endl;
  ss.getline( line, 20, delim );
  cout << line << endl;
  nRow = atoi(line);
  ss << "";
  ss.clear();
  line = "";
*/
//int rows[5] = {6, 358, 3793, 14320, 78300}; 
//nRow = rows[level];
    //store the buffer into a two dimensional array
    
//    mat <<"";
//    mat.clear();
    
//    mat << buffer;
    //delete[] buffer; buffer = NULL;
    SafeDelete(buffer);
    
    int **myMat = new int*[nRow];
    int nCol = 100;
    
	cout << "************************\n";
    for (int row = 0; row < nRow; row++ ) {
//cout << row << endl;
        mat.getline( line,  BUFFSIZE );
//cout << line << endl;
//	if(row == 0)
//		cout << line << endl;
        // copy the entire buffered line into the stringstream

        ss << line;
        myMat[row] = new int[nCol];        
        for (int col = 0; ss.getline( line, 20, delim ); col++ ){
			//	if(row == 0)
			//		    cout << line << "\n";
            myMat[row][col] = atoi(line);
//	   cout << myMat[row][col] << endl;
		}        
        ss << "";
        ss.clear();
    }
    
    SafeDelete(line);
    //delete[] line; line = NULL;
    
	//  for (int i = 0 ; i < nRow; i++) {
	//     cout << myMat[i][0] << "\t" << myMat[i][1] << "\t" <<myMat[i][2] << endl;
	// }
    
    
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
    //cout << "There are exactly " << nRow << " rows and " << nCol << " columns in the matrix" << endl;

//	printMat(myMat, nRow, 2);

//    for (int i = 0; i < 3; i++)
//	cout << myMat[i][0] << "\t" << myMat[i][1] << endl;

    return myMat;
    
}



/*
  double dlmread3(char *filePath, char delim)
 
  Read a matrix stored in a text file and save it into one dimensional array.
  Matrix dimension is determined.
  A simple loop is used for splitting the string.
*/
/*
  void *dlmread3(char *filePath, char delim, int &nRow, int &nCol, int flag ){
    
  //Read the file into a buffer
  ifstream pFile;
  if (!pFile.is_open()) //check is file has been opened
  {
  pFile.open (filePath, ios::in | ios::out);
        
  if (!pFile)
  {
  cerr << "Failed to open '" << filePath << "' " << endl;
  exit(EXIT_FAILURE);  //abort program
  }
  }
    
  int fileLen;
  char * buffer;
  string lines;
    
  pFile.seekg(0, ios::end);
  fileLen = pFile.tellg();
  pFile.seekg (0, ios::beg);
    
  buffer = new char [fileLen];
    
  pFile.read (buffer, fileLen);
  pFile.close();
    
    
  //Save the file into 1 dimensional array
  int arrSize = floor(fileLen/3);
  float *myArr = new float[arrSize];
  char *cell;
  int oldInd = 0, ii = 0, ind = 0, tmpCol = 1;
  nCol = 1; nRow = 0;
    
  for ( ; (ind < fileLen - 1) && (buffer[ind] != '\n'); ind++) {
        
  oldInd = ind;
  int sizeCell = 0;
  for ( ; (buffer[ind] != '\n') && buffer[ind] != delim ; ind++) {
  sizeCell++;
  }
        
  cell = new char[sizeCell];
  memcpy(cell, buffer+oldInd, sizeCell*sizeof(char));
  if(buffer[ind] == delim){
  tmpCol++;
  }
  else if ((ind < fileLen - 1) && (buffer[ind] == '\r') || (buffer[ind] == '\n')){
  nRow++;
            
  if (nRow == 1)  nCol = tmpCol;
  if (tmpCol != nCol){
  cerr << "Wrong Number of Filds in line " << nRow << endl;
  exit(1);
  }
  tmpCol = 1;
  }
        
  //????????????????????????????????????????
  if (ii < arrSize){
  myArr[ii++] = atof(cell);
  }
  else{ //reallocate
  cerr << "size err" << endl;
  exit(1);
  }
  }
    
  delete[] buffer; buffer = NULL;
    
  cout << "There are exactly " << nRow << " rows and " << nCol << " columns in the matrix" << endl;
  if (flag == 1)
  return myArr
  ;
  else{
  cout << ii << endl;
  cout << nRow << "\t" << nCol << endl;
        
  float **myMat = new float* [nRow];
        
  for (int jj = 0; jj < nCol; jj++) {
  myMat[jj] = new float[nCol];
  myMat[jj] = myArr + jj*nCol;
  }
        
  for (int jj = 0; jj < nCol; jj++) {
  cout << myMat[jj][2] << "\t";
  }
  cout << endl;
        
  delete[] myArr;
  myArr = NULL;
        
  for (int jj = 0; jj < nCol; jj++) {
  cout << myMat[jj][2] << "\t";
  }
  cout << endl;
        
  return myMat;
  }
    
  }


  double dlmread3(char *filePath, char delim)
 
  Read a matrix stored in a text file and save it into a two dimensional array.
  Matrix dimension is determined.
  A simple loop is used for splitting the string.
*/

/*
  float **dlmread3(char *filePath, char delim, int &nRow, int &nCol){
    
  //Read the file into a buffer
  ifstream pFile;
  if (!pFile.is_open()) //check is file has been opened
  {
  pFile.open (filePath, ios::in | ios::out);
        
  if (!pFile)
  {
  cerr << "Failed to open '" << filePath << "' " << endl;
  exit(EXIT_FAILURE);  //abort program
  }
  }
    
  int fileLen;
  char * buffer;
  string lines;
    
  pFile.seekg(0, ios::end);
  fileLen = pFile.tellg();
  pFile.seekg (0, ios::beg);
    
  buffer = new char [fileLen];
    
  pFile.read (buffer, fileLen);
  pFile.close();
    
    
  //First save the file into 1 dimensional array (myArr) and find the matrix dimension
  float *myArr = new float[fileLen/2];
  char *cell;
  int oldInd = 0, ii = 0, ind = 0, tmpCol = 1;
  nCol = 1; nRow = 0; //matrix dimension
    
  for ( ; (ind < fileLen - 1) && (buffer[ind] != '\n'); ind++) {
        
  oldInd = ind;
  int sizeCell = 0;
  for ( ; (buffer[ind] != '\n') && buffer[ind] != delim ; ind++) {
  sizeCell++;
  }
        
  cell = new char[sizeCell];
  memcpy(cell, buffer+oldInd, sizeCell*sizeof(char));
  if(buffer[ind] == delim){
  tmpCol++;
  }
  else if ((ind < fileLen - 1) && (buffer[ind] == '\r') || (buffer[ind] == '\n')){
  nRow++;
            
  if (nRow == 1)  nCol = tmpCol;
  if (tmpCol != nCol){
  cerr << "Wrong Number of Filds in line " << nRow << endl;
  exit(1);
  }
  tmpCol = 1;
  }
        
  //????????????????????????????????????????
  if (ii < fileLen/2){
  myArr[ii++] = atof(cell);
  }
  else{
  cerr << "size err" << endl;
  exit(1);
  }
  }
    
  delete[] buffer; buffer = NULL;
    
  cout << ii << endl;
  cout << nRow << "\t" << nCol << endl;
    
  float **myMat = new float* [nRow];
    
  for (int jj = 0; jj < nCol; jj++) {
  myMat[jj] = new float[nCol];
  myMat[jj] = myArr + jj*nCol;
  }
    
  for (int jj = 0; jj < nCol; jj++) {
  cout << myMat[jj][2] << "\t";
  }
  cout << endl;
    
  delete[] myArr;
  myArr = NULL;
    
  for (int jj = 0; jj < nCol; jj++) {
  cout << myMat[jj][2] << "\t";
  }
  cout << endl;
    
    
  cout << "There are exactly " << nRow << " rows and " << nCol << " columns in the matrix" << endl;
  return myMat;
    
  }
*/
