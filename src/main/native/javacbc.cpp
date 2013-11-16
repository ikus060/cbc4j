// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>
#include <string>
using namespace std;



// For Branch and bound
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "OsiClpSolverInterface.hpp"
// Preprocessing
#include "CglPreProcess.hpp"

#include  "CoinTime.hpp"
#include "javacbc.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcHeuristicFPump.hpp"

// Included for solver
#include "CbcOrClpParam.hpp"

CoinPackedMatrix CPM(int numberRows, int numberColumns, int numberElements, int *startint, int nbelems, int *length, int *rows, double *elements){
    if(nbelems>0){
        CoinBigIndex *start;
        start = new CoinBigIndex[nbelems];
        int i=0;
        for(i=0;i<nbelems;i++){
            start[i] = startint[i];
        }
        CoinPackedMatrix matrix(true,numberRows,numberColumns,numberElements,elements,rows,start,length);
        return matrix;
    }
    CoinPackedMatrix matrix;
    return matrix;
}

CoinPackedMatrix CPM(int *rowIndices, int *colIndices, double *elements, int nbelems){
    if(nbelems>0){
        CoinPackedMatrix matrix(true,rowIndices,colIndices,elements,nbelems);
        return matrix;
    }
    CoinPackedMatrix matrix;
    return matrix;
}

// Solve initial LP relaxation.
void initialSolve(OsiClpSolverInterface& solver){
  solver.initialSolve();
  // // Flush the stdout.
  std::cout.flush();
}

// Solve initial LP relaxation.
void initialSolve(OsiSolverInterface& solver){
  solver.initialSolve();
  // Flush the stdout.
  std::cout.flush();
}

void setInt(OsiClpSolverInterface& solver, int i){
  solver.setInteger(i);
}

// Default Constructor. 
OsiClpSolverInterface* newOsiClpSolverInterface(){
  return new OsiClpSolverInterface();
}

// Copy constructor.
OsiClpSolverInterface* newOsiClpSolverInterface(OsiClpSolverInterface & s){
  return new OsiClpSolverInterface(s);
}

// Delete OsiClpSolverInterface
void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver) {
  delete solver;
}

CoinBigIndex CBI(int start[]){
  int length = sizeof(start)/sizeof(int);
  CoinBigIndex cbi[length];
  for(int i=0;i<length;i++){
    cbi[i] = start[i];
  }
  return *cbi;
}

// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiClpSolverInterface& s){
    return new CbcModel(s);
}
// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiSolverInterface& s){
    return new CbcModel(s);
}
void deleteCbcModel(CbcModel* c) {
    delete c;
}
void initialSolve(CbcModel& c) {
    c.initialSolve();
}
void branchAndBound(CbcModel& c, int doStatistics){
    c.branchAndBound(doStatistics);
    std::cout.flush();
}
int status(CbcModel& c) {
    return c.status();
}
int secondaryStatus(CbcModel& c) {
    return c.secondaryStatus();
}
bool isProvenOptimal(CbcModel& c){
    return c.isProvenOptimal();
}
bool isProvenInfeasible(CbcModel& c){
    return c.isProvenInfeasible();
}
bool isContinuousUnbounded(CbcModel& c){
    return c.isContinuousUnbounded();
}
bool isNodeLimitReached(CbcModel& c){
    return c.isNodeLimitReached();
}
bool isSecondsLimitReached(CbcModel& c){
    return c.isSecondsLimitReached();
}
int getNumCols(CbcModel& c){
    return c.getNumCols();
}
int getNumRows(CbcModel& c){
    return c.getNumRows();
}
CoinBigIndex getNumElements (CbcModel& c){
    return c.getNumElements();
}
int numberIntegers(CbcModel& c){
    return c.numberIntegers();
}
char integerType(CbcModel& c, int i){
    return c.integerType(i);
}
double getColLower(CbcModel& c, int colIndex){
    return c.getColLower()[colIndex];
}
double getColUpper(CbcModel& c, int colIndex){
    return c.getColUpper()[colIndex];
}
const char* getRowSense(CbcModel& c){
    return c.getRowSense();
}
double getRightHandSide(CbcModel& c, int rowIndex){
    return c.getRightHandSide()[rowIndex];
}
double getRowRange(CbcModel& c, int rowIndex){
    return c.getRowRange()[rowIndex];
}
double getRowLower(CbcModel& c, int rowIndex){
    return c.getRowLower()[rowIndex];
}
double getRowUpper(CbcModel& c, int rowIndex){
    return c.getRowUpper()[rowIndex];
}
double getObjCoefficients(CbcModel& c, int colIndex){
    return c.getObjCoefficients()[colIndex];
}
double getObjSense(CbcModel& c){
    return c.getObjSense();
}
void setObjSense(CbcModel& c, double s){
    c.setObjSense(s);
}
bool isContinuous(CbcModel& c, int colIndex){
return c.isContinuous(colIndex);
}
bool isBinary(CbcModel& c, int colIndex){
return c.isBinary(colIndex);
}
bool isInteger(CbcModel& c, int colIndex){
    return c.isInteger(colIndex);
}
bool isIntegerNonBinary(CbcModel& c, int colIndex){
    return c.isIntegerNonBinary(colIndex);
}
bool isFreeBinary(CbcModel& c, int colIndex){
    return c.isFreeBinary(colIndex);
}
const CoinPackedMatrix* getMatrixByRow(CbcModel& c){
    return c.getMatrixByRow();
}
const CoinPackedMatrix* getMatrixByCol(CbcModel& c){
    return c.getMatrixByCol();
}
double getInfinity(CbcModel& c){
    return c.getInfinity();
}
double getCbcColLower(CbcModel& c, int colIndex){
    return c.getCbcColLower()[colIndex];
}
double getCbcColUpper(CbcModel& c, int colIndex){
    return c.getCbcColUpper()[colIndex];
}
double getCbcRowLower(CbcModel& c, int rowIndex){
    return c.getCbcRowLower()[rowIndex];
}
double getCbcRowUpper(CbcModel& c, int rowIndex){
    return c.getCbcRowUpper()[rowIndex];
}
double getCbcColSolution(CbcModel& c, int colIndex){
    return c.getCbcColSolution()[colIndex];
}
double getCbcRowPrice(CbcModel& c, int rowIndex){
    return c.getCbcRowPrice()[rowIndex];
}
double getObjValue(CbcModel& c){
    return c.getObjValue();
}
double bestSolution(CbcModel& c, int colIndex){
    return c.bestSolution()[colIndex];
}
int getNumberThreads(CbcModel& c){
    return c.getNumberThreads();
}
void setNumberThreads(CbcModel& c, int value){
    c.setNumberThreads(value);
}
void setLogLevel(CbcModel& c, int n){
    c.setLogLevel(n);
}
double getCurrentObjValue(CbcModel& c){
    return c.getCurrentObjValue();
}




//OsiClpInterface proxy methods
//Methods returning info on how the solution process terminated
bool isAbandoned(OsiClpSolverInterface& osi){ //Are there a numerical difficulties?
    return osi.isAbandoned();
}
bool isProvenOptimal(OsiClpSolverInterface& osi){ //Is optimality proven?
    return osi.isProvenOptimal();
}
bool isProvenPrimalInfeasible(OsiClpSolverInterface& osi){ // Is primal infeasiblity proven?
    return osi.isProvenPrimalInfeasible();
}
bool isProvenDualInfeasible(OsiClpSolverInterface& osi){ // Is dual infeasiblity proven?
    return osi.isProvenDualInfeasible();
 }
bool isPrimalObjectiveLimitReached(OsiClpSolverInterface& osi){ // Is the given primal objective limit reached?
    return osi.isPrimalObjectiveLimitReached();
 }
bool isDualObjectiveLimitReached(OsiClpSolverInterface& osi){ // Is the given dual objective limit reached?
    return osi.isDualObjectiveLimitReached();
}
bool isIterationLimitReached(OsiClpSolverInterface& osi){ // Iteration limit reached?
    return osi.isIterationLimitReached();
}

//Methods related to querying the input data
int getNumCols(OsiClpSolverInterface& osi){ // Get number of columns.
    return osi.getNumCols();
}
int getNumRows(OsiClpSolverInterface& osi){ // Get number of rows.
    return osi.getNumRows();
}
int getNumElements(OsiClpSolverInterface& osi){ // Get number of nonzero elements.
    return osi.getNumElements();
}
std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen){ // Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
    return osi.getRowName(rowIndex,maxLen);
}
std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen){ // Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
    return osi.getColName(colIndex,maxLen);
}
double getColLower(OsiClpSolverInterface& osi, int colIndex){ // Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of column lower bounds.
    return osi.getColLower()[colIndex];
}
double getColUpper(OsiClpSolverInterface& osi, int colIndex){ // Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of column upper bounds.
    return osi.getColUpper()[colIndex];
}
double getObjCoefficients(OsiClpSolverInterface& osi, int colIndex){ // Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of objective function coefficients.
    return osi.getObjCoefficients()[colIndex];
}
const char * getRowSense(OsiClpSolverInterface& osi){ // Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row constraint senses.
    return osi.getRowSense();
}
double getRightHandSide(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of rows right-hand sides.
    return osi.getRightHandSide()[rowIndex];
}
double getRowRange(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row ranges.
    return osi.getRowRange()[rowIndex];
}
double getRowLower(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row lower bounds.
    return osi.getRowLower()[rowIndex];
}
double getRowUpper(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row upper bounds.
    return osi.getRowUpper()[rowIndex];
}
double getObjSense(OsiClpSolverInterface& osi){ // Get objective function sense (1 for min (default), -1 for max).
    return osi.getObjSense();
}
void setObjSense(OsiClpSolverInterface& osi, double s){
    osi.setObjSense(s);
}
bool isContinuous(OsiClpSolverInterface& osi, int colNumber){ // Return true if column is continuous.
    return osi.isContinuous(colNumber);
}
bool isBinary(OsiClpSolverInterface& osi, int colIndex){ // Return true if variable is binary.
    return osi.isBinary(colIndex);
}
bool isInteger(OsiClpSolverInterface& osi, int colIndex){ // Return true if column is integer.
    return osi.isInteger(colIndex);
}
bool isIntegerNonBinary(OsiClpSolverInterface& osi, int colIndex){ // Return true if variable is general integer.
    return osi.isIntegerNonBinary(colIndex);
}
bool isFreeBinary(OsiClpSolverInterface& osi, int colIndex){ // Return true if variable is binary and not fixed at either bound.
    return osi.isFreeBinary(colIndex);
}
const char * getColType(OsiClpSolverInterface& osi, bool refresh){ // Return array of column length 0 - continuous 1 - binary (may get fixed later) 2 - general integer (may get fixed later).
    return osi.getColType(refresh);
}
const CoinPackedMatrix * getMatrixByRow(OsiClpSolverInterface& osi){ // Get pointer to row-wise copy of matrix.
    return osi.getMatrixByRow();
}
const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface& osi){ // Get pointer to column-wise copy of matrix.
    return osi.getMatrixByCol();
}
CoinPackedMatrix *getMutableMatrixByCol(OsiClpSolverInterface& osi){ // Get pointer to mutable column-wise copy of matrix.
    return osi.getMutableMatrixByCol();
}
double getInfinity(OsiClpSolverInterface& osi){ // Get solver's value for infinity.
    return osi.getInfinity();
}

//Methods related to querying the solution
double getColSolution(OsiClpSolverInterface& osi, int colIndex){ // Get pointer to array[getNumCols()] of primal solution vector.
    return osi.getColSolution()[colIndex];
}
double getRowPrice(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows()] of dual prices.
    return osi.getRowPrice()[rowIndex];
}
double getReducedCost(OsiClpSolverInterface& osi, int colIndex){ // Get a pointer to array[getNumCols()] of reduced costs.
    return osi.getReducedCost()[colIndex];
}
double getRowActivity(OsiClpSolverInterface& osi, int rowIndex){ // Get pointer to array[getNumRows()] of row activity levels (constraint matrix times the solution vector.
    return osi.getRowActivity()[rowIndex];
}
double getObjValue(OsiClpSolverInterface& osi){ // Get objective function value.
    return osi.getObjValue();
}
int getIterationCount(OsiClpSolverInterface& osi){ // Get how many iterations it took to solve the problem (whatever "iteration" mean to the solver.
    return osi.getIterationCount();
}
std::vector< double * > getDualRays(OsiClpSolverInterface& osi, int maxNumRays){ // Get as many dual rays as the solver can provide.
    return osi.getDualRays(maxNumRays);
}
std::vector< double * > getPrimalRays(OsiClpSolverInterface& osi, int maxNumRays){ // Get as many primal rays as the solver can provide.
    return osi.getPrimalRays(maxNumRays);
}



//Changing bounds on variables and constraints
void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set an objective function setObjCoeff.
    osi.setObjCoeff(elementIndex,elementValue);
}
void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single column lower bound, Use -DBL_MAX for -infinity.
    osi.setColLower(elementIndex,elementValue);
}
void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single column upper bound, Use DBL_MAX for infinity.
    osi.setColUpper(elementIndex,elementValue);
}
void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper){// Set a single column lower and upper bound.
    osi.setColBounds(elementIndex,lower,upper);
}
void setColSetBounds(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const double *boundList){// Set the bounds on a number of columns simultaneously, The default implementation just invokes setColLower() and setColUpper() over and over again.
    osi.setColSetBounds(indexFirst,indexLast,boundList);
}
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single row lower bound, Use -DBL_MAX for -infinity.
    osi.setRowLower(elementIndex,elementValue);
}
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single row upper bound, Use DBL_MAX for infinity.
    osi.setRowUpper(elementIndex,elementValue);
}
void setRowBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper){// Set a single row lower and upper bound.
    osi.setRowBounds(elementIndex,lower,upper);
}
void setRowType(OsiClpSolverInterface& osi, int index, char sense, double rightHandSide, double range){// Set the type of a single row.
    osi.setRowType(index,sense,rightHandSide,range);
}
void setRowSetBounds(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const double *boundList){// Set the bounds on a number of rows simultaneously. The default implementation just invokes setRowLower() and setRowUpper() over and over again.
    osi.setRowSetBounds(indexFirst,indexLast,boundList);
}
void setRowSetTypes(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const char *senseList, const double *rhsList, const double *rangeList){// Set the type of a number of rows simultaneously. The default implementation just invokes setRowType() over and over again.
    osi.setRowSetTypes(indexFirst,indexLast,senseList,rhsList,rangeList);
}
void setObjective(OsiClpSolverInterface& osi, const double *array){// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
    osi.setObjective(array);
}
void setColLower(OsiClpSolverInterface& osi, const double *array){// Set the lower bounds for all columns array [getNumCols()] is an array of values for the objective.
    osi.setColLower(array);
}
void setColUpper(OsiClpSolverInterface& osi, const double *array){// Set the upper bounds for all columns array [getNumCols()] is an array of values for the objective.
    osi.setColUpper(array);
}
void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name){// Set name of row.
    osi.setRowName(rowIndex,name);
}
void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name){// Set name of column.
    osi.setColName(colIndex,name);
}


//Integrality related changing methods
void setContinuous(OsiClpSolverInterface& osi, int index){// Set the index-th variable to be a continuous variable.
    osi.setContinuous(index);
}
void setInteger(OsiClpSolverInterface& osi, int index){// Set the index-th variable to be an integer variable.
    osi.setInteger(index);
}
void setContinuous(OsiClpSolverInterface& osi, const int *indices, int len){// Set the variables listed in indices (which is of length len) to be continuous variables.
    osi.setContinuous(indices,len);
}
void setInteger(OsiClpSolverInterface& osi, const int *indices, int len){// Set the variables listed in indices (which is of length len) to be integer variables.
    osi.setInteger(indices,len);
}



/*Methods to expand a problem.
  Note that if a column is added then by default
  it will correspond to a continuous variable.  */

void addCol(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const double collb, const double colub, const double obj){// Add a column (primal variable) to the problem.
    osi.addCol(vec,collb,colub,obj);
}
void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj){// Add a column (primal variable) to the problem.
    osi.addCol(numberElements,rows,elements,collb,colub,obj);
}
void addCols(OsiClpSolverInterface& osi, const int numcols, const CoinPackedVectorBase *const *cols, const double *collb, const double *colub, const double *obj){// Add a set of columns (primal variables) to the problem.
    osi.addCols (numcols,cols,collb,colub,obj);
}
void addCols(OsiClpSolverInterface& osi, const int numcols, const int *columnStarts, const int *rows, const double *elements, const double *collb, const double *colub, const double *obj){// Add a set of columns (primal variables) to the problem.
    osi.addCols (numcols,columnStarts,rows,elements,collb,colub,obj);
}
void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices){// Remove a set of columns (primal variables) from the problem.
    osi.deleteCols(num,colIndices);
}
void addRow(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const double rowlb, const double rowub){// Add a row (constraint) to the problem.
    osi.addRow(vec,rowlb,rowub);
}
void addRow(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const char rowsen, const double rowrhs, const double rowrng){// Add a row (constraint) to the problem.
    osi.addRow(vec,rowsen,rowrhs,rowrng);
}
void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub){// Add a row (constraint) to the problem.
    osi.addRow(numberElements,columns,element,rowlb,rowub);
}
void addRows(OsiClpSolverInterface& osi, const int numrows, const CoinPackedVectorBase *const *rows, const double *rowlb, const double *rowub){
    osi.addRows(numrows,rows,rowlb,rowub);
}
void addRows(OsiClpSolverInterface& osi, const int numrows, const CoinPackedVectorBase *const *rows, const char *rowsen, const double *rowrhs, const double *rowrng){// Add a set of rows (constraints) to the problem.
    osi.addRows(numrows,rows,rowsen,rowrhs,rowrng);
}
void addRows(OsiClpSolverInterface& osi, const int numrows, const int *rowStarts, const int *columns, const double *element, const double *rowlb, const double *rowub){// Add a set of rows (constraints) to the problem.
    osi.addRows(numrows,rowStarts,columns,element,rowlb,rowub);
}
void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero){
    osi.modifyCoefficient(row,column,newElement,keepZero);
}
void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices){// Delete a set of rows (constraints) from the problem.
    osi.deleteRows(num,rowIndices);
}
void saveBaseModel(OsiClpSolverInterface& osi){// If solver wants it can save a copy of "base" (continuous) model here.
    osi.saveBaseModel();
}
void restoreBaseModel(OsiClpSolverInterface& osi, int numberRows){// Strip off rows to get to this number of rows.
    osi.restoreBaseModel(numberRows);
}
void applyRowCuts(OsiClpSolverInterface& osi, int numberCuts, const OsiRowCut *cuts){// Apply a collection of row cuts which are all effective.
    osi.applyRowCuts(numberCuts,cuts);
}
void applyRowCuts(OsiClpSolverInterface& osi, int numberCuts, const OsiRowCut **cuts){// Apply a collection of row cuts which are all effective.
    osi.applyRowCuts(numberCuts,cuts);
}



//Methods to input a problem
void loadProblem(OsiClpSolverInterface& osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub){// Load in an problem by copying the arguments (the constraints on the rows are given by lower and upper bounds).
    osi.loadProblem(matrix,collb,colub,obj,rowlb,rowub);
}
void assignProblem(OsiClpSolverInterface& osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, double *&rowlb, double *&rowub){// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by lower and upper bounds).
    osi.assignProblem(matrix,collb,colub,obj,rowlb,rowub);
}
void loadProblem(OsiClpSolverInterface& osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng){// Load in an problem by copying the arguments (the constraints on the rows are given by sense/rhs/range triplets).
    osi.loadProblem(matrix,collb,colub,obj,rowsen,rowrhs,rowrng);
}
void assignProblem(OsiClpSolverInterface& osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, char *&rowsen, double *&rowrhs, double *&rowrng){// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by sense/rhs/range triplets).
    osi.assignProblem(matrix,collb,colub,obj,rowsen,rowrhs,rowrng);
}
void loadProblem(OsiClpSolverInterface& osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub){// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
    osi.loadProblem(numcols,numrows,start,index,value,collb,colub,obj,rowlb,rowub);
}
void loadProblem(OsiClpSolverInterface& osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng){// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
    osi.loadProblem(numcols,numrows,start,index,value,collb,colub,obj,rowsen,rowrhs,rowrng);
}
int loadFromCoinModel(OsiClpSolverInterface& osi, CoinModel &modelObject, bool keepSolution){// This loads a model from a coinModel object - returns number of errors.
    osi.loadFromCoinModel(modelObject,keepSolution);
}
int readMps(OsiClpSolverInterface& osi, const char *filename, const char *extension){// Read an mps file from the given filename (defaults to Osi reader) - returns number of errors (see OsiMpsReader class).
    osi.readMps(filename,extension);
}
int readMps(OsiClpSolverInterface& osi, const char *filename, bool keepNames, bool allowErrors){// Read an mps file from the given filename returns number of errors (see OsiMpsReader class).
    osi.readMps(filename,keepNames,allowErrors);
}
void writeMps(OsiClpSolverInterface& osi, const char *filename, const char *extension, double objSense){// Write the problem into an mps file of the given filename.
    osi.writeMps(filename,extension,objSense);
}
int writeMpsNative(OsiClpSolverInterface& osi, const char *filename, const char **rowNames, const char **columnNames, int formatType, int numberAcross, double objSense){// Write the problem into an mps file of the given filename, names may be null.
    osi.writeMpsNative(filename,rowNames,columnNames,formatType,numberAcross,objSense);
}
int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon){// Read file in LP format (with names).
    osi.readLp(filename,epsilon);
}
void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension, double epsilon, int numberAcross, int decimals, double objSense, bool useRowNames){// Write the problem into an Lp file of the given filename.
    osi.writeLp(filename,extension,epsilon,numberAcross,decimals,objSense,useRowNames);
}
void writeLp(OsiClpSolverInterface& osi, FILE *fp, double epsilon, int numberAcross, int decimals, double objSense, bool useRowNames){// Write the problem into the file pointed to by the parameter fp.
    osi.writeLp(fp,epsilon,numberAcross,decimals,objSense,useRowNames);
}

// Create a new pre-processor
CglPreProcess* newCglPreProcess() {
    return new CglPreProcess();
}

// Delete the pre-processor
void deleteCglPreProcess(CglPreProcess* p) {
    delete p;
}

// Run the preprocessor
OsiSolverInterface* preProcess(CglPreProcess& p, OsiClpSolverInterface& osi, bool makeEquality, int numberPasses) {
    OsiSolverInterface* o = p.preProcess(osi,makeEquality, numberPasses);
    // Flush the stdout.
    std::cout.flush();
    return o;
}

// Create solution in original model
void postProcess(CglPreProcess& p, OsiSolverInterface& osi) {
    p.postProcess(osi);
    // Flush the stdout.
    std::cout.flush();
}

// Resolve an LP relaxation after problem modification.
void resolve(OsiSolverInterface& osi) {
    osi.resolve();
}

// Returns solver - has current state. 
OsiClpSolverInterface* solver(CbcModel& m) {
    OsiSolverInterface *solver = m.solver();
    return dynamic_cast<OsiClpSolverInterface *>(solver);
}

// Create a new CbcHeuristicFPump
CbcHeuristicFPump* newCbcHeuristicFPump() {
    return new CbcHeuristicFPump();
}
// Create a new CbcHeuristicFPump
CbcHeuristicFPump* newCbcHeuristicFPump(CbcModel& m, double downValue, bool roundExpensive) {
    return new CbcHeuristicFPump(m, downValue, roundExpensive);
}
// Delete CbcHeuristicFPump
void deleteCbcHeuristicFPump(CbcHeuristicFPump * h) {
    delete h;
}
// Add one heuristic - up to user to delete. 
void addHeuristic(CbcModel& m, CbcHeuristicFPump* h) {
    m.addHeuristic(h);
}
// Delete OsiSolverInterface
void deleteOsiSolverInterface(OsiSolverInterface* s) {
    delete s;
}
// Set current log (detail) level. 
void setLogLevel(CglPreProcess& p, int n) {
    p.messageHandler()->setLogLevel(n);
}
// Set current log (detail) level.
void setLogLevel(OsiSolverInterface& s, int n) {
    s.messageHandler()->setLogLevel(n);
}
// Set current log (detail) level.
void setLogLevel(OsiClpSolverInterface& s, int n) {
    s.messageHandler()->setLogLevel(n);
}

// Get pointer to array[getNumCols()] of objective function coefficients. 
const double* getObjCoefficients(OsiClpSolverInterface& osi) {
    return osi.getObjCoefficients();
}

// Sets the objective coefficient
void setObjCoefficients(OsiClpSolverInterface &osi, int numberElements, const int *cols, const double *coefs) {
    int colCount = osi.getNumCols();
    for(int col=0;col<colCount;col++) {
        osi.setObjCoeff(col, 0);
    }
    for(int col=0;col<numberElements;col++) {
        osi.setObjCoeff(cols[col], coefs[col]);
    }
}

// Sets the coefficients for a row
void setCoefficients(OsiClpSolverInterface &osi, int row, int numberElements, const int *cols, const double *coefs) {
    int colCount = osi.getNumCols();
    for(int col=0;col<colCount;col++) {
        osi.modifyCoefficient(row, col, 0.0, false);
    }
    for(int col=0;col<numberElements;col++) {
        osi.modifyCoefficient(row, cols[col], coefs[col], false);
    }
}

// Return the coefficient.
double getCoefficient(OsiClpSolverInterface &osi, int row, int col) {
    return osi.getMatrixByCol()->getCoefficient(row, col);
}

// Run the solver using CbcMain1
void callCbc0 (CbcModel &m) {
    CbcMain0(m);
    std::cout.flush();
}

// Run solver using CbcMain1
int callCbc1 (int argc, const char *argv[], CbcModel &m) {
    int val = CbcMain1(argc, argv, m);
    std::cout.flush();
    return val;
}

// Get objective function value.
double getObjValue(OsiSolverInterface& c) {
    return c.getObjValue();
}

// Get pointer to array[getNumCols()] of primal solution vector.
double getColSolution(OsiSolverInterface& osi, int colIndex){ 
    return osi.getColSolution()[colIndex];
}

// Simple method to solve the model with default settings.
void solve(CbcModel &model) {

	int logLevel = 0;
	bool preProcess = true;

	// Reduce printout
	model.setLogLevel(logLevel);
	model.messageHandler()->setLogLevel(logLevel);
    model.solver()->messageHandler()->setLogLevel(logLevel);
	if (logLevel <= 1) {
		model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
	} else {
		model.solver()->setHintParam(OsiDoReducePrint, false, OsiHintTry);
	}
	
	// TODO initialSolve on Clp
	
	model.initialSolve();
	
	// Could tune more
	double objValue = model.solver()->getObjSense() * model.solver()->getObjValue();
	double minimumDropA=CoinMin(1.0,fabs(objValue)*1.0e-3+1.0e-4);
	double minimumDrop= fabs(objValue)*1.0e-4+1.0e-4;
	model.setMinimumDrop(minimumDrop);

	if (model.getNumCols()<500)
		model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
	else if (model.getNumCols()<5000)
		model.setMaximumCutPassesAtRoot(100); // use minimum drop
	else
		model.setMaximumCutPassesAtRoot(20);
		model.setMaximumCutPasses(10);

	// Switch off strong branching if wanted
	// model.setNumberStrong(0);
	// Do more strong branching if small
	if (model.getNumCols()<5000)
		model.setNumberStrong(10);
	else
		model.setNumberStrong(20);
		//model.setNumberStrong(5);
	model.setNumberBeforeTrust(5);

	model.solver()->setIntParam(OsiMaxNumIterationHotStart,100);
	
	// Default strategy will leave cut generators as they exist already
	// so cutsOnlyAtRoot (1) ignored
	// numberStrong (2) is 5 (default)
	// numberBeforeTrust (3) is 5 (default is 0)
	// printLevel (4) defaults (0)
	CbcStrategyDefault strategy(true,5,5);
	// Set up pre-processing to find sos if wanted
	if (preProcess) {
		strategy.setupPreProcessing(2);
		model.setStrategy(strategy);
	}
	
	// Do complete search
	model.branchAndBound();
	
}
