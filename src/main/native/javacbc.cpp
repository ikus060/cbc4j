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

#include  "CoinTime.hpp"
#include "javacbc.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcHeuristicFPump.hpp"

// Included for solver
#include "CbcOrClpParam.hpp"



// Solve initial LP relaxation.
void initialSolve(OsiClpSolverInterface& solver){
  solver.initialSolve();
  // // Flush the stdout.
  std::cout.flush();
}

void setInt(OsiClpSolverInterface& solver, int i){
  solver.setInteger(i);
}

// Default Constructor.
OsiClpSolverInterface* newOsiClpSolverInterface(){
  return new OsiClpSolverInterface();
}

// Delete OsiClpSolverInterface
void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver) {
  delete solver;
}

// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiClpSolverInterface& s){
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
bool isProvenOptimal(CbcModel& c){
    return c.isProvenOptimal();
}
bool isProvenInfeasible(CbcModel& c){
    return c.isProvenInfeasible();
}
int getNumCols(CbcModel& c){
    return c.getNumCols();
}
int getNumRows(CbcModel& c){
    return c.getNumRows();
}
double getObjSense(CbcModel& c){
    return c.getObjSense();
}
double getColLower(CbcModel& c, int colIndex){
    return c.getColLower()[colIndex];
}
double getColUpper(CbcModel& c, int colIndex){
    return c.getColUpper()[colIndex];
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
const CoinPackedMatrix* getMatrixByCol(CbcModel& c){
    return c.getMatrixByCol();
}
double getInfinity(CbcModel& c){
    return c.getInfinity();
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
double getObjValue(CbcModel& c){
    return c.getObjValue();
}
void setLogLevel(CbcModel& c, int n){
    c.setLogLevel(n);
}

//OsiClpInterface proxy methods
//Methods returning info on how the solution process terminated
bool isProvenOptimal(OsiClpSolverInterface& osi){ //Is optimality proven?
    return osi.isProvenOptimal();
}

//Methods related to querying the input data
int getNumCols(OsiClpSolverInterface& osi){ // Get number of columns.
    return osi.getNumCols();
}
int getNumRows(OsiClpSolverInterface& osi){ // Get number of rows.
    return osi.getNumRows();
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
const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface& osi){ // Get pointer to column-wise copy of matrix.
    return osi.getMatrixByCol();
}
double getInfinity(OsiClpSolverInterface& osi){ // Get solver's value for infinity.
    return osi.getInfinity();
}

//Methods related to querying the solution
double getColSolution(OsiClpSolverInterface& osi, int colIndex){ // Get pointer to array[getNumCols()] of primal solution vector.
    return osi.getColSolution()[colIndex];
}
double getObjValue(OsiClpSolverInterface& osi){ // Get objective function value.
    return osi.getObjValue();
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
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single row lower bound, Use -DBL_MAX for -infinity.
    osi.setRowLower(elementIndex,elementValue);
}
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue){// Set a single row upper bound, Use DBL_MAX for infinity.
    osi.setRowUpper(elementIndex,elementValue);
}
void setObjective(OsiClpSolverInterface& osi, const double *array){// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
    osi.setObjective(array);
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

/*Methods to expand a problem.
  Note that if a column is added then by default
  it will correspond to a continuous variable.  */

void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj){// Add a column (primal variable) to the problem.
    osi.addCol(numberElements,rows,elements,collb,colub,obj);
}
void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices){// Remove a set of columns (primal variables) from the problem.
    osi.deleteCols(num,colIndices);
}
void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub){// Add a row (constraint) to the problem.
    osi.addRow(numberElements,columns,element,rowlb,rowub);
}
void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero){
    osi.modifyCoefficient(row,column,newElement,keepZero);
}
void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices){// Delete a set of rows (constraints) from the problem.
    osi.deleteRows(num,rowIndices);
}

//Methods to input a problem
int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon){// Read file in LP format (with names).
    osi.readLp(filename,epsilon);
}
void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension, double epsilon, int numberAcross, int decimals, double objSense, bool useRowNames){// Write the problem into an Lp file of the given filename.
    osi.writeLp(filename,extension,epsilon,numberAcross,decimals,objSense,useRowNames);
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

double bestSolution(CbcModel& c, int colIndex){
return c.bestSolution()[colIndex];
}

// Simple method to solve the model with default settings.
void solve(CbcModel &model) {

	int logLevel = 0;

	// Reduce printout
	model.setLogLevel(logLevel);
	model.messageHandler()->setLogLevel(logLevel);
    model.solver()->messageHandler()->setLogLevel(logLevel);
		model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);

	// TODO initialSolve on Clp

	model.initialSolve();

	// Could tune more
	double objValue = model.solver()->getObjSense() * model.solver()->getObjValue();
	double minimumDropA=CoinMin(1.0,fabs(objValue)*1.0e-3+1.0e-4);
	double minimumDrop= fabs(objValue)*1.0e-4+1.0e-4;
	model.setMinimumDrop(minimumDrop);



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
	CbcStrategyDefault strategy(true,0,0);

	// Do complete search
	model.branchAndBound();

}
