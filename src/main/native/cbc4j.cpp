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
#include "cbc4j.hpp"
#include "CoinPackedMatrix.hpp"
#include "CbcHeuristicFPump.hpp"

// Included for solver
#include "CbcOrClpParam.hpp"



void addCol(OsiClpSolverInterface& osi, int numberElements, const int* rows, const double* elements, const double collb, const double colub, const double obj){// Add a column (primal variable) to the problem.
    osi.addCol(numberElements,rows,elements,collb,colub,obj);
}

void addRow(OsiClpSolverInterface& osi, int numberElements, const int* columns, const double* element, const double rowlb, const double rowub){// Add a row (constraint) to the problem.
    osi.addRow(numberElements,columns,element,rowlb,rowub);
}

double bestSolution(CbcModel& c, int colIndex){
    return c.bestSolution()[colIndex];
}

void branchAndBound(CbcModel& c, int doStatistics){
    c.branchAndBound(doStatistics);
    std::cout.flush();
}

void callCbc0 (CbcModel &m) {
    CbcMain0(m);
    std::cout.flush();
}

int callCbc1 (int argc, const char* argv[], CbcModel &m) {
    int val = CbcMain1(argc, argv, m);
    std::cout.flush();
    return val;
}

void deleteCbcModel(CbcModel* c) {
    delete c;
}

// Remove a set of columns (primal variables) from the problem.
void deleteCols(OsiClpSolverInterface& osi, const int num, const int* colIndices){
    osi.deleteCols(num,colIndices);
}

void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver) {
    delete solver;
}

// Delete a set of rows (constraints) from the problem.
void deleteRows(OsiClpSolverInterface& osi, const int num, const int* rowIndices){
    osi.deleteRows(num,rowIndices);
}

double getCbcColSolution(CbcModel& c, int colIndex){
    return c.getCbcColSolution()[colIndex];
}

double getCbcRowLower(CbcModel& c, int rowIndex){
    return c.getCbcRowLower()[rowIndex];
}

double getCbcRowUpper(CbcModel& c, int rowIndex){
    return c.getCbcRowUpper()[rowIndex];
}

double getCoefficient(OsiClpSolverInterface &osi, int row, int col) {
    return osi.getMatrixByCol()->getCoefficient(row, col);
}

double getColLower(CbcModel& c, int colIndex){
    return c.getColLower()[colIndex];
}

// Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of column lower bounds.
double getColLower(OsiClpSolverInterface& osi, int colIndex){ 
    return osi.getColLower()[colIndex];
}

// Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen){ 
    return osi.getColName(colIndex,maxLen);
}

// Get pointer to array[getNumCols()] of primal solution vector.
double getColSolution(OsiClpSolverInterface& osi, int colIndex){ 
    return osi.getColSolution()[colIndex];
}

double getColSolution(OsiSolverInterface& osi, int colIndex){
    return osi.getColSolution()[colIndex];
}

double getColUpper(CbcModel& c, int colIndex){
    return c.getColUpper()[colIndex];
}

// Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of column upper bounds.
double getColUpper(OsiClpSolverInterface& osi, int colIndex){ 
    return osi.getColUpper()[colIndex];
}

double getInfinity(CbcModel& c){
    return c.getInfinity();
}

// Get solver's value for infinity.
double getInfinity(OsiClpSolverInterface& osi){ 
    return osi.getInfinity();
}

int getNumCols(CbcModel& c){
    return c.getNumCols();
}

// Get number of columns.
int getNumCols(OsiClpSolverInterface& osi){ 
    return osi.getNumCols();
}

int getNumRows(CbcModel& c){
    return c.getNumRows();
}

// Get number of rows.
int getNumRows(OsiClpSolverInterface& osi){
    return osi.getNumRows();
}


double getObjCoefficients(CbcModel& c, int colIndex){
    return c.getObjCoefficients()[colIndex];
}

// Get pointer to array[getNumCols(OsiClpSolverInterface& osi)] of objective function coefficients.
double getObjCoefficients(OsiClpSolverInterface& osi, int colIndex){ 
    return osi.getObjCoefficients()[colIndex];
}

double getObjSense(CbcModel& c){
    return c.getObjSense();
}

// Get objective function sense (1 for min (default), -1 for max).
double getObjSense(OsiClpSolverInterface& osi){ 
    return osi.getObjSense();
}

double getObjValue(CbcModel& c){
    return c.getObjValue();
}

// Get objective function value.
double getObjValue(OsiClpSolverInterface& osi){ 
    return osi.getObjValue();
}

double getObjValue(OsiSolverInterface& c) {
    return c.getObjValue();
}

double getRowLower(CbcModel& c, int rowIndex){
    return c.getRowLower()[rowIndex];
}

// Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row lower bounds.
double getRowLower(OsiClpSolverInterface& osi, int rowIndex){ 
    return osi.getRowLower()[rowIndex];
}

// Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen){ 
    return osi.getRowName(rowIndex,maxLen);
}

double getRowUpper(CbcModel& c, int rowIndex){
    return c.getRowUpper()[rowIndex];
}

// Get pointer to array[getNumRows(OsiClpSolverInterface& osi)] of row upper bounds.
double getRowUpper(OsiClpSolverInterface& osi, int rowIndex){ 
    return osi.getRowUpper()[rowIndex];
}

void initialSolve(CbcModel& c) {
    c.initialSolve();
}

void initialSolve(OsiClpSolverInterface& solver){
    solver.initialSolve();
    std::cout.flush();
}

bool isBinary(CbcModel& c, int colIndex){
    return c.isBinary(colIndex);
}

// Return true if variable is binary.
bool isBinary(OsiClpSolverInterface& osi, int colIndex){
    return osi.isBinary(colIndex);
}

bool isContinuous(CbcModel& c, int colIndex){
    return c.isContinuous(colIndex);
}

// Return true if column is continuous.
bool isContinuous(OsiClpSolverInterface& osi, int colNumber){ 
    return osi.isContinuous(colNumber);
}

bool isInteger(CbcModel& c, int colIndex){
    return c.isInteger(colIndex);
}

// Return true if column is integer.
bool isInteger(OsiClpSolverInterface& osi, int colIndex){ 
    return osi.isInteger(colIndex);
}

bool isProvenInfeasible(CbcModel& c){
    return c.isProvenInfeasible();
}

// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiClpSolverInterface& s){
    return new CbcModel(s);
}

// Default Constructor.
bool isProvenOptimal(CbcModel& c){
    return c.isProvenOptimal();
}

//Is optimality proven?
bool isProvenOptimal(OsiClpSolverInterface& osi){ 
    return osi.isProvenOptimal();
}

void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero){
    osi.modifyCoefficient(row,column,newElement,keepZero);
}

OsiClpSolverInterface* newOsiClpSolverInterface(){
    return new OsiClpSolverInterface();
}

// Read file in LP format (with names).
int readLp(OsiClpSolverInterface& osi, const char* filename, const double epsilon){
    osi.readLp(filename,epsilon);
}

//Integrality related changing methods
void resolve(OsiSolverInterface& osi) {
    osi.resolve();
}

// Returns solver - has current state.
OsiClpSolverInterface* solver(CbcModel& m) {
    OsiSolverInterface* solver = m.solver();
    return dynamic_cast<OsiClpSolverInterface*>(solver);
}

void setCoefficients(OsiClpSolverInterface &osi, int row, int numberElements, const int* cols, const double* coefs) {
    int colCount = osi.getNumCols();
    for(int col=0;col<colCount;col++) {
        osi.modifyCoefficient(row, col, 0.0, false);
    }
    for(int col=0;col<numberElements;col++) {
        osi.modifyCoefficient(row, cols[col], coefs[col], false);
    }
}

// Set a single column lower and upper bound.
void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper){
    osi.setColBounds(elementIndex,lower,upper);
}

// Set a single column lower bound, Use -DBL_MAX for -infinity.
void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue){
    osi.setColLower(elementIndex,elementValue);
}

// Set name of column.
void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name){
    osi.setColName(colIndex,name);
}

// Set a single column upper bound, Use DBL_MAX for infinity.
void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue){
    osi.setColUpper(elementIndex,elementValue);
}

// Set the index-th variable to be a continuous variable.
void setContinuous(OsiClpSolverInterface& osi, int index){
    osi.setContinuous(index);
}

// Set the index-th variable to be an integer variable.
void setInteger(OsiClpSolverInterface& osi, int index){
    osi.setInteger(index);
}
void setLogLevel(CbcModel& c, int n){
    c.setLogLevel(n);
}

// Resolve an LP relaxation after problem modification.
void setLogLevel(OsiClpSolverInterface& s, int n) {
    s.messageHandler()->setLogLevel(n);
}

// Get pointer to array[getNumCols()] of objective function coefficients.
const double* getObjCoefficients(OsiClpSolverInterface& osi) {
    return osi.getObjCoefficients();
}

// Set an objective function setObjCoeff.
void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue){
    osi.setObjCoeff(elementIndex,elementValue);
}

// Sets the objective coefficient
void setObjCoefficients(OsiClpSolverInterface &osi, int numberElements, const int* cols, const double* coefs) {
    int colCount = osi.getNumCols();
    for(int col=0;col<colCount;col++) {
        osi.setObjCoeff(col, 0);
    }
    for(int col=0;col<numberElements;col++) {
        osi.setObjCoeff(cols[col], coefs[col]);
    }
}

// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
void setObjective(OsiClpSolverInterface& osi, const double* array){
    osi.setObjective(array);
}

// Return the coefficient.
void setObjSense(CbcModel& c, double s){
    c.setObjSense(s);
}

// Run the solver using CbcMain1
void setObjSense(OsiClpSolverInterface& osi, double s){
    osi.setObjSense(s);
}

// Set a single row lower bound, Use -DBL_MAX for -infinity.
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue){
    osi.setRowLower(elementIndex,elementValue);
}

// Set name of row.
void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name){
    osi.setRowName(rowIndex,name);
}

// Set a single row upper bound, Use DBL_MAX for infinity.
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue){
    osi.setRowUpper(elementIndex,elementValue);
}

int status(CbcModel& c) {
    return c.status();
}

void writeLp(OsiClpSolverInterface& osi, const char* filename, const char* extension, double epsilon, int numberAcross, int decimals, double objSense, bool useRowNames){// Write the problem into an Lp file of the given filename.
    osi.writeLp(filename,extension,epsilon,numberAcross,decimals,objSense,useRowNames);
}
