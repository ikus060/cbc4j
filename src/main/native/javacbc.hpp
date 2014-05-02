#ifndef JAVACBC_H
#define JAVACBC_H

#include <cassert>
#include <iomanip>
#include <string>
using namespace std;
// For Branch and bound
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "OsiClpSolverInterface.hpp"

#include  "CoinTime.hpp"
#include "CbcHeuristicFPump.hpp"

// Solve initial LP relaxation.
void initialSolve(OsiClpSolverInterface& solver);

void setInt(OsiClpSolverInterface& solver, int i);

// Default Constructor.
OsiClpSolverInterface* newOsiClpSolverInterface();

// Delete OsiClpSolverInterface
void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);

// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiClpSolverInterface& s);
void deleteCbcModel(CbcModel* c);

void initialSolve(CbcModel& c);
void branchAndBound(CbcModel& c, int doStatistics = 0);
int status(CbcModel& c);
bool isProvenOptimal(CbcModel& c);
bool isProvenInfeasible(CbcModel& c);
int getNumCols(CbcModel& c);
int getNumRows(CbcModel& c) ;
double getColLower(CbcModel& c, int colIndex);
double getColUpper(CbcModel& c, int colIndex);
double getRowLower(CbcModel& c, int rowIndex);
double getRowUpper(CbcModel& c, int rowIndex);
double getObjCoefficients(CbcModel& c, int colIndex);
double getObjSense(CbcModel& c);
void setObjSense(CbcModel& c, double s);
bool isContinuous(CbcModel& c, int colIndex);
bool isBinary(CbcModel& c, int colIndex);
bool isInteger(CbcModel& c, int colIndex);
const CoinPackedMatrix* getMatrixByCol(CbcModel& c);
double getInfinity(CbcModel& c);
double getCbcRowLower(CbcModel& c, int rowIndex);
double getCbcRowUpper(CbcModel& c, int rowIndex);
double getCbcColSolution(CbcModel& c, int colIndex);
double getObjValue(CbcModel& c);
void setLogLevel(CbcModel& c, int n);

//OsiClpInterface proxy methods
//Methods returning info on how the solution process terminated
bool isProvenOptimal(OsiClpSolverInterface& osi); //Is optimality proven?

//Methods related to querying the input data
int getNumCols(OsiClpSolverInterface& osi);// Get number of columns.
int getNumRows(OsiClpSolverInterface& osi);// Get number of rows.
std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen=std::string::npos);// Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen=std::string::npos);// Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
double getColLower(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of column lower bounds.
double getColUpper(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of column upper bounds.
double getObjCoefficients(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of objective function coefficients.
double getRowLower(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row lower bounds.
double getRowUpper(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row upper bounds.
double getObjSense(OsiClpSolverInterface& osi);// Get objective function sense (1 for min (default), -1 for max).
void setObjSense(OsiClpSolverInterface& osi, double s); // Sets objective function sense (1 for min (default), -1 for max).
bool isContinuous(OsiClpSolverInterface& osi, int colNumber);// Return true if column is continuous.
bool isBinary(OsiClpSolverInterface& osi, int colIndex);// Return true if variable is binary.
bool isInteger(OsiClpSolverInterface& osi, int colIndex);// Return true if column is integer.
const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface& osi);// Get pointer to column-wise copy of matrix.
double getInfinity(OsiClpSolverInterface& osi);// Get solver's value for infinity.

//Methods related to querying the solution
double getColSolution(OsiClpSolverInterface& osi, int colIndex);// Get pointer to array[getNumCols()] of primal solution vector.
double getObjValue(OsiClpSolverInterface& osi);// Get objective function value.

//Changing bounds on variables and constraints
void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set an objective function coefficient.
void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single column lower bound, Use -DBL_MAX for -infinity.
void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single column upper bound, Use DBL_MAX for infinity.
void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper);// Set a single column lower and upper bound.
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single row lower bound, Use -DBL_MAX for -infinity.
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single row upper bound, Use DBL_MAX for infinity.
void setObjective(OsiClpSolverInterface& osi, const double *array);// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name);// Set name of row.
void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name);// Set name of column.


//Integrality related changing methods
void setContinuous(OsiClpSolverInterface& osi, int index);// Set the index-th variable to be a continuous variable.
void setInteger(OsiClpSolverInterface& osi, int index);// Set the index-th variable to be an integer variable.

/*Methods to expand a problem.
  Note that if a column is added then by default
  it will correspond to a continuous variable.  */

void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices);// Remove a set of columns (primal variables) from the problem.
void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero=false);
void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices);// Delete a set of rows (constraints) from the problem.

//Methods to input a problem
int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon=1e-5);// Read file in LP format (with names).
void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into an Lp file of the given filename.

// Resolve an LP relaxation after problem modification.
void resolve(OsiSolverInterface& osi);

// Returns solver - has current state.
OsiClpSolverInterface* solver(CbcModel& m);
// Set current log (detail) level.
void setLogLevel(OsiClpSolverInterface& s, int n);
// Get pointer to array[getNumCols()] of objective function coefficients.
const double* getObjCoefficients(OsiClpSolverInterface& osi);
// Sets the objective coefficient
void setObjCoefficients(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *coefs);
// Sets the coefficients for a row
void setCoefficients(OsiClpSolverInterface& osi, int row, int numberElements, const int *cols, const double *coefs);
// Return the coefficient.
double getCoefficient(OsiClpSolverInterface& osi, int row, int col);
// Sun the solver using CbcMain1
void callCbc0 (CbcModel &m);
int callCbc1 (int argc, const char *argv[], CbcModel & babSolver);
// Get objective function value.
double getObjValue(OsiSolverInterface& c);
// Get pointer to array[getNumCols()] of primal solution vector.
double getColSolution(OsiSolverInterface& osi, int colIndex);

double bestSolution(CbcModel& c, int colIndex);
#endif
