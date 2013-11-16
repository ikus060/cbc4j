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
// Preprocessing
#include "CglPreProcess.hpp"
#include  "CoinTime.hpp"
#include "CbcHeuristicFPump.hpp"

// Solve initial LP relaxation.
void initialSolve(OsiClpSolverInterface& solver);

// Solve initial LP relaxation.
void initialSolve(OsiSolverInterface& solver);

void setInt(OsiClpSolverInterface& solver, int i);

// Default Constructor. 
OsiClpSolverInterface* newOsiClpSolverInterface();

// Copy constructor.
OsiClpSolverInterface* newOsiClpSolverInterface(OsiClpSolverInterface & s);

// Delete OsiClpSolverInterface
void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);

CoinBigIndex CBI(int start[]);

CoinPackedMatrix CPM(int numberRows, int numberColumns, int numberElements, int *startint, int nbelems, int *length, int *rows, double *elements);
CoinPackedMatrix CPM(int *rowIndices, int *colIndices, double *elements, int nbelems);

// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiClpSolverInterface& s);
// Create a new instance of CbcModel
CbcModel* newCbcModel(OsiSolverInterface& s);
void deleteCbcModel(CbcModel* c);

void initialSolve(CbcModel& c);
void branchAndBound(CbcModel& c, int doStatistics = 0);
int status(CbcModel& c);
int secondaryStatus(CbcModel& c);
bool isProvenOptimal(CbcModel& c);
bool isProvenInfeasible(CbcModel& c);
bool isContinuousUnbounded(CbcModel& c);
bool isNodeLimitReached(CbcModel& c);
bool isSecondsLimitReached(CbcModel& c);
int getNumCols(CbcModel& c);
int getNumRows(CbcModel& c) ;
CoinBigIndex getNumElements (CbcModel& c);
int numberIntegers(CbcModel& c) ;
char integerType(CbcModel& c, int i);
double getColLower(CbcModel& c, int colIndex);
double getColUpper(CbcModel& c, int colIndex);
const char* getRowSense(CbcModel& c);
double getRightHandSide(CbcModel& c, int rowIndex);
double getRowRange(CbcModel& c, int rowIndex);
double getRowLower(CbcModel& c, int rowIndex);
double getRowUpper(CbcModel& c, int rowIndex);
double getObjCoefficients(CbcModel& c, int colIndex);
double getObjSense(CbcModel& c);
void setObjSense(CbcModel& c, double s);
bool isContinuous(CbcModel& c, int colIndex);
bool isBinary(CbcModel& c, int colIndex);
bool isInteger(CbcModel& c, int colIndex);
bool isIntegerNonBinary(CbcModel& c, int colIndex);
bool isFreeBinary(CbcModel& c, int colIndex);
const CoinPackedMatrix* getMatrixByRow(CbcModel& c);
const CoinPackedMatrix* getMatrixByCol(CbcModel& c);
double getInfinity(CbcModel& c);
double getCbcColLower(CbcModel& c, int colIndex);
double getCbcColUpper(CbcModel& c, int colIndex);
double getCbcRowLower(CbcModel& c, int rowIndex);
double getCbcRowUpper(CbcModel& c, int rowIndex);
double getCbcColSolution(CbcModel& c, int colIndex);
double getCbcRowPrice(CbcModel& c, int colIndex);
double getObjValue(CbcModel& c);
double bestSolution(CbcModel& c, int colIndex);
int getNumberThreads(CbcModel& c);
void setNumberThreads(CbcModel& c, int value);
void setLogLevel(CbcModel& c, int n);
double getCurrentObjValue(CbcModel& c);

//OsiClpInterface proxy methods
//Methods returning info on how the solution process terminated
bool isAbandoned(OsiClpSolverInterface& osi); //Are there a numerical difficulties?
bool isProvenOptimal(OsiClpSolverInterface& osi); //Is optimality proven?
bool isProvenPrimalInfeasible(OsiClpSolverInterface& osi);// Is primal infeasiblity proven?
bool isProvenDualInfeasible(OsiClpSolverInterface& osi);// Is dual infeasiblity proven?
bool isPrimalObjectiveLimitReached(OsiClpSolverInterface& osi);// Is the given primal objective limit reached?
bool isDualObjectiveLimitReached(OsiClpSolverInterface& osi);// Is the given dual objective limit reached?
bool isIterationLimitReached(OsiClpSolverInterface& osi);// Iteration limit reached?

//Methods related to querying the input data
int getNumCols(OsiClpSolverInterface& osi);// Get number of columns.
int getNumRows(OsiClpSolverInterface& osi);// Get number of rows.
int getNumElements(OsiClpSolverInterface& osi);// Get number of nonzero elements.
std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen=std::string::npos);// Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen=std::string::npos);// Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
double getColLower(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of column lower bounds.
double getColUpper(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of column upper bounds.
double getObjCoefficients(OsiClpSolverInterface& osi, int elementIndex);// Get pointer to array[getNumCols(OsiClpSolverInterface osi)] of objective function coefficients.
const char * getRowSense(OsiClpSolverInterface& osi);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row constraint senses.
double getRightHandSide(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of rows right-hand sides.
double getRowRange(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row ranges.
double getRowLower(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row lower bounds.
double getRowUpper(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows(OsiClpSolverInterface osi)] of row upper bounds.
double getObjSense(OsiClpSolverInterface& osi);// Get objective function sense (1 for min (default), -1 for max).
void setObjSense(OsiClpSolverInterface& osi, double s); // Sets objective function sense (1 for min (default), -1 for max).
bool isContinuous(OsiClpSolverInterface& osi, int colNumber);// Return true if column is continuous.
bool isBinary(OsiClpSolverInterface& osi, int colIndex);// Return true if variable is binary.
bool isInteger(OsiClpSolverInterface& osi, int colIndex);// Return true if column is integer.
bool isIntegerNonBinary(OsiClpSolverInterface& osi, int colIndex);// Return true if variable is general integer.
bool isFreeBinary(OsiClpSolverInterface& osi, int colIndex);// Return true if variable is binary and not fixed at either bound.
const char * getColType(OsiClpSolverInterface& osi, bool refresh=false);// Return array of column length 0 - continuous 1 - binary (may get fixed later) 2 - general integer (may get fixed later).
const CoinPackedMatrix * getMatrixByRow(OsiClpSolverInterface& osi);// Get pointer to row-wise copy of matrix.
const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface& osi);// Get pointer to column-wise copy of matrix.
CoinPackedMatrix *getMutableMatrixByCol(OsiClpSolverInterface& osi);// Get pointer to mutable column-wise copy of matrix.
double getInfinity(OsiClpSolverInterface& osi);// Get solver's value for infinity.

//Methods related to querying the solution
double getColSolution(OsiClpSolverInterface& osi, int colIndex);// Get pointer to array[getNumCols()] of primal solution vector.
double getRowPrice(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows()] of dual prices.
double getReducedCost(OsiClpSolverInterface& osi, int colIndex);// Get a pointer to array[getNumCols()] of reduced costs.
double getRowActivity(OsiClpSolverInterface& osi, int rowIndex);// Get pointer to array[getNumRows()] of row activity levels (constraint matrix times the solution vector.
double getObjValue(OsiClpSolverInterface& osi);// Get objective function value.
int getIterationCount(OsiClpSolverInterface& osi);// Get how many iterations it took to solve the problem (whatever "iteration" mean to the solver.
std::vector< double * > getDualRays(OsiClpSolverInterface& osi, int maxNumRays);// Get as many dual rays as the solver can provide.
std::vector< double * > getPrimalRays(OsiClpSolverInterface& osi, int maxNumRays);// Get as many primal rays as the solver can provide.

//Changing bounds on variables and constraints
void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set an objective function coefficient.
void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single column lower bound, Use -DBL_MAX for -infinity.
void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single column upper bound, Use DBL_MAX for infinity.
void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper);// Set a single column lower and upper bound.
void setColSetBounds(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const double *boundList);// Set the bounds on a number of columns simultaneously, The default implementation just invokes setColLower() and setColUpper() over and over again.
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single row lower bound, Use -DBL_MAX for -infinity.
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);// Set a single row upper bound, Use DBL_MAX for infinity.
void setRowBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper);// Set a single row lower and upper bound.
void setRowType(OsiClpSolverInterface& osi, int index, char sense, double rightHandSide, double range);// Set the type of a single row.
void setRowSetBounds(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const double *boundList);// Set the bounds on a number of rows simultaneously. The default implementation just invokes setRowLower() and setRowUpper() over and over again.
void setRowSetTypes(OsiClpSolverInterface& osi, const int *indexFirst, const int *indexLast, const char *senseList, const double *rhsList, const double *rangeList);// Set the type of a number of rows simultaneously. The default implementation just invokes setRowType() over and over again.
void setObjective(OsiClpSolverInterface& osi, const double *array);// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
void setColLower(OsiClpSolverInterface& osi, const double *array);// Set the lower bounds for all columns array [getNumCols()] is an array of values for the objective.
void setColUpper(OsiClpSolverInterface& osi, const double *array);// Set the upper bounds for all columns array [getNumCols()] is an array of values for the objective.
void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name);// Set name of row.
void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name);// Set name of column.


//Integrality related changing methods
void setContinuous(OsiClpSolverInterface& osi, int index);// Set the index-th variable to be a continuous variable.
void setInteger(OsiClpSolverInterface& osi, int index);// Set the index-th variable to be an integer variable.
void setContinuous(OsiClpSolverInterface& osi, const int *indices, int len);// Set the variables listed in indices (which is of length len) to be continuous variables.
void setInteger(OsiClpSolverInterface& osi, const int *indices, int len);// Set the variables listed in indices (which is of length len) to be integer variables.


/*Methods to expand a problem.
  Note that if a column is added then by default
  it will correspond to a continuous variable.  */

void addCol(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
void addCols(OsiClpSolverInterface& osi, const int numcols, const CoinPackedVectorBase *const *cols, const double *collb, const double *colub, const double *obj);// Add a set of columns (primal variables) to the problem.
void addCols(OsiClpSolverInterface& osi, const int numcols, const int *columnStarts, const int *rows, const double *elements, const double *collb, const double *colub, const double *obj);// Add a set of columns (primal variables) to the problem.
void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices);// Remove a set of columns (primal variables) from the problem.
void addRow(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
void addRow(OsiClpSolverInterface& osi, const CoinPackedVectorBase &vec, const char rowsen, const double rowrhs, const double rowrng);// Add a row (constraint) to the problem.
void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
void addRows(OsiClpSolverInterface& osi, const int numrows, const CoinPackedVectorBase *const *rows, const double *rowlb, const double *rowub);
void addRows(OsiClpSolverInterface& osi, const int numrows, const CoinPackedVectorBase *const *rows, const char *rowsen, const double *rowrhs, const double *rowrng);// Add a set of rows (constraints) to the problem.
void addRows(OsiClpSolverInterface& osi, const int numrows, const int *rowStarts, const int *columns, const double *element, const double *rowlb, const double *rowub);// Add a set of rows (constraints) to the problem.
void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero=false);
void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices);// Delete a set of rows (constraints) from the problem.
void saveBaseModel(OsiClpSolverInterface& osi);// If solver wants it can save a copy of "base" (continuous) model here.
void restoreBaseModel(OsiClpSolverInterface& osi, int numberRows);// Strip off rows to get to this number of rows.
void applyRowCuts(OsiClpSolverInterface& osi, int numberCuts, const OsiRowCut *cuts);// Apply a collection of row cuts which are all effective.
void applyRowCuts(OsiClpSolverInterface& osi, int numberCuts, const OsiRowCut **cuts);// Apply a collection of row cuts which are all effective.


//Methods to input a problem
void loadProblem(OsiClpSolverInterface& osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub);// Load in an problem by copying the arguments (the constraints on the rows are given by lower and upper bounds).
void assignProblem(OsiClpSolverInterface& osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, double *&rowlb, double *&rowub);// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by lower and upper bounds).
void loadProblem(OsiClpSolverInterface& osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng);// Load in an problem by copying the arguments (the constraints on the rows are given by sense/rhs/range triplets).
void assignProblem(OsiClpSolverInterface& osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, char *&rowsen, double *&rowrhs, double *&rowrng);// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by sense/rhs/range triplets).
void loadProblem(OsiClpSolverInterface& osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub);// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
void loadProblem(OsiClpSolverInterface& osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng);// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
int loadFromCoinModel(OsiClpSolverInterface& osi, CoinModel &modelObject, bool keepSolution=false);// This loads a model from a coinModel object - returns number of errors.
int readMps(OsiClpSolverInterface& osi, const char *filename, const char *extension="mps");// Read an mps file from the given filename (defaults to Osi reader) - returns number of errors (see OsiMpsReader class).
int readMps(OsiClpSolverInterface& osi, const char *filename, bool keepNames, bool allowErrors);// Read an mps file from the given filename returns number of errors (see OsiMpsReader class).
void writeMps(OsiClpSolverInterface& osi, const char *filename, const char *extension="mps", double objSense=0.0);// Write the problem into an mps file of the given filename.
int writeMpsNative(OsiClpSolverInterface& osi, const char *filename, const char **rowNames, const char **columnNames, int formatType=0, int numberAcross=2, double objSense=0.0);// Write the problem into an mps file of the given filename, names may be null.
int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon=1e-5);// Read file in LP format (with names).
void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into an Lp file of the given filename.
void writeLp(OsiClpSolverInterface& osi, FILE *fp, double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into the file pointed to by the parameter fp.

// Create a new pre-processor
CglPreProcess* newCglPreProcess();

// Delete the pre-processor
void deleteCglPreProcess(CglPreProcess* p);

// Run the preprocessor
OsiSolverInterface* preProcess(CglPreProcess& p, OsiClpSolverInterface& osi, bool makeEquality, int numberPasses);

// Create solution in original model
void postProcess(CglPreProcess& p, OsiSolverInterface& osi);

// Resolve an LP relaxation after problem modification.
void resolve(OsiSolverInterface& osi);

// Returns solver - has current state. 
OsiClpSolverInterface* solver(CbcModel& m);
// Create a new CbcHeuristicFPump
CbcHeuristicFPump* newCbcHeuristicFPump();
// Create a new CbcHeuristicFPump
CbcHeuristicFPump* newCbcHeuristicFPump(CbcModel& m, double downValue, bool roundExpensive );
// Delete CbcHeuristicFPump
void deleteCbcHeuristicFPump(CbcHeuristicFPump * h);
// Add one heuristic - up to user to delete. 
void addHeuristic(CbcModel& m, CbcHeuristicFPump* h);
// Delete OsiSolverInterface
void deleteOsiSolverInterface(OsiSolverInterface* s);
// Set current log (detail) level. 
void setLogLevel(CglPreProcess& p, int n);
// Set current log (detail) level.
void setLogLevel(OsiSolverInterface& s, int n);
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
#endif
