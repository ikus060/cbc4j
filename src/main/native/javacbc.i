%module javacbc
%{
#include "javacbc.hpp"
%}

%include "arrays_java.i";
%include "std_string.i";
%include "various.i"

%apply double[] {double *};
%apply int[] {int *};
%apply char **STRING_ARRAY { const char *argv[] }

using namespace std;
typedef std::string String;

// Solve initial LP relaxation.
extern void initialSolve(OsiClpSolverInterface &solver);

// Solve initial LP relaxation.
extern void initialSolve(OsiSolverInterface& solver);

extern void setInt(OsiClpSolverInterface &solver, int i);

// Default Constructor.
extern OsiClpSolverInterface* newOsiClpSolverInterface();

// Copy constructor.
extern OsiClpSolverInterface* newOsiClpSolverInterface(OsiClpSolverInterface & s);

// Delete OsiClpSolverInterface
extern void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);

extern CoinBigIndex CBI(int start[]);
extern CoinPackedMatrix CPM(int numberRows, int numberColumns, int numberElements, int *startint, int nbelems, int *length, int *rows, double *elements);
extern CoinPackedMatrix CPM(int *rowIndices, int *colIndices, double *elements, int nbelems);

// Create a new instance of CbcModel
extern CbcModel* newCbcModel(OsiClpSolverInterface &s);
// Create a new instance of CbcModel
extern CbcModel* newCbcModel(OsiSolverInterface& s);
extern void deleteCbcModel(CbcModel* c);
extern void initialSolve(CbcModel &c);
extern void branchAndBound(CbcModel &c, int doStatistics = 0);
extern int status(CbcModel &c);
extern int secondaryStatus(CbcModel &c);
extern bool isProvenOptimal(CbcModel &c);
extern bool isProvenInfeasible(CbcModel &c);
extern bool isContinuousUnbounded(CbcModel &c);
extern bool isNodeLimitReached(CbcModel &c);
extern bool isSecondsLimitReached(CbcModel &c);
extern int getNumCols(CbcModel &c);
extern int getNumRows(CbcModel &c) ;
extern CoinBigIndex getNumElements(CbcModel &c);
extern int numberIntegers(CbcModel &c) ;
extern char integerType(CbcModel &c,int i);

extern double getColLower(CbcModel &c, int colIndex);
extern double getColUpper(CbcModel &c, int colIndex);
extern double getObjCoefficients(CbcModel &c, int colIndex);
extern double getCbcColLower(CbcModel &c, int colIndex);
extern double getCbcColUpper(CbcModel &c, int colIndex);
extern double getCbcColSolution(CbcModel &c, int colIndex);

extern double getRightHandSide(CbcModel &c, int rowIndex);
extern double getRowRange(CbcModel &c, int rowIndex);
extern double getRowLower(CbcModel &c, int rowIndex);
extern double getRowUpper(CbcModel &c, int rowIndex);
extern double getCbcRowLower(CbcModel &c, int rowIndex);
extern double getCbcRowUpper(CbcModel &c, int rowIndex);
extern double getCbcRowPrice(CbcModel &c, int rowIndex);

extern double bestSolution(CbcModel &c, int colIndex);
double getCurrentObjValue(CbcModel &c);

extern const char* getRowSense(CbcModel &c);
extern double getObjSense(CbcModel &c);
extern void setObjSense(CbcModel &c, double s);
extern bool isContinuous(CbcModel &c, int colIndex);
extern bool isBinary(CbcModel &c, int colIndex);
extern bool isInteger(CbcModel &c, int colIndex);
extern bool isIntegerNonBinary(CbcModel &c, int colIndex);
extern bool isFreeBinary(CbcModel &c, int colIndex);
extern const CoinPackedMatrix* getMatrixByRow(CbcModel &c);
extern const CoinPackedMatrix* getMatrixByCol(CbcModel &c);
extern double getInfinity(CbcModel &c);
extern double getObjValue(CbcModel &c);
extern int getNumberThreads(CbcModel &c);
extern void setNumberThreads(CbcModel &c, int value);
extern void setLogLevel(CbcModel &c, int n);

//OsiClpSolverInterface proxy functions and constructor
//Methods returning info on how the solution process terminated
extern bool isAbandoned(OsiClpSolverInterface &osi); //Are there a numerical difficulties?
extern bool isProvenOptimal(OsiClpSolverInterface &osi); //Is optimality proven?
extern bool isProvenPrimalInfeasible(OsiClpSolverInterface &osi);// Is primal infeasiblity proven?
extern bool isProvenDualInfeasible(OsiClpSolverInterface &osi);// Is dual infeasiblity proven?
extern bool isPrimalObjectiveLimitReached(OsiClpSolverInterface &osi);// Is the given primal objective limit reached?
extern bool isDualObjectiveLimitReached(OsiClpSolverInterface &osi);// Is the given dual objective limit reached?
extern bool isIterationLimitReached(OsiClpSolverInterface &osi);// Iteration limit reached?

//Methods related to querying the input data
extern int getNumCols(OsiClpSolverInterface &osi);// Get number of columns.
extern int getNumRows(OsiClpSolverInterface &osi);// Get number of rows.
extern int getNumElements(OsiClpSolverInterface &osi);// Get number of nonzero elements.
extern std::string getRowName(OsiClpSolverInterface &osi, int rowIndex, unsigned maxLen=std::string::npos);// Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
extern std::string getColName(OsiClpSolverInterface &osi, int colIndex, unsigned maxLen=std::string::npos);// Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!

extern double getColLower(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of column lower bounds.
extern double getColUpper(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of column upper bounds.
extern double getObjCoefficients(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of objective function coefficients.

extern double getRightHandSide(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of rows right-hand sides.
extern double getRowRange(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row ranges.
extern double getRowLower(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row lower bounds.
extern double getRowUpper(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row upper bounds.

extern const char * getRowSense(OsiClpSolverInterface &osi);// Get pointer to array[getNumRows()] of row constraint senses.

extern double getObjSense(OsiClpSolverInterface &osi);// Get objective function sense (1 for min (default), -1 for max).
extern void setObjSense(OsiClpSolverInterface &osi, double s);// Sets objective function sense (1 for min (default), -1 for max).
extern bool isContinuous(OsiClpSolverInterface &osi, int colNumber);// Return true if column is continuous.
extern bool isBinary(OsiClpSolverInterface &osi, int colIndex);// Return true if variable is binary.
extern bool isInteger(OsiClpSolverInterface &osi, int colIndex);// Return true if column is integer.
extern bool isIntegerNonBinary(OsiClpSolverInterface &osi, int colIndex);// Return true if variable is general integer.
extern bool isFreeBinary(OsiClpSolverInterface &osi, int colIndex);// Return true if variable is binary and not fixed at either bound.
extern const char * getColType(OsiClpSolverInterface &osi, bool refresh=false);// Return array of column length 0 - continuous 1 - binary (may get fixed later) 2 - general integer (may get fixed later).
extern const CoinPackedMatrix * getMatrixByRow(OsiClpSolverInterface &osi);// Get pointer to row-wise copy of matrix.
extern const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface &osi);// Get pointer to column-wise copy of matrix.
extern CoinPackedMatrix *getMutableMatrixByCol(OsiClpSolverInterface &osi);// Get pointer to mutable column-wise copy of matrix.
extern double getInfinity(OsiClpSolverInterface &osi);// Get solver's value for infinity.

//Methods related to querying the solution
extern double getColSolution(OsiClpSolverInterface &osi, int colIndex);// Get pointer to array[getNumCols()] of primal solution vector.
extern double getReducedCost(OsiClpSolverInterface &osi, int colIndex);// Get a pointer to array[getNumCols()] of reduced costs.

extern double getRowPrice(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of dual prices.
extern double getRowActivity(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row activity levels (constraint matrix times the solution vector.

extern double getObjValue(OsiClpSolverInterface &osi);// Get objective function value.
extern int getIterationCount(OsiClpSolverInterface &osi);// Get how many iterations it took to solve the problem (whatever "iteration" mean to the solver.
extern std::vector< double * > getDualRays(OsiClpSolverInterface &osi, int maxNumRays);// Get as many dual rays as the solver can provide.
extern std::vector< double * > getPrimalRays(OsiClpSolverInterface &osi, int maxNumRays);// Get as many primal rays as the solver can provide.


//Changing bounds on variables and constraints
void setObjCoeff(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set an objective function coefficient.
void setColLower(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single column lower bound, Use -DBL_MAX for -infinity.
void setColUpper(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single column upper bound, Use DBL_MAX for infinity.
void setColBounds(OsiClpSolverInterface &osi, int elementIndex, double lower, double upper);// Set a single column lower and upper bound.
void setColSetBounds(OsiClpSolverInterface &osi, const int *indexFirst, const int *indexLast, const double *boundList);// Set the bounds on a number of columns simultaneously, The default implementation just invokes setColLower() and setColUpper() over and over again.
void setRowLower(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single row lower bound, Use -DBL_MAX for -infinity.
void setRowUpper(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single row upper bound, Use DBL_MAX for infinity.
void setRowBounds(OsiClpSolverInterface &osi, int elementIndex, double lower, double upper);// Set a single row lower and upper bound.
void setRowType(OsiClpSolverInterface &osi, int index, char sense, double rightHandSide, double range);// Set the type of a single row.
void setRowSetBounds(OsiClpSolverInterface &osi, const int *indexFirst, const int *indexLast, const double *boundList);// Set the bounds on a number of rows simultaneously. The default implementation just invokes setRowLower() and setRowUpper() over and over again.
void setRowSetTypes(OsiClpSolverInterface &osi, const int *indexFirst, const int *indexLast, const char *senseList, const double *rhsList, const double *rangeList);// Set the type of a number of rows simultaneously. The default implementation just invokes setRowType() over and over again.
void setObjective(OsiClpSolverInterface &osi, const double *array);// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
void setColLower(OsiClpSolverInterface &osi, const double *array);// Set the lower bounds for all columns array [getNumCols()] is an array of values for the objective.
void setColUpper(OsiClpSolverInterface &osi, const double *array);// Set the upper bounds for all columns array [getNumCols()] is an array of values for the objective.
void setRowName(OsiClpSolverInterface &osi, int rowIndex, std::string name);// Set name of row.
void setColName(OsiClpSolverInterface &osi, int colIndex, std::string name);// Set name of column.


//Integrality related changing methods
void setContinuous(OsiClpSolverInterface &osi, int index);// Set the index-th variable to be a continuous variable.
void setInteger(OsiClpSolverInterface &osi, int index);// Set the index-th variable to be an integer variable.
void setContinuous(OsiClpSolverInterface &osi, const int *indices, int len);// Set the variables listed in indices (which is of length len) to be continuous variables.
void setInteger(OsiClpSolverInterface &osi, const int *indices, int len);// Set the variables listed in indices (which is of length len) to be integer variables.



//Methods to expand a problem.
//Note that if a column is added then by default
//it will correspond to a continuous variable.

extern void addCol(OsiClpSolverInterface &osi, const CoinPackedVectorBase &vec, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
extern void addCol(OsiClpSolverInterface &osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
extern void addCols(OsiClpSolverInterface &osi, const int numcols, const CoinPackedVectorBase *const *cols, const double *collb, const double *colub, const double *obj);// Add a set of columns (primal variables) to the problem.
extern void addCols(OsiClpSolverInterface &osi, const int numcols, const int *columnStarts, const int *rows, const double *elements, const double *collb, const double *colub, const double *obj);// Add a set of columns (primal variables) to the problem.
extern void deleteCols(OsiClpSolverInterface &osi, const int num, const int *colIndices);// Remove a set of columns (primal variables) from the problem.
extern void addRow(OsiClpSolverInterface &osi, const CoinPackedVectorBase &vec, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
extern void addRow(OsiClpSolverInterface &osi, const CoinPackedVectorBase &vec, const char rowsen, const double rowrhs, const double rowrng);// Add a row (constraint) to the problem.
extern void addRow(OsiClpSolverInterface &osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
extern void addRows(OsiClpSolverInterface &osi, const int numrows, const CoinPackedVectorBase *const *rows, const double *rowlb, const double *rowub);
extern void addRows(OsiClpSolverInterface &osi, const int numrows, const CoinPackedVectorBase *const *rows, const char *rowsen, const double *rowrhs, const double *rowrng);// Add a set of rows (constraints) to the problem.
extern void addRows(OsiClpSolverInterface &osi, const int numrows, const int *rowStarts, const int *columns, const double *element, const double *rowlb, const double *rowub);// Add a set of rows (constraints) to the problem.
extern void modifyCoefficient(OsiClpSolverInterface &osi, int row, int column, double newElement, bool keepZero=false);
extern void deleteRows(OsiClpSolverInterface &osi, const int num, const int *rowIndices);// Delete a set of rows (constraints) from the problem.
extern void saveBaseModel(OsiClpSolverInterface &osi);// If solver wants it can save a copy of "base" (continuous) model here.
extern void restoreBaseModel(OsiClpSolverInterface &osi, int numberRows);// Strip off rows to get to this number of rows.
extern void applyRowCuts(OsiClpSolverInterface &osi, int numberCuts, const OsiRowCut *cuts);// Apply a collection of row cuts which are all effective.
extern void applyRowCuts(OsiClpSolverInterface &osi, int numberCuts, const OsiRowCut **cuts);// Apply a collection of row cuts which are all effective.

//Methods to input a problem
extern void loadProblem(OsiClpSolverInterface &osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub);// Load in an problem by copying the arguments (the constraints on the rows are given by lower and upper bounds).
extern void assignProblem(OsiClpSolverInterface &osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, double *&rowlb, double *&rowub);// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by lower and upper bounds).
extern void loadProblem(OsiClpSolverInterface &osi, const CoinPackedMatrix &matrix, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng);// Load in an problem by copying the arguments (the constraints on the rows are given by sense/rhs/range triplets).
extern void assignProblem(OsiClpSolverInterface &osi, CoinPackedMatrix *&matrix, double *&collb, double *&colub, double *&obj, char *&rowsen, double *&rowrhs, double *&rowrng);// Load in an problem by assuming ownership of the arguments (the constraints on the rows are given by sense/rhs/range triplets).
extern void loadProblem(OsiClpSolverInterface &osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const double *rowlb, const double *rowub);// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
extern void loadProblem(OsiClpSolverInterface &osi, const int numcols, const int numrows, const CoinBigIndex *start, const int *index, const double *value, const double *collb, const double *colub, const double *obj, const char *rowsen, const double *rowrhs, const double *rowrng);// Just like the other loadProblem() methods except that the matrix is given in a standard column major ordered format (without gaps).
extern int loadFromCoinModel(OsiClpSolverInterface &osi, CoinModel &modelObject, bool keepSolution=false);// This loads a model from a coinModel object - returns number of errors.
extern int readMps(OsiClpSolverInterface &osi, const char *filename, const char *extension="mps");// Read an mps file from the given filename (defaults to Osi reader) - returns number of errors (see OsiMpsReader class).
extern int readMps(OsiClpSolverInterface &osi, const char *filename, bool keepNames, bool allowErrors);// Read an mps file from the given filename returns number of errors (see OsiMpsReader class).
extern void writeMps(OsiClpSolverInterface &osi, const char *filename, const char *extension="mps", double objSense=0.0);// Write the problem into an mps file of the given filename.
extern int writeMpsNative(OsiClpSolverInterface &osi, const char *filename, const char **rowNames, const char **columnNames, int formatType=0, int numberAcross=2, double objSense=0.0);// Write the problem into an mps file of the given filename, names may be null.
extern int readLp(OsiClpSolverInterface &osi, const char *filename, const double epsilon=1e-5);// Read file in LP format (with names).
extern void writeLp(OsiClpSolverInterface &osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into an Lp file of the given filename.
extern void writeLp(OsiClpSolverInterface &osi, FILE *fp, double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into the file pointed to by the parameter fp.

// Create a new pre-processor
extern CglPreProcess* newCglPreProcess();

// Delete the pre-processor
extern void deleteCglPreProcess(CglPreProcess* p);

// Run the preprocessor
extern OsiSolverInterface* preProcess(CglPreProcess& p, OsiClpSolverInterface &osi, bool makeEquality, int numberPasses);

// Create solution in original model
extern void postProcess(CglPreProcess& p, OsiSolverInterface& osi);

// Resolve an LP relaxation after problem modification.
extern void resolve(OsiSolverInterface& osi);

// Returns solver - has current state. 
extern OsiClpSolverInterface* solver(CbcModel &m);
// Create a new CbcHeuristicFPump
extern CbcHeuristicFPump* newCbcHeuristicFPump();
// Create a new CbcHeuristicFPump
extern CbcHeuristicFPump* newCbcHeuristicFPump(CbcModel &m, double downValue, bool roundExpensive );
// Delete CbcHeuristicFPump
extern void deleteCbcHeuristicFPump(CbcHeuristicFPump * h);
// Add one heuristic - up to user to delete. 
extern void addHeuristic(CbcModel &m, CbcHeuristicFPump* h);
// Delete OsiSolverInterface
extern void deleteOsiSolverInterface(OsiSolverInterface* s);
// Set current log (detail) level. 
extern void setLogLevel(CglPreProcess& p, int n);
// Set current log (detail) level.
extern void setLogLevel(OsiSolverInterface& s, int n);
// Set current log (detail) level.
extern void setLogLevel(OsiClpSolverInterface &s, int n);
// Get pointer to array[getNumCols()] of objective function coefficients. 
%typemap(out) const double* {
$result = SWIG_JavaArrayOutDouble(jenv, $1, arg1->getNumCols());
}
extern const double* getObjCoefficients(OsiClpSolverInterface &osi);
// Sets the objective coefficient
extern void setObjCoefficients(OsiClpSolverInterface &osi, int numberElements, const int *rows, const double *coefs);
// Sets the coefficients for a row
extern void setCoefficients(OsiClpSolverInterface &osi, int row, int numberElements, const int *cols, const double *coefs);
// Return the coefficient.
extern double getCoefficient(OsiClpSolverInterface &osi, int row, int col);
// Run the solver using CbcMain0
extern void callCbc0 (CbcModel &babSolver);
// Run the solver using CbcMain1
extern int callCbc1 (int argc, const char *argv[], CbcModel &babSolver);
// Get objective function value.
extern double getObjValue(OsiSolverInterface& c);
// Get pointer to array[getNumCols()] of primal solution vector.
extern double getColSolution(OsiSolverInterface& osi, int colIndex);
