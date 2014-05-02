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

extern void setInt(OsiClpSolverInterface &solver, int i);

// Default Constructor.
extern OsiClpSolverInterface* newOsiClpSolverInterface();

// Delete OsiClpSolverInterface
extern void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);

// Create a new instance of CbcModel
extern CbcModel* newCbcModel(OsiClpSolverInterface &s);
extern void deleteCbcModel(CbcModel* c);
extern void initialSolve(CbcModel &c);
extern void branchAndBound(CbcModel &c, int doStatistics = 0);
extern int status(CbcModel &c);
extern bool isProvenOptimal(CbcModel &c);
extern bool isProvenInfeasible(CbcModel &c);
extern int getNumCols(CbcModel &c);
extern int getNumRows(CbcModel &c) ;

extern double getColLower(CbcModel &c, int colIndex);
extern double getColUpper(CbcModel &c, int colIndex);
extern double getObjCoefficients(CbcModel &c, int colIndex);
extern double getCbcColSolution(CbcModel &c, int colIndex);

extern double getRowLower(CbcModel &c, int rowIndex);
extern double getRowUpper(CbcModel &c, int rowIndex);
extern double getCbcRowLower(CbcModel &c, int rowIndex);
extern double getCbcRowUpper(CbcModel &c, int rowIndex);

extern double getObjSense(CbcModel &c);
extern void setObjSense(CbcModel &c, double s);
extern bool isContinuous(CbcModel &c, int colIndex);
extern bool isBinary(CbcModel &c, int colIndex);
extern bool isInteger(CbcModel &c, int colIndex);
extern const CoinPackedMatrix* getMatrixByCol(CbcModel &c);
extern double getInfinity(CbcModel &c);
extern double getObjValue(CbcModel &c);
extern void setLogLevel(CbcModel &c, int n);

//OsiClpSolverInterface proxy functions and constructor
//Methods returning info on how the solution process terminated
extern bool isProvenOptimal(OsiClpSolverInterface &osi); //Is optimality proven?

//Methods related to querying the input data
extern int getNumCols(OsiClpSolverInterface &osi);// Get number of columns.
extern int getNumRows(OsiClpSolverInterface &osi);// Get number of rows.
extern std::string getRowName(OsiClpSolverInterface &osi, int rowIndex, unsigned maxLen=std::string::npos);// Return name of row if one exists or Rnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!
extern std::string getColName(OsiClpSolverInterface &osi, int colIndex, unsigned maxLen=std::string::npos);// Return name of column if one exists or Cnnnnnnn maxLen is currently ignored and only there to match the signature from the base class!

extern double getColLower(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of column lower bounds.
extern double getColUpper(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of column upper bounds.
extern double getObjCoefficients(OsiClpSolverInterface &osi, int elementIndex); // Get pointer to array[getNumCols()] of objective function coefficients.

extern double getRowLower(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row lower bounds.
extern double getRowUpper(OsiClpSolverInterface &osi, int rowIndex);// Get pointer to array[getNumRows()] of row upper bounds.

extern double getObjSense(OsiClpSolverInterface &osi);// Get objective function sense (1 for min (default), -1 for max).
extern void setObjSense(OsiClpSolverInterface &osi, double s);// Sets objective function sense (1 for min (default), -1 for max).
extern bool isContinuous(OsiClpSolverInterface &osi, int colNumber);// Return true if column is continuous.
extern bool isBinary(OsiClpSolverInterface &osi, int colIndex);// Return true if variable is binary.
extern bool isInteger(OsiClpSolverInterface &osi, int colIndex);// Return true if column is integer.
extern const CoinPackedMatrix * getMatrixByCol(OsiClpSolverInterface &osi);// Get pointer to column-wise copy of matrix.
extern double getInfinity(OsiClpSolverInterface &osi);// Get solver's value for infinity.

//Methods related to querying the solution
extern double getColSolution(OsiClpSolverInterface &osi, int colIndex);// Get pointer to array[getNumCols()] of primal solution vector.
extern double getObjValue(OsiClpSolverInterface &osi);// Get objective function value.

//Changing bounds on variables and constraints
void setObjCoeff(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set an objective function coefficient.
void setColLower(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single column lower bound, Use -DBL_MAX for -infinity.
void setColUpper(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single column upper bound, Use DBL_MAX for infinity.
void setColBounds(OsiClpSolverInterface &osi, int elementIndex, double lower, double upper);// Set a single column lower and upper bound.
void setRowLower(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single row lower bound, Use -DBL_MAX for -infinity.
void setRowUpper(OsiClpSolverInterface &osi, int elementIndex, double elementValue);// Set a single row upper bound, Use DBL_MAX for infinity.
void setObjective(OsiClpSolverInterface &osi, const double *array);// Set the objective coefficients for all columns array [getNumCols()] is an array of values for the objective.
void setRowName(OsiClpSolverInterface &osi, int rowIndex, std::string name);// Set name of row.
void setColName(OsiClpSolverInterface &osi, int colIndex, std::string name);// Set name of column.


//Integrality related changing methods
void setContinuous(OsiClpSolverInterface &osi, int index);// Set the index-th variable to be a continuous variable.
void setInteger(OsiClpSolverInterface &osi, int index);// Set the index-th variable to be an integer variable.

//Methods to expand a problem.
//Note that if a column is added then by default
//it will correspond to a continuous variable.

extern void addCol(OsiClpSolverInterface &osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);// Add a column (primal variable) to the problem.
extern void deleteCols(OsiClpSolverInterface &osi, const int num, const int *colIndices);// Remove a set of columns (primal variables) from the problem.
extern void addRow(OsiClpSolverInterface &osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);// Add a row (constraint) to the problem.
extern void modifyCoefficient(OsiClpSolverInterface &osi, int row, int column, double newElement, bool keepZero=false);
extern void deleteRows(OsiClpSolverInterface &osi, const int num, const int *rowIndices);// Delete a set of rows (constraints) from the problem.

//Methods to input a problem
extern int readLp(OsiClpSolverInterface &osi, const char *filename, const double epsilon=1e-5);// Read file in LP format (with names).
extern void writeLp(OsiClpSolverInterface &osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);// Write the problem into an Lp file of the given filename.

// Resolve an LP relaxation after problem modification.
extern void resolve(OsiSolverInterface& osi);

// Returns solver - has current state. 
extern OsiClpSolverInterface* solver(CbcModel &m);
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

extern double bestSolution(CbcModel& c, int colIndex);
