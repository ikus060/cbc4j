%module cbc4j
%{
#include "cbc4j.hpp"
%}
%include "arrays_java.i";
%include "std_string.i";
%include "various.i"
%apply double[] {double *};
%apply int[] {int *};
%apply char **STRING_ARRAY { const char *argv[] }
using namespace std;
typedef std::string String;
%typemap(out) const double* {
$result = SWIG_JavaArrayOutDouble(jenv, $1, arg1->getNumCols());
}

extern void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);
extern void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);
extern double bestSolution(CbcModel& c, int colIndex);
extern void branchAndBound(CbcModel& c, int doStatistics = 0);
extern void callCbc0 (CbcModel &m);
extern int callCbc1 (int argc, const char *argv[], CbcModel & babSolver);
extern void deleteCbcModel(CbcModel* c);
extern void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices);
extern void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices);
extern void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);
extern double getCbcColSolution(CbcModel& c, int colIndex);
extern double getCbcRowLower(CbcModel& c, int rowIndex);
extern double getCbcRowUpper(CbcModel& c, int rowIndex);
extern double getCoefficient(OsiClpSolverInterface& osi, int row, int col);
extern double getColLower(CbcModel& c, int colIndex);
extern double getColLower(OsiClpSolverInterface& osi, int elementIndex);
extern std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen=std::string::npos);
extern double getColSolution(OsiClpSolverInterface& osi, int colIndex);
extern double getColSolution(OsiSolverInterface& osi, int colIndex);
extern double getColUpper(CbcModel& c, int colIndex);
extern double getColUpper(OsiClpSolverInterface& osi, int elementIndex);
extern double getInfinity(CbcModel& c);
extern double getInfinity(OsiClpSolverInterface& osi);
extern double getObjCoefficients(CbcModel& c, int colIndex);
extern double getObjSense(CbcModel& c);
extern int getNumCols(CbcModel& c);
extern int getNumCols(OsiClpSolverInterface& osi);
extern int getNumRows(CbcModel& c);
extern int getNumRows(OsiClpSolverInterface& osi);
extern double getObjCoefficients(OsiClpSolverInterface& osi, int elementIndex);
extern const double* getObjCoefficients(OsiClpSolverInterface& osi);
extern double getObjValue(CbcModel& c);
extern double getObjValue(OsiClpSolverInterface& osi);
extern double getObjValue(OsiSolverInterface& c);
extern double getObjSense(OsiClpSolverInterface& osi);
extern double getRowLower(CbcModel& c, int rowIndex);
extern double getRowLower(OsiClpSolverInterface& osi, int rowIndex);
extern std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen=std::string::npos);
extern double getRowUpper(CbcModel& c, int rowIndex);
extern double getRowUpper(OsiClpSolverInterface& osi, int rowIndex);
extern void initialSolve(CbcModel& c);
extern void initialSolve(OsiClpSolverInterface& solver);
extern bool isBinary(CbcModel& c, int colIndex);
extern bool isBinary(OsiClpSolverInterface& osi, int colIndex);
extern bool isContinuous(CbcModel& c, int colIndex);
extern bool isContinuous(OsiClpSolverInterface& osi, int colNumber);
extern bool isInteger(CbcModel& c, int colIndex);
extern bool isInteger(OsiClpSolverInterface& osi, int colIndex);
extern bool isProvenOptimal(CbcModel& c);
extern bool isProvenOptimal(OsiClpSolverInterface& osi);
extern bool isProvenInfeasible(CbcModel& c);
extern CbcModel* newCbcModel(OsiClpSolverInterface& s);
extern OsiClpSolverInterface* newOsiClpSolverInterface();
extern void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero=false);
extern int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon=1e-5);
extern void resolve(OsiSolverInterface& osi);
extern void setCoefficients(OsiClpSolverInterface& osi, int row, int numberElements, const int *cols, const double *coefs);
extern void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper);
extern void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
extern void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name);
extern void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
extern void setContinuous(OsiClpSolverInterface& osi, int index);
extern void setInteger(OsiClpSolverInterface& osi, int index);
extern void setLogLevel(CbcModel& c, int n);
extern void setLogLevel(OsiClpSolverInterface& s, int n);
extern void setObjective(OsiClpSolverInterface& osi, const double *array);
extern void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
extern void setObjCoefficients(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *coefs);
extern void setObjSense(CbcModel& c, double s);
extern void setObjSense(OsiClpSolverInterface& osi, double s); 
extern void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
extern void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name);
extern void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
extern OsiClpSolverInterface* solver(CbcModel& m);
extern int status(CbcModel& c);
extern void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);
