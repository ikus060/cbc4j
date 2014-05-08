#ifndef CBC4J_H
#define CBC4J_H
#include <cassert>
#include <iomanip>
#include <string>
using namespace std;
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "OsiClpSolverInterface.hpp"
#include  "CoinTime.hpp"
#include "CbcHeuristicFPump.hpp"

void addCol(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *elements, const double collb, const double colub, const double obj);
void addRow(OsiClpSolverInterface& osi, int numberElements, const int *columns, const double *element, const double rowlb, const double rowub);
double bestSolution(CbcModel& c, int colIndex);
void branchAndBound(CbcModel& c, int doStatistics = 0);
void callCbc0 (CbcModel &m);
int callCbc1 (int argc, const char *argv[], CbcModel & babSolver);
void deleteCbcModel(CbcModel* c);
void deleteCols(OsiClpSolverInterface& osi, const int num, const int *colIndices);
void deleteRows(OsiClpSolverInterface& osi, const int num, const int *rowIndices);
void deleteOsiClpSolverInterface(OsiClpSolverInterface* solver);
double getCbcColSolution(CbcModel& c, int colIndex);
double getCbcRowLower(CbcModel& c, int rowIndex);
double getCbcRowUpper(CbcModel& c, int rowIndex);
double getCoefficient(OsiClpSolverInterface& osi, int row, int col);
double getColLower(CbcModel& c, int colIndex);
double getColLower(OsiClpSolverInterface& osi, int elementIndex);
std::string getColName(OsiClpSolverInterface& osi, int colIndex, unsigned maxLen=std::string::npos);
double getColSolution(OsiClpSolverInterface& osi, int colIndex);
double getColSolution(OsiSolverInterface& osi, int colIndex);
double getColUpper(CbcModel& c, int colIndex);
double getColUpper(OsiClpSolverInterface& osi, int elementIndex);
double getInfinity(CbcModel& c);
double getInfinity(OsiClpSolverInterface& osi);
double getObjCoefficients(CbcModel& c, int colIndex);
double getObjSense(CbcModel& c);
int getNumCols(CbcModel& c);
int getNumCols(OsiClpSolverInterface& osi);
int getNumRows(CbcModel& c);
int getNumRows(OsiClpSolverInterface& osi);
double getObjCoefficients(OsiClpSolverInterface& osi, int elementIndex);
const double* getObjCoefficients(OsiClpSolverInterface& osi);
double getObjValue(CbcModel& c);
double getObjValue(OsiClpSolverInterface& osi);
double getObjValue(OsiSolverInterface& c);
double getObjSense(OsiClpSolverInterface& osi);
double getRowLower(CbcModel& c, int rowIndex);
double getRowLower(OsiClpSolverInterface& osi, int rowIndex);
std::string getRowName(OsiClpSolverInterface& osi, int rowIndex, unsigned maxLen=std::string::npos);
double getRowUpper(CbcModel& c, int rowIndex);
double getRowUpper(OsiClpSolverInterface& osi, int rowIndex);
void initialSolve(CbcModel& c);
void initialSolve(OsiClpSolverInterface& solver);
bool isBinary(CbcModel& c, int colIndex);
bool isBinary(OsiClpSolverInterface& osi, int colIndex);
bool isContinuous(CbcModel& c, int colIndex);
bool isContinuous(OsiClpSolverInterface& osi, int colNumber);
bool isInteger(CbcModel& c, int colIndex);
bool isInteger(OsiClpSolverInterface& osi, int colIndex);
bool isProvenOptimal(CbcModel& c);
bool isProvenOptimal(OsiClpSolverInterface& osi);
bool isProvenInfeasible(CbcModel& c);
CbcModel* newCbcModel(OsiClpSolverInterface& s);
OsiClpSolverInterface* newOsiClpSolverInterface();
void modifyCoefficient(OsiClpSolverInterface& osi, int row, int column, double newElement, bool keepZero=false);
int readLp(OsiClpSolverInterface& osi, const char *filename, const double epsilon=1e-5);
void resolve(OsiSolverInterface& osi);
void setCoefficients(OsiClpSolverInterface& osi, int row, int numberElements, const int *cols, const double *coefs);
void setColBounds(OsiClpSolverInterface& osi, int elementIndex, double lower, double upper);
void setColLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
void setColName(OsiClpSolverInterface& osi, int colIndex, std::string name);
void setColUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
void setContinuous(OsiClpSolverInterface& osi, int index);
void setInteger(OsiClpSolverInterface& osi, int index);
void setLogLevel(CbcModel& c, int n);
void setLogLevel(OsiClpSolverInterface& s, int n);
void setObjective(OsiClpSolverInterface& osi, const double *array);
void setObjCoeff(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
void setObjCoefficients(OsiClpSolverInterface& osi, int numberElements, const int *rows, const double *coefs);
void setObjSense(CbcModel& c, double s);
void setObjSense(OsiClpSolverInterface& osi, double s); 
void setRowLower(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
void setRowName(OsiClpSolverInterface& osi, int rowIndex, std::string name);
void setRowUpper(OsiClpSolverInterface& osi, int elementIndex, double elementValue);
OsiClpSolverInterface* solver(CbcModel& m);
int status(CbcModel& c);
void writeLp(OsiClpSolverInterface& osi, const char *filename, const char *extension="lp", double epsilon=1e-5, int numberAcross=10, int decimals=5, double objSense=0.0, bool useRowNames=true);
#endif
