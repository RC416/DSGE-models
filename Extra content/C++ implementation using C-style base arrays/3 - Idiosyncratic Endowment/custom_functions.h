#pragma once
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <vector>

// Data structure to store parameters.
struct Parameters
{
	double sigma, beta;
	const int number_of_a_values, number_of_e_values;
	double* a_grid;
	double* e_grid;
	double* e_probs;
};

// Functions to solve the household's problem, Value Function, and Population Distribution.
void Solve_HH_Problem(double q, int a_start_index, int e_start_index,
	double* Value_Function, Parameters& params,
	double& V, double& g, int& g_index);

void Solve_Value_Function(double q, Parameters& params,
	double* Value_Function, double* Policy_Function, int* Policy_Function_Index,
	double* Value_Function_New, double* Policy_Function_New, int* Policy_Function_Index_New);

void Get_Population_Distribution(int* Policy_Function_Index, Parameters& params,
	double* Population_Distribution, double* Population_Distribution_New);

double Calculate_Market_Clearing_Condition(
	double* Population_Distribution,
	double* Policy_Function,
	int n_rows, int n_cols);

double GetDistanceArray2D(
	double* Array1,
	double* Array2,
	int n_rows, int n_cols);

// Functions to copy and write arrays.
template<class number_type>
void DeepCopyArray2D(
	number_type* Array_Original,
	number_type* Array_Copy,
	int n_rows, int n_cols);

template<class number_type> void WriteArrayToCSV(number_type* Array, const char* file_name, int n_rows, int n_cols);
template<class number_type> void WriteArrayToCSV(number_type* Array, const char* file_name, int n_rows, int n_cols);