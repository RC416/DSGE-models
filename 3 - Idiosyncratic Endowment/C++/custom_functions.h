#pragma once
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <vector>

// Containers for dynamic arrays.
using std::vector;
template<class number_type> using Array2D = vector<vector<number_type>>;

// Data structure to store parameters.
struct Parameters
{
	double sigma, beta;
	const int number_of_a_values, number_of_e_values;
	vector<double> a_grid;
	vector<double> e_grid;
	Array2D<double> e_probs;
};

// Functions to solve the household's problem, Value Function, and Population Distribution.
void Solve_HH_Problem(double q, int a_start_index, int e_start_index, Array2D<double>& Value_Function,
	Parameters& params, double& V, double& g, int& g_index);
void Solve_Value_Function(double q, Parameters& params,
	Array2D<double>& Value_Function, Array2D<double>& Policy_Function, Array2D<int>& Policy_Function_Index);
void Get_Population_Distribution(Array2D<int>& Policy_Function_Index, Parameters& params,
	Array2D<double>& Population_Distribution);
double Calculate_Market_Clearing_Condition(Array2D<double>& Population_Distribution, Array2D<double>& Policy_Function);

// Functions to get distance, copy, and write arrays.
double GetDistanceArray2D(Array2D<double>& Array1, Array2D<double>& Array2);
template<class number_type> void DeepCopyArray2D(Array2D<number_type>& Array_Original, Array2D<number_type>& Array_Copy);
template<class number_type> void WriteArrayToCSV(Array2D<number_type>& Array, const char* file_name);