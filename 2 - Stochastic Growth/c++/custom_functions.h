#pragma once
#include<string>
#include<vector>

// Define vector-based type aliases for arrays.
using std::vector;
template<class number_type>
using Array2D = vector<vector<number_type>>;
template<class number_type>
using Array3D = vector<vector<vector<number_type>>>;

// Data structure to store parameters.
struct Parameters
{
	double alpha, beta, delta;
	const int number_of_k_values, number_of_z_values;
	vector<double> k_values;
	vector<double> z_values;
	Array2D<double> z_probs;
};

// Function to solve the household's problem for a given starting state.
void Solve_HH_Problem(Array3D<double>& Value_Function, Array3D<double>& Policy_Function, 
	int iteration, int kt0_index, int zt_index, Parameters& params);

// Functions to read and write arrays.
void WriteArrayToCSV(Array2D<double> Array, int n_rows, int n_cols, const char* file_name);
vector<double> ReadVectorFromCSV(std::string path);
Array2D<double> ReadArrayFromCSV(std::string path);