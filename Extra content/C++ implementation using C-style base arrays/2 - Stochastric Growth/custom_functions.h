#pragma once
#include<string>
#include<vector>

// Data structure to store parameters.
struct Parameters
{
	double alpha, beta, delta;
	const int number_of_k_values, number_of_z_values;
	double* k_values;
	double* z_values;
	double** z_probs;
};

// Function to solve Household's problem for given starting states.
void Solve_HH_Problem(double*** Value_Function, double*** Policy_Function, int iteration, int kt0_index, int zt_index, Parameters params);

// Functions to initialize, deallocate, read and write arrays.
double*** InitializeArray3D(int size_dim_1, int size_dim_2, int size_dim_3);
double** InitializeArray2D(int n_rows, int n_cols);
void DeleteArray3D(double*** Array3D, int size_dim_1, int size_dim_2, int size_dim_3);
void DeleteArray2D(double** Array2D, int n_rows, int n_cols);
void WriteArrayToCSV(double** Array2D, int n_rows, int n_cols, const char* file_name);
double** ReadArrayFromCSV(std::string path);
double* ReadVectorFromCSV(std::string path);