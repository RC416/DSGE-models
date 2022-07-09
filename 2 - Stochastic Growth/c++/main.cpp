/* 
Stochastic Growth Model implemented in c++.

Steps:
	1 - Define utility parameters, grids, and parameter struct.
	2 - Perform value function iteration.
	3 - Display results, save results to file, clean up.
	
Relies on Iterate_Value_Function and helper functions.
*/

#include "iteration_functions.h"
#include <iostream>
#include <string>
#include <vector>

int main()
{
	// Assign parameter values.
	double alpha = 0.400;
	double beta  = 0.987;
	double delta = 0.012;
	const int number_of_iterations = 1000;

	// Calculate the steady-state level of capital.
	double k_steady = pow(((1 - beta * (1 - delta)) / (alpha * beta)), (1 / (alpha - 1)));

	// Create a grid of capital values around steady-state (+/- 2%).
	const int number_of_k_values = 201;
	double k_low_pct  = 0.98;
	double k_high_pct = 1.02;
	double k_values[number_of_k_values];

	for (int i = 0; i < number_of_k_values; i++)
	{
		k_values[i] = k_low_pct * k_steady + (double(i) / (double(number_of_k_values) - 1)) * ((k_high_pct - k_low_pct) * k_steady);
	}

	// Get productivity levels and transition probabilities from csv files (as vector arrays).
	std::string directory = "C:/Users/Ray/Documents/GitHub/DSGE-models/2 - Stochastic Growth/c++/Inputs/";
	std::vector<std::vector<double>> z_probs_csv = ReadArrayFromCSV(directory + "z_probs.csv");
	std::vector<std::vector<double>> z_values_csv = ReadArrayFromCSV(directory + "z_values.csv");

	// Convert std::vector arrays to base arrays.
	const int number_of_z_values = static_cast<int>(z_values_csv.size());	// = 11, cast from unsigned int to signed int
	double** z_probs = InitializeArray2D(number_of_z_values, number_of_z_values);
	double* z_values = new double[number_of_z_values];

	for (int z1_index = 0; z1_index < number_of_z_values; z1_index++)
	{
		for (int z2_index = 0; z2_index < number_of_z_values; z2_index++)
		{
			z_probs[z1_index][z2_index] = z_probs_csv[z1_index][z2_index];	// copy z_probs to new array
		}
		z_values[z1_index] = z_values_csv[z1_index][0];				// copy z_values to new vector from column 0 of vector array
	}

	// Initialize Value Function and Policy Function (as arrays).
	double*** Value_Function  = InitializeArray3D(number_of_iterations, number_of_k_values, number_of_z_values);
	double*** Policy_Function = InitializeArray3D(number_of_iterations, number_of_k_values, number_of_z_values);
	
	// Assign value of 0 to first iteration.
	for (int i = 0; i < number_of_k_values; i++)
	{
		for (int j = 0; j < number_of_z_values; j++)
		{
			Value_Function[0][i][j] = 0.0;
			Policy_Function[0][i][j] = 0.0;
		}
	}

	// Store utility parameters and capital/productivity grids in a struct for passing to a function.
	Parameters params = { alpha, beta, delta, number_of_k_values, number_of_z_values, k_values, z_values, z_probs };

	// Perform value function iteration.
	for (int iteration = 1; iteration < number_of_iterations; iteration++)
	{
		// Loop over all possible starting states.
		for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
		{
			for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
			{
				// Solve Value Function and Policy Function.
				Iterate_Value_Function(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params);
			}
		}
	}

	// Write final Value Function and Policy Function to csv files.
	double** Final_Value_Function  = InitializeArray2D(number_of_k_values, number_of_z_values);
	double** Final_Policy_Function = InitializeArray2D(number_of_k_values, number_of_z_values);

	for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
	{
		for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
		{
			Final_Value_Function[kt0_index][zt_index] = Value_Function[number_of_iterations - 1][kt0_index][zt_index];
			Final_Policy_Function[kt0_index][zt_index] = Policy_Function[number_of_iterations - 1][kt0_index][zt_index];
		}
	}

	WriteArrayToCSV(Final_Value_Function, number_of_k_values, number_of_z_values, "Value_Function.csv");
	WriteArrayToCSV(Final_Policy_Function, number_of_k_values, number_of_z_values, "Policy_Function.csv");

	// Display subset of results: final Value Function for certain capital and prodcutivity values.
	for (int kt0_index = 0; kt0_index < 10; kt0_index++)
	{
		for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
		{
			std::cout << Value_Function[number_of_iterations-1][kt0_index][zt_index] << " ";
		}
		std::cout << "\n";
	}

	// Remove arrays from memory.
	DeleteArray3D(Value_Function, number_of_iterations, number_of_k_values, number_of_z_values);
	DeleteArray3D(Policy_Function, number_of_iterations, number_of_k_values, number_of_z_values);
	DeleteArray2D(Final_Value_Function, number_of_k_values, number_of_z_values);
	DeleteArray2D(Final_Policy_Function, number_of_k_values, number_of_z_values);
	std::cin.get();
	return 0;
}
