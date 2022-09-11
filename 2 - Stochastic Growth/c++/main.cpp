/*
Stochastic Growth Model implemented in c++.

Steps:
	1 - Define utility parameters, grids, and parameter struct.
	2 - Perform value function iteration.
	3 - Display results, save results to file, clean up.

Relies on Solve_Value_Function and helper functions in custom_funcions file.
*/

#include "custom_functions.h"
#include <iostream>
#include <string>
#include <vector>

int main()
{
	// Assign parameter values.
	double alpha = 0.400;
	double beta = 0.987;
	double delta = 0.012;
	int const number_of_iterations = 1000;

	// Calculate the steady-state level of capital.
	double k_steady = pow(((1 - beta * (1 - delta)) / (alpha * beta)), (1 / (alpha - 1)));

	// Create a grid of capital values around steady-state (+/- 2%).
	int const number_of_k_values = 201;
	double k_low_pct = 0.98;
	double k_high_pct = 1.02;
	double k_values[number_of_k_values];

	for (int i = 0; i < number_of_k_values; i++)
	{
		k_values[i] = k_low_pct * k_steady + (double(i) / (double(number_of_k_values) - 1)) * ((k_high_pct - k_low_pct) * k_steady);
	}

	// Get productivity levels and transition probabilities from csv files (as vector arrays).
	std::string directory = "C:/Users/Ray/Documents/GitHub/DSGE-models/2 - Stochastic Growth/c++/Inputs/";
	double* z_values = ReadVectorFromCSV(directory + "z_values.csv");
	double** z_probs = ReadArrayFromCSV(directory + "z_probs.csv");
	int const number_of_z_values = 11;

	// Initialize Value Function and Policy Function (as arrays).
	double*** Value_Function = InitializeArray3D(number_of_iterations, number_of_k_values, number_of_z_values);
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
		for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
		{
			for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
			{
				// Solve Value Function and Policy Function.
				Solve_HH_Problem(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params);
			}
		}
	}

	// Get slice of Value Function and Policy Function for final iteration.
	double** Final_Value_Function = InitializeArray2D(number_of_k_values, number_of_z_values);
	double** Final_Policy_Function = InitializeArray2D(number_of_k_values, number_of_z_values);

	for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
	{
		for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
		{
			Final_Value_Function[kt0_index][zt_index] = Value_Function[number_of_iterations - 1][kt0_index][zt_index];
			Final_Policy_Function[kt0_index][zt_index] = Policy_Function[number_of_iterations - 1][kt0_index][zt_index];
		}
	}

	// Write final Value Function and Policy Function to csv files.
	WriteArrayToCSV(Final_Value_Function, number_of_k_values, number_of_z_values, "Value_Function.csv");
	WriteArrayToCSV(Final_Policy_Function, number_of_k_values, number_of_z_values, "Policy_Function.csv");

	// Display subset of results: final Value Function for certain capital and prodcutivity values.
	for (int kt0_index = 0; kt0_index < 10; kt0_index++)
	{
		for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
		{
			std::cout << Value_Function[number_of_iterations - 1][kt0_index][zt_index] << "\t";
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