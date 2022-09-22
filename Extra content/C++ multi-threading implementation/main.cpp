/*
Stochastic Growth Model implemented in C++.

This is the basic multi-threaded implementation using only the standard library.
The changes for this multi-threaded implementations are on lines 76-120.

Steps:
	1 - Define utility parameters, grids, and parameter struct.
	2 - Perform value function iteration.
	3 - Display results, save results to file, clean up.

Notes on multi-threaded implementation
	- Introduces a lambda function which applies the Solve_HH_Problem over all capital values 
	for a fixed productivity value (i.e. solves one column of the value function).
	- The lambda function is run in parallel threads for all productivity values (i.e. the columns
	of the value function are solved in parallel).
	- Threads are stored in a list and must be terminated before passing to the next iteration.
	- Threads are created and destroyed in each iteration. This avoids the need for additional code
	or packages to manage threads, but introduces additional overhead.
	- Net improvement is ~4x (uses up to 85% of 6-core 12-thread CPU)
*/

#include "custom_functions.h"
#include <iostream>
#include <string>
#include <vector>
#include <thread>
using std::thread;

int main()
{
	for (int m = 0; m < 10; m++) {
		// Assign parameter values.
		double alpha = 0.400;
		double beta = 0.987;
		double delta = 0.012;
		int const number_of_iterations = 1000;

		// Calculate the steady-state level of capital.
		double k_steady = pow(((1 - beta * (1 - delta)) / (alpha * beta)), (1 / (alpha - 1)));

		// Create a grid of capital values around the steady-state (+/- 2%).
		int const number_of_k_values = 201;
		double k_low_pct = 0.98;
		double k_high_pct = 1.02;
		vector<double> k_values(number_of_k_values);

		for (int i = 0; i < number_of_k_values; i++)
		{
			k_values[i] = k_low_pct * k_steady + (double(i) / (double(number_of_k_values) - 1)) * ((k_high_pct - k_low_pct) * k_steady);
		}

		// Get productivity levels and transition probabilities from csv files (as vector arrays).
		std::string directory = "C:/Users/Ray/Documents/GitHub/DSGE-models/2 - Stochastic Growth/c++/Inputs/";
		Array2D<double> z_probs = ReadArrayFromCSV(directory + "z_probs.csv");
		vector<double> z_values = ReadVectorFromCSV(directory + "z_values.csv");
		int number_of_z_values = z_values.size();

		// Initialize the Value Function and Policy Function (as arrays).
		Array3D<double> Value_Function(number_of_iterations, vector<vector<double>>(number_of_k_values, vector<double>(number_of_z_values)));
		Array3D<double> Policy_Function(number_of_iterations, vector<vector<double>>(number_of_k_values, vector<double>(number_of_z_values)));

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
			/* // Original loop (inner 2 loops) that gets converted into multi-threaded implementation.
			for (int iteration = 1; iteration < number_of_iterations; iteration++)
			{
				// Loop over all possible starting states.
				for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
				{
					for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
					{
						// Solve the Value Function and Policy Function.
						Solve_HH_Problem(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params);
					}
				}
			} */

			// Lambda function which takes a zt_index and solves the HH problem for all k_values at this zt_index.
			// The function is applied in parallel over the 11 zt_index values.
			// See documentation: https://learn.microsoft.com/en-us/cpp/cpp/lambda-expressions-in-cpp?view=msvc-170
			auto Solve_HH_Lambda = [&, iteration](int zt_index) {

				// Loop over all kt0_index starting states.
				for (int kt0_index = 0; kt0_index < number_of_k_values; kt0_index++)
				{
					// Solve the Value Function and Policy Function.
					Solve_HH_Problem(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params);
				}
			};

			// Vector to store thread objects.
			vector<thread> Threads;

			// Loop over all possible zt_index starting states.
			for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
			{
				// Create a thread for this zt_index and store in vector of threads.
				thread thread_object = thread(Solve_HH_Lambda, zt_index);
				Threads.push_back(std::move(thread_object));

				// Alternatively: 
				//Threads.push_back(thread(Solve_HH_Lambda, zt_index));
			}

			// Wait for all threads to terminate.
			for (thread& thread_object : Threads)
			{
				// .join() method forces the program to wait for each thread to finish.
				if (thread_object.joinable()) { thread_object.join(); }
			}
		}

		// Get slice of Value Function and Policy Function for final iteration.
		Array2D<double> Final_Value_Function(number_of_k_values, vector<double>(number_of_z_values));
		Array2D<double> Final_Policy_Function(number_of_k_values, vector<double>(number_of_z_values));

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

		// Display a subset of results: the final Value Function for certain capital and prodcutivity values.
		for (int kt0_index = 0; kt0_index < 10; kt0_index++)
		{
			for (int zt_index = 0; zt_index < number_of_z_values; zt_index++)
			{
				std::cout << Value_Function[number_of_iterations - 1][kt0_index][zt_index] << "\t";
			}
			std::cout << "\n";
		}


	}
	std::cin.get();
	return 0;
}