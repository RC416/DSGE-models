! Stochastic Growth Model implemented in Fortran (2008).

! Steps:
! 	1 - Define utility parameters, grids, and parameter struct.
!	2 - Perform value function iteration.
!	3 - Display results, save results to file, clean up.
	
program main
    
    use custom_functions     
    implicit none

    ! ---------------------------------------------------------------------------------------------------
    ! Declaration section.

    ! Double precision parameter values.
    real(8), parameter :: alpha = 0.400
    real(8), parameter :: beta  = 0.987
    real(8), parameter :: delta = 0.012
    real(8) :: k_steady, k_pct_low, k_pct_high

    ! Integer variables.
    integer, parameter :: number_of_iterations = 1000
    integer, parameter :: number_of_k_values = 201
    integer, parameter :: number_of_z_values = 11
    integer ::  i, iteration, kt0_index, zt_index

    ! Value Function, Policy Function, and other arrays.
    real(8), allocatable :: Value_Function(:,:,:), Policy_Function(:,:,:)
    real(8), dimension(number_of_k_values) :: k_values
    real(8), dimension(number_of_z_values) :: z_values
    real(8), dimension(number_of_z_values, number_of_z_values) :: z_probs

    ! Struct to store utility parameters and capital/productivity grids for passing to a function.   
    type(Parameters) :: params

    ! Allocate memory for allocatable arrays.
    allocate(Value_Function(number_of_z_values, number_of_k_values, number_of_iterations))
    allocate(Policy_Function(number_of_z_values, number_of_k_values, number_of_iterations))

    ! ---------------------------------------------------------------------------------------------------
    ! Execution section.

    ! Calculate the steady-state level of capital.
    k_steady = ((1 - beta * (1 - delta)) / (alpha * beta)) ** (1 / (alpha - 1))

    ! Create a grid of capital values around the steady-state (+/- 2%).
    k_pct_low = 0.98
    k_pct_high = 1.02

    do i = 1, number_of_k_values
        k_values(i) = k_steady * (k_pct_low + ((real(i, 8) - 1) / real(number_of_k_values, 8) * (k_pct_high - k_pct_low)))
    end do
    
    ! Assign value of 0 to first value function iteration.
    where (Value_Function /= 0.0) Value_Function = 0
    where (Policy_Function /= 0.0) Policy_Function = 0

    ! Get productivity levels and transition probabilities from csv files.
    call ReadVectorFromCSV("Inputs\z_values.csv", number_of_z_values, z_values)
    call ReadArrayFromCSV("Inputs\z_probs.csv", number_of_z_values, number_of_z_values, z_probs)

    ! Store utility parameters and capital/productivity grids in a struct for passing to a function.
    params = Parameters(alpha, beta, delta, number_of_k_values, number_of_z_values, k_values, z_values, z_probs)

    ! Perform value function iteration.
    do iteration = 2, number_of_iterations
        
        ! Loop over all possible starting states.
        do kt0_index = 1, number_of_k_values
            do zt_index = 1, number_of_z_values

                ! Solve the Value Function and Policy Function.
                call Solve_HH_Problem(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params)
            end do 
        end do
    end do

    ! Write the final Value Function and Policy Function to csv files.
    call WriteArrayToCSV(transpose(Value_Function(:,:,number_of_iterations)), & 
    number_of_k_values, number_of_z_values, "Value_Function.csv")
    call WriteArrayToCSV(transpose(Policy_Function(:,:,number_of_iterations)), &
    number_of_k_values, number_of_z_values, "Policy_Function.csv")

    ! Display a subset of results: the final Value Function for certain capital and productivity values.
    do kt0_index = 1, 10
        do zt_index = 1, 10
            write(*, '(4x, f7.3)', advance='no') Value_Function(zt_index, kt0_index, number_of_iterations)
        end do
        write(*,*) NEW_LINE('A')
    end do
    call sleep(100)

end program

! Compilier instructions.
! debug mode: gfortran src/custom_functions.f08 src/main.f08 -g -O0 -Wall -o start build.exe
! release mode: gfortran src/custom_functions.f08 src/main.f08 -O3 -o start build.exe
! run program: start build