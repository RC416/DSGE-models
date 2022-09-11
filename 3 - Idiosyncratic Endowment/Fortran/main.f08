! Idiosyncratic Endowment model implemented in Fortran (2008).
! Steps:
!    1 - Define utility parameters, grids, and parameter struct.
!    2 - Solve model and get market bond price using binary search.
!        - Guess bond price
!        - Solve model and get market bond price using binary search.
!        - Get the distribution of credit and productivity levels
!        - Check the market clearing condition and update bond price guess
!    3 - Print results and save to file.

program main
    
    use custom_functions     
    implicit none

    ! ---------------------------------------------------------------------------------------------------
    ! Declaration section.

    ! Endowment parameters.
    real(8), parameter :: e_high = 1.0
    real(8), parameter :: e_low = 0.1
    integer, parameter :: number_of_e_values = 2
    real(8), dimension(number_of_e_values) :: e_grid = [e_low, e_high]
    real(8), dimension(number_of_e_values, number_of_e_values) :: e_probs = &
        reshape([0.500, 0.075, 0.500, 0.925], [number_of_e_values, number_of_e_values])
    
    ! Utility parameters.
    real(8) :: sigma = 1.5
    real(8) :: beta = 0.99322

    ! Credit parameters.
    real(8) :: a_high = 4
    real(8) :: a_low = -2
    integer, parameter :: number_of_a_values = 100
    real(8), dimension(number_of_a_values) :: a_grid

    ! Bonds prices and iteration parameters.
    real(8) :: q, q_min, q_max, mcc
    real(8) :: dist, tolerance
    integer :: iteration_count, max_iterations, a_index, e_index, i

    ! Value Function, Policy Functions, and Population Distribution.
    real(8), dimension(number_of_e_values, number_of_a_values) :: Value_Function = 0
    real(8), dimension(number_of_e_values, number_of_a_values) :: Policy_Function, &
    Population_Distribution, Value_Function_New, Policy_Function_New, Population_Distribution_New
    integer, dimension(number_of_e_values, number_of_a_values) :: Policy_Function_Index, &
    Policy_Function_Index_New

    ! Struct to store utility parameters and capital/productivity grids for passing to a function.   
    type(Parameters) :: params

    ! ---------------------------------------------------------------------------------------------------
    ! Execution section.

    ! Create grid of credit values.
    do i = 1, number_of_a_values
        a_grid(i) = a_low + (a_high - a_low) * real(i - 1, 8) / real(number_of_a_values - 1, 8)
    end do

    ! Store parameters in a struct for passing to a function.
    params = Parameters(sigma, beta, number_of_a_values, number_of_e_values, a_grid, e_grid, e_probs)

    ! Range of bond values to search.
    q_min = 0.985
    q_max = 1.100

    ! Optional: floor for q_min such that those at credit limit can still afford positive consumption.
    q_min = (a_low + e_low) / a_low

    ! Placeholder for the market clearing condition.
    mcc = huge(mcc);

    ! Iteration parameters.
    dist = huge(dist);
    iteration_count = 0;
    max_iterations = 20;
    tolerance = 1e-3;

    ! Solve for market price q.
    do while ((dist > tolerance) .and. (iteration_count < max_iterations))

        ! Get value of q from middle of range.
        q = (q_min + q_max) / 2

        ! Solve for the Value Function and Policy Function.
        call Solve_Value_Function(q, params, Value_Function, Policy_Function, Policy_Function_Index, &
        Value_Function_New, Policy_Function_New, Policy_Function_Index_New)

        ! Get the Population Distribution.
        call Get_Population_Distribution(Policy_Function_Index, params, & 
        Population_Distribution, Population_Distribution_New)

        ! Check the market clearing condition.
        call Calculate_Market_Clearing_Condition(mcc, Population_Distribution, Policy_Function)
        
        ! Update search parameters.
        dist = abs(mcc)
        iteration_count = iteration_count + 1

        ! Update range of q according to the sign of the market clearing condition.
        if (mcc > 0) q_min = q
        if (mcc < 0) q_max = q

        ! Print results.
        write (*,*) "Iteration", iteration_count, "   ", "q=", q, "mcc=", mcc;
    end do

    ! Write final Value Function, Policy Function, and Population Distribution to csv files.
    call WriteArrayToCSV(transpose(Value_Function), number_of_a_values, number_of_e_values, "Value_Function.csv")
    call WriteArrayToCSV(transpose(Policy_Function), number_of_a_values, number_of_e_values, "Policy_Function.csv")
    call WriteArrayToCSV(transpose(Population_Distribution), number_of_a_values, number_of_e_values, "Population_Distribution.csv")

    ! Display subset of results: Value Function for certain capital and productivity terms.
    write(*,*) NEW_LINE('A'), "Value Function", NEW_LINE('A')
    do a_index = 1, 10
        do e_index = 1, 2
            write(*, '(4x, f8.3)', advance='no') Value_Function(e_index, a_index)
        end do
        write(*,*) NEW_LINE('A')
    end do
    call sleep(100)

end program

! Compilier instructions.
! debug: gfortran src/custom_functions.f08 src/main.f08 -g -O0 -Wall -o build.exe
! release: gfortran src/custom_functions.f08 src/main.f08 -O3 -o build.exe

! Debug and profiling:
! gfortran -g -O0 -Wall -pg -o obj/custom_functions.o -c src/custom_functions.f08
! gfortran -g -O0 -Wall -pg -o obj/main.o -c src/main.f08
! gfortran -g -O0 -Wall -pg -o build obj/custom_functions.o obj/main.o 
! run program: start build
! profile: gprof build.exe

! Release:
! gfortran -O3 -o obj/custom_functions.o -c src/custom_functions.f08
! gfortran -O3 -o obj/main.o -c src/main.f08
! gfortran -O3 -o build obj/custom_functions.o obj/main.o 
! rm obj/custom_functions.o obj/main.o custom_functions.d main.d build