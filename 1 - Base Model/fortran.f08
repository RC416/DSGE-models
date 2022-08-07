! Implementation of Base Model in Fortran (2008).

program main
    implicit none

    ! ---------------------------------------------------------------------------------------------------
    ! Declaration Section.

    ! Double precision variables.
    real(8), parameter  :: alpha    = 0.400
    real(8), parameter  :: beta     = 0.987
    real(8), parameter  :: delta    = 1.000
    real(8) :: k_steady, k_pct_low, k_pct_high
    real(8) :: v_max, kt1_optimal, new_value_function_value
    
    ! Integer variables.
    integer, parameter  :: number_of_iterations = 1000
    integer, parameter  :: number_of_k_values   = 201
    integer :: i
    integer :: iteration, kt0_index, kt1_index

    ! Value Function, Policy Function, and other arrays.
    real(8), dimension(number_of_iterations, number_of_k_values) :: Value_Function, Policy_Function
    real(8), dimension(number_of_k_values)   :: k_values

    ! ---------------------------------------------------------------------------------------------------
    ! Execution Section.

    ! Calculate the steady-state level of capital.
    k_steady = ((1-beta*(1-delta))/(alpha*beta*1)) ** (1/(alpha-1))

    ! Create grid of capital values around steady-state (+/- 50%).
    k_pct_low = 0.50
    k_pct_high = 1.50

    do i = 1, number_of_k_values
        k_values(i) = k_steady * (k_pct_low + ((real(i,8) - 1) / real(number_of_k_values,8) * (k_pct_high - k_pct_low)))
    end do

    ! Assign value of 0 to first value function iteration.
    where (Value_Function /= 0.0) Value_Function = 0
    where (Policy_Function /= 0.0) Policy_Function = 0

    ! Perform value function iteration.
    do iteration = 2, number_of_iterations      
        do kt0_index = 1, number_of_k_values        ! for each level of starting capital...
            
            ! Initialize variables to store candiate optimal values.
            v_max = -huge(k_steady)
            kt1_optimal = 0.0
            new_value_function_value = 0.0

            do kt1_index = 1, number_of_k_values    ! ... check all next period capital choices

                ! Calculate Value Function for given starting capital and next period capital choice.
                new_value_function_value = &
                    log(k_values(kt0_index)**alpha + (1-delta)*k_values(kt0_index) - k_values(kt1_index)) &
                    + beta*Value_Function(iteration - 1, kt1_index)

                ! Check if this capital choice gives highest Value Function value.
                if (new_value_function_value > v_max) then

                    ! Update candidate values.
                    v_max = new_value_function_value
                    kt1_optimal = k_values(kt1_index)
                end if
            end do
            
            ! Update Value Function and Policy Function with optimal values.
            Value_Function(iteration, kt0_index) = v_max
            Policy_Function(iteration, kt0_index) = kt1_optimal

        end do
    end do

    ! Write Value Function and Policy Function to csv files.
    call write_to_csv_file(Value_Function, number_of_iterations, number_of_k_values, "Value_Function.csv")
    call write_to_csv_file(Policy_Function, number_of_iterations, number_of_k_values, "Policy_Function.csv")

    ! Display the first and last 5 values of the Value Function and Policy Function.
    write(*,*) "Value Function V(k)"
    do i = 1,5
        print '(2 (A, G0.4) )', "V(", k_values(i), ") = ", Value_Function(number_of_iterations, i)
    end do
    write(*,*) "..."
    do i = number_of_k_values-4,number_of_k_values
        print '(2 (A, G0.4) )', "V(", k_values(i), ") = ", Value_Function(number_of_iterations, i)
    end do

    write(*,*)
    write(*,*) "Policy Function g(k)"
    do i = 1,5
        print '(2 (A, G0.4) )', "g(", k_values(i), ") = ", Policy_Function(number_of_iterations, i)
    end do
    write(*,*) "..."
    do i = number_of_k_values-4,number_of_k_values
        print '(2 (A, G0.4) )', "g(", k_values(i), ") = ", Policy_Function(number_of_iterations, i)
    end do

    ! Leave window open after program terminates.
    call sleep(100)

end program


! ---------------------------------------------------------------------------------------------------
! Subroutines.
! 1 - Write 2-dimensional array to .csv file.

! 1 - Writes a 2-dimensional array of real(8) numbers (double floating point precision) to a csv file.
subroutine write_to_csv_file(Array2D, n_rows, n_cols, file_name)

    ! -----------------------------------------------------------------------------------------------
    ! Declaration section.

    ! Input variables.
    integer, intent(in) :: n_rows, n_cols
    character(len=*), intent(in) :: file_name
    real(8), dimension(n_rows, n_cols), intent(in) :: Array2D

    ! Local variables. 
    integer :: i

    ! -----------------------------------------------------------------------------------------------
    ! Execution section.
    
    ! Formatting for CSV (label = 101)
    101 format(1x, *(g0, ", "))

    ! Create and open file.
    open(unit = 10, access = "sequential", action = "write", &
         status = "replace", file = file_name, form = "formatted")

    ! Write each row to file.
    do i = 1, n_rows
        write(10, 101) Array2D(i,:)
    end do

    ! Close file.
    close(10)
 
end subroutine write_to_csv_file