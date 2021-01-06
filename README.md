# DSGE-models
Implementation of a simple DSGE model with Python, Julia and Matlab

I wanted to compare Python and Julia for solving dynamic macroeconomic models. I wanted to know:
- Is Julia meaningfully faster than Python for this task?
- Is Julia a suitable alternative to something like c++ for writing fast code?
- Are there other notable advantages/disadvantages to Julia?

Unfortunately, Julia did not shine in this application. It was not faster than Python and was not otherwise better. I liked some of the available features and plan to try it in other projects.

The task of solving this model breaks down to three steps:
1. Iterating over many combinations of variables
2. Calculating a scalar value for each combination of variables.
3. Appending that scalar to an n-dimensional array.

This would be an excellent candidate for parallel computing if not for the dynamic nature of the problem; iteration 2 depends on the results from iteration 1, so iterations cannot be done independently and in parallel.

In my testing, I found that Julia was not faster than Python (slower if anything). I expected Julia to be faster since I was using multiple for loops in 1.

Steps 2. and 3. used the math and numpy packages. I hypothesize that steps 2. and 3. are the bottleneck and using the precompiled packages let Python be just as fast as Julia. 

# Overall pros/cons of the languages for this applicaton.

## Python:
Pros
- Great support since the language is mature with large user base.
- Excellent selection of packages for nearly any task.

Cons
- Requires coding in a second language (such as c or c++) or a precompiled package to acheive top speeds.
- Old langauge which may carry some poor design choices.


## Julia
Pros
- Simple syntax which is similar to Python.
- Speed approaching the fast languages (c, c++, fortran).
- Some nice frills like unicode characters and having useful constants/functions loaded automatically.
- Can use Python packages, including plotting.

Cons
- VSCode is the only supported IDE (as of Jan 2021). My initial imporession is that it is less user frendly than R Studio and Spyder (for R and Python, respectively) for doing work with data.
- Less support as user base is small and language is still new (version 1.5 as of Jan 2021).
