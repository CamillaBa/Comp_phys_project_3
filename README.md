# Comp_phys_project_4

There are 3 code files.

- Main.cpp
- PythonApplication1.py
- module1.py

The file "Main.cpp" contains the main C++ program,
and all other functions that we have used in our code.
That includes the implementation of the metropolis algorithm.
The code snippets that control the various functions and produce
outputs are separated by rows of comments, for example

//=========================================================
// b)
//=========================================================

The code snippets right under answer some part of the project,
in the sense that they generate the corresponding data needed to make
a plot.

When the code is run as it stands, it will print:

"
Energy function works as intended.
Magnetization function works as intended.
completed cycle: 10000
completed cycle: 20000
completed cycle: 30000
completed cycle: 40000
completed cycle: 50000
completed cycle: 60000
completed cycle: 70000
completed cycle: 80000
completed cycle: 90000
completed cycle: 100000
Success!"

A file is created labeled "small_lattice_2_2__cycles_100000_T_1.000000.txt".
This contains CSV formated data for the expectation values of the quanteties
asked about in exercise b).


Similar files are created for all other such code snippets.





As for PythonApplication1.py, it is formatted similarly with respect to the exercises,
but it only creates plots based on the data that each exercise requires.
When it is run as is (after plot data for b) is generated), it will print the plot for exercise b).
