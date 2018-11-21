#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <ctime>
#include <string>
#include <omp.h>

//=============================================================================================
// Random numbers
//=============================================================================================

int rand_int(int min, int max) {
	/*Function to quickly get random integer between from {min,..,max-1}.
	Using uniform distribution from <random>.
	*/
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());


	std::uniform_int_distribution<int>  distr(min, max);
	return distr(generator);
}

double rand_double(double min, double max) {
	/*Function to quickly get random double between from {min,..,max-1}.
	Using uniform distribution from <random>.
	*/
	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());

	std::uniform_real_distribution<>    distr(min, max);
	return distr(generator);
}

//=============================================================================================
// simple vector class
//=============================================================================================

template <class T>
class vector {
private:
	// class variables
	T* entries;

public:
	int n;

	// constructor
	vector(int n) {
		entries = new T[n]{};
		vector::n = n;
	}

	// copy constructor
	vector(vector const & copy) {
		matrix::m = copy.m;

		// create new entries
		entries = new T[n]{};
		*entries = *copy.entries;
	}

	// destructor
	~vector() {
		delete[] entries;
	}

	// access i'th entry in vector
	T & operator[] (int i) {
		return entries[i];
	}
};

//=============================================================================================

template<class T>
void print_vector_to_file(vector<T> &vec, std::ofstream &file) {
	/* Function to print the contents of an instance of vector<T>
	as a line in the given file, formated as CSV.
	*/
	int n = vec.n; // get length of vector
	for (int i = 0; i < n - 1; i++) {
		file << vec[i] << ',';
	}
	file << vec[n - 1] << "\n";
}

//=============================================================================================
// simple matrix class
//=============================================================================================

template <class T>
class matrix {
private:
	// class variables
	T** entries;

public:
	// shape
	int m;
	int n;

	// constructor
	matrix(int m, int n) {
		entries = new T *[m];
		for (int i = 0; i < m; i++) {
			entries[i] = new T[n]();
		}
		matrix::m = m;
		matrix::n = n;
	}

	// copy constructor
	matrix(matrix const & copy) {
		matrix::m = copy.m;
		matrix::n = copy.n;

		// create new entries
		entries = new T *[m];
		for (int i = 0; i < m; i++) {
			entries[i] = new T[n]();
		}
		** entries = **copy.entries;
	}

	// access i'th row in matrix
	T* & operator[] (int i) {
		return entries[i];
	}

	// destructor
	~matrix() {
		for (int i = 0; i < matrix::m; i++) {
			delete[] entries[i];
		}
		delete[] entries;
	}
};

//=============================================================================================

template<class T>
void print(matrix<T> &A) {
	/* Function to print the entries of a given matrix<int> instance to terminal.
	*/

	int m = A.m; int n = A.n;
	std::cout << '(';
	for (int i = 0; i < m; i++) {
		std::cout << '(';
		for (int j = 0; j < n; j++) {
			if (j < n - 1) {
				std::cout << std::setw(4) << A[i][j] << ", ";
			}
			else {
				std::cout << std::setw(4) << A[i][j];
			}

		}
		if (i < m - 1) {
			std::cout << ')' << std::endl;
		}
		else {
			std::cout << ')';
		}

	}
	std::cout << ')' << std::endl;
}

template<class T>
void save_matrix_file(matrix<T>&s, std::string filename) {
	/* Simple function to save a matrix to a file.*/

	int m = s.m; int n = s.n;

	// open file
	std::ofstream myfile(filename);


	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (j < n - 1) {
				myfile << s[i][j] << ",";
			}
			else {
				myfile << s[i][j];
			}

		}
		if (i < m - 1) {
			myfile << std::endl;
		}
	}

	myfile.close();
}


//=============================================================================================
// Physics
//=============================================================================================

int energy(matrix<int> &s) {
	/* Function that returns the energy with periodic boundary conditions
	of a given spin configuration s in a rectangular grid and J=1.*/


	int  m = s.m; int  n = s.n;  // get size of matrix
	int E = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i < m - 1) {
				// interaction of neighbours
				E -= s[i][j] * s[i + 1][j];
			}
			else {
				// periodic boundary conditions
				E -= s[i][j] * s[0][j];
			}

			if (j < n - 1) {
				// interaction of neighbours
				E -= s[i][j] * s[i][j + 1];
			}
			else {
				// boundary conditions
				E -= s[i][j] * s[i][0];
			}
		};
	};
	return E;
}

int magnetization(matrix<int> &s) {
	/* Function that returns the magnetiation
	of a given spin configuration s in a rectangular grid.*/

	int  m = s.m; int  n = s.n;
	int M = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			M += s[i][j];
		}
	}
	return M;
}

//=============================================================================================
// Unit tests
//=============================================================================================

void unit_test_energy_magnetization() {
	/*Computing energies with boundary conditions of
	2 by 2 lattice for known spin configurations.*/

	// energy 8, M = 0
	matrix<int> s1(2, 2);
	s1[0][0] = 1; s1[0][1] = -1;
	s1[1][0] = -1; s1[1][1] = 1;

	// energy -8, M = 4
	matrix<int> s2(2, 2);
	s2[0][0] = 1; s2[0][1] = 1;
	s2[1][0] = 1; s2[1][1] = 1;

	if (energy(s1) == 8 && energy(s2) == -8) {
		std::cout << "Energy function works as intended." << std::endl;
	}

	if (magnetization(s1) == 0 && magnetization(s2) == 4) {
		std::cout << "Magnetization function works as intended." << std::endl;
	}
}

//=============================================================================================
// Metropolis Algorithm (algorithm & and control functions)
//=============================================================================================




void get_random_matrix(matrix <int> &input, int L) {
	/* Function to fill a matrix with random entries from {0,..., L-1}*/

	int m = input.m; int n = input.n;

	std::random_device                  rand_dev;
	std::mt19937                        generator(rand_dev());
	std::uniform_int_distribution<int>  distr(0, L-1);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			input[i][j] = distr(generator);
		}
	}
}


void initialize_random_spin_configuration(matrix<int> &s) {
	/* Function to initialize the given matrix<int> s as a spin configuration,
	i.e., as a matrix<int> with random entries -1 or 1.*/

	// get size of s
	int m = s.m; int n = s.n;
	// numbers that we pick from
	static int range[] = { -1, 1 };

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {

			// load random number 0 or 1
			int index = rand_int(0, 1);

			// update s[i][j] to random entry from range
			s[i][j] = range[index];
		}
	}
}

int mod(int n, int m) {
	/* Function that returns n (mod m) if we know
	n to lie in the interval {-1,0,...,m-1,m}. Does not work for n
	outside that range.*/

	if (n == -1) {
		return m - 1;
	};
	if (n == m) {
		return 0;
	}
	return n;
}

void monte_carlo_cycle(matrix<int> &s, int &E, int &M, double *w) {
	/* Iterating over all entries in the spin configuration s,
	and flipping a random entry each time and checking whether the
	energy difference is acceptable. If yes, then we update the spin configuration,
	and we update the energy E and the magnetization M.*/

	// get size of s
	int n = s.n; int m = s.m;
	int nn = n * n;

	// get_random_matrix
	matrix<int> random(nn, 2);
	get_random_matrix(random,n);

	for (int i = 0; i < nn; i++) {
			
		//int k = rand_int(0, m - 1);
		//int l = rand_int(0, n - 1);

		int k = random[i][0];
		int l = random[i][1];


		// calculating the change in energy based on flipping s[k][l]
		int delta_E = 2 * s[k][l] * (s[mod((k - 1), m)][l] + s[mod((k + 1), m)][l] + s[k][mod((l - 1), n)] + s[k][mod((l + 1), n)]);

		// if the energy difference is non-positive, then we accept the new spin configuration, 
		// or if the energy difference positive, then check if it satisfies the metropolois condition

		if (delta_E <= 0) {
			// update spin configuration
			s[k][l] *= -1;
			// adding new energy and magnetization
			E += delta_E;
			M += 2 * s[k][l];
		}
		else {
			if (rand_double(0, 1) <= w[delta_E + 8]) {
					// update spin configuration
				s[k][l] *= -1;

				// adding new energy and magnetization
				E += delta_E;
				M += 2 * s[k][l];
			}
		}

	}	
}

void warm_up(matrix<int> &s, double T, int warm_up_cycles) {
	/*Function to warm up a spin configuration so that it gets closer to steady state.*/

	// precalculate exp(delta E) for the 5 possible energy changes
	double w[17];
	for (int delta_E = -8; delta_E <= 8; delta_E++) w[delta_E + 8] = 0;
	for (int delta_E = -8; delta_E <= 8; delta_E += 4) w[delta_E + 8] = exp(-delta_E / T);

	// calculate energy and magnetization of s
	int E = energy(s); int M = magnetization(s);

	// warming up
	for (int cycle = 1; cycle <= warm_up_cycles; cycle++) {
		monte_carlo_cycle(s, E, M, w);

		if (cycle % 10000 == 0) {
			std::cout << "completed warm up cycle: " << cycle << std::endl;
		}
	}
}


void write_quantities_to_file(matrix<int> &s, int cycles, std::string filename, double T) {
	/* Function to run a metropolis algorithm on the incoming spin configuration s
	for the given number of cycles. Then the quantities: <E>, <M>, (sigma_E)^2, (sigma_M)^2, and <|M|>,
	are calculated (exectation values of energy and magnetization, specific heat, susceptibility,
	and expectation value of magnetizatio respectively). After each sycle these quantities
	are saved in a CSV format in the given order to a file labeled, for example,
	"filename_lattice_20_20__cycles_10_T_2.400000.txt", s is 20 by 20, cycles=10 and T=2.4.*/

	// get initial energy and magnetization 
	int E = energy(s); int M = magnetization(s);

	// precalculate exp(delta E) for the 5 possible energy changes
	double w[17];
	for (int delta_E = -8; delta_E <= 8; delta_E++) w[delta_E + 8] = 0;
	for (int delta_E = -8; delta_E <= 8; delta_E += 4) w[delta_E + 8] = exp(-delta_E / T);

	// opening file and formatting filename
	std::ofstream myfile(filename
		+ "_lattice_" + std::to_string(s.m) + "_" + std::to_string(s.n) + "_" +
		"_cycles_" + std::to_string(cycles) +
		"_T_" + std::to_string(T) + ".txt");

	// place to store quantities each cycle
	vector<double> quantities(5);

	// place to store output each cycle
	vector<double> output(6);

	// go through cycles and save quantities to file
	for (int cycle = 1; cycle <= cycles; cycle++) {
		monte_carlo_cycle(s, E, M, w);

		// add newest term in average
		quantities[0] += E; quantities[1] += E * E; // energy
		quantities[2] += M; quantities[3] += M * M; quantities[4] += fabs(M); // magnetization

		// update various expectation values
		double expect_E = quantities[0]/ cycle;
		double expect_E2 = quantities[1] / cycle;
		double expect_M = quantities[2] / cycle;
		double expect_M2 = quantities[3] / cycle;
		double expect_absM = quantities[4] / cycle;

		// store output values into vector<double> output
		output[0] = expect_E;
		output[1] = expect_M;
		output[2] = expect_E2 - expect_E * expect_E;
		output[3] = expect_M2 - expect_absM * expect_absM; // some of the results were made without this line
		output[4] = expect_M2 - expect_M * expect_M;
		output[5] = expect_absM;

		// save content of output as line formated as CSV file
		print_vector_to_file(output, myfile);

		if (cycle % 10000 == 0) {
			std::cout << "completed cycle: " << cycle << std::endl;
		}
	};

	// close file
	myfile.close();
}

//=======================================================================================================

void write_energy_to_file(matrix<int> &s, int cycles, std::string filename, double T) {
	/* Function to run a metropolis algorithm on the incoming spin configuration s
	for the given number of cycles. */

	// get initial energy and magnetization
	int E = energy(s); int M = magnetization(s);

	// precalculate exp(delta E) for the 5 possible energy changes
	double w[17];
	for (int delta_E = -8; delta_E <= 8; delta_E++) w[delta_E + 8] = 0;
	for (int delta_E = -8; delta_E <= 8; delta_E += 4) w[delta_E + 8] = exp(-delta_E / T);

	// open file
	std::ofstream myfile(filename
		+ "energy_lattice_" + std::to_string(s.m) + "_" + std::to_string(s.n) + "_" +
		"_cycles_" + std::to_string(cycles) +
		"_T_" + std::to_string(T) + ".txt");

	// place to store output each cycle
	vector<double> output(2);

	// go through cycles and save quantities to file
	for (int cycle = 1; cycle <= cycles; cycle++) {
		monte_carlo_cycle(s, E, M, w);


		// store output values into vector<double> output
		output[0] = E;

		// save content of output as line formated as CSV file
		print_vector_to_file(output, myfile);

		if (cycle % 100 == 0) {
			std::cout << "completed cycle: " << cycle << std::endl;
		}

	};

	// close file
	myfile.close();
}

//=======================================================================================================

void monte_carlo_cycle_pass_counter(matrix<int> &s, int & counter, double *w) {
	/* Iterating over all entries in the spin configuration s,
	and flipping a random entry each time to check whether the
	energy difference is acceptable. If yes, then we update the spin configuration,
	and we update the counter*/

	// get size of s
	int m = s.m; int n = s.n;
	int iterations = m * n;

	for (int iteration = 0; iteration < iterations; iteration++) {
		int k = rand_int(0, m - 1);
		int l = rand_int(0, n - 1);

		// calculating the change in energy based on flipping s[k][l]
		int delta_E = 2 * s[k][l] * (s[mod((k - 1), m)][l] + s[mod((k + 1), m)][l] + s[k][mod((l - 1), n)] + s[k][mod((l + 1), n)]);

		// or if the energy difference positive, then check if it satisfies the metropolis condition
		if (delta_E <= 0) {
			// update spin configuration
			s[k][l] *= -1;
			counter += 1;

		}
		else {
			if (rand_double(0, 1) <= w[delta_E + 8]) {
				// update spin configuration
				s[k][l] *= -1;
				counter += 1;
			}
		}
	}
}

void write_passes_vs_cycles_to_file(matrix<int> &s, int cycles, std::string filename, double T) {

	// precalculate exp(delta E) for the 5 possible energy changes
	double w[17];
	for (int delta_E = -8; delta_E <= 8; delta_E++) w[delta_E + 8] = 0;
	for (int delta_E = -8; delta_E <= 8; delta_E += 4) w[delta_E + 8] = exp(-delta_E / T);

	// open file
	std::ofstream myfile("passes_" + filename +
		"_lattice_" + std::to_string(s.m) + "_" + std::to_string(s.n) + "_" +
		"_cycles_" + std::to_string(cycles) +
		"_T_" + std::to_string(T) + ".txt");

	int counter = 0;
	vector<int> output(1);

	// go through cycles and save quantities to file
	for (int cycle = 1; cycle <= cycles; cycle++) {
		monte_carlo_cycle_pass_counter(s, counter, w);
		output[0] = counter;
		// save content of output as line formated as CSV file
		print_vector_to_file(output, myfile);

		if (cycle % 100 == 0) {
			std::cout << "completed cycle: " << cycle << std::endl;
		}
	};

	// close file
	myfile.close();
}


//===============================================================================================
// Main program
//===============================================================================================

int main() {
	int L, cycles;

	// run unit tests
	unit_test_energy_magnetization();

	//===========================================================================================
	// Exercise b)
	//===========================================================================================

	// save quantities to file (L=2)
	L = 2;
	cycles = 100000;
	matrix<int> s(L, L);
	initialize_random_spin_configuration(s);
	write_quantities_to_file(s, cycles, "small", 1.0);

	//===========================================================================================
	// Exercise c)
	//===========================================================================================

	//// generate data for different temperatures
	//vector<double> temperatures(2);
	//temperatures[0] = 1.0;
	//temperatures[1] = 2.4;

	//for (int index = 0; index < 2; index++) {
	//	// save quantities to file (L=20) random
	//	int L = 20;
	//	int cycles = 100000;
	//	double T = temperatures[index];
	//		
	//	matrix<int> s(L, L);
	//	//random initial condition
	//	initialize_random_spin_configuration(s);

	//	write_quantities_to_file(s, cycles, "random", T);

	//	//random initial condition
	//	initialize_random_spin_configuration(s);

	//	write_passes_vs_cycles_to_file(s, cycles, "random", T);

	//	// save quantities to file (L=20) oriented

	//	//orient s
	//	for (int i = 0; i < L; i++) {
	//		for (int j = 0; j < L; j++) {
	//			s[i][j] = 1;
	//		}
	//	}

	//	write_quantities_to_file(s, cycles, "oriented", T);
	//	for (int i = 0; i < L; i++) {
	//		for (int j = 0; j < L; j++) {
	//			s[i][j] = 1;
	//		}
	//	}
	//	write_passes_vs_cycles_to_file(s, cycles, "oriented", T);
	//}

	//===========================================================================================
	// Exercise d)
	//===========================================================================================

	//// generate data for different temperatures
	//vector<double> temperatures(2);
	//temperatures[0] = 1.0;
	//temperatures[1] = 2.4;

	//for (int index = 0; index < 2; index++) {
	//	// save quantities to file (L=20) random
	//	int L = 20;
	//	int cycles = 1000000;
	//	double T = temperatures[index];

	//	matrix<int> s(L, L);
	//	//random initial condition
	//	initialize_random_spin_configuration(s);

	//	write_energy_to_file(s, cycles, "random_", T);
	//}

	//===========================================================================================
	// Exercise e)
	//===========================================================================================

	//// generate vector<double> containing different temperatures [2.00,2.05,...,2.35]
	//vector<double> temperatures(8);
	//for (int index = 0; index < 8; index++) {
	//	temperatures[index] = 2.0 + 0.05*index;
	//}

	//// generate vector<int> containing different lattice sizes
	//vector<int> lattice_sizes(4);
	//lattice_sizes[0] = 80;
	//lattice_sizes[1] = 60;
	//lattice_sizes[2] = 80;
	//lattice_sizes[2] = 100;

	//// iterate over lattice sizes
	//for (int i = 0; i < 1; i++) {
	//	int L = lattice_sizes[i];
	//	// parallelize calculation for different temperatures
	//    #pragma omp parallel for
	//	for (int j = 0; j < 8; j++) {
	//		// save quantities to file
	//		int warm_up_cycles = 200000;
	//		int cycles = 1500000;
	//		double T = temperatures[j];

	//		matrix<int> s(L, L);

	//		// random initial condition
	//		initialize_random_spin_configuration(s);

	//		// warm up spin configuration
	//		warm_up(s, T, warm_up_cycles);

	//		// save warmed up spin configuration (for later use?)
	//		std::string filename = "warmed_up_matrix_" + std::to_string(L) + "_" + std::to_string(L) + "_" +
	//			"_cycles_" + std::to_string(warm_up_cycles) +
	//			"_T_" + std::to_string(T) + ".txt";
	//		save_matrix_file(s, filename);

	//		// write quantities to file
	//		write_quantities_to_file(s, cycles, "random", T);
	//	}
	//}

	//// timing of parallelization for rectangular lattice L=10
	//int N = 10;
	//vector<double> times(N);
	//double average_time = 0;
	//for (int lap = 0; lap < N; lap++){
	//	using namespace std::chrono;
	//	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//	#pragma omp parallel for //commented/uncomment to compare speed
	//	for (int j = 0; j < 6; j++) {
	//		// save quantities to file
	//		int cycles = 10000;
	//		double T = temperatures[j];

	//		matrix<int> s(10, 10);

	//		// random initial condition
	//		initialize_random_spin_configuration(s);

	//		// write quantities to file
	//		write_quantities_to_file(s, cycles, "random", T);
	//	}
	//	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	//	times[lap] = time_span.count();
	//	average_time += time_span.count();
	//}
	//std::cout << average_time/N << std::endl;

	//===========================================================================================
	// Success message
	//===========================================================================================

	std::cout << "Success!";
	std::cin.get();
}