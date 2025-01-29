#include <iostream>
#include <random>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <fstream>
#include <array>
#include <cmath>
using namespace std;

// TODO : McKean-Vlasov

double my_unif(){
	double number = rand();
	return number / RAND_MAX;
};

double my_gaussian(){
	double first, second, radius;
	first = my_unif();
	second = my_unif();

	radius = sqrt(-2.0 * log(first));
	return radius * sin(2.0 * 3.14 * second);
};

std::vector<double> two_gaussians(){
	double first, second, radius;
	first = my_unif();
	second = my_unif();

	radius = sqrt(-2.0 * log(first));

	std::vector<double> gaus_array;
	gaus_array.push_back(radius * sin(2.0 * 3.14 * second));
	gaus_array.push_back(radius * cos(2.0 * 3.14 * second));
	return gaus_array;
};

double identity(double x){
	return x;
};

// necessary?
double my_div(int k, double m){
	return (double) k / m;
};

double my_div(double m, int k){
	return m / (double) k;
};

double my_div(int m, int k){
	return (double) m / (double) k;
}

double my_abs(double x){
	if (x > 0)
		return x;
	return -x;
}

// here we are going to store matrices. We intend to keep them as 1d objects that will pretend to be 2D or 3D
class Array{
public:
	std::vector<double> data;
	int n_dim;
	std::vector<int> dimensions;

	// possible constructors:
	Array(int);
	Array(int, int);
	Array(int, int, int);
	Array(int, std::vector<double>);
	Array(int, int, std::vector<double>);
	Array(int, int, int, std::vector<double>);

	// general number of dimensions, although not sure why would I need it
	Array(std::vector<int>);
	Array(std::vector<int>, std::vector<double>);

	// basic operations with overloading go here
	double get_value(int);
	double get_value(int, int);
	double get_value(int, int, int);
	double get_value(std::vector<int>);

	// filling in with data goes here
	void insert_value(int, double);
	void insert_value(int, int, double);
	void insert_value(int, std::vector<double>);

	// debugging and output
	void dump_to_textfile(std::string);
	void print_first_values(int);

	// operators
	Array& operator*=(const double& mult){
		for (int i = 0; i < data.size(); i++)
			data[i] *= mult;
		return *this;
	}
};

Array::Array(int size){
	n_dim = 1;

	this->dimensions.push_back(size);
	for (int i = 0; i < size; i++){
		this->data.push_back(0);
	};
};

Array::Array(int sizex, int sizey){
	this->n_dim = 2;
	this->dimensions.push_back(sizex);
	this->dimensions.push_back(sizey);
	for (int i = 0; i < sizex * sizey; i++)
		this->data.push_back(0);
};

double Array::get_value(int i){
	assert(i >= 0 && i < this->dimensions[0]);
	return this->data[i];
};

double Array::get_value(int line, int column){
	return this->data[ this->dimensions[1] * line + column ];
}

void Array::insert_value(int i, double val){
	assert(i >= 0 && i < this->dimensions[0]);
	this->data[i] = val;
}

void Array::insert_value(int i, int j, double val){
	assert(this->n_dim == 2);
	this->data[ this->dimensions[1] * i + j ] = val;
};

void Array::insert_value(int line, std::vector<double> values){
	assert(this->n_dim == 2 and values.size() == this->dimensions[1]);

	for (int j = 0; j < this->dimensions[1]; j++){
		this->data[ this->dimensions[1] * line + j ] = values[j];
	};
};

void Array::dump_to_textfile(std::string filename){

	ofstream im_file;
	im_file.open(filename);
	if (this->n_dim == 1){
		for (int i = 0; i < this->data.size(); i++){
			im_file << this->data[i] << " ";
		};
	};
	if (this->n_dim == 2) {
		for (int i = 0; i < this->dimensions[0]; i++){
			for (int j = 0; j < this->dimensions[1]; j++){
				im_file << this->data[ i * this->dimensions[1] + j ] << " ";
			};
			im_file << endl;
		};
	};
	im_file.close();
};

double mean_tot(Array a){
	int tot = a.data.size();
	double m(0);
	for (int i = 0; i < tot; i++){
		m += a.data[i];
	};
	return my_div(m, tot);
};

double mean_tot(std::vector<double> a){
	double m(0);
	int tot = a.size();
	for (int i = 0; i < tot; i ++){
		m += a[i];
	};
	return my_div(m, tot);
};

double my_min(double x, double y){
	if (x<y)
		return x;
	return y;
};

void Array::print_first_values(int number){

	if (this->n_dim == 1){
		for (int i = 0; i < my_min(number, this->data.size()); i++){
			cout << this->data[i] << " ";
		};
	};
	if (this->n_dim == 2) {
		for (int i = 0; i < my_min(number, this->dimensions[0]); i++){
			for (int j = 0; j < this->dimensions[1]; j++){
				cout << this->data[ i * this->dimensions[1] + j ] << " ";
			};
			cout << endl;
		};
	};
};

// here we will simulate fBm in a couple of different ways. Possibly we will move it out to a separate file
Array sample_brownian_motion(int N){
	Array values(N);
	double timestep;
	timestep = 1/(double) N;
	timestep = sqrt(timestep);
	for (int i = 1; i < N; i++){
		values.insert_value(i, values.get_value(i-1) + timestep * my_gaussian());
	};

	return values;
};

// samples of fBm following Picard's "Representation formulae for fBm". Note that this can only work on the interval of [0,1]
// TODO : write a simulation for an arbitrary interval so that you can simulate long-time behaviour
Array sample_fbm_fourier(double hurst, int N, int cutoff){
	Array values(N);
	Array gaussians(cutoff, 2);
	double timestep, power;
	timestep = 1/(double) N;
	power = hurst + 0.5;

	// By Remark 6.7 in Picard's we can compute much simpler process than (77)
	for (int i = 0; i < cutoff; i++){
		gaussians.insert_value(i, two_gaussians());
	};

	// cout << gaussians.get_value(1, 0) << " " << gaussians.get_value(1, 1) << endl;

	for (int i = 0; i < N; i++){
		double loc_val(0), time;
		time = timestep * i;
		loc_val = 0;
		for (int j = 1; j < cutoff+1; j++){
			double contribution(0);
			contribution = gaussians.get_value(j-1, 0) * sin(3.14 * j * time) + gaussians.get_value(j-1, 1) * (cos(3.14 * j * time) - 1);
			contribution = contribution / pow((3.14 * j), power);
			loc_val += contribution;
		};
		values.insert_value(i, loc_val);
	};

	return values;
};

// create iterated integral, if possible maybe this could be a method in subclass Path
// also in Path you could have sth like "get increment"
// Array generate_roughpath(){

// };

// this function will take a function for the drift, the function for the noise and a path
Array integrate_simple_sde(std::function<double(double)> drift, double volat, int N, Array path, double init_val = 0){
	double timestep = 1/ (double) N;
	double movement(0);
	Array results(N);
	results.insert_value(0, init_val);

	for (int i = 1; i < N; i++){
		movement = drift(results.get_value(i-1)) * timestep + volat * (path.get_value(i)-path.get_value(i-1));
		results.insert_value(i, results.get_value(i-1) + movement);
	};

	return results;
};

Array integrate_sde_splitting(std::function<double(double, double)> flow, double volat, int N, Array path, double init_val = 0){
	double timestep = my_div(1., N);
	double movement;
	Array results(N);
	results.insert_value(0, init_val);

	for (int i = 1; i < N; i++){
		movement = volat * (path.get_value(i)-path.get_value(i-1));
		results.insert_value(i, flow(results.get_value(i-1) + movement, timestep));
	};
	return results;
};

// to be rewritten with new noises ? we can easily add drift as well 
Array integrate_ips_sde(std::function<double(double, double)> interaction, int Ninter, int Nsteps, double mult = 1, int report_progress = 100){
	double timestep = 1 / (double) Nsteps;

	// first we need to generate noises
	Array noises(Ninter, Nsteps);
	Array values(Ninter, Nsteps);
	for (int i = 0; i < Ninter; i++){
		noises.insert_value(i, sample_brownian_motion(Nsteps).data);
	};

	noises *= mult;

	// inside the loop we have another loop which is responsible for encoding interactions
	for (int t = 1; t < Nsteps; t++){
		for (int i = 0; i < Ninter; i++){

			double movement(0);
			for (int j = 0; j < Ninter; j++){
				if (j != i){
					movement += interaction(values.get_value(i, t-1), values.get_value(j, t-1)) * timestep;
				};
			};
			movement = my_div(movement, Ninter-1);
			movement += noises.get_value(i, t) - noises.get_value(i, t-1);
			// movement += drift(values.get_value(i-1, t-1)) * timestep;

			values.insert_value(i, t, values.get_value(i, t-1) + movement);
		};

		if (t % report_progress == 0)
			cout << "Timestep " << t << ", which is " << 100 * my_div(t, Nsteps) << "% of progress" << endl;
	};

	return values;
};

// idea: P( X \in [x_i, x_{i+1}]) = \int_{x_i}^{x_{i+1}} g(y) dy ~ g(x_i) (x_{i+1} - x_i)
// which means that we will divide the area in bins of length binsize and use:
// P( X \in [x_i, x_{i+1}]) = # nr of X in [x_i, x_{i+1}] / # total nr of X ~ g(x_i) (x_{i+1} - x_i)
// so that g(x_i) ~ 1/(x_{i+1} - x_i) * # ( X \in [x_i, x_{i+1}]) / # total nr of X
// in the end this is just histogram, it's just quicker to write a new implementation than look for it on stackoverflow and deal with adapting it

Array compute_density_profile(std::vector<double> data, double binsize, bool replace_borders = 0, double minval = -1, double maxval = 1){
	// we pass by copy on purpose not to ruin the order of the original data

	std::sort(data.begin(), data.end());

	double miv, mav;
	int data_length = data.size();
	if (replace_borders){
		miv = minval;
		mav = maxval;
	} else {
		// at this point they are already sorted
		miv = data[0];
		mav = data.back();
	};

	int N = (mav - miv)/ binsize + 1;

	if (N > data_length){
		cout << "number of bins is higher than number of datapoints" << endl;
	};

	Array histogram(N);
	int curr_index(0), curr_hist_index(0), count(0);
	double curr_lhs;
	curr_lhs = miv;
	// note that we will put the last bin by hand because otherwise we add one count too much and we run into a segfault

	while (curr_hist_index < N-1 and curr_index < data_length){
		while (data[curr_index + count] < curr_lhs + binsize)
			count++;

		curr_index = curr_index + count;
		histogram.insert_value(curr_hist_index, my_div(count, data_length));

		curr_hist_index++;
		curr_lhs += binsize;
		count = 0;
	};

	histogram.insert_value(curr_hist_index, my_div(data_length - curr_index, data_length));

	return histogram;
};

double integrate_function(std::function<double(double)> f, int N){
	double val(0);
	double timestep;
	timestep = 1/(double) N;
	for (int i = 0; i < N; i++){
		val += f(i * timestep) * timestep;
	};

	return val;
};

// run integration N times and take the average, would be good to implement changing noise.
std::vector<double> integrate_montecarlo_return_all(int N, int Nint, std::function<double(Array)> f, function<double(double)> drift){
	std::vector<double> val;
	for (int i = 0; i < N; i++){
		Array path = sample_brownian_motion(Nint);
		Array results = integrate_simple_sde(drift, 1, Nint, path);
		val.push_back(f(results));
	};
	return val;
};

double integrate_montecarlo(int N, int Nint, function<double(Array)> f, function<double(double)> drift){
	return mean_tot(integrate_montecarlo_return_all(N, Nint, f, drift));
};

// 

// dumb functions created to be passed as arguments

// value of -0.5 is hardcoded, not sure how to make it passable atm, will check
double power_fn(double x){
	if (my_abs(x) < 0.001)
		return 0;
	if (x > 0)
		return my_div(1, pow(my_abs(x), 0.5));
	return -my_div(1, pow(my_abs(x), 0.5));
};

double interaction_power(double x, double y){
	return power_fn(x-y);
}

// it stems from x_t = ( x_0^(1-g) + t)^(1/(1-g)), this one corresponds to g = -1/2, but I plan on extending it to variable g
double flow_fn(double x, double h){
	return pow( pow(x, 1.5) + h, 0.67 );
};

int main(){
	srand((unsigned) time(NULL));

	int Ninter(10000), Nsteps(1000);
	Array values = integrate_ips_sde(interaction_power, Ninter, Nsteps, 0.01);

	values.dump_to_textfile("result_ips.txt");



	// Array u = sample_brownian_motion(1000);
	// Array t = sample_fbm_fourier(0.25, 1000, 1000);
	// Array result(1000);
	// t.dump_to_textfile("fbm.txt");
	// u.dump_to_textfile("bm.txt");
	// double epsilon = 1;
	// std::string filename;

	// for (int e = 1; e < 4; e++){
	// 	epsilon = my_div(1, pow(10, e-1));
	// 	cout << "considering " << epsilon << endl;
	// 	cout << "power_fn test" << power_fn(0) << endl;
	// 	result = integrate_simple_sde(power_fn, epsilon, 1000, u); // double volat, int N, Array path, double init_val = 0
	// 	filename = "integr_simple_" + std::to_string(e) + ".txt";
	// 	cout << endl;
	// 	result.dump_to_textfile(filename);

	// 	result = integrate_sde_splitting(flow_fn, epsilon, 1000, u);
	// 	filename = "integr_split_" + std::to_string(e) + ".txt";
	// 	result.dump_to_textfile(filename);
	
};