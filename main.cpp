#include <iostream>
#include <random>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <fstream>
#include <array>
#include <cmath>
#include <unistd.h>
using namespace std;

// TODO : Davies-Harty dla fBm (bo obecne nie zadziala dla T > 1) ALBO przerobic Fouriera przez skalowanie T^H B_s =_{\PP} B_t dla s\in[0,1] i t \in[0,T] 
// TODO : symulacje miary niezmienniczej w dlugim czasie (zrobione?)

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

// necessary? is it possible to use templates or something like it ?
double my_div(int k, double m){
	return (double) k / m;
};

double my_div(double m, int k){
	return m / (double) k;
};

double my_div(int m, int k){
	return (double) m / (double) k;
}

template <typename T>
T my_abs(T x){
	if (x > 0)
		return x;
	return -x;
}

template <typename T>
T my_prod(std::vector<T> x){
	int N, r(1);
	N = x.size();
	for (int i = 0; i <N; i++)
		r *= x[i];
	return r;
};

// here we are going to store matrices. We intend to keep them as 1d objects that will pretend to be 2D or 3D
class Array{
public:
	std::vector<double> data;
	int n_dim;
	std::vector<int> dimensions;

	// possible constructors:
	Array();
	Array(int);
	Array(int, int);
	Array(int, int, int);
	Array(int, std::vector<double>);
	Array(int, int, std::vector<double>);
	Array(int, int, int, std::vector<double>);

	// general number of dimensions
	Array(std::vector<int>);
	Array(std::vector<int>, std::vector<double>);

	// getters
	double get_value(int);
	double get_value(int, int);
	double get_value(int, int, int);
	double get_value(std::vector<int>);

	// filling in with data goes here
	void insert_value(int, double);
	void insert_value(int, int, double);
	void insert_value(int, std::vector<double>);

	// debugging and output
	void dump_to_textfile(std::string, bool);
	void print_first_values(int);

	// operators
	Array& operator*=(const double& mult){
		for (int i = 0; i < data.size(); i++)
			data[i] *= mult;
		return *this;
	}
};

Array::Array(){
	n_dim = 0;
	this->dimensions = std::vector<int>();
	this->data = std::vector<double>();
}

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

Array::Array(std::vector<int> dim_vec){
	if (dim_vec.size() == 1){
		this->dimensions.push_back(dim_vec[0]);
	} else if (dim_vec.size() == 2){
		this->dimensions.push_back(dim_vec[0]);
		this->dimensions.push_back(dim_vec[1]);
	} else if (dim_vec.size() == 3){
		this->dimensions.push_back(dim_vec[0]);
		this->dimensions.push_back(dim_vec[1]);
		this->dimensions.push_back(dim_vec[2]);
	} else {
		throw std::invalid_argument("dimension vector too large, currently we accept only up to three dimensions");
	};
};

double Array::get_value(int i){
	assert(i >= 0 && i < this->dimensions[0]);
	return this->data[i];
};

double Array::get_value(int line, int column){
	return this->data[ this->dimensions[1] * line + column ];
};

void Array::insert_value(int i, double val){
	assert(i >= 0 && i < this->dimensions[0] && this->n_dim == 1);
	this->data[i] = val;
};

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

void Array::dump_to_textfile(std::string filename, bool append){

	ofstream im_file;

	if (append)
		im_file.open(filename, std::ios_base::app);
	else
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
	im_file << endl;
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

template <typename T>
T my_min(T x, T y){
	if (x<y)
		return x;
	return y;
};

void Array::print_first_values(int number){

	if (this->n_dim == 1){
		for (int i = 0; i < my_min<int>(number, this->data.size()); i++){
			cout << this->data[i] << " ";
		};
	};
	if (this->n_dim == 2) {
		for (int i = 0; i < my_min<int>(number, this->dimensions[0]); i++){
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
Array sample_fbm_fourier(double hurst, double timestep, int cutoff){
	int N;
	N = (int) 1/timestep;
	Array values(N);
	Array gaussians(cutoff, 2);
	double power;
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

// this function will take a function for the drift, the function for the noise and a path
Array integrate_simple_sde(std::function<double(double)> drift, double volat, int N, Array path, double init_val = 0){
	double timestep = 1/ (double) N;
	double movement(0);
	if (path.n_dim == 1){

		Array results(N);
		results.insert_value(0, init_val);

		for (int i = 1; i < N; i++){
			movement = drift(results.get_value(i-1)) * timestep + volat * (path.get_value(i)-path.get_value(i-1));
			results.insert_value(i, results.get_value(i-1) + movement);
		};

		return results;
	} else if (path.n_dim == 2){

		Array results(path.dimensions);

	} else {
		throw std::invalid_argument("path array has incorrect number of dimensions, expected 1 or 2");
	}
};

Array integrate_sde_splitting(std::function<double(double, double, double, double)> flow, double volat, double timestep, 
					double gamma, double drift_multip, Array &path, double init_val = 0){
	double movement;
	double val;
	int N = (int)(1/timestep);
	Array results(N);
	results.insert_value(0, init_val);

	for (int i = 1; i < N; i++){
		movement = volat * (path.get_value(i)-path.get_value(i-1));
		val = results.get_value(i-1) + movement;
		if (val > 0)
			results.insert_value(i, flow(val, timestep, gamma, drift_multip));
		else
			results.insert_value(i, -flow(-val, timestep, gamma, 1));
	};
	return results;
};

// to be rewritten with new noises ? we can easily add drift as well 
Array integrate_ips_sde(std::function<double(double, double, double)> interaction, int Ninter, 
	double timestep, double mult = 1, int report_progress = 100, double g = -0.5, double hurst=0.5, bool bm=1, int cutoff = 1000){
	int Nsteps;

	// first we need to generate noises
	Array noises(Ninter, Nsteps);
	Array values(Ninter, Nsteps);
	if (bm){
		for (int i = 0; i < Ninter; i++){
			noises.insert_value(i, sample_brownian_motion(Nsteps).data);
		};
	} else{
		for (int i = 0; i < Ninter; i++){
			noises.insert_value(i, sample_fbm_fourier(hurst, timestep, cutoff));
		}
	}

	noises *= mult;

	// inside the loop we have another loop which is responsible for encoding interactions
	for (int t = 1; t < Nsteps; t++){
		for (int i = 0; i < Ninter; i++){

			double movement(0);
			for (int j = 0; j < Ninter; j++){
				if (j != i){
					movement += interaction(values.get_value(i, t-1), values.get_value(j, t-1), g) * timestep;
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

// value of -0.5 is hardcoded, not sure how to make it passable atm, will check
double power_fn(double x, double g){
	if (my_abs(x) < 0.001)
		return 0;
	if (x > 0)
		return my_div(1, pow(my_abs(x), g));
	return -my_div(1, pow(my_abs(x), g));
};

double interaction_power(double x, double y, double g){
	return power_fn(x-y, g);
};

// it stems from x_t = ( x_0^(1-g) + t)^(1/(1-g)), this one corresponds to g = -1/2, but I plan on extending it to variable g
double flow_fn(double x, double h, double g, double multip){
	return pow( pow(x, 1-g) + h, multip * my_div(1, 1-g) );
};

void run_zero_noise(int N, double H, double g, double timestep, double eps, double Amult, int cutoff, std::string result_fname){
	Array noise, sde;
	int count_positive;
	int length = (int)(1/timestep);
	for (int i = 0; i < N; i++){
		noise = sample_fbm_fourier(H, timestep, cutoff);
		sde = integrate_sde_splitting(flow_fn, eps, timestep, 
					g, Amult, noise);
		sde.dump_to_textfile(result_fname, 1);
		if (sde.get_value(length-1) > 0)
			count_positive++;
	};
};

int main(int argc, char** argv){
	srand((unsigned) time(NULL));

	double g, timestep, eps(1), H, Amult;
	int N(100), cutoff(1000);
	bool interacting_system(0);
	std::string fname;

	for (int i = 1; i < argc; i++){
		if (std::string(argv[i]) == "-g")
			g = atof( argv[i+1]);
		if (std::string(argv[i]) == "-eps")
			eps = atof( argv[i+1]);
		if (std::string(argv[i]) == "-ts")
			timestep = atof(argv[i+1]);
		if (std::string(argv[i]) == "-H")
			H = atof(argv[i+1]);
		if (std::string(argv[i]) == "-N")
			N = atoi(argv[i+1]);
		if (std::string(argv[i]) == "--ips")
			interacting_system = 1;
		if (std::string(argv[i]) == "-A")
			Amult = atof(argv[i+1]);
		if (std::string(argv[i]) == "-cutoff")
			cutoff = atoi(argv[i+1]);
		if (std::string(argv[i]) == "-fname")
			fname = argv[i+1];
	};

	cout << "gamma " << g << ", timestep " << timestep << ", H " << H << ", eps " << eps << ", N " << N << ", cutoff " << cutoff << endl;

	if (interacting_system){
		Array interactions = integrate_ips_sde(interaction_power, N, 
			timestep, eps, N / 10, g, H, cutoff);
		interactions.dump_to_textfile(fname, 0);
	} else {
		run_zero_noise(N, H, g, timestep, eps, Amult, cutoff, fname);
	}
};