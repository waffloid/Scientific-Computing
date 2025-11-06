#include<functional>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include <cassert>
#include<vector>
#include<filesystem>
#include<cmath>

#include "IterativeRule.h"


// changes todo:
// create a separate file for the initial conditions
// create a folder for each of the "experiments"
//	and have each file be an experiment in there
//


using StateVector = std::vector<double>;
// Forward declaration

struct NumericalSchemeParams {
	int N;
	double dt = 0.01;
	std::function<double(double)> init_conditions;
	IterativeRule& iterative_rule;
	std::string& output_dir;
	bool transmissive_boundary = true;
	int log_freq=1;
};

class StateData {
    int N;
    StateVector velocity;
	bool transmissive_boundary;
	
    public: 
        StateData(int N, bool transmissive) : N(N), transmissive_boundary(transmissive) {
			if (transmissive_boundary) {
				velocity = std::vector<double>(N+2);
			} else {
				velocity = std::vector<double>(N);
			}
		}

        void set_velocity(int i, double vel) {
            if (transmissive_boundary) {
                velocity[i+1] = vel;

                // ghost boundary -- if we set velocity at boundary, make the ghost cells match
                if (i==0) {
                    velocity[0] = vel;
                } else if (i==N-1) {
                    velocity[N+1] = vel;
                }
            } else {
                velocity[i % N] = vel;
            }
        }

        double get_velocity(int i) const {
			if (transmissive_boundary) {
            	return velocity.at(i+1);
			} else {
				return velocity.at(((i % N) + N) % N);
			}
        }
};

// pass in dx, dt, N, initial value as a function of x(?)
// not our responsibility to 

class NumericalScheme {
    StateData data;
	IterativeRule& iterative_rule;
    
	
	int iter_count = 0;
	int log_freq;
	int N;
	double dx;
	double dt;
	std::string& output_dir;

	bool transmissive_boundary;

	public:

   		std::function<double(double, double)> analytic_solution;

		NumericalScheme(NumericalSchemeParams params) 
			: data(StateData(params.N, params.transmissive_boundary)), 
			  dx(1/double(params.N)),
              N(params.N),
			  log_freq(params.log_freq), 
			  dt(params.dt), 
			  iterative_rule(params.iterative_rule),
			  output_dir(params.output_dir),
              transmissive_boundary(params.transmissive_boundary)
        {

			std::filesystem::remove_all("Output/" + output_dir);
			std::filesystem::create_directories("Output/" + output_dir + "/" + "simulation");
			std::filesystem::create_directories("Output/" + output_dir + "/" + "error");
			std::filesystem::create_directories("Output/" + output_dir + "/" + "error_l1");

			for(int i=0;i<N;++i){
				double x = double(i) / double(N);
				double velocity = params.init_conditions(x);

				data.set_velocity(i, velocity);
			}
			log_init_condition();
		}

		void run_timestep() {
			StateData next_data = StateData(N, transmissive_boundary);

			for (int i = 0; i < N; ++i) {
				// Build stencil vector
				std::vector<double> stencil;
				for (int offset = iterative_rule.left; offset <= iterative_rule.right; ++offset) {
					stencil.push_back(data.get_velocity(i + offset));
				}

				double new_velocity = iterative_rule.apply(stencil, dx, dt);
				next_data.set_velocity(i, new_velocity);
			}
			log_velocity();
			//log_error();
			//log_l1_error();

			data = next_data;
			iter_count++;
		}

		double error_L1() const {
			double error = 0.0;

			for (int i = 0; i<N; ++i) {
				double x = i * dx;
				double t = iter_count * dt;
				double u_exact = analytic_solution(x, t);
				error += abs(data.get_velocity(i) - u_exact);
			}

			return error * dx;
		}

	private:
		void log_data(const std::string& dir,  const std::vector<double>& data) const {
			std::ofstream file(
				"Output/" + output_dir + "/" + dir  + "/" + std::to_string(iter_count)  + ".dat"
			);

			for (const double & value : data){
				file << value << "\n";
			}
			file.close();
		}
		
		void log_error() const {
			const std::string directory_name = "error";
			std::vector<double> error;

			for (int i = 0; i < N; ++i) {
				double x = i * dx;
				double t = iter_count * dt;
				double u_exact = analytic_solution(x, t);

				error.emplace_back(u_exact - data.get_velocity(i));
			}

			log_data(directory_name, error);


		}
		
		void log_l1_error() const {
			const std::string directory_name = "error_l1";
			const std::vector<double> error_vector = {error_L1()};
			log_data(directory_name, error_vector);
		}

		void log_init_condition() const {
			std::ofstream file("Output/" + output_dir + "/" + "initial_condition.dat");
	
			for (int i = 0; i < N; ++i) {
				file << data.get_velocity(i) << "\n";
			}
			file.close();
		}

		void log_velocity() const {
			std::ofstream file(
				"Output/" + output_dir + "/simulation/" 
				+ std::to_string(iter_count) + ".dat"
			);


			for (int i = 0; i < N; ++i){
				file <<  data.get_velocity(i) << "\n";
			}
			file.close();
		}
};


//// INIT CONDITIONS

double IndicatorInitCondition(double x) {
	if ((0.25 < x) && (x < 0.75)) {
		return 1.0;
	}

	return 0.0;
};

double GaussianInitialCondition(double x) {
	return exp(-16 * (x-0.5) * (x-0.5));
};

std::function<double(double, double)> AdvectionSolution(
		std::function<double(double)> init_condition, 
		double speed=1.0
) {
	return [init_condition, speed](double x, double t) {
		double characteristic_source = x - speed*t;
		double interval_contained = floor(x - speed*t);
		double location = characteristic_source - interval_contained;

		return init_condition(location);
	};
}

////

// NONSENSE FUNCTIONS


void run_burger_lax_friedrichs() {
	Burger_Lax_Friedrichs iter_rule;
	std::string output_dir = "Burger Lax Friedrichs";

	NumericalSchemeParams params = {
		.N = 100,
		.dt = 0.01,
		.init_conditions = GaussianInitialCondition,
		.iterative_rule = iter_rule,
		.output_dir = output_dir,
		.transmissive_boundary = true,
		.log_freq=1
	};

	NumericalScheme scheme = NumericalScheme(params);
	
	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}

}

void run_burger_richtimeyer() {
    Burger_Richtimeyer iter_rule;
    std::string output_dir = "Burger Richtimeyer";

	NumericalSchemeParams params = {
		.N = 100,
		.dt = 0.01,
		.init_conditions = GaussianInitialCondition,
		.iterative_rule = iter_rule,
		.output_dir = output_dir,
		.transmissive_boundary = true,
		.log_freq=1
	};

	NumericalScheme scheme = NumericalScheme(params);
	
	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}
}


void run_burger_FORCE() {
    Burger_FORCE iter_rule;
    std::string output_dir = "Burger FORCE";

	NumericalSchemeParams params = {
		.N = 100,
		.dt = 0.01,
		.init_conditions = GaussianInitialCondition,
		.iterative_rule = iter_rule,
		.output_dir = output_dir,
		.transmissive_boundary = true,
		.log_freq=1
	};

	NumericalScheme scheme = NumericalScheme(params);
	
	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}

}


int main() {
	run_burger_lax_friedrichs();
    run_burger_richtimeyer();
    run_burger_FORCE();

	return 0;
};
