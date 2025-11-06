#include<functional>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include <cassert>
#include<vector>
#include<filesystem>
#include<cmath>

#include "IterativeRule.h"


using StateVector = std::vector<double>;
// Forward declaration

class StateData {
    int N;
    StateVector velocity;
	
    public: 
        StateData(int N) : N(N), velocity(N) {}

        void set_velocity(int i, double vel) {
            velocity[i] = vel;
        }

        double get_velocity(int i) const {
            return velocity.at(i);
        }

        StencilView get_stencil(int left, int index, int right) const {
            StencilView stencil = StencilView(velocity, left, index, right, N);

            return stencil;
        }
};

// pass in dx, dt, N, initial value as a function of x(?)
// not our responsibility to 

class NumericalScheme {
    StateData data;
	IterativeRule& rule;
	std::string simulation_name;
	
	int N;
	double dx;
	double dt;

	int iter_count = 0;
	int log_freq = 1;

	public:

   		std::function<double(double, double)> analytic_solution;
		
		NumericalScheme(int N, std::function<double(double)> init_conditions,
					IterativeRule& rule, std::string& name, double dt=0.01, int log_freq = 1)
			: data(StateData(N)), N(N), dx(1/double(N)), log_freq(log_freq), dt(dt), rule(rule), simulation_name(name) {

			std::filesystem::remove_all(name);
			std::filesystem::create_directories(name + "/" + "simulation");
			std::filesystem::create_directories(name + "/" + "metadata");

			for(int i=0;i<N;++i){
				double x = double(i) / double(N);
				double velocity = init_conditions(x);

				data.set_velocity(i, velocity);
			}
			log_init_condition();
		}

		void run_timestep() {
			StateData next_data = StateData(N);

			for (int i = 0; i < N; ++i) {
				StencilView stencil = data.get_stencil(i - rule.left, i, i + rule.right);
				double new_velocity = rule.apply(stencil, dx, dt);

				next_data.set_velocity(i, new_velocity);
			}

			log_velocity();
			log_error();

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
		void log_data(const std::string& filename, const std::vector<double>& data) const {
			std::ofstream file(
				simulation_name + "/" + filename + ".dat"
			);

			for (const double & value : data){
				file << value << "\n";
			}
			file.close();
		}
		
		void log_error() const {
			const std::string directory_name = "/metadata/error_L1_" + std::to_string(iter_count);
			const std::vector<double> error_vector = {error_L1()};
			log_data(directory_name, error_vector);
		}

		void log_init_condition() const {
			std::ofstream file(simulation_name + "/" + "initial_condition.dat");
	
			for (int i = 0; i < N; ++i) {
				file << data.get_velocity(i) << "\n";
			}
			file.close();
		}

		void log_velocity() const {
			std::ofstream file(
				simulation_name + "/simulation/" 
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

std::function<double(double, double)> AdvectionSolution(
		std::function<double(double)> init_condition, 
		double speed=1.0
) {
	return [init_condition, speed](double x, double t) {
		return init_condition(x - speed * t);
	};
}

////

// NONSENSE FUNCTIONS

void Run_Warming_Beam() {
	Advection_Warming_Beam iter_rule(1.0);

	std::string output_name = "Warming_Beam";
	NumericalScheme scheme = NumericalScheme(
		100,
		IndicatorInitCondition,
		iter_rule,
		output_name,
		0.005 // dx
	);

	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}

}

void Run_Lax_Friedrichs() {
	Advection_Lax_Friedrichs iter_rule(1.0);
	std::string output_name = "Lax Friedrichs";

	NumericalScheme scheme = NumericalScheme(
		100,
		IndicatorInitCondition,
		iter_rule,
		output_name,
		0.008
	);

	scheme.analytic_solution = AdvectionSolution(IndicatorInitCondition);


	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}
}



void Run_Lax_Wendroff() {
	Advection_Lax_Wendroff iter_rule(1.0);
	std::string output_name = "Lax Wendroff";

	NumericalScheme scheme = NumericalScheme(
		100,
		IndicatorInitCondition,
		iter_rule,
		output_name,
		0.01
	);

	scheme.analytic_solution = AdvectionSolution(IndicatorInitCondition);

	for (int i = 0; i<100; ++i){
		scheme.run_timestep();
	}

}

int main() {
	// Advection_ForwardDifference iter_rule(-1.0);
	
	//Diffusion iter_rule;
	Run_Lax_Friedrichs();
	Run_Lax_Wendroff();

	return 0;
};
