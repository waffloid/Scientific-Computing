#include<vector>
#include<iostream>
#include<cmath>

using StateVector = std::vector<double>;

class NumericalScheme {
    public:
        StateVector initialConditions;
        std::vector<StateVector> velocityLog;

        double t = 0;
        double dx = 0.01;
        double dt = 0.01;
        int N;
        int index = 0;

        NumericalScheme(double dx = 0.01, double dt = 0.01, int N = 100) 
            : dx(dx), dt(dt), N(N) {}

        void setInitConditions(StateVector initialConditions) {
            if(index != 0) {
                std::cout << "Bad! Only initialize ICs before running simulation!" << std::endl;
            }

            velocityLog.emplace_back(initialConditions);
        }
};


class IterativeRule {
    public:
        int stencilLeft;
        int stencilRight;
        NumericalScheme *parent;

        IterativeRule(int left, int right, NumericalScheme *scheme)
            : stencilLeft(left), stencilRight(right), parent(scheme) {}

        virtual double apply(const StateVector& data) = 0;
        virtual ~IterativeRule() = default;

};

class Advection_ForwardDifference : public IterativeRule {
    public:
        double speed = 1;


        double apply(const StateVector& data) {
            double u = data.at(0);
            double uPlus = data.at(1);
            double dx = parent->dx;
            double dt = parent->dt;
            double C = speed * (dt / dx);

            return u - C * (uPlus - u);
        }

        Advection_ForwardDifference(NumericalScheme *scheme, double speed)
            : IterativeRule(0, 1, scheme), speed(speed) {}

};

void ApplySineInitCondition(NumericalScheme& scheme) {
    StateVector initialConditions;

    for (int i = 0; i < scheme.N; i++) {
        double value = std::sin(2 * M_PI * i / scheme.N);
        initialConditions.emplace_back(value);
    }

    scheme.setInitConditions(initialConditions);
}


int main() {
    NumericalScheme scheme = NumericalScheme();
    ApplySineInitCondition(scheme);

    std::cout << "Successful!" << std::endl;

    return 0;
}

/*


Things to be able to do:

* 1D numerical simulations
    * u -> u data, (1d for now, but future should be able to template higher dim u), 1 spatial dim
    * flux data! conservative schemes.
    * able to pass in different "solvers", via enum or via fn.
    * methods:
        * batch experiment -- pass in a family of simulations.
            maybe have an experiment class?

design goals:
    * framework to easily express different schemes -- conservative, non-conservative, half-step, full-step.
        * different boundary conditions
    * easy to get graphics & output. in  particular,
        * variations in input:
            * choice of scheme,
            * parameters:
                * N
                * speed of sound
                * ...
        * variations in output:
            * (final) error
                * L1? ...
            * solution indexed by time
            * final solution

*/

