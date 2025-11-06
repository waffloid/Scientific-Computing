#ifndef ITERATIVE_RULE_H
#define ITERATIVE_RULE_H

#include <vector>
#include <tuple>

using StateVector = std::vector<double>;

class IterativeRule {
    public:
        int left;
        int right;

        IterativeRule(int left, int right)
            : left(left), right(right) {}

        virtual double apply(const std::vector<double>& stencil, double dx, double dt) = 0;
        virtual ~IterativeRule() = default;

};

// ADVECTION EQUATION //

class Advection_ForwardDifference : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double u = stencil[0];
            double uPlus = stencil[1];
            double C = speed * (dt / dx);

            return u - C * (uPlus - u);
        }

        Advection_ForwardDifference(double speed)
            : IterativeRule(0, 1), speed(speed) {}

};

class Advection_Lax_Friedrichs : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double u_minus = stencil[0];
            double u_plus = stencil[1];

            double C = speed * (dt / dx);

            return 0.5 * (1 + C) * u_minus + 0.5 * (1 - C) * u_plus;
        }

        Advection_Lax_Friedrichs(double speed)
            : IterativeRule(-1, 1), speed(speed) {}
};

class Advection_Lax_Wendroff : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double u_minus = stencil[0];
            double u = stencil[1];
            double u_plus = stencil[2];

            double C = speed * (dt / dx);

            return 0.5 * C * (1+C) * u_minus + (1 - C*C) * u - 0.5 * C * (1 - C) * u_plus;
        }

        Advection_Lax_Wendroff(double speed)
            : IterativeRule(-1, 1), speed(speed) {}
};

// nb, only handles positive speed!
class Advection_Warming_Beam : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double u_mm = stencil[0];
            double u_m = stencil[1];
            double u = stencil[2];

            double C = speed * (dt / dx);

            return - 0.5 * C * (1 - C) * u_mm + C * (2 - C) * u_m + 0.5 * (C-1) * (C-2) * u;
        }

        Advection_Warming_Beam(double speed)
            : IterativeRule(-2, 0), speed(speed) {}

};

// DIFFUSION EQUATION //

class Diffusion : public IterativeRule {
    public:
        double alpha = 0.3;

        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double u_left = stencil[0];
            double u = stencil[1];
            double u_right = stencil[2];

            return (alpha / 2) * (u_left + u_right) + (1 - alpha) * u;
        }

        Diffusion() : IterativeRule(-1, 1) {}

};


//////////////////////////////////////////
///

class Burger_Lax_Friedrichs : public IterativeRule {
	private:
		double get_flux(double u, double u_plus, double dx, double dt){
		return 0.25 * (u * u + u_plus * u_plus) - 0.5 * (dx / dt) * (u_plus - u);
        }

	public:
        double get_flux_balance(const std::vector<double>& stencil, double dx, double dt) {
            auto [u_minus, u, u_plus] = std::tie(stencil[0], stencil[1], stencil[2]);

            double flux_l = get_flux(u_minus, u, dx, dt);
            double flux_r = get_flux(u, u_plus, dx, dt);

            return flux_r - flux_l;
        }

		double apply(const std::vector<double>& stencil, double dx, double dt) override {
			double flux_balance = get_flux_balance(stencil, dx, dt);
            double u = stencil[1];

            return u - (dt/dx) * flux_balance;
		}

		Burger_Lax_Friedrichs() : IterativeRule(-1, 1) {}
};

class Burger_Richtimeyer : public IterativeRule {
	private:
		double analytic_flux(double u){
			return 0.5 * u * u;
		}

		double get_flux(double u, double u_plus, double dx, double dt) {
			double avg_u = 0.5 * (u + u_plus);
			double u_balanced = avg_u - 0.5 * (dt / dx) * (analytic_flux(u_plus) - analytic_flux(u));

			return analytic_flux(u_balanced);
		}

	public:
		double get_flux_balance(const std::vector<double>& stencil, double dx, double dt) {
			auto [u_minus, u, u_plus] = std::tie(stencil[0], stencil[1], stencil[2]);
		    
            double flux_r = get_flux(u, u_plus, dx, dt);
            double flux_l = get_flux(u_minus, u, dx, dt);

            return flux_r - flux_l;
        }
        
        double apply(const std::vector<double>& stencil, double dx, double dt) override {
            double flux_balance = get_flux_balance(stencil, dx, dt);
            double u = stencil[1];

			return u - (dt/dx) * flux_balance;
		}

		Burger_Richtimeyer() : IterativeRule(-1, 1) {}
};

class Burger_FORCE : public IterativeRule {
	private:
		Burger_Lax_Friedrichs LF_iterrule;
		Burger_Richtimeyer R_iterrule;
	public:
		double apply(const std::vector<double>& stencil, double dx, double dt) override {
			double flux_balance_lax_friedrichs = LF_iterrule.get_flux_balance(stencil, dx, dt);
			double flux_balance_richtimeyer = R_iterrule.get_flux_balance(stencil, dx, dt);
            double flux_balance = flux_balance_lax_friedrichs;//(flux_balance_lax_friedrichs + flux_balance_richtimeyer) / 2;
            
            double u = stencil[1];

			return u - (dt / dx) * flux_balance;
		}

		Burger_FORCE() : IterativeRule(-1, 1), LF_iterrule(Burger_Lax_Friedrichs()), R_iterrule(Burger_Richtimeyer()) {}
};


class Burger_Godunov : public IterativeRule {
	private:
		double get_u_i_half(double u, double u_plus) {
			double S = (u + u_plus) / 2;
			if (S > 0) {
				if (u > u_plus) {
					return u;
				} else {
					return u_plus;
				}
			} else {
				if ((u < 0) & (u_plus > 0)) {
					return 0;
				} else if (u < u_plus) {
					return u;
				} else {
					return u_plus;
				}
			}
		}
		double get_flux(double u, double u_plus) {
			double u_i_half = get_u_i_half(u, u_plus);
			return 0.5 * u_i_half * u_i_half;
		}
	public:
		double apply(const std::vector<double>& stencil, double dx, double dt) override {
			double u_minus = stencil[0];
			double u = stencil[1];
			double u_plus = stencil[2];

			double flux_L = get_flux(u_minus, u);
			double flux_R = get_flux(u, u_plus);

			return u - (dt/dx) * (flux_R - flux_L);
		}

		Burger_Godunov() : IterativeRule(-1, 1) {}
};


#endif 
