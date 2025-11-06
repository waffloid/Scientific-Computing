#ifndef ITERATIVE_RULE_H
#define ITERATIVE_RULE_H

using StateVector = std::vector<double>;

class StencilView {
    const StateVector& velocity;
    int left;
    int center;
    int right;
    int N;  

    public:
        StencilView(const StateVector& velocity, int l, int c, int r, int N)
            : velocity(velocity), left(l), center(c), right(r), N(N) {}

    double get_velocity(int i) const {
        int idx = (((center + i) % N) + N) % N;

        return velocity[idx];
    }
};

class IterativeRule {
    public:
        int left;
        int right;

        IterativeRule(int left, int right)
            : left(left), right(right) {}

        virtual double apply(const StencilView& stencil, double dx, double dt) = 0;
        virtual ~IterativeRule() = default;

};

// ADVECTION EQUATION //

class Advection_ForwardDifference : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const StencilView& stencil, double dx, double dt) override {
            double u = stencil.get_velocity(0);
            double uPlus = stencil.get_velocity(1);
            double C = speed * (dt / dx);

            return u - C * (uPlus - u);
        }

        Advection_ForwardDifference(double speed)
            : IterativeRule(0, 1), speed(speed) {}

};

class Advection_Lax_Friedrichs : public IterativeRule {
    public:
        double speed = 1.0;
        
        double apply(const StencilView& stencil, double dx, double dt) override {
            double u_minus = stencil.get_velocity(-1);
            double u_plus = stencil.get_velocity(1);

            double C = speed * (dt / dx);

            return 0.5 * (1 + C) * u_minus + 0.5 * (1 - C) * u_plus;
        }

        Advection_Lax_Friedrichs(double speed)
            : IterativeRule(-1, 1), speed(speed) {}
};

class Advection_Lax_Wendroff : public IterativeRule {
    public:
        double speed = 1.0;

        double apply(const StencilView& stencil, double dx, double dt) override {
            double u_minus = stencil.get_velocity(-1);
            double u = stencil.get_velocity(0);
            double u_plus = stencil.get_velocity(1);

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

        double apply(const StencilView& stencil, double dx, double dt) override {
            double u_mm = stencil.get_velocity(-2);
            double u_m = stencil.get_velocity(-1);
            double u = stencil.get_velocity(0);

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

        double apply(const StencilView& stencil, double dx, double dt) override {
            double u_left = stencil.get_velocity(-1);
            double u = stencil.get_velocity(0);
            double u_right = stencil.get_velocity(1);

            return (alpha / 2) * (u_left + u_right) + (1 - alpha) * u;
        }

        Diffusion() : IterativeRule(-1, 1) {}

};

#endif 
