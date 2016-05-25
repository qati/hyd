#ifndef _SOLVERWA_HPP_
#define _SOLVERWA_HPP_


#include "solver.hpp"


class solverwa : public solver{
    protected:
        /**
         * @a: analitic solution
         */
        hyd a;

        /**
         * Relative and absolute difference.
         */
        std::vector<dtype> rd, ad;

        /**
         * @XYVW: scale variables.
         */
        dtype X, Y, V, W;

        /**
         * Speed of sound.
         */
        dtype cs2;

        /**
         * Init values.
         */
        struct iv_t {dtype rho0, p0;} iv;

        /**
         * Analitic solution.
         */
        virtual void analitic(const uint&, const uint&);

        /**
         * Calulate errors.
         */
        void error(const uint&, const uint&);

        /**
         * Start/End time step calc function.
         */
        virtual void start_step();
        virtual void end_step();

        /**
         * Loop index.
         */
        uint ik;

        /**
         * Energy correction.
         */
        dtype energy_s_corr;

    public:
        solverwa(io*);
        virtual ~solverwa() = 0;

        /**
         * Initialization function.
         */
        void init();

         /**
         * Forwarder functions.
         */
        friend void fw_calc(void *, const uint&, const uint&);
};




/**
 * Symmetric exact solution.
 * @rho = rho0*(X0/X)^2*exp(-(x^2+y^2)/X^2)
 * @v   = dX/dt * (x,y)
 * @eps = eps0*(X0/X)^(2+2cs2)*exp(-(x^2+y^2)/X^2)
 * @X   = X0*sqrt(1+t^2)
 * This module contain shocktube test! (initvalue_type1)
 */
class symsol : public solverwa{
    dtype X0;

    void start_step();
    void end_step();
    void analitic(const uint&, const uint&);

    public:
        symsol(io*);
};




/**
 * Asymmetric exact solution.
 * @rho = rho0*(V0/V)*exp(-s)
 * @v   = (X'*x/X,Y'*y/Y)
 * @eps = eps0*(V0/V)^(1+cs2)*exp(-s)
 * @V   = X*Y
 * @s   = x^2/X^2+y^2/Y^2
 * @XY  : X''X=Y''Y=(T0/m)*(V0/V)^(1/kappa)
 */
class asymsol : public solverwa{
    dtype X0, Y0;
    dtype k1, k2, l1, l2, p1, p2, q1, q2;
    dtype alpha;

    void start_step();
    void end_step();
    void numeric_XY();
    void analitic(const uint&, const uint&);

    public:
        asymsol(io*);
};




/**
 * New asymmetric exact solution with const preasure.
 * @rho = rho0*(V0/V)*exp(-s)
 * @v   = (V0/X * x, W0/Y * y)
 * @p   = p0*(V0/V)^(1+cs2)
 * @s   = r^2/R^2 * (1+eps_N*cos(N*fi))
 * @XY  : X=V0*t, Y=W0*t
 */
class nassol : public solverwa {
    dtype X0, Y0, s;
    dtype R2i, E2, E3, E4;

    void start_step();
    void end_step();
    void analitic(const uint&, const uint&);

    public:
        nassol(io*);
};

#endif // _SOLVERWA_HPP_
