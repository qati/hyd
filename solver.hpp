#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_


#include "mff.hpp"
#include "io.hpp"


class solver{
protected:
    /**
     * Numeric solver.
     */
    mff nsolver;

    /**
     * Solution container.
     */
    hyd p;

    /**
     * IO module.
     */
    io * fs;

    /**
     * Speed of sound.
     */
    dtype kappa;

    /**
     * Mass: rho = n * m.
     */
    dtype m;

    /**
     * eps parameters.
     */
    std::vector<dtype> eps;
    uint LAST_V_ID;


    /**
     * Integrals to calculate eps parameters.
     */
    std::vector<dtype> eint, pint, lint;

    /**
     * Integral radius.
     */
    dtype ir;

    /**
     * x, y point.
     */
    dtype x, y, r, fi;

    /**
     * Space.
     */
    struct st_t {dtype x1, x2, y1, y2, dx, dy;} st;

    /**
     * Index of grid point where rho/eps is max.
     */
    uint i0, j0;

    /**
     * Time resolution.
     */
    dtype dt;

    /**
     * Time.
     */
    dtype t;

    /**
     * Calculate eps3 parameter.
     */
    void calc_integral(const uint&, const uint&);
    void null_integral();
    void calc_eps();

    /**
     * Start/end step.
     */
    virtual void start_step();
    virtual void end_step();

    /**
     * Functions.
     */
    std::vector<void (*)(void*, const uint&, const uint&)> fptr;

    /**
     * Initial distributions
     */
    void gauss();
    void woods_saxon();

    /**
     * Hadronized cells.
     */
    hyd sol;

    /**
     * Hadronization temperature.
     */
    dtype T0, T;


    public:
        solver(io*);
        virtual ~solver();

        /**
         * Init function.
         * Sets up initial value.
         */
        virtual void init();

        /**
         * Loop function.
         */
        void loop();

        /**
         * Get solution.
         */
        hyd* get_sol();

        /**
         * Forwarder functions.
         */
        friend void fw_calc_integral(void*, const uint&, const uint&);
};


#endif // _SOLVER_HPP_

