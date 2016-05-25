#ifndef _HYDRO_HPP_
#define _HYDRO_HPP_


#include "matrix.hpp"
#include "eos.hpp"
#include <cmath>


class hydro{
    protected:
        /**
         * Spacetime resolution.
         */
        dtype dt, dx, dy;

        /**
         * Square of speed of sound.
         */
        dtype kappa;

        /**
         * Viscosity parameters.
         */
        dtype u1, u2;

        /**
         * Kinematic viscosity.
         */
        dtype visc1, visc2;

        /**
         * @p: physical quantity: n, vx, vy, eps! eps contains M! (eps'+rho*v^2/2) => (eps+n*v^2/2)
         * @q: transformed quantity
         */
        hyd *p, q;

        dtype gamma;
        
        /**
         * Temperature.
         */
        dtype T;

        /**
         * EoS
         */
        EoS *eos;

        /**
         * Check settings.
         */
        void check();

        /**
         * Calculate fluxes from q1, q2, q3, q4.
         */
        void flux_f(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&);
        void flux_g(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&);

        /**
         * Calculate viscosus fluxes from v, derivates of v and entropy.
         */
        void vflux_f(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&);
        void vflux_g(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&, const dtype&);

        uint r, c;

        /**
         * Boundary condition.
         */
        virtual void boundary(const uint&, const uint&) = 0;

        /**
         * Set variables.
         */
        void set(const dtype&, const dtype&,const dtype&, hyd*);

        /**
         * Number of nan in rho.
         * Number of zero in rho.
         */
         uint non;

    public:
        hydro();
        hydro(const dtype&, const dtype&,const dtype&, hyd*);

        /**
         * Set variables.
         */
        void set(const dtype&);
        void set(const dtype&, const dtype&);

         /**
         * Makes p to q transformation.
         */
        void pq();
        /**
         * Makes q to p transformation.
         */
        void qp(const uint&, const uint&);

        dtype& get_dt(){return dt;}

        virtual void step() = 0;

        /**
         * Get non.
         * Get noz.
         * Get num of points.
         */
        uint& get_non();
        uint get_noz();
        uint get_nop();
};

#endif // _HYDRO_HPP_
