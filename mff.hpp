#ifndef _MFF_HPP_
#define _MFF_HPP_


#include "hydro.hpp"


#define MUSTA_K 12

const dtype CFLN      = 0.4 ;
const dtype CFLN_WV   = 0.9 ;
const dtype CFLN_VISC = 1.5 ;

const bool FLUX_X = 1 ;
const bool FLUX_Y = 0 ;


class mff : public hydro{
    /**
     * f_{i+1/2} flux.
     */
	//hyd fi_x, fi_y, vfi_x, vfi_y;

    hyd k1, k2, k3, k1y, k2y, k3y;
    dtype k4, k4y;


    std::vector<dtype> flx, flx2, flx_y, flx2_y;

    /**
     * Matrix of speed,entropy and energy fields.
     */
    dtype u, v, s, e; //!!!!!!!!!!!!!!! s is not the entropy!!! s = rho!!!!!

    /**
     * Derivates of v.
     */
    std::vector<dtype> v_deriv;

    /**
     * TMP variables for MUSTA method.
     */
    std::vector<dtype> ql, qm, qr, fl, fm, fr, QL, QR;


    /**
     * Calculate fi with MUSTA FORCE K method.
     */
    void iflux(void (hydro::*)(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&),bool);


    void CalcThings();

    /**
     * Function forwarder and obj.
     */
    std::vector<void*> call_context;
    std::vector<void (*)(void*, const uint&, const uint&)> call_func;

    /**
     * Boundary condition.
     */
    void boundary(const uint&, const uint&);

    /**
     * Find maximums for stability criterions.
     */
    dtype vmax, cmax, findcmax_c;
    dtype Dmax;
    bool mdirx;

    void init_max();
    void find_max(const uint&, const uint&);

    /**
     * Change dt to satisfy stability criterions.
     */
    void update_dt();

    public:
        mff();
        mff(const dtype&, const dtype&,const dtype&, hyd*);

        /**
         * Set variables.
         */
        void set(const dtype& kappa){hydro::set(kappa);}
        void set(const dtype& U1, const dtype& U2){hydro::set(U1, U2);}
        void set(const dtype&, const dtype&,const dtype&, hyd*);

        void add_callfunc(void*, void (*)(void*, const uint&, const uint&));

        void step();

        dtype get_flux(bool, const uint&, const uint&, const uint&){return 1.0;};
};


#endif // _MFF_HPP_
