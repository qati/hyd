#ifndef _EOS_HPP_
#define _EOS_HPP_


#include "matrix.hpp"

#define EOS_CONST_KAPPA 0
#define EOS_QCD_KAPPA 1

const double  h0 = 0.1396;
const double  h1 = -0.18;
const double  h2 = 0.035;
const double  f0 = 2.76;
const double  f1 = 6.79;
const double  f2 = -5.29;
const double  g1 = -0.47;
const double  g2 = 1.04;

const double dT  = 1;

class EoS{
private:
    dtype * kappa;
    double t, I, T, p;
    std::vector<dtype> v_kappa;
    uint index;
    dtype x1,x2,y1,y2, m;
    bool use_QCD;
public:
    EoS(dtype*);
    void ck(const dtype&);
    void set_mode(bool);
};





#endif
