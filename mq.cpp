#include "mq.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


mq::mq(io * FS)
{
    fs     = FS;

    ps.pt0 = fs->get("pt1");
    ps.pt1 = fs->get("pt2");
    ps.dpt = fs->get("dpt");
    ps.fi0 = fs->get("fi1");
    ps.fi1 = fs->get("fi2");
    ps.dfi = fs->get("dfi");

    x1     = fs->get("x1");
    y1     = fs->get("y1");
    dx     = fs->get("dx");
    dy     = fs->get("dy");

    cs2    = fs->get("cs2");

    LAST_V_ID = fs->get("LAST_V_ID");

    pt.resize((ps.pt0+ps.pt1)/ps.dpt);
    vn.resize(pt.size(), LAST_V_ID);
    fi.resize((fabs(ps.fi0)+fabs(ps.fi1))/ps.dfi);
    Nfi.resize(fi.size(), 1);
    Npt.resize(pt.size(), 1);

    fs->add_m(&pt, &vn);
    fs->add_m(&fi, &Nfi);
    fs->add_m(&pt, &Npt);

    N.resize(pt.size(), fi.size());

    m = fs->get("m");
}


void mq::add_sol(hyd* P)
{
    p = P;

    r = p->row();
    c = p->col();

    return;
}


void mq::calcN(const dtype& _pt, const dtype& _fi)
{
    tN = 0;
    for(k=1;k<=r;++k){
        for(l=1;l<=c;++l){

            T = cs2*p->at(3, k, l)/p->at(0, k, l);

            pu = pow( (m*p->at(1, k, l)-_pt*cos(_fi)) , 2) + pow( (m*p->at(2, k, l)-_pt*sin(_fi)) , 2);

            tN += p->at(0, k, l)*exp(-pu/2/m/T) * dx*dy;

            if (tN!=tN || tN>1e+10){
                cout<<k<<" "<<l<<" p="<<cs2*p->at(3, k, l)<<" n="<<p->at(0, k, l)<<"; exp("<<(m*p->at(1, k, l)-_pt*cos(_fi))*(m*p->at(2, k, l)-_pt*sin(_fi))/2/m/T<<") "<<endl;

            }
        }
    }
    return;
}


void mq::calc()
{
    dtype _fi = ps.fi0,
          _pt = ps.pt0;
    uint i, j, q;

    dtype iN;

    for(i=0;i<pt.size();++i){
        iN = 0;
        _fi = ps.fi0;
        for(q=0;q<LAST_V_ID;q++){
            vn.at(i, q) = 0;
        }
        pt.at(i) = _pt;
        for(j=0;j<fi.size();++j){
            calcN(_pt, _fi);
            for(q=0;q<LAST_V_ID;q++){
                vn.at(i, q) += tN*cos((q+1)*_fi)*ps.dfi;
            }
            iN += tN*ps.dfi;
            N.at(i, j) = tN;
            _fi += ps.dfi;
        }
        Npt.at(i, 0) = iN;
        for(q=0;q<LAST_V_ID;q++){
            vn.at(i, q) /= iN;
        }
        _pt += ps.dpt;
    }

    _fi = ps.fi0;
    _pt = ps.pt0;
    for(i=0;i<fi.size();++i){
        iN  = 0;
        _pt = ps.pt0;
        fi.at(i)=_fi;
        for(j=0;j<pt.size();++j){
            iN  += N.at(j, i)*ps.dpt;
            _pt += ps.dpt;
        }
        fi.at(i)  = _fi;
        Nfi.at(i, 0) = iN;
        _fi += ps.dfi;
    }
    fs->write();
    return;
}
