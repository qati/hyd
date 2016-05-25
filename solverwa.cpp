#include "solverwa.hpp"
#include <cassert>
#include <iostream>
#include <cstdlib>


using namespace std;


void fw_calc(void * context, const uint& i, const uint& j)
{
    static_cast<solverwa*>(context)->analitic(i, j);
    static_cast<solverwa*>(context)->error(i, j);

    return;
}


solverwa::solverwa(io * FS) : solver(FS)
{
    fptr.push_back(&fw_calc);

    a.resize(p.row(), p.col());

    rd.resize(4);
    ad.resize(4);

    vector<vector<dtype>*> tmp;
    tmp.push_back(&ad);
    tmp.push_back(&rd);

    fs->add_hyd(&a);
    fs->add_vec(tmp);

    cs2     = 1/fs->get("kappa");
    iv.rho0 = fs->get("rho0");
    iv.p0   = fs->get("p0");

    X       = fs->get("X0");
    Y       = fs->get("Y0");
    V       = fs->get("V0");
    W       = fs->get("W0");

    energy_s_corr = fs->get("p_s_corr");

}

solverwa::~solverwa()
{
}


void solverwa::analitic(const uint& i, const uint& j)
{
    cout<<"Calling virtual function testmod::analitic!"<<endl;
    assert(0);
    return;
}


void solverwa::error(const uint& i, const uint& j)
{
    for(ik=0;ik<4;ik++){
        rd.at(ik)     += fabs(p.at(ik, i, j) - a.at(ik, i, j));
        ad.at(ik)     += fabs(p.at(ik, i, j) - a.at(ik, i, j));
    }

    return;
}


void solverwa::start_step()
{
    cout<<"Calling virtual function testmod::start_step!"<<endl;
    assert(0);
}


void solverwa::end_step()
{
    cout<<"Calling virtual function testmod::end_step!"<<endl;
    assert(0);
}


void solverwa::init()
{
    uint j, k;

    start_step();

    p.at(0, i0, j0) = iv.rho0;
    p.at(3, i0, j0) = iv.p0/cs2;
    for(uint i=1;i<=a.row();i++){
        for(j=1;j<=a.col();j++){
            analitic(i, j);
            for(k=0;k<4;k++) p.at(k, i, j) = a.at(k, i, j);
            calc_integral(i, j);
        }
    }
    calc_eps();

    nsolver.pq();

    fs->write(t, 0);

    return;
}




/**
 * SYMSOL MODULE.
 */
symsol::symsol(io * FS) : solverwa(FS)
{
    X0 = fs->get("X0");
    cout << X0 <<endl;
}


void symsol::start_step()
{
    X = X0*sqrt(1+t*t),
    V = X0*t/sqrt(1+t*t);
    for(ik=0;ik<4;ik++){
        ad.at(ik) = 0;
        rd.at(ik) = 0;
    }

    return;
}


void symsol::end_step()
{
    for(ik=0;ik<4;ik++){
        rd.at(ik) /= a.row()*a.col();
        ad.at(ik) /= a.row()*a.col();
    }

    rd.at(0) /= iv.rho0*pow(X0, 2)*M_PI*(erf(st.x2/X)-erf(st.x1/X))*(erf(st.y2/X)-erf(st.y1/X))/4;
    rd.at(1) /= V/X*(st.y2-st.y1)*(st.x2*st.x2+st.x1*st.x1)/2;
    rd.at(2) /= V/X*(st.x2-st.x1)*(st.y2*st.y2+st.y1*st.y1)/2;
    rd.at(3) /= (iv.p0)*pow(X0, 2+2*cs2)*(M_PI/4)*(erf(st.x2/X)-erf(st.x1/X))*(erf(st.y2/X)-erf(st.y1/X))/pow(X, 2*cs2);

    return;
}


void symsol::analitic(const uint& i, const uint& j)
{
    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    a.at(0, i, j) = iv.rho0 * pow(X0/X, 2) * exp(-(x*x+y*y)/X/X);
    a.at(1, i, j) = x * V / X;
    a.at(2, i, j) = y * V / X;
    a.at(3, i, j) = (iv.p0) * pow(X0/X, 2*(1+cs2)) * exp(-energy_s_corr*(x*x+y*y)/X/X);

    return;
}




/**
 * ASYMSOL MODULE.
 */
asymsol::asymsol(io * FS) : solverwa(FS)
{
    X0 = fs->get("X0");
    Y0 = fs->get("Y0");
}


void asymsol::start_step()
{
    numeric_XY();

    for(ik=0;ik<4;ik++){
        ad.at(ik) = 0;
        rd.at(ik) = 0;
    }

    return;
}


void asymsol::end_step()
{
    for(ik=0;ik<4;ik++){
        rd.at(ik) /= a.row()*a.col();
        ad.at(ik) /= a.row()*a.col();
    }

    rd.at(0) /= iv.rho0*X0*Y0*M_PI*(erf(st.x2/X)-erf(st.x1/X))*(erf(st.y2/Y)-erf(st.y1/Y))/4;
    rd.at(1) /= V/X*(st.y2-st.y1)*(st.x2*st.x2+st.x1*st.x1)/2;
    rd.at(2) /= W/Y*(st.x2-st.x1)*(st.y2*st.y2+st.y1*st.y1)/2;
    rd.at(3) /= (iv.p0)*pow(X0*Y0, 1+cs2)*(M_PI/4)*(erf(st.x2/X)-erf(st.x1/X))*(erf(st.y2/Y)-erf(st.y1/Y))/pow(X*Y, cs2);

    return;
}


void asymsol::numeric_XY()
{
    alpha = iv.p0*pow( X0*Y0 , cs2)/iv.rho0*2;  //miert kell *4

    k1 = dt*alpha* pow(Y, -cs2) * pow(X, -cs2-1);
    l1 = dt*alpha* pow(X, -cs2) * pow(Y, -cs2-1);
    p1 = dt*V;
    q1 = dt*W;

    k2 = dt*alpha* pow(Y+q1/2, -cs2) * pow(X+p1/2, -cs2-1);
    l2 = dt*alpha* pow(X+p1/2, -cs2) * pow(Y+q1/2, -cs2-1);
    p2 = dt*(V+k1/2);
    q2 = dt*(W+l1/2);

    V += k2;
    W += l2;
    X += p2;
    Y += q2;

    return;
}


void asymsol::analitic(const uint& i, const uint& j)
{
    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    a.at(0, i, j) = iv.rho0 * (X0*Y0)/(X*Y) * exp(-pow(x, 2)/pow(X, 2) - pow(y, 2)/pow(Y, 2));
    a.at(1, i, j) = (V / X) * x;
    a.at(2, i, j) = (W / Y) * y;
    a.at(3, i, j) = (iv.p0) * pow((X0*Y0)/(X*Y), cs2 + 1) * exp(-energy_s_corr*(pow(x, 2)/pow(X, 2) + pow(y, 2)/pow(Y, 2)));

    return;
}




/**
 * New asymmetric exact solution.
 */
nassol::nassol(io * FS) : solverwa(FS)
{
    X0  = FS->get("X0");
    Y0  = FS->get("Y0");
    R2i = fs->get("R2i");
    E2  = fs->get("E2");
    E3  = fs->get("E3");
    E4  = fs->get("E4");

    if (t<0.1){
        cerr<<"Start time must be greather then 0.1!"<<endl;
        cout<<"Please enter a new t0: ";
        cin>>t;
        if (t<0.1) exit(1);
    }
}


void nassol::start_step()
{
    for(ik=0;ik<4;ik++){
        ad.at(ik) = 0;
        rd.at(ik) = 0;
    }
    return;
}


void nassol::end_step()
{
    for(ik=0;ik<4;ik++){
        rd.at(ik) /= a.row()*a.col();
        ad.at(ik) /= a.row()*a.col();
    }

    rd.at(0) /= iv.rho0*X0*Y0*M_PI*(erf(st.x2/X)-erf(st.x1/X))*(erf(st.y2/Y)-erf(st.y1/Y))/4;
    rd.at(1) /= V/X*(st.y2-st.y1)*(st.x2*st.x2+st.x1*st.x1)/2;
    rd.at(2) /= W/Y*(st.x2-st.x1)*(st.y2*st.y2+st.y1*st.y1)/2;
    rd.at(3) /= (iv.p0/cs2)*pow(X0*Y0/X/Y, 1+cs2)*(st.x2-st.x1)*(st.y2-st.y1);

    return;
}


void nassol::analitic(const uint& i, const uint& j)
{
    /**
     * If p0>rho0 then space must be: [0, x]x[0,y]
     */
    X = V*t;
    Y = W*t;

    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    if (x==0 && y==0){
        s = 0;
    } else {
        fi = atan2(y, x);
        s = pow(x*x+y*y, 1)  * pow(R2i, 1) * (1+E2*cos(2*fi));//+E4*cos(4*fi)+E3*cos(3*fi));
    }

    s=x*x/X/X+y*y/Y/Y;

    a.at(0, i, j) = iv.rho0 * (X0*Y0)/(X*Y) * exp(-s);
    a.at(1, i, j) = x/t;
    a.at(2, i, j) = y/t;
    a.at(3, i, j) = (iv.p0/cs2) * pow(X0*Y0/X/Y, 1+cs2);

    return;
}
