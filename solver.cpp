#include "solver.hpp"
#include <iostream>
#include <cassert>
#include <cstdlib>


using namespace std;


void fw_calc_integral(void* context, const uint& i, const uint& j)
{
    static_cast<solver*>(context)->calc_integral(i, j);
    return;
}


solver::solver(io * FS)
{
    st.x1 = FS->get("x1");
    st.x2 = FS->get("x2");
    st.y1 = FS->get("y1");
    st.y2 = FS->get("y2");
    st.dx = FS->get("dx");
    st.dy = FS->get("dy");

    kappa   = FS->get("kappa");

    m     = FS->get("m");

    p.resize( (fabs(st.x1)+fabs(st.x2))/st.dx+1, (fabs(st.y1)+fabs(st.y2))/st.dy+1 );

    nsolver.set(FS->get("dt"), FS->get("dx"), FS->get("dy"), &p);
    nsolver.set(FS->get("kappa"));
    nsolver.set(FS->get("visc_u1"), FS->get("visc_u2"));

    LAST_V_ID = FS->get("LAST_V_ID");
    eps.resize(LAST_V_ID*3);
    eint.resize((LAST_V_ID+1)*3);
    pint.resize(4);
    lint.resize(4);

    ir = FS->get("int_rad");

    fs = FS;
    fs->set_sgrid(p.row(), p.col());
    fs->add_hyd(&p);
    fs->add_vec(vector<vector<dtype>*>(1, &eps));
    vector<vector<dtype>*> tmp(2);
    tmp.at(0) = &pint;
    tmp.at(1) = &lint;
    fs->add_vec(tmp);

    i0 = -st.x1/st.dx+1,
    j0 = -st.y1/st.dy+1;

    dt = fs->get("dt");
    t  = fs->get("t0");

    fptr.push_back(&fw_calc_integral);

    T0 = fs->get("T0");

    sol.resize(p.row(), p.col());

}

solver::~solver()
{
    fs = NULL;
}


void solver::gauss()
{
    dtype kappa  = fs->get("kappa"),
          rho0 = fs->get("rho0"),
          p0   = fs->get("p0");

    dtype t0 = fs->get("t0");

    dtype R  = fs->get("R"),
          E2 = fs->get("E2"),
          E3 = fs->get("E3"),
          E4 = fs->get("E4");

    dtype psi2 = fs->get("psi2"),
          psi3 = fs->get("psi3"),
          psi4 = fs->get("psi4");

    dtype p_s_corr = fs->get("p_s_corr");

    dtype s;

    uint i, j;
    uint r = p.row(),
         c = p.col();

    p.at(0, i0, j0) = rho0;
    p.at(3, i0, j0) = p0;
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            x = st.x1 + (i-1)*st.dx;
            y = st.y1 + (j-1)*st.dy;

            if (x==0 && y==0){
                s = 0;
            } else {
                fi = atan2(y, x);
                s = (x*x+y*y) * pow(R, -2) * (  1+E2*cos( 2*(fi-psi2) )+E4*cos( 4*(fi-psi4) )+E3*cos( 3*(fi-psi3) )  );
            }
            //if (s>-R*0.25 && s<R*0.25){
                p.at(0, i, j) = rho0 * exp(-s/2);
                p.at(1, i, j) = x/t0;
                p.at(2, i, j) = y/t0;
                p.at(3, i, j) = p0 * exp(-p_s_corr*s/2);
            /*} else {
                p.at(0,i,j) = 0.0;
                p.at(1,i,j) = 0.0;
                p.at(2,i,j) = 0.0;
                p.at(3,i,j) = 0.0;
            }*/
            calc_integral(i, j);
        }
    }
    calc_eps();

    return;
}


void solver::woods_saxon()
{
     dtype kappa  = fs->get("kappa"),
           rho0 = fs->get("rho0"),
           p0   = fs->get("p0");

    dtype E2  = fs->get("E2"),
          E3  = fs->get("E3"),
          E4  = fs->get("E4");

    dtype R = fs->get("WS_R"),
          a = fs->get("WS_a");

    dtype s;

    uint i, j;
    uint r = p.row(),
         c = p.col();

    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            x = st.x1 + (i-1)*st.dx;
            y = st.y1 + (j-1)*st.dy;

            if (x==0 && y==0){
                s = 0;
            } else {
                fi = atan2(y, x);
                s = sqrt(x*x+y*y)*(1+E2*cos(2*fi)+E3*cos(3*fi)+E4*cos(4*fi));
            }
            p.at(0, i, j) =  rho0 / (1+exp((s-R)/a));
            p.at(1, i, j) = 0;
            p.at(2, i, j) = 0;
            p.at(3, i, j) =  (p0*kappa) / (1+exp((s-R)/a));
            calc_integral(i, j);
        }
    }
    calc_eps();

    return;
}


void solver::init()
{
    gauss();

    fs->write(-1);
    
    string call = "python make_smth.py "+fs->get_dfn()+" "+to_string(t);
    cout << "Call: "<<call<<endl;
    //system(call.c_str());
    fs->read_datfile_to_p(fs->get_dfn()); 
    cout << "End of data smoothing "<<endl;
    nsolver.pq();
    
    fs->write(0, t);

    return;
}


void solver::calc_integral(const uint& i, const uint& j)
{
    uint k;

    x = st.x1 + (i-1)*st.dx;
    y = st.y1 + (j-1)*st.dy;

    pint.at(0) += p.at(0, i, j)                                                                 * st.dx*st.dy;
    pint.at(1) += p.at(0, i, j)*p.at(1, i, j)                                                   * st.dx*st.dy;
    pint.at(2) += p.at(0, i, j)*p.at(2, i, j)                                                   * st.dx*st.dy;
    pint.at(3) += (p.at(3, i, j)+p.at(0, i, j)*(pow(p.at(1, i, j), 2)+pow(p.at(2, i, j), 2))/2) * st.dx*st.dy;

    if (i==1 && (j>1 && j<p.col())){
        for (k=0;k<4;k++){
            lint.at(k) -= (nsolver.get_flux(FLUX_X, k, p.row(), j) - nsolver.get_flux(FLUX_X, k, 0, j)) * st.dy * dt;
        }

    }
    if (j==1){
        for (k=0;k<4;k++){
            lint.at(k) -= (nsolver.get_flux(FLUX_Y, k, i, p.col()) - nsolver.get_flux(FLUX_Y, k, i, 0)) * st.dx * dt;
        }
    }


    if ((x*x+y*y)>ir*(st.x2*st.x2+st.y2*st.y2)) return;

    if (x==0 && y==0){
        for(vector<dtype>::iterator it=eint.begin();it!=eint.end();++it){
            (*it) -= 1;
        }
    } else {
        fi = atan2(y, x);
        for(k=0;k<=LAST_V_ID;k++){
            eint.at(k)               += p.at(0, i, j)                                       * cos(k*fi);
            eint.at(LAST_V_ID+1+k)   += exp(-sqrt(1+pow(p.at(1, i, j), 2)+pow(p.at(2, i, j), 2))) * cos(k*fi);
            eint.at(2*LAST_V_ID+2+k) += (p.at(3, i, j)+p.at(0,i,j))                                       * cos(k*fi);
        }
    }
    return;
}


void solver::null_integral()
{
    eint.clear();
    eint.resize((LAST_V_ID+1)*3);
    pint.clear();
    pint.resize(4);
    return;
}


void solver::calc_eps()
{
    for(uint i=0;i<LAST_V_ID;i++){
        if (eint.at(i+1)==eint.at(0) && eint.at(i+1)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+1)/eint.at(0);
    }
    for(uint i=LAST_V_ID;i<2*LAST_V_ID;i++){
        if (eint.at(i+2)==eint.at(LAST_V_ID+1) && eint.at(i+2)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+2)/eint.at(LAST_V_ID+1);
    }
    for(uint i=2*LAST_V_ID;i<3*LAST_V_ID;i++){
        if (eint.at(i+3)==eint.at(2*LAST_V_ID+2) && eint.at(i+3)==0) eps.at(i) = 0;
        else eps.at(i) = eint.at(i+3)/eint.at(2*LAST_V_ID+2);

    }

    return;
}


void solver::start_step()
{
    return;
}


void solver::end_step()
{
    return;
}


void solver::loop()
{
    uint N = 0;
    uint non;
    dtype tmax = t + fs->get("RTT");
    dtype tws  = (dtype)fs->get("RTT")/fs->get("NOPTBW"),
          ct   = 0;
    uint ws = fs->get("WS");

    for(vector<void (*)(void*, const uint&, const uint&)>::iterator it=fptr.begin();it!=fptr.end();++it){
        nsolver.add_callfunc(this, *it);
    }
    cout <<"STARTING LOOP, WS="<<ws<<endl;
    for(uint i=1;t<tmax;i++){
        dt = nsolver.get_dt();
        t += dt;

        null_integral();
        start_step();
        nsolver.step();
        end_step();

        ct += dt;
        if (!(i%ws)){
            ct = 0;
            calc_eps();
            fs->write(t);
            fs->write_nos(N++);
            fs->write_dt(dt);
            non = nsolver.get_non();
            cout<<"param: dt="<<dt<<", p0="<<fs->get("p0")<<", pc="<<fs->get("p_s_corr")
                   <<", visc1="<<fs->get("visc_u1")<<"; "<<"time: "<<t<<"fm/c; nop: "<<nsolver.get_nop()<<"; non: "<<non<<endl;
            non = 0;
        }

    }
    for(uint k=1;k<=p.row();++k){
        for(uint l=1;l<=p.col();++l){
            sol.at(0, k, l) = p.at(0, k, l);
            sol.at(1, k, l) = p.at(1, k, l);
            sol.at(2, k, l) = p.at(2, k, l);
            sol.at(3, k, l) = p.at(3, k, l);
        }
    }
    cout<<"HADRONIZED"<<endl;
    return;
}

hyd* solver::get_sol()
{
    return &sol;
}
