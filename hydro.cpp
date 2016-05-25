#include "hydro.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>

using namespace std;


hydro::hydro()
{
    p   = NULL;

    kappa = 0;
    u1  = 0;
    u2  = 0;

    eos = new EoS(&kappa);
}

hydro::hydro(const dtype& DT, const dtype& DX, const dtype& DY, hyd* h)
{
    set(DT, DX, DY, h);

    p   = NULL;

    kappa = 0;
    u1  = 0;
    u2  = 0;

    eos = new EoS(&kappa);
}


void hydro::set(const dtype& DT, const dtype& DX, const dtype& DY, hyd* h)
{
    p  = h;

    dt = DT;
    dx = DX;
    dy = DY;

    q.resize(p->row(), p->col());

    r = p->row();
    c = p->col();

    return;
}


void hydro::set(const dtype& KAPPA)
{
    if (KAPPA <= 0) eos->set_mode(EOS_QCD_KAPPA);
    else eos->set_mode(EOS_CONST_KAPPA);
    kappa = KAPPA;
    return;
}


void hydro::set(const dtype& U1, const dtype& U2)
{
    u1 = U1;
    u2 = U2;
    return;
}


void hydro::check()
{
    if (p==NULL){
        cerr<<"hydro::check FATAL ERROR!"<<endl<<"No pointer to phys variables set!"<<endl;
        exit(1);
    }
    if (kappa==0){
        cerr<<"hydro::check FATAL ERROR!"<<endl<<"Not set speed of sound!"<<endl;
        exit(1);
    }
    return;
}


void hydro::pq()
{
    check();
    uint i, j;
    try {
    for(i=1;i<=r;i++){
        for(j=1;j<=c;j++){
            //gamma = (dtype)sqrt(1+pow(p->at(1, i, j),2)+pow(p->at(2, i, j),2));

            q.at(0, i, j) = (dtype)(p->at(0, i, j));
            q.at(1, i, j) = (dtype)(p->at(0, i, j) * p->at(1, i, j)) ;
            q.at(2, i, j) = (dtype)(p->at(0, i, j) * p->at(2, i, j)) ;
            //T = p->at(3, i, j)/p->at(0, i, j);
            //eos->ck(T);
            q.at(3, i, j) = (dtype)(p->at(3, i, j)*kappa - p->at(0,i,j) + p->at(0, i, j) * (pow(p->at(1,i,j),2)+pow(p->at(2,i,j),2))/2);
        }
    }
    } catch(logic_error& e){
        string w(e.what());
        if (w=="nan") non++;
        //else if (w=="0") cerr<<"0 in rho!"<<endl;
    }

    return;
}

void hydro::qp(const uint& i, const uint& j)
{
    try {
        //if (q.at(0, i, j)==0)             throw logic_error("0");
        if (q.at(0, i, j)!=q.at(0, i, j)) throw logic_error("nan");
        if (q.at(0, i, j)<DEF_NULL || q.at(0,i,j)==0.0){
            p->at(0, i, j) = 0.0;
            p->at(1, i, j) = 0.0;
            p->at(2, i, j) = 0.0;
            p->at(3, i, j) = 0.0;
        }
        //gamma = (dtype)sqrt(1+pow(q.at(1, i, j)/q.at(0, i, j),2)+pow(q.at(2, i, j)/q.at(0, i, j),2));
        //if (gamma<DEF_NULL){
         //   gamma = 1;
       // }
        p->at(0, i, j) = (dtype)(q.at(0, i, j));
        p->at(1, i, j) = (dtype)(q.at(1, i, j)   / q.at(0, i, j));
        p->at(2, i, j) = (dtype)(q.at(2, i, j)  / q.at(0, i, j));
        p->at(3, i, j) = (dtype)((q.at(3, i, j) + q.at(0,i,j) - (pow(q.at(1, i, j), 2)+pow(q.at(2, i, j), 2))/q.at(0, i, j)/2)/kappa);
        
       if (p->at(0,i,j)!=p->at(0,i,j)) throw logic_error("nan");
       if (p->at(0,i,j)<DEF_NULL || p->at(3,i, j)<DEF_NULL){
        p->at(0,i,j)=0.0;
        p->at(1,i,j)=0.0;
        p->at(2,i,j)=0.0;
        p->at(3,i,j)=0.0;

       }
        // T = p->at(3, i, j)/p->at(0, i, j);
       // eos->ck(T);

    } catch(logic_error& e){
        string w(e.what());
        if (w=="nan") non++;
        //else if (w=="0") cerr<<"0 in rho!"<<endl;
    }
    return;
}


void hydro::flux_f(vector<dtype>& f, const dtype& q1, const dtype& q2, const dtype& q3, const dtype& q4)
{
    //T = (q4-(q2*q2+q3*q3)/q1/2)/q1/kappa;
   // eos->ck(T);
    if (q1==0.0 || q1<DEF_NULL){
        f[0]=0;
        f[1]=0;
        f[2]=0;
        f[3]=0;
        return;
    }
       
     f.at(0) = (dtype)(q2);
   // f.at(1) = (dtype)(pow(q2, 2)/q1 + (  q4 - (pow(q2,2)+pow(q3,2))/q1/2  )*kappa +q1*kappa);
    f.at(1) = (dtype)(pow(q2, 2)/q1 + (  q4 + q1 - (pow(q2,2)+pow(q3,2))/q1/2  )/kappa);   // - (q2+q3)*(q2/q1)/q1);
    f.at(2) = (dtype)(q2*q3/q1);// - (q2+q3)*(q2/q1)/q1/gamma);
   // f.at(3) = (dtype)((  (1+kappa)*q4 - kappa*( pow(q2,2)+pow(q3,2) )/q1/2 + q1/kappa )*q2/q1);
    f.at(3) = (dtype)((  (1+kappa)*(  q4 - (pow(q2,2)+pow(q3,2))/q1/2  )/kappa + (pow(q2,2)+pow(q3,2))/q1/2 + q1/kappa )*q2/q1);
    //if (f[0]!=f[0] || f[1]!=f[1] || f[2]!=f[2] || f[3]!=f[3]){
    //    cout<<"flux_f: " << q1<<" "<<q2<<" " <<q3 <<" "<<q4<<endl;
   // }
    return;
}

void hydro::flux_g(vector<dtype>& f, const dtype& q1, const dtype& q2, const dtype& q3, const dtype& q4)
{
   // T = (q4-(q2*q2+q3*q3)/q1/2)/q1/kappa;
   if (q1==0.0 || q1<DEF_NULL){
       f[0]=0;
       f[1]=0;
       f[2]=0;
       f[3]=0;
       return;
   }
   // eos->ck(T);
     //gamma = (dtype)sqrt(1+pow(q2/q1,2)+pow(q3/q1,2));
     //if (gamma<DEF_NULL) gamma=1;
    f.at(0) = (dtype)(q3);
    f.at(1) = (dtype)(q2*q3/q1);//-(q2+q3)*(q2/q1)/q1/gamma);
    //f.at(2) = (dtype)(pow(q3, 2)/q1 + (  q4 - (pow(q2,2)+pow(q3,2))/q1/2  )*kappa + q1*kappa);
    f.at(2) = (dtype)(pow(q3, 2)/q1 + (  q4 + q1 - (pow(q2,2)+pow(q3,2))/q1/2  )/kappa); //- (q2+q3)*(q3/q1)/g1/gamma);
    //f.at(3) = (dtype)((  (1+kappa)*q4 - kappa*( pow(q2,2)+pow(q3,2) )/q1/2  +q1*kappa)*q3/q1);
    f.at(3) = (dtype)((  (1+kappa)*(  q4 - (pow(q2,2)+pow(q3,2))/q1/2  )/kappa + (pow(q2,2)+pow(q3,2))/q1/2 + q1/kappa)*q3/q1);

    return;
}


void hydro::vflux_f(std::vector<dtype>& vf, const dtype& s, const dtype& u, const dtype& v,
    const dtype& dudx, const dtype& dvdx, const dtype& dudy, const dtype& dvdy)
{
    //eos->ck(T);
    visc1 = u1 * s;
    visc2 = (u2-2*u1/3) * s * (1+1/kappa);

    vf.at(0) = 0;
    vf.at(1) = - (2*visc1*dudx + visc2*(dudx+dvdy));
    vf.at(2) = -(visc1*(dvdx+dudy));
    vf.at(3) = -(vf.at(1)*u+vf.at(2)*v);

    return;
}

void hydro::vflux_g(std::vector<dtype>& vf, const dtype& s, const dtype& u, const dtype& v,
    const dtype& dudx, const dtype& dvdx, const dtype& dudy, const dtype& dvdy)
{
    //eos->ck(T);
    visc1 = u1 * s;
    visc2 = (u2-2*u1/3) * s * (1+1/kappa);

    vf.at(0) = 0;
    vf.at(1) = -(visc1*(dvdx+dudy));
    vf.at(2) = -(2*visc1*dvdy + visc2*(dudx+dvdy));
    vf.at(3) = -(vf.at(1)*u+vf.at(2)*v);

    return;
}



uint& hydro::get_non()
{
    return non;
}

uint hydro::get_nop()
{
    return r*c;
}
