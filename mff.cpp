#include "mff.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;


mff::mff()
{
}

mff::mff(const dtype& DT, const dtype& DX,const dtype& DY, hyd* P)
{
    mff::set(DT, DX, DY, P);
}


void mff::set(const dtype& DT, const dtype& DX,const dtype& DY, hyd* P)
{
    hydro::set(DT, DX, DY, P);

    flx.resize(4);
    flx2.resize(4);
    flx_y.resize(4);
    flx2_y.resize(4);

    k1.resize(P->row(), P->col());
    k2.resize(P->row(), P->col());
    k3.resize(P->row(), P->col());
    k1y.resize(P->row(), P->col());
    k2y.resize(P->row(), P->col());
    k3y.resize(P->row(), P->col());

    v_deriv.resize(4);
    ql.resize(4);
    QL.resize(4);
    qm.resize(4);
    qr.resize(4);
    QR.resize(4);
    fl.resize(4);
    fm.resize(4);
    fr.resize(4);


    return;
}


void mff::boundary(const uint& i, const uint& j)
{
    for(uint k=0;k<4;k++){
        q.at(k, i, 0)   = q.at(k, i, 1);
        q.at(k, i, c+1) = q.at(k, i, c);
        q.at(k, 0, j)   = q.at(k, 1, j);
        q.at(k, r+1, j) = q.at(k, r, j);

        if (i==0 && j==0) q.at(k, 0, 0)     = q.at(k, 1, 1);
        if (i==0)         q.at(k, 0, c+1)   = q.at(k, 1, c);
        if (j==0)         q.at(k, r+1, 0)   = q.at(k, r, 1);
        if (j==c)         q.at(k, r+1, c+1) = q.at(k, r, c);
    }
    return;
}


void mff::add_callfunc(void *context, void(*func)(void*, const uint&, const uint&))
{
    call_context.push_back(context);
    call_func.push_back(func);

    return;
}


void mff::iflux(void (hydro::*flux)(std::vector<dtype>&, const dtype&, const dtype&, const dtype&, const dtype&), bool hs = false)
{
  //  uint i, j, k, l, ri, rj;
	uint k,l;
  //  uint v_deriv_ind;
    dtype ds, DT;
  //  dtype tmp;
  //  bool dir_x;

    if (hs) DT = dt*0.5;
    else DT = dt;

    if (flux==&mff::flux_f){
        ds          = dx;
       // dir_x       = true;
       // v_deriv_ind = 0;
    } else {
        ds          = dy;
     //   dir_x       = false;
       // v_deriv_ind = 2;
    }


    for(l=0;l<MUSTA_K;l++){

        (this->*flux)(fl, ql.at(0), ql.at(1), ql.at(2), ql.at(3));
        (this->*flux)(fr, qr.at(0), qr.at(1), qr.at(2), qr.at(3));

        for(k=0;k<4;k++){
            qm.at(k) = (ql.at(k)+qr.at(k))/2 - (fr.at(k)-fl.at(k))*DT/ds/2;
        }
        (this->*flux)(fm, qm.at(0), qm.at(1), qm.at(2), qm.at(3));

        for(k=0;k<4;k++){
            flx.at(k) = (fl.at(k)+2*fm.at(k)+fr.at(k) - (qr.at(k)-ql.at(k))*ds/DT )/4;
        }

        for(k=0;k<4 && (l+1)<MUSTA_K;k++){
            ql.at(k) -= (flx.at(k)- fl.at(k)  ) * DT/ds;
            qr.at(k) -= (fr.at(k) - flx.at(k) ) * DT/ds;
        }
    }
    return;
}



void mff::init_max()
{
    vmax  = 0;
    cmax  = 0;
    mdirx = true;
    Dmax  = 0;
    return;
}


void mff::find_max(const uint& i, const uint& j)
{
    if (vmax<fabs(p->at(1, i, j))){
        vmax  = fabs(p->at(1, i, j));
        mdirx = true;
    }
    if (vmax<fabs(p->at(2, i, j))){
        vmax  = fabs(p->at(2, i, j));
        mdirx = false;
    }
    if (p->at(0,i,j)<DEF_NULL) return;
    eos->ck(p->at(3, i, j)/p->at(0, i, j));
    findcmax_c = sqrt( (1+kappa) * p->at(3, i, j) / p->at(0, i, j));
    
    if (cmax<findcmax_c){
        cmax = findcmax_c;
    }

    return;
}

void mff::CalcThings()
{
	v_deriv[0] = (qr[1]/qr[0]-ql[1]/ql[0])/dx;
	v_deriv[1] = (qr[2]/qr[0]-ql[2]/ql[0])/dx;
	v_deriv[2] = (QR[1]/QR[0]-QL[1]/QL[0])/dy;
	v_deriv[3] = (QR[2]/QR[0]-QL[2]/QL[0])/dy;
	s = (1+kappa) * 0.25*(ql[0]+qr[0]+QR[0]+QL[0]);
	e = 0.25*(ql[3]+qr[3]+QR[3]+QL[3])/kappa;
	u = 0.25*(qr[1]/qr[0]+ql[1]/ql[0]+QR[1]/QR[0]+QL[1]/QL[0]);
	v = 0.25*(qr[2]/qr[0]+ql[2]/ql[0]+QR[2]/QR[0]+QL[2]/QL[0]);
	return;
}


void mff::update_dt()
{
    if (cmax==0 && vmax==0){
        return;
    }
    //cout << "vmax="<<vmax<<"; cmax="<<cmax<<endl<<endl;
    dt = CFLN_WV * ( mdirx ? dx : dy ) / (vmax+cmax);
    if (u1>0){
        dtype T1 = CFLN * ( (mdirx) ? dx : dy ) / (vmax+cmax);
        dtype T2 = CFLN_VISC  * pow((dx<dy ? dx : dy), 2)/(Dmax==0 ? 1 : Dmax);
        dt = T1<T2 ? T1 : T2;

        /*if (pow((vmax+cmax)*T1/dx,2)>T2*Dmax/dx/dx){
            cerr<<"CFL BIG! STABILITY CRITERIA CAN'T BE SATISFIED!"<<endl;
        }
        if (T2*Dmax/dx/dx>1){
            cerr<<"STABILITY ERROR! 2s>1"<<endl;
        }*/
    }
    return;
}


void mff::step()
{
    uint i, j, k;
    init_max();
    for(i=1;i<=(r+3);i++){
    	for(j=1;j<=c;j++){
    		if (i<=r){
    			boundary(i, j);
    			boundary(i+1, j+1);
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i-1, j);
    				qr[k] = q.at(k, i, j);
    			}
    			iflux(&mff::flux_f);
    			flx2 = flx;
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i, j);
    				qr[k] = q.at(k, i+1, j);
    			}
    			iflux(&mff::flux_f);
    			for(k=0;k<4;k++){
    				k1.at(k, i, j)  = (flx2[k]-flx[k])/dx;
    				k1.at(k, 0, j)  = k1.at(k, 1, j);
    			}
    		}

    		if (i==1) continue;

    		if (i<=(r+1)){
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i-2, j) + 0.5*dt*k1.at(k, i-2, j);
    				qr[k] = q.at(k, i-1, j) + 0.5*dt*k1.at(k, i-1, j);
    			}
				iflux(&mff::flux_f, true);
				flx2 = flx;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-1, j) + 0.5*dt*k1.at(k, i-1, j);
					qr[k] = q.at(k, i, j)   + 0.5*dt*k1.at(k, i, j);
				}
				iflux(&mff::flux_f, true);
				for(k=0;k<4;k++){
					k2.at(k, i-1, j) = (flx2[k]-flx[k])/dx;
					k2.at(k, 0,   j) =  k2.at(k, 1, j);
				}
    		}

    		if (i==2) continue;

			if (i<=(r+2)){
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-3, j) + 0.5*dt*k2.at(k, i-3, j);
					qr[k] = q.at(k, i-2, j) + 0.5*dt*k2.at(k, i-2, j);
				}
				iflux(&mff::flux_f, true);
				flx2 = flx;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-2, j) + 0.5*dt*k2.at(k, i-2, j);
					qr[k] = q.at(k, i-1, j) + 0.5*dt*k2.at(k, i-1, j);
				}
				iflux(&mff::flux_f, true);
				for(k=0;k<4;k++){
					k3.at(k, i-2, j) = (flx2[k]-flx[k])/dx;
					k3.at(k, 0, j)   = k3.at(k, 1, j);
				}
			}

			if (i==3) continue;

			for(k=0;k<4;k++){
				ql[k] = q.at(k, i-4, j) + dt*k3.at(k, i-4, j);
				qr[k] = q.at(k, i-3, j) + dt*k3.at(k, i-3, j);
			}
			iflux(&mff::flux_f);
			flx2 = flx;
			for(k=0;k<4;k++){
				ql[k] = q.at(k, i-3, j) + dt*k3.at(k, i-3, j);
				qr[k] = q.at(k, i-2, j) + dt*k3.at(k, i-2, j);
			}
			iflux(&mff::flux_f);
			for(k=0;k<4;k++){
				k4  = (flx2[k]-flx[k])/dx;
				q.at(k,  i-3, j) += (dt/6)*(k1.at(k, i-3, j)+2*k2.at(k, i-3, j)+2*k3.at(k, i-3, j)+k4);
			}
    	}
    }


    for(i=1;i<=r;i++){
    	for(j=1;j<=(c+3);j++){
    		if (j<=c){
    			boundary(i, j);
    			boundary(i+1, j+1);
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i, j-1);
    				qr[k] = q.at(k, i, j);

    			}
    			iflux(&mff::flux_g);
    			flx2 = flx;
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i, j);
    				qr[k] = q.at(k, i, j+1);
    			}
    			iflux(&mff::flux_g);
    			for(k=0;k<4;k++){
    				k1.at(k, i, j)  = (flx2[k]-flx[k])/dx;
    				k1.at(k, i, 0)  = k1.at(k, i, 1);
    			}
    		}

    		if (j==1) continue;

    		if (j<=(c+1)){
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i, j-2) + 0.5*dt*k1.at(k, i, j-1);
    				qr[k] = q.at(k, i, j-1) + 0.5*dt*k1.at(k, i, j-2);
    			}
				iflux(&mff::flux_g, true);
				flx2 = flx;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i, j-1) + 0.5*dt*k1.at(k, i, j-1);
					qr[k] = q.at(k, i, j)   + 0.5*dt*k1.at(k, i, j);
				}
				iflux(&mff::flux_g, true);
				for(k=0;k<4;k++){
					k2.at(k, i, j-1) = (flx2[k]-flx[k])/dx;
					k2.at(k, i,   0) =  k2.at(k, i, 1);
				}
    		}

    		if (j==2) continue;

			if (j<=(c+2)){
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i, j-3) + 0.5*dt*k2.at(k, i, j-3);
					qr[k] = q.at(k, i, j-2) + 0.5*dt*k2.at(k, i, j-2);
				}
				iflux(&mff::flux_g, true);
				flx2 = flx;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i, j-2) + 0.5*dt*k2.at(k, i, j-2);
					qr[k] = q.at(k, i, j-1) + 0.5*dt*k2.at(k, i, j-1);
				}
				iflux(&mff::flux_g, true);
				for(k=0;k<4;k++){
					k3.at(k, i, j-2) = (flx2[k]-flx[k])/dx;
					k3.at(k, i, 0)   = k3.at(k, i, 1);
				}
			}

			if (j==3) continue;

			for(k=0;k<4;k++){
				ql[k] = q.at(k, i, j-4) + dt*k3.at(k, i, j-4);
				qr[k] = q.at(k, i, j-3) + dt*k3.at(k, i, j-3);
			}
			iflux(&mff::flux_g);
			flx2 = flx;
			for(k=0;k<4;k++){
				ql[k] = q.at(k, i, j-3) + dt*k3.at(k, i, j-3);
				qr[k] = q.at(k, i, j-2) + dt*k3.at(k, i, j-2);
			}
			iflux(&mff::flux_g);
			for(k=0;k<4;k++){
				k4  = (flx2[k]-flx[k])/dx;
				q.at(k,  i, j-3) += (dt/6)*(k1.at(k, i, j-3)+2*k2.at(k, i, j-3)+2*k3.at(k, i, j-3)+k4);
			}
			if (u1<=0){
				qp(i, j-3);
				find_max(i, j-3);
				if (!call_func.empty() || call_func.size()!=call_context.size()){
					for(k=0;k<call_func.size();k++){
						call_func.at(k)(call_context.at(k), i, j-3);
					}
				}
			}
    	}

    }
    

    if (u1>0){
    for(i=1;i<=(r+3);i++){
    	for(j=1;j<=(c+3);j++){
    		if (i<=r && j<=c){
    			boundary(i, j);
    			boundary(i+1, j+1);
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i-1, j);
    				qr[k] = q.at(k, i, j);
    				QL[k] = q.at(k, i, j-1);
    				QR[k] = q.at(k, i, j);
    			}
    			CalcThings();
    	        vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    	        vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    			flx2   = flx;
    			flx2_y = flx_y;
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i, j);
    				qr[k] = q.at(k, i+1, j);
    				QL[k] = q.at(k, i, j);
    				QR[k] = q.at(k, i, j+1);
    			}
    			CalcThings();
    			vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    	        vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    			for(k=0;k<4;k++){
    				k1.at(k, i, j)   = (flx2[k]-flx[k])/dx;
    				k1y.at(k, i, j)  = (flx2_y[k]-flx_y[k])/dy;
    				k1.at(k, 0, j)   = k1.at(k, 1, j);
    				k1y.at(k, i, 0)  = k1.at(k, 1, 0);
    			}
    		}

    		if (j==1 || i==1) continue;

    		if (i<=(r+1) && j<=(c+1)){
    			for(k=0;k<4;k++){
    				ql[k] = q.at(k, i-2, j-1) + 0.5*dt*k1.at(k, i-2, j-1);
    				qr[k] = q.at(k, i-1, j-1) + 0.5*dt*k1.at(k, i-1, j-1);
    				QL[k] = q.at(k, i-1, j-2) + 0.5*dt*k1y.at(k, i-1, j-2);
    				QR[k] = q.at(k, i-1, j-1) + 0.5*dt*k1y.at(k, i-1, j-1);
    			}
    			CalcThings();
    			vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    			vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
    			flx2   = flx;
    			flx2_y =flx_y;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-1, j-1) + 0.5*dt*k1.at(k, i-1, j-1);
					qr[k] = q.at(k, i, j-1)   + 0.5*dt*k1.at(k, i, j-1);
					QL[k] = q.at(k, i-1, j-1) + 0.5*dt*k1y.at(k, i-1, j-1);
					QR[k] = q.at(k, i-1, j)   + 0.5*dt*k1y.at(k, i-1, j);
				}
				CalcThings();
				vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				for(k=0;k<4;k++){
					k2.at(k, i-1, j-1)  = (flx2[k]-flx[k])/dx;
					k2y.at(k, i-1, j-1) = (flx2_y[k]-flx_y[k])/dy;
					k2.at(k, 0,   j-1) =  k2.at(k, 1, j-1);
					k2y.at(k, i-1,  0) =   k2y.at(k, i-1, 0);
				}
    		}

    		if (j==2 || i==2) continue;

			if (i<=(r+2) && j<=(c+1)){
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-3, j-2) + 0.5*dt*k2.at(k, i-3, j-2);
					qr[k] = q.at(k, i-2, j-2) + 0.5*dt*k2.at(k, i-2, j-2);
					QL[k] = q.at(k, i-2, j-3) + 0.5*dt*k2y.at(k, i-2, j-3);
					QR[k] = q.at(k, i-2, j-2) + 0.5*dt*k2y.at(k, i-2, j-2);
				}
				CalcThings();
				vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				flx2   = flx;
				flx2_y = flx_y;
				for(k=0;k<4;k++){
					ql[k] = q.at(k, i-2, j-2) + 0.5*dt*k2.at(k, i-2, j-2);
					qr[k] = q.at(k, i-1, j-2) + 0.5*dt*k2.at(k, i-1, j-2);
					QL[k] = q.at(k, i-2, j-2) + 0.5*dt*k2y.at(k, i-2, j-2);
					QR[k] = q.at(k, i-2, j-1) + 0.5*dt*k2y.at(k, i-2, j-1);
				}
				CalcThings();
				vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
				for(k=0;k<4;k++){
					k3.at(k, i-2, j-2) = (flx2[k]-flx[k])/dx;
					k3y.at(k, i-2, j-2)= (flx2_y[k]-flx_y[k])/dy;
					k3.at(k, 0, j-2)   = k3.at(k, 1, j-2);
					k3y.at(k, i-2, 0)   = k3y.at(k, i-2, 0);
				}
			}

			if (i==3 || j==3) continue;

			for(k=0;k<4;k++){
				ql[k] = q.at(k, i-4, j-3) + dt*k3.at(k, i-4, j-3);
				qr[k] = q.at(k, i-3, j-3) + dt*k3.at(k, i-3, j-3);
				QL[k] = q.at(k, i-3, j-4) + dt*k3y.at(k, i-3, j-4);
				QR[k] = q.at(k, i-3, j-3) + dt*k3y.at(k, i-3, j-3);
			}
			CalcThings();
			vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
			vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
			flx2   = flx;
			flx2_y = flx_y;
			for(k=0;k<4;k++){
				ql[k] = q.at(k, i-3, j-3) + dt*k3.at(k, i-3, j-3);
				qr[k] = q.at(k, i-2, j-3) + dt*k3.at(k, i-2, j-3);
				QL[k] = q.at(k, i-3, j-3) + dt*k3y.at(k, i-3, j-3);
				QR[k] = q.at(k, i-3, j-2) + dt*k3y.at(k, i-3, j-2);
			}
			CalcThings();
			vflux_f(flx,   s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
			vflux_g(flx_y, s, u, v, v_deriv.at(0), v_deriv.at(1), v_deriv.at(2), v_deriv.at(3));
			for(k=0;k<4;k++){
				k4  = (flx2[k]-flx[k])/dx;
				k4y = (flx2_y[k]-flx_y[k])/dy;
				q.at(k,  i-3, j-3) += (dt/6)*(k1.at(k, i-3, j-3)+2*k2.at(k, i-3, j-3)+2*k3.at(k, i-3, j-3)+k4);
				q.at(k,  i-3, j-3) += (dt/6)*(k1y.at(k, i-3, j-3)+2*k2y.at(k, i-3, j-3)+2*k3y.at(k, i-3, j-3)+k4y);
			}
			if (Dmax<fabs((u*u+v*v)*4*u1*s/e*2)){
				Dmax  = fabs((u*u+v*v)*4*u1*s/e*2);
			}
			qp(i-3, j-3);
			find_max(i-3, j-3);
			if (!call_func.empty() || call_func.size()!=call_context.size()){
				for(k=0;k<call_func.size();k++){
					call_func.at(k)(call_context.at(k), i-3, j-3);
				}
			}
    	}
    }}
    update_dt();
    return;
}

