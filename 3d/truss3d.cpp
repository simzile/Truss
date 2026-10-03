#include <iostream>
#include <math.h>
#include <fstream>
#include "json.hpp"
#include <string>
#define space  "       |       "
class Read{
    public:
        unsigned int n=0;
        unsigned int nn=0;
        int nc=0;
        unsigned int free=n-nc;
        unsigned int m=0;
        unsigned int mf=0;
        unsigned int numof_sections=0;
        unsigned int ndof=0;
        unsigned int* loca_loads;
        double* x_loads;
        double* y_loads;
        double* z_loads;
        double* x=nullptr;
        double* y=nullptr;
        double* z=nullptr;
        unsigned int* ihl=nullptr;
        unsigned int* ihr=nullptr;
        double* EA=nullptr;
        Read(const std::string& path){
            load(path);
        }
        void load(const std::string& path){
            std::ifstream file(path);
            if(!file)  {std::cerr<<"cannot open input.json in this path\n"<<std::endl;}
            nlohmann::json data=nlohmann::json::parse(file);
            n=data.at("number_of_total_nodes").get<int>();
            nc=data.at("number_of_fixed_nodes").get<int>();
            nn=3*(n-nc);
            m=data.at("number_of_rods").get<int>();
            x=new double[n];
            y=new double[n];
            z=new double[n];
            for(unsigned int i=0;i<n;i++){
                x[i]=data.at("x_for_nodes").at(i).get<double>();
                y[i]=data.at("y_for_nodes").at(i).get<double>();
                z[i]=data.at("z_for_nodes").at(i).get<double>();
            }
            ihl=new unsigned int[m];
            ihr=new unsigned int[m];
            EA=new double[m];
            mf=data.at("location_of_loads_by_rods").size();
            loca_loads=new unsigned int[mf];
            x_loads=new double[mf];
            y_loads=new double[mf];
            z_loads=new double[mf];
            for(unsigned int i=0;i<m;i++){
                ihl[i]=data.at("left_node_of_rod").at(i).get<unsigned int>();
                ihr[i]=data.at("right_node_of_rod").at(i).get<unsigned int>();
                EA[i]=data.at("EA").at(i).get<double>();
            }
            for(unsigned int i=0;i<mf;i++){
                loca_loads[i]=data.at("location_of_loads_by_rods").at(i).get<unsigned int>();
                x_loads[i]=data.at("x_of_loads").at(i).get<double>();
                y_loads[i]=data.at("y_of_loads").at(i).get<double>();
                z_loads[i]=data.at("z_of_loads").at(i).get<double>();
            }
        }
};
class Kmatrix{
    private:
    public:
        std::array<double,4> ch(unsigned int k){
            Read rod("input.json");
            unsigned int ihl=rod.ihl[k-1];
            unsigned int ihr=rod.ihr[k-1];
            double xi=rod.x[ihl-1],yi=rod.y[ihl-1],zi=rod.z[ihl-1],xj=rod.x[ihr-1],yj=rod.y[ihr-1],zj=rod.z[ihr-1];
            double l=sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi));
            double lx=(xj-xi)/l;
            double ly=(yj-yi)/l;
            double lz=(zj-zi)/l;
            return {l,lx,ly,lz};
        }
        std::array<int,2> i0j0(unsigned int k){
            Read rod("input.json");
            int i=rod.ihl[k-1];
            int j=rod.ihr[k-1];
            int i0=3*(i-rod.nc-1);
            int j0=3*(j-rod.nc-1);
            return {i0,j0};
        }
        void stif(unsigned int k,double* ke){
            Read rod("input.json");
            double EA=rod.EA[k-1];
            double lk=ch(k)[0];
            double lx=ch(k)[1];
            double ly=ch(k)[2];
            double lz=ch(k)[3];
            double stiffness=EA/lk;
            double l[3]={lx,ly,lz};
            for(unsigned int i=0;i<3;i++){
                for(unsigned int j=0;j<3;j++){
                    ke[i*3+j]=stiffness*l[i]*l[j];
                }
            }
        }
        void dor(double* r,size_t nn){
            Read rod("input.json");
            for(unsigned int k=1;k<=rod.m;k++){
                int i0=i0j0(k)[0];
                int j0=i0j0(k)[1];
                double rd[9]{};
                unsigned int truss3d[3]={0,1,2};
                stif(k,rd);
                if(i0>=0){
                    for(unsigned int i:truss3d){
                        for(unsigned int j:truss3d){
                            r[(i0+i)*nn+(i0+j)]+=rd[i*3+j];
                            r[(i0+i)*nn+(j0+j)]-=rd[i*3+j];
                            r[(j0+i)*nn+(i0+j)]-=rd[i*3+j];
                            r[(j0+i)*nn+(j0+j)]+=rd[i*3+j];
                        }
                    }
                }
                else{
                    for(unsigned int i:truss3d){
                        for(unsigned int j:truss3d){
                            r[(j0+i)*nn+(j0+j)]+=rd[i*3+j];
                        }
                    }
                }
            }
        }
        double doh(unsigned int i,const double* dp,unsigned int m){
            Read rod("input.json");
            Kmatrix K;
            double EA=rod.EA[i-1];
            double l=K.ch(i)[0];
            double lx=K.ch(i)[1];
            double ly=K.ch(i)[2];
            double lz=K.ch(i)[3];
            double internal_force=0;
            int i0=K.i0j0(i)[0];
            int j0=K.i0j0(i)[1];
            if(i0>=0){
                internal_force=(EA/l)*((dp[j0]*lx+dp[j0+1]*ly+dp[j0+2]*lz)-(dp[i0]*lx+dp[i0+1]*ly+dp[i0+2]*lz));
            }
            else{
                internal_force=(EA/l)*(dp[j0]*lx+dp[j0+1]*ly+dp[j0+2]*lz);
            }
            return internal_force;
        }
};
class solver{
    private:
        double alpha=0;
        double beta=0;
        void precondition(double* z,const double* A,double* r,size_t n){
            for(unsigned int i=0;i<n;i++){
                z[i]=r[i]/(A[i*n+i]);
            }
        }
        double Tol(const double* r,size_t n){
            double norm=0;
            for(unsigned int i=0;i<n;i++){
                norm+=(r[i]*r[i]);
            }
            return norm;
        }
        void res(const double*b,const double* A,const double* x,double* r,size_t n){
            double* ax=new double[n]{};
            for(unsigned int j=0;j<n;j++){
                for(unsigned int i=0;i<n;i++){
                    ax[j]+=A[j*n+i]*x[i];
                }
                r[j]=b[j]-ax[j];
            }
            delete[] ax;
        }
        double Alpha(const double* r,const double* z,const double* p,const double* A,size_t n){
            double rr=0,pAp=0;
            for(unsigned int i=0;i<n;i++){
                rr+=(r[i]*z[i]);
                double* Ap=new double[n]{};
                for(unsigned int j=0;j<n;j++){
                    Ap[i]+=A[i*n+j]*p[j];
                }
                pAp+=(Ap[i]*p[i]);
                delete[] Ap;
            }
            double alpha=rr/pAp;
            return alpha;
        };
        void X(const double* p,double* x,const double alpha,size_t n){
            for(unsigned int i=0;i<n;i++){
                x[i]+=(alpha*p[i]);
            }
        }
        double Beta(const double* r,const double* R,const double* z,double* Z,size_t n){
            double RR=0,rr=0;
            for(unsigned int i=0;i<n;i++){
                RR+=R[i]*Z[i];
                rr+=r[i]*z[i];
            }
            double beta=RR/rr;
            return beta;
        };
        void P(const double* R,const double beta,double* p,size_t n){
            for(unsigned int i=0;i<n;i++){
                p[i]=R[i]+beta*p[i];
            }
        }
        void equ(const double* R,double* r,size_t n){
            for(unsigned int i=0;i<n;i++){
                r[i]=R[i];
            }
        }
    public:
        void CG(const double* K,double*dp,const double* f,size_t nn){
            double* p=new double[nn]{};
            double* r=new double[nn]{};
            double* z=new double[nn]{};
            res(f,K,dp,r,nn);
            precondition(z,K,r,nn);
            equ(z,p,nn);
            double* R=new double[nn]{};
            double* Z=new double[nn]{};
            double tol=1;
            unsigned int k=0;
            while(tol>1e-2){
                alpha=Alpha(r,z,p,K,nn);
                X(p,dp,alpha,nn);
                res(f,K,dp,R,nn);
                precondition(Z,K,R,nn);
                beta=Beta(r,R,z,Z,nn);
                P(Z,beta,p,nn);
                equ(R,r,nn);
                equ(Z,z,nn);
                tol=Tol(r,nn);
                k++;
                printf("%d th iteration,tolerance is %.9e\n",k,tol);
            }
            delete[] p;
            delete[] r;
            delete[] z;
            delete[] R;
            delete[] Z;
        }
};
class fout{
    public:
        void write_disp(const double* dp,int n,int nc){
            Read rod("input.json");
            std::ofstream displace("displacements.txt");
            if(!displace) return;
            displace<<"num of nodes | dx | dy | dz\n";
            for(unsigned int i=0;i<n;i++){
                if(i<rod.nc){
                    displace<<i+1<<space<<0<<space<<0<<space<<0<<"\n";
                }
                else{
                    int i0=3*(i-nc);
                    displace<<i+1<<space<<dp[i0]<<space<<dp[i0+1]<<space<<dp[i0+2]<<"\n";
                }
            }
            displace.close();
        }
        void write_stress(const double* dp){
            Read rod("input.json");
            Kmatrix K;
            std::ofstream stress("stress.txt");
            if(!stress) return;
            stress<<"num of rods"<<space<<"length"<<space<<"lx"<<space<<"ly"<<space<<"lz"<<space<<"axial force"<<'\n';
            double internal_force=0;
            for(unsigned int i=1;i<=rod.m;i++){
                double length=K.ch(i)[0];
                internal_force=K.doh(i,dp,rod.m);
                stress<<i<<space<<length<<space<<K.ch(i)[1]<<space<<K.ch(i)[2]<<space<<K.ch(i)[3]<<space<<internal_force<<"\n";
            }
            stress.close();
        }
};
void F(double* f){
    Read rod("input.json");
    Kmatrix K;
    for(unsigned int i=0;i<rod.mf;i++){
        unsigned int j=rod.loca_loads[i];
        int i0=3*(j-rod.nc-1);
        f[i0]+=rod.x_loads[i];
        f[i0+1]+=rod.y_loads[i];
        f[i0+2]+=rod.z_loads[i];
    }
}
void slveq(const double* r,double* dp,const double* f){
    Read rod("input.json");
    solver solve;
    solve.CG(r,dp,f,rod.nn);
    fout write;
    write.write_disp(dp,rod.n,rod.nc);
}
int main(){
    std::string path="input.json";
    Read rod(path);
    size_t nn=rod.nn;
    double* r=new double[nn*nn]{};
    double* dp=new double[nn]{};
    double* f=new double[nn]{};
    F(f);
    Kmatrix K;
    K.dor(r,nn);
    slveq(r,dp,f);
    fout stress;
    stress.write_stress(dp);
    return 0;
}
