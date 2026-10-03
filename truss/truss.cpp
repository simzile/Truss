#include <iostream>
#include <math.h>
#include <fstream>
#include "json.hpp"
#include <string>
#define space  "       |       "
//float doh(unsigned int i,const float* dp,unsigned int m);
class Read{
    public: 
        // total number of nodes
        unsigned int n=0;
        unsigned int nn=0;
        // number of fixed nodes
        int nc=0;
        // number of free nodes
        unsigned int free=n-nc;
        // number of rods
        unsigned int m=0;
        // number of loads
        unsigned int mf=0;
        // number of sections
        unsigned int numof_sections=0;

        unsigned int ndof=0;

        unsigned int* loca_loads;
        double* x_loads;
        double* y_loads;
        // x coordinate of nodes
        double* x=nullptr;
        
        // y coordinate of nodes
        double* y=nullptr;
        
        // beginning node of nodes
        unsigned int* ihl=nullptr;
        // end node of nodes
        unsigned int* ihr=nullptr;

        //material of nodes
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
            nn=2*(n-nc);
            //need bios
            m=data.at("number_of_rods").get<int>();
            x=new double[n];
            y=new double[n];
            for(unsigned int i=0;i<n;i++){
                x[i]=data.at("x_for_nodes").at(i).get<double>();
                y[i]=data.at("y_for_nodes").at(i).get<double>();
            }
            ihl=new unsigned int[m];
            ihr=new unsigned int[m];
            EA=new double[m];

            mf=data.at("location_of_loads_by_rods").size();
            loca_loads=new unsigned int[mf];
            x_loads=new double[mf];
            y_loads=new double[mf];

            for(unsigned int i=0;i<m;i++){
                ihl[i]=data.at("left_node_of_rod").at(i).get<unsigned int>();
                ihr[i]=data.at("right_node_of_rod").at(i).get<unsigned int>();
                EA[i]=data.at("EA").at(i).get<double>();
            }

            for(unsigned int i=0;i<mf;i++){
                loca_loads[i]=data.at("location_of_loads_by_rods").at(i).get<unsigned int>();
                x_loads[i]=data.at("x_of_loads").at(i).get<double>();
                y_loads[i]=data.at("y_of_loads").at(i).get<double>();
            }
        }
};
class Kmatrix{
    private:
    public:        
        std::array<double,3> ch(unsigned int k){
            /*
                compute the length and theta
            */
            Read rod("input.json");
            unsigned int ihl=rod.ihl[k-1];
            unsigned int ihr=rod.ihr[k-1];
            // xi yi xj yj
            double xi=rod.x[ihl-1],yi=rod.y[ihl-1],xj=rod.x[ihr-1],yj=rod.y[ihr-1];
            double l=sqrt((xj-xi)*(xj-xi)+(yj-yi)*(yj-yi));
            double cos=(xj-xi)/l;
            double sin=(yj-yi)/l;
            return {l,cos,sin};   
        }
        std::array<int,2> i0j0(unsigned int k){
            /*
                input kth --- i j
                output i0 j0
                i0=2*(i0-nc-1)
            */
            Read rod("input.json");
            int i=rod.ihl[k-1];
            int j=rod.ihr[k-1];
            int i0=2*(i-rod.nc-1);
            int j0=2*(j-rod.nc-1);
            return {i0,j0};
        }    
        void stif(unsigned int k,double* ke){
            /*
                compute for elemental stiffness matrix in the global coordinate
                k=EA/l*[1,0;0,1],  kk=T*k*TT
            */
            //for kth rod
            Read rod("input.json");
            double EA=rod.EA[k-1];
            double lk=ch(k)[0];
            double cos=ch(k)[1];
            double sin=ch(k)[2];
            double stiffness=EA/lk;
            *ke=stiffness*cos*cos;
            *(ke+1)=stiffness*cos*sin;
            *(ke+2)=stiffness*cos*sin;
            *(ke+3)=stiffness*sin*sin;
            //clxy[3]={l,cos,sin};
        }
        void dor(double* r,size_t nn){
            Read rod("input.json");
            // m is the num of rods
            for(unsigned int k=1;k<=rod.m;k++){   
                int i0=i0j0(k)[0];
                int j0=i0j0(k)[1];
                //unsigned int nn=rod.ndof;
                //float* r=new float[nn*nn]{};
                //store K matrix by column K[i*ndof+j]=K[i][j]
                double rd[4]{};

                //ke(i*2+j)=k(i,j)
                unsigned int truss2d[2]={0,1};
                //both are free
                stif(k,rd);
                if(i0>=0){      
                    for(unsigned int i:truss2d){
                        for(unsigned int j:truss2d){
                            //(i0,i0)
                            r[(i0+i)*nn+(i0+j)]+=rd[i*2+j];
                            //(i0,j0),(j0,i0)
                            r[(i0+i)*nn+(j0+j)]-=rd[i*2+j];
                            r[(j0+i)*nn+(i0+j)]-=rd[i*2+j];
                            //(j0,j0)
                            r[(j0+i)*nn+(j0+j)]+=rd[i*2+j];
                            j++;
                        }
                        i++; 
                    }
                    /*
                    values[]=
                    for(int x : values)
                    */
                }
                else{
                    // only for (j0,j0)
                    for(unsigned int i:truss2d){
                        for(unsigned int j:truss2d){
                            r[(j0+i)*nn+(j0+j)]+=rd[i*2+j];
                            j++;
                        }
                        i++; 
                    }
                } 
            }
            //we get the global stiffness matrix K[ndof,ndof]
        }
        double doh(unsigned int i,const double* dp,unsigned int m){
            /*
                compute the axial force in the section of rods
            */
            Read rod("input.json");
            Kmatrix K;    
            /*
            sx    k11 k12 xi
            sy    k21 k22 yi
            */
            // i=1,2,3,4,5,6
            double EA=rod.EA[i-1];
            double l=K.ch(i)[0];
            double cos=K.ch(i)[1];
            double sin=K.ch(i)[2];
            double internal_force=0;
            int i0=K.i0j0(i)[0];
            //int i0=2*(i-rod.nc);
            int j0=K.i0j0(i)[1];
            if(i0>=0){
                internal_force=(EA/l)*((dp[j0]*cos+dp[j0+1]*sin)-(dp[i0]*cos+dp[i0+1]*sin));
            } 
            else{
                internal_force=(EA/l)*(dp[j0]*cos+dp[j0+1]*sin);
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
                //z=M_inv*r
                //1/A[i,i]*r[i]
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
            //x=x+alpha*p
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
        };
        void equ(const double* R,double* r,size_t n){
            for(unsigned int i=0;i<n;i++){
                r[i]=R[i];
            }
        };

    public:
        void CG(const double* K,double*dp,const double* f,size_t nn){
            double* p=new double[nn]{};
            double* r=new double[nn]{};
            double* z=new double[nn]{};
            res(f,K,dp,r,nn);
            precondition(z,K,r,nn);
            equ(z,p,nn);// p=r
            double* R=new double[nn]{};
            double* Z=new double[nn]{};
            double tol=1;
            
            unsigned int k=0;
            while(tol>1e-2){

                alpha=Alpha(r,z,p,K,nn); //compute for searching distance
                X(p,dp,alpha,nn); //update dp
                res(f,K,dp,R,nn); //compute the new bias
                precondition(Z,K,R,nn);
                
                beta=Beta(r,R,z,Z,nn); //compute beta
                P(Z,beta,p,nn); //compute for the searching direction
                equ(R,r,nn); // update r
                equ(Z,z,nn);
                //beta=R*R/r*r                
                //p=r+beta*p
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
            //std::ios::app
            std::ofstream displace("displacements.txt");
            if(!displace) return;
            displace<<"num of nodes | dx | dy | theta\n";
            for(unsigned int i=0;i<n;i++){
                if(i<rod.nc){
                    //fixed   
                    displace<<i+1<<space<<0<<space<<0<<space<<0<<"\n";
                }
                else{
                    int i0=2*(i-nc);
                    double dx=0,dy=0;
                    dx=dp[i0];
                    dy=dp[i0+1];
                    displace<<i+1<<space<<dx<<space<<dy<<space<<0<<"\n";
                }

            }
            
            displace.close();
        }
        void write_stress(const double* dp){

            Read rod("input.json");
            Kmatrix K;
            std::ofstream stress("stress.txt");
            if(!stress) return;
            stress<<"num of rods"<<space<<"length"<<space<<"theta"<<'\n';
            double internal_force=0;
            for(unsigned int i=1;i<=rod.m;i++){
                    double length=K.ch(i)[0];
                    internal_force=K.doh(i,dp,rod.m);
                    stress<<i+1<<space<<length<<space<<internal_force<<space<<0<<"\n";
            }
            stress.close();
        }
        
};
void F(double* f){
    Read rod("input.json");
    Kmatrix K;

    for(unsigned int i=0;i<rod.mf;i++){
        unsigned int j=rod.loca_loads[i];
        int i0=2*(j-rod.nc-1);
        f[i0]+=rod.x_loads[i];
        f[i0+1]+=rod.y_loads[i];
    }
}
void slveq(const double* r,double* dp,const double* f){
    Read rod("input.json");
    solver solve;
    solve.CG(r,dp,f,rod.nn);
    fout write;
    // displacement of all nodes dp[nn]
    write.write_disp(dp,rod.n,rod.nc);
}
int main(){    
    std::string path="input.json";
    //std::cin>>path;
    Read rod(path);
    size_t nn=rod.nn;
    double* r=new double[nn*nn]{};
    double* dp=new double[nn]{};
    double* f=new double[nn]{};
    F(f);
    Kmatrix K;
    K.dor(r,nn);// get K
    slveq(r,dp,f);// get dp
    fout stress;
    stress.write_stress(dp);
    return 0;
}