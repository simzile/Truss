#include <iostream>
#include <math.h>
#include <fstream>
#include "json.hpp"
#include <string>
#define space  "       |       "
//float doh(unsigned int i,const float* dp,unsigned int m);
void force(unsigned int ind,double* rx,double* ry,double* m,double dist,double l,double g);

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

        unsigned int* loca_loads=nullptr;
        unsigned int* category=nullptr;
        double* magnitude=nullptr;
        double* dist_load=nullptr;
        double* moments=nullptr;

        unsigned int* loca_sect=nullptr;
        double* dist_sect=nullptr;
        // x coordinate of nodes
        double* x=nullptr;
        
        // y coordinate of nodes
        double* y=nullptr;

        double* theta=nullptr;
        
        // beginning node of nodes
        unsigned int* ihl=nullptr;
        // end node of nodes
        unsigned int* ihr=nullptr;

        unsigned int nqq=0;
        unsigned int npp=0;
        //material of nodes
        double* E=nullptr;
        double* A=nullptr;
        double* I=nullptr;
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
            E=new double[m];
            A=new double[m];
            I=new double[m];

            mf=data.at("location_of_loads_by_rods").size();
            loca_loads=new unsigned int[mf];
            category=new unsigned int[mf];
            magnitude=new double[mf];
            dist_load=new double[mf];

            for(unsigned int i=0;i<m;i++){
                ihl[i]=data.at("left_node_of_rod").at(i).get<unsigned int>();
                ihr[i]=data.at("right_node_of_rod").at(i).get<unsigned int>();
                E[i]=data.at("E").at(i).get<double>();
                A[i]=data.at("A").at(i).get<double>();
                I[i]=data.at("I").at(i).get<double>();
            }
            
            nqq=data.at("number_of_loads").get<unsigned int>();
            npp=data.at("number_of_section").get<unsigned int>();
            for(unsigned int ild=0;ild<nqq;ild++){
                loca_loads[ild]=data.at("location_of_loads_by_rods").at(ild).get<unsigned int>();
                category[ild]=data.at("category_of_loads").at(ild).get<unsigned int>();
                magnitude[ild]=data.at("magnitude_of_loads").at(ild).get<double>();
                dist_load[ild]=data.at("distance_between_loads_and_leftend").at(ild).get<double>();
            }
            loca_sect=new unsigned int[npp]{};
            dist_sect=new double[npp]{};
            for(unsigned int j=0;j<npp;j++){
                loca_sect[j]=data.at("location_of_section_by_rods").at(j).get<unsigned int>();
                dist_sect[j]=data.at("distance_between_section_and_leftend").at(j).get<double>();
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
            int i0=3*(i-rod.nc-1);
            int j0=3*(j-rod.nc-1);
            return {i0,j0};
        }
        void stif(unsigned int k,double* rd){
            /*
                compute for elemental stiffness matrix in the global coordinate
                k=EA/l*[1,0;0,1],  kk=T*k*TT
            */
            //for kth rod
            Read rod("input.json");
            double E=rod.E[k-1];
            double A=rod.A[k-1];
            double I=rod.I[k-1];
            double lk=ch(k)[0];
            double cos=ch(k)[1];
            double sin=ch(k)[2];
            double k11=E*A/lk,k22=12*E*I/(lk*lk*lk),k23=-6*E*I/(lk*lk),k33=4*E*I/lk;
            double k32=k23;
            double ke[9]{};
            *ke=cos*cos*k11+sin*sin*k22;
            ke[1]=sin*cos*(k11-k22);
            ke[2]=-sin*k32;
            ke[5]=cos*k32;
            ke[4]=sin*sin*k11+cos*cos*k22;
            ke[8]=k33;
            ke[3]=*(ke+1);
            ke[6]=*(ke+2);
            ke[7]=*(ke+5);
            
            double kk[36]{};
            //1
            kk[0]=ke[0];
            kk[1]=ke[1];
            kk[2]=ke[2];
            kk[3]=-ke[0];
            kk[4]=-ke[1];
            kk[5]=ke[2];
            //2
            kk[6]=kk[1];
            kk[7]=ke[4];
            kk[8]=ke[5];
            kk[9]=-ke[1];
            kk[10]=-ke[4];
            kk[11]=ke[5];
            //3
            kk[12]=kk[2];
            kk[13]=kk[8];
            kk[14]=ke[8];
            kk[15]=-ke[2];
            kk[16]=-ke[5];
            kk[17]=0.5*ke[8];
            //4
            kk[18]=kk[3];
            kk[19]=kk[9];
            kk[20]=kk[15];
            kk[21]=ke[0];
            kk[22]=ke[1];
            kk[23]=-ke[2];
            //5
            kk[24]=kk[4];
            kk[25]=kk[10];
            kk[26]=kk[16];
            kk[27]=kk[22];
            kk[28]=ke[4];
            kk[29]=-ke[5];
            //6
            kk[30]=kk[5];
            kk[31]=kk[11];
            kk[32]=kk[17];
            kk[33]=kk[23];
            kk[34]=kk[29];
            kk[35]=ke[8];
            for(unsigned int i=0;i<36;i++){
                rd[i]=kk[i];
            }
            
            //clxy[3]={l,cos,sin};
        }
        void dor(double* r,size_t nn){
            /*
                assemble the global stiffness matrix
            */
            Read rod("input.json");
            // m is the num of rods
            for(unsigned int k=1;k<=rod.m;k++){
                int i0=i0j0(k)[0];
                int j0=i0j0(k)[1];
                //unsigned int nn=rod.ndof;
                //float* r=new float[nn*nn]{};
                //store K matrix by column K[i*ndof+j]=K[i][j]
                double rd[36]{};
                //ke(i*2+j)=k(i,j)
                unsigned int truss2d[3]={0,1,2};
                //both are free
                stif(k,rd);
                if(i0>=0){
                    for(unsigned int i:truss2d){
                        for(unsigned int j:truss2d){
                            //(i0,i0)
                            r[(i0+i)*nn+(i0+j)]+=rd[i*6+j];
                            //(i0,j0),(j0,i0)
                            r[(i0+i)*nn+(j0+j)]+=rd[i*6+j+3];
                            r[(j0+i)*nn+(i0+j)]+=rd[i*6+j+18];
                            //(j0,j0)
                            r[(j0+i)*nn+(j0+j)]+=rd[i*6+j+21];
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
                            r[(j0+i)*nn+(j0+j)]+=rd[i*6+j+21];
                            j++;
                        }
                        i++;
                    }
                }
            }
            //we get the global stiffness matrix K[ndof,ndof]
        }
        void doh(unsigned int im,const double* dp,double* qmn,double dist_sect){
            /*
                compute the axial force in the section of rods
            */
            Read rod("input.json");
            Kmatrix K;
            

            double l=K.ch(im)[0];
            double cos=K.ch(im)[1];
            double sin=K.ch(im)[2];
            double E=rod.E[im-1];
            double A=rod.A[im-1];
            double J=rod.I[im-1];
            double rd[36]{};
            stif(im,rd);

            double l_fx=0,r_fx=0.0,l_fy=0.0,r_fy=0.0,l_m=0.0,r_m=0.0;

            double Neq=0,Veq=0,Meq=0;
            double ui=0,vi=0,ti=0,uj=0,vj=0,tj=0;
            // n=Sx*c+Sy*s
            // q=-(Sx*s-Sy*c)
            // m=m+(Sx*s-Sy*c)*dist
            int i0=i0j0(im)[0];
            int j0=i0j0(im)[1];

            if(i0>=0){
                ui=cos*dp[i0]+sin*dp[i0+1];
                vi=-sin*dp[i0]+cos*dp[i0+1];
                ti=dp[i0+2];
            }
            if(j0>=0){
                uj=cos*dp[j0]+sin*dp[j0+1];
                vj=-sin*dp[j0]+cos*dp[j0+1];
                tj=dp[j0+2];
            }
            double Nd=E*A/l*(ui-uj);
            double Vd=(12*E*J/(l*l*l))*(vi-vj)+(-6*E*J/(l*l))*(ti+tj);
            double Md=(-6*E*J/(l*l))*(vi-vj)+(4*E*J/l)*ti+(2*E*J/l)*tj;
            
            double addN=0,addV=0,addM=0;
            double x=dist_sect;
            //for sec 1
            for(unsigned int ild=0;ild<rod.nqq;ild++){
                //every load
                if(rod.loca_loads[ild]!=im) continue;
                double rx[2]{},ry[2]{},mm[2]{};
                double a=rod.dist_load[ild];
                double g=rod.magnitude[ild];
                force(rod.category[ild],rx,ry,mm,rod.dist_load[ild],l,g);
                Neq+=rx[0];
                Veq+=ry[0];
                Meq+=mm[0];
                unsigned int type=rod.category[ild];
                if (type!= 1 && type != 2 && type != 6) {
                    throw std::runtime_error("doh currently supports load types 1, 2, 6 only");
                }
                if(type==1&&x>=a){
                    addV+=g;
                    addM+=g*(x-a);
                }
                else if(type==2){
                    double b=(x<a)?x:a;
                    addV+=g*b;
                    addM+=g*b*(x-b/2.0);
                }
                else if(type==6&&x>=a){
                    addM+=g;
                }
            }
            double N=Nd-Neq;
            double V=Vd-Veq;
            double M=Md-Meq;
            qmn[0]=-N+addN;    //rbx
            qmn[1]=V+addV; 
            qmn[2]=M+V*dist_sect+addM;
            // qmn[0]=-fx[0];    //rbx
            // qmn[1]=fy[0]; 
            // qmn[2]=m[0]+fy[0]*dist_sect;
            
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
            while(tol>1e-9){

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
            displace<<"num of nodes | dx | dy | theta \n";
            for(unsigned int i=0;i<n;i++){
                if(i<rod.nc){
                    //fixed
                    displace<<i+1<<space<<0<<space<<0<<space<<0<<"\n";
                }
                else{
                    int i0=3*(i-nc);
                    double dx=0,dy=0,theta=0;
                    dx=dp[i0];
                    dy=dp[i0+1];
                    theta=dp[i0+2];
                    displace<<i+1<<space<<dx<<space<<dy<<space<<theta<<"\n";
                }
            }
            displace.close();
        }
        void write_stress(const double* dp){

            Read rod("input.json");
            Kmatrix K;
            std::ofstream stress("stress.txt");
            if(!stress) return;
            stress<<"num of sections"<<space<<"interfx"<<space<<"interfy"<<space<<"moment"<<space<<"middle"<<'\n';
            double inter_fx=0,inter_fy=0,moment=0;
            double qmn[3]{};
            
            for(unsigned int i=0;i<rod.npp;i++){
                    //traverse by section

                    unsigned int loca_sect=rod.loca_sect[i];  //rod 
                    double dist_sect=rod.dist_sect[i];

                    K.doh(loca_sect,dp,qmn,dist_sect);
                    inter_fx=qmn[0];
                    inter_fy=qmn[1];
                    moment=qmn[2];
                    stress<<i+1<<space<<inter_fx<<space<<inter_fy<<space<<moment<<"\n";
            }
            stress.close();
        }
        
};
// f vec
//force(ind,rx,ry,m,dist,l,)
void force(unsigned int ind,double* rx,double* ry,double* m,double dist,double l,double g){
    double rax=rx[0],rbx=rx[1],ray=ry[0],rby=ry[1],ma=m[0],mb=m[1];
    double xq=dist;
    switch (ind){
    case 1:
        // concentrate force
        ray=g*(1+2*xq/l)*(1-xq/l)*(1-xq/l);
        rby=g-ray;
        ma=-g*xq*(1-xq/l)*(1-xq/l);
        mb=g*xq*xq*(l-xq)/(l*l);
        break;
    case 2:
        //uniform force
        ray=g*xq*(1-(xq*xq)/(l*l)+(xq*xq*xq)/(2*l*l*l));
        rby=g*xq-ray;
        ma=-g*xq*xq*(6-8*xq/l+3*xq*xq/(l*l))/12;
        mb=g*xq*xq*xq*(4-3*xq/l)/(12*l);
        break;
    case 3:
        //transverse concentrate force
        rax=g*(1-xq/l);
        rbx=g-rax;
        break;
    case 4:
        //uniform force
        rbx=g*xq*xq/(2*l);
        rax=g*xq-rbx;
        break;
    case 5:
        ray=g*xq*(2-3*xq*xq/(l*l)+1.6*xq*xq*xq/(l*l*l))/4;
        rby=g*xq/2-ray;
        ma=-g*xq*xq*(2-3*xq/l+1.2*xq*xq/(l*l))/6;
        mb=g*xq*xq*xq*(1-0.8*xq/l)/(4*l);
        break;
    case 6:
        ray=-6*g*(1-xq/l)*xq/(l*l);
        rby=-ray;
        ma=g*(1-xq/l)*(2-3*(1-xq/l));
        mb=g*(xq/l)*(2.0-3*xq/l);
        break;
    case 7:
        // xq =a
        //xq=E*a
        rax=g*xq;
        rbx=-rax;

        break;
    case 8:
        ma=g*xq;
        mb=-ma;
        break;
    default:
        std::cerr<<"error\n"<<std::endl;
    }
    rx[0]=rax;
    rx[1]=rbx;
    ry[0]=ray;
    ry[1]=rby;
    m[0]=ma;
    m[1]=mb;
}
void F(double* f){
    Read rod("input.json");
    Kmatrix K;
    for(unsigned int ild=0;ild<rod.nqq;ild++){
        //traverse all loads
        double rax=0,ray=0,rbx=0,rby=0,ma=0,mb=0;
        double rx[2]={rax,rbx};
        double ry[2]={ray,rby};
        double m[2]={ma,mb};

        unsigned int j=rod.loca_loads[ild]; //by rods
        unsigned int a=rod.ihl[j-1],b=rod.ihr[j-1];
        double l=K.ch(j)[0];
        double cos=K.ch(j)[1];
        double sin=K.ch(j)[2];

        unsigned int ind=rod.category[ild]; //4
        double G=rod.magnitude[ild];
        double dist=rod.dist_load[ild];

        force(ind,rx,ry,m,dist,l,G);
        //局部
        int i0=3*(a-rod.nc-1);
        int j0=3*(b-rod.nc-1);
        if(i0>=0){
            f[i0]+=cos*rx[0]-sin*ry[0];
            f[i0+1]+=sin*rx[0]+cos*ry[0];
            f[i0+2]+=m[0];
        }
        f[j0]+=cos*rx[1]-sin*ry[1];
        f[j0+1]+=sin*rx[1]+cos*ry[1];
        f[j0+2]+=m[1];
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
