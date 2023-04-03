#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;
double g = 1.4;
void mul_matrix(double *A,double *B,double *ret,int m=3,int o=3,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = 0.0;
            for(int k=0;k<o;k++){
                ret[n*i+j] += A[i*o+k]*B[k*n+j];
            }
        }
    }
}
void add_matrix(double *A,double *B,double *ret,int m=3,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = A[n*i+j]+B[n*i+j];
        }
    }
}
void print_matrix(double *p,int m=3,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cout << p[i*m+j] << "\t";
        }
        cout << endl;
    }
}
void Su(double *U,double *S){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    S[0] = 0.5*u*u - c*c/(g-1);
    S[1] = -u;
    S[2] = 1.;
    S[3] = -u -(g-1)/c*0.5*u*u;
    S[4] = 1+(g-1)/c*u;
    S[5] = -(g-1)/c;
    S[6] = -u+(g-1)/c*0.5*u*u;
    S[7] = 1-(g-1)/c*u;
    S[8] = (g-1)/c;
}
void invSu(double *U,double *invS){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    invS[0] = -(g-1)/(c*c);
    invS[1] = -1./(2*c);
    invS[2] = 1./(2*c);
    invS[3] = -(g-1)/(c*c)*u;
    invS[4] = -(u-c)/(2*c);
    invS[5] = (u+c)/(2*c);
    invS[6] = -(g-1)/(c*c)*0.5*u*u;
    invS[7] = -1./(2*c)*(h-u*c);
    invS[8] = 1./(2*c)*(h+u*c);
}
inline double lp(double lambda){
    return 0.5*(lambda+sqrt(lambda*lambda + 1e-12));
}
inline double lm(double lambda){
    return 0.5*(lambda-sqrt(lambda*lambda + 1e-12));
}
void fU(double *lbd,double *f,double rho,double u,double c){
    f[0] = rho/(2*g)*(2*(g-1)*lbd[0] + lbd[1] +lbd[2]);
    f[1] = rho/(2*g)*(2*(g-1)*lbd[0]*u + lbd[1]*(u-c)+lbd[2]*(u+c));
    double w = (3-g)*(lbd[1]+lbd[2])*c*c/(2*(g-1));
    f[2] = rho/(2*g)*((g-1)*lbd[0]*u*u + 0.5*lbd[1]*(u-c)*(u-c) + 0.5*lbd[2]*(u+c)*(u+c) + w);
}
void SWf(double *U,int Nx,double *fp,double*fm){
    //通量矢量分裂，SW
    double rho,u,p,h,c,w;
    double lbd[3];
    double lbdp[3];
    double ftmp[3];
    for(int j=0;j<Nx;j++){
        rho = U[j];u = U[j+Nx]/U[j];p =(g-1)*(U[j+Nx*2]-0.5*rho*u*u);
        c = sqrt(g*p/rho);
        lbd[0] = u ; lbd[1] = u-c ;lbd[2]=u+c;
        //lambda+
        for(int i=0;i<3;i++) lbdp[i] = lp(lbd[i]);
        fU(lbdp,ftmp,rho,u,c);
        for(int i=0;i<3;i++) fp[i*Nx+j] = ftmp[i];
        //lambda-
        for(int i=0;i<3;i++) lbdp[i] = lm(lbd[i]);
        fU(lbdp,ftmp,rho,u,c);
        for(int i=0;i<3;i++) fm[i*Nx+j] = ftmp[i];
    }
}
//差分格式f^tilde_j+1/2 = f(...f^tilde....)
void GVC2_p(double *ft,double &ftm){
    int j=1;
    if(abs(ft[j]-ft[j-1]) < abs(ft[j+1]-ft[j])){
        ftm = 0.5*(3*ft[j]-ft[j-1]);
    }
    else{
        ftm = 0.5*(ft[j] + ft[j+1]);
    }
}
void GVC2_m(double *ft,double &ftm){
    int j=1;
    if(abs(ft[j]-ft[j+1]) < abs(ft[j-1]-ft[j])){
        ftm = 0.5*(3*ft[j]-ft[j+1]);
    }
    else{
        ftm = 0.5*(ft[j] + ft[j-1]);
    }
}
void Flux_e(double *U,int Nx,double *fo,int bpn,int sn,int bopn){
    double thisU[3];
    double S[9];
    double invS[9];
    double fhp[3][bpn];  //特征空间fhat^+
    double fhm[3][bpn]; //特征空间fhat^-
    double fh[3];//特征空间fhat^+_j+1/2 + fhat^-_j+1/2
    double ft1[3][bpn];
    double ft2[3][bpn];
    int i,j,k;
    double fp[3][Nx]; //原空间f^+
    double fm[3][Nx]; //原空间f^-
    double fhpm[3]; //fhat^+_j+1/2
    double fhmm[3]; //fhat^-_j+1/2
    double fom[3]; //原空间 f_j+1/2
    //先求f_j+1/2
    SWf(U,Nx,fp[0],fm[0]); //通量矢量分裂SW
    for(j=bopn;j<Nx-1-bopn;j++){
        for(i=0;i<3;i++) thisU[i] = 0.5*(U[j+i*Nx]+U[j+1+i*Nx]);
        Su(thisU,S); //更新S_j+1/2
        invSu(thisU,invS); //更新S^-1_j+1/2
        
        //正矢量通量
        for(k=j-sn;k<=j-sn+bpn-1;k++){
            for(i=0;i<3;i++) ft1[i][k-j+sn] = fp[i][k];
        }
        mul_matrix(S,ft1[0],ft2[0],3,3,bpn);  //求fhat^+_k = S_j+1/2*f^+_k
        for(k=0;k<bpn;k++){
            for(i=0;i<3;i++) fhp[i][k] = ft2[i][k];
        }
        
        //负矢量通量
        for(k=j-(bpn-1-sn)+1;k<=j+sn+1;k++){
            for(i=0;i<3;i++) ft1[i][k-(j-(bpn-1-sn)+1)] = fm[i][k];
        }
        mul_matrix(S,ft1[0],ft2[0],3,3,bpn);  //求fhat^-_k = S_j+1/2*f^-_k
        for(k=0;k<bpn;k++){
            for(i=0;i<3;i++) fhm[i][k] = ft2[i][k];
        }

        for(i=0;i<3;i++) {GVC2_p(fhp[i],fhpm[i]);GVC2_m(fhm[i],fhmm[i]);} //ft1 = fhat^+
        add_matrix(fhpm,fhmm,fh,1,3);  //fhat_j+1/2 = fhat^+_j+1/2+fhat^-_j+1/2
        mul_matrix(invS,fh,fom,3,3,1); //f_j+1/2 = S^-1_j+1/2 * fhat^j+1/2
        for(i=0;i<3;i++) fo[i*(Nx-2*bopn-1)+j-bopn] = fom[i];
    }
}
void updateU(double *U,int Nx,double dx,double dt,int Nt,int bpn,int sn){
    //bpn ： 基架点数量 ； sn（shift num)：格式向左偏移数
    int bopn = (sn>bpn-1-sn ? sn : bpn-1-sn); //边界点数量
    double fo[3][Nx-2*bopn-1];  //原空间f_j+1/2
    double U1[3][Nx],U2[3][Nx],Unext[3][Nx];
    double cfl = dt/dx;
    int i,j,N;
    for(N = 1;N<=Nt;N++){
        //时间步推进
        Flux_e(U,Nx,fo[0],bpn,sn,bopn);
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                U1[i][j] = U[j+i*Nx];
                U1[i][Nx-j-1] = U[i*Nx+Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                U1[i][j] = U[j+i*Nx] - cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]);
            }    
        }
        Flux_e(U1[0],Nx,fo[0],bpn,sn,bopn);
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                U2[i][j] = U1[i][j];
                U2[i][Nx-j-1] = U1[i][Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                U2[i][j] = 0.75*U[j+i*Nx] + 0.25*(U1[i][j]- cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]));
            }
        }
        Flux_e(U2[0],Nx,fo[0],bpn,sn,bopn);
        for(i=0;i<3;i++) {
            for(j=0;j<=bopn;j++) {
                Unext[i][j] = U2[i][j];
                Unext[i][Nx-j-1] = U2[i][Nx-j-1];
            }
        }
        for(j=bopn+1;j<Nx-bopn-1;j++){
            for(i=0;i<3;i++){
                Unext[i][j] = 1.0*U[j+i*Nx]/3.0 + 2.0/3.0*(U2[i][j]- cfl*(fo[i][j-bopn]-fo[i][j-bopn-1]));
            }
        }
        for(i=0;i<3;i++){
            for(j=0;j<Nx;j++){
                U[j+i*Nx] = Unext[i][j];
            }
        }
    }    
}
int main(){
    //开始计算
    int Nx = 101;
    double dx = 1./(Nx-1),dt = 0.001;
    double t_end = 0.14,Nt = round(t_end/dt);
    cout << "Nt=" << Nt << endl;
    double U[3][Nx];
    
    for(int i=0;i<Nx;i++){
        //cout << i*dx << endl;
        U[0][i] = (i>Nx/2 ? 0.125 : 1);  //rho
        U[1][i] = 0.;  // rho*u
        U[2][i] = (i>Nx/2 ? 0.1 : 1)/(g-1); //E = 1/2*rho*u^2 + p/(g-1)
    }
    updateU(U[0],Nx,dx,dt,Nt,3,1);  //使用GVC2格式，3个基架点，第一个点为j-1

    //输出
    ofstream csvfile;
    csvfile.open("General_SW-GVC2_t=0.14.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    //system("pause");
    return 0;
}