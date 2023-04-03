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
//差分格式WENO5
void WENO5p(double *f,double &fm,double eps=1e-6,int p=2){
    int j=2;
    double fpi[3] = {f[j-2]/3.0 - 7.0*f[j-1]/6.0 + 11.0/6.0*f[j],
                -f[j-1]/6.0 + 5.0*f[j]/6.0 + 1.0/3.0*f[j+1],
                f[j]/3.0 + 5.0*f[j+1]/6.0 - 1.0/6.0*f[j+2]};
    double IS[3] = {pow(f[j-2]-4.0*f[j-1]+3*f[j],2)/4.0 + pow(f[j-2]-2*f[j-1]+f[j],2)*13.0/12.0,
                    pow(f[j-1]-f[j+1],2)/4.0 + pow(f[j-1]-2*f[j]+f[j+1],2)*13.0/12.0,
                    pow(f[j]-4.0*f[j+1]+3*f[j+2],2)/4.0 + pow(f[j]-2*f[j+1]+f[j+2],2)*13.0/12.0};
    double C[3] = {1.0/10.0,6.0/10.0,3.0/10.0};
    double alp[3] = {C[0]/pow(eps+IS[0],p),C[1]/pow(eps+IS[1],p),C[2]/pow(eps+IS[2],p)};
    double alps = alp[0] + alp[1] + alp[2];
    double w[3] = {alp[0]/alps,alp[1]/alps,alp[2]/alps};
    fm = w[0]*fpi[0] + w[1]*fpi[1] + w[2]*fpi[2];     
}
void WENO5m(double *f,double &fm,double eps=1e-6,int p=2){
    int j=2;
    double fmi[3] = {f[j+2]/3.0 - 7.0*f[j+1]/6.0 + 11.0/6.0*f[j],
                -f[j+1]/6.0 + 5.0*f[j]/6.0 + 1.0/3.0*f[j-1],
                f[j]/3.0 + 5.0*f[j-1]/6.0 - 1.0/6.0*f[j-2]};
    double IS[3] = {pow(f[j+2]-4.0*f[j+1]+3*f[j],2)/4.0 + pow(f[j+2]-2*f[j+1]+f[j],2)*13.0/12.0,
                    pow(f[j+1]-f[j-1],2)/4.0 + pow(f[j+1]-2*f[j]+f[j-1],2)*13.0/12.0,
                    pow(f[j]-4.0*f[j-1]+3*f[j-2],2)/4.0 + pow(f[j]-2*f[j-1]+f[j-2],2)*13.0/12.0};
    double C[3] = {1.0/10.0,6.0/10.0,3.0/10.0};
    double alp[3] = {C[0]/pow(eps+IS[0],p),C[1]/pow(eps+IS[1],p),C[2]/pow(eps+IS[2],p)};
    double alps = alp[0] + alp[1] + alp[2];
    double w[3] = {alp[0]/alps,alp[1]/alps,alp[2]/alps};
    fm = w[0]*fmi[0] + w[1]*fmi[1] + w[2]*fmi[2];     
}
void Flux_e(double *U,int Nx,double *fo){
    double thisU[3];
    double S[9];
    double invS[9];
    double fhp[3][5];  //特征空间fhat^+
    double fhm[3][5]; //特征空间fhat^-
    double fh[3];//特征空间fhat^+_j+1/2 + fhat^-_j+1/2
    double ft1[3][5];
    double ft2[3][5];
    int i,j,k;
    double fp[3][Nx]; //原空间f^+
    double fm[3][Nx]; //原空间f^-
    double fhpm[3]; //fhat^+_j+1/2
    double fhmm[3]; //fhat^-_j+1/2
    double fom[3]; //原空间 f_j+1/2
    //先求f_j+1/2
    SWf(U,Nx,fp[0],fm[0]); //通量矢量分裂SW
    for(j=2;j<Nx-3;j++){
        for(i=0;i<3;i++) thisU[i] = 0.5*(U[j+i*Nx]+U[j+1+i*Nx]);
        Su(thisU,S); //更新S_j+1/2
        invSu(thisU,invS); //更新S^-1_j+1/2
        
        //正矢量通量
        for(k=j-2;k<=j+2;k++){
            for(i=0;i<3;i++) ft1[i][k-j+2] = fp[i][k];
        }
        mul_matrix(S,ft1[0],ft2[0],3,3,5);  //求fhat^+_k = S_j+1/2*f^+_k
        for(k=0;k<5;k++){
            for(i=0;i<3;i++) fhp[i][k] = ft2[i][k];
        }
        
        //负矢量通量
        for(k=j-1;k<=j+3;k++){
            for(i=0;i<3;i++) ft1[i][k-j+1] = fm[i][k];
        }
        mul_matrix(S,ft1[0],ft2[0],3,3,5);  //求fhat^-_k = S_j+1/2*f^-_k
        for(k=0;k<5;k++){
            for(i=0;i<3;i++) fhm[i][k] = ft2[i][k];
        }

        for(i=0;i<3;i++) {WENO5p(fhp[i],fhpm[i]);WENO5m(fhm[i],fhmm[i]);} //ft1 = fhat^+
        add_matrix(fhpm,fhmm,fh,3,1);  //fhat_j+1/2 = fhat^+_j+1/2+fhat^-_j+1/2
        mul_matrix(invS,fh,fom,3,3,1); //f_j+1/2 = S^-1_j+1/2 * fhat^j+1/2
        for(i=0;i<3;i++) fo[i*(Nx-5)+j-2] = fom[i];
    }
}
void updateU(double *U,int Nx,double dx,double dt,int Nt){
    double fo[3][Nx-5];  //原空间f_j+1/2
    double U1[3][Nx],U2[3][Nx],Unext[3][Nx];
    double cfl = dt/dx;
    int i,j,N;
    for(N = 1;N<=Nt;N++){
        //时间步推进
        Flux_e(U,Nx,fo[0]);
        for(i=0;i<3;i++) {
            U1[i][0] = U[i*Nx];U1[i][1] = U[1+i*Nx];U1[i][2] = U[2+i*Nx];
            U1[i][Nx-3] = U[i*Nx+Nx-3];U1[i][Nx-2] = U[i*Nx+Nx-2];U1[i][Nx-1] = U[i*Nx+Nx-1];
        }
        for(j=3;j<Nx-3;j++){
            for(i=0;i<3;i++){
                U1[i][j] = U[j+i*Nx] - cfl*(fo[i][j-2]-fo[i][j-3]);
            }    
        }
        Flux_e(U1[0],Nx,fo[0]);
        for(i=0;i<3;i++) {
            U2[i][0] = U1[i][0];U2[i][1] = U1[i][1];U2[i][2] = U1[i][2];
            U2[i][Nx-3] = U1[i][Nx-3];U2[i][Nx-2] = U1[i][Nx-2];U2[i][Nx-1] = U1[i][Nx-1];
        }
        for(j=3;j<Nx-3;j++){
            for(i=0;i<3;i++){
                U2[i][j] = 0.75*U[j+i*Nx] + 0.25*(U1[i][j]- cfl*(fo[i][j-2]-fo[i][j-3]));
            }
        }
        Flux_e(U2[0],Nx,fo[0]);
        for(i=0;i<3;i++) {
            Unext[i][0] = U2[i][0];Unext[i][1] = U2[i][1];Unext[i][2] = U2[i][2];
            Unext[i][Nx-3] = U2[i][Nx-3];Unext[i][Nx-2] = U2[i][Nx-2];Unext[i][Nx-1] = U2[i][Nx-1];
        }
        for(j=3;j<Nx-3;j++){
            for(i=0;i<3;i++){
                Unext[i][j] = 1.0*U[j+i*Nx]/3.0 + 2.0/3.0*(U2[i][j]- cfl*(fo[i][j-2]-fo[i][j-3]));
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
    updateU(U[0],Nx,dx,dt,Nt);

    //输出
    ofstream csvfile;
    csvfile.open("SW-weno5_t=0.14.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    system("pause");
    return 0;
}