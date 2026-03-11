#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <time.h>
#include "trex.h"
#include <string.h>


// The following is to define the variables declard globally in trex.h
// This only needs to happen once in nz.c to avoid linking errors when compiling

const gsl_rng_type* T;
gsl_rng* r;
double ubound, lbound;
int NK;
gsl_matrix_int* Ktable;
double epsilon;
int HL;
gsl_matrix* efl;
gsl_matrix* gc;
gsl_matrix* map;
gsl_matrix* ZZ;
gsl_matrix* iZ;
gsl_matrix* IH;
gsl_matrix* gamma_hat;
gsl_vector* dd;
gsl_vector* loglambda;

double ml, mgc;


double MVNsigmax, MVNsigmay, MVNrho;
double nI2;



int ncoeffs;
double prior;
double anchor;
int nbreak;
double CEIL;

gsl_bspline_workspace* Bw;
gsl_matrix* Bx;
gsl_vector* B;
gsl_vector* Cx;
gsl_matrix* Cov;

double MaxKnot;

void adapt_jump_proposals(gsl_matrix *osu,
                    double acc_b1, double acc_X, double acc_U,
                    double target, int adapt_t_b1, int adapt_t_X, int adapt_t_U)
{
    const double c0   = 0.25;   // gain factor (0.01–0.05 typical)
    const double c1   = 0.1;   // gain factor (0.01–0.05 typical)
    const double c2   = 0.1;   // gain factor (0.01–0.05 typical)
    const double minv = 1e-6;
    const double maxv = 10.0;

    // --- beta1 ---
    {
        double jump     = gsl_matrix_get(osu, 0, 0);
        double log_jump = log(jump);
        double eta      = c0 / sqrt((double)adapt_t_b1);

        log_jump += eta * (acc_b1 - target);
        if (log_jump < log(minv)) log_jump = log(minv);
        if (log_jump > log(maxv)) log_jump = log(maxv);

        gsl_matrix_set(osu, 0, 0, exp(log_jump));
    }

    // --- X ---
    {
        double jump     = gsl_matrix_get(osu, 1, 0);
        double log_jump = log(jump);
        double eta      = c1 / sqrt((double)adapt_t_X);

        log_jump += eta * (acc_X - target);
        if (log_jump < log(minv)) log_jump = log(minv);
        if (log_jump > log(maxv)) log_jump = log(maxv);

        gsl_matrix_set(osu, 1, 0, exp(log_jump));
    }

    // --- U ---
    {
        double jump     = gsl_matrix_get(osu, 2, 0);
        double log_jump = log(jump);
        double eta      = c2 / sqrt((double)adapt_t_U);

        log_jump += eta * (acc_U - target);
        if (log_jump < log(minv)) log_jump = log(minv);
        if (log_jump > log(maxv)) log_jump = log(maxv);

        gsl_matrix_set(osu, 2, 0, exp(log_jump));
    }
}

void dual_averaging_step_size_update(double accept_stat, double target,
                           double *log_eps, double *log_eps_avg,
                           double *Hbar, double mu,
                           double gamma, double t0, double kappa,
                           int t)
{
    // 1. Running average of (target - accept)
    double eta = 1.0 / (t + t0);
    *Hbar = (1 - eta) * (*Hbar) + eta * (target - accept_stat);

    // 2. Update log epsilon
    *log_eps = mu - sqrt((double) t) / gamma * (*Hbar);

    // 3. Smooth log epsilon (Nesterov averaging)
    double w = pow(t, -kappa);
    *log_eps_avg = w * (*log_eps) + (1 - w) * (*log_eps_avg);
}


void adapt_step_num(int *HL, double acc_S)
{
    int hl = *HL;

    if (acc_S < 0.60) {
        hl -= 3;             
    }
    else if (acc_S < 0.80) {
        hl -= 2;             
    }
    else if (acc_S < 0.93) {
        hl -= 1;            
    }
    else if (acc_S < 0.95) {
        hl += 1;            
    }
    else if (acc_S <= 0.96) {
        hl += 2;
    }
    else {
        hl += 3;
    }

    if (hl < 5)  hl = 5;    
    if (hl > 40) hl = 40;    

    *HL = hl;
}



void enforce_identifiability(gsl_matrix *S, int ns) {
    double x1 = gsl_matrix_get(S, 0, 0);
    double y1 = gsl_matrix_get(S, 0, 1);
    double z1 = gsl_matrix_get(S, 0, 2);

    // 1) Translate so first point is (0,0,0)
    for (int i = 0; i < ns; i++) {
        gsl_matrix_set(S, i, 0, gsl_matrix_get(S, i, 0) - x1);
        gsl_matrix_set(S, i, 1, gsl_matrix_get(S, i, 1) - y1);
        gsl_matrix_set(S, i, 2, gsl_matrix_get(S, i, 2) - z1);
    }

    // Point ns-1 is last locus
    double xn = gsl_matrix_get(S, ns-1, 0);
    double yn = gsl_matrix_get(S, ns-1, 1);
    double zn = gsl_matrix_get(S, ns-1, 2);

    // 2) Rotate about z-axis to make yn = 0
    double d = sqrt(xn*xn + yn*yn);
    if (d > 1e-12) {
        double Rz00 = xn/d, Rz01 = yn/d;
        double Rz10 = -yn/d, Rz11 = xn/d;

        for (int i = 0; i < ns; i++) {
            double x = gsl_matrix_get(S, i, 0);
            double y = gsl_matrix_get(S, i, 1);
            double z = gsl_matrix_get(S, i, 2);
            gsl_matrix_set(S, i, 0, Rz00*x + Rz01*y);
            gsl_matrix_set(S, i, 1, Rz10*x + Rz11*y);
            gsl_matrix_set(S, i, 2, z);
        }
    }

    // 3) Rotate about y-axis to make zn = 0 and xn > 0
    xn = gsl_matrix_get(S, ns-1, 0);
    yn = gsl_matrix_get(S, ns-1, 1);
    zn = gsl_matrix_get(S, ns-1, 2);
    d = sqrt(xn*xn + zn*zn);
    if (d > 1e-12) {
        double Ry00 = xn/d, Ry02 = zn/d;
        double Ry20 = -zn/d, Ry22 = xn/d;

        for (int i = 0; i < ns; i++) {
            double x = gsl_matrix_get(S, i, 0);
            double y = gsl_matrix_get(S, i, 1);
            double z = gsl_matrix_get(S, i, 2);
            gsl_matrix_set(S, i, 0, Ry00*x + Ry02*z);
            gsl_matrix_set(S, i, 1, y);
            gsl_matrix_set(S, i, 2, Ry20*x + Ry22*z);
        }
    }

    // 4) Rotate about x-axis to make second point have z > 0
    double x2 = gsl_matrix_get(S, 1, 0);
    double y2 = gsl_matrix_get(S, 1, 1);
    double z2 = gsl_matrix_get(S, 1, 2);
    d = sqrt(y2*y2 + z2*z2);
    if (d > 1e-12) {
        double Rx11 = z2/d, Rx12 = -y2/d;
        double Rx21 = y2/d, Rx22 = z2/d;
        for (int i = 0; i < ns; i++) {
            double x = gsl_matrix_get(S, i, 0);
            double y = gsl_matrix_get(S, i, 1);
            double z = gsl_matrix_get(S, i, 2);
            gsl_matrix_set(S, i, 0, x);
            gsl_matrix_set(S, i, 1, Rx11*y + Rx12*z);
            gsl_matrix_set(S, i, 2, Rx21*y + Rx22*z);
        }
    }

    // 5) Mirror so the y-coordinate of the 3rd point is positive
    double y3 = gsl_matrix_get(S, 2, 1);
    if (y3 < 0) {
        for (int i = 0; i < ns; i++) {
            gsl_matrix_set(S, i, 1, -gsl_matrix_get(S, i, 1));
        }
    }
}



void zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double x=0.0,y=0.0,z=0.0,xj=0.0,yj=0.0,zj=0.0,d=0.0,v=0.0;
    int i,j,k,howmany,kk;
    for(i=0;i<ns;i++){
        x=gsl_matrix_get(S,i,0);
        y=gsl_matrix_get(S,i,1);
        z=gsl_matrix_get(S,i,2);
        howmany=gsl_matrix_int_get(Lookup,i,0);
        //printf("i howmany %d %d\n",i,howmany);
        for(k=1;k<howmany+1;k++){
            //printf("k %d \n",k);
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            if(i>=j) continue;
            //printf("kk %d j %d \n",kk,j);
            xj=gsl_matrix_get(S,j,0);
            yj=gsl_matrix_get(S,j,1);
            zj=gsl_matrix_get(S,j,2);
            d=(x-xj)*(x-xj)+(y-yj)*(y-yj)+(z-zj)*(z-zj);
            v=sqrt(d);
            gsl_matrix_set(Metric,kk,0,v);
            gsl_matrix_set(logMet,kk,0,log(v));
            //printf("kk %d i %d j %d %f\n",kk,i,j,gsl_matrix_get(Metric,kk,0));
        }
    }
    return;
}

double zgetlkI(I,logMet,X,U,beta,ns,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta;
int ns;
gsl_matrix_int *Lookup,*ijTable;
{
    int i,j,k,kk,howmany;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    double retv=0.0,lmet=0.0,iv=0.0,part1=0.0,lambda=0.0;
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0;
    
    
    for(i=0;i<ns;i++){
        howmany=gsl_matrix_int_get(Lookup,i,0);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            if(i>=j) continue;
            
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            part1=b0+b1*lmet+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            lambda=exp(part1);
            
            if(lambda<10.0) retv+=iv*part1-log(exp(lambda)-1.0);
            else retv+=iv*part1-lambda;
        }
    }
    return retv;
}


double zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta;
int ns;
gsl_matrix_int *Lookup,*ijTable;
{
    int i,j,k,kk,howmany;
    double b1=gsl_matrix_get(beta,1,0);
    double xi=0.0,xj=0.0,lmet=0.0,iv=0.0,retv=0.0;
    double part1=0.0,uij=0.0,lambda=0.0;
    
    for(i=0;i<ns;i++){
        howmany=gsl_matrix_int_get(Lookup,i,0);
        xi=gsl_matrix_get(X,i,0);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            if(i>=j) continue;
            
            xj=gsl_matrix_get(X,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            part1=b1*lmet+xi+xj+uij;//printf("part1 %f ",part1);
            
            lambda=exp(part1);
            
            if(lambda<10.0) retv+=iv*part1-log(exp(lambda)-1.0);
            else retv+=iv*part1-lambda;
        }
    }
    return retv;
}

double zgetlkX(X,thetaX,ns,iSigma,beta)
gsl_matrix *X,*thetaX,*iSigma,*beta;
int ns;
{
    double retv,v;
    double sigma=gsl_matrix_get(thetaX,0,0);
    double beta0=gsl_matrix_get(beta,0,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    double efi=0.0,gci=0.0,mapi=0.0;
    double bias=0.0;
    double tXX=0.0;
    
    int i;
    for(i=0;i<ns;i++){
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        bias=b2*efi+b3*gci+b4*mapi;
        v=gsl_matrix_get(X,i,0)-beta0-bias;
        tXX+=v*v;
    }
    
    retv=-0.5*ns*log(sigma)-0.5*tXX/sigma;
    return retv;
}


double zgetlkU(U,thetaX,nI)
gsl_matrix *U,*thetaX;
int nI;
{
    double retv=0.0,vu=0.0;
    
    gsl_matrix *Ut=gsl_matrix_calloc(1,nI);
    gsl_matrix *UtU=gsl_matrix_calloc(1,1);
    
    double sigma=gsl_matrix_get(thetaX,2,0);
    retv=-0.5*nI*log(sigma);
    
    gsl_matrix_transpose_memcpy(Ut,U);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Ut, U,0.0, UtU);
    
    vu=gsl_matrix_get(UtU,0,0);
    retv-=0.5*vu/sigma;
    
    gsl_matrix_free(Ut);
    gsl_matrix_free(UtU);
    
    return retv;
}



double zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric)
gsl_matrix *X,*U,*I,*logMet,*beta,*lk,*osu,*thetaX,*iSigma,*S,*Metric;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    double delta=0.0,deltaI=0.0,deltaX=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0,iv=0.0,v=0.0,uij=0.0;
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double xi=0.0,xj=0.0,xni=0.0,lmet=0.0;
    double rxi=0.0;
    
    int i,j,k,kk,howmany;
    double accepted=0.0;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    double sigma=gsl_matrix_get(thetaX,0,0);
    double osuv=gsl_matrix_get(osu,1,0);
    
    double efi=0.0,gci=0.0,mapi=0.0,bias=0.0;
    
    
    double mx=0.0;
    double check=0.0;
    
    for(i=0;i<ns;i++){
        xi=gsl_matrix_get(X,i,0);
        rxi=gsl_ran_gaussian (r, osuv);
        
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        
        bias=b2*efi+b3*gci+b4*mapi;
        xni=xi+rxi;//Modified on March 10th
        //xni=rxi+b0+bias;//Modified on March 10th
        howmany=gsl_matrix_int_get(Lookup,i,0);
        lkIv=0.0;lk0=0.0;//printf("%d %d\n",i,howmany);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            npart1=b1*lmet+xni+xj+uij;
            part1=b1*lmet+xi+xj+uij;
            
            nlambda=exp(npart1);
            lambda=exp(part1);
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
            
            if(lambda<10.0) lk0+=iv*part1-log(exp(lambda)-1.0);
            else lk0+=iv*part1-lambda;
            
        }
        deltaI=lkIv-lk0;
        

        deltaX=(xni-bias-b0)*(xni-bias-b0)-(xi-bias-b0)*(xi-bias-b0);
        deltaX*=-0.5/sigma;
        delta=deltaI+deltaX;
        if(delta>=0.0){
            gsl_matrix_set(X,i,0,xni);
            v=gsl_matrix_get(lk,0,0);
            v+=deltaI;gsl_matrix_set(lk,0,0,v);
            accepted+=1.0;
            check+=xni;
            mx+=(xni-bias);
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(X,i,0,xni);
                v=gsl_matrix_get(lk,0,0);
                v+=deltaI;gsl_matrix_set(lk,0,0,v);
                accepted+=1.0;
                check+=xni;
                mx+=(xni-bias);
            }else {
                check+=xi;
                mx+=(xi-bias);
            }
        }
        //Rprintf("update %f %f %f X (newX=%f oldX=%f) (deltaI=%f deltaX=%f) accepted %f\n",bias,b0,sigma,xni,xi,deltaI,deltaX,accepted);
    }
    
    mx/=ns;
    check/=ns;
    //printf("%f %f \n",mx,check);
    double C=(mx-b0);
    //C=0.0;
    double scale=exp(2.0*C/b1);
    //scale=1.0;
    double s0,s1,s2;
    //double lkv=getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    //printf("lkv %f %f \n",lkv,gsl_matrix_get(lk,0,0));
    //printf("C scale %f %f %f-> ",C,scale,lkv);
    for(i=0;i<ns;i++){
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        xi=gsl_matrix_get(X,i,0);
        s0*=scale;s1*=scale;s2*=scale;
        gsl_matrix_set(S,i,0,s0);gsl_matrix_set(S,i,1,s1);gsl_matrix_set(S,i,2,s2);
        xi-=C;
        gsl_matrix_set(X,i,0,xi);
    }
    //getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    gsl_matrix_scale (Metric,scale);
    double lscale=log(scale);
    gsl_matrix_add_constant (logMet, lscale);
    //lkv=getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    //printf("%f\n",lkv);
    v=zgetlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,v);
    
    return accepted/ns;
}



double zupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    int i,j,k,kk,howmany;
    double met=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    double iv=0.0,v=0.0,uij=0.0;
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    
    
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0,bias=0.0;
    double accepted=0.0,delta=0.0;
    double osuv=gsl_matrix_get(osu,3,0);
    
    for(i=0;i<ns;i++){
        
        s0=gsl_matrix_get(S,i,0);
        s1=gsl_matrix_get(S,i,1);
        s2=gsl_matrix_get(S,i,2);
        
        sn0=gsl_ran_gaussian (r, osuv);sn0+=s0;
        sn1=gsl_ran_gaussian (r, osuv);sn1+=s1;
        sn2=gsl_ran_gaussian (r, osuv);sn2+=s2;
        
        if(sn0<lbound || sn0>ubound) continue;
        if(sn1<lbound || sn1>ubound) continue;
        if(sn2<lbound || sn2>ubound) continue;
        
        howmany=gsl_matrix_int_get(Lookup,i,0);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        
        lkIv=0.0;lk0=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            
            met=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            nlmet=log(met);
            gsl_matrix_set(Metricnew,k,0,met);
            gsl_matrix_set(logMetnew,k,0,nlmet);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            bias=b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            npart1=b0+b1*nlmet+bias;
            part1=b0+b1*lmet+bias;
            
            nlambda=exp(npart1);
            lambda=exp(part1);
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
            
            if(lambda<10.0) lk0+=iv*part1-log(exp(lambda)-1.0);
            else lk0+=iv*part1-lambda;
            
        }
        delta=lkIv-lk0;
        //printf("%d (%f %f %f) (%f %f %f) %f\n",i,sn0,sn1,sn2,s0,s1,s2,delta);
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            v=gsl_matrix_get(lk,0,0);
            v+=delta;gsl_matrix_set(lk,0,0,v);
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                met=gsl_matrix_get(Metricnew,k,0);
                nlmet=gsl_matrix_get(logMetnew,k,0);
                gsl_matrix_set(Metric,kk,0,met);
                gsl_matrix_set(logMet,kk,0,nlmet);
            }
            accepted+=1.0;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
                v=gsl_matrix_get(lk,0,0);
                v+=delta;gsl_matrix_set(lk,0,0,v);
                for(k=1;k<howmany+1;k++){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    met=gsl_matrix_get(Metricnew,k,0);
                    nlmet=gsl_matrix_get(logMetnew,k,0);
                    gsl_matrix_set(Metric,kk,0,met);
                    gsl_matrix_set(logMet,kk,0,nlmet);
                }
                accepted+=1.0;
            }
        }
        //printf("update S (newS=%f oldS=%f) (deltaI=%f deltaX=%f) accepted %f\n",sn1,s1,deltaI,deltaX,accepted);
    }
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double zupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b1=gsl_matrix_get(beta,1,0);
    int i,j,k,kk,howmany;
    double met=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    double sr0=0.0,sr1=0.0,sr2=0.0;
    
    
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0;
    double accepted=0.0,delta=0.0;
    double osuv=gsl_matrix_get(osu,3,0);
    double xi=0.0,xj=0.0,uij=0.0;
    double v=0.0,iv=0.0;
    
    for(i=0;i<ns;i++){
        
        s0=gsl_matrix_get(S,i,0);
        s1=gsl_matrix_get(S,i,1);
        s2=gsl_matrix_get(S,i,2);
        
        sr0=gsl_ran_gaussian (r, osuv);sn0=s0+sr0;
        sr1=gsl_ran_gaussian (r, osuv);sn1=s1+sr1;
        sr2=gsl_ran_gaussian (r, osuv);sn2=s2+sr2;
        
        if(sn0<lbound || sn0>ubound) continue;
        if(sn1<lbound || sn1>ubound) continue;
        if(sn2<lbound || sn2>ubound) continue;
        
        howmany=gsl_matrix_int_get(Lookup,i,0);
        xi=gsl_matrix_get(X,i,0);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        
        lkIv=0.0;lk0=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            
            met=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            nlmet=log(met);
            gsl_matrix_set(Metricnew,k,0,met);
            gsl_matrix_set(logMetnew,k,0,nlmet);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            
            npart1=b1*nlmet+xi+xj+uij; //+bias;
            part1=b1*lmet+xi+xj+uij; //+bias;
            
            nlambda=exp(npart1);
            lambda=exp(part1);
            
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
            
            if(lambda<10.0) lk0+=iv*part1-log(exp(lambda)-1.0);
            else lk0+=iv*part1-lambda;
            
        }
        delta=lkIv-lk0;
        //printf("%d (%f %f %f) (%f %f %f) %f\n",i,sn0,sn1,sn2,s0,s1,s2,delta);
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            v=gsl_matrix_get(lk,0,0);
            v+=delta;gsl_matrix_set(lk,0,0,v);
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                met=gsl_matrix_get(Metricnew,k,0);
                nlmet=gsl_matrix_get(logMetnew,k,0);
                gsl_matrix_set(Metric,kk,0,met);
                gsl_matrix_set(logMet,kk,0,nlmet);
            }
            accepted+=1.0;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
                v=gsl_matrix_get(lk,0,0);
                v+=delta;gsl_matrix_set(lk,0,0,v);
                for(k=1;k<howmany+1;k++){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    met=gsl_matrix_get(Metricnew,k,0);
                    nlmet=gsl_matrix_get(logMetnew,k,0);
                    gsl_matrix_set(Metric,kk,0,met);
                    gsl_matrix_set(logMet,kk,0,nlmet);
                }
                accepted+=1.0;
            }
        }
        //printf("update S (newS=%f oldS=%f) (deltaI=%f deltaX=%f) accepted %f\n",sn1,s1,deltaI,deltaX,accepted);
    }
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double zHMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    int i,j,k,kk,howmany,jL;
    double met=0.0,met2=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    double factor=1.0;
    
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    
    double accepted=0.0,delta=0.0;
    double iv=0.0;
    
    double pv0=0.0,pv1=0.0,pv2=0.0;
    double dw0=0.0,dw1=0.0,dw2=0.0;
    double current_K=0.0,new_K=0.0;
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0;
    
    for(i=0;i<ns;i++){
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        pv0=gsl_ran_gaussian (r, 1.0);pv1=gsl_ran_gaussian (r, 1.0);pv2=gsl_ran_gaussian (r, 1.0);
        current_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        lk0=0.0;
        dw0=0.0;dw1=0.0;dw2=0.0;
        howmany=gsl_matrix_int_get(Lookup,i,0);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            met2=(s0-sj0)*(s0-sj0)+(s1-sj1)*(s1-sj1)+(s2-sj2)*(s2-sj2);
            lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            part1=b0+b1*lmet+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            
            lambda=exp(part1);
            if(lambda<10.0) {
                lk0+=iv*part1-log(exp(lambda)-1.0);
                factor=exp(lambda)/(exp(lambda)-1.0);
            }else {
                lk0+=iv*part1-lambda;
                factor=1.0;
            }
            
            dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/met2;
            dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/met2;
            dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/met2;
            
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        sn0=s0;sn1=s1;sn2=s2;
        for(jL=1;jL<HL;jL++){
            sn0+=epsilon*pv0;sn1+=epsilon*pv1;sn2+=epsilon*pv2;
            dw0=0.0;dw1=0.0;dw2=0.0;
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                j=gsl_matrix_int_get(ijTable,i,k);
                
                sj0=gsl_matrix_get(S,j,0);
                sj1=gsl_matrix_get(S,j,1);
                sj2=gsl_matrix_get(S,j,2);
                efj=gsl_matrix_get(efl,j,0);
                gcj=gsl_matrix_get(gc,j,0);
                mapj=gsl_matrix_get(map,j,0);
                
                met2=(sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2);
                nlmet=0.5*log(met2);
                iv=gsl_matrix_get(I,kk,0);
                part1=b0+b1*nlmet+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
                lambda=exp(part1);
                if(lambda<10.0) factor=exp(lambda)/(exp(lambda)-1.0);
                else factor=1.0;
                dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/met2;
                dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/met2;
                dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/met2;
                
            }
            pv0+=epsilon*dw0;pv1+=epsilon*dw1;pv2+=epsilon*dw2;
            
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        new_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        
        lkIv=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            met=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            nlmet=log(met);
            gsl_matrix_set(Metricnew,k,0,met);
            gsl_matrix_set(logMetnew,k,0,nlmet);
            
            iv=gsl_matrix_get(I,kk,0);
            npart1=b0+b1*nlmet+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            //lkIv+=iv*npart1-exp(npart1);
            
            nlambda=exp(npart1);
            //lkIv+=iv*npart1-log(exp(nlambda)-1);
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
            
            
        }
        
        delta=lkIv+new_K-lk0-current_K;
        //printf("%d (%f %f %f)<-(%f %f %f) %f\n",i,sn0,sn1,sn2,s0,s1,s2,delta);
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                met=gsl_matrix_get(Metricnew,k,0);
                nlmet=gsl_matrix_get(logMetnew,k,0);
                gsl_matrix_set(Metric,kk,0,met);
                gsl_matrix_set(logMet,kk,0,nlmet);
            }
            accepted+=1.0;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
                for(k=1;k<howmany+1;k++){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    met=gsl_matrix_get(Metricnew,k,0);
                    nlmet=gsl_matrix_get(logMetnew,k,0);
                    gsl_matrix_set(Metric,kk,0,met);
                    gsl_matrix_set(logMet,kk,0,nlmet);
                }
                accepted+=1.0;
            }
        }
    }
    lk0=zgetlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,lk0);
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b1=gsl_matrix_get(beta,1,0);
    
    int i,j,k,kk,howmany,jL;
    double met=0.0,met2=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    double xi=0.0,xj=0.0,uij=0.0;
    
    double factor=1.0;
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    
    double accepted=0.0,delta=0.0;
    double iv=0.0;
    
    
    double pv0=0.0,pv1=0.0,pv2=0.0;
    double dw0=0.0,dw1=0.0,dw2=0.0,dw3=0.0,dw4=0.0;
    double current_K=0.0,new_K=0.0;
    
    
    for(i=0;i<ns;i++){
        
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        xi=gsl_matrix_get(X,i,0);
        pv0=gsl_ran_gaussian (r, 1.0);pv1=gsl_ran_gaussian (r, 1.0);pv2=gsl_ran_gaussian (r, 1.0);
        current_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        lk0=0.0;
        dw0=0.0;dw1=0.0;dw2=0.0;dw3=0.0;dw4=0.0;
        
        howmany=gsl_matrix_int_get(Lookup,i,0);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            
            met2=(s0-sj0)*(s0-sj0)+(s1-sj1)*(s1-sj1)+(s2-sj2)*(s2-sj2);
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            part1=b1*lmet+xi+xj+uij; //+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            
            lambda=exp(part1);
            if(lambda<10.0) {
                lk0+=iv*part1-log(exp(lambda)-1.0);
                factor=exp(lambda)/(exp(lambda)-1.0);
            }else {
                lk0+=iv*part1-lambda;
                factor=1.0;
            }
            dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/met2;
            dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/met2;
            dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/met2;
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        sn0=s0;sn1=s1;sn2=s2;
        for(jL=1;jL<HL;jL++){
            sn0+=epsilon*pv0;sn1+=epsilon*pv1;sn2+=epsilon*pv2;
            dw0=0.0;dw1=0.0;dw2=0.0;dw3=0.0;
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                j=gsl_matrix_int_get(ijTable,i,k);
                xj=gsl_matrix_get(X,j,0);
                sj0=gsl_matrix_get(S,j,0);
                sj1=gsl_matrix_get(S,j,1);
                sj2=gsl_matrix_get(S,j,2);
                
                met2=(sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2);
                nlmet=0.5*log(met2);
                uij=gsl_matrix_get(U,kk,0);
                iv=gsl_matrix_get(I,kk,0);
                part1=b1*nlmet+xi+xj+uij;//+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
                
                lambda=exp(part1);
                if(lambda<10.0) factor=exp(lambda)/(exp(lambda)-1.0);
                else factor=1.0;
                
                dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/met2;
                dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/met2;
                dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/met2;
            }
            pv0+=epsilon*dw0;pv1+=epsilon*dw1;pv2+=epsilon*dw2;
            
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        new_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        
        lkIv=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            sj0=gsl_matrix_get(S,j,0);
            sj1=gsl_matrix_get(S,j,1);
            sj2=gsl_matrix_get(S,j,2);
            
            met=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            nlmet=log(met);
            gsl_matrix_set(Metricnew,k,0,met);
            gsl_matrix_set(logMetnew,k,0,nlmet);
            
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            npart1=b1*nlmet+xi+xj+uij; //+b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            nlambda=exp(npart1);
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
            
        }
        
        delta=lkIv+new_K-lk0-current_K;
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                met=gsl_matrix_get(Metricnew,k,0);
                nlmet=gsl_matrix_get(logMetnew,k,0);
                gsl_matrix_set(Metric,kk,0,met);
                gsl_matrix_set(logMet,kk,0,nlmet);
            }
            accepted+=1.0;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
                for(k=1;k<howmany+1;k++){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    met=gsl_matrix_get(Metricnew,k,0);
                    nlmet=gsl_matrix_get(logMetnew,k,0);
                    gsl_matrix_set(Metric,kk,0,met);
                    gsl_matrix_set(logMet,kk,0,nlmet);
                }
                accepted+=1.0;
            }
        }
    }
    
    lk0=zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,lk0);
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    //printf("accepted %f\n",accepted/ns);
    
    return accepted/ns;
}



double zupdatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b0=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b0=gsl_matrix_get(beta,0,0);
    
    
    osuv=gsl_matrix_get(osu,4,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b0;
    
    gsl_matrix_set(betanew,0,0,betan);
    
    lkv=zgetlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;
    //printf("beta0 lkv delta (%f %f) %f %f\n",b0,betan,lkv,delta);
    
    if(delta>=0){
        gsl_matrix_set(beta,0,0,betan);
        gsl_matrix_set(lk,0,0,lkv);
        accepted=1.0;
    }else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            gsl_matrix_set(beta,0,0,betan);
            gsl_matrix_set(lk,0,0,lkv);
            accepted=1.0;
        }
    }
    gsl_matrix_free(betanew);
    return accepted;
}


double zupdatebetaZero2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu,*thetaX,*iSigma;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double tOneX=0.0;
    double sigma=gsl_matrix_get(thetaX,0,0);
    double rho=gsl_matrix_get(thetaX,1,0);
    double beta0=gsl_matrix_get(beta,0,0);
    
    int i;
    for(i=0;i<ns;i++){
        tOneX+=gsl_matrix_get(X,i,0);
    }
    
    double part0=ns-rho*ns*ns/(1.0-rho+ns*rho);part0/=(1.0-rho);
    double part1=part0/sigma;
    double part2=tOneX-rho*tOneX*ns/(1.0-rho+ns*rho);part2/=(1.0-rho);
    
    part1+=0.001;
    double betan=gsl_ran_gaussian (r, sqrt(1.0/part1));betan+=part2/(sigma*part1);
    gsl_matrix_set(beta,0,0,betan);
    double delta=(betan*betan-beta0*beta0)*part0-2.0*(betan-beta0)*part2;
    delta*=-0.5/sigma;
    double v= gsl_matrix_get(lk,1,0);
    v+=delta;
    //    printf("betazero2 %f %f\n",delta,lkX-v);
    gsl_matrix_set(lk,1,0,v);
    //    printf("%f %f\n",sqrt(1.0/part1),part2/(sigma*part1));
    
    return 0.0;
}


double zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b1=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b1=gsl_matrix_get(beta,1,0);
    
    osuv=gsl_matrix_get(osu,0,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b1;
    if(betan>0.0) return 0.0;
    //if(betan < -1.5) betan=-1.5;
    gsl_matrix_set(betanew,1,0,betan);
    
    
    lkv=zgetlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;//printf("lkv delta (%f %f) %f %f\n",b1,betan,lkv,delta);
    //printf("\nUpdate One (%f %f) %f %f %f\n",b1,betan,lkv,lk0,delta);
    
    if(delta>=0){
        gsl_matrix_set(beta,1,0,betan);
        gsl_matrix_set(lk,0,0,lkv);
        accepted=1.0;
    }else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            gsl_matrix_set(beta,1,0,betan);
            gsl_matrix_set(lk,0,0,lkv);
            accepted=1.0;
        }
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}


double zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b1=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    b1=gsl_matrix_get(beta,1,0);
    osuv=gsl_matrix_get(osu,0,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b1;
    if(betan>0.0) return 0.0;
    //if(betan < -1.5) betan=-1.5;
    
    gsl_matrix_set(betanew,1,0,betan);
    
    
    lkv=zgetlkI2(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;//printf("lkv delta (%f %f) %f %f\n",b1,betan,lkv,delta);
    //printf("\nUpdate One (%f %f) %f %f %f\n",b1,betan,lkv,lk0,delta);
    
    if(delta>=0){
        gsl_matrix_set(beta,1,0,betan);
        gsl_matrix_set(lk,0,0,lkv);
        accepted=1.0;
    }else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            gsl_matrix_set(beta,1,0,betan);
            gsl_matrix_set(lk,0,0,lkv);
            accepted=1.0;
        }
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}


double zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta)
gsl_matrix *X,*lk,*osu,*thetaX,*iSigma,*beta;
int ns,nI;
gsl_matrix_int *Lookup;
{
    double xi=0.0,efi=0.0,gci=0.0,mapi=0.0,bias=0.0;
    int i;
    
    double b0=gsl_matrix_get(beta,0,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    gsl_matrix *Xstar=gsl_matrix_calloc(ns,1);
    gsl_matrix *Xstart=gsl_matrix_calloc(1,ns);
    gsl_matrix *IHX=gsl_matrix_calloc(ns,1);
    gsl_matrix *vx=gsl_matrix_calloc(1,1);
    gsl_matrix *mu=gsl_matrix_calloc(2,1);
    
    for(i=0;i<ns;i++){
        xi=gsl_matrix_get(X,i,0);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        bias=b2*efi+b3*gci+b4*mapi;
        xi-=b0+bias;
        gsl_matrix_set(Xstar,i,0,xi);
        //Rprintf("%f %f %f\n",bias,b0,xi);
    }
    
    gsl_matrix_transpose_memcpy(Xstart,Xstar);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, IH,Xstar,0.0, IHX);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Xstart, IHX,0.0, vx);
    
    double shape=0.001+0.5*(ns-2);
    double SSE=gsl_matrix_get(vx,0,0);
    double scale=0.5*SSE+0.001;
    double iscale=1.0/scale;
    double isigma=gsl_ran_gamma (r, shape, iscale);
    double sigmanew=1.0/isigma;
    //if(sigmanew>0.1 || sigmanew<0.00000001) return 0.0;
    gsl_matrix_set(thetaX,0,0,sigmanew);
    
    //double lkX=-0.5*ns*log(sigmanew)-0.5*SSE/sigmanew;
    double lkX=zgetlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,lkX);
    
    gsl_matrix_free(Xstar);
    gsl_matrix_free(Xstart);
    gsl_matrix_free(mu);
    gsl_matrix_free(IHX);
    gsl_matrix_free(vx);
    return 1.0;
    
}

double zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S)
gsl_matrix *U,*X,*I,*logMet,*beta,*osu,*thetaX,*lk,*S;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    double lknv=0.0,lkv=0.0,deltaI=0.0,deltaU=0.0,delta=0.0;
    double accept=0.0,accepted=0.0;
    double npart1=0.0,part1=0.0,iv=0.0,xi=0.0,xj=0.0;
    double uij=0.0,unij=0.0,lmet=0.0,v=0.0,nlambda=0.0,lambda=0.0;
    
    double sigmau=gsl_matrix_get(thetaX,2,0);
    double osuv=gsl_matrix_get(osu,2,0);
    double b1=gsl_matrix_get(beta,1,0);
    //double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0;
    
    int i,j,k,howmany,kk;
    k=0;
    //double mv=0.0;
    for(i=0;i<ns;i++){
        howmany=gsl_matrix_int_get(Lookup,i,0);
        xi=gsl_matrix_get(X,i,0);
        //efi=gsl_matrix_get(efl,i,0);
        //gci=gsl_matrix_get(gc,i,0);
        //mapi=gsl_matrix_get(map,i,0);
        lknv=0.0;
        lkv=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            if(i>=j) continue;
            
            xj=gsl_matrix_get(X,j,0);
            //efj=gsl_matrix_get(efl,j,0);
            //gcj=gsl_matrix_get(gc,j,0);
            //mapj=gsl_matrix_get(map,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            unij=gsl_ran_gaussian (r, osuv);
            unij+=uij;
            
            npart1=b1*lmet+xi+xj+unij; //+bias;
            part1=b1*lmet+xi+xj+uij; //+bias;
            
            nlambda=exp(npart1);
            lambda=exp(part1);
            
            if(nlambda<10.0) lknv+=iv*npart1-log(exp(nlambda)-1.0);
            else lknv+=iv*npart1-nlambda;
            
            if(lambda<10.0) lkv+=iv*part1-log(exp(lambda)-1.0);
            else lkv+=iv*part1-lambda;
            
            
            deltaI=lknv-lkv;
            deltaU=0.5/sigmau*(uij*uij-unij*unij);
            delta=deltaI+deltaU;
            
            //Rprintf("delta %f",exp(delta));
            //Rprintf("deltaI U %f %f \n",deltaI,deltaU);
            
            if(delta>=0.0){
                gsl_matrix_set(U,kk,0,unij);
                v=gsl_matrix_get(lk,0,0);
                v+=deltaI;gsl_matrix_set(lk,0,0,v);
                v=gsl_matrix_get(lk,2,0);
                v+=deltaU;gsl_matrix_set(lk,2,0,v);
                accepted+=1.0;
                //mv+=unij;
            }
            else{
                accept=gsl_rng_uniform (r);
                if(accept<exp(delta)){
                    gsl_matrix_set(U,kk,0,unij);
                    v=gsl_matrix_get(lk,0,0);
                    v+=deltaI;gsl_matrix_set(lk,0,0,v);
                    v=gsl_matrix_get(lk,2,0);
                    v+=deltaU;gsl_matrix_set(lk,2,0,v);
                    accepted+=1.0;
                    //mv+=unij;
                }//else mv+=uij;
            }
            
        }
    }
/*
    mv/=nI;
    double vu=0.0;
    for(i=0;i<nI;i++){
        v=gsl_matrix_get(U,i,0);
        v-=mv;
        gsl_matrix_set(U,i,0,v);
        vu+=v*v;
    }

    v=gsl_matrix_get(lk,2,0);
    v+=0.5*mv*mv*nI2/sigmau;

    
    Rprintf("updateU1 %f ",v);
    v=-0.5*nI*log(sigmau);
    v-=0.5*vu/sigmau;
    gsl_matrix_set(lk,2,0,v);
    Rprintf("updateU2 %f ",v);
    
    v=exp(mv/b1);
    Rprintf("adjust %f %f %f %f\n",v,mv,b1,sigmau);
    
    gsl_matrix_scale (S, v);

*/
    //double vu=gsl_matrix_get(lk,2,0);
    //v=zgetlkU(U,thetaX,nI);
    //gsl_matrix_set(lk,2,0,v);
    
    //Rprintf("updateU %f right %f \n",vu,v);
 
    return accepted/nI;
}

double zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI)
gsl_matrix *U,*lk,*osu,*thetaX;
int ns,nI;
gsl_matrix_int *Lookup;
{
    double vu=0.0,lkU=0.0;
    double sigma=gsl_matrix_get(thetaX,2,0);
    
    vu=gsl_matrix_get(lk,2,0);
    vu+=nI*0.5*log(sigma);
    vu*=-2.0*sigma;

/*
    double mv=0.0;int i=0;
    double v=0.0;
    for(i=0;i<nI;i++){
	v=gsl_matrix_get(U,i,0);
	mv+=v*v;
    } 
    mv/=nI2;
*/
    double shape=0.001+0.5*nI;
    double scale=0.5*vu+0.001;
    double iscale=1.0/scale;
    double isigma=gsl_ran_gamma (r, shape, iscale);

    double sigmanew=1.0/isigma;
    
    lkU=gsl_matrix_get(lk,2,0);
    //Rprintf("raw estimate iscale isigma %f %f %f %f curr %f",shape,scale,scale/(shape+1.0),sigmanew,lkU);
    //	if(sigmanew>0.5 || sigmanew<0.00000001) return 0.0;
    gsl_matrix_set(thetaX,2,0,sigmanew);
    
    //vu+=mv*mv/nI;
    lkU=-0.5*nI*log(sigmanew);
    lkU-=0.5*vu/sigmanew;
    gsl_matrix_set(lk,2,0,lkU);
    //Rprintf(" new %f\n",lkU);
    return 0.0;
}



double zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu,*thetaX,*iSigma;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    int i;
    double xi=0.0,mapi=0.0;
    double b0=gsl_matrix_get(beta,0,0);
    double b4=gsl_matrix_get(beta,4,0);
    double sigma=gsl_matrix_get(thetaX,0,0);
    
    
    gsl_matrix *Xstar=gsl_matrix_calloc(ns,1);
    gsl_matrix *mu=gsl_matrix_calloc(2,1);
    
    for(i=0;i<ns;i++){
        xi=gsl_matrix_get(X,i,0);
        mapi=gsl_matrix_get(map,i,0);
        xi-=b0+b4*mapi;
        gsl_matrix_set(Xstar,i,0,xi);
    }
    
    double vx,vy;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, gamma_hat, Xstar,0.0, mu);
    //printf("ran_bivariate %f %f\n",gsl_matrix_get(mu,0,0),gsl_matrix_get(mu,1,0));
    //printf("ran_bivariate %f %f %f\n",MVNsigmax,MVNsigmay,MVNrho);
    double usigmax=MVNsigmax*sqrt(sigma);
    double usigmay=MVNsigmay*sqrt(sigma);
    //double urho=MVNrho*sigma/(usigmax*usigmay);
    gsl_ran_bivariate_gaussian (r, usigmax,usigmay,MVNrho,&vx,&vy);
    
    //printf("vx,vy %f %f \n",usigmax,usigmay);
    
    
    gsl_matrix_set(beta,2,0,vx);
    gsl_matrix_set(beta,3,0,vy);
    
    double lkX=zgetlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,lkX);
    
    gsl_matrix_free(Xstar);
    gsl_matrix_free(mu);
    
    return 1.0;
}


double zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b2=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    
    b2=gsl_matrix_get(beta,2,0);
    osuv=gsl_matrix_get(osu,5,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b2;
    gsl_matrix_set(betanew,2,0,betan);
    
    lkv=zgetlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;
    double priordelta=0.5*b2*b2*prior-0.5*betan*betan*prior;
    delta+=priordelta;
    
    if(delta>=0){
        gsl_matrix_set(beta,2,0,betan);
        gsl_matrix_set(lk,0,0,lkv);
        accepted=1.0;
    }else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            gsl_matrix_set(beta,2,0,betan);
            gsl_matrix_set(lk,0,0,lkv);
            accepted=1.0;
        }
        else gsl_matrix_set(lk,0,0,lk0);
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}


double zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b3=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b3=gsl_matrix_get(beta,3,0);
    osuv=gsl_matrix_get(osu,6,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b3;
    gsl_matrix_set(betanew,3,0,betan);
    
    lkv=zgetlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;
    double priordelta=0.5*b3*b3*prior-0.5*betan*betan*prior;
    delta+=priordelta;
    
    if(delta>=0){
        gsl_matrix_set(beta,3,0,betan);
        gsl_matrix_set(lk,0,0,lkv);
        accepted=1.0;
    }else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            gsl_matrix_set(beta,3,0,betan);
            gsl_matrix_set(lk,0,0,lkv);
            accepted=1.0;
        }
        else gsl_matrix_set(lk,0,0,lk0);
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}



void nz(
        double *Contact,
        double *bias,
        int *n, 
        int *repn, 
        int *repb, 
        double *argv0, 
        double *argv1,
        double *argv2,
        int *nHL,
        double *vepsilon,
        int *thinning,
        int *gear,
        double *result)
{
    int ns=*n;
    int nI=0;
    nI2=ns*(ns-1)*0.5;
    int nu2r=(int)ns*(ns-1)*0.5;
    gsl_matrix *U2R = gsl_matrix_calloc(1,nu2r);
    
    int rep_b=*repb;
    int rep_n=*repn;
    
    int i,j,k,howmany,ii,jj,kk,uk;
    double v,Mi,line_t;
    double max_count;
    
    prior=0.01;
    
    gsl_matrix_int *Lookup=gsl_matrix_int_alloc(ns,ns);
    gsl_matrix_int *ijTable=gsl_matrix_int_alloc(ns,ns);
    
    gsl_matrix *Jc=gsl_matrix_calloc(ns,ns);
    gsl_matrix *X=gsl_matrix_calloc(ns,1);
    gsl_matrix *S=gsl_matrix_calloc(ns,3);
    efl=gsl_matrix_calloc(ns,1);
    gc=gsl_matrix_calloc(ns,1);
    map=gsl_matrix_calloc(ns,1);
    
    gsl_matrix *Imat=gsl_matrix_calloc(ns,ns);
    gsl_matrix *Jmat=gsl_matrix_calloc(ns,ns);
    gsl_matrix *iSigma=gsl_matrix_calloc(ns,ns);
    gsl_matrix *beta=gsl_matrix_calloc(5,1);
    gsl_matrix *thetaX=gsl_matrix_calloc(3,1);
    gsl_matrix *osu=gsl_matrix_calloc(7,1);
    gsl_matrix *lk=gsl_matrix_calloc(3,1);
    
    gsl_matrix *Z=gsl_matrix_calloc(ns,2);
    gsl_matrix *Zt=gsl_matrix_calloc(2,ns);
    gsl_matrix *Zgamma=gsl_matrix_calloc(ns,ns);
    
    gsl_vector *collectorXs=gsl_vector_calloc(ns);
    gsl_vector *collectorb=gsl_vector_calloc(5);
    gsl_vector *collectorthetaXs=gsl_vector_calloc(3);
    gsl_vector *collectorLK=gsl_vector_calloc(3);
    
    int  filter_factor=*thinning;
    int  filter_n=rep_n/filter_factor;
    int ufilter_factor=filter_factor*10;
    int ufilter_n=rep_n/ufilter_factor;
    
	gsl_matrix *accept=gsl_matrix_calloc(filter_n,4);
	gsl_matrix *collectorS=gsl_matrix_calloc(filter_n,3*ns);
	gsl_matrix *collectorX=gsl_matrix_calloc(filter_n,ns);
	gsl_matrix *collectorbeta=gsl_matrix_calloc(filter_n,5);
	gsl_matrix *collectorthetaX=gsl_matrix_calloc(filter_n,3);
	gsl_matrix *collectorlk=gsl_matrix_calloc(filter_n,3);
    
	gsl_matrix_set(osu,0,0,*argv0); //betaOne
	gsl_matrix_set(osu,1,0,*argv1); //X
	gsl_matrix_set(osu,2,0,*argv2); //U

	gsl_matrix_set(osu,3,0,0.07); //S
	gsl_matrix_set(osu,4,0,0.03); //betaZero

	gsl_matrix_set(osu,5,0,0.03); //EFL
	gsl_matrix_set(osu,6,0,0.03); //GC
    
    if(*gear == 0){
        gsl_matrix_set(osu,5,0,*argv1);
        gsl_matrix_set(osu,6,0,*argv2);
    }

	
	gsl_rng_env_setup();
 	T = gsl_rng_default;
 	r = gsl_rng_alloc (T);
 	gsl_rng_set (r, 3);
 	
	double initial=gsl_rng_uniform (r);

	ubound=1.0; //atof(argv[29]);		  
	lbound=-1.0*ubound;
	
	HL=*nHL;
	epsilon=*vepsilon;
	double log_eps = log(epsilon);

	k=0;
	for(i=0;i<ns;i++) {
		for(j=i+1;j<ns;j++){
			v=Contact[k];
			if(v != 0 ) nI++;
			gsl_matrix_set(Jc,i,j,v);
	  		gsl_matrix_set(Jc,j,i,v);
			k++;
		}
	}
			
	
	for(i=0;i<ns;i++){
	   gsl_matrix_int_set(Lookup,i,0,0);
	}
	//printf("Number of non-zero counts %d \n",nI);
	
	gsl_matrix *I=gsl_matrix_calloc(nI,1);
	gsl_matrix *U=gsl_matrix_calloc(nI,1);

	gsl_matrix *Metric=gsl_matrix_calloc(nI,1);
	gsl_matrix *logMet=gsl_matrix_calloc(nI,1);

    gsl_vector *collectorUs=gsl_vector_calloc(nI);
    gsl_matrix *collectorU=gsl_matrix_calloc(ufilter_n,nI);

	max_count=0;
	
	k=0;Mi=0.0;
	line_t=ubound*initial;


	ZZ=gsl_matrix_calloc(2,2);
	IH=gsl_matrix_calloc(ns,ns);
	gamma_hat=gsl_matrix_calloc(2,ns);
	iZ=gsl_matrix_calloc(2,2);

	k=0;
	for(i=0;i<ns;i++){
		v=bias[k];
		gsl_matrix_set(efl,i,0,v);gsl_matrix_set(Z,i,0,v);
		k++;v=bias[k];
		gsl_matrix_set(gc,i,0,v);gsl_matrix_set(Z,i,1,v);

		k++;v=bias[k];
		gsl_matrix_set(map,i,0,v);
	}

	k=0;max_count=0.0;
	for(i=0;i<ns;i++){	
		for(j=i+1;j<ns;j++){
			v=gsl_matrix_get(Jc,i,j);
			if(max_count<v) max_count=v;
			if(v != 0.0){
				Mi+=v;
				howmany=gsl_matrix_int_get(Lookup,i,0)+1;
				gsl_matrix_int_set(Lookup,i,howmany,k);
				gsl_matrix_int_set(ijTable,i,howmany,j);
				gsl_matrix_int_set(Lookup,i,0,howmany);
				gsl_matrix_int_set(ijTable,i,0,howmany);
				//printf("%d %d %d %d \n",i,j,howmany,k);
				howmany=gsl_matrix_int_get(Lookup,j,0)+1;
				gsl_matrix_int_set(Lookup,j,howmany,k);
				gsl_matrix_int_set(ijTable,j,howmany,i);
				gsl_matrix_int_set(Lookup,j,0,howmany);
				gsl_matrix_int_set(ijTable,j,0,howmany);
				
				gsl_matrix_set(I,k,0,v);
				//printf("%d %d %d %d %f\n",i,j,howmany,k,v);
				k++;
			}
		}	
		gsl_matrix_set(S,i,0,line_t*i/(ns-1));
		gsl_matrix_set(S,i,1,0.0);
		gsl_matrix_set(S,i,2,0.0);
        gsl_matrix_set(X,i,0,0.0);
	}

    
	gsl_matrix_transpose_memcpy(Zt,Z);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Zt, Z,0.0, ZZ);

	double a,b,c,d;
	a=gsl_matrix_get(ZZ,0,0);
	b=gsl_matrix_get(ZZ,0,1);
	c=gsl_matrix_get(ZZ,1,0);
	d=gsl_matrix_get(ZZ,1,1);

	gsl_matrix_set(iZ,0,0,d/(a*d-b*c));
	gsl_matrix_set(iZ,0,1,-1.0*c/(a*d-b*c));
	gsl_matrix_set(iZ,1,0,-1.0*b/(a*d-b*c));
	gsl_matrix_set(iZ,1,1,a/(a*d-b*c));
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, iZ, Zt,0.0, gamma_hat);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Z, gamma_hat,0.0, Zgamma);
	gsl_matrix_set_identity (IH);
	gsl_matrix_sub (IH, Zgamma);
	MVNsigmax=sqrt(d/(a*d-b*c));
	MVNsigmay=sqrt(a/(a*d-b*c));
	MVNrho= -1.0*c/(a*d-b*c)/(MVNsigmax*MVNsigmay);
	
    enforce_identifiability(S,ns);

	zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
	
	gsl_matrix_set(beta,0,0,initial*log(max_count));
	gsl_matrix_set(beta,1,0,-2.0*initial);
	
	//printf("Initial beta0 beta1 %f %f X S %f\n",initial*log(max_count),-2.0*initial,line_t);
	
	gsl_matrix_set(thetaX,0,0,0.1);
	gsl_matrix_set(thetaX,1,0,0.0);
	gsl_matrix_set(thetaX,2,0,0.1);
	
	gsl_matrix_set_identity(iSigma);
	gsl_matrix_set_identity(Imat);
	gsl_matrix_set_all(Jmat,1.0); 
	
	gsl_matrix *tOne=gsl_matrix_calloc(1,ns);
	gsl_matrix *One=gsl_matrix_calloc(ns,1);
	gsl_matrix_set_all(tOne,1.0);
	gsl_matrix_set_all(One,1.0);
	
	gsl_matrix_set_all(Jmat,1.0); 
	
	gsl_matrix_set(lk,0,0,0.0);
	gsl_matrix_set(lk,1,0,0.0);
	gsl_matrix_set(lk,2,0,0.0);
	
	double sv=0.0;
	double s1,s2,s3;
	
	double diff = 0.0;
    time_t start;
    time_t stop;
    time(&start);



    
    double v0,v1,v2,v3,mv;
    double scale,b1;
    int fi=0;int ufi=0;
    double adjust=0.1;
    
    int nn=0;
    int howoften=0;


	v=gsl_matrix_get(beta,0,0);
	gsl_matrix_set(beta,0,0,0.5*v);
	for(i=0;i<ns;i++){
		line_t=gsl_ran_gaussian (r, 0.1*ubound);
        gsl_matrix_set(X,i,0,0.0);
	}
//===========================================
    
	v=gsl_matrix_get(beta,0,0);

	gsl_matrix_set(beta,1,0,0.0);
	gsl_matrix_set(beta,0,0,0.5*v);
	sv= zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
	gsl_matrix_set(lk,0,0,sv);
	sv= zgetlkX(X,thetaX,ns,iSigma,beta);
	gsl_matrix_set(lk,1,0,sv);
	sv= zgetlkU(U,thetaX,nI);
    gsl_matrix_set(lk,2,0,sv);
    
    gsl_matrix_set(beta,4,0,1.0);//printf("MAP is set off\n");
    	//printf("%f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
//	printf("INITIALIZATION \n");
    //=====================================
    gsl_matrix_set(beta,1,0,0.0);
	for(i=0;i<10000;i++){
		zupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		zupdatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
	}

    
    v=gsl_matrix_get(beta,0,0);
    v*=0.5;
	gsl_matrix_set(beta,0,0,1.5);
    b1=gsl_matrix_get(beta,1,0);
	scale=exp(2.0*(v-1.5)/b1);
    
	for(i=0;i<ns;i++){
	 	s1=gsl_matrix_get(S,i,0);s2=gsl_matrix_get(S,i,1);s3=gsl_matrix_get(S,i,2);
	 	s1*=scale;s2*=scale;s3*=scale;
	 	gsl_matrix_set(S,i,0,s1);gsl_matrix_set(S,i,1,s2);gsl_matrix_set(S,i,2,s3);
        gsl_matrix_set(X,i,0,1.5);
	}
    
    enforce_identifiability(S,ns);
    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    v=zgetlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,v);
    v=zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);

    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
    int iter1 = 5000;
    int iter2 = 2000;
    int iter3 = 5000;
    int iter4 = 2000;
    int iter5 = 10000;
    
    int block1_size = 500;
    int block2_size = 500;
    int block3_size = 500;
    int block4_size = 500;
    int block5_size = 500;
    
    double kappa_init = 0.5; // learning rate for averaging
    
    double lr_mass = 0.1;   // learning rate for EMA (tune 0.05–0.2) in cov update
    
    
    double mu;        // log of the initial epsilon * 10 (Stan uses log(10*epsilon_init))
    double log_eps_avg; // running average of log eps (smoothed)
    double Hbar;      // running average of (target - accept_stat)
    double gamma;     // controls shrinkage (Stan uses 0.05)
    double t0;        // stabilizes early iterations (Stan uses 10)
    int adapt_iter;   // counter for adaptation iterations
    
    gamma = 0.05;
    t0 = 10.0;
    
    // initialize dual-averaging
    mu = log(10*epsilon);   // Stan uses log(10 * epsilon_init)
    log_eps_avg = log_eps;
    Hbar = 0.0;
    adapt_iter = 1;
    
    int adapt_t_b1 = 1;
    int adapt_t_X  = 1;
    int adapt_t_U  = 1;
    
    Rprintf("\n====================================\n");
    int check=0;
if(*gear == 1){
	Rprintf("Jumping rules in Proposals\n");
	Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
	Rprintf("beta1 : osu %f \n",gsl_matrix_get(osu,0,0));
	Rprintf("X : osu %f \n",gsl_matrix_get(osu,1,0));
	Rprintf("U : osu %f \n",gsl_matrix_get(osu,2,0));
    Rprintf("Looking for jumping rules...\n");
    

    int max_cycle = 10;
        
        for (int cycle = 0; cycle < max_cycle; cycle++) {
            
            // Parameters for adaptation
            double target = 0.95;    // target acceptance rate
            double kappa = kappa_init;    // learning rate
            double kappa_decay = 0.95; //learning rate decay
            
            // reset log-epsilon?
            
            // log_eps = log(epsilon_init);
            // epsilon  = epsilon_init;
            
            Rprintf("\n==================== Adaptation Cycle %d ====================\n", cycle+1);
        
            /* -------------------- Phase 1: Step-size adaptation -------------------- */
            Rprintf("Phase 1: Step-size adaptation\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            int n_iter = iter1; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            int block_size = block1_size;
            double acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            double acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;

            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                
                // every 'block_size' iterations, update epsilon and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    dual_averaging_step_size_update(acc_S, target,
                      &log_eps, &log_eps_avg,
                      &Hbar, mu,
                      gamma, t0, kappa,
                      adapt_iter);

                    epsilon = exp(log_eps);
                    adapt_iter++;
                    
                    adapt_step_num(&HL, acc_S);
                    
                    adapt_jump_proposals(osu, acc_b1, acc_X, acc_U, 0.25, adapt_t_b1, adapt_t_X, adapt_t_U);
                    
                    adapt_t_b1 += 1;
                    adapt_t_X += 1;
                    adapt_t_U += 1;
                    
                    // Report diagnostics
                    Rprintf("\n[Phase 1 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
            }
            
            /* ---------- Phase summary ---------- */
            double mean_acc_S  = acc_phase_S  / n_iter;
            double mean_acc_b1 = acc_phase_b1 / n_iter;
            double mean_acc_X  = acc_phase_X  / n_iter;
            double mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 1 completed | Cycle %d ===\n", cycle + 1);
            
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
        
        
            /* -------------------- Phase 2: Step-size + Mass adaptation -------------------- */
            Rprintf("Phase 2: Mass + Step-size adaptation\n");
            
            kappa = kappa_init;
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter2; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block2_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
        
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                // every 'block_size' iterations, update epsilon and mass matrix and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    // Report diagnostics
                    Rprintf("\n[Phase 2 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
            }
            
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 2 completed | Cycle %d ===\n", cycle + 1);
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 3: Freeze mass, stabilize epsilon -------------------- */
            Rprintf("Phase 3: Freeze Mass, tune step size \n");
            
            kappa = kappa_init;
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter3; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block3_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                // every 'block_size' iterations, update epsilon and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    // Step-size adaptation for S
                    // if (acc_S < target){
                    //     adapt_step_size(&epsilon, &log_eps, acc_S, target, &kappa, kappa_decay);
                    //     adapt_step_num(&HL, acc_S);
                    // }
                    
                    dual_averaging_step_size_update(acc_S, target,
                      &log_eps, &log_eps_avg,
                      &Hbar, mu,
                      gamma, t0, kappa,
                      adapt_iter);

                    epsilon = exp(log_eps);
                    adapt_iter++;
                    
            
                    // Report diagnostics
                    Rprintf("\n[Phase 3 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 3 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 4: Stability run (no adaptation) -------------------- */
            Rprintf("Phase 4: Stability run\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter4; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block4_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                
                // every 'block_size' iterations report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
                    
                    
            
                    // Report diagnostics
                    Rprintf("\n[Phase 4 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 4 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 5: Tune jump proposals -------------------- */
            Rprintf("Phase 5: Tuning beta1, X, and U proposals\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter5; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block5_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
            
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                // every 'block_size' iterations report acceptance and adapt jump proposals
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    // === Adapt jump proposals ===
                    adapt_jump_proposals(osu, acc_b1, acc_X, acc_U, 0.25, adapt_t_b1, adapt_t_X, adapt_t_U);
                    
                    adapt_t_b1 += 1;
                    adapt_t_X += 1;
                    adapt_t_U += 1;
            
                    // Report diagnostics
                    Rprintf("\n[Phase 5 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 5 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
                    
            /* -------------------- Early stopping check -------------------- */
            int good_S  = (mean_acc_S  > 0.80 && mean_acc_S  < 0.85);
            int good_b1 = (mean_acc_b1 > 0.20 && mean_acc_b1 < 0.35);
            int good_X  = (mean_acc_X  > 0.20 && mean_acc_X  < 0.35);
            int good_U  = (mean_acc_U  > 0.20 && mean_acc_U  < 0.35);
            
            if (good_S && good_b1 && good_X && good_U) {
                Rprintf("\n*** Early stopping: acceptance rates are good at Cycle %d ***\n",
                        cycle + 1);
            
                // free matrices before break
                if (accept != NULL) {
                    gsl_matrix_free(accept);
                    accept = NULL;
                }
            
                // break out of cycle loop and skip resetting structure
                break;
            }
        
            if(cycle < max_cycle-1){
                /* -------------------- Reset between cycles -------------------- */
                Rprintf("Reset structure for next cycle\n");
                
                    // ---- Phase accumulators and counters ----
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                    acc_phase_S = acc_phase_b1 = acc_phase_X = acc_phase_U = 0.0;
                    fi = 0;
                    ufi = 0;
                    
                    Hbar        = 0.0;           // running average of (target - accept)
                    log_eps_avg = 0.0;           // smoothed average of log-epsilon
                    adapt_iter  = 1;             // dual-averaging iteration counter
                    mu = log(10.0 * *vepsilon);
                    adapt_t_b1 = 1;
                    adapt_t_X  = 1;
                    adapt_t_U  = 1;
                    
                    // ---- Free any leftover matrices from previous phases ----
                    if (accept != NULL) {
                        gsl_matrix_free(accept);
                        accept = NULL;
                    }
                
                    gsl_matrix_set(beta,0,0,initial*log(max_count));
                    v=gsl_matrix_get(beta,0,0);
                    for(i=0;i<ns;i++){
                        gsl_matrix_set(S,i,0,line_t*i/(ns-1));
                        gsl_matrix_set(S,i,1,0.0);
                        gsl_matrix_set(S,i,2,0.0);
                        //line_t=gsl_ran_gaussian (r, 0.1*ubound);
    
                    }
                    
                    gsl_matrix_set(thetaX,0,0,0.1); 
                    gsl_matrix_set(thetaX,1,0,0.0);
                    gsl_matrix_set(thetaX,2,0,0.1);
                    
                    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                    gsl_matrix_set(beta,1,0,-2.0*initial);
                    sv= zgetlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
                    gsl_matrix_set(lk,0,0,sv);
                    for(i=0;i<10000;i++){
                        zupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                        zupdatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                        zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                    }
                    
                    
                    v=gsl_matrix_get(beta,0,0);
                    v*=0.5;
                    gsl_matrix_set(beta,0,0,1.5);
                    b1=gsl_matrix_get(beta,1,0);
                    scale=exp(2.0*(v-1.5)/b1);
                    
                    for(i=0;i<ns;i++){
                        s1=gsl_matrix_get(S,i,0);s2=gsl_matrix_get(S,i,1);s3=gsl_matrix_get(S,i,2);
                        s1*=scale;s2*=scale;s3*=scale;
                        gsl_matrix_set(S,i,0,s1);gsl_matrix_set(S,i,1,s2);gsl_matrix_set(S,i,2,s3);
                        gsl_matrix_set(X,i,0,1.5);
                    }
                    
                    
                    enforce_identifiability(S, ns);
                    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                    
                    sv= zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
                    gsl_matrix_set(lk,0,0,sv);
                    sv= zgetlkX(X,thetaX,ns,iSigma,beta);
                    gsl_matrix_set(lk,1,0,sv);
                    sv= zgetlkU(U,thetaX,nI);
                    gsl_matrix_set(lk,2,0,sv);
                    gsl_matrix_set(beta,4,0,1.0);
            
            }
        }
    
    
    accept = gsl_matrix_calloc(rep_b, 4);
        
    double acc_rate = 0.0;
    
    fi=0;ufi=0;
	Rprintf("\nBurning has started\n");
	for(i=0;i<rep_b;i++){
		v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
		v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
		v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		
		enforce_identifiability(S, ns);
        zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
		
		v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
		v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
		v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
		
		gsl_matrix_set(accept,fi,0,v0);
        gsl_matrix_set(accept,fi,1,v1);
        gsl_matrix_set(accept,fi,2,v2);
        gsl_matrix_set(accept,fi,3,v3);
            
        fi += 1;    
		
		if ((i+1)==rep_b) {
			Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
			Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
			//Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0)+gsl_matrix_get(lk,1,0)+gsl_matrix_get(lk,2,0));
		}
	}

	Rprintf("Burning has ended\n");
	
	for(j=0;j<4;j++) {
            acc_rate=0.0;
            for(i=0;i<rep_b;i++){
                acc_rate+=gsl_matrix_get(accept,i,j);
            }
            acc_rate/=rep_b;
            result[k]=acc_rate;k++;
            if(j==0) Rprintf("Burn-in S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL, epsilon, acc_rate);
            if(j==1) Rprintf("Burn-in beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),acc_rate);
            if(j==2) Rprintf("Burn-in X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),acc_rate);
            if(j==3) Rprintf("Burn-in U : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),acc_rate);
    }
        
        
    if (accept != NULL) {
                    gsl_matrix_free(accept);
                    accept = NULL;
                }
                
    accept = gsl_matrix_calloc(filter_n, 4);
    
    fi = 0;

	Rprintf("\nPosterior sampling has started\n");
	for(i=0;i<rep_n;i++){
		v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
		v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
		v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		
        enforce_identifiability(S, ns);
        zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    
		v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
		v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
		v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
		
		if(i%filter_factor==0){	
			gsl_matrix_set(accept,fi,0,v0);
			gsl_matrix_set(accept,fi,1,v1);
			gsl_matrix_set(accept,fi,2,v2);
			gsl_matrix_set(accept,fi,3,v3);

				
			for(j=0;j<ns;j++){
					s1=gsl_matrix_get(S,j,0);
					s2=gsl_matrix_get(S,j,1);
					s3=gsl_matrix_get(S,j,2);
					gsl_matrix_set(collectorS,fi,3*j,s1);
					gsl_matrix_set(collectorS,fi,3*j+1,s2);
					gsl_matrix_set(collectorS,fi,3*j+2,s3);
			}
			
			gsl_matrix_get_col (collectorb, beta, 0);
			gsl_matrix_set_row (collectorbeta, fi, collectorb);
					
			gsl_matrix_get_col (collectorthetaXs, thetaX, 0);
			gsl_matrix_set_row (collectorthetaX, fi, collectorthetaXs);
			
            gsl_matrix_get_col (collectorXs, X, 0);
            gsl_matrix_set_row (collectorX, fi, collectorXs);
            
            if(i%ufilter_factor==0){
                gsl_matrix_get_col (collectorUs, U, 0);
                gsl_matrix_set_row (collectorU, ufi, collectorUs);
                ufi+=1;
            }
            
            gsl_matrix_get_col (collectorLK, lk, 0);
            gsl_matrix_set_row (collectorlk, fi, collectorLK);
            fi+=1;
		}		
		if ((i+1)%ith==0) {
			Rprintf("At iteration %d beta1 %f ",i+1,gsl_matrix_get(beta,1,0));
			Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
			//Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0)+gsl_matrix_get(lk,1,0)+gsl_matrix_get(lk,2,0));
				
		}
	}
}else if(*gear == 0){
	Rprintf("Jumping rules in Proposals\n");
	Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
	Rprintf("beta1 : %f \n",gsl_matrix_get(osu,0,0));
	Rprintf("cov1 : %f \n",gsl_matrix_get(osu,5,0));
	Rprintf("cov2 : %f \n",gsl_matrix_get(osu,6,0));

    fi=0;ufi=0;
	Rprintf("\nBurning has started\n");
	for(i=0;i<rep_b;i++){		
		v0=zHMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		
		enforce_identifiability(S,ns);
        zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
		
		
		v1=zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        v2=zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        v3=zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		if ((i+1)==rep_b) {
			Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
			Rprintf("cov1 %f cov2 %f \n",gsl_matrix_get(beta,2,0),gsl_matrix_get(beta,3,0));
			//Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0));
		}
	}

	Rprintf("Burning has ended\n");

	Rprintf("\nPosterior sampling has started\n");
	for(i=0;i<rep_n;i++){		
		v0=zHMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		
		enforce_identifiability(S,ns);
        zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
		
		
		v1=zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        v2=zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        v3=zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		
		if(i%filter_factor==0){	
			gsl_matrix_set(accept,fi,0,v0);
			gsl_matrix_set(accept,fi,1,v1);
			gsl_matrix_set(accept,fi,2,v2);
			gsl_matrix_set(accept,fi,3,v3);

				
			for(j=0;j<ns;j++){
					s1=gsl_matrix_get(S,j,0);
					s2=gsl_matrix_get(S,j,1);
					s3=gsl_matrix_get(S,j,2);
					gsl_matrix_set(collectorS,fi,3*j,s1);
					gsl_matrix_set(collectorS,fi,3*j+1,s2);
					gsl_matrix_set(collectorS,fi,3*j+2,s3);
			}
			
			gsl_matrix_get_col (collectorb, beta, 0);
			gsl_matrix_set_row (collectorbeta, fi, collectorb);
					
			gsl_matrix_get_col (collectorthetaXs, thetaX, 0);
			gsl_matrix_set_row (collectorthetaX, fi, collectorthetaXs);
			
			gsl_matrix_get_col (collectorLK, lk, 0);
			gsl_matrix_set_row (collectorlk, fi, collectorLK);
			fi+=1;
		}		
		if ((i+1)%ith==0) {
			Rprintf("At iteration %d beta1 %f ",i+1,gsl_matrix_get(beta,1,0));
			Rprintf("cov1 %f cov2 %f \n",gsl_matrix_get(beta,2,0),gsl_matrix_get(beta,3,0));
			//Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0));
				
		}
	}
}

	
	time(&stop);
    diff = difftime(stop, start);
 
    Rprintf("==========Summary Report================\n");
    Rprintf("CPU TIME IS %f\n",diff/60.0);
    
    double mean;
    if(*gear==1){
	    
    k=0;
    for(i=0;i<filter_n;i++){
        for(j=0;j<ns;j++){
            result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
            result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
            result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
        }
    }
    for(i=0;i<filter_n;i++){
        result[k]=gsl_matrix_get(collectorbeta,i,1);k++;
        result[k]=gsl_matrix_get(collectorbeta,i,2);k++;
        result[k]=gsl_matrix_get(collectorbeta,i,3);k++;
        result[k]=gsl_matrix_get(collectorthetaX,i,0);k++;
        result[k]=gsl_matrix_get(collectorthetaX,i,2);k++;
    }
    
    for(j=0;j<ns;j++){
        mean=0.0;
        for(i=0;i<filter_n;i++){
            mean+=gsl_matrix_get(collectorX,i,j);
        }
        mean/=filter_n;
        result[k]=mean;k++;
    }
    
    for(i=0;i<ns;i++){
        howmany=gsl_matrix_int_get(Lookup,i,0);
        for(jj=1;jj<howmany+1;jj++){
            kk=gsl_matrix_int_get(Lookup,i,jj);
            j=gsl_matrix_int_get(ijTable,i,jj);
            if(i>=j) continue;
            mean=0.0;
            for(ii=0;ii<ufilter_n;ii++){
                mean+=gsl_matrix_get(collectorU,ii,kk);
            }
            mean/=ufilter_n;
            uk=ns*i-(i+1)*i*0.5+j-i-1;
            gsl_matrix_set(U2R,0,uk,mean);
            //Rprintf("(i,j,k,kk,uij) %d %d %d %d %f\n",i,j,kk,uk,mean);
        }
    }
    
    for(j=0;j<nu2r;j++){
        v=gsl_matrix_get(U2R,0,j);//Rprintf("(j,uij) %d %f\n",j,v);
        result[k]=v;k++;
    }
    
    for(j=0;j<4;j++) {
        mean=0.0;
        for(i=0;i<filter_n;i++){
            mean+=gsl_matrix_get(accept,i,j);
        }
        mean/=filter_n;
        result[k]=mean;k++;
        if(j==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
        if(j==1) Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
        if(j==2) Rprintf("X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),mean);
        if(j==3) Rprintf("U : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),mean);
    }
    
}else if(*gear==0){
        for(i=0;i<4;i++) {
            mean=0.0;
            for(j=0;j<filter_n;j++){
                mean+=gsl_matrix_get(accept,j,i);
            }
            mean/=filter_n;
            if(i==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
            if(i==1) Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
            if(i==2) Rprintf("cov1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,5,0),mean);
            if(i==3) Rprintf("cov2 : proposal jump %f acceptance rate%f\n",gsl_matrix_get(osu,6,0),mean);
        }
        k=0;
        for(i=0;i<filter_n;i++){
            for(j=0;j<ns;j++){
                result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
            }
        }
        for(i=0;i<filter_n;i++){
            result[k]=gsl_matrix_get(collectorbeta,i,1);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,2);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,3);k++;
        }
    }

    for(i=0;i<filter_n;i++){
	v=gsl_matrix_get(collectorlk,i,0);
	v+=gsl_matrix_get(collectorlk,i,1);
	v+=gsl_matrix_get(collectorlk,i,2);
        result[k]=v;k++;
    }
 

	gsl_matrix_int_free(Lookup);
	gsl_matrix_int_free(ijTable);
	gsl_matrix_free(Jc);
	gsl_matrix_free(X);
	gsl_matrix_free(S);
	gsl_matrix_free(efl);
	gsl_matrix_free(gc);
	gsl_matrix_free(map);
	gsl_matrix_free(Imat);
	gsl_matrix_free(Jmat);
	gsl_matrix_free(iSigma);
	gsl_matrix_free(beta);	
	gsl_matrix_free(thetaX);
	gsl_matrix_free(osu);
	gsl_matrix_free(lk);	
	gsl_matrix_free(Z);
	gsl_matrix_free(Zt);
	gsl_matrix_free(Zgamma);
	gsl_matrix_free(I);
	gsl_matrix_free(U);
	gsl_matrix_free(Metric);
	gsl_matrix_free(logMet);
	gsl_matrix_free(tOne);
	gsl_matrix_free(One);
	gsl_matrix_free(accept);
	gsl_matrix_free(collectorS);
	gsl_matrix_free(collectorbeta);
	gsl_matrix_free(collectorthetaX);
	gsl_matrix_free(collectorlk);
    
    gsl_matrix_free(collectorX);
    gsl_matrix_free(collectorU);
    gsl_matrix_free(U2R);
    
    gsl_vector_free(collectorb);
    gsl_vector_free(collectorthetaXs);
    gsl_vector_free(collectorLK);
    gsl_vector_free(collectorXs);
    gsl_vector_free(collectorUs);


	return;

}



void nz2(
        double *Contact,
        int *n,
        int *repn,
        int *repb,
        double *argv0,
        double *argv1,
        double *argv2,
        int *nHL,
        double *vepsilon,
        int *thinning,
        int *gear,
        double *result)
{
    int ns=*n;
    int nI=0;
    nI2=ns*(ns-1)*0.5;
    int nu2r=(int)ns*(ns-1)*0.5;
    gsl_matrix *U2R = gsl_matrix_calloc(1,nu2r);
    
    int rep_b=*repb;
    int rep_n=*repn;
    
    int i,j,k,howmany,ii,jj,kk,uk;
    double v,Mi,line_t;
    double max_count;
    prior=0.01;
    
    gsl_matrix_int *Lookup=gsl_matrix_int_alloc(ns,ns);
    gsl_matrix_int *ijTable=gsl_matrix_int_alloc(ns,ns);
    
    gsl_matrix *Jc=gsl_matrix_calloc(ns,ns);
    gsl_matrix *X=gsl_matrix_calloc(ns,1);
    gsl_matrix *S=gsl_matrix_calloc(ns,3);
    efl=gsl_matrix_calloc(ns,1);
    gc=gsl_matrix_calloc(ns,1);
    map=gsl_matrix_calloc(ns,1);
    
    gsl_matrix *Imat=gsl_matrix_calloc(ns,ns);
    gsl_matrix *Jmat=gsl_matrix_calloc(ns,ns);
    gsl_matrix *iSigma=gsl_matrix_calloc(ns,ns);
    gsl_matrix *beta=gsl_matrix_calloc(5,1);
    gsl_matrix *thetaX=gsl_matrix_calloc(3,1);
    gsl_matrix *osu=gsl_matrix_calloc(7,1);
    gsl_matrix *lk=gsl_matrix_calloc(3,1);
    
    gsl_matrix *Z=gsl_matrix_calloc(ns,2);
    gsl_matrix *Zt=gsl_matrix_calloc(2,ns);
    gsl_matrix *Zgamma=gsl_matrix_calloc(ns,ns);
    
    IH=gsl_matrix_calloc(ns,ns);
    gsl_matrix_set_identity (IH);
    
    
    gsl_vector *collectorXs=gsl_vector_calloc(ns);
    gsl_vector *collectorb=gsl_vector_calloc(5);
    gsl_vector *collectorthetaXs=gsl_vector_calloc(3);
    gsl_vector *collectorLK=gsl_vector_calloc(3);
    
    int  filter_factor=*thinning;
    int  filter_n=rep_n/filter_factor;
    int ufilter_factor=filter_factor*10;
    int ufilter_n=rep_n/ufilter_factor;
    
    gsl_matrix *accept=gsl_matrix_calloc(filter_n,4);
    gsl_matrix *collectorS=gsl_matrix_calloc(filter_n,3*ns);
    gsl_matrix *collectorX=gsl_matrix_calloc(filter_n,ns);

    gsl_matrix *collectorbeta=gsl_matrix_calloc(filter_n,5);
    gsl_matrix *collectorthetaX=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorlk=gsl_matrix_calloc(filter_n,3);
    
    
    
    gsl_matrix_set(osu,0,0,*argv0); //betaOne
    gsl_matrix_set(osu,1,0,*argv1); //X
    gsl_matrix_set(osu,2,0,*argv2); //U
    
    gsl_matrix_set(osu,3,0,0.07); //S
    gsl_matrix_set(osu,4,0,0.03); //betaZero
    
    gsl_matrix_set(osu,5,0,0.03); //EFL
    gsl_matrix_set(osu,6,0,0.03); //GC
    
    if(*gear == 0){
        gsl_matrix_set(osu,5,0,*argv1);
        gsl_matrix_set(osu,6,0,*argv2);
    }
    
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, 3);
    
    double initial=gsl_rng_uniform (r);
    
    ubound=1.0;
    lbound=-1.0*ubound;
    
    line_t=ubound*initial;
    
    HL=*nHL;
	epsilon=*vepsilon;
	double log_eps = log(epsilon);
    
    k=0;
    for(i=0;i<ns;i++) {
        for(j=i+1;j<ns;j++){
            v=Contact[k];
            if(v != 0 ) nI++;
            gsl_matrix_set(Jc,i,j,v);
            gsl_matrix_set(Jc,j,i,v);
            k++;
        }
    }
    
    gsl_vector *collectorUs=gsl_vector_calloc(nI);
    gsl_matrix *collectorU=gsl_matrix_calloc(ufilter_n,nI);
    
    
    for(i=0;i<ns;i++){
        gsl_matrix_int_set(Lookup,i,0,0);
    }
    //printf("Number of non-zero counts %d \n",nI);
    
    gsl_matrix *I=gsl_matrix_calloc(nI,1);
    gsl_matrix *U=gsl_matrix_calloc(nI,1);
    
    gsl_matrix *Metric=gsl_matrix_calloc(nI,1);
    gsl_matrix *logMet=gsl_matrix_calloc(nI,1);
    
    ZZ=gsl_matrix_calloc(2,2);
    gamma_hat=gsl_matrix_calloc(2,ns);
    iZ=gsl_matrix_calloc(2,2);
    
    k=0;
    for(i=0;i<ns;i++){
        gsl_matrix_set(efl,i,0,0.0);
        gsl_matrix_set(gc,i,0,0.0);
        gsl_matrix_set(map,i,0,0.0);
    }
    
    k=0;max_count=0.0;Mi=0.0;
    for(i=0;i<ns;i++){
        for(j=i+1;j<ns;j++){
            v=gsl_matrix_get(Jc,i,j);
            if(max_count<v) max_count=v;
            if(v != 0.0){
                Mi+=v;
                howmany=gsl_matrix_int_get(Lookup,i,0)+1;
                gsl_matrix_int_set(Lookup,i,howmany,k);
                gsl_matrix_int_set(ijTable,i,howmany,j);
                gsl_matrix_int_set(Lookup,i,0,howmany);
                gsl_matrix_int_set(ijTable,i,0,howmany);
                //printf("%d %d %d %d \n",i,j,howmany,k);
                howmany=gsl_matrix_int_get(Lookup,j,0)+1;
                gsl_matrix_int_set(Lookup,j,howmany,k);
                gsl_matrix_int_set(ijTable,j,howmany,i);
                gsl_matrix_int_set(Lookup,j,0,howmany);
                gsl_matrix_int_set(ijTable,j,0,howmany);
                
                gsl_matrix_set(I,k,0,v);
                //printf("%d %d %d %d %f\n",i,j,howmany,k,v);
                k++;
            }
        }
        gsl_matrix_set(S,i,0,line_t*i/(ns-1));
        gsl_matrix_set(S,i,1,0.0);
        gsl_matrix_set(S,i,2,0.0);
        //Rprintf("Initial beta0 beta1 %f %f X S %f\n",initial*log(max_count),-2.0*initial,line_t);
    }
    
    enforce_identifiability(S, ns);
    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    
    gsl_matrix_set(beta,0,0,initial*log(max_count));
    //gsl_matrix_set(beta,1,0,-2.0*initial);
    
    //Rprintf("Initial beta0 beta1 %f %f X S %f\n",initial*log(max_count),-2.0*initial,line_t);
    
    gsl_matrix_set(thetaX,0,0,1.0);
    gsl_matrix_set(thetaX,1,0,0.0);
    gsl_matrix_set(thetaX,2,0,0.1);
    
    gsl_matrix_set_identity(iSigma);
    gsl_matrix_set_identity(Imat);
    gsl_matrix_set_all(Jmat,1.0);
    
    gsl_matrix *tOne=gsl_matrix_calloc(1,ns);
    gsl_matrix *One=gsl_matrix_calloc(ns,1);
    gsl_matrix_set_all(tOne,1.0);
    gsl_matrix_set_all(One,1.0);
    
    gsl_matrix_set_all(Jmat,1.0);
    
    gsl_matrix_set(lk,0,0,0.0);
    gsl_matrix_set(lk,1,0,0.0);
    gsl_matrix_set(lk,2,0,0.0);
    
    double sv=0.0;
    double s1,s2,s3;
    double diff = 0.0;
    time_t start;
    time_t stop;
    time(&start);
    
    
    double v0,v1,v2,v3,mv;
    double scale,b1;
    int fi=0;int ufi=0;
    double adjust=0.1;
    
    int nn=0;
    int howoften=0;
    
/*    v=gsl_matrix_get(beta,0,0);
    for(i=0;i<ns;i++){
        line_t=gsl_ran_gaussian (r, 0.1*ubound);
        gsl_matrix_set(X,i,0,line_t+0.5*v);
    }
*/
    //===========================================
    

    //printf("MAP is set off\n");
    //printf("%f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
    //	printf("INITIALIZATION \n");
    //=====================================
    gsl_matrix_set(beta,1,0,-2.0*initial);
    
    sv= zgetlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    for(i=0;i<10000;i++){
        zupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
        zupdatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
    }

    Rprintf("Initial1 : beta1 %f %f ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0));
    Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    v=gsl_matrix_get(beta,0,0);
    v*=0.5;
    gsl_matrix_set(beta,0,0,1.5);
    b1=gsl_matrix_get(beta,1,0);
    scale=exp(2.0*(v-1.5)/b1);
    
    for(i=0;i<ns;i++){
        s1=gsl_matrix_get(S,i,0);s2=gsl_matrix_get(S,i,1);s3=gsl_matrix_get(S,i,2);
        s1*=scale;s2*=scale;s3*=scale;
        gsl_matrix_set(S,i,0,s1);gsl_matrix_set(S,i,1,s2);gsl_matrix_set(S,i,2,s3);
        gsl_matrix_set(X,i,0,1.5);
    }
    
    enforce_identifiability(S, ns);
    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    
    sv= zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    sv= zgetlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,sv);
    sv= zgetlkU(U,thetaX,nI);
    gsl_matrix_set(lk,2,0,sv);
    gsl_matrix_set(beta,4,0,1.0);
    
    Rprintf("Initial2 : beta1 %f %f ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0));
    Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    

    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
    int iter1 = 5000;
    int iter2 = 2000;
    int iter3 = 5000;
    int iter4 = 2000;
    int iter5 = 10000;
    
    int block1_size = 500;
    int block2_size = 500;
    int block3_size = 500;
    int block4_size = 500;
    int block5_size = 500;
    
    double kappa_init = 0.5;
    
    double mu;        // log of the initial epsilon * 10 (Stan uses log(10*epsilon_init))
    double log_eps_avg; // running average of log eps (smoothed)
    double Hbar;      // running average of (target - accept_stat)
    double gamma;     // controls shrinkage (Stan uses 0.05)
    double t0;        // stabilizes early iterations (Stan uses 10)
    int adapt_iter;   // counter for adaptation iterations
    
    gamma = 0.05;
    t0 = 10.0;
    
    // initialize dual-averaging
    mu = log(10 * epsilon);   // Stan uses log(10 * epsilon_init)
    log_eps_avg = log_eps;
    Hbar = 0.0;
    adapt_iter = 1;
    
    int adapt_t_b1 = 1;
    int adapt_t_X  = 1;
    int adapt_t_U  = 1;
    
    Rprintf("\n====================================\n");
    int check=0;
    if(*gear == 1){
        Rprintf("Jumping rules in Proposals\n");
        Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
        Rprintf("beta1 : osu %f \n",gsl_matrix_get(osu,0,0));
        Rprintf("X : osu %f \n",gsl_matrix_get(osu,1,0));
        Rprintf("U : osu %f \n",gsl_matrix_get(osu,2,0));
        
        Rprintf("Looking for jumping rules...\n");

        int max_cycle = 10;
        
        for (int cycle = 0; cycle < max_cycle; cycle++) {
            
            // Parameters for adaptation
            double target = 0.95;    // target acceptance rate
            double kappa = kappa_init;    // learning rate
            double kappa_decay = 0.95; //learning rate decay
            
            // reset log-epsilon?
            
            // log_eps = log(epsilon_init);
            // epsilon  = epsilon_init;
            
            Rprintf("\n==================== Adaptation Cycle %d ====================\n", cycle+1);
        
            /* -------------------- Phase 1: Step-size adaptation -------------------- */
            Rprintf("Phase 1: Step-size adaptation\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            int n_iter = iter1; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            int block_size = block1_size;
            double acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            double acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;

            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                // every 'block_size' iterations, update epsilon and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    
                    dual_averaging_step_size_update(acc_S, target,
                      &log_eps, &log_eps_avg,
                      &Hbar, mu,
                      gamma, t0, kappa,
                      adapt_iter);

                    epsilon = exp(log_eps);
                    adapt_iter++;
                    
                    adapt_step_num(&HL, acc_S);
                    
                    adapt_jump_proposals(osu, acc_b1, acc_X, acc_U, 0.25, adapt_t_b1, adapt_t_X, adapt_t_U);
                    
                    adapt_t_b1 += 1;
                    adapt_t_X += 1;
                    adapt_t_U += 1;
                    
                    
                    // Report diagnostics
                    Rprintf("\n[Phase 1 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
            }
            
            /* ---------- Phase summary ---------- */
            double mean_acc_S  = acc_phase_S  / n_iter;
            double mean_acc_b1 = acc_phase_b1 / n_iter;
            double mean_acc_X  = acc_phase_X  / n_iter;
            double mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 1 completed | Cycle %d ===\n", cycle + 1);
            
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
        
        
            /* -------------------- Phase 2: Step-size + Mass adaptation -------------------- */
            Rprintf("Phase 2: Mass + Step-size adaptation\n");
            
            kappa = kappa_init;
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter2; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block2_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
        
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                

                
                // every 'block_size' iterations, update epsilon and mass matrix and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
            
                    // Report diagnostics
                    Rprintf("\n[Phase 2 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
            }
            
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 2 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 3: Freeze mass, stabilize epsilon -------------------- */
            Rprintf("Phase 3: Freeze Mass, tune step size \n");
            
            kappa = kappa_init;
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter3; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block3_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                // every 'block_size' iterations, update epsilon and report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    dual_averaging_step_size_update(acc_S, target,
                      &log_eps, &log_eps_avg,
                      &Hbar, mu,
                      gamma, t0, kappa,
                      adapt_iter);

                    epsilon = exp(log_eps);
                    adapt_iter++;
                    
                    
            
                    // Report diagnostics
                    Rprintf("\n[Phase 3 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 3 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 4: Stability run (no adaptation) -------------------- */
            Rprintf("Phase 4: Stability run\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter4; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block4_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                
                // every 'block_size' iterations report acceptance
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
                    
            
                    // Report diagnostics
                    Rprintf("\n[Phase 4 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 4 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
            
        
            /* -------------------- Phase 5: Tune jump proposals -------------------- */
            Rprintf("Phase 5: Tuning beta1, X, and U proposals\n");
            
            if (accept != NULL)
            gsl_matrix_free(accept);
        
            n_iter = iter5; 
            accept = gsl_matrix_calloc(n_iter, 4);
            
            block_size = block5_size;
            acc_block_S = 0.0, acc_block_b1 = 0.0, acc_block_X = 0.0, acc_block_U = 0.0;
            acc_phase_S = 0.0, acc_phase_b1 = 0.0, acc_phase_X = 0.0, acc_phase_U = 0.0;
            
            fi=0;ufi=0;howoften=1;
            
            for (i = 0; i < n_iter; i++) {
                v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
                // v0=zHMCupdateS2_mass(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
                enforce_identifiability(S, ns);
                zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
                v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
            
                
                // accumulate for block and overall means
                
                acc_block_S  += v0;  acc_phase_S  += v0;
                acc_block_b1 += v1;  acc_phase_b1 += v1;
                acc_block_X  += v2;  acc_phase_X  += v2;
                acc_block_U  += v3;  acc_phase_U  += v3;
                
                
                // every 'block_size' iterations report acceptance and adapt jump proposals
                if ((i + 1) % block_size == 0) {
                    double acc_S  = acc_block_S  / block_size;
                    double acc_b1 = acc_block_b1 / block_size;
                    double acc_X  = acc_block_X  / block_size;
                    double acc_U  = acc_block_U  / block_size;
            
                    // === Adapt jump proposals ===
                    
                    adapt_jump_proposals(osu, acc_b1, acc_X, acc_U, 0.25, adapt_t_b1, adapt_t_X, adapt_t_U);
                    
                    adapt_t_b1 += 1;
                    adapt_t_X += 1;
                    adapt_t_U += 1;
                    
                    
                    // Report diagnostics
                    Rprintf("\n[Phase 5 | Cycle %d | Iter %5d]\n", cycle + 1, i + 1);
                    Rprintf("  Adjusting : beta1 %.6f | sigma2_x %.6f | sigma2_u %.6f\n",
                            gsl_matrix_get(beta, 1, 0),
                            gsl_matrix_get(thetaX, 0, 0),
                            gsl_matrix_get(thetaX, 2, 0));
                    Rprintf("  log-likelihoods: %.6f %.6f %.6f\n",
                            gsl_matrix_get(lk, 0, 0),
                            gsl_matrix_get(lk, 1, 0),
                            gsl_matrix_get(lk, 2, 0));
            
                    Rprintf("  Block acceptance:\n");
                    Rprintf("    S     : acc=%.3f | ε=%.6f | HL=%d\n", acc_S, epsilon, HL);
                    Rprintf("    beta1 : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,0,0), acc_b1);
                    Rprintf("    X     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,1,0), acc_X);
                    Rprintf("    U     : jump=%.5f | acc=%.3f\n", gsl_matrix_get(osu,2,0), acc_U);
            
                    // Reset block accumulators
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                }
                
            }
            
            
            /* ---------- Phase summary ---------- */
            mean_acc_S  = acc_phase_S  / n_iter;
            mean_acc_b1 = acc_phase_b1 / n_iter;
            mean_acc_X  = acc_phase_X  / n_iter;
            mean_acc_U  = acc_phase_U  / n_iter;
            
            Rprintf("\n=== Phase 5 completed | Cycle %d ===\n", cycle + 1);
            
            
            Rprintf("  Mean acceptance: S=%.3f, beta1=%.3f, X=%.3f, U=%.3f\n",
                    mean_acc_S, mean_acc_b1, mean_acc_X, mean_acc_U);
            Rprintf("  Final epsilon=%.6f | log_eps=%.3f | HL=%d | kappa=%.4f\n ",
                    epsilon, log_eps, HL, kappa);
                    
            Rprintf("  Final params: beta1=%.6f, sigma2_x=%.6f, sigma2_u=%.6f\n",
                    gsl_matrix_get(beta, 1, 0),
                    gsl_matrix_get(thetaX, 0, 0),
                    gsl_matrix_get(thetaX, 2, 0));
            Rprintf("  Final log-likelihoods: %.6f %.6f %.6f\n",
                    gsl_matrix_get(lk, 0, 0),
                    gsl_matrix_get(lk, 1, 0),
                    gsl_matrix_get(lk, 2, 0));
                    
                    
            /* -------------------- Early stopping check -------------------- */
            int good_S  = (mean_acc_S  > 0.80 && mean_acc_S  < 0.85);
            int good_b1 = (mean_acc_b1 > 0.20 && mean_acc_b1 < 0.35);
            int good_X  = (mean_acc_X  > 0.20 && mean_acc_X  < 0.35);
            int good_U  = (mean_acc_U  > 0.20 && mean_acc_U  < 0.35);
            
            if (good_S && good_b1 && good_X && good_U) {
                Rprintf("\n*** Early stopping: acceptance rates are good at Cycle %d ***\n",
                        cycle + 1);
            
                // free matrices before break
                if (accept != NULL) {
                    gsl_matrix_free(accept);
                    accept = NULL;
                }
            
                // break out of cycle loop and skip resetting structure
                break;
            }

        
            if(cycle < max_cycle-1){
                /* -------------------- Reset between cycles -------------------- */
                Rprintf("Reset structure for next cycle\n");
                
                    // ---- Phase accumulators and counters ----
                    acc_block_S = acc_block_b1 = acc_block_X = acc_block_U = 0.0;
                    acc_phase_S = acc_phase_b1 = acc_phase_X = acc_phase_U = 0.0;
                    fi = 0;
                    ufi = 0;
                    
                    Hbar        = 0.0;           // running average of (target - accept)
                    log_eps_avg = 0.0;           // smoothed average of log-epsilon
                    adapt_iter  = 1;             // dual-averaging iteration counter
                    mu = log(10.0 * *vepsilon);
                    adapt_t_b1 = 1;
                    adapt_t_X  = 1;
                    adapt_t_U  = 1;
                    
                    // ---- Free any leftover matrices from previous phases ----
                    if (accept != NULL) {
                        gsl_matrix_free(accept);
                        accept = NULL;
                    }
                
                    gsl_matrix_set(beta,0,0,initial*log(max_count));
                    v=gsl_matrix_get(beta,0,0);
                    for(i=0;i<ns;i++){
                        gsl_matrix_set(S,i,0,line_t*i/(ns-1));
                        gsl_matrix_set(S,i,1,0.0);
                        gsl_matrix_set(S,i,2,0.0);
                        //line_t=gsl_ran_gaussian (r, 0.1*ubound);
    
                    }
                    
                    gsl_matrix_set(thetaX,0,0,0.1); 
                    gsl_matrix_set(thetaX,1,0,0.0);
                    gsl_matrix_set(thetaX,2,0,0.1);
                    
                    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                    gsl_matrix_set(beta,1,0,-2.0*initial);
                    sv= zgetlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
                    gsl_matrix_set(lk,0,0,sv);
                    for(i=0;i<10000;i++){
                        zupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                        zupdatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                        zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                    }
                    
                    
                    v=gsl_matrix_get(beta,0,0);
                    v*=0.5;
                    gsl_matrix_set(beta,0,0,1.5);
                    b1=gsl_matrix_get(beta,1,0);
                    scale=exp(2.0*(v-1.5)/b1);
                    
                    for(i=0;i<ns;i++){
                        s1=gsl_matrix_get(S,i,0);s2=gsl_matrix_get(S,i,1);s3=gsl_matrix_get(S,i,2);
                        s1*=scale;s2*=scale;s3*=scale;
                        gsl_matrix_set(S,i,0,s1);gsl_matrix_set(S,i,1,s2);gsl_matrix_set(S,i,2,s3);
                        gsl_matrix_set(X,i,0,1.5);
                    }
                    
                    
                    enforce_identifiability(S, ns);
                    zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                    
                    sv= zgetlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
                    gsl_matrix_set(lk,0,0,sv);
                    sv= zgetlkX(X,thetaX,ns,iSigma,beta);
                    gsl_matrix_set(lk,1,0,sv);
                    sv= zgetlkU(U,thetaX,nI);
                    gsl_matrix_set(lk,2,0,sv);
                    gsl_matrix_set(beta,4,0,1.0);
            
            }
        }

        
        
        fi=0;ufi=0;
        
        // set accept there
        
        accept = gsl_matrix_calloc(rep_b, 4);
        
        double acc_rate = 0.0;
        Rprintf("\nBurning has started\n");
        for(i=0;i<rep_b;i++){
            v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
            v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
            v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
            enforce_identifiability(S, ns);
            zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
        
            v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
            v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
//            v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
            
            gsl_matrix_set(accept,fi,0,v0);
            gsl_matrix_set(accept,fi,1,v1);
            gsl_matrix_set(accept,fi,2,v2);
            gsl_matrix_set(accept,fi,3,v3);
            
            fi += 1;    

            if ((i+1)==rep_b) {
                Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
                Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
                Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
            }
        }
        
        Rprintf("Burning has ended\n");
        
        for(j=0;j<4;j++) {
            acc_rate=0.0;
            for(i=0;i<rep_b;i++){
                acc_rate+=gsl_matrix_get(accept,i,j);
            }
            acc_rate/=rep_b;
            result[k]=acc_rate;k++;
            if(j==0) Rprintf("Burn-in S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL, epsilon, acc_rate);
            if(j==1) Rprintf("Burn-in beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),acc_rate);
            if(j==2) Rprintf("Burn-in X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),acc_rate);
            if(j==3) Rprintf("Burn-in U : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),acc_rate);
        }
        
        Rprintf("\nPosterior sampling has started\n");
        
        if (accept != NULL) {
                    gsl_matrix_free(accept);
                    accept = NULL;
                }
                
        accept = gsl_matrix_calloc(filter_n, 4);
        
        fi = 0;
        
        for(i=0;i<rep_n;i++){
            v2=zupdateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
            v3=zupdateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable,S);
            v0=zHMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
            enforce_identifiability(S, ns);
            zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
        
            v1=zupdatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            v=zupdateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
            v=zupdateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
//            v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
            
            if(i%filter_factor==0){
                gsl_matrix_set(accept,fi,0,v0);
                gsl_matrix_set(accept,fi,1,v1);
                gsl_matrix_set(accept,fi,2,v2);
                gsl_matrix_set(accept,fi,3,v3);
                
                
                for(j=0;j<ns;j++){
                    s1=gsl_matrix_get(S,j,0);
                    s2=gsl_matrix_get(S,j,1);
                    s3=gsl_matrix_get(S,j,2);
                    gsl_matrix_set(collectorS,fi,3*j,s1);
                    gsl_matrix_set(collectorS,fi,3*j+1,s2);
                    gsl_matrix_set(collectorS,fi,3*j+2,s3);
                }
                
                gsl_matrix_get_col (collectorb, beta, 0);
                gsl_matrix_set_row (collectorbeta, fi, collectorb);
                
                gsl_matrix_get_col (collectorthetaXs, thetaX, 0);
                gsl_matrix_set_row (collectorthetaX, fi, collectorthetaXs);
                
                gsl_matrix_get_col (collectorXs, X, 0);
                gsl_matrix_set_row (collectorX, fi, collectorXs);
                
                if(i%ufilter_factor==0){
                    gsl_matrix_get_col (collectorUs, U, 0);
                    gsl_matrix_set_row (collectorU, ufi, collectorUs);
                    ufi+=1;
                }
                
                gsl_matrix_get_col (collectorLK, lk, 0);
                gsl_matrix_set_row (collectorlk, fi, collectorLK);
                fi+=1;
            }
            if ((i+1)%ith==0) {
                Rprintf("At iteration %d beta1 %f ",i+1,gsl_matrix_get(beta,1,0));
                Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
                //Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0)+gsl_matrix_get(lk,1,0)+gsl_matrix_get(lk,2,0));
                
            }
        }
    }else if(*gear == 0){
        Rprintf("Jumping rules in Proposals\n");
        Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
        Rprintf("beta1 : %f \n",gsl_matrix_get(osu,0,0));
//        Rprintf("cov1 : %f \n",gsl_matrix_get(osu,5,0));
//        Rprintf("cov2 : %f \n",gsl_matrix_get(osu,6,0));
        
        fi=0;
        Rprintf("\nBurning has started\n");
        for(i=0;i<rep_b;i++){
            v0=zHMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
            enforce_identifiability(S, ns);
            zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
            v1=zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
//            v2=zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
//            v3=zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            if ((i+1)==rep_b) {
                Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
                Rprintf("cov1 %f cov2 %f \n",gsl_matrix_get(beta,2,0),gsl_matrix_get(beta,3,0));
                //Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0));
            }
        }
        
        Rprintf("Burning has ended\n");
        
        Rprintf("\nPosterior sampling has started\n");
        for(i=0;i<rep_n;i++){
            v0=zHMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            
            enforce_identifiability(S, ns);
            zgetDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
            
            v1=zupdatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
//            v2=zupdatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
//            v3=zupdatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            
            if(i%filter_factor==0){
                gsl_matrix_set(accept,fi,0,v0);
                gsl_matrix_set(accept,fi,1,v1);
                gsl_matrix_set(accept,fi,2,v2);
                gsl_matrix_set(accept,fi,3,v3);
                
                
                for(j=0;j<ns;j++){
                    s1=gsl_matrix_get(S,j,0);
                    s2=gsl_matrix_get(S,j,1);
                    s3=gsl_matrix_get(S,j,2);
                    gsl_matrix_set(collectorS,fi,3*j,s1);
                    gsl_matrix_set(collectorS,fi,3*j+1,s2);
                    gsl_matrix_set(collectorS,fi,3*j+2,s3);
                }
                
                gsl_matrix_get_col (collectorb, beta, 0);
                gsl_matrix_set_row (collectorbeta, fi, collectorb);
                
                gsl_matrix_get_col (collectorthetaXs, thetaX, 0);
                gsl_matrix_set_row (collectorthetaX, fi, collectorthetaXs);
                
                gsl_matrix_get_col (collectorLK, lk, 0);
                gsl_matrix_set_row (collectorlk, fi, collectorLK);
                fi+=1;
            }
            if ((i+1)%ith==0) {
                Rprintf("At iteration %d beta1 %f ",i+1,gsl_matrix_get(beta,1,0));
                //Rprintf("cov1 %f cov2 %f \n",gsl_matrix_get(beta,2,0),gsl_matrix_get(beta,3,0));
                Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0));
                
            }
            
        }
        
        
    }
    
    
    time(&stop);
    diff = difftime(stop, start);
    
    Rprintf("==========Summary Report================\n");
    Rprintf("CPU TIME IS %f\n",diff/60.0);
    
    double mean;
    if(*gear==1){
        
        k=0;
        for(i=0;i<filter_n;i++){
            for(j=0;j<ns;j++){
                result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
            }
        }
        for(i=0;i<filter_n;i++){
            result[k]=gsl_matrix_get(collectorbeta,i,1);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,2);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,3);k++;
            result[k]=gsl_matrix_get(collectorthetaX,i,0);k++;
            result[k]=gsl_matrix_get(collectorthetaX,i,2);k++;
        }
        
        for(j=0;j<ns;j++){
            mean=0.0;
            for(i=0;i<filter_n;i++){
                mean+=gsl_matrix_get(collectorX,i,j);
            }
            mean/=filter_n;
            result[k]=mean;k++;
        }
        
        for(i=0;i<ns;i++){
            howmany=gsl_matrix_int_get(Lookup,i,0);
            for(jj=1;jj<howmany+1;jj++){
                kk=gsl_matrix_int_get(Lookup,i,jj);
                j=gsl_matrix_int_get(ijTable,i,jj);
                if(i>=j) continue;
                mean=0.0;
                for(ii=0;ii<ufilter_n;ii++){
                    mean+=gsl_matrix_get(collectorU,ii,kk);
                }
                mean/=ufilter_n;
                uk=ns*i-(i+1)*i*0.5+j-i-1;
                gsl_matrix_set(U2R,0,uk,mean);
                //Rprintf("(i,j,k,kk,uij) %d %d %d %d %f\n",i,j,kk,uk,mean);
            }
        }
                
        for(j=0;j<nu2r;j++){
            v=gsl_matrix_get(U2R,0,j);//Rprintf("(j,uij) %d %f\n",j,v);
            result[k]=v;k++;
        }
        
        for(j=0;j<4;j++) {
            mean=0.0;
            for(i=0;i<filter_n;i++){
                mean+=gsl_matrix_get(accept,i,j);
            }
            mean/=filter_n;
            result[k]=mean;k++;
            if(j==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
            if(j==1) Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
            if(j==2) Rprintf("X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),mean);
            if(j==3) Rprintf("U : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),mean);
        }
        
    }else if(*gear==0){
        for(i=0;i<4;i++) {
            mean=0.0;
            for(j=0;j<filter_n;j++){
                mean+=gsl_matrix_get(accept,j,i);
            }
            mean/=filter_n;
            if(i==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
            if(i==1) Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
            //if(i==2) Rprintf("cov1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,5,0),mean);
            //if(i==3) Rprintf("cov2 : proposal jump %f acceptance rate%f\n",gsl_matrix_get(osu,6,0),mean);
        }
        k=0;
        for(i=0;i<filter_n;i++){
            for(j=0;j<ns;j++){
                result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
            }
        }
        for(i=0;i<filter_n;i++){
            result[k]=gsl_matrix_get(collectorbeta,i,1);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,2);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,3);k++;
        }
    }
    
    for(i=0;i<filter_n;i++){
	v=gsl_matrix_get(collectorlk,i,0);
	v+=gsl_matrix_get(collectorlk,i,1);
	v+=gsl_matrix_get(collectorlk,i,2);
        result[k]=v;k++;
    }

    gsl_matrix_int_free(Lookup);
    gsl_matrix_int_free(ijTable);
    gsl_matrix_free(Jc);
    gsl_matrix_free(X);
    gsl_matrix_free(S);
    gsl_matrix_free(efl);
    gsl_matrix_free(gc);
    gsl_matrix_free(map);
    gsl_matrix_free(Imat);	
    gsl_matrix_free(Jmat);	
    gsl_matrix_free(iSigma);
    gsl_matrix_free(beta);	
    gsl_matrix_free(thetaX);
    gsl_matrix_free(osu);
    gsl_matrix_free(lk);	
    gsl_matrix_free(Z);
    gsl_matrix_free(Zt);
    gsl_matrix_free(Zgamma);
    gsl_matrix_free(IH);
    gsl_matrix_free(I);
    gsl_matrix_free(U);
    gsl_matrix_free(Metric);
    gsl_matrix_free(logMet);
    gsl_matrix_free(tOne);
    gsl_matrix_free(One);
    gsl_matrix_free(accept);
    gsl_matrix_free(collectorS);
    gsl_matrix_free(collectorbeta);
    gsl_matrix_free(collectorthetaX);
    gsl_matrix_free(collectorlk);
    gsl_matrix_free(collectorX);
    gsl_matrix_free(collectorU);
    gsl_matrix_free(U2R);
    
    gsl_vector_free(collectorb);
    gsl_vector_free(collectorthetaXs);
    gsl_vector_free(collectorLK);
    gsl_vector_free(collectorXs);
    gsl_vector_free(collectorUs);
    
    return;
    
}

