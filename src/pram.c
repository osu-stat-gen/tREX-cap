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



void getDist(S,Metric,logMet,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double x,y,z,xj,yj,zj,d,v;
    int i,j,k,howmany,kk;
    kk=0;
    for(i=0;i<ns;i++){
        x=gsl_matrix_get(S,i,0);
        y=gsl_matrix_get(S,i,1);
        z=gsl_matrix_get(S,i,2);
        howmany=gsl_matrix_int_get(Lookup,i,0);
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(S,j,0);
            yj=gsl_matrix_get(S,j,1);
            zj=gsl_matrix_get(S,j,2);
            d=(x-xj)*(x-xj)+(y-yj)*(y-yj)+(z-zj)*(z-zj);
            v=sqrt(d);
            gsl_matrix_set(Metric,kk,0,v);
            gsl_matrix_set(logMet,kk,0,log(v));
        }
    }
    return;
}

double getlkI(I,logMet,X,U,beta,ns,Lookup,ijTable)
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
    double retv=0.0;
    double lmet,iv,part1;
    double efi,efj,gci,gcj,mapi,mapj;
    
    
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
            retv+=iv*part1-exp(part1);
            
        }
    }
    return retv;
}


double getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta;
int ns;
gsl_matrix_int *Lookup,*ijTable;
{
    int i,j,k,kk,howmany;
    double b1=gsl_matrix_get(beta,1,0);
    double retv=0.0;
    double xi,xj,lmet,iv,part1,uij;
    
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
            part1=b1*lmet+xi+xj+uij;
            retv+=iv*part1-exp(part1);
            
        }
    }
    return retv;
}

double getlkX(X,thetaX,ns,iSigma,beta)
gsl_matrix *X,*thetaX,*iSigma,*beta;
int ns;
{
    double retv,v;
    double sigma=gsl_matrix_get(thetaX,0,0);
    double beta0=gsl_matrix_get(beta,0,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    double efi,gci,mapi;
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


double getlkU(U,thetaX,nI)
gsl_matrix *U,*thetaX;
int nI;
{
    double retv;
    double vu;
    
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

double updateX(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S)
gsl_matrix *X,*U,*I,*logMet,*beta,*lk,*osu,*thetaX,*iSigma,*S;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    double delta,deltaI,deltaX,accept;
    double lkIv,lk0,xni,xi,xj,lmet,iv,v,uij;
    double npart1, part1;
    
    int i,j,k,kk,howmany;
    double accepted=0.0;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    double sigma=gsl_matrix_get(thetaX,0,0);
    double osuv=gsl_matrix_get(osu,1,0);
    double rho=gsl_matrix_get(thetaX,1,0);
    
    double tOneX=0.0;
    double ntOneX=0.0;
    double tXX=0.0;
    double deltatXX=0.0;
    double efi,efj,gci,gcj,mapi,mapj,bias;
    
    for(i=0;i<ns;i++){
        v=gsl_matrix_get(X,i,0)-b0;
        tOneX+=v;
        tXX+=v*v;
    }
    double mx=0.0;
    
    for(i=0;i<ns;i++){
        ntOneX=tOneX;
        xi=gsl_matrix_get(X,i,0);
        xni=gsl_ran_gaussian (r, osuv);
        xni+=xi;
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        
        howmany=gsl_matrix_int_get(Lookup,i,0);
        lkIv=0.0;lk0=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            bias=b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            npart1=b1*lmet+xni+xj+uij+bias;
            
            part1=b1*lmet+xi+xj+uij+bias;
            lkIv+=iv*npart1-exp(npart1);
            lk0+=iv*part1-exp(part1);
        }
        deltaI=lkIv-lk0;
        
        ntOneX=tOneX-xi+xni;
        deltatXX=(xni-b0)*(xni-b0)-(xi-b0)*(xi-b0);
        deltaX=1.0/(1.0-rho)*deltatXX;
        
        deltaX-=rho/((1.0-rho)*(1.0-rho+ns*rho))*(ntOneX*ntOneX-tOneX*tOneX);
        deltaX*=-0.5/sigma;
        
        delta=deltaI+deltaX;
        if(delta>=0.0){
            gsl_matrix_set(X,i,0,xni);
            v=gsl_matrix_get(lk,0,0);
            v+=deltaI;gsl_matrix_set(lk,0,0,v);
            v=gsl_matrix_get(lk,1,0);
            v+=deltaX;gsl_matrix_set(lk,1,0,v);
            tOneX+=-xi+xni;
            tXX+=deltatXX;
            accepted+=1.0;
            mx+=xni;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(X,i,0,xni);
                v=gsl_matrix_get(lk,0,0);
                v+=deltaI;gsl_matrix_set(lk,0,0,v);
                v=gsl_matrix_get(lk,1,0);
                v+=deltaX;gsl_matrix_set(lk,1,0,v);
                tOneX+=-xi+xni;
                tXX+=deltatXX;
                accepted+=1.0;
                mx+=xni;
            }else mx+=xi;
        }
    }
    
    return accepted/ns;
}



double updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric)
gsl_matrix *X,*U,*I,*logMet,*beta,*lk,*osu,*thetaX,*iSigma,*S,*Metric;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    double delta,deltaI,deltaX,accept;
    double lkIv,lk0,xni,xi,xj,lmet,iv,v,uij;
    double npart1, part1;
    
    int i,j,k,kk,howmany;
    double accepted=0.0;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    double b3=gsl_matrix_get(beta,3,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    double sigma=gsl_matrix_get(thetaX,0,0);
    double osuv=gsl_matrix_get(osu,1,0);
    double efi,gci,mapi,bias;
    
    
    double mx=0.0;
    double check=0.0;
    
    for(i=0;i<ns;i++){
        xi=gsl_matrix_get(X,i,0);
        xni=gsl_ran_gaussian (r, osuv);
        xni+=xi;
        howmany=gsl_matrix_int_get(Lookup,i,0);
        lkIv=0.0;lk0=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            xj=gsl_matrix_get(X,j,0);
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            npart1=b1*lmet+xni+xj+uij;
            part1=b1*lmet+xi+xj+uij;
            lkIv+=iv*npart1-exp(npart1);
            lk0+=iv*part1-exp(part1);
        }
        deltaI=lkIv-lk0;
        
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        
        bias=b2*efi+b3*gci+b4*mapi;
        
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
    }
    
    mx/=ns;
    check/=ns;
    //printf("%f %f \n",mx,check);
    double C=(mx-b0);
    //C=0.0;
    double scale=exp(2.0*C/b1);
    //scale=1.0;
    double s0,s1,s2;
    for(i=0;i<ns;i++){
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        xi=gsl_matrix_get(X,i,0);
        s0*=scale;s1*=scale;s2*=scale;
        gsl_matrix_set(S,i,0,s0);gsl_matrix_set(S,i,1,s1);gsl_matrix_set(S,i,2,s2);
        xi-=C;
        gsl_matrix_set(X,i,0,xi);
    }
    
    gsl_matrix_scale (Metric,scale);
    double lscale=log(scale);
    gsl_matrix_add_constant (logMet, lscale);
    
    v=getlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,v);
    
    return accepted/ns;
}



double updateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
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
    double delta,met,accept;
    double lkIv,lk0;
    double s0,s1,s2,sn0,sn1,sn2,sj0,sj1,sj2;
    double nlmet,lmet,iv,v,npart1,part1;//,uij;
    double efi,efj,gci,gcj,mapi,mapj,bias;
    double accepted=0.0;
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
            //uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            bias=b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            npart1=b0+b1*nlmet+bias;
            part1=b0+b1*lmet+bias;
            lkIv+=iv*npart1-exp(npart1);
            lk0+=iv*part1-exp(part1);
        }
        delta=lkIv-lk0;
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
    }
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double updateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    
    double b1=gsl_matrix_get(beta,1,0);
    
    
    int i,j,k,kk,howmany;
    double delta,met,accept;
    double lkIv,lk0;
    double s0,s1,s2,sn0,sn1,sn2,sj0,sj1,sj2;
    double xi,xj,nlmet,lmet,iv,v,npart1,part1,uij;
    double accepted=0.0;
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
        xi=gsl_matrix_get(X,i,0);
        
        
        lkIv=0.0;lk0=0.0;
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
            
            lmet=gsl_matrix_get(logMet,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            npart1=b1*nlmet+xi+xj+uij;
            part1=b1*lmet+xi+xj+uij;
            lkIv+=iv*npart1-exp(npart1);
            lk0+=iv*part1-exp(part1);
        }
        delta=lkIv-lk0;
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
    }
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double HMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
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
    
    int i,j,k,kk,jL;
    double delta,met,met2,accept;
    double lkIv,lk0;
    double s0,s1,s2,sn0,sn1,sn2,sj0,sj1,sj2;
    double nlmet,lmet,iv,npart1,part1;
    double accepted=0.0;
    double pv0,pv1,pv2;
    double dw0,dw1,dw2;
    double current_K=0.0;double new_K=0.0;
    double efi,efj,gci,gcj,mapi,mapj;
    
    for(i=0;i<ns;i++){
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        pv0=gsl_ran_gaussian (r, 1.0);pv1=gsl_ran_gaussian (r, 1.0);pv2=gsl_ran_gaussian (r, 1.0);
        current_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        lk0=0.0;
        dw0=0.0;dw1=0.0;dw2=0.0;
        for(k=1;k<ns;k++){
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
            dw0+=(iv-exp(part1))*b1*(s0-sj0)/met2;
            dw1+=(iv-exp(part1))*b1*(s1-sj1)/met2;
            dw2+=(iv-exp(part1))*b1*(s2-sj2)/met2;
            lk0+=iv*part1-exp(part1);
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        sn0=s0;sn1=s1;sn2=s2;
        for(jL=1;jL<HL;jL++){
            sn0+=epsilon*pv0;sn1+=epsilon*pv1;sn2+=epsilon*pv2;
            dw0=0.0;dw1=0.0;dw2=0.0;
            for(k=1;k<ns;k++){
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
                dw0+=(iv-exp(part1))*b1*(sn0-sj0)/met2;
                dw1+=(iv-exp(part1))*b1*(sn1-sj1)/met2;
                dw2+=(iv-exp(part1))*b1*(sn2-sj2)/met2;
                
            }
            pv0+=epsilon*dw0;pv1+=epsilon*dw1;pv2+=epsilon*dw2;
            
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        new_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        
        lkIv=0.0;
        for(k=1;k<ns;k++){
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
            lkIv+=iv*npart1-exp(npart1);
        }
        
        delta=lkIv+new_K-lk0-current_K;
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            for(k=1;k<ns;k++){
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
                for(k=1;k<ns;k++){
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
    lk0=getlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,lk0);
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    
    return accepted/ns;
}


double HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*Metric,*logMet,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_matrix *Metricnew=gsl_matrix_calloc(ns,1);
    gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b1=gsl_matrix_get(beta,1,0);
    
    int i,j,k,kk,howmany,jL;
    double delta,met,met2,accept;
    double lkIv,lk0;
    double s0,s1,s2,sn0,sn1,sn2,sj0,sj1,sj2,xi,xj;
    double nlmet,lmet,iv,npart1,part1,uij;
    double accepted=0.0;
    double pv0,pv1,pv2;
    double dw0,dw1,dw2;
    double current_K=0.0;double new_K=0.0;
    
    
    for(i=0;i<ns;i++){
        
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        xi=gsl_matrix_get(X,i,0);
        pv0=gsl_ran_gaussian (r, 1.0);pv1=gsl_ran_gaussian (r, 1.0);pv2=gsl_ran_gaussian (r, 1.0);
        current_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        lk0=0.0;
        dw0=0.0;dw1=0.0;dw2=0.0;
        
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
            part1=b1*lmet+xi+xj+uij;
            dw0+=(iv-exp(part1))*b1*(s0-sj0)/met2;
            dw1+=(iv-exp(part1))*b1*(s1-sj1)/met2;
            dw2+=(iv-exp(part1))*b1*(s2-sj2)/met2;
            lk0+=iv*part1-exp(part1);
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        sn0=s0;sn1=s1;sn2=s2;
        for(jL=1;jL<HL;jL++){
            sn0+=epsilon*pv0;sn1+=epsilon*pv1;sn2+=epsilon*pv2;
            dw0=0.0;dw1=0.0;dw2=0.0;
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
                part1=b1*nlmet+xi+xj+uij;
                dw0+=(iv-exp(part1))*b1*(sn0-sj0)/met2;
                dw1+=(iv-exp(part1))*b1*(sn1-sj1)/met2;
                dw2+=(iv-exp(part1))*b1*(sn2-sj2)/met2;
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
            npart1=b1*nlmet+xi+xj+uij;
            lkIv+=iv*npart1-exp(npart1);
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
    
    lk0=getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,lk0);
    gsl_matrix_free(Metricnew);
    gsl_matrix_free(logMetnew);
    //printf("accepted %f\n",accepted/ns);
    
    return accepted/ns;
}



double updatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta,accept,accepted;
    double b0,betan,osuv,lkv,lk0;
    
    accepted=0.0;
    
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b0=gsl_matrix_get(beta,0,0);
    
    
    osuv=gsl_matrix_get(osu,4,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b0;
    
    gsl_matrix_set(betanew,0,0,betan);
    
    lkv=getlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
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


double updatebetaZero2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma)
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
    //    double lkX=getlkX(X,thetaX,ns,iSigma,beta);
    double delta=(betan*betan-beta0*beta0)*part0-2.0*(betan-beta0)*part2;
    delta*=-0.5/sigma;
    double v= gsl_matrix_get(lk,1,0);
    v+=delta;
    //    printf("betazero2 %f %f\n",delta,lkX-v);
    gsl_matrix_set(lk,1,0,v);
    //    printf("%f %f\n",sqrt(1.0/part1),part2/(sigma*part1));
    
    return 0.0;
}


double updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta,accept,accepted;
    double b1,betan,osuv,lkv,lk0;
    accepted=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b1=gsl_matrix_get(beta,1,0);
    
    osuv=gsl_matrix_get(osu,0,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b1;
    if(betan>0.0) return 0.0;
    
    gsl_matrix_set(betanew,1,0,betan);
    
    
    lkv=getlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
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


double updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta,accept,accepted;
    double b1,betan,osuv,lkv,lk0;
    
    accepted=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    b1=gsl_matrix_get(beta,1,0);
    osuv=gsl_matrix_get(osu,0,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b1;
    if(betan>0.0) return 0.0;
    
    gsl_matrix_set(betanew,1,0,betan);
    
    
    lkv=getlkI2(I,logMet,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;//printf("lkv delta (%f %f) %f %f\n",b1,betan,lkv,delta);
    
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


double updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta)
gsl_matrix *X,*lk,*osu,*thetaX,*iSigma,*beta;
int ns,nI;
gsl_matrix_int *Lookup;
{
    double xi,mapi;
    int i;
    
    double b0=gsl_matrix_get(beta,0,0);
    double b4=gsl_matrix_get(beta,4,0);
    
    gsl_matrix *Xstar=gsl_matrix_calloc(ns,1);
    gsl_matrix *Xstart=gsl_matrix_calloc(1,ns);
    
    for(i=0;i<ns;i++){
        xi=gsl_matrix_get(X,i,0);
        mapi=gsl_matrix_get(map,i,0);
        xi-=b0+b4*mapi;
        gsl_matrix_set(Xstar,i,0,xi);
    }
    
    gsl_matrix_transpose_memcpy(Xstart,Xstar);
    gsl_matrix *IHX=gsl_matrix_calloc(ns,1);
    gsl_matrix *vx=gsl_matrix_calloc(1,1);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, IH,Xstar,0.0, IHX);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Xstart, IHX,0.0, vx);
    
    double shape=0.001+0.5*(ns-2);
    double SSE=gsl_matrix_get(vx,0,0);
    double scale=0.5*SSE+0.001;
    double iscale=1.0/scale;
    double isigma=gsl_ran_gamma (r, shape, iscale);
    double sigmanew=1.0/isigma;
    //printf("%f %f\n",scale/(shape+1.0),sigmanew);
    //	if(sigmanew>0.5 || sigmanew<0.00000001) return 0.0;
    sigmanew=scale/(shape+1.0);
    gsl_matrix_set(thetaX,0,0,sigmanew);
    
    //double lkX2=-0.5*ns*log(sigmanew)-0.5*SSE/sigmanew;
    double lkX=getlkX(X,thetaX,ns,iSigma,beta);
    //printf("%f %f\n",lkX2,lkX);
    gsl_matrix_set(lk,1,0,lkX);
    
    gsl_matrix_free(IHX);
    gsl_matrix_free(vx);
    gsl_matrix_free(Xstar);
    gsl_matrix_free(Xstart);
    return 1.0;
    
}

double updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable)
gsl_matrix *U,*X,*I,*logMet,*beta,*osu,*thetaX,*lk;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    double lknv,lkv,deltaI,deltaU,delta,accept,accepted;
    double part1n,part1,iv,xi,xj,uij,unij,lmet,v;
    
    double sigmau=gsl_matrix_get(thetaX,2,0);
    double osuv=gsl_matrix_get(osu,2,0);
    double b1=gsl_matrix_get(beta,1,0);
    
    accepted=0;
    int i,j,k,howmany,kk;
    k=0;
    
    for(i=0;i<ns;i++){
        howmany=gsl_matrix_int_get(Lookup,i,0);
        xi=gsl_matrix_get(X,i,0);
        lknv=0.0;
        lkv=0.0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            if(i>=j) continue;
            
            xj=gsl_matrix_get(X,j,0);
            lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            uij=gsl_matrix_get(U,kk,0);
            unij=gsl_ran_gaussian (r, osuv);
            unij+=uij;
            
            part1n=b1*lmet+xi+xj+unij; //+bias;
            part1=b1*lmet+xi+xj+uij; //+bias;
            
            lknv=iv*part1n-exp(part1n);
            lkv=iv*part1-exp(part1);
            
            deltaI=lknv-lkv;
            deltaU=0.5/sigmau*(uij*uij-unij*unij);
            delta=deltaI+deltaU;
            
            if(delta>=0.0){
                gsl_matrix_set(U,kk,0,unij);
                v=gsl_matrix_get(lk,0,0);
                v+=deltaI;gsl_matrix_set(lk,0,0,v);
                v=gsl_matrix_get(lk,2,0);
                v+=deltaU;gsl_matrix_set(lk,2,0,v);
                accepted+=1.0;
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
                }
            }
            
        }
    }
    return accepted/nI;
}

double updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI)
gsl_matrix *U,*lk,*osu,*thetaX;
int ns,nI;
gsl_matrix_int *Lookup;
{
    double vu,lkU;
    
    
    double sigma=gsl_matrix_get(thetaX,2,0);	
    
    vu=gsl_matrix_get(lk,2,0);
    vu+=nI*0.5*log(sigma);
    vu*=-2.0*sigma;
    double shape=0.001+0.5*nI;
    double scale=0.5*vu+0.001;
    double iscale=1.0/scale;
    double isigma=gsl_ran_gamma (r, shape, iscale);
    double sigmanew=1.0/isigma;
    gsl_matrix_set(thetaX,2,0,sigmanew);
    
    lkU=-0.5*nI*log(sigmanew);
    lkU-=0.5*vu/sigmanew;
    gsl_matrix_set(lk,2,0,lkU);
    return 0.0;
}



double updateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu,*thetaX,*iSigma;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    int i;
    double b0=gsl_matrix_get(beta,0,0);
    double b4=gsl_matrix_get(beta,4,0);
    double xi; 
    double sigma=gsl_matrix_get(thetaX,0,0);
    double mapi;
    
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
    double usx=MVNsigmax*sqrt(sigma);
    double usy=MVNsigmay*sqrt(sigma);
    gsl_ran_bivariate_gaussian (r, usx,usy,MVNrho,&vx,&vy);
    
    vx=gsl_matrix_get(mu,0,0);
    vy=gsl_matrix_get(mu,1,0);
    
    gsl_matrix_set(beta,2,0,vx);
    gsl_matrix_set(beta,3,0,vy);
    
    double lkX=getlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,lkX);
    gsl_matrix_free(Xstar);
    gsl_matrix_free(mu);
    
    return 1.0;
}


double updatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta,accept,accepted;
    double b2,betan,osuv,lkv,lk0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    accepted=0.0;
    b2=gsl_matrix_get(beta,2,0);
    
    osuv=gsl_matrix_get(osu,5,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b2;
    gsl_matrix_set(betanew,2,0,betan);
    
    lkv=getlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
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



double updatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*logMet,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta,accept,accepted;
    double b3,betan,osuv,lkv,lk0;
    gsl_matrix *betanew=gsl_matrix_calloc(5,1);
    gsl_matrix_memcpy(betanew,beta);
    
    accepted=0.0;
    b3=gsl_matrix_get(beta,3,0);
    
    osuv=gsl_matrix_get(osu,6,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b3;
    gsl_matrix_set(betanew,3,0,betan);
    
    lkv=getlkI(I,logMet,X,U,betanew,ns,Lookup,ijTable);
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



void pram(
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
          double *result)
{
    int ns=*n;
    int nI=0;
    
    int rep_b=*repb;
    int rep_n=*repn;
    
    int i,j,k,howmany;
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
    
    
    //gsl_vector *collectorXs=gsl_vector_calloc(ns);
    //gsl_vector *collectorUs=gsl_vector_calloc(nI);
    gsl_vector *collectorb=gsl_vector_calloc(5);
    gsl_vector *collectorthetaXs=gsl_vector_calloc(3);
    gsl_vector *collectorLK=gsl_vector_calloc(3);
    
    int filter_factor=*thinning;
	int filter_n=rep_n/filter_factor;
	gsl_matrix *accept=gsl_matrix_calloc(filter_n,4);
	gsl_matrix *collectorS=gsl_matrix_calloc(filter_n,3*ns);
	//gsl_matrix *collectorX=gsl_matrix_calloc(filter_n,ns);
	//gsl_matrix *collectorU=gsl_matrix_calloc(filter_n,nI);
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

	
	gsl_rng_env_setup();
 	T = gsl_rng_default;
 	r = gsl_rng_alloc (T);
 	gsl_rng_set (r, 3);
 	
	double initial=gsl_rng_uniform (r);

	ubound=1.0; //atof(argv[29]);		  
	lbound=-1.0*ubound;
	
	HL=*nHL;
	epsilon=*vepsilon;

	k=0;
	for(i=0;i<ns;i++) {
		for(j=i+1;j<ns;j++){
			v=Contact[k];
            nI++;
			//if(v != 0 ) nI++;
			gsl_matrix_set(Jc,i,j,v);
	  		gsl_matrix_set(Jc,j,i,v);
            k++;
		}
	}
			
	for(i=0;i<ns;i++){
	   gsl_matrix_int_set(Lookup,i,0,0);
	}
	//Rprintf("Number of non-zero counts %d \n",nI);
	
	gsl_matrix *I=gsl_matrix_calloc(nI,1);
	gsl_matrix *U=gsl_matrix_calloc(nI,1);

	gsl_matrix *Metric=gsl_matrix_calloc(nI,1);
	gsl_matrix *logMet=gsl_matrix_calloc(nI,1);


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
        k++;
	}

	k=0;max_count=0.0;
	for(i=0;i<ns;i++){	
		for(j=i+1;j<ns;j++){
			v=gsl_matrix_get(Jc,i,j);
			if(max_count<v) max_count=v;
			//if(v != 0.0){
			//	Mi+=v;
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
			//}
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


	
	getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
	
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


double v0,v1,v2,v3;	
double scale,b1;
int fi=0;
//double adjust=0.1;

//int nn=0;
//int howoften=0;



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
	sv=getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);	
	gsl_matrix_set(lk,0,0,sv);
	sv= getlkX(X,thetaX,ns,iSigma,beta);
	gsl_matrix_set(lk,1,0,sv);
	sv=getlkU(U,thetaX,nI);
    gsl_matrix_set(lk,2,0,sv);
    
    gsl_matrix_set(beta,4,0,1.0);//printf("MAP is set off\n");
    	//printf("%f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
	//Rprintf("INITIALIZATION \n");
    //=====================================
    gsl_matrix_set(beta,1,0,0.0);
	for(i=0;i<10000;i++){
		updateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
		updatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		updatebetaEFL(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		updatebetaGC(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		/*if (i%1000==0) {
			Rprintf("After initialization %d beta0 %f beta1 %f ",i,gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0));
            Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
			Rprintf("log-likelihood %f\n",gsl_matrix_get(lk,0,0));
		}*/

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

    getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    v=getlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,v);
    v=getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);

    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
	Rprintf("\n====================================\n");
	Rprintf("Jumpinp rules in Proposals\n");
	Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
	Rprintf("beta1 : %f \n",gsl_matrix_get(osu,0,0));
	Rprintf("X : %f \n",gsl_matrix_get(osu,1,0));
	Rprintf("U : %f \n",gsl_matrix_get(osu,2,0));

    /*
    int nn;
    int howoften;
for(j=0;j<3;j++){
	fi=0;
	if(j<10) {nn=10000;howoften=1;}
	else {nn =10000;howoften=1;}
	for(i=0;i<nn;i++){		
		v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);		
		v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);		
		v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
		v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
		v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
		v=updateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);

		if(i%howoften==0){
			gsl_matrix_set(accept,fi,0,v0);
			gsl_matrix_set(accept,fi,1,v1);
			gsl_matrix_set(accept,fi,2,v2);
			gsl_matrix_set(accept,fi,3,v3);
			fi+=1;	
		}
	}
	
	v0=0.0;v1=0.0;v2=0.0;v3=0.0;
	for(i=0;i<10000;i++){
		v0+=gsl_matrix_get(accept,i,0);
		v1+=gsl_matrix_get(accept,i,1);
		v2+=gsl_matrix_get(accept,i,2);
		v3+=gsl_matrix_get(accept,i,3);

	}

	Rprintf("Acceptance rate after %d iterations\n",(j+1)*10000);
	Rprintf("S : %f\n",v0/10000.0);
	Rprintf("beta1 : %f \n",v1/10000.0);
	Rprintf("X : %f\n",v2/10000.0);
	Rprintf("U : %f\n",v3/10000.0);	
}
	Rprintf("========Evaluation has ended\n");
     */

	fi=0;
	Rprintf("\nBurning has started\n");
	for(i=0;i<rep_b;i++){		
		v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);		
		v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);		
		v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
		v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
		v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
		v=updateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
		if ((i+1)%rep_b==0) {
			Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
			Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
			Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0)+gsl_matrix_get(lk,1,0)+gsl_matrix_get(lk,2,0));
            Rprintf("X %f\n",gsl_matrix_get(X,0,0));
		}
	}

	Rprintf("Burning has ended\n");

	Rprintf("\nPosterior sampling has started\n");
	for(i=0;i<rep_n;i++){		
		v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);		
		v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);		
		v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
		v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
		v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
		v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
		v=updateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
		
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
		
			//gsl_matrix_get_col (collectorXs, X, 0);
			//gsl_matrix_set_row (collectorX, fi, collectorXs);
			
			//gsl_matrix_get_col (collectorUs, U, 0);
			//gsl_matrix_set_row (collectorU, fi, collectorUs);
			
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
			Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
			//Rprintf("log-likelihood %f \n",gsl_matrix_get(lk,0,0)+gsl_matrix_get(lk,1,0)+gsl_matrix_get(lk,2,0));
				
		}
 
}

	
	time(&stop);
    diff = difftime(stop, start);
 
    Rprintf("==========Summary Report================\n");
    Rprintf("CPU TIME IS %f\n",diff/60.0);
    
    double mean;
	for(i=0;i<4;i++) {
		mean=0.0;
		for(j=0;j<filter_n;j++){
			mean+=gsl_matrix_get(accept,j,i);
		}
		mean/=filter_n;
		if(i==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
		if(i==1) Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
		if(i==2) Rprintf("X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),mean);
		if(i==3) Rprintf("U : proposal jump %f acceptance rate%f\n",gsl_matrix_get(osu,2,0),mean);

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
			result[k]=gsl_matrix_get(collectorthetaX,i,0);k++;
			result[k]=gsl_matrix_get(collectorthetaX,i,2);k++;
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
	gsl_vector_free(collectorb);
	gsl_vector_free(collectorthetaXs);
	gsl_vector_free(collectorLK);


	return;

}





void pram2(
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
    
    getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    
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
    
    sv= getlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    for(i=0;i<10000;i++){
        updateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
        updatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
        updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
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
    
    getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
    
    sv= getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    sv= getlkX(X,thetaX,ns,iSigma,beta);
    gsl_matrix_set(lk,1,0,sv);
    sv= getlkU(U,thetaX,nI);
    gsl_matrix_set(lk,2,0,sv);
    gsl_matrix_set(beta,4,0,1.0);
    
    Rprintf("Initial2 : beta1 %f %f ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0));
    Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
    
    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
    Rprintf("\n====================================\n");
    int check=0;
    if(*gear == 1){
        Rprintf("Jumping rules in Proposals\n");
        Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
        Rprintf("beta1 : osu %f \n",gsl_matrix_get(osu,0,0));
        Rprintf("X : osu %f \n",gsl_matrix_get(osu,1,0));
        Rprintf("U : osu %f \n",gsl_matrix_get(osu,2,0));
        
        Rprintf("Looking for jumping rules...\n");
        
        //for(j=0;j<7;j++){
        j=0;
        while(check < 4){
            fi=0;ufi=0;
            check=0;j++;
            if(j==30) break;
            nn=10000;howoften=1;
            //if(j<10) {nn=10000;howoften=1;adjust=0.03;printf("adjust %f\n",adjust);}
            //else {nn =10000;howoften=1;adjust=0.01;}
            for(i=0;i<nn;i++){
                v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
                v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
                v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
                v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
                
                if (i%10000==0) {
                    Rprintf("Adjusting : beta1 %f ",gsl_matrix_get(beta,1,0));
                    Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
                    Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
                }
                if(i%howoften==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    gsl_matrix_set(accept,fi,2,v2);
                    gsl_matrix_set(accept,fi,3,v3);
                    fi+=1;
                }
            }
            
            v0=0.0;v1=0.0;v2=0.0;v3=0.0;
            for(i=0;i<10000;i++){
                v0+=gsl_matrix_get(accept,i,0);
                v1+=gsl_matrix_get(accept,i,1);
                v2+=gsl_matrix_get(accept,i,2);
                v3+=gsl_matrix_get(accept,i,3);
            }
            
            Rprintf("j-th adjustment %d\n",j);
            Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,v0/10000.0);
            Rprintf("beta1 : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),v1/10000.0);
            Rprintf("X : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),v2/10000.0);
            Rprintf("U : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),v3/10000.0);
            
            if(v0 < 5000) {
                epsilon*=0.1;
            }else if(v0 < 7000 & v0 >= 5000){
                epsilon*=0.3;
            }else if(v0 < 9000 & v0 >= 7000){
                epsilon*=0.5;
            }else if(v0 < 9800 & v0 >= 9000){
                epsilon*=0.8;
            }else if(v0 > 9900){
                epsilon*=1.1;
            }else check++;
            
            if(v1 < 1200) {
                mv=gsl_matrix_get(osu,0,0);v=0.3*mv;gsl_matrix_set(osu,0,0,v);
            }else if(v1 > 7000){
                mv=gsl_matrix_get(osu,0,0);v=3.0*mv;gsl_matrix_set(osu,0,0,v);
            }else if(v1 > 4000 & v1 <= 7000){
                mv=gsl_matrix_get(osu,0,0);v=2.0*mv;gsl_matrix_set(osu,0,0,v);
            }else if(v1 > 3000 & v1 <= 4000){
                mv=gsl_matrix_get(osu,0,0);v=1.5*mv;gsl_matrix_set(osu,0,0,v);
            }else if(v1 > 2200 & v1 <= 3000){
                mv=gsl_matrix_get(osu,0,0);v=1.1*mv;gsl_matrix_set(osu,0,0,v);
            }else if(v1 >= 1000 & v1 < 1900){
                mv=gsl_matrix_get(osu,0,0);v=0.9*mv;gsl_matrix_set(osu,0,0,v);
            }else check++;
            
            
            if(v2 < 1200) {
                mv=gsl_matrix_get(osu,1,0);v=0.7*mv;gsl_matrix_set(osu,1,0,v);
            }else if(v2 > 7000){
                mv=gsl_matrix_get(osu,1,0);v=3.0*mv;gsl_matrix_set(osu,1,0,v);
            }else if(v2 > 4000 & v2 <= 7000){
                mv=gsl_matrix_get(osu,1,0);v=2.0*mv;gsl_matrix_set(osu,1,0,v);
            }else if(v2 > 3000 & v2 <= 4000){
                mv=gsl_matrix_get(osu,1,0);v=1.5*mv;gsl_matrix_set(osu,1,0,v);
            }else if(v2 > 2000 & v2 <= 3000){
                mv=gsl_matrix_get(osu,1,0);v=1.1*mv;gsl_matrix_set(osu,1,0,v);
            }else if(v2 >= 1200 & v2 < 1900){
                mv=gsl_matrix_get(osu,1,0);v=0.9*mv;gsl_matrix_set(osu,1,0,v);
            }else check++;
            
            
            if(v3 < 1200) {
                mv=gsl_matrix_get(osu,2,0);v=0.7*mv;gsl_matrix_set(osu,2,0,v);
            }else if(v3 > 7000){
                mv=gsl_matrix_get(osu,2,0);v=3.0*mv;gsl_matrix_set(osu,2,0,v);
            }else if(v3 > 4000 & v3 <= 7000){
                mv=gsl_matrix_get(osu,2,0);v=2.0*mv;gsl_matrix_set(osu,2,0,v);
            }else if(v3 > 3000 & v3 <= 4000){
                mv=gsl_matrix_get(osu,2,0);v=1.5*mv;gsl_matrix_set(osu,2,0,v);
            }else if(v3 > 2000 & v3 <= 3000){
                mv=gsl_matrix_get(osu,2,0);v=1.1*mv;gsl_matrix_set(osu,2,0,v);
            }else if(v3 >= 1000 & v3 < 1900){
                mv=gsl_matrix_get(osu,2,0);v=0.9*mv;gsl_matrix_set(osu,2,0,v);
            }else check++;
            
            //--------------MODIFIED on March 10.
            /*
             if(v3 < 1000) {
             mv=gsl_matrix_get(osu,2,0);v=0.9*mv;gsl_matrix_set(osu,2,0,v);
             }else if(v3 > 7000){
             mv=gsl_matrix_get(osu,2,0);v=3.0*mv;gsl_matrix_set(osu,2,0,v);
             }else if(v3 > 4000 & v3 <= 7000){
             mv=gsl_matrix_get(osu,2,0);v=2.0*mv;gsl_matrix_set(osu,2,0,v);
             }else if(v3 > 3000 & v3 <= 4000){
             mv=gsl_matrix_get(osu,2,0);v=1.5*mv;gsl_matrix_set(osu,2,0,v);
             }else if(v3 > 2000 & v3 <= 3000){
             mv=gsl_matrix_get(osu,2,0);v=1.2*mv;gsl_matrix_set(osu,2,0,v);
             }else if(v3 >= 1200 & v3 <= 2000){
             mv=gsl_matrix_get(osu,2,0);v=1.1*mv;gsl_matrix_set(osu,2,0,v);
             }else check++;
             */
            
            Rprintf("%d %d\n",j,check);
            //=================================================================================
            if(check < 3){
                gsl_matrix_set(beta,0,0,initial*log(max_count));
                v=gsl_matrix_get(beta,0,0);
                for(i=0;i<ns;i++){
                    gsl_matrix_set(S,i,0,line_t*i/(ns-1));
                    gsl_matrix_set(S,i,1,0.0);
                    gsl_matrix_set(S,i,2,0.0);
                    //line_t=gsl_ran_gaussian (r, 0.1*ubound);
                    
                }
                getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                gsl_matrix_set(beta,1,0,-2.0*initial);
                sv= getlkI(I,logMet,X,U,beta,ns,Lookup,ijTable);
                gsl_matrix_set(lk,0,0,sv);
                for(i=0;i<10000;i++){
                    updateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                    updatebetaZero(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                    updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
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
                
                getDist(S,Metric,logMet,ns,nI,Lookup,ijTable);
                
                sv= getlkI2(I,logMet,X,U,beta,ns,Lookup,ijTable);
                gsl_matrix_set(lk,0,0,sv);
                sv= getlkX(X,thetaX,ns,iSigma,beta);
                gsl_matrix_set(lk,1,0,sv);
                sv= getlkU(U,thetaX,nI);
                gsl_matrix_set(lk,2,0,sv);
                gsl_matrix_set(beta,4,0,1.0);
            }
        }
        
        fi=0;ufi=0;
        Rprintf("\nBurning has started\n");
        for(i=0;i<rep_b;i++){
            v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
            v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
            v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
            v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
            //            v=zupdateBias(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,thetaX,iSigma);
            if ((i+1)==rep_b) {
                Rprintf("After burn-in : beta1 %f ",gsl_matrix_get(beta,1,0));
                Rprintf("sigma2_x %f sigma2_u %f \n",gsl_matrix_get(thetaX,0,0),gsl_matrix_get(thetaX,2,0));
                Rprintf("log-likelihood %f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
            }
        }
        
        Rprintf("Burning has ended\n");
        
        Rprintf("\nPosterior sampling has started\n");
        for(i=0;i<rep_n;i++){
            v2=updateX2(X,U,I,logMet,beta,lk,ns,nI,osu,thetaX,Lookup,ijTable,iSigma,S,Metric);
            v0=HMCupdateS2(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            v1=updatebetaOne2(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            v=updateSigma(X,lk,osu,thetaX,Lookup,ns,nI,iSigma,beta);
            v3=updateU(U,X,I,logMet,beta,osu,thetaX,lk,Lookup,ns,nI,ijTable);
            v=updateSigmaU(U,lk,osu,thetaX,Lookup,ns,nI);
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
            v0=HMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            v1=updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
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
            v0=HMCupdateS(S,Metric,logMet,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
            v1=updatebetaOne(I,logMet,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
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


