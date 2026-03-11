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
#include <gsl/gsl_multifit.h>


void bngetDist(S,ns,nI,Lookup,ijTable)
gsl_matrix *S;
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
            gsl_vector_set(dd,kk,v);
            //gsl_vector_set(logMet,kk,log(v));
            //printf("kk %d i %d j %d %f\n",kk,i,j,gsl_vector_get(dd,kk));
        }
    }
    return;
}

double bngetlkI(I,X,U,beta,ns,Lookup,ijTable)
gsl_matrix *I,*X,*U,*beta;
int ns;
gsl_matrix_int *Lookup,*ijTable;
{
    int i,j,k,kk,howmany;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    //double b3=gsl_matrix_get(beta,3,0);
    //double b4=gsl_matrix_get(beta,4,0);
    double retv=0.0,lmet=0.0,iv=0.0,part1=0.0,lambda=0.0;
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0,bias=0.0;
    
    
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
            
            //lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
            //printf("bias %f\n",bias);
            part1=anchor+gsl_vector_get(loglambda,kk)+bias;
            lambda=exp(part1);
            //retv+=iv*part1-log(exp(lambda)-1);
            
            if(lambda<10.0) retv+=iv*part1-log(exp(lambda)-1.0);
            else retv+=iv*part1-lambda;
            
            
        }
    }
    return retv;
}


double bngetlkI3(I,X,U,beta,ns,Lookup,ijTable,loglambdastar)
gsl_matrix *I,*X,*U,*beta;
int ns;
gsl_matrix_int *Lookup,*ijTable;
gsl_vector *loglambdastar;
{
    int i,j,k,kk,howmany;
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    //double b3=gsl_matrix_get(beta,3,0);
    //double b4=gsl_matrix_get(beta,4,0);
    double retv=0.0,lmet=0.0,iv=0.0,part1=0.0,lambda=0.0;
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0,bias=0.0;
    
    
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
            
            //lmet=gsl_matrix_get(logMet,kk,0);
            iv=gsl_matrix_get(I,kk,0);
            
            bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
            //printf("bias %f\n",bias);
            part1=anchor+gsl_vector_get(loglambdastar,kk)+bias;
            lambda=exp(part1);
            //retv+=iv*part1-log(exp(lambda)-1);
            
            if(lambda<10.0) retv+=iv*part1-log(exp(lambda)-1.0);
            else retv+=iv*part1-lambda;
            
        }
    }
    return retv;
}


double bnupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_vector *newdd=gsl_vector_calloc(ns);
    gsl_vector *newloglambda=gsl_vector_calloc(ns);
    //gsl_matrix *logMetnew=gsl_matrix_calloc(ns,1);
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    //double b3=gsl_matrix_get(beta,3,0);
    //double b4=gsl_matrix_get(beta,4,0);
    
    int i,j,k,kk,howmany;
    double dv=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    double iv=0.0,v=0.0,uij=0.0;
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    
    
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0,bias=0.0;
    double accepted=0.0;double delta=0.0;
    double osuv=0.003;//gsl_matrix_get(osu,1,0);
    double xi,yi,yerr,x0,y0;
    
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
            
            dv=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            gsl_vector_set(newdd,k,dv);
            
            iv=gsl_matrix_get(I,kk,0);
            gsl_bspline_eval(dv, B, Bw);
            gsl_multifit_linear_est(B, Cx, Cov, &yi, &yerr);
            
            gsl_vector_set(newloglambda,k,yi);
            y0=gsl_vector_get(loglambda,kk);
            
            bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
            npart1=anchor+yi+bias;
            part1=anchor+y0+bias;
            
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
                dv=gsl_vector_get(newdd,k);
                gsl_vector_set(dd,kk,dv);
                yi=gsl_vector_get(newloglambda,k);
                gsl_vector_set(loglambda,kk,yi);
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
                    dv=gsl_vector_get(newdd,k);
                    gsl_vector_set(dd,kk,dv);
                    yi=gsl_vector_get(newloglambda,k);
                    gsl_vector_set(loglambda,kk,yi);
                }
                accepted+=1.0;
            }
        }
        //printf("update S (newS=%f oldS=%f) delta=%f accepted %f\n",sn1,s1,delta,accepted);
    }
    gsl_vector_free(newdd);
    gsl_vector_free(newloglambda);
    return accepted/ns;
}



double bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *S,*X,*U,*I,*thetaX;
gsl_matrix *beta,*lk,*osu;
int ns,nI;
gsl_matrix_int *Lookup,*ijTable;
{
    gsl_vector *newdd=gsl_vector_calloc(ns);
    gsl_vector *newloglambda=gsl_vector_calloc(ns);
    
    double b0=gsl_matrix_get(beta,0,0);
    double b1=gsl_matrix_get(beta,1,0);
    double b2=gsl_matrix_get(beta,2,0);
    //double b3=gsl_matrix_get(beta,3,0);
    //double b4=gsl_matrix_get(beta,4,0);
    
    int i,j,k,kk,howmany,jL;
    double met=0.0,met2=0.0,accept=0.0;
    double lkIv=0.0,lk0=0.0;
    double s0=0.0,s1=0.0,s2=0.0,sn0=0.0,sn1=0.0,sn2=0.0,sj0=0.0,sj1=0.0,sj2=0.0;
    double factor=1.0;
    
    
    double npart1=0.0,part1=0.0;
    double nlambda=1.0,lambda=0.0;
    double lmet=0.0,nlmet=0.0;
    
    double accepted=0.0,delta=0.0;
    double iv=0.0,dv=0.0;
    
    double pv0=0.0,pv1=0.0,pv2=0.0;
    double dw0=0.0,dw1=0.0,dw2=0.0;
    double current_K=0.0,new_K=0.0;
    double efi=0.0,efj=0.0,gci=0.0,gcj=0.0,mapi=0.0,mapj=0.0,bias=0.0;
    double xi,yi,yerr,x0,y0,v;
    int check=0;
    
    for(i=0;i<ns;i++){
        s0=gsl_matrix_get(S,i,0);s1=gsl_matrix_get(S,i,1);s2=gsl_matrix_get(S,i,2);
        pv0=gsl_ran_gaussian (r, 1.0);pv1=gsl_ran_gaussian (r, 1.0);pv2=gsl_ran_gaussian (r, 1.0);
        //printf("pv %f %f %f\n",pv0,pv1,pv2);
        current_K=-0.5*(pv0*pv0+pv1*pv1+pv2*pv2);
        efi=gsl_matrix_get(efl,i,0);
        gci=gsl_matrix_get(gc,i,0);
        mapi=gsl_matrix_get(map,i,0);
        lk0=0.0;
        dw0=0.0;dw1=0.0;dw2=0.0;
        howmany=gsl_matrix_int_get(Lookup,i,0);
        check=0;
        for(k=1;k<howmany+1;k++){
            kk=gsl_matrix_int_get(Lookup,i,k);
            j=gsl_matrix_int_get(ijTable,i,k);
            
            efj=gsl_matrix_get(efl,j,0);
            gcj=gsl_matrix_get(gc,j,0);
            mapj=gsl_matrix_get(map,j,0);
            
            dv=gsl_vector_get(dd,kk);
            iv=gsl_matrix_get(I,kk,0);
            bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
            y0=gsl_vector_get(loglambda,kk);
            
            
            //bias=b2*(efi+efj)+b3*(gci+gcj)+b4*(mapi+mapj);
            part1=anchor+y0+bias;
            
            
            lambda=exp(part1);
            if(lambda<10.0) {
                lk0+=iv*part1-log(exp(lambda)-1.0);
                factor=exp(lambda)/(exp(lambda)-1.0);
            }else {
                lk0+=iv*part1-lambda;
                factor=1.0;
            }
            
            dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/dv;
            dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/dv;
            dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/dv;
            
        }
        
        pv0+=0.5*epsilon*dw0;pv1+=0.5*epsilon*dw1;pv2+=0.5*epsilon*dw2;
        //printf("pv %f %f %f\n",pv0,pv1,pv2);
        sn0=s0;sn1=s1;sn2=s2;
        for(jL=1;jL<HL;jL++){
            sn0+=epsilon*pv0;sn1+=epsilon*pv1;sn2+=epsilon*pv2;
            dw0=0.0;dw1=0.0;dw2=0.0;
            for(k=1;k<howmany+1;k++){
                if(check==0){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    j=gsl_matrix_int_get(ijTable,i,k);
                    
                    sj0=gsl_matrix_get(S,j,0);
                    sj1=gsl_matrix_get(S,j,1);
                    sj2=gsl_matrix_get(S,j,2);
                    efj=gsl_matrix_get(efl,j,0);
                    gcj=gsl_matrix_get(gc,j,0);
                    mapj=gsl_matrix_get(map,j,0);
                    
                    dv=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
                    iv=gsl_matrix_get(I,kk,0);
                    bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
                    if(dv > MaxKnot) {
                        check=1;
                        continue;
                    }
                    
                    //printf("dv %f sn %f %f %f\n",dv,sn0,sn1,sn2);
                    //printf("MX KNOT DD %f %f\n",MaxKnot,dv);
                    gsl_bspline_eval(dv, B, Bw);
                    gsl_multifit_linear_est(B, Cx, Cov, &yi, &yerr);
                    
                    part1=anchor+yi+bias;
                    
                    lambda=exp(part1);
                    if(lambda<10.0) factor=exp(lambda)/(exp(lambda)-1.0);
                    else factor=1.0;
                    dw0+=(iv-factor*lambda)*b1*(sn0-sj0)/dv;
                    dw1+=(iv-factor*lambda)*b1*(sn1-sj1)/dv;
                    dw2+=(iv-factor*lambda)*b1*(sn2-sj2)/dv;
                }
                
            }
            if(check==0) pv0+=epsilon*dw0;pv1+=epsilon*dw1;pv2+=epsilon*dw2;
            
        }
        if(check==1) continue;
        
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
            
            dv=sqrt((sn0-sj0)*(sn0-sj0)+(sn1-sj1)*(sn1-sj1)+(sn2-sj2)*(sn2-sj2));
            
            gsl_vector_set(newdd,k,dv);
            iv=gsl_matrix_get(I,kk,0);
            bias=b0*(efi+efj)+b1*(gci+gcj)+b2*(mapi+mapj);
            //printf("MX KNOT DD %f %f\n",MaxKnot,dv);
            gsl_bspline_eval(dv, B, Bw);
            gsl_multifit_linear_est(B, Cx, Cov, &yi, &yerr);
            gsl_vector_set(newloglambda,k,yi);
            
            npart1=anchor+yi+bias;
            nlambda=exp(npart1);
            
            if(nlambda<10.0) lkIv+=iv*npart1-log(exp(nlambda)-1.0);
            else lkIv+=iv*npart1-nlambda;
        }
        
        delta=lkIv+new_K-lk0-current_K;
        //printf("HMCUpdate %d (%f %f %f)<-(%f %f %f) %f\n",i,sn0,sn1,sn2,s0,s1,s2,delta);
        if(delta>=0.0){
            gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
            for(k=1;k<howmany+1;k++){
                kk=gsl_matrix_int_get(Lookup,i,k);
                dv=gsl_vector_get(newdd,k);
                gsl_vector_set(dd,kk,dv);
                yi=gsl_vector_get(newloglambda,k);
                gsl_vector_set(loglambda,kk,yi);
            }
            accepted+=1.0;
        }
        else{
            accept=gsl_rng_uniform (r);
            if(accept<exp(delta)){
                gsl_matrix_set(S,i,0,sn0);gsl_matrix_set(S,i,1,sn1);gsl_matrix_set(S,i,2,sn2);
                for(k=1;k<howmany+1;k++){
                    kk=gsl_matrix_int_get(Lookup,i,k);
                    dv=gsl_vector_get(newdd,k);
                    gsl_vector_set(dd,kk,dv);
                    yi=gsl_vector_get(newloglambda,k);
                    gsl_vector_set(loglambda,kk,yi);
                }
                accepted+=1.0;
            }
        }
    }
    lk0=bngetlkI(I,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,lk0);
    //printf("lk0 %f\n",lk0);
    gsl_vector_free(newdd);
    gsl_vector_free(newloglambda);
    return accepted/ns;
}



double bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S)
gsl_matrix *I,*X,*U,*beta,*lk,*osu,*S;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    //gsl_multifit_linear_workspace *mw;
    gsl_vector *Bxi=gsl_vector_alloc(ncoeffs);
    gsl_vector *Cxstar=gsl_vector_alloc(ncoeffs);
    gsl_vector *loglambdastar=gsl_vector_alloc(nI);
    //gsl_vector *Cxsort=gsl_vector_alloc(ncoeffs);
    gsl_vector *breaks=gsl_vector_alloc(nbreak);
    
    double accepted=0.0,xi,v;
    double accept=0.0;
    double chisq,delta1,delta2;
    double jump=gsl_matrix_get(osu,0,0);
    double ddsort[nI];
    
    int i,j;
    
    //mw = gsl_multifit_linear_alloc(nI, ncoeffs);
    
    MaxKnot  = 0.0;
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        if(v > MaxKnot) MaxKnot=v;
        
    }
    MaxKnot*=1.1;
    
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        ddsort[i]=v;
        
    }
    
    gsl_sort (ddsort, 1, nI);
    
    double bk1= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.05);
    double bk2= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.10);
    double bk3= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.15);
    double bk4= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.35);
    double bk5= gsl_stats_median_from_sorted_data (ddsort,1, nI);
    double bk6= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.65);
    double bk7= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.75);
    
    //printf("%f %f %f\n",lowerq,median,upperq);
    
    gsl_vector_set(breaks,0,0.0);
    gsl_vector_set(breaks,1,bk1);
    gsl_vector_set(breaks,2,bk2);
    gsl_vector_set(breaks,3,bk3);
    gsl_vector_set(breaks,4,bk4);
    gsl_vector_set(breaks,5,bk5);
    gsl_vector_set(breaks,6,bk6);
    gsl_vector_set(breaks,7,bk7);
    gsl_vector_set(breaks,8,MaxKnot);
    
    gsl_bspline_knots (breaks, Bw);
    //gsl_bspline_knots_uniform(0.0, MaxKnot, Bw);
    
    /* construct the fit matrix X */
    for (i = 0; i < nI; ++i)
    {
        double xi = gsl_vector_get(dd, i);
        gsl_bspline_eval(xi, B, Bw);
        for (j = 0; j < ncoeffs; ++j)
        {
            double Bj = gsl_vector_get(B, j);
            gsl_matrix_set(Bx, i, j, Bj);
        }
    }

    double v0,v1,v2;
    gsl_vector_set(Cxstar,0,CEIL);
    for(i=1;i<ncoeffs;i++){
        v1=gsl_ran_gaussian (r, jump);
        v1+=gsl_vector_get(Cx,i);
        if(v1 > CEIL){
            //gsl_vector_free(Bxi);gsl_vector_free(Cxstar);gsl_vector_free(loglambdastar);
            //return 0.0;
            //gsl_vector_set(Cxstar,i,CEIL);
            v1=CEIL-0.03*gsl_rng_uniform(r);
        }
        if(i==1){
            v2=gsl_vector_get(Cx,(i+1));
            if(v1 > v2) gsl_vector_set(Cxstar,i,v1);
            else gsl_vector_set(Cxstar,i,v2);
            
        }else if(i==(ncoeffs-1)){
            v0=gsl_vector_get(Cx,(i-1));
            if(v0 > v1) gsl_vector_set(Cxstar,i,v1);
            else gsl_vector_set(Cxstar,i,v0);
            
        }else{
            v0=gsl_vector_get(Cx,(i-1));
            v2=gsl_vector_get(Cx,(i+1));
            if((v0 >= v1) && (v1>= v2)) gsl_vector_set(Cxstar,i,v1);
            else{
                if(v0 < v1) gsl_vector_set(Cxstar,i,v0);
                else if(v2 > v1) gsl_vector_set(Cxstar,i,v2);
            }
        }
    }
    
    //gsl_multifit_linear(Bx, loglambda, Cxstar, Cov, &chisq, mw);
    
    
    //v0=gsl_vector_get(Cxstar,0);v1=gsl_vector_get(Cxstar,1);v2=gsl_vector_get(Cxstar,2);
    //v3=gsl_vector_get(Cxstar,3);v4=gsl_vector_get(Cxstar,4);v5=gsl_vector_get(Cxstar,5);
    //v5=gsl_vector_get(Cxstar,6);
    //printf("proposed spline %f %f %f %f %f %f %f\n",v0,v1,v2,v3,v4,v5,v6);
    //printf("current spline %f %f %f %f %f\n",gsl_vector_get(Cx,0),gsl_vector_get(Cx,1),gsl_vector_get(Cx,2),gsl_vector_get(Cx,3),gsl_vector_get(Cx,4));
    
    
    for(i=0;i<nI;i++){
        v=0.0;
        gsl_matrix_get_row (Bxi, Bx, i);
        gsl_vector_mul (Bxi, Cxstar);
        for(j=0;j<ncoeffs;j++){
            v+=gsl_vector_get(Bxi,j);
        }
        gsl_vector_set(loglambdastar,i,v);
    }
    
    delta2=bngetlkI3(I,X,U,beta,ns,Lookup,ijTable,loglambdastar);
    delta1=gsl_matrix_get(lk,0,0);
    //getlkI3(I,X,U,beta,ns,Lookup,ijTable,loglambda);
    double delta=delta2-delta1;
    
    //printf("delta %f (%f %f) %f\n",delta2,delta1,gsl_matrix_get(lk,0,0),delta);
    if(delta>=0.0){
        accepted=1.0;
    }
    else{
        accept=gsl_rng_uniform (r);
        if(accept<exp(delta)){
            accepted=1.0;
        }
    }
    
    if(accepted==1.0){
        //printf("accepted spline %f %f %f %f %f MaxKnot %f\n",v0,v1,v2,v3,v4,MaxKnot);
        //printf("S %f %f %f\n",gsl_matrix_get(S,0,0),gsl_matrix_get(S,0,1),gsl_matrix_get(S,0,2));
        gsl_vector_memcpy(Cx,Cxstar);
        gsl_vector_memcpy(loglambda,loglambdastar);
        gsl_matrix_set(lk,0,0,delta2);
    }
    
    gsl_vector_free(Bxi);
    gsl_vector_free(Cxstar);
    gsl_vector_free(loglambdastar);
    gsl_vector_free(breaks);
    
    return accepted;
}



double bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b2=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(3,1);
    gsl_matrix_memcpy(betanew,beta);
    
    
    b2=gsl_matrix_get(beta,0,0);
    osuv=gsl_matrix_get(osu,1,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b2;
    gsl_matrix_set(betanew,0,0,betan);
    
    lkv=bngetlkI(I,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;
    double priordelta=0.5*b2*b2*prior-0.5*betan*betan*prior;
    delta+=priordelta;
    
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
        else gsl_matrix_set(lk,0,0,lk0);
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}



double bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable)
gsl_matrix *I,*X,*U,*beta,*lk,*osu;
gsl_matrix_int *Lookup,*ijTable;
int ns,nI;
{
    double delta=0.0,accept=0.0,accepted=0.0;
    double b3=0.0,betan=0.0,osuv=0.0,lkv=0.0,lk0=0.0;
    gsl_matrix *betanew=gsl_matrix_calloc(3,1);
    gsl_matrix_memcpy(betanew,beta);
    
    b3=gsl_matrix_get(beta,1,0);
    osuv=gsl_matrix_get(osu,2,0);
    betan=gsl_ran_gaussian (r, osuv);betan+=b3;
    gsl_matrix_set(betanew,1,0,betan);
    
    lkv=bngetlkI(I,X,U,betanew,ns,Lookup,ijTable);
    lk0=gsl_matrix_get(lk,0,0);
    delta=lkv-lk0;
    double priordelta=0.5*b3*b3*prior-0.5*betan*betan*prior;
    delta+=priordelta;
    
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
        else gsl_matrix_set(lk,0,0,lk0);
    }
    
    gsl_matrix_free(betanew);
    return accepted;
}



void bn(
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
    //gsl_matrix *U2R = gsl_matrix_calloc(1,nu2r);
    
    int rep_b=*repb;
    int rep_n=*repn;
    
    int i,j,k,howmany,ii,jj,kk,uk;
    double v,Mi,line_t;
    double max_count;
    
    //MaxKnot=1.0;
    CEIL=3.0;
    anchor=0.0;
    
    ncoeffs = 11;
    nbreak = ncoeffs-2;
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
    gsl_matrix *beta=gsl_matrix_calloc(3,1);
    gsl_matrix *thetaX=gsl_matrix_calloc(3,1);
    gsl_matrix *osu=gsl_matrix_calloc(7,1);
    gsl_matrix *lk=gsl_matrix_calloc(3,1);

    gsl_matrix *Z=gsl_matrix_calloc(ns,2);
    gsl_matrix *Zt=gsl_matrix_calloc(2,ns);
    gsl_matrix *Zgamma=gsl_matrix_calloc(ns,ns);
    
    gsl_vector *collectorXs=gsl_vector_calloc(ns);
    gsl_vector *collectorthetaXs=gsl_vector_calloc(3);
    gsl_vector *collectorbs=gsl_vector_calloc(3);
    gsl_vector *collectorLK=gsl_vector_calloc(3);
    gsl_vector *collectorCxs=gsl_vector_calloc(ncoeffs);


    int  filter_factor=*thinning;
    int  filter_n=rep_n/filter_factor;
    int ufilter_factor=filter_factor*10;
    int ufilter_n=rep_n/ufilter_factor;
    
    
    gsl_matrix_set(osu,0,0,*argv0); //Cx
    gsl_matrix_set(osu,1,0,*argv1); //EFL
    gsl_matrix_set(osu,2,0,*argv2); //GC
    
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
            if(v != 0 ) nI++;
            gsl_matrix_set(Jc,i,j,v);
            gsl_matrix_set(Jc,j,i,v);
            k++;
        }
    }

    gsl_matrix *accept=gsl_matrix_calloc(filter_n,4);
    gsl_matrix *collectorS=gsl_matrix_calloc(filter_n,3*ns);
    gsl_matrix *collectorX=gsl_matrix_calloc(filter_n,ns);
    gsl_matrix *collectorbeta=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorthetaX=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorlk=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorCx=gsl_matrix_calloc(filter_n,ncoeffs);
    gsl_matrix *collectorLamb=gsl_matrix_calloc(filter_n,nI);
    gsl_matrix *collectorDD=gsl_matrix_calloc(filter_n,nI);
    
    
    for(i=0;i<ns;i++){
        gsl_matrix_int_set(Lookup,i,0,0);
    }
    //printf("Number of non-zero counts %d \n",nI);
    
    gsl_matrix *I=gsl_matrix_calloc(nI,1);
    gsl_matrix *U=gsl_matrix_calloc(nI,1);
    
    
    dd=gsl_vector_calloc(nI);
    loglambda=gsl_vector_calloc(nI);

    
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
    

    bngetDist(S,ns,nI,Lookup,ijTable);
    /* allocate a cubic bspline workspace (k = 4) */
    Bw = gsl_bspline_alloc(4, nbreak);
    Bx = gsl_matrix_alloc(nI, ncoeffs);
    B = gsl_vector_alloc(ncoeffs);
    Cx = gsl_vector_alloc(ncoeffs);
    Cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    //mw = gsl_multifit_linear_alloc(nI, ncoeffs);
    gsl_vector *Bxi=gsl_vector_alloc(ncoeffs);
    
    MaxKnot  = 0.0;
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        if(v > MaxKnot) MaxKnot=v;
        
    }
    MaxKnot*=1.1;
    
    double ddsort[nI];
    gsl_vector *breaks=gsl_vector_alloc(nbreak);
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        ddsort[i]=v;
    }
    
    gsl_sort (ddsort, 1, nI);
    
    double bk1= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.05);
    double bk2= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.10);
    double bk3= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.15);
    double bk4= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.35);
    double bk5= gsl_stats_median_from_sorted_data (ddsort,1, nI);
    double bk6= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.65);
    double bk7= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.75);
    
    gsl_vector_set(breaks,0,0.0);
    gsl_vector_set(breaks,1,bk1);
    gsl_vector_set(breaks,2,bk2);
    gsl_vector_set(breaks,3,bk3);
    gsl_vector_set(breaks,4,bk4);
    gsl_vector_set(breaks,5,bk5);
    gsl_vector_set(breaks,6,bk6);
    gsl_vector_set(breaks,7,bk7);
    gsl_vector_set(breaks,8,MaxKnot);
    
    gsl_bspline_knots (breaks, Bw);

    /* construct the fit matrix X */
    double xi=0.0;double Bj=0.0;
    for (i = 0; i < nI; ++i){
        xi = gsl_vector_get(dd, i);
        
        /* compute B_j(xi) for all j */
        gsl_bspline_eval(xi, B, Bw);
        
        /* fill in row i of X */
        for (jj = 0; jj < ncoeffs; ++jj)
        {
            Bj = gsl_vector_get(B, jj);
            gsl_matrix_set(Bx, i, jj, Bj);
        }
    }
    
   // Rprintf("ncoeffs %d \n",ncoeffs);
    gsl_vector_set(Cx,0,3.0);gsl_vector_set(Cx,1,2.7);gsl_vector_set(Cx,2,2.5);
    gsl_vector_set(Cx,3,2.0);gsl_vector_set(Cx,4,1.8);gsl_vector_set(Cx,5,1.5);
    gsl_vector_set(Cx,6,1.2);gsl_vector_set(Cx,7,0.8);gsl_vector_set(Cx,8,0.5);
    gsl_vector_set(Cx,9,0.3);gsl_vector_set(Cx,10,0.1);
    
    
    // Rprintf("Cx:");
    for(j=0;j<ncoeffs;j++){
        v=gsl_vector_get(Cx,j);
        //Rprintf("%f ",v);
    }
    
    for(i=0;i<nI;i++){
        v=0.0;
        gsl_matrix_get_row (Bxi, Bx, i);
        gsl_vector_mul (Bxi, Cx);
        for(jj=0;jj<ncoeffs;jj++){
            v+=gsl_vector_get(Bxi,jj);
        }
        gsl_vector_set(loglambda,i,v);
        //printf("%f ",gsl_vector_get(loglambda,i));
    }
    
    
    /* do the fit */
    //gsl_vector_set_all (w, 1.0);
    //gsl_multifit_linear(Bx, loglambda, Cx, Cov, &chisq, mw);

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
    
    gsl_matrix_set(beta,0,0,0.0);gsl_matrix_set(beta,1,0,0.0);gsl_matrix_set(beta,2,0,1.0);
    Rprintf("S %f %f %f\n",gsl_matrix_get(S,0,1),gsl_matrix_get(S,1,0),gsl_matrix_get(S,2,0));
    
    sv= bngetlkI(I,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    
    Rprintf("MAP is set off\n");
    //printf("%f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
    //	printf("INITIALIZATION \n");
    //=====================================

    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
    Rprintf("\n====================================\n");
    int check=0;

        Rprintf("Jumping rules in Proposals\n");
        Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
        Rprintf("Cx : osu %f \n",gsl_matrix_get(osu,0,0));
        Rprintf("EFL : osu %f \n",gsl_matrix_get(osu,1,0));
        Rprintf("GC : osu %f \n",gsl_matrix_get(osu,2,0));
        Rprintf("Looking for jumping rules...\n");
        
        //for(j=0;j<7;j++){
        j=0;
        while(check < 4){
            fi=0;ufi=0;
            check=0;j++;
            if(j==20) break;
            nn=10000;howoften=1;
            //if(j<10) {nn=10000;howoften=1;adjust=0.03;printf("adjust %f\n",adjust);}
            //else {nn =10000;howoften=1;adjust=0.01;}
            for(i=0;i<nn;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                
                if (i%5000==0) {
                    Rprintf("i %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
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
            Rprintf("Cx : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),v1/10000.0);
            Rprintf("EFL : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),v2/10000.0);
            Rprintf("GC : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),v3/10000.0);
            
            
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
            
            Rprintf("%d %d\n",j,check);
            //=================================================================================
            if(check < 4){
                for(i=0;i<ns;i++){
                    gsl_matrix_set(S,i,0,line_t*i/(ns-1));
                    gsl_matrix_set(S,i,1,0.0);
                    gsl_matrix_set(S,i,2,0.0);
                }
                
                bngetDist(S,ns,nI,Lookup,ijTable);
                Bw = gsl_bspline_alloc(4, nbreak);
                Bx = gsl_matrix_alloc(nI, ncoeffs);
                B = gsl_vector_alloc(ncoeffs);
                Cx = gsl_vector_alloc(ncoeffs);
                Cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
                //mw = gsl_multifit_linear_alloc(nI, ncoeffs);
                //gsl_vector *Bxi=gsl_vector_alloc(ncoeffs);
                
                MaxKnot  = 0.0;
                for(i=0;i<nI;i++){
                    v=gsl_vector_get(dd,i);
                    if(v > MaxKnot) MaxKnot=v;
                    
                }
                MaxKnot*=1.1;

                for(i=0;i<nI;i++){
                    v=gsl_vector_get(dd,i);
                    ddsort[i]=v;
                }
                
                gsl_sort (ddsort, 1, nI);
                
                bk1= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.05);
                bk2= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.10);
                bk3= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.15);
                bk4= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.35);
                bk5= gsl_stats_median_from_sorted_data (ddsort,1, nI);
                bk6= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.65);
                bk7= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.75);
                
                gsl_vector_set(breaks,0,0.0);
                gsl_vector_set(breaks,1,bk1);
                gsl_vector_set(breaks,2,bk2);
                gsl_vector_set(breaks,3,bk3);
                gsl_vector_set(breaks,4,bk4);
                gsl_vector_set(breaks,5,bk5);
                gsl_vector_set(breaks,6,bk6);
                gsl_vector_set(breaks,7,bk7);
                gsl_vector_set(breaks,8,MaxKnot);
                
                gsl_bspline_knots (breaks, Bw);
                
                /* construct the fit matrix X */
                for (i = 0; i < nI; ++i){
                    xi = gsl_vector_get(dd, i);
                    
                    /* compute B_j(xi) for all j */
                    gsl_bspline_eval(xi, B, Bw);
                    
                    /* fill in row i of X */
                    for (jj = 0; jj < ncoeffs; ++jj)
                    {
                        Bj = gsl_vector_get(B, jj);
                        gsl_matrix_set(Bx, i, jj, Bj);
                    }
                }

                Rprintf("ncoeffs %d \n",ncoeffs);
                gsl_vector_set(Cx,0,3.0);gsl_vector_set(Cx,1,2.7);gsl_vector_set(Cx,2,2.5);
                gsl_vector_set(Cx,3,2.0);gsl_vector_set(Cx,4,1.8);gsl_vector_set(Cx,5,1.5);
                gsl_vector_set(Cx,6,1.2);gsl_vector_set(Cx,7,0.8);gsl_vector_set(Cx,8,0.5);
                gsl_vector_set(Cx,9,0.3);gsl_vector_set(Cx,10,0.1);
                
                for(i=0;i<nI;i++){
                    v=0.0;
                    gsl_matrix_get_row (Bxi, Bx, i);
                    gsl_vector_mul (Bxi, Cx);
                    for(jj=0;jj<ncoeffs;jj++){
                        v+=gsl_vector_get(Bxi,jj);
                    }
                    gsl_vector_set(loglambda,i,v);
                }
                
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
            
                gsl_matrix_set(beta,0,0,0.0);gsl_matrix_set(beta,1,0,0.0);gsl_matrix_set(beta,2,0,1.0);
                Rprintf("S %f %f %f\n",gsl_matrix_get(S,0,1),gsl_matrix_get(S,1,0),gsl_matrix_get(S,2,0));
            
                sv= bngetlkI(I,X,U,beta,ns,Lookup,ijTable);
                gsl_matrix_set(lk,0,0,sv);
            }
        }
    
            fi=0;ufi=0;
            Rprintf("\nBurning has started\n");
            for(i=0;i<rep_b;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                if ((i+1)==rep_b) {
                    Rprintf("i %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
                }
            }
        
            Rprintf("Burning has ended\n");
        
            Rprintf("\nPosterior sampling has started\n");
            for(i=0;i<rep_n;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
            
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
                
                    gsl_matrix_get_col (collectorbs, beta, 0);
                    gsl_matrix_set_row (collectorbeta, fi, collectorbs);
                
                    gsl_matrix_set_row (collectorCx, fi, Cx);
                    gsl_matrix_set_row (collectorLamb, fi, loglambda);
                    gsl_matrix_set_row (collectorDD, fi, dd);
                    
                    gsl_matrix_get_col (collectorLK, lk, 0);
                    gsl_matrix_set_row (collectorlk, fi, collectorLK);
                
                
                    fi+=1;
                }
                if ((i+1)%ith==0) {
                    Rprintf("At iteration %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
                
                }
            }
    
    time(&stop);
    diff = difftime(stop, start);
    
    Rprintf("==========Summary Report================\n");
    Rprintf("CPU TIME IS %f\n",diff/60.0);
    
    double mean;

        k=0;
        for(i=0;i<filter_n;i++){
            for(j=0;j<ns;j++){
                result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
            }
        }
        for(i=0;i<filter_n;i++){
            result[k]=gsl_matrix_get(collectorbeta,i,0);k++;
            result[k]=gsl_matrix_get(collectorbeta,i,1);k++;
        }
        for(i=0;i<filter_n;i++){
            for(j=0;j<ncoeffs;j++){
                result[k]=gsl_matrix_get(collectorCx,i,j);k++;
            }
        }

        for(j=0;j<4;j++) {
            mean=0.0;
            for(i=0;i<filter_n;i++){
                mean+=gsl_matrix_get(accept,i,j);
            }
            mean/=filter_n;
            result[k]=mean;k++;
            if(j==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
            if(j==1) Rprintf("Cx : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
            if(j==2) Rprintf("EFL : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),mean);
            if(j==3) Rprintf("GC : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),mean);
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
    gsl_matrix_free(tOne);
    gsl_matrix_free(One);
    gsl_matrix_free(accept);
    gsl_matrix_free(collectorS);
    gsl_matrix_free(collectorbeta);
    gsl_matrix_free(collectorthetaX);
    gsl_matrix_free(collectorlk);
    gsl_matrix_free(collectorCx);
    gsl_matrix_free(collectorLamb);
    gsl_matrix_free(collectorDD);
            
    gsl_matrix_free(collectorX);
    gsl_matrix_free(collectorU);
    //gsl_matrix_free(U2R);
    
    gsl_vector_free(collectorbs);
    gsl_vector_free(collectorthetaXs);
    gsl_vector_free(collectorLK);
    gsl_vector_free(collectorXs);
    gsl_vector_free(collectorUs);
    
    
    return;
    
}



void bn2(
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
    //gsl_matrix *U2R = gsl_matrix_calloc(1,nu2r);
    
    int rep_b=*repb;
    int rep_n=*repn;
    
    Rprintf("%d %d %d %d\n",ns,nI,rep_b,rep_n);
    
    int i,j,k,howmany,ii,jj,kk,uk;
    double v,Mi,line_t;
    double max_count;
    
    CEIL=3.0;
    anchor=0.0;
    
    ncoeffs = 11;
    nbreak = ncoeffs-2;
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
    gsl_matrix *beta=gsl_matrix_calloc(3,1);
    gsl_matrix *thetaX=gsl_matrix_calloc(3,1);
    gsl_matrix *osu=gsl_matrix_calloc(7,1);
    gsl_matrix *lk=gsl_matrix_calloc(3,1);
    
    gsl_matrix *Z=gsl_matrix_calloc(ns,2);
    gsl_matrix *Zt=gsl_matrix_calloc(2,ns);
    gsl_matrix *Zgamma=gsl_matrix_calloc(ns,ns);
    
    gsl_vector *collectorXs=gsl_vector_calloc(ns);
    gsl_vector *collectorthetaXs=gsl_vector_calloc(3);
    gsl_vector *collectorbs=gsl_vector_calloc(3);
    gsl_vector *collectorLK=gsl_vector_calloc(3);
    gsl_vector *collectorCxs=gsl_vector_calloc(ncoeffs);
    
    
    int  filter_factor=*thinning;
    int  filter_n=rep_n/filter_factor;
    int ufilter_factor=filter_factor*10;
    int ufilter_n=rep_n/ufilter_factor;
    
    Rprintf(" %d %d %d %d\n",filter_factor,filter_n,ufilter_factor,ufilter_n);
    
    gsl_matrix_set(osu,0,0,*argv0); //Cx
    gsl_matrix_set(osu,1,0,*argv1); //EFL
    gsl_matrix_set(osu,2,0,*argv2); //GC
    

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
            if(v != 0 ) nI++;
            gsl_matrix_set(Jc,i,j,v);
            gsl_matrix_set(Jc,j,i,v);
            k++;
        }
    }
    
    gsl_matrix *accept=gsl_matrix_calloc(filter_n,4);
    gsl_matrix *collectorS=gsl_matrix_calloc(filter_n,3*ns);
    gsl_matrix *collectorX=gsl_matrix_calloc(filter_n,ns);
    gsl_matrix *collectorbeta=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorthetaX=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorlk=gsl_matrix_calloc(filter_n,3);
    gsl_matrix *collectorCx=gsl_matrix_calloc(filter_n,ncoeffs);
    gsl_matrix *collectorLamb=gsl_matrix_calloc(filter_n,nI);
    gsl_matrix *collectorDD=gsl_matrix_calloc(filter_n,nI);
    
    for(i=0;i<ns;i++){
        gsl_matrix_int_set(Lookup,i,0,0);
    }
    Rprintf("Number of non-zero counts %d \n",nI);
    
    gsl_matrix *I=gsl_matrix_calloc(nI,1);
    gsl_matrix *U=gsl_matrix_calloc(nI,1);
    
    
    dd=gsl_vector_calloc(nI);
    loglambda=gsl_vector_calloc(nI);
    
    
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
        v=0.0;
        gsl_matrix_set(efl,i,0,v);gsl_matrix_set(Z,i,0,v);
        gsl_matrix_set(gc,i,0,v);gsl_matrix_set(Z,i,1,v);
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
    
    bngetDist(S,ns,nI,Lookup,ijTable);
    /* allocate a cubic bspline workspace (k = 4) */
    Bw = gsl_bspline_alloc(4, nbreak);
    Bx = gsl_matrix_alloc(nI, ncoeffs);
    B = gsl_vector_alloc(ncoeffs);
    Cx = gsl_vector_alloc(ncoeffs);
    Cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    //mw = gsl_multifit_linear_alloc(nI, ncoeffs);
    gsl_vector *Bxi=gsl_vector_alloc(ncoeffs);
    
    MaxKnot  = 0.0;
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        if(v > MaxKnot) MaxKnot=v;
        
    }
    MaxKnot*=1.1;
    
    double ddsort[nI];
    gsl_vector *breaks=gsl_vector_alloc(nbreak);
    for(i=0;i<nI;i++){
        v=gsl_vector_get(dd,i);
        ddsort[i]=v;
    }
    
    gsl_sort (ddsort, 1, nI);
    
    double bk1= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.05);
    double bk2= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.10);
    double bk3= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.15);
    double bk4= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.35);
    double bk5= gsl_stats_median_from_sorted_data (ddsort,1, nI);
    double bk6= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.65);
    double bk7= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.75);
    
    gsl_vector_set(breaks,0,0.0);
    gsl_vector_set(breaks,1,bk1);
    gsl_vector_set(breaks,2,bk2);
    gsl_vector_set(breaks,3,bk3);
    gsl_vector_set(breaks,4,bk4);
    gsl_vector_set(breaks,5,bk5);
    gsl_vector_set(breaks,6,bk6);
    gsl_vector_set(breaks,7,bk7);
    gsl_vector_set(breaks,8,MaxKnot);
    
    gsl_bspline_knots (breaks, Bw);
    
    /* construct the fit matrix X */
    double xi=0.0;double Bj=0.0;
    for (i = 0; i < nI; ++i){
        xi = gsl_vector_get(dd, i);
        
        /* compute B_j(xi) for all j */
        gsl_bspline_eval(xi, B, Bw);
        
        /* fill in row i of X */
        for (jj = 0; jj < ncoeffs; ++jj)
        {
            Bj = gsl_vector_get(B, jj);
            gsl_matrix_set(Bx, i, jj, Bj);
        }
    }
    
    Rprintf("ncoeffs %d \n",ncoeffs);
    gsl_vector_set(Cx,0,3.0);gsl_vector_set(Cx,1,2.7);gsl_vector_set(Cx,2,2.5);
    gsl_vector_set(Cx,3,2.0);gsl_vector_set(Cx,4,1.8);gsl_vector_set(Cx,5,1.5);
    gsl_vector_set(Cx,6,1.2);gsl_vector_set(Cx,7,0.8);gsl_vector_set(Cx,8,0.5);
    gsl_vector_set(Cx,9,0.3);gsl_vector_set(Cx,10,0.1);
    
    
    Rprintf("Cx:");
    for(jj=0;jj<ncoeffs;jj++){
        v=gsl_vector_get(Cx,jj);
        Rprintf("%f ",v);
    }
    
    for(i=0;i<nI;i++){
        v=0.0;
        gsl_matrix_get_row (Bxi, Bx, i);
        gsl_vector_mul (Bxi, Cx);
        for(jj=0;jj<ncoeffs;jj++){
            v+=gsl_vector_get(Bxi,jj);
        }
        gsl_vector_set(loglambda,i,v);
        //printf("%f ",gsl_vector_get(loglambda,i));
    }
    
    
    /* do the fit */
    //gsl_vector_set_all (w, 1.0);
    //gsl_multifit_linear(Bx, loglambda, Cx, Cov, &chisq, mw);
    
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
    
    gsl_matrix_set(beta,0,0,0.0);gsl_matrix_set(beta,1,0,0.0);gsl_matrix_set(beta,2,0,1.0);
    Rprintf("S %f %f %f\n",gsl_matrix_get(S,0,1),gsl_matrix_get(S,1,0),gsl_matrix_get(S,2,0));
    
    sv= bngetlkI(I,X,U,beta,ns,Lookup,ijTable);
    gsl_matrix_set(lk,0,0,sv);
    
    Rprintf("MAP is set off\n");
    //printf("%f %f %f\n",gsl_matrix_get(lk,0,0),gsl_matrix_get(lk,1,0),gsl_matrix_get(lk,2,0));
    
    //	printf("INITIALIZATION \n");
    //=====================================
    
    int ith=10000;
    if(rep_n<100000) ith=rep_n;
    if(rep_n>=100000 & rep_n<=1000000) ith=100000;
    if(rep_n>1000000) ith=500000;
    
    Rprintf("\n====================================\n");
    int check=0;

        Rprintf("Jumping rules in Proposals\n");
        Rprintf("S : Leap.L Leap.e %d %f \n",HL,epsilon);
        Rprintf("Cx : osu %f \n",gsl_matrix_get(osu,0,0));
        //Rprintf("EFL : osu %f \n",gsl_matrix_get(osu,1,0));
        //Rprintf("GC : osu %f \n",gsl_matrix_get(osu,2,0));
        Rprintf("Looking for jumping rules...\n");
        
        //for(j=0;j<7;j++){
        j=0;
        while(check < 2){
            fi=0;ufi=0;
            check=0;j++;
            if(j==20) break;
            nn=10000;howoften=1;
            //if(j<10) {nn=10000;howoften=1;adjust=0.03;printf("adjust %f\n",adjust);}
            //else {nn =10000;howoften=1;adjust=0.01;}
            for(i=0;i<nn;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                //v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                //v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                
                if (i%5000==0) {
                    Rprintf("i %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    //Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("\n lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
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
            Rprintf("Cx : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),v1/10000.0);
            //Rprintf("EFL : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,1,0),v2/10000.0);
            //Rprintf("GC : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,2,0),v3/10000.0);
            
            
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
            
            
            Rprintf("%d %d\n",j,check);
            //=================================================================================
            if(check < 2){
                for(i=0;i<ns;i++){
                    gsl_matrix_set(S,i,0,line_t*i/(ns-1));
                    gsl_matrix_set(S,i,1,0.0);
                    gsl_matrix_set(S,i,2,0.0);
                }
                
                bngetDist(S,ns,nI,Lookup,ijTable);
                Bw = gsl_bspline_alloc(4, nbreak);
                Bx = gsl_matrix_alloc(nI, ncoeffs);
                B = gsl_vector_alloc(ncoeffs);
                Cx = gsl_vector_alloc(ncoeffs);
                Cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
                //mw = gsl_multifit_linear_alloc(nI, ncoeffs);
                //gsl_vector *Bxi=gsl_vector_alloc(ncoeffs);
                
                MaxKnot  = 0.0;
                for(i=0;i<nI;i++){
                    v=gsl_vector_get(dd,i);
                    if(v > MaxKnot) MaxKnot=v;
                    
                }
                MaxKnot*=1.1;
                
                for(i=0;i<nI;i++){
                    v=gsl_vector_get(dd,i);
                    ddsort[i]=v;
                }
                
                gsl_sort (ddsort, 1, nI);
                
                bk1= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.05);
                bk2= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.10);
                bk3= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.15);
                bk4= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.35);
                bk5= gsl_stats_median_from_sorted_data (ddsort,1, nI);
                bk6= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.65);
                bk7= gsl_stats_quantile_from_sorted_data (ddsort,1, nI,0.75);
                
                gsl_vector_set(breaks,0,0.0);
                gsl_vector_set(breaks,1,bk1);
                gsl_vector_set(breaks,2,bk2);
                gsl_vector_set(breaks,3,bk3);
                gsl_vector_set(breaks,4,bk4);
                gsl_vector_set(breaks,5,bk5);
                gsl_vector_set(breaks,6,bk6);
                gsl_vector_set(breaks,7,bk7);
                gsl_vector_set(breaks,8,MaxKnot);
                
                gsl_bspline_knots (breaks, Bw);
                
                /* construct the fit matrix X */
                for (i = 0; i < nI; ++i){
                    xi = gsl_vector_get(dd, i);
                    
                    /* compute B_j(xi) for all j */
                    gsl_bspline_eval(xi, B, Bw);
                    
                    /* fill in row i of X */
                    for (jj = 0; jj < ncoeffs; ++jj)
                    {
                        Bj = gsl_vector_get(B, jj);
                        gsl_matrix_set(Bx, i, jj, Bj);
                    }
                }
                
                Rprintf("ncoeffs %d \n",ncoeffs);
                gsl_vector_set(Cx,0,3.0);gsl_vector_set(Cx,1,2.7);gsl_vector_set(Cx,2,2.5);
                gsl_vector_set(Cx,3,2.0);gsl_vector_set(Cx,4,1.8);gsl_vector_set(Cx,5,1.5);
                gsl_vector_set(Cx,6,1.2);gsl_vector_set(Cx,7,0.8);gsl_vector_set(Cx,8,0.5);
                gsl_vector_set(Cx,9,0.3);gsl_vector_set(Cx,10,0.1);

                
                for(i=0;i<nI;i++){
                    v=0.0;
                    gsl_matrix_get_row (Bxi, Bx, i);
                    gsl_vector_mul (Bxi, Cx);
                    for(jj=0;jj<ncoeffs;jj++){
                        v+=gsl_vector_get(Bxi,jj);
                    }
                    gsl_vector_set(loglambda,i,v);
                }
                
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
                
                gsl_matrix_set(beta,0,0,0.0);gsl_matrix_set(beta,1,0,0.0);gsl_matrix_set(beta,2,0,1.0);
                Rprintf("S %f %f %f\n",gsl_matrix_get(S,0,1),gsl_matrix_get(S,1,0),gsl_matrix_get(S,2,0));
                
                sv= bngetlkI(I,X,U,beta,ns,Lookup,ijTable);
                gsl_matrix_set(lk,0,0,sv);
            }
            Rprintf("Check %d\n",check);
        }
    
            fi=0;ufi=0;
            Rprintf("\nBurning has started\n");
            for(i=0;i<rep_b;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                //v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                //v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                if ((i+1)==rep_b) {
                    Rprintf("i %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    //Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
                }
            }
            
            Rprintf("Burning has ended\n");
            
            Rprintf("\nPosterior sampling has started\n");
            for(i=0;i<rep_n;i++){
                v0=bnHMCupdateS(S,X,U,I,thetaX,beta,lk,osu,ns,nI,Lookup,ijTable);
                v1=bnupdateCx(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable,S);
                //v2=bnupdatebetaEFL(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                //v3=bnupdatebetaGC(I,X,U,beta,lk,osu,ns,nI,Lookup,ijTable);
                
                if(i%filter_factor==0){
                    gsl_matrix_set(accept,fi,0,v0);
                    gsl_matrix_set(accept,fi,1,v1);
                    //gsl_matrix_set(accept,fi,2,v2);
                    //gsl_matrix_set(accept,fi,3,v3);
                    
                    for(j=0;j<ns;j++){
                        s1=gsl_matrix_get(S,j,0);
                        s2=gsl_matrix_get(S,j,1);
                        s3=gsl_matrix_get(S,j,2);
                        gsl_matrix_set(collectorS,fi,3*j,s1);
                        gsl_matrix_set(collectorS,fi,3*j+1,s2);
                        gsl_matrix_set(collectorS,fi,3*j+2,s3);
                    }
                    
                    //gsl_matrix_get_col (collectorbs, beta, 0);
                    //gsl_matrix_set_row (collectorbeta, fi, collectorbs);
                    
                    //gsl_matrix_get_col (collectorthetaXs, thetaX, 0);
                    //gsl_matrix_set_row (collectorthetaX, fi, collectorthetaXs);
                    
                    gsl_matrix_set_row (collectorCx, fi, Cx);
                    gsl_matrix_set_row (collectorLamb, fi, loglambda);
                    gsl_matrix_set_row (collectorDD, fi, dd);
                    
                    gsl_matrix_get_col (collectorLK, lk, 0);
                    gsl_matrix_set_row (collectorlk, fi, collectorLK);
                    
                    
                    fi+=1;
                }
                if ((i+1)%ith==0) {
                    Rprintf("At iteration %d Cx ",i);
                    for(k=0;k<ncoeffs;k++){
                        v=gsl_vector_get(Cx,k);
                        Rprintf("%f ",v);
                    }
                    //Rprintf("beta (%f %f %f) ",gsl_matrix_get(beta,0,0),gsl_matrix_get(beta,1,0),gsl_matrix_get(beta,2,0));
                    Rprintf("lk %f ",gsl_matrix_get(lk,0,0));
                    Rprintf("MaxKnot %f\n",MaxKnot);
                    
                }
            }

    
    time(&stop);
    diff = difftime(stop, start);
    
    Rprintf("==========Summary Report================\n");
    Rprintf("CPU TIME IS %f\n",diff/60.0);
    
    double mean;

        k=0;
        for(i=0;i<filter_n;i++){
            for(j=0;j<ns;j++){
                result[k]=gsl_matrix_get(collectorS,i,3*j);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+1);k++;
                result[k]=gsl_matrix_get(collectorS,i,3*j+2);k++;
            }
        }
    
        for(i=0;i<filter_n;i++){
            for(j=0;j<ncoeffs;j++){
                result[k]=gsl_matrix_get(collectorCx,i,j);k++;
            }
        }
        
        for(j=0;j<2;j++) {
            mean=0.0;
            for(i=0;i<filter_n;i++){
                mean+=gsl_matrix_get(accept,i,j);
            }
            mean/=filter_n;
            result[k]=mean;k++;
            if(j==0) Rprintf("S : Leap.L Leap.e %d %f  acceptance rate %f\n",HL,epsilon,mean);
            if(j==1) Rprintf("Cx : proposal jump %f acceptance rate %f\n",gsl_matrix_get(osu,0,0),mean);
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
    gsl_matrix_free(tOne);
    gsl_matrix_free(One);
    gsl_matrix_free(accept);
    gsl_matrix_free(collectorS);
    gsl_matrix_free(collectorbeta);
    gsl_matrix_free(collectorthetaX);
    gsl_matrix_free(collectorlk);
    gsl_matrix_free(collectorCx);
    gsl_matrix_free(collectorLamb);
    gsl_matrix_free(collectorDD);
    
    gsl_matrix_free(collectorX);
    gsl_matrix_free(collectorU);
    //gsl_matrix_free(U2R);
    
    gsl_vector_free(collectorbs);
    gsl_vector_free(collectorthetaXs);
    gsl_vector_free(collectorLK);
    gsl_vector_free(collectorXs);
    gsl_vector_free(collectorUs);
    
    
    return;
    
}

