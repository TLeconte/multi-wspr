/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.
 
 File name: wsprd.c
 
 Copyright 2001-2018, Joe Taylor, K1JT
 
 Much of the present code is based on work by Steven Franke, K9AN,
 which in turn was based on earlier work by K1JT.
 
 Copyright 2014-2018, Steven Franke, K9AN

 Modified for use in multi_wspr 
 Copyright 2019, Thierry Leconte, F4DWV
 
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "fano.h"
#include "jelinek.h"
#include "nhash.h"
#include "wsprd_utils.h"
#include "wsprsim_utils.h"
#include "metric_tables.h"

#include "wsprd.h"
#include "wsprnet.h"

#define max(x,y) ((x) > (y) ? (x) : (y))

#ifdef OSDWSPR
extern void osdwspr_ (float [], unsigned char [], int *, unsigned char [], int *, float *);
#endif

static const unsigned char pr3[162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};

static struct snode *stack;


//***************************************************************************
static void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode)
{
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/
    
    static __thread float fplast=-10000.0;
    const float dt=1.0/375.0, df=375.0/256.0;
    const float twopidt=2*M_PI*dt, df15=df*1.5, df05=df*0.5;

    int i, j, k, lag;
    float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
    float p0,p1,p2,p3,cmet,totp,syncmax,fac;
    float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float f0=0.0, fp, ss, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
    int best_shift = 0, ifreq;

    syncmax=-1e30;
    if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1;}
    if( mode == 1 ) {lagmin=*shift1;lagmax=*shift1;f0=*f1;}
    if( mode == 2 ) {lagmin=*shift1;lagmax=*shift1;ifmin=0;ifmax=0;f0=*f1;}
    
    for(ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
            ss=0.0;
            totp=0.0;
            for (i=0; i<162; i++) {
                fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=twopidt*(fp-df15);
                    cdphi0=cos(dphi0);
                    sdphi0=sin(dphi0);
                    
                    dphi1=twopidt*(fp-df05);
                    cdphi1=cos(dphi1);
                    sdphi1=sin(dphi1);
                    
                    dphi2=twopidt*(fp+df05);
                    cdphi2=cos(dphi2);
                    sdphi2=sin(dphi2);
                    
                    dphi3=twopidt*(fp+df15);
                    cdphi3=cos(dphi3);
                    sdphi3=sin(dphi3);
                    
                    c0[0]=1; s0[0]=0;
                    c1[0]=1; s1[0]=0;
                    c2[0]=1; s2[0]=0;
                    c3[0]=1; s3[0]=0;
                    
                    for (j=1; j<256; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
                    }
                    fplast = fp;
                }
                
                i0[i]=0.0; q0[i]=0.0;
                i1[i]=0.0; q1[i]=0.0;
                i2[i]=0.0; q2[i]=0.0;
                i3[i]=0.0; q3[i]=0.0;
                
                for (j=0; j<256; j++) {
                    k=lag+i*256+j;
                    if( (k>0) && (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];

                p0=sqrt(p0);
                p1=sqrt(p1);
                p2=sqrt(p2);
                p3=sqrt(p3);
                
                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss = (pr3[i] == 1) ? ss+cmet : ss-cmet;
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]==1) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
            }
            ss=ss/totp;
            if( ss > syncmax ) {          //Save best parameters
                syncmax=ss;
                best_shift=lag;
                fbest=f0;
            }
        } // lag loop
    } //freq loop
    
    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }
    
    if( mode == 2 ) {
        *sync=syncmax;
        for (i=0; i<162; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/162.0;
            f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
        }
        fac=sqrt(f2sum-fsum*fsum);
        for (i=0; i<162; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i]=fsymb[i] + 128;
        }
        return;
    }
    return;
}

static void noncoherent_sequence_detection(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1,  int *shift1,
                         float *drift1, int symfac, int *nblocksize)
{
    /************************************************************************
     *  Noncoherent sequence detection for wspr.                            *
     *  Allowed block lengths are nblock=1,2,3,6, or 9 symbols.             * 
     *  Longer block lengths require longer channel coherence time.         *
     *  The whole block is estimated at once.                               *
     *  nblock=1 corresponds to noncoherent detection of individual symbols *
     *     like the original wsprd symbol demodulator.                      *
     ************************************************************************/
    static __thread float fplast=-10000.0;
    const float dt=1.0/375.0, df=375.0/256.0;
    const float twopidt=2*M_PI*dt, df15=df*1.5, df05=df*0.5;
    
    int i, j, k, lag, itone, ib, b, nblock, nseq, imask;
    float xi[512],xq[512];
    float is[4][162],qs[4][162],cf[4][162],sf[4][162],cm,sm,cmp,smp;
    float p[512],fac,xm1,xm0;
    float c0[257],s0[257],c1[257],s1[257],c2[257],s2[257],c3[257],s3[257];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float f0, fp, fsum=0.0, f2sum=0.0, fsymb[162];
    
    f0=*f1;
    lag=*shift1;
    nblock=*nblocksize;
    nseq=1<<nblock;

    for (i=0; i<162; i++) {
        fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
        if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
            dphi0=twopidt*(fp-df15);
            cdphi0=cos(dphi0);
            sdphi0=sin(dphi0);
                    
            dphi1=twopidt*(fp-df05);
            cdphi1=cos(dphi1);
            sdphi1=sin(dphi1);
                    
            dphi2=twopidt*(fp+df05);
            cdphi2=cos(dphi2);
            sdphi2=sin(dphi2);
                    
            dphi3=twopidt*(fp+df15);
            cdphi3=cos(dphi3);
            sdphi3=sin(dphi3);
                    
            c0[0]=1; s0[0]=0;
            c1[0]=1; s1[0]=0;
            c2[0]=1; s2[0]=0;
            c3[0]=1; s3[0]=0;
                    
            for (j=1; j<257; j++) {
                c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
            }

            fplast = fp;
        }

        cf[0][i]=c0[256]; sf[0][i]=s0[256];
        cf[1][i]=c1[256]; sf[1][i]=s1[256];
        cf[2][i]=c2[256]; sf[2][i]=s2[256];
        cf[3][i]=c3[256]; sf[3][i]=s3[256];

        is[0][i]=0.0; qs[0][i]=0.0;
        is[1][i]=0.0; qs[1][i]=0.0;
        is[2][i]=0.0; qs[2][i]=0.0;
        is[3][i]=0.0; qs[3][i]=0.0;
                
        for (j=0; j<256; j++) {
            k=lag+i*256+j;
            if( (k>0) && (k<np) ) {
                is[0][i]=is[0][i] + id[k]*c0[j] + qd[k]*s0[j];
                qs[0][i]=qs[0][i] - id[k]*s0[j] + qd[k]*c0[j];
                is[1][i]=is[1][i] + id[k]*c1[j] + qd[k]*s1[j];
                qs[1][i]=qs[1][i] - id[k]*s1[j] + qd[k]*c1[j];
                is[2][i]=is[2][i] + id[k]*c2[j] + qd[k]*s2[j];
                qs[2][i]=qs[2][i] - id[k]*s2[j] + qd[k]*c2[j];
                is[3][i]=is[3][i] + id[k]*c3[j] + qd[k]*s3[j];
                qs[3][i]=qs[3][i] - id[k]*s3[j] + qd[k]*c3[j];
            }
        }
    }

    for (i=0; i<162; i=i+nblock) {
        for (j=0;j<nseq;j++) {
            xi[j]=0.0; xq[j]=0.0;
            cm=1; sm=0;
            for (ib=0; ib<nblock; ib++) {
                b=(j&(1<<(nblock-1-ib)))>>(nblock-1-ib);
                itone=pr3[i+ib]+2*b;
                xi[j]=xi[j]+is[itone][i+ib]*cm + qs[itone][i+ib]*sm;
                xq[j]=xq[j]+qs[itone][i+ib]*cm - is[itone][i+ib]*sm;
                cmp=cf[itone][i+ib]*cm - sf[itone][i+ib]*sm;
                smp=sf[itone][i+ib]*cm + cf[itone][i+ib]*sm;
                cm=cmp; sm=smp;
            }
            p[j]=xi[j]*xi[j]+xq[j]*xq[j];
            p[j]=sqrt(p[j]);
        }
        for (ib=0; ib<nblock; ib++) {
            imask=1<<(nblock-1-ib);
            xm1=0.0; xm0=0.0;
            for (j=0; j<nseq; j++) {
                if((j & imask)!=0) {
                    if(p[j] > xm1) xm1=p[j];
                }
                if((j & imask)==0) {
                    if(p[j]>xm0) xm0=p[j];
                }
            }
            fsymb[i+ib]=xm1-xm0;
        }
    }
    for (i=0; i<162; i++) {              //Normalize the soft symbols
        fsum=fsum+fsymb[i]/162.0;
        f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
    }
    fac=sqrt(f2sum-fsum*fsum);
    for (i=0; i<162; i++) {
        fsymb[i]=symfac*fsymb[i]/fac;
        if( fsymb[i] > 127) fsymb[i]=127.0;
        if( fsymb[i] < -128 ) fsymb[i]=-128.0;
        symbols[i]=fsymb[i] + 128;
    }
    return;
}


/******************************************************************************
 Fully coherent signal subtraction
 *******************************************************************************/
static void subtract_signal(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
    const float dt=1.0/375.0, df=375.0/256.0, twopidt=2.0*M_PI*dt;
    float phi=0, dphi, cs;
    int i, j, k, ii, nsym=162, nspersym=256,  nfilt=256; //nfilt must be even number.
    int nsig=nsym*nspersym;
    int nc2=45000;
    
    float *refi, *refq, *ci, *cq, *cfi, *cfq;

    refi=calloc(nc2,sizeof(float));
    refq=calloc(nc2,sizeof(float));
    ci=calloc(nc2,sizeof(float));
    cq=calloc(nc2,sizeof(float));
    cfi=calloc(nc2,sizeof(float));
    cfq=calloc(nc2,sizeof(float));
   
    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/
    
    // create reference wspr signal vector, centered on f0.
    //
    for (i=0; i<nsym; i++) {
        
        cs=(float)channel_symbols[i];
        
        dphi=twopidt*
        (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
         );
        
        for ( j=0; j<nspersym; j++ ) {
            ii=nspersym*i+j;
            refi[ii]=cos(phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=sin(phi);
            phi=phi+dphi;
        }
    }
    
    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (i=0; i<nsym*nspersym; i++) {
        k=shift0+i;
        if( (k>0) && (k<np) ) {
            ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }
    
    //lowpass filter and remove startup transient
    float w[nfilt], norm=0, partialsum[nfilt];
    for (i=0; i<nfilt; i++) partialsum[i]=0.0;
    for (i=0; i<nfilt; i++) {
        w[i]=sin(M_PI*(float)i/(float)(nfilt-1));
        norm=norm+w[i];
    }
    for (i=0; i<nfilt; i++) {
        w[i]=w[i]/norm;
    }
    for (i=1; i<nfilt; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }
    
    // LPF
    for (i=nfilt/2; i<45000-nfilt/2; i++) {
        cfi[i]=0.0; cfq[i]=0.0;
        for (j=0; j<nfilt; j++) {
            cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
            cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
        }
    }
    
    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift0.
    for (i=0; i<nsig; i++) {
        if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[nfilt/2+i];
        } else if( i>(nsig-1-nfilt/2) ) {
            norm=partialsum[nfilt/2+nsig-1-i];
        } else {
            norm=1.0;
        }
        k=shift0+i;
        j=i+nfilt;
        if( (k>0) && (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }
    
    free(refi);
    free(refq);
    free(ci);
    free(cq);
    free(cfi);
    free(cfq);

    return;
}

typedef struct {
  fftwf_complex *fftin, *fftout;
  fftwf_plan PLAN;
  hashtelt_t hashtab[32768];
} chndata_t;
static chndata_t chndata[4];

void loadHashtable(uint32_t n, uint32_t fr)
{
        FILE *fhash;
        char line[80];
        char filename[256];

	memset(chndata[n].hashtab,0,sizeof(hashtelt_t)*32768);

	sprintf(filename,"/tmp/hash_%d.txt",fr/1000);

        if((fhash=fopen(filename,"r"))==NULL) return ;

        while (fgets(line, sizeof(line), fhash) != NULL) {
  		int32_t nh;
		time_t t;
        	char hcall[13];
                if(sscanf(line,"%d %ld %12s",&nh,&t,hcall)<3) continue;
                strcpy(chndata[n].hashtab[nh].call,hcall);
                chndata[n].hashtab[nh].t=t;
        }
        fclose(fhash);
}

void saveHashtable(uint32_t n, uint32_t fr) 
{
	FILE *fhash;
        char filename[256];

	sprintf(filename,"/tmp/hash_%d.txt",fr/1000);
        fhash=fopen(filename,"w");
        for (uint32_t i=0; i<32768; i++) {
            if( chndata[n].hashtab[i].t) {
                fprintf(fhash,"%d %ld %s\n",i,chndata[n].hashtab[i].t,chndata[n].hashtab[i].call);
            }
        }
        fclose(fhash);
}

void insHashtable(hashtelt_t *hashtab,char *call)
{
  uint32_t ihash;
  time_t t;

  t=time(NULL);
  ihash=nhash(call,strlen(call),(uint32_t)146);
  strcpy(hashtab[ihash].call,call);
  hashtab[ihash].t=t;
	
}

char *getHashtable(hashtelt_t *hashtab,uint32_t ihash)
{
  time_t t;

  t=time(NULL);
  if(hashtab[ihash].t<t-7200) {
 	// too old
	return NULL;
  }
  return hashtab[ihash].call;
}

// Parameters used for performance-tuning:
static const unsigned int maxcycles=10000;            //Decoder timeout limit
static const float minsync1=0.10;                     //First sync limit
static const int iifac=8;                             //Step size in final DT peakup
static const int symfac=50;                           //Soft-symbol normalizing factor
static const int subtraction=1;
static const int npasses=2;
static const int delta=60;                            //Fano threshold step
static const float bias=0.45;                        //Fano metric bias (used for both Fano and stack algorithms)
static const int more_candidates=1, stackdecoder=0;

static int mettab[2][256];
static float fftwindow[512];

//***************************************************************************
void initwsprd(uint32_t nbc)
{
    uint32_t n,i;

    // setup metric table
    for(i=0; i<256; i++) {
        mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
        mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
    }

    for(i=0; i<512; i++) {
        fftwindow[i]=sin(0.006147931*i);
    }
    for(n=0;n<nbc;n++) {

     chndata[n].fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
     chndata[n].fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);

     // Possible options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT,
     // FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
     chndata[n].PLAN = fftwf_plan_dft_1d(512, chndata[n].fftin, chndata[n].fftout, FFTW_FORWARD, FFTW_EXHAUSTIVE);
  }
}
    
void wspr_decode(float *idat, float *qdat, uint32_t npoints,uint32_t fr,uint32_t chn, char *p_date, char *p_uttime)
{

    const float minrms=52.0 * (symfac/64.0);      //Final test for plausible decoding
    const float dt=1.0/375.0 , df=375.0/256.0/2;

    int i,j,k;
    unsigned char *symbols, *decdata, *channel_symbols, *apmask, *cw;
    signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
    int ipass, nblocksize;
    int maxdrift;
    int shift1, lagmin, lagmax, lagstep, ifmin, ifmax, worth_a_try, not_decoded;
    unsigned int nbits=81, stacksize=200000;
    unsigned int metric, cycles, maxnp;
    float minsync2=0.12;                     //Second sync limit
    float freq0[200],snr0[200],drift0[200],sync0[200];
    int shift0[200];
    double freq_print;
    float f1, fstep, sync1, drift1;
    float psavg[512];
    float dt_print;
    float allfreqs[100];
    char allcalls[100][13];

    symbols=calloc(nbits*2,sizeof(unsigned char));
    apmask=calloc(162,sizeof(unsigned char));
    cw=calloc(162,sizeof(unsigned char));
    decdata=calloc(11,sizeof(unsigned char));
    channel_symbols=calloc(nbits*2,sizeof(unsigned char));
    if( stackdecoder ) {
        stack=calloc(stacksize,sizeof(struct snode));
    }

    for (i=0; i<100; i++) allfreqs[i]=0.0;
    memset(allcalls,0,sizeof(char)*100*13);
    
    int uniques=0, noprint=0, ndecodes_pass=0;
    
    struct result { float sync; float snr; float dt; double freq; 
	 	    char loc[7]; char pwr[4]; char call[13]; float drift;
                    unsigned int cycles; int jitter; int blocksize; unsigned int metric; 
                    };
    struct result decodes[50];
 
    // Do windowed ffts over 2 symbols, stepped by half symbols
    int nffts=4*floor(npoints/512)-1;
    float ps[512][nffts];

    //*************** main loop starts here *****************
    for (ipass=0; ipass<npasses; ipass++) {
        if(ipass == 0) {
            nblocksize=1;
            maxdrift=2;
            minsync2=0.12;
        } else {
           nblocksize=3;  // try all blocksizes up to 3
           maxdrift=1;    // no drift for smaller frequency estimator variance
           minsync2=0.10;
        }
        ndecodes_pass=0;   // still needed?
        
        for (i=0; i<nffts; i++) {
            for(j=0; j<512; j++ ) {
                k=i*128+j;
                chndata[chn].fftin[j][0]=idat[k] * fftwindow[j];
                chndata[chn].fftin[j][1]=qdat[k] * fftwindow[j];
            }
            fftwf_execute(chndata[chn].PLAN);
            for (j=0; j<512; j++ ) {
                k=j+256;
                if( k>511 )
                    k=k-512;
                ps[j][i]=chndata[chn].fftout[k][0]*chndata[chn].fftout[k][0]+chndata[chn].fftout[k][1]*chndata[chn].fftout[k][1];
            }
        }
        
        // Compute average spectrum
        for (i=0; i<512; i++) psavg[i]=0.0;
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                psavg[j]=psavg[j]+ps[j][i];
            }
        }
        
        // Smooth with 7-point window 
        float smspec[506];
        for (i=0; i<506; i++) {
            smspec[i]=0.0;
            for(j=-3; j<=3; j++) {
                k=3+i+j;
                smspec[i]=smspec[i]+psavg[k];
            }
        }
        
        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[506];
        for (j=0; j<506; j++) {
            tmpsort[j]=smspec[j];
        }
        qsort(tmpsort, 506, sizeof(float), floatcomp);
        
        // Noise level of spectrum is estimated as 152/506= 30'th percentile
        float noise_level = tmpsort[152];
        
        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */
        
        float min_snr, snr_scaling_factor;
        min_snr = pow(10.0,-8.0/10.0); //this is min snr in wspr bw
        snr_scaling_factor=26.3;
        for (j=0; j<506; j++) {
            smspec[j]=smspec[j]/noise_level - 1.0;
            if( smspec[j] < min_snr) smspec[j]=0.1*min_snr;
            continue;
        }
        
        // Find all local maxima in smoothed spectrum.
        for (i=0; i<200; i++) {
            freq0[i]=0.0;
            snr0[i]=0.0;
            drift0[i]=0.0;
            shift0[i]=0;
            sync0[i]=0.0;
        }
        
        int npk=0;
        unsigned char candidate;
        if( more_candidates ) {
            for(j=0; j<506; j=j+2) {
                candidate = (smspec[j]>min_snr) && (npk<200);
                if ( candidate ) {
                    freq0[npk]=(j-253)*df;
                    snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
                    npk++;
                }
            }
        } else {
            for(j=1; j<505; j++) {
                candidate = (smspec[j]>smspec[j-1]) &&
                            (smspec[j]>smspec[j+1]) &&
                            (npk<200);
                if ( candidate ) {
                    freq0[npk]=(j-253)*df;
                    snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
                    npk++;
                }
            }
        }

        // bubble sort on snr, bringing freq along for the ride
        int pass;
        float tmp;
        for (pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass ; k++) {
                if (snr0[k] < snr0[k+1]) {
                    tmp = snr0[k];
                    snr0[k] = snr0[k+1];
                    snr0[k+1] = tmp;
                    tmp = freq0[k];
                    freq0[k] = freq0[k+1];
                    freq0[k+1] = tmp;
                }
            }
        }
        
        /* Make coarse estimates of shift (DT), freq, and drift
         
         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
         to nominal start time, which is 2 seconds into the file
         
         * Calculates shift relative to the beginning of the file
         
         * Negative shifts mean that signal started before start of file
         
         * The program prints DT = shift-2 s
         
         * Shifts that cause sync vector to fall off of either end of the data
         vector are accommodated by "partial decoding", such that missing
         symbols produce a soft-decision symbol value of 128
         
         * The frequency drift model is linear, deviation of +/- drift/2 over the
         span of 162 symbols, with deviation equal to 0 at the center of the
         signal vector.
         */

        int idrift,ifr,if0,ifd,k0;
        int kindex;
        float smax,ss,pow,p0,p1,p2,p3;
        for(j=0; j<npk; j++) {                              //For each candidate...
            smax=-1e30;
            if0=freq0[j]/df+256;
            for (ifr=if0-2; ifr<=if0+2; ifr++) {                      //Freq search
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
                        ss=0.0;
                        pow=0.0;
                        for (k=0; k<162; k++) {                             //Sum over symbols
                            ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
                            kindex=k0+2*k;
                            if( kindex >= 0 && kindex < nffts ) {
                                p0=ps[ifd-3][kindex];
                                p1=ps[ifd-1][kindex];
                                p2=ps[ifd+1][kindex];
                                p3=ps[ifd+3][kindex];
                                
                                p0=sqrt(p0);
                                p1=sqrt(p1);
                                p2=sqrt(p2);
                                p3=sqrt(p3);
                                
                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                pow=pow+p0+p1+p2+p3;
                            }
                        }
                        sync1=ss/pow;
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
                            shift0[j]=128*(k0+1);
                            drift0[j]=idrift;
                            freq0[j]=(ifr-256)*df;
                            sync0[j]=sync1;
                        }
                    }
                }
            }
        }

        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].
         
         Function sync_and_demodulate has three modes of operation
         mode is the last argument:
         
         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.
         
         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */
        for (j=0; j<npk; j++) {
            memset(symbols,0,sizeof(char)*nbits*2);

            f1=freq0[j];
            drift1=drift0[j];
            shift1=shift0[j];
            sync1=sync0[j];

            // coarse-grid lag and freq search, then if sync>minsync1 continue
            fstep=0.0; ifmin=0; ifmax=0;
            lagmin=shift1-128;
            lagmax=shift1+128;
            lagstep=64;
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

            fstep=0.25; ifmin=-2; ifmax=2;
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

            if(ipass == 0) {
                // refine drift estimate
                fstep=0.0; ifmin=0; ifmax=0;
                float driftp,driftm,syncp,syncm;
                driftp=drift1+0.5;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &driftp, symfac, &syncp, 1);
            
                driftm=drift1-0.5;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &driftm, symfac, &syncm, 1);
            
                if(syncp>sync1) {
                    drift1=driftp;
                    sync1=syncp;
                } else if (syncm>sync1) {
                    drift1=driftm;
                    sync1=syncm;
                }
            }

            // fine-grid lag and freq search
            if( sync1 > minsync1 ) {
        
                lagmin=shift1-32; lagmax=shift1+32; lagstep=16;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
            
                // fine search over frequency
                fstep=0.05; ifmin=-2; ifmax=2;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

                worth_a_try = 1;
            } else {
                worth_a_try = 0;
            }
            
            int idt, ii, jittered_shift;
            float y,sq,rms;
            not_decoded=1;
            int ib=1, blocksize;

            while( ib <= nblocksize && not_decoded ) {
                blocksize=ib;
                idt=0; ii=0;
                while ( worth_a_try && not_decoded && idt<=(128/iifac)) {
                    ii=(idt+1)/2;
                    if( idt%2 == 1 ) ii=-ii;
                    ii=iifac*ii;
                    jittered_shift=shift1+ii;

                // Use mode 2 to get soft-decision symbols
                    noncoherent_sequence_detection(idat, qdat, npoints, symbols, &f1,
                                    &jittered_shift, &drift1, symfac, &blocksize);
                
                    sq=0.0;
                    for(i=0; i<162; i++) {
                        y=(float)symbols[i] - 128.0;
                        sq += y*y;
                    }
                    rms=sqrt(sq/162.0);

                    if((sync1 > minsync2) && (rms > minrms)) {
                        deinterleave(symbols);
                    
                        if ( stackdecoder ) {
                            not_decoded = jelinek(&metric, &cycles, decdata, symbols, nbits,
                                              stacksize, stack, mettab,maxcycles);
                        } else {
                            not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                                           mettab,delta,maxcycles);
                        }
                    }
                    idt++;
                }
                ib++;
            } 

            if( worth_a_try && !not_decoded ) {
    		char call_loc_pow[23];
	 	char callsign[13];
		char call[13];
		char loc[7];
		char pwr[4];

                ndecodes_pass++;
                
                for(i=0; i<11; i++) {
                    
                    if( decdata[i]>127 ) {
                        message[i]=decdata[i]-256;
                    } else {
                        message[i]=decdata[i];
                    }
                    
                }

                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
		noprint=unpk_(message,chndata[chn].hashtab,callsign,call_loc_pow,call,loc,pwr);

		if(noprint) continue;

                // subtract even on last pass
                if( subtraction && (ipass < npasses ) ) {
                    if( get_wspr_channel_symbols(call_loc_pow, channel_symbols) ) {
                        subtract_signal(idat, qdat, npoints, f1, shift1, drift1, channel_symbols);
                    } else {
                        break;
                    }
                }

                // Remove dupes (same callsign and freq within 3 Hz)
                int dupe=0;
                for (i=0; i<uniques; i++) {
                    if(!strcmp(callsign,allcalls[i]) &&
                       (fabs(f1-allfreqs[i]) <3.0)) dupe=1;
                }
                if( !dupe) {
                    strcpy(allcalls[uniques],callsign);
                    allfreqs[uniques]=f1;
                    uniques++;
                    
                    // Add an extra space at the end of each line so that wspr-x doesn't
                    // truncate the power (TNX to DL8FCL!)
                    
                    freq_print=(fr+f1)/1e6;
                    dt_print=shift1*dt-2;
                    
                    decodes[uniques-1].sync=sync1;
                    decodes[uniques-1].snr=snr0[j];
                    decodes[uniques-1].dt=dt_print;
                    decodes[uniques-1].freq=freq_print;
                    decodes[uniques-1].drift=drift1;
		    strcpy(decodes[uniques-1].call,call);
                    strcpy(decodes[uniques-1].loc,loc);
                    strcpy(decodes[uniques-1].pwr,pwr);
                    decodes[uniques-1].cycles=cycles;
                    decodes[uniques-1].jitter=ii;
                    decodes[uniques-1].blocksize=blocksize;
                    decodes[uniques-1].metric=metric;
                }
            }
        }
    }

    // sort the result in order of increasing frequency
    struct result temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j ; k++) {
            if (decodes[k].freq > decodes[k+1].freq) {
                temp = decodes[k];
                decodes[k]=decodes[k+1];;
                decodes[k+1] = temp;
            }
        }
    }
    
    for (i=0; i<uniques; i++) {
	postSpot(p_date, p_uttime, decodes[i].freq , decodes[i].sync, decodes[i].snr, decodes[i].dt, 
		 decodes[i].call, decodes[i].loc, decodes[i].pwr,decodes[i].drift, decodes[i].cycles,decodes[i].jitter);
    }

   if(uniques ==0 ) 
   	postNospot(p_date, p_uttime, fr/1e6);

    free(apmask);
    free(cw);
    free(decdata);
    free(channel_symbols);
    free(symbols);
    if( stackdecoder ) {
        free(stack);
    }
}

