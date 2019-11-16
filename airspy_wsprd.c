/*
 * FreeBSD License
 * Copyright (c) 2016, Guenael
 * All rights reserved.
 *
 * This file is based on AirSpy project & HackRF project
 *   Copyright 2012 Jared Boone <jared@sharebrained.com>
 *   Copyright 2014-2015 Benjamin Vernoux <bvernoux@airspy.com>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>

#include <libairspy/airspy.h>

#include "airspy_wsprd.h"
#include "wsprd.h"
#include "wsprnet.h"
#include "filter.h"

#define SIGNAL_LENGHT      116
#define SIGNAL_SAMPLE_RATE 625
#define AIRSPY_SAMPLE_RATE 5000000
#define CICDOWNSAMPLE (AIRSPY_SAMPLE_RATE/SIGNAL_SAMPLE_RATE)


/* Global declaration for these structs */
struct decoder_options  dec_options;

static struct receiver_options rx_options;
static struct receiver_state   rx_state;

static struct airspy_device*   device = NULL;
static airspy_read_partid_serialno_t readSerial;

/* Thread stuff for separate decoding */
struct decoder_state {
    pthread_t        thread;

    pthread_cond_t   ready_cond;
    pthread_mutex_t  ready_mutex;
};
struct decoder_state dec;


#define N 4
#define FIRLEN 31

/* Callback for each buffer received */
int rx_callback(airspy_transfer_t* transfer) {
    int16_t *sigIn = (int16_t*) transfer->samples;
    uint32_t sigLenght = transfer->sample_count;
    uint32_t len=0;

    static uint32_t decimationIndex=0,mixerphase=0;

    /* CIC integrator buffers */
    static int64_t  Ix[N],Qx[N];
    
    if (rx_state.decode_flag == true) 
		return 0; 

    //printf("1:rx_callback %d %d\n",rx_state.iqIndex, rx_state.decode_flag);

    for(int32_t i=0; i<sigLenght; i++) {
	int32_t st;

    	/* mixer + 1st CIC integrator */
	switch(mixerphase) {
		case 0 :
   			Ix[0] -= (int64_t)(sigIn[i]-2048);
			break;
		case 1 :
   			Qx[0] -= (int64_t)(sigIn[i]-2048);
			break;
		case 2 :
   			Ix[0] += (int64_t)(sigIn[i]-2048);
			break;
		case 3 :
   			Qx[0] += (int64_t)(sigIn[i]-2048);
			break;
	}
	mixerphase=(mixerphase+1) & 3;

	/* N-1 cic integrators */
	for(st=1;st<N;st++)  {
        	Ix[st] += Ix[st-1]; Qx[st] += Qx[st-1];
	}

        /* Decimation by CICDOWNSAMPLE */
        decimationIndex++;
        if (decimationIndex < CICDOWNSAMPLE) continue;
        decimationIndex = 0;

        if ((rx_state.iqIndex+len) < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE) && rx_state.decode_flag == false) {
            rx_state.iSamples[rx_state.iqIndex+len] = Ix[N-1];
            rx_state.qSamples[rx_state.iqIndex+len] = Qx[N-1];
	    len++;
        } 
    }

    /* Send a signal to the decoding threads */
    pthread_mutex_lock(&dec.ready_mutex);
    rx_state.iqIndex+=len;
    if (rx_state.iqIndex == (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE)) {
                rx_state.decode_flag = true;
    }
    //printf("2:rx_callback %d %d\n",rx_state.iqIndex, rx_state.decode_flag);
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
    return 0;
}


static void *wsprDecoder(void *arg) {

    static uint32_t polyphase=0;

    /* CIC comb filter buffers */    
    static int64_t  Itz[N],Qtz[N];

    /* FIR compensation filter buffers */
    static float    firI[FIRLEN], firQ[FIRLEN];

    /* WSPR decoder use buffers of 45000 samples (hardcoded)
       (120 sec max @ 375sps = 45000 samples)
    */
    static float iSamples[45000]= {0};
    static float qSamples[45000]= {0};
    static uint32_t read_len=0,write_len=0;

    while (!dec_options.exit_flag) {

	uint32_t len;
	bool ct;

        pthread_mutex_lock(&dec.ready_mutex);
        while(read_len>=rx_state.iqIndex && !dec_options.exit_flag)
		pthread_cond_wait(&dec.ready_cond, &dec.ready_mutex);

        if(dec_options.exit_flag) {
        	pthread_mutex_unlock(&dec.ready_mutex);
		pthread_exit(NULL);
	}

	len=rx_state.iqIndex-read_len;
        ct = rx_state.decode_flag ;
        pthread_mutex_unlock(&dec.ready_mutex);

        //printf("1:wsprdecode %d %d\n",len, ct);

        if (len) {
    	 for(int32_t i=0; i<len; i++) {
		uint32_t j,k,st;
   		float Isum,Qsum;
        	int64_t Iy,Ity,Qy,Qty;
		int64_t Iin, Qin;

		Iin=rx_state.iSamples[read_len+i];
		Qin=rx_state.qSamples[read_len+i];

        	/* CIC Comb stages */
       		Iy  = Iin - Itz[0]; Itz[0] = Iin; Ity=Iy;
       		Qy  = Qin - Qtz[0]; Qtz[0] = Qin; Qty=Qy;
		for(st=1;st<N;st++)  {
        		Iy  = Ity - Itz[st]; Itz[st] = Ity; Ity = Iy;
        		Qy  = Qty - Qtz[st]; Qtz[st] = Qty; Qty = Qy;
		}

        	/* FIR compensation filter  + polphase 5/3 downsampler */
        	for (j=0; j<FIRLEN-1; j++) {
                	firI[j] = firI[j+1];
                	firQ[j] = firQ[j+1];
		}
       		firI[FIRLEN-1] = (float)Iy;
       		firQ[FIRLEN-1] = (float)Qy;

		/* 5/3 downsampling */
		polyphase+=3;
		if(polyphase<5) continue ;
		polyphase-=5;

        	Isum=Qsum=0.0;
        	for (j=0, k=polyphase; k<91; j++,k+=3) {
            		Isum += firI[j]*zCoef[k];
            		Qsum += firQ[j]*zCoef[k];
        	}

        	/* Save the result in the buffer */
        	iSamples[write_len]=Isum;
        	qSamples[write_len]=Qsum;
		write_len++;
	 }
	}
	read_len+=len;

        if (!ct) continue;

        printf("2:wsprdecode %d %d\n",read_len,write_len);

        /* Search & decode the signal */
        wspr_decode(iSamples, qSamples, write_len);
	read_len=write_len=0;

        pthread_mutex_lock(&dec.ready_mutex);
        while(rx_state.decode_flag == true && !dec_options.exit_flag)
		pthread_cond_wait(&dec.ready_cond, &dec.ready_mutex);
        //printf("3:wsprdecode %d\n",rx_state.decode_flag);
        pthread_mutex_unlock(&dec.ready_mutex);

        if(dec_options.exit_flag) 
		pthread_exit(NULL);
    }
    
    pthread_exit(NULL);
}


double atofs(char *s) {
    /* standard suffixes */
    char last;
    uint32_t len;
    double suff = 1.0;
    len = strlen(s);
    last = s[len-1];
    s[len-1] = '\0';
    switch (last) {
    case 'g':
    case 'G':
        suff *= 1e3;
    case 'm':
    case 'M':
        suff *= 1e3;
    case 'k':
    case 'K':
        suff *= 1e3;
        suff *= atof(s);
        s[len-1] = last;
        return suff;
    }
    s[len-1] = last;
    return atof(s);
}


int32_t parse_u64(char* s, uint64_t* const value) {
    uint_fast8_t base = 10;
    char* s_end;
    uint64_t u64_value;

    if( strlen(s) > 2 ) {
        if( s[0] == '0' ) {
            if( (s[1] == 'x') || (s[1] == 'X') ) {
                base = 16;
                s += 2;
            } else if( (s[1] == 'b') || (s[1] == 'B') ) {
                base = 2;
                s += 2;
            }
        }
    }

    s_end = s;
    u64_value = strtoull(s, &s_end, base);
    if( (s != s_end) && (*s_end == 0) ) {
        *value = u64_value;
        return AIRSPY_SUCCESS;
    } else {
        return AIRSPY_ERROR_INVALID_PARAM;
    }
}


/* Reset flow control variable & decimation variables */
void initSampleStorage() {
    pthread_mutex_lock(&dec.ready_mutex);
    rx_state.decode_flag = false;
    rx_state.iqIndex=0;
    //printf("initSampleStorage\n");
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
}


/* Default options for the decoder */
void initDecoder_options() {
    dec_options.usehashtable = 1;
    dec_options.npasses = 2;
    dec_options.subtraction = 1;
    dec_options.quickmode = 0;
}


/* Default options for the receiver */
void initrx_options() {
    rx_options.linearitygain= 12;
    rx_options.bias = 0;       // No bias
    rx_options.shift = 0;
    rx_options.serialnumber = 0;
    rx_options.packing = 1;
}


void sigint_callback_handler(int signum) {
    fprintf(stdout, "Caught signal %d\n", signum);
    pthread_mutex_lock(&dec.ready_mutex);
    dec_options.exit_flag = true;
    pthread_cond_broadcast(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
}


void usage(void) {
    fprintf(stderr,
            "airspy_wsprd, a simple WSPR daemon for AirSpy receivers\n\n"
            "Use:\tairspy_wsprd -f frequency -c callsign -g locator [options]\n"
            "\t-f dial frequency [(,k,M) Hz], check http://wsprnet.org/ for freq.\n"
            "\t-c your callsign (12 chars max)\n"
            "\t-g your locator grid (6 chars max)\n"
            "Receiver extra options:\n"
            "\t-l linearity gain [0-21] (default: 12)\n"
            "\t-b set Bias Tee [0-1], (default: 0 disabled)\n"
            "\t-r sampling rate [2.5M, 3M, 6M, 10M], (default: 2.5M)\n"
            "\t-p frequency correction (default: 0)\n"
            "\t-u upconverter (default: 0, example: 125M)\n"
            "\t-s S/N: Open device with specified 64bits serial number\n"
            "\t-k packing: Set packing for samples, \n"
            "\t   1=enabled(12bits packed), 0=disabled(default 16bits not packed)\n"
            "Decoder extra options:\n"
            "\t-H do not use (or update) the hash table\n"
            "\t-Q quick mode, doesn't dig deep for weak signals\n"
            "\t-S single pass mode, no subtraction (same as original wsprd)\n"
            "Example:\n"
            "\tairspy_wsprd -f 144.489M -r 2.5M -c A1XYZ -g AB12cd -l 10 -m 7 -v 7\n");
    exit(1);
}


int main(int argc, char** argv) {
    uint32_t opt;
    uint32_t result;
    uint32_t exit_code = EXIT_SUCCESS;

    initrx_options();
    initDecoder_options();

    /* RX buffer allocation */
    rx_state.iSamples=malloc(sizeof(int64_t)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);
    rx_state.qSamples=malloc(sizeof(int64_t)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);

    /* Stop condition setup */
    dec_options.exit_flag   = false;
    rx_state.decode_flag = true;

    if (argc <= 1)
        usage();

    while ((opt = getopt(argc, argv, "f:c:g:r:l:b:s:p:u:k:H:Q:S")) != -1) {
        switch (opt) {
        case 'f': // Frequency
            rx_options.dialfreq = (uint32_t)atofs(optarg);
            break;
        case 'c': // Callsign
            sprintf(dec_options.rcall, "%.12s", optarg);
            break;
        case 'g': // Locator / Grid
            sprintf(dec_options.rloc, "%.6s", optarg);
            break;
        case 'l': // LNA gain
            rx_options.linearitygain = (uint32_t)atoi(optarg);
            if (rx_options.linearitygain < 0) rx_options.linearitygain = 0;
            if (rx_options.linearitygain > 21 ) rx_options.linearitygain = 21;
            break;
        case 'b': // Bias setting
            rx_options.bias = (uint32_t)atoi(optarg);
            if (rx_options.bias < 0) rx_options.bias = 0;
            if (rx_options.bias > 1) rx_options.bias = 1;
            break;
        case 's': // Serial number
            parse_u64(optarg, &rx_options.serialnumber);
            break;
        case 'p': // Fine frequency correction
            rx_options.shift = (int32_t)atoi(optarg);
            break;
        case 'u': // Upconverter frequency
            rx_options.upconverter = (uint32_t)atofs(optarg);
            break;
        case 'k': // Bit packing
            rx_options.packing = (uint32_t)atoi(optarg);
            if (rx_options.packing < 0) rx_options.packing = 0;
            if (rx_options.packing > 1) rx_options.packing = 1;
            break;
        case 'H': // Decoder option, use a hastable
            dec_options.usehashtable = 0;
            break;
        case 'Q': // Decoder option, faster
            dec_options.quickmode = 1;
            break;
        case 'S': // Decoder option, single pass mode (same as original wsprd)
            dec_options.subtraction = 0;
            dec_options.npasses = 1;
            break;
        default:
            usage();
            break;
        }
    }

    if (rx_options.dialfreq == 0) {
        fprintf(stderr, "Please specify a dial frequency.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rcall[0] == 0) {
        fprintf(stderr, "Please specify your callsign.\n");
    }

    if (dec_options.rloc[0] == 0) {
        fprintf(stderr, "Please specify your locator.\n");
    }

    /* Calcule decimation rate & frequency offset for fs/4 shift */
    rx_options.realfreq = rx_options.dialfreq + rx_options.shift + rx_options.upconverter;

    /* Store the frequency used for the decoder */
    dec_options.freq = rx_options.dialfreq;

    /* If something goes wrong... */
    signal(SIGINT, &sigint_callback_handler);
    signal(SIGILL, &sigint_callback_handler);
    signal(SIGFPE, &sigint_callback_handler);
    signal(SIGSEGV, &sigint_callback_handler);
    signal(SIGTERM, &sigint_callback_handler);
    signal(SIGABRT, &sigint_callback_handler);

    if(dec_options.usehashtable)
	loadHashtable();

    pthread_cond_init(&dec.ready_cond, NULL);
    pthread_mutex_init(&dec.ready_mutex, NULL);
    pthread_create(&dec.thread, NULL, wsprDecoder, NULL);

    result = airspy_init();
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_init() failed: %s (%d)\n", airspy_error_name(result), result);
        return EXIT_FAILURE;
    }

    if( rx_options.serialnumber ) {
        result = airspy_open_sn(&device, rx_options.serialnumber);
        if( result != AIRSPY_SUCCESS ) {
            printf("airspy_open_sn() failed: %s (%d)\n", airspy_error_name(result), result);
            airspy_exit();
            return EXIT_FAILURE;
        }
    } else {
        result = airspy_open(&device);
        if( result != AIRSPY_SUCCESS ) {
            printf("airspy_open() failed: %s (%d)\n", airspy_error_name(result), result);
            airspy_exit();
            return EXIT_FAILURE;
        }
    }

    result = airspy_set_sample_type(device, AIRSPY_SAMPLE_UINT16_REAL);
    if (result != AIRSPY_SUCCESS) {
        printf("airspy_set_sample_type() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    result = airspy_set_samplerate(device, AIRSPY_SAMPLE_RATE/2);
    if (result != AIRSPY_SUCCESS) {
        printf("airspy_set_samplerate() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    if(rx_options.packing) {
        result = airspy_set_packing(device, 1);
        if( result != AIRSPY_SUCCESS ) {
            printf("airspy_set_packing() failed: %s (%d)\n", airspy_error_name(result), result);
            airspy_close(device);
            airspy_exit();
            return EXIT_FAILURE;
        }
    }

    result = airspy_set_rf_bias(device, rx_options.bias);
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_set_rf_bias() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    result = airspy_set_linearity_gain(device, rx_options.linearitygain);
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_set_linearity_gain() failed: %s (%d)\n", airspy_error_name(result), result);
    }

    result = airspy_set_freq(device, rx_options.realfreq + 1500); 
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_set_freq() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    result = airspy_board_partid_serialno_read(device, &readSerial);
    if (result != AIRSPY_SUCCESS) {
        fprintf(stderr, "airspy_board_partid_serialno_read() failed: %s (%d)\n",
                airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    /* reduce if bandwidth */
    airspy_r820t_write(device, 11, 0xE0 | 11);

    /* Sampling run non-stop, for stability and sample are dropped or stored */
    result = airspy_start_rx(device, rx_callback, NULL);
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_start_rx() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    /* Print used parameter */
    time_t rawtime;
    time ( &rawtime );
    struct tm *gtm = gmtime(&rawtime);
    printf("\nStarting airspy-wsprd (%04d-%02d-%02d, %02d:%02dz) -- Version 0.2\n",
           gtm->tm_year + 1900, gtm->tm_mon + 1, gtm->tm_mday, gtm->tm_hour, gtm->tm_min);
    printf("  Callsign     : %s\n", dec_options.rcall);
    printf("  Locator      : %s\n", dec_options.rloc);
    printf("  Dial freq.   : %d Hz\n", rx_options.dialfreq);
    printf("  Real freq.   : %d Hz\n", rx_options.realfreq);
    printf("  Gain         : %d\n", rx_options.linearitygain);
    printf("  Bias         : %s\n", rx_options.bias ? "yes" : "no");
    printf("  Bits packing : %s\n", rx_options.packing ? "yes" : "no");
    printf("  S/N          : 0x%08X%08X\n", readSerial.serial_no[2], readSerial.serial_no[3]);

    initWsprNet();

    /* Main loop : Wait, read, decode */
    while (!dec_options.exit_flag) {
	uint64_t usec,uwait;
	struct timespec tp;

        /* Wait for time Sync on 2 mins */
	clock_gettime(CLOCK_REALTIME, &tp);
        usec  = tp.tv_sec * 1000000 + tp.tv_nsec/1000;
        uwait = 118000000 - usec % 120000000 ;

        usleep(uwait);
        //printf("SYNC! RX started\n");

        /* Start to store the samples */
        initSampleStorage();

        /* Use the Store the date at the begin of the frame */
        rawtime=time(NULL)+2;
	gtm = gmtime(&rawtime);
        sprintf(dec_options.date,"%02d%02d%02d", gtm->tm_year - 100, gtm->tm_mon + 1, gtm->tm_mday);
        sprintf(dec_options.uttime,"%02d%02d", gtm->tm_hour, gtm->tm_min);

        usleep(100000000);
    }

    result = airspy_stop_rx(device);
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_stop_rx() failed: %s (%d)\n", airspy_error_name(result), result);
    }

    if(device != NULL) {
        result = airspy_close(device);
        if( result != AIRSPY_SUCCESS ) {
            printf("airspy_close() failed: %s (%d)\n", airspy_error_name(result), result);
        }
        airspy_exit();
    }

    if(dec_options.usehashtable)
	saveHashtable();

    printf("Bye!\n");

    stopWrprNet();

    /* Wait the thread join (send a signal before to terminate the job) */
    pthread_mutex_lock(&dec.ready_mutex);
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
    pthread_join(dec.thread, NULL);

    /* Destroy the lock/cond/thread */
    pthread_cond_destroy(&dec.ready_cond);
    pthread_mutex_destroy(&dec.ready_mutex);
    pthread_exit(NULL);

    return exit_code;
}
