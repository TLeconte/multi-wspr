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
#include <curl/curl.h>

#include <libairspy/airspy.h>

#include "airspy_wsprd.h"
#include "wsprd.h"
#include "filter.h"

/* TODO
 - BUG : bit packing not working
 - clean/fix serial number section
 - multispot report in one post
 - type fix (uint32_t etc..)
 - verbose option
*/


#define SIGNAL_LENGHT      116
#define SIGNAL_SAMPLE_RATE 375
#define AIRSPY_SAMPLE_RATE 5000000
#define CICDOWNSAMPLE (3*AIRSPY_SAMPLE_RATE/SIGNAL_SAMPLE_RATE/5)

/* Global declaration for these structs */
struct receiver_state   rx_state;
struct receiver_options rx_options;
struct decoder_options  dec_options;
struct decoder_results  dec_results[50];
struct airspy_device*   device = NULL;
airspy_read_partid_serialno_t readSerial;


/* Thread stuff for separate decoding */
struct decoder_state {
    pthread_t        thread;

    pthread_rwlock_t rw;
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

    static uint32_t decimationIndex=0,mixerphase=0,polyphase=0;

    /* CIC buffers */
    static int64_t  Ix[N],Qx[N];
    static int64_t  Itz[N],Qtz[N];

    /* FIR compensation filter buffers */
    static float    firI[FIRLEN], firQ[FIRLEN];
    
    for(int32_t i=0; i<sigLenght; i++) {
	int32_t st,j,k;
    	float Isum,Qsum;
        int64_t  Iy,Ity,Qy,Qty;

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

        /* CIC Comb stages */
       	Iy  = Ix[N-1] - Itz[0]; Itz[0] = Ix[N-1]; Ity=Iy;
       	Qy  = Qx[N-1] - Qtz[0]; Qtz[0] = Qx[N-1]; Qty=Qy;
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

        Isum=0.0, Qsum=0.0;
        for (j=0, k=polyphase; k<91; j++,k+=3) {
            Isum += firI[j]*zCoef[k];
            Qsum += firQ[j]*zCoef[k];
        }

        /* Save the result in the buffer */
        if (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE)) {
            /* Lock the buffer during writing */     // Overkill ?!
            pthread_rwlock_wrlock(&dec.rw);
            rx_state.iSamples[rx_state.iqIndex] = Isum;
            rx_state.qSamples[rx_state.iqIndex] = Qsum;
            pthread_rwlock_unlock(&dec.rw);
            rx_state.iqIndex++;
        } else {
            if (rx_state.decode_flag == false) {
                /* Send a signal to the other thread to start the decoding */
                pthread_mutex_lock(&dec.ready_mutex);
                pthread_cond_signal(&dec.ready_cond);
                pthread_mutex_unlock(&dec.ready_mutex);
                rx_state.decode_flag = true;
                //printf("RX done! [Buffer size: %d]\n", rx_state.iqIndex);
            }
        }
    }
    return 0;
}


void postSpots(uint32_t n_results) {
    CURL *curl;
    CURLcode res;
    char url[1024]; 

    if (n_results == 0) {
        printf("%s No spot\n",dec_options.uttime);
        fflush(stdout);
        return ;
    }

    curl = curl_easy_init();
    if(curl) {
            curl_easy_setopt(curl, CURLOPT_NOBODY, 1);
            curl_easy_setopt(curl, CURLOPT_TIMEOUT, 15);
    } else {
           fprintf(stderr, "curl_easy_init() failed\n");
    }

    for (uint32_t i=0; i<n_results; i++) {

        sprintf(url,"http://wsprnet.org/post?function=wspr&rcall=%s&rgrid=%s&rqrg=%.6f&date=%s&time=%s&sig=%.0f&dt=%.1f&tqrg=%.6f&tcall=%s&tgrid=%s&dbm=%s&version=0.1_wsprd&mode=2",
                dec_options.rcall, dec_options.rloc, dec_results[i].freq, dec_options.date, dec_options.uttime,
                dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
                dec_results[i].call, dec_results[i].loc, dec_results[i].pwr);

        printf("%s %3.2f %4.2f %10.6f %2d  %-s\n",
               dec_options.uttime,dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
               (int)dec_results[i].drift, dec_results[i].message);

        if(curl) {
	    int try=0;
	    do {
		usleep(500000);
		try++;

                curl_easy_setopt(curl, CURLOPT_URL, url);
            	res = curl_easy_perform(curl);

	    } while(res != CURLE_OK && try<3 ); 
            if(res != CURLE_OK)
               	fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
        }
    }
    if(curl) curl_easy_cleanup(curl);
}


static void *wsprDecoder(void *arg) {
    /* WSPR decoder use buffers of 45000 samples (hardcoded)
       (120 sec max @ 375sps = 45000 samples)
    */
    static float iSamples[45000]= {0};
    static float qSamples[45000]= {0};
    static uint32_t samples_len;
    int32_t n_results=0;

    while (!rx_state.exit_flag) {
        pthread_mutex_lock(&dec.ready_mutex);
        pthread_cond_wait(&dec.ready_cond, &dec.ready_mutex);
        pthread_mutex_unlock(&dec.ready_mutex);

        if(rx_state.exit_flag)  // Abord case, final sig
            break;

        /* Lock the buffer access and make a local copy */
        pthread_rwlock_wrlock(&dec.rw);
        memcpy(iSamples, rx_state.iSamples, rx_state.iqIndex * sizeof(float));
        memcpy(qSamples, rx_state.qSamples, rx_state.iqIndex * sizeof(float));
        samples_len = rx_state.iqIndex;  // Overkill ?
        pthread_rwlock_unlock(&dec.rw);

        /* Date and time will be updated/overload during the search & decoding process
           Make a simple copy
        */
        memcpy(dec_options.date, rx_options.date, sizeof(rx_options.date));
        memcpy(dec_options.uttime, rx_options.uttime, sizeof(rx_options.uttime));

        /* DEBUG -- Save samples
        printf("Writing file\n");
        FILE* fd = NULL;
        fd = fopen("samples.bin", "wb");
        int r=fwrite(rx_state.iSamples, sizeof(float), samples_len, fd);
        printf("%d samples written file\n", r);
        fclose(fd);
        */

        /* Search & decode the signal */
        wspr_decode(iSamples, qSamples, samples_len, dec_options, dec_results, &n_results);
        postSpots(n_results);

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
    rx_state.decode_flag = false;
    rx_state.iqIndex=0;
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
    rx_options.packing = 0;
}


void sigint_callback_handler(int signum) {
    fprintf(stdout, "Caught signal %d\n", signum);
    rx_state.exit_flag = true;
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
    rx_state.iSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);
    rx_state.qSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);

    /* Stop condition setup */
    rx_state.exit_flag   = false;
    rx_state.decode_flag = false;

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
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rloc[0] == 0) {
        fprintf(stderr, "Please specify your locator.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
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

    if(rx_options.packing) {
        result = airspy_set_packing(device, 1);
        if( result != AIRSPY_SUCCESS ) {
            printf("airspy_set_packing() failed: %s (%d)\n", airspy_error_name(result), result);
            airspy_close(device);
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

    /* Time alignment stuff */
    struct timeval lTime;
    gettimeofday(&lTime, NULL);
    uint32_t sec   = lTime.tv_sec % 120;
    uint32_t usec  = sec * 1000000 + lTime.tv_usec;
    uint32_t uwait = 120000000 - usec;
    printf("Wait for time sync (start in %d sec)\n\n", uwait/1000000);

    /* Create a thread and stuff for separate decoding
       Info : https://computing.llnl.gov/tutorials/pthreads/
    */
    pthread_rwlock_init(&dec.rw, NULL);
    pthread_cond_init(&dec.ready_cond, NULL);
    pthread_mutex_init(&dec.ready_mutex, NULL);
    pthread_create(&dec.thread, NULL, wsprDecoder, NULL);

    /* Main loop : Wait, read, decode */
    while (!rx_state.exit_flag) {
        /* Wait for time Sync on 2 mins */
        gettimeofday(&lTime, NULL);
        sec   = lTime.tv_sec % 120;
        usec  = sec * 1000000 + lTime.tv_usec;
        uwait = 120000000 - usec + 10000;  // Adding 10ms, to be sure to reach this next minute
        usleep(uwait);
        //printf("SYNC! RX started\n");

        /* Use the Store the date at the begin of the frame */
        time ( &rawtime );
        gtm = gmtime(&rawtime);
        sprintf(rx_options.date,"%02d%02d%02d", gtm->tm_year - 100, gtm->tm_mon + 1, gtm->tm_mday);
        sprintf(rx_options.uttime,"%02d%02d", gtm->tm_hour, gtm->tm_min);

        /* Start to store the samples */
        initSampleStorage();

        while( (airspy_is_streaming(device) == AIRSPY_TRUE) &&
                (rx_state.exit_flag == false) &&
                (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE) ) ) {
            usleep(250000);
        }
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

    /* Wait the thread join (send a signal before to terminate the job) */
    pthread_mutex_lock(&dec.ready_mutex);
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
    pthread_join(dec.thread, NULL);

    /* Destroy the lock/cond/thread */
    pthread_rwlock_destroy(&dec.rw);
    pthread_cond_destroy(&dec.ready_cond);
    pthread_mutex_destroy(&dec.ready_mutex);
    pthread_exit(NULL);

    return exit_code;
}
