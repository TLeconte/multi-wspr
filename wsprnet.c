#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <curl/curl.h>
#include "wsprd.h"

extern struct receiver_options rx_options;

typedef struct spot_s spot_t;
struct spot_s {
	spot_t *next;
	char *date;
	char *uttime;
	double freq;
	float sync;
	float snr;
	float dt;
	float drift;
	int32_t jitter;
	char *message;
	char *call;
	char *loc;
	char *pwr;
};

static spot_t *spot_head=NULL;
static spot_t *spot_tail=NULL;
static pthread_cond_t   spot_cond;
static pthread_mutex_t  spot_mutex;
static pthread_t  spot_thread;


void postSpot(char *date, char *uttime, double freq, float sync, float snr, float dt, float drift, int32_t jitter,char *message, char *call, char *loc, char *pwr)
{
  spot_t *spot;

  spot=malloc(sizeof(spot_t));

  spot->next=NULL;
  spot->date=strdup(date);
  spot->uttime=strdup(uttime);
  spot->freq=freq;
  spot->sync=sync;
  spot->snr=snr;
  spot->dt=dt;
  spot->drift=drift;
  spot->jitter=jitter;
  spot->message=strdup(message);
  spot->call=strdup(call);
  spot->loc=strdup(loc);
  spot->pwr=strdup(pwr);

  /* put in list at tail */
  pthread_mutex_lock(&spot_mutex);
  if(spot_tail) 
	spot_tail->next=spot;
   else
	spot_head=spot;
  spot_tail=spot;
  pthread_cond_signal(&spot_cond);
  pthread_mutex_unlock(&spot_mutex);

}


static void *sendSpots(void * args) {
    CURL *curl;
    CURLcode res;
    char url[1024];
    spot_t *spot;


   curl=NULL;
   while(!dec_options.exit_flag) {
    pthread_mutex_lock(&spot_mutex);
    while (spot_head == NULL)  {
        if(dec_options.exit_flag) return 0;
	if(curl) { 
    		curl_easy_cleanup(curl); curl=NULL;
	}
	pthread_cond_wait(&spot_cond, &spot_mutex);
    }


    spot=spot_head;
    spot_head=spot_head->next;
    if(spot_head == NULL) spot_tail=NULL;

    pthread_mutex_unlock(&spot_mutex);

    printf("%s %3.2f %4.2f %10.6f %2d  %-s\n",
               spot->uttime,spot->snr,spot->dt,spot->freq,(int)spot->drift,spot->message);

    url[0]=0;
    if(dec_options.rcall[0] && dec_options.rloc[0])
	sprintf(url,"http://wsprnet.org/post?function=wspr&rcall=%s&rgrid=%s&rqrg=%.6f&date=%s&time=%s&sig=%.0f&dt=%.1f&tqrg=%.6f&tcall=%s&tgrid=%s&dbm=%s&version=0.1_wsprd&mode=2",
                dec_options.rcall, dec_options.rloc, spot->freq, spot->date, spot->uttime,
                spot->snr, spot->dt, spot->freq,
                spot->call, spot->loc, spot->pwr);

    free(spot->date);
    free(spot->uttime);
    free(spot->message);
    free(spot->call);
    free(spot->loc);
    free(spot->pwr);
    free(spot); 

    if(url[0]) 
      do {
    	if(curl==NULL)  {
    		curl = curl_easy_init();
    		if(curl==NULL)  {
           		fprintf(stderr, "curl_easy_init() failed\n");
	   		break;
    		}
        	curl_easy_setopt(curl, CURLOPT_NOBODY, 1);
        	curl_easy_setopt(curl, CURLOPT_TIMEOUT, 20);
      }

    curl_easy_setopt(curl, CURLOPT_URL, url);
    res = curl_easy_perform(curl);
    if(res != CURLE_OK) {
        fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));
    	curl_easy_cleanup(curl); curl=NULL;
	usleep(500000);
    }
    } while (res != CURLE_OK);

   }
   return 0;
}

void initWsprNet(void)
{
    pthread_cond_init(&spot_cond, NULL);
    pthread_mutex_init(&spot_mutex, NULL);
    pthread_create(&spot_thread, NULL, sendSpots, NULL);
}

void stopWrprNet(void)
{

  pthread_mutex_lock(&spot_mutex);
  pthread_cond_signal(&spot_cond);
  pthread_mutex_unlock(&spot_mutex);

}

