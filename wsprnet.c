/*
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

 wsprnet.c
 Copyright (c) 2019-2020, Thierry Leconte, F4DWV

 extracted from github/Guenael/airspy-wsprd/airspy_wsprd.c
 Copyright (c) 2016, Guenael

*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <curl/curl.h>
#include "multi_wspr.h"

typedef struct spot_s spot_t;
struct spot_s {
	spot_t *next;
	char *date;
	char *uttime;
	double freq;
	float sync;
	float snr;
	float dt;
	int metric;
	char *call;
	char *loc;
	char *pwr;
	float drift;
	int cycles;
	int jitter;
};

static spot_t *spot_head = NULL;
static spot_t *spot_tail = NULL;
static pthread_cond_t spot_cond;
static pthread_mutex_t spot_mutex;
static pthread_t spot_thread;

void postSpot(char *date, char *uttime, double freq, float sync , float snr, float dt,
              char *call, char *loc, char *pwr, float drift, int cycles, int jitter)
{
	spot_t *spot;

	spot = malloc(sizeof(spot_t));

	spot->next = NULL;
	spot->date = strdup(date);
	spot->uttime = strdup(uttime);
	spot->freq = freq;
	spot->sync = sync;
	spot->snr = snr;
	spot->dt = dt;
	spot->call = strdup(call);
	spot->loc = strdup(loc);
	spot->pwr = strdup(pwr);
	spot->drift = drift;
	spot->cycles = cycles;
	spot->jitter = jitter;

	printf("%s %s %2d %+3.2f %+4.2f %10.6f %s %s %s %2d %2d %2d\n",
	       spot->date,spot->uttime, (int)(10*(spot->sync)), spot->snr, spot->dt, spot->freq,
	       spot->call, spot->loc, spot->pwr,(int)(spot->drift),spot->cycles/81,spot->jitter);
	fflush(stdout);

	/* put in list at tail */
	pthread_mutex_lock(&spot_mutex);
	if (spot_tail)
		spot_tail->next = spot;
	else
		spot_head = spot;
	spot_tail = spot;
	pthread_cond_signal(&spot_cond);
	pthread_mutex_unlock(&spot_mutex);

}

void postNospot(char *date, char *uttime, double freq)
{
	printf("%s %10.6f NoSpot\n", uttime, freq);
	fflush(stdout);
}

static void *sendSpots(void *args)
{
	CURL *curl;
	CURLcode res;
	char url[2048];
	spot_t *spot;

	curl = NULL;
	while (1) {
		pthread_mutex_lock(&spot_mutex);
		while (spot_head == NULL) {
			if (curl) {
				curl_easy_cleanup(curl);
				curl = NULL;
			}
			pthread_cond_wait(&spot_cond, &spot_mutex);
		}

		spot = spot_head;
		spot_head = spot_head->next;
		if (spot_head == NULL)
			spot_tail = NULL;
		pthread_mutex_unlock(&spot_mutex);

		url[0] = 0;
		if (dec_options.wsprnet)
			sprintf(url,
				"http://wsprnet.org/post?function=wspr&rcall=%s&rgrid=%s&rqrg=%.6f&date=%s&time=%s&sig=%.0f&dt=%.1f&tqrg=%.6f&tcall=%s&tgrid=%s&dbm=%s&version=1.0_multi_wspr&mode=2",
				dec_options.rcall, dec_options.rloc, spot->freq,
				spot->date, spot->uttime, spot->snr, spot->dt,
				spot->freq, spot->call, spot->loc, spot->pwr);

		free(spot->date);
		free(spot->uttime);
		free(spot->call);
		free(spot->loc);
		free(spot->pwr);
		free(spot);

		if (url[0])
			do {
				if (curl == NULL) {
					usleep(7000000);
					curl = curl_easy_init();
					if (curl == NULL) {
						fprintf(stderr, "curl_easy_init() failed\n");
						break;
					}
					curl_easy_setopt(curl, CURLOPT_NOBODY,
							 1);
					curl_easy_setopt(curl, CURLOPT_TIMEOUT,
							 30);
				}

				curl_easy_setopt(curl, CURLOPT_URL, url);
				res = curl_easy_perform(curl);
				if (res != CURLE_OK) {
					fprintf(stderr, "curl_easy_perform() failed: %s\n",
						curl_easy_strerror(res));
					curl_easy_cleanup(curl);
					curl = NULL;
					if (res == CURLE_OPERATION_TIMEDOUT)
						usleep(5000000);
					else
						usleep(20000000);
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
