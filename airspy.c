#include <stdlib.h>
#include <stdio.h>
#include <libairspy/airspy.h>
#include "multi_wspr.h"

static struct airspy_device*   device = NULL;

extern struct receiver_options rx_options;
extern int rx_callback(airspy_transfer_t* transfer);

int startairspy(void)
{
    uint32_t result;

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

    result = airspy_set_sample_type(device, AIRSPY_SAMPLE_FLOAT32_REAL);
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

    result = airspy_set_freq(device, rx_options.realfreq); 
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_set_freq() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }

    /* Sampling run non-stop, for stability and sample are dropped or stored */
    result = airspy_start_rx(device, rx_callback, NULL);
    if( result != AIRSPY_SUCCESS ) {
        printf("airspy_start_rx() failed: %s (%d)\n", airspy_error_name(result), result);
        airspy_close(device);
        airspy_exit();
        return EXIT_FAILURE;
    }
    return 0;
}

void stopairspy(void) 
{
    uint32_t result;

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

}
