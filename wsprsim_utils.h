#ifndef WSPRSIM_UTILS_H
#define WSPRSIM_UTILS_H

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

int get_wspr_channel_symbols(char* message, unsigned char* symbols);

#endif
