# MULTI_WSPR
multi_wspr, is a Multi-simultaneous frequency WSPR receiver for AirSpy SDR

## Features :

 * up to 4 bands received simultaneously

## Usage
> multi_wspr -f frequency_set -c callsign -g locator [options]

 -f frequency set :
		0 : 137.5Khz
		1 : 1.8381Mhz, 3.5701Mhz, 5.3647Mhz, 7.0401Mhz
		2 : 3.5701Mhz, 5.3647Mhz, 7.0401Mhz, 10.1402Mhz
		3 : 7.0401Mhz, 10.1402Mhz, 14.0971Mhz
		4 : 10.1402Mhz, 14.0971Mhz, 18.1061Mhz
		5 : 14.0971Mhz, 18.1061Mhz, 21.0961Mhz
		6 : 18.1061Mhz, 21.0961Mhz, 24.9261Mhz
		7 : 21.0961Mhz, 24.9261Mhz, 28.1261Mhz

  -c your callsign (12 chars max)
  -g your locator grid (6 chars max)

Receiver extra options:

  -l linearity gain [0-21] (default: 12)
  -b set Bias Tee [0-1], (default: 0 disabled)
  -p frequency correction (default: 0)
  -u upconverter (default: 0, example: 120M for SpyVerter)
  -k packing: Set packing for samples, 
	   1=enabled(12bits packed), 0=disabled(default 16bits not packed)

## Example

> multi_wspr -f 1 -u 120M -b 1 -c F4DWV -g IN98bc -l 18 

## Compilation
multi_wspr  must compile directly on any modern Linux distrib by doing :

> make

It depends on some external libraries :
 * libusb
 * libairspy
 * libfftw3f
 * libcurl

## Copyrights 
multi_wspr is Copyright Thierry Leconte 2019

This started as a fork of https://github.com/Guenael/airspy-wsprd, but the frontend have been rewritten to receive multi-simultaneous frequencies, and the decoder part upgraded to the last code from  https://sourceforge.net/p/wsjt/wsjtx/ci/master/tree/lib/wsprd/

See the sources for the others copyrights 


 
