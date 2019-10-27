# airspy-wsprd -- WSPR daemon for AirSpy receivers

This is a modified version of https://github.com/Guenael/airspy-wsprd : a non-interactive application allows automatic reporting of WSPR spots on WSPRnet.  

List of mods :
	- Use linearity gain instead of Lna/mixer/vga gains
	- increase order of cic filter to 4
	- Use a polyphase FIR filter for cic correction in order to have integer downsampling factor
	- limit input sampling rate to 2.5Ms/s (to ease polyphase filer implementation)
	- combine Fs/4 mixer with 1st cic filter integrator
	- Load/save hashtable only once and don't load/save fftwisdom data.
 
