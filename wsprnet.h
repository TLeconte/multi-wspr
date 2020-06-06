void postSpot(char *date, char *uttime, double freq, float sync , float snr, float dt,
	      char *call, char *loc, char *pwr, float drift, int cycles, int jitter);
void postNospot(char *date, char *uttime, double freq);

void initWsprNet(void);
void stopWrprNet(void);
