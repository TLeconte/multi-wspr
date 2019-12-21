void postSpot(char *date, char *uttime, double freq, float snr, float dt,
	      float drift, char *call, char *loc, char *pwr);
void postNospot(char *date, char *uttime, double freq);

void initWsprNet(void);
void stopWrprNet(void);
