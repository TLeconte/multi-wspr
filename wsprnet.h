void postSpot(char *date, char *uttime, double freq, float sync, float snr, float dt, float drift, int32_t jitter,char *message, char *call, char *loc, char *pwr);
void postNospot(char *date, char *uttime, double freq, char *call, char *loc, char *pwr);

void initWsprNet(void);
void stopWrprNet(void);

