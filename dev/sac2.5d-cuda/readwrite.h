#include "iotypes.h"
#include <cstdio>

int createlog(char *logfile);
int appendlog(char *logfile, params p, state s);
int writeconfig(char *name,int n,params p, meta md, real *w);
int readconfig(char *cfgfile, params p, meta md, real *w);
