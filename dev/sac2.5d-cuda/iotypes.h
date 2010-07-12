#ifndef TYPES_H_
#define TYPES_H_

struct Constants {
       float g;
       float u0; 
       float v0;
       float b0;
       float h0;
     };


struct Domain {
	int ni;
        float xmax;
 	int nj;
        float ymax;
	int nt;
        float tmax;
        float step;
        int steeringenabled;
        int finishsteering;
};

struct Source {
  float freq;
  float amp;
  float xloc;
  float yloc;
};

struct Meta {
   char *directory ;
   char *author;
   char *sdate;
   char *platform;
   char *desc;
   char *name;
};

struct Iome {
    char *server;
    int port;
    int id;
};

struct params {
	int ni;
 	int nj;     
	float dt;
        float dx;
        float dy;
        float g;       
};


typedef struct Source source;
typedef struct Constants constants;
typedef struct Domain domain;
typedef struct Iome iome;
typedef struct Meta meta;

#endif

