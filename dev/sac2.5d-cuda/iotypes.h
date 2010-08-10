#ifndef TYPES_H_
#define TYPES_H_

#define DEFINE_PRECISION(T) \
  typedef T real;


#ifdef USE_REAL
DEFINE_PRECISION(float)
#else
DEFINE_PRECISION(double)
#endif




#undef DEFINE_PRECISION



struct Meta {
   char *directory ;
   char *author;
   char *sdate;
   char *platform;
   char *desc;
   char *name;
   char *ini_file;
   char *log_file;
   char *out_file;
};

struct Iome {
    char *server;
    int port;
    int id;
};

struct params {
	int ni;
 	int nj;

        float xmax;
        float ymax;
	int nt;
        float tmax;

        double *boundxu;
        double *boundxl;
        double *boundyu;
        double *boundyl;

        float cmax;
        int steeringenabled;
        int finishsteering;     
	float dt;
        float dx;
        float dy;
        float g;
        float gamma;
/*constant used for adiabatic hydrodynamics*/
         #ifdef ADIABHYDRO
            float adiab;
        #endif
        float mu;
        float eta;
        float g1;
        float g2;
        float g3;
	int sodifon;
        int rkon;
        int moddton;
        int divbon;
        int divbfix;
        int cfgsavefrequency; 

        int readini;        
};

//it   t   dt    rho m1 m2 e bx by
struct state{
	int it;
	float t;
	float dt;
	float rho;
        float m1;
	float m2;
	float m3;
	float e;
        float b1;
        float b2;
        float b3;
};

struct hydrovars{
    int numvars; //variables each vector component
	int num;   //total number of dimensions including any ghost variables
	float *w;

};



typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3} CEV;
typedef enum dvars {current1,current2,current3,pressuret,pressurek,bdotv,soundspeed,divb} DEV;

typedef struct Source source;
typedef struct Constants constants;
typedef struct Domain domain;
typedef struct Iome iome;
typedef struct Meta meta;
typedef struct Stateinfo stateinfo;

#endif

