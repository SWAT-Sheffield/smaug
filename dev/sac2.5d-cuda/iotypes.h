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

        real xmax;
        real ymax;
	int nt;
        real tmax;

        real *boundxu;
        real *boundxl;
        real *boundyu;
        real *boundyl;

        real cmax;
        int steeringenabled;
        int finishsteering;     
	real dt;
        real dx;
        real dy;
        real g;
        real gamma;
/*constant used for adiabatic hydrodynamics*/
         #ifdef ADIABHYDRO
            real adiab;
        #endif
        real mu;
        real eta;
        real g1;
        real g2;
        real g3;
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
	real t;
	real dt;
	real rho;
        real m1;
	real m2;
	real m3;
	real e;
        real b1;
        real b2;
        real b3;
};

struct hydrovars{
    int numvars; //variables each vector component
	int num;   //total number of dimensions including any ghost variables
	real *w;

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

