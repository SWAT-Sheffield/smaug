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

#define NDIM 2
#define NVECDIM 3
#ifdef USE_SAC
   #define NVAR 13
 #else
   #define NVAR 8
 #endif

#define NDERV 11


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
	int n[NDIM];

        real xmax[NDIM];
	int nt;
        real tmax;

        real boundu[NDIM][NVAR];
        real boundl[NDIM][NVAR];

        real cmax;
        int steeringenabled;
        int finishsteering;     
	real dt;
        real dx[NDIM];

        real gamma;
/*constant used for adiabatic hydrodynamics*/
         #ifdef ADIABHYDRO
            real adiab;
        #endif
        real mu;
        real eta;
        real g[NDIM];

	int sodifon;
        int rkon;
        int moddton;
        int divbon;
        int divbfix;
        int cfgsavefrequency;
        int hyperdifmom; 

        int readini;

        real maxviscoef;
        real chyp;
        real chyp3;      
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

/*         #ifdef USE_SAC

         #else

         #endif*/


        #ifdef USE_SAC
           typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,rhob,energyb,b1b,b2b,b3b} CEV;
         #else
           typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3} CEV;
         #endif


typedef enum dvars {current1,current2,current3,pressuret,pressurek,bdotv,soundspeed,divb,cfast,hdnur,hdnul} DEV;
typedef enum tempvars {tmp1, tmp2, tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9, tmprhol, tmprhor } TEV;

typedef struct Source source;
typedef struct Constants constants;
typedef struct Domain domain;
typedef struct Iome iome;
typedef struct Meta meta;
typedef struct Stateinfo stateinfo;
typedef struct params Params;
#endif

