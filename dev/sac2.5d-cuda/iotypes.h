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
   //#define NVAR 13
   #define NVAR 10
   //#define NDERV 19
   #define NDERV 16
   #define NTEMP 11
   #define NDIM 2
   #define NVECDIM 2
#endif
 #ifdef USE_VAC
   //#define NVAR 8
   #define NVAR 6
   //#define NDERV 17
   #define NDERV 14
   #define NTEMP 11
   #define NDIM 2
   #define NVECDIM 2
 #endif
 #ifdef ADIABHYDRO
   //#define NVAR 4
   #define NVAR 3
   //#define NDERV 9
   #define NDERV 6
   #define NTEMP 1
   #define NDIM 2
   #define NVECDIM 3
 #endif

#define PI 3.14159265358979
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

        int mnthreads;   
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

typedef enum oldvars {mom3, b3,b3b} CEVOLD;
#ifdef USE_SAC
   //typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,rhob,energyb,b1b,b2b,b3b} CEV;
   typedef enum vars {rho, mom1, mom2, energy, b1, b2,rhob,energyb,b1b,b2b} CEV;
#else
   //typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3} CEV;
   typedef enum vars {rho, mom1, mom2, energy, b1, b2} CEV;
#endif

#ifdef USE_SAC
//	typedef enum dvars {vel1,vel2,vel3,f1,f2,f3,soundspeed,pressuret,pressurek,current1,current2,current3,bdotv,divb,cfast,hdnur,hdnul,ptb,pkb} DEV;

	typedef enum dvars {vel1,vel2,f1,f2,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul,ptb,pkb} DEV;
#else
//	typedef enum dvars {vel1,vel2,vel3,f1,f2,f3,soundspeed,pressuret,pressurek,current1,current2,current3,bdotv,divb,cfast,hdnur,hdnul} DEV;
	typedef enum dvars {vel1,vel2,f1,f2,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul} DEV;
#endif


typedef enum tempvars {tmp1, tmp2, tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9, tmprhol, tmprhor } TEV;



typedef struct Source source;
typedef struct Constants constants;
typedef struct Domain domain;
typedef struct Iome iome;
typedef struct Meta meta;
typedef struct Stateinfo stateinfo;
typedef struct params Params;
#endif

