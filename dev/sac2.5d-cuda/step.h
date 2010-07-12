int stepfunc();
int cuinit(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h);
int cuprop(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h);
int cufinish(struct params **p, float **u, float **v, float **b, float **h,struct params **d_p, float **d_u, float **d_v, float **d_b, float **d_h);

