int stepfunc();
int cuinit(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b);
int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b);
int cufinish(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b);

