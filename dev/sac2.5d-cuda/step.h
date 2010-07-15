
int cuinit(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd);
//int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b);
int cufinish(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd);
int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b, float **d_wmod, float **d_dwn1, float **d_wd);

