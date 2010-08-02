



int cuinit(struct params **p, float **w, float **wnew,   struct state **state, struct params **d_p, float **d_w, float **d_wnew, float **d_wmod, float **d_dwn1, float **d_wd, struct state **d_state);
//int cuprop(struct params **p, float **w, float **wnew, float **b,struct params **d_p, float **d_w, float **d_wnew, float **d_b);
int cufinish(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd);
int cuprop(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd);
int cuboundary(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd);
int cuupdate(struct params **p, float **w, float **wnew,  struct state **state,struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd, struct state **d_state);
int cuderivcurrent(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd, int order);
int cuderivsource(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd);
int cuadvance(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd);
int cupredictor(struct params **p, float **w, float **wnew, struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd, int order);
int cudivb(struct params **p, float **w, float **wnew,  struct state **state,struct params **d_p, float **d_w, float **d_wnew,  float **d_wmod, float **d_dwn1, float **d_wd, struct state **d_state);

