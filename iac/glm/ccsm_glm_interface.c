void initglm_(int *restart, int *curryear) {
    initialize(restart, curryear);
}
void runglm_() {
    run_model();
}
void stepglm_(int *curryear,double *glmi,int* glmi_fdim1, int* glmi_fdim2, double *glmi_wh, int *glmi_wh_fdim1,double *glmo, int* glmo_fdim1, int* glmo_fdim2) {
  stepglm_ccsm(curryear,glmi,glmi_fdim1,glmi_fdim2,glmi_wh,glmi_wh_fdim1,glmo,glmo_fdim1,glmo_fdim2);
}
void finalizeglm_() {
    finalize();
}
