#define GLMONFLDS 9
#define PLONFLDS 23

/*  this declaration matches the one in updateannuallanduse_v2.c  */
void updateannuallanduse_main(double glmo[][GLMONFLDS], double plodata[][PLONFLDS], int *inyear,
			int *crop_addtreeonly, double *crop_setherbfracrem, double *crop_setavailtreefracrem,
                int *pasture_addtreeonly, double *pasture_setherbfracrem, double *pasture_setavailtreefracrem);

void updateannuallanduse_(double *glmo[][GLMONFLDS], double *plodata[][PLONFLDS], int *myear,
		int *crop_addtreeonly, double *crop_setherbfracrem, double *crop_setavailtreefracrem,
		int *pasture_addtreeonly, double *pasture_setherbfracrem, double *pasture_setavailtreefracrem) {
    updateannuallanduse_main(glmo, plodata, myear, crop_addtreeonly, crop_setherbfracrem, crop_setavailtreefracrem,
			pasture_addtreeonly, pasture_setherbfracrem, pasture_setavailtreefracrem);
}
