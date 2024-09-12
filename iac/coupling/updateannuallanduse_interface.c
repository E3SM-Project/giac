void updateannuallanduse_(double *glmo[][9], double *plodata[][23], int *myear,
		int *crop_addtreeonly, double *crop_setherbfracrem, double *crop_setavailtreefracrem,
		int *pasture_addtreeonly, double *pasture_setherbfracrem, double *pasture_setavailtreefracrem) {
    updateannuallanduse_main(glmo, plodata, myear, crop_addtreeonly, crop_setherbfracrem, crop_setavailtreefracrem,
			pasture_addtreeonly, pasture_setherbfracrem, pasture_setavailtreefracrem);
}
