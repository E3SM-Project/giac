#define GLMONFLDS 9
#define PLONFLDS 23

/*  this declaration matches the one in updateannuallanduse_v2.c  */
/*  since C is row major (last index varies fastest) the glmo/plodata pointers do not need to know the size
 *  of the first (row) dimension here, as long at the code knows it for proper indexing. they need to know only
 *  the size of the second dimension to know how long each row is */
void updateannuallanduse_main(double glmo[][GLMONFLDS], double plodata[][PLONFLDS], int *inyear,
			int *crop_addtreeonly, double *crop_setherbfracrem, double *crop_setavailtreefracrem,
                int *pasture_addtreeonly, double *pasture_setherbfracrem, double *pasture_setavailtreefracrem);

/*  fortran passes a flat pointer to contiguous memory. it is column major (first index varies fastest) so
 *  the fields vary fastest on disk because they are the first dim in the fortran code. the c code has fixed
 *  2d arrays for glmo and plodata, with a pointer to the contiguous memory. this syntax should be more clear
 *  and robust */
void updateannuallanduse_(double *glmo, double *plodata, int *myear,
		int *crop_addtreeonly, double *crop_setherbfracrem, double *crop_setavailtreefracrem,
		int *pasture_addtreeonly, double *pasture_setherbfracrem, double *pasture_setavailtreefracrem) {
    updateannuallanduse_main((double (*)[GLMONFLDS]) glmo, (double (*)[PLONFLDS]) plodata, myear, crop_addtreeonly,
		    crop_setherbfracrem, crop_setavailtreefracrem, pasture_addtreeonly, pasture_setherbfracrem,
		    pasture_setavailtreefracrem);
}
