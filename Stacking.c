#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592

/*Evaluate the volume of a sphere of radius r*/
float Vol(float r){
	float resp;
	if(r<=0.0) resp = 0.0;
	else resp = 4.0*PI*r*r*r/3.0;
	return resp;
}

int main(int argc,char *argv[])
{
FILE *profiles, *out;
char profilesfile[100], outfile[100];
int np, nh, i, j, k, n, n_temp, bins, f, label, cont, j_temp, bin, h_vel, Do_log, flag;
float n_m, radius, len_b, len, r_min, r_max, R_min, R_max, r, log_min, r_times, r_times_v, v, pos_n, log_b, r_mean;
float **density, **velocity;
double *den,  *vel, *err, *errv;

if (argc != 7){
	printf("Wrong number of arguments.\n");
	printf("arg1: Name of the void/halos's profile file\n");
	printf("arg2: Name of the output file\n");
	printf("arg3: Smaller void/halo's radius in stack\n");
	printf("arg4: Bigest void/halo's radius in stack\n");
	printf("arg5: Use the voids/halos that grows on others? Yes (1) or No (0)\n");
	printf("arg6: The particle catalog also have velocities? Yes (1) or No (0).\n\n");
	exit(0);
}

sprintf(profilesfile,"%s", argv[1]);
sprintf(outfile,"%s", argv[2]);
R_min = atof(argv[3]);
R_max = atof(argv[4]);
flag = atoi(argv[5]);
h_vel = atoi(argv[6]);	/*Catalog have velocities?*/

/*Check the velocity option*/
if(h_vel != 1 && h_vel != 0){
	printf("You need to say if your particle catalog have velocities! Yes (1) or No (0)?\n");
	exit(0);
}

/*Check if use voids/halos that grows over others*/
if(flag != 1 && flag != 0){
	printf("You need to say if you use halo that grows over other! Yes (1) or No (0)?\n");
	exit(0);
}

/*open the Profiles file*/
profiles = fopen(profilesfile,"r");
if (profiles == NULL) {
	printf("Unable to open %s\n", profilesfile);
	exit(0);
}

fscanf(profiles, "%d", &nh);	/*Total number of voids/halos*/
fscanf(profiles, "%f", &n_m);	/*Mean number of particle in (Mpc/h)^-3*/
fscanf(profiles, "%d", &bins);	/*Number of bins in input file*/
fscanf(profiles, "%d", &Do_log);/*information about the spacing. log (1) or linear (0) */
fscanf(profiles, "%f", &r_times); /*Catch as many times the radius of the voids/halos*/

/*Check the bin spacing option*/
if(Do_log != 1 && Do_log != 0){
	printf("You need to say if your bins are log spaced! Yes (1) or No (0)?\n");
	exit(0);
}

density = (float **)malloc(nh*sizeof(float *));		/*Array with the density array for each halo*/
velocity = (float **)malloc(nh*sizeof(float *));

for(i=0;i<nh;i++){
	density[i] = (float *)malloc(bins*sizeof(float));	/*Array with density in each bin*/
	velocity[i] = (float *)malloc(bins*sizeof(float));
}

for(i=0;i<nh;i++)
	for(j=0;j<bins;j++){
		density[i][j] = 0.0;
		velocity[i][j] = 0.0;
	}

/*Allocating the final density array*/
den = (double *)malloc(bins*sizeof(double));
vel = (double *)malloc(bins*sizeof(double));

for(j=0;j<bins;j++){
	den[j] = 0.0;
	vel[j] = 0.0;
}

/*Allocating the error in density and velocity*/
err = (double *)malloc(bins*sizeof(double));
errv = (double *)malloc(bins*sizeof(double));

for(j=0;j<bins;j++){
	err[j] = 0.0;
	errv[j] = 0.0;
}

/*Allocating the temporary in density*/
r_mean = 0.0;		/*The mean radius in the stack*/
cont = 0; 		/*Number of voids/halos in stack*/

len_b = r_times/bins;
log_min = log10(len_b);
log_b = (log10(r_times) - log_min)/(bins-1);

for(i=0;i<nh;i++){

	fscanf(profiles, "%d %f %d %f", &label, &radius, &f, &r_times_v);	/*The label, radius  and flag of halo*/

	if((f==0 || (f==1 && flag==1 )) && radius>=R_min && radius<=R_max){
		for(j=0;j<bins;j++){
			
			if(h_vel==1)	fscanf(profiles, "%f %d %f", &r, &n, &v);
			else{
				fscanf(profiles, "%f %d", &r, &n);
				v = 0.0; 
			}

			if(j==0){
				if(Do_log == 1) density[cont][j] += ((float) n)/Vol(pow(10, log_min)*radius);
				else density[cont][j] += ((float) n)/Vol(len_b*radius);	
			}

			else{
				if(Do_log == 1) density[cont][j] += ((float) n)/(Vol(pow(10, log_min + j*log_b)*radius) - Vol(pow(10, log_min + (j-1)*log_b)*radius));
				else density[cont][j] += ((float) n)/(Vol((j+1)*len_b*radius) - Vol(j*len_b*radius));
			}

			velocity[cont][j] += v;
		}
		cont ++;
		r_mean += radius;
	}
	else{
		if(h_vel==1)
			for(j=0;j<bins;j++)
				fscanf(profiles, "%f %d %f", &r, &n, &v);
		else
			for(j=0;j<bins;j++)
				fscanf(profiles, "%f %d", &r, &n);
	}
}

fclose(profiles);

/*Normalize the density*/
for(i=0;i<cont;i++)
	for(j=0;j<bins;j++)
		density[i][j] = density[i][j]/(n_m);

/*Stacking of all halos*/
for(j=0;j<bins;j++)
	for(i=0;i<cont;i++){
		den[j] += ((double)density[i][j])/cont;
		vel[j] += ((double)velocity[i][j])/cont;
	}

/*Evaluating the error in profiles*/
for(i=0;i<bins;i++)	
	for(j=0;j<cont;j++){
		err[i] = (double)pow(den[i] - density[j][i], 2)/(cont-1.0);
		errv[i] = (double)pow(vel[i] - velocity[j][i], 2)/(cont-1.0);
}

/*Evaluate the error of the mean*/
for(i=0;i<bins;i++){
	err[i] = err[i]/cont;
	errv[i] = errv[i]/cont;
}

/*Open the output file*/
out = fopen(outfile,"w");
if (out == NULL) {
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/*Print the stacked profile*/
fprintf(out,"%d\n", cont);
fprintf(out,"%f %f %f\n", R_min, R_max, r_mean/cont);
for(j=0;j<bins;j++){
	if(Do_log == 1){
		if(h_vel==1)	fprintf(out,"%f %lf %lf %lf %lf\n", pow(10, log_min + j*log_b), den[j], sqrt(err[j]), vel[j], sqrt(errv[j]));
		else	fprintf(out,"%f %lf %lf\n", pow(10, log_min + j*log_b), den[j], sqrt(err[j]));
	}
	else{
		if(h_vel == 1)	fprintf(out,"%f %lf %lf %lf %lf\n", (j+1)*len_b, den[j], sqrt(err[j]), vel[j], sqrt(errv[j]));
		else	fprintf(out,"%f %lf %lf\n", (j+1)*len_b, den[j], sqrt(err[j]));
	}
}
fclose(out);

return 0;
}
