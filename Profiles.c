#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define PI 3.141592

/*Structure for the particles in the original catalog*/
typedef struct Particle {
	int label; 		/*Label of the particle*/
	float pos[3];		/*Array with the position of the particle*/
	float vel[3];		/*Array with the velocity of the particle*/
} PARTICLE;

/*Structure for the voids/halos*/
typedef struct void_halo {
	int label; 		/*Label of the structure*/
	float pos[3];		/*Array with the position of the center of structure*/
	float rad;		/*Radius of structure*/
	int n;			/*Number of particles in structure*/
	float den;		/*Density of structure*/
} VH;

/*Structure for the particles around each void*/
typedef struct Radial{
	float dist;		/*Distance to the structure center*/
	float vr;		/*Radial component of the velocity particle*/
} RAD;

/*Returns The mean of vector x*/
float mean(float *x, float sum_vol, int N){
	int i;
	float resp;
	if(N==0)
		resp = 0.0;
	else{
		resp = 0.0;
		for(i=0;i<N;i++)
			resp += x[i];
		resp = resp/sum_vol;
	}

	return resp;
}

/*Find correct grid box of a particle*/
void indices(float x[], int xt[], float Lb, int nd){
	xt[0] = (int) x[0]/Lb;
	xt[1] = (int) x[1]/Lb;
	xt[2] = (int) x[2]/Lb;
	if(xt[0]==nd)	xt[0] -= nd;
	if(xt[1]==nd)	xt[1] -= nd;
	if(xt[2]==nd)	xt[2] -= nd;
}

/*Free the counts*/
void free_cont(int ***cont, int nd){
	int i, j;
	
	for(i=0;i<nd;i++)
		for(j=0;j<nd;j++)
			free(cont[i][j]);

	for(i=0;i<nd;i++)
		free(cont[i]);

	free(cont);
}

/*Free the grid of particles*/
void free_m(PARTICLE ****m, int nd){
	int i, j, k;
	
	for(i=0;i<nd;i++)
		for(j=0;j<nd;j++)
			for(k=0;k<nd;k++)
				free(m[i][j][k]);

	for(i=0;i<nd;i++)
		for(j=0;j<nd;j++)
			free(m[i][j]);
	for(i=0;i<nd;i++)
		free(m[i]);

	free(m);
}

int main(int argc,char *argv[])
{
FILE *part, *structures, *out;
char partfile[100], structurefile[100], outfile[100];
int np, np2, nv, i, j, k, a, b , c, n, flag, bin, res, nd, nb, h_vel, Do_log;
int ***cont1, ***cont2, xt[3], N_cores, percent, div;
float dist, n_m, r_times, len_b, L, L_c, Lb, trash, Delta;
float volume;
PARTICLE *p, ****m;
VH *v;

if (argc != 11){
	printf("Wrong number of arguments.\n");	
	printf("arg1: Name of the particles file\n");
	printf("arg2: Name of the voids/halos file\n");
	printf("arg3: Name of the output file\n");
	printf("arg4: Box size\n");
	printf("arg5: Catch as many times the radius of the voids/halos\n");
	printf("arg6: The number of bins in the range r_times*radius\n");
	printf("arg7: The particle catalog also have velocities? Yes (1) or No (0).\n");
	printf("arg8: Use the log spacing in radial bins? Yes (1) or No (0).\n");
	printf("arg9: Give the number of cores for the parallel computation.\n");
	printf("arg10. The mean density contrast of the voids/halos in the catalog.\n\n");	
	exit(0);
}

sprintf(partfile,"%s", argv[1]);
sprintf(structurefile,"%s", argv[2]);
sprintf(outfile,"%s", argv[3]);

/*read the parameters given by the user*/
r_times = atof(argv[5]);	/*How much grows the radius for profile*/
L = atof(argv[4]);		/*Box size*/
res = atof(argv[6]);		/*The number of bins in the range r_times*radius*/
h_vel = atoi(argv[7]);		/*Catalog have velocities?*/
Do_log = atoi(argv[8]);		/*Use the radial bins log spaced*/
N_cores = atoi(argv[9]);	/*Number of cores for the parallelization*/
Delta = atof(argv[10]);		/*Density contrast of the voids/halos*/

/*Check the velocity option*/
if(h_vel != 1 && h_vel !=0){
	printf("You need to say if your particle catalog have velocities! Yes (1) or No (0)?\n");
	exit(0);
}

/*Check the bin spacing option*/
if(Do_log != 1 && Do_log != 0){
	printf("You need to say if your bins are log spaced! Yes (1) or No (0)?\n");
	exit(0);
}

/*Open the particles file*/
part = fopen(partfile,"r");
if (part == NULL) {
	printf("Unable to open %s\n",partfile);
	exit(0);
}

/*Read the total number of particles and alocate they*/
fscanf(part,"%d", &np);	
p = (PARTICLE *)malloc(np*sizeof(PARTICLE));

/*Some parameters of simulation and grid*/
nd = (int)floor(pow(np,1.0/3.0)/2+0.9);	/*Number of divisons for the allocation*/
Lb = L/nd;				/*Size of each sub-box*/		
nb = (int)floor(L_c/Lb + 0.99);		/*Number of grids that contain L_c*/
n_m = np/(L*L*L);			/*The mean density in (Mpc/h)^-3*/

/*Open the voids/halos file*/
structures = fopen(structurefile,"r");
if (structures == NULL) {
	printf("Unable to open %s\n",structurefile);
	exit(0);
}

/*Allocating the counters*/
cont1 = (int ***)malloc(nd*sizeof(int **));
cont2 = (int ***)malloc(nd*sizeof(int **));

for(i=0;i<nd;i++){
	cont1[i] = (int **)malloc(nd*sizeof(int *));
	cont2[i] = (int **)malloc(nd*sizeof(int *));
}

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++){
		cont1[i][j] = (int *)malloc(nd*sizeof(int));
		cont2[i][j] = (int *)malloc(nd*sizeof(int));
}

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++){
			cont1[i][j][k] = 0;
			cont2[i][j][k] = 0;
}

/*Read the position of all particles and determines yours sub-box*/
for (i=0;i<np;i++){
	p[i].label = i;
	for (j=0;j<3;j++)	
		fscanf(part,"%f", &p[i].pos[j]);
	for (j=0;j<3 && (h_vel == 1);j++)	
		fscanf(part,"%f", &p[i].vel[j]);
	indices(p[i].pos, xt, Lb, nd);	
	cont1[xt[0]][xt[1]][xt[2]] += 1;
}
fclose(part);

/*Allocating m*/
m = (PARTICLE ****)malloc(nd*sizeof(PARTICLE ***));

for(i=0;i<nd;i++)
	m[i] = (PARTICLE ***)malloc(nd*sizeof(PARTICLE **));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		m[i][j] = (PARTICLE **)malloc(nd*sizeof(PARTICLE *));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++)
			m[i][j][k] = (PARTICLE *)malloc(cont1[i][j][k]*sizeof(PARTICLE));

/*Saving particles in m*/
for(i=0;i<np;i++){
	indices(p[i].pos, xt, Lb, nd);
	m[xt[0]][xt[1]][xt[2]][cont2[xt[0]][xt[1]][xt[2]]].label = p[i].label;
	for(j=0;j<3;j++){
		m[xt[0]][xt[1]][xt[2]][cont2[xt[0]][xt[1]][xt[2]]].pos[j] = p[i].pos[j];
		if(h_vel == 1) m[xt[0]][xt[1]][xt[2]][cont2[xt[0]][xt[1]][xt[2]]].vel[j] = p[i].vel[j];
	}
	cont2[xt[0]][xt[1]][xt[2]] +=1;
}
free(p);

/*Read the total number of voids/halos and alocate they*/
fscanf(structures,"%d", &nv);
div = nv/N_cores;	
v = (VH *)malloc(nv*sizeof(VH));
for(i=0;i<nv;i++){
	fscanf(structures,"%d", &v[i].label);
	for (j=0;j<3;j++)	
		fscanf(structures,"%f", &v[i].pos[j]); 
	fscanf(structures,"%f", &v[i].rad);	
	fscanf(structures,"%f", &trash);
	fscanf(structures,"%d", &v[i].n);
	fscanf(structures,"%f", &v[i].den);	
}	
fclose(structures);	

/*Open the output file*/
out = fopen(outfile,"w");
if (out == NULL) {
	printf("Unable to open out\n");
	exit(0);
}

fprintf(out,"%d\n", nv);	/*Total number of voids/halos*/
fprintf(out,"%f\n", n_m);	/*The mean density of catalog*/
fprintf(out,"%d\n", res);	/*Number of bins for each void*/
fprintf(out,"%d\n", Do_log);	/*information about the spacing*/
fprintf(out,"%f\n", r_times);	/*How much the radius of the voids/halos are get*/

printf("There are %d halos/voids to compute the profile.\n", nv);
percent = 0;
/**************************************/
/*Start to parallelize from this point*/
/**************************************/
omp_set_num_threads(N_cores);
#pragma omp parallel for private(i, j, k, a, b, c, L_c, flag, len_b, nb, percent)
/*Calculate the desnsity profile of each halo*/
for (i=0;i<nv;i++){

	/*variables used in the main loop*/
	float *x_c, *dx, *vels, *vols, r_times_v, distance, log_min, log_dist, log_b;
	int x, y, z, *nbins, *x_t, bin, cont;
	RAD *r;

	/*Jump this void/halo if it have radius zero*/
	if(v[i].rad == 0.0){
		#pragma omp critical
		{
		/*Print the relevant void/halos information*/
		fprintf(out,"%d %f %d %f\n", v[i].label, v[i].rad, 0, 0.0);
	
		/*Prints the number of particle up to each bin radius*/
		for(j=0;j<res;j++){
			if(Do_log == 1){
				if(h_vel==1)	fprintf(out,"%f %d %f\n", 0.0, 0, 0.0);
				else	fprintf(out,"%f %d\n", 0.0, 0);
			}
			else{
				if(h_vel == 1)	fprintf(out,"%f %d %f\n", 0.0, 0, 0.0);
				else	fprintf(out,"%f %d\n", 0.0, 0);
			}
		}
		}
 		continue;
	}

	/*Jump this void/halo if it has a bigger radius*/
	if(r_times*v[i].rad > L/2){
		#pragma omp critical
		{
		/*Print the relevant void/halos information*/
		fprintf(out,"%d %f %d %f\n", v[i].label, v[i].rad, 0, 0.0);
	
		/*Prints the number of particle up to each bin radius*/
		for(j=0;j<res;j++){
			if(Do_log == 1){
				if(h_vel==1)	fprintf(out,"%f %d %f\n", 0.0, 0, 0.0);
				else	fprintf(out,"%f %d\n", 0.0, 0);
			}
			else{
				if(h_vel == 1)	fprintf(out,"%f %d %f\n", 0.0, 0, 0.0);
				else	fprintf(out,"%f %d\n", 0.0, 0);
			}
		}
		}
 		continue;
	}

	/*Array with the possible particles in the profile*/
	r = (RAD *)malloc(np*sizeof(RAD));

	/*Allocate the variables*/
	x_c = (float *)malloc(3*sizeof(float));
	dx = (float *)malloc(3*sizeof(float));
	nbins = (int *)malloc(res*sizeof(int));
	vols = (float *)malloc(res*sizeof(float));
	vels = (float *)malloc(res*sizeof(float));
	x_t = (int *)malloc(3*sizeof(int));

	/*Bin lenght in normalized units*/
	r_times_v = r_times;
	len_b = r_times_v/res;			

	/*Values needed for the log spacing*/
	log_min = log10(len_b);
	log_b = (log10(r_times_v) - log_min)/(res-1);

	/*Initialize the number, volume and radial velocity of particles in the profile */
	for(j=0;j<res;j++){
		nbins[j] = 0;	
		vols[j] = 0.0;
		vels[j] = 0.0;
	}

	/*The center of void/halo*/
	for(k=0;k<3;k++)	
		x_c[k] = v[i].pos[k];

	/*Determine the bin of center of the void/halo*/
	cont = 0;
	indices(x_c, x_t, Lb, nd);
	nb = (int)floor(r_times_v*v[i].rad/Lb + 0.99);

	for(a=-nb;a<=nb;a++)
		for(b=-nb;b<=nb;b++)
			for(c=-nb;c<=nb;c++){
				x = x_t[0] + a;
				y = x_t[1] + b;
				z = x_t[2] + c;
				if(x<0)	x += nd; if(x>=nd) x -= nd;
				if(y<0)	y += nd; if(y>=nd) y -= nd;
				if(z<0)	z += nd; if(z>=nd) z -= nd;

				for(j=0;j<cont1[x][y][z];j++){

					for (k=0;k<3;k++){
						dx[k] = m[x][y][z][j].pos[k] - x_c[k];	
						if(dx[k] < -L/2.0)
							dx[k] = dx[k] + L;
						if(dx[k] > L/2.0)
							dx[k] = dx[k] - L;				
					}
					distance = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]); /*Evaluate the distance beetween the particle and the halo center*/
					if(distance < r_times_v*v[i].rad){ /*If the distance is small than r_times times the halo radius get the main information about the particle*/
						r[cont].dist = distance;
						if(h_vel == 1 && distance > 0.0) r[cont].vr = (dx[0]*m[x][y][z][j].vel[0] + dx[1]*m[x][y][z][j].vel[1] + dx[2]*m[x][y][z][j].vel[2])/distance;
						cont ++;
					}
				}
			}

	/*Determines the bin of each particle (in log or linear spacing) and cont this particle*/
	for(j=0;j<cont;j++){
		if(Do_log == 1){
			log_dist = log10(r[j].dist/v[i].rad);	
			bin = (int)floor((log_dist - log_min)/log_b+0.99);
			if(bin<0) bin = 0;
		}
		else{
			bin = (int)floor(r[j].dist/v[i].rad/len_b);
		}
		nbins[bin] ++;
		if(h_vel==1) vels[bin] += r[j].vr;
		vols[bin] += 1.0;
	}
	
	/*Normalize the weighted velocities*/
	for(j=0;j<res && (h_vel==1);j++){
		if(vols[j] == 0.0) vels[j] = 0.0;
		else vels[j] = vels[j]/vols[j];
	}

	/*This flag denotes if this void/halos grows on other halo*/
	if(Delta > 1.0 && v[i].den > Delta)	flag = 1;
	else if(Delta > 1.0 && v[i].den <= Delta)	flag = 0;
	else if(Delta < 1.0 && v[i].den < Delta)	flag = 1;
	else if(Delta < 1.0 && v[i].den >= Delta)	flag = 0;

	#pragma omp critical
	{
	/*Print the relevant void/halos information*/
	fprintf(out,"%d %f %d %f\n", v[i].label, v[i].rad, flag, r_times_v);
	
	/*Prints the number of particle up to each bin radius*/
	for(j=0;j<res;j++){
		if(Do_log == 1){
			if(h_vel==1)	fprintf(out,"%f %d %f\n", pow(10, log_min + j*log_b), nbins[j], vels[j]);
			else	fprintf(out,"%f %d\n", pow(10, log_min + j*log_b), nbins[j]);
		}
		else{
			if(h_vel == 1)	fprintf(out,"%f %d %f\n", (j+1)*len_b, nbins[j], vels[j]);
			else	fprintf(out,"%f %d\n", (j+1)*len_b, nbins[j]);
		}
	}
	}

	free(dx);
	free(x_c);
	free(x_t);
	free(nbins);
	free(vels);
	free(vols);
	free(r);

	if(i<div && percent < i*100/div){
		percent = i*100/div;
		printf("%d%% of the halos/voids profiles computed\n", percent);
	}
}

/**************************/
/*Stop the parallelization*/
/**************************/
fclose(out);

/*Free the memory*/
free(v);
free_cont(cont1, nd);
free_cont(cont2, nd);
free_m(m, nd);

return 0;
}
