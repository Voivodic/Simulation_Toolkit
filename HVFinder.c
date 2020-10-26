#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define PI 3.141592
#define SUB 1e+10

/*Structure for the particles in the original catalog*/
typedef struct Particle {
	float pos[3];		/*Array with the position of the particle*/
} PARTICLE;

/*List with the particles around the center of the Spherical void*/
typedef struct Radius {
	int label;		/*Label of teh particle*/
	float dist;		/*Distance for the center*/
} RAD;

/*Structure used to allocate the particles in the grid*/
typedef struct Grid_Particles {
	int label;		/*Label of the particle*/
	float pos[3];		/*Array with the position of the particle*/
} GRID;

/*Structure used to temporary allocate the centers*/
typedef struct Center {
	float pos[3];		/*Array with the position of the center*/
	float den;		/*Density of the central particle in the halo*/
	int label;		/*Label of the central particle*/
} CENTER;

/*Structure used to allocate the centers*/
typedef struct Center2 {
	float pos[3];		/*Array with the position of the center*/
	int label;		/*Label of the central particle*/
} CENTER2;

/*Structure with the final Halos/Voids*/
typedef struct Halo_Void {
	int label;
	float pos[3];
	float R;
	float errR;
	int np;
	float den;
} HV;

/*Give a indice for the partricle*/
void indices(float x[], int xt[], float Lb, int nd){
	xt[0] = (int) x[0]/Lb;
	xt[1] = (int) x[1]/Lb;
	xt[2] = (int) x[2]/Lb;
	if(xt[0]==nd)	xt[0] -=1;
	if(xt[1]==nd)	xt[1] -=1;
	if(xt[2]==nd)	xt[2] -=1;
}

/*Partition function used in the radius quicksort*/
int partition( RAD a[], int l, int r) {
   	int i, j;
	RAD pivot, t;
  	pivot.label = a[l].label;
  	pivot.dist = a[l].dist;	
   	i = l; j = r+1;
		
   	while( 1){
   	do ++i; while( a[i].dist <= pivot.dist && i <= r );
   	do --j; while( a[j].dist > pivot.dist );
   	if( i >= j ) break;
   	t.dist = a[i].dist; a[i].dist = a[j].dist; a[j].dist = t.dist;
   	t.label = a[i].label; a[i].label = a[j].label; a[j].label = t.label;	
   	}
   	t.dist = a[l].dist; a[l].dist = a[j].dist; a[j].dist= t.dist;
	t.label = a[l].label; a[l].label = a[j].label; a[j].label = t.label;
   	return j;
}

/*The quicksort algorithm to sort the RAD list*/
void quickSort( RAD a[], int l, int r){
	int j;

   	if( l < r ){
   	// divide and conquer
        j = partition( a, l, r);
       	quickSort( a, l, j-1);
      	quickSort( a, j+1, r);
   	}	
}

/*Partition function used in the radius quicksort*/
int partition_center(CENTER a[], int l, int r) {
   	int i, j, k;
	CENTER pivot, t;
  	pivot.label = a[l].label;
  	pivot.den = a[l].den;	
   	i = l; j = r+1;
		
   	while( 1){
   	do ++i; while( a[i].den <= pivot.den && i <= r );
   	do --j; while( a[j].den > pivot.den );
   	if( i >= j ) break;
   	t.den = a[i].den; a[i].den = a[j].den; a[j].den = t.den;
   	t.label = a[i].label; a[i].label = a[j].label; a[j].label = t.label;
	for(k=0;k<3;k++){ 
		t.pos[k] = a[i].pos[k]; a[i].pos[k] = a[j].pos[k]; a[j].pos[k] = t.pos[k];
	}
   	}
   	t.den = a[l].den; a[l].den = a[j].den; a[j].den= t.den;
	t.label = a[l].label; a[l].label = a[j].label; a[j].label = t.label;
	for(k=0;k<3;k++){
		t.pos[k] = a[l].pos[k]; a[l].pos[k] = a[j].pos[k]; a[j].pos[k] = t.pos[k];
	}
   	return j;
}

/*The quicksort algorithm to sort the center list*/
void quickSort_center(CENTER a[], int l, int r){
	int j;

   	if( l < r ){
   	// divide and conquer
        j = partition_center( a, l, r);
       	quickSort_center( a, l, j-1);
      	quickSort_center( a, j+1, r);
   	}	
}

int main(int argc,char *argv[])
{
FILE *part, *outh, *outv, *crit, *outp;
char partfile[100], outhfile[100], outvfile[100], critfile[100], propfile[100];
int i, j, k, l, a, b, c, cont, h_vel, N_cores, Ns, np, np2, nbh, nbv, nd, nh, nv, DOhv, n, div, rest, conth, contv;
int ***cont1, ***cont2, xt[3], *flag, percent;
float L, Lch, Lcv, Lb, Lm, n_m, denTOT, Dh, Dv, xc[3], den;
PARTICLE *p;
GRID ****m;
HV *hv;
CENTER *centerh, *centerv;
CENTER2 **cen;

/*Check all inputs*/
if (argc != 10){
	printf("Wrong number of arguments.\n");
	printf("arg1: The preffix of the input file\n");
	printf("arg2: The preffix of the output file\n");
	printf("arg3: Find the Halos (0), Voids (1) or both (2)?\n");
	printf("arg4: Halo's caracteristic lenght.\n");
	printf("arg5: Void's caracteristic lenght.\n");
	printf("arg6: Halo's contrast density.\n");
	printf("arg7: Void's contrast density.\n");
	printf("arg8: Number of critical files.\n");
	printf("arg9: Give the number of cores for the parallel computation.\n\n");
	exit(0);
}

/*Get the name of all files*/
sprintf(partfile,"%s", argv[1]);
sprintf(outhfile,"%s_halos.dat", argv[2]);
sprintf(outvfile,"%s_voids.dat", argv[2]);
sprintf(propfile,"%s_props.dat", argv[2]);

/*Read the parameters given by the user*/
DOhv = atoi(argv[3]);		/*Finde halos, voids or both parameter*/
Lch = atof(argv[4]);		/*Caracteristic lenght of halos*/
Lcv = atof(argv[5]);		/*Caracteristic lenght of halos*/
Dh = atof(argv[6]);		/*Halo's contrast density*/
Dv = atof(argv[7]);		/*Void's contrast density*/
Ns = atoi(argv[8]);		/*Number critical files to be read*/
N_cores = atoi(argv[9]);	/*Number of cores for the parallelization*/

/*Check if the kind of finder is ok*/
if(DOhv != 0 && DOhv !=1 && DOhv != 2){
	printf("You need to say the kind of structure that you need to find, Halos (0), Voids (1) or both (2)?\n");
	exit(0);
}

/*Open the particles file*/
part = fopen(partfile,"rb");
if (part == NULL) {
	printf("Unable to open %s\n",partfile);
	exit(0);
}

/*Read the total number of particles and the size of teh box*/
fread(&np, sizeof(int), 1, part);	
fread(&L, sizeof(float), 1, part);	

/*Some parameters of simulation and grid*/
nd = (int)floor(pow(np,1.0/3.0)/2+0.5);	/*Number of divisons for the allocation*/
Lb = L/nd;				/*Size of each sub-box*/		
nbh = (int)floor(Lch/Lb + 0.9);		/*Number of grids that contain Lch*/
nbv = (int)floor(Lcv/Lb + 0.9);		/*Number of grids that contain Lcv*/
n_m = np/(L*L*L);			/*The mean density in (Mpc/h)^-3*/
Lm = L/pow(np,1.0/3.0); 		/*The mean separation between particles*/

printf("The mean density of the catalog is %f\n", n_m);
printf("np = %d and nd = %d\n", np, nd);

/*Allocating the counts*/
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
p = (PARTICLE *)malloc(np*sizeof(PARTICLE));
for (i=0;i<np;i++){
	for(j=0;j<3;j++)	
		fread(&p[i].pos[j], sizeof(float), 1, part);

	indices(p[i].pos, xt, Lb, nd);	
	cont1[xt[0]][xt[1]][xt[2]] += 1;
}
fclose(part);

/*Allocating the grid estructure*/
m = (GRID ****)malloc(nd*sizeof(GRID ***));

for(i=0;i<nd;i++)
	m[i] = (GRID ***)malloc(nd*sizeof(GRID **));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		m[i][j] = (GRID **)malloc(nd*sizeof(GRID *));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++)
			m[i][j][k] = (GRID *)malloc(cont1[i][j][k]*sizeof(GRID));

/*Array with the information if the particle is in a halo or a void*/
flag = (int *)malloc(np*sizeof(int));	

/*Saving particles in the grid structure m (used to find the nearst particles faster)*/
for(i=0;i<np;i++){
	flag[i] = 0;
	indices(p[i].pos, xt, Lb, nd);
	m[xt[0]][xt[1]][xt[2]][cont2[xt[0]][xt[1]][xt[2]]].label = i;
	for(j=0;j<3;j++)
		m[xt[0]][xt[1]][xt[2]][cont2[xt[0]][xt[1]][xt[2]]].pos[j] = p[i].pos[j];
	cont2[xt[0]][xt[1]][xt[2]] += 1;
}

/*Free p*/
free(p);

/*Free the count2*/
for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		free(cont2[i][j]);

for(i=0;i<nd;i++)
	free(cont2[i]);

free(cont2);

/*Read the total number of halos and voids and allocate them*/
nh = 0;
nv = 0;
printf("Ns = %d\n", Ns);
for(i=0;i<Ns;i++)
	for(j=0;j<Ns;j++)
		for(k=0;k<Ns;k++){
			sprintf(critfile, "%s_%d_%d_%d.crit", argv[1], i, j, k);

			crit = fopen(critfile, "rb");
			if (crit == NULL){
				printf("Unable to open critical file %s.\n\n", critfile);
				exit(0);
			}

			fread(&a, sizeof(int), 1, crit);
			fread(&b, sizeof(int), 1, crit);

			fclose(crit);

			nh += a;
			nv += b;
		}

printf("There are %d halos and %d voids\n", nh, nv);

/*Alloc the arrays with the centers*/
if(DOhv == 0 || DOhv == 2) centerh = (CENTER *)malloc(nh*sizeof(CENTER));	/*Array with the halo's centers*/
if(DOhv == 1 || DOhv == 2) centerv = (CENTER *)malloc(nv*sizeof(CENTER));	/*Array with the void's centers*/

/*Read the information about peaks and troughs*/
conth = 0;
contv = 0;
for(i=0;i<Ns;i++)
	for(j=0;j<Ns;j++)
		for(k=0;k<Ns;k++){
			sprintf(critfile, "%s_%d_%d_%d.crit", argv[1], i, j, k);

			crit = fopen(critfile, "rb");
			if (crit == NULL){
				printf("Unable to open critical file %s.\n\n", critfile);
				exit(0);
			}

			fread(&a, sizeof(int), 1, crit);
			fread(&b, sizeof(int), 1, crit);
			
			for(c=0;c<(a + b);c++){
				for(l=0;l<3;l++)
					fread(&xc[l], sizeof(float), 1, crit);
				fread(&den, sizeof(float), 1, crit);

				if(den >= Dh && (DOhv == 0 || DOhv == 2)){
					centerh[conth].label = conth;
					for(l=0;l<3;l++)
						centerh[conth].pos[l] = xc[l];
					centerh[conth].den = den;	

					conth ++;
				}

				if(den <= Dv && (DOhv == 1 || DOhv == 2)){
					centerv[contv].label = contv;
					for(l=0;l<3;l++)
						centerv[contv].pos[l] = xc[l];
					centerv[contv].den = den;	

					contv ++;
				}
			}
		}
fclose(crit);

/*Check the number of halos*/
if((DOhv == 0 || DOhv == 2) && nh != conth){
	printf("There are some problem with the halo number! %d != %d\n", nh, conth);
	exit(0);
}

/*Check the number of voids*/
if((DOhv == 1 || DOhv == 2) && nv != contv){
	printf("There are some problem with the void number! %d != %d\n", nv, contv);
	exit(0);
}

/*************************/
/*Start to grow the halos*/
/*************************/
if(DOhv == 0 || DOhv == 2){

/*Sort the peaks using its density*/
quickSort_center(centerh, 0, nh-1);

/*Number of halos for each core*/
div = nh/N_cores;
rest = nh%N_cores;

/*Alloc the array with the final centers*/
cen = (CENTER2 **)malloc(N_cores*sizeof(CENTER2 *));

for(i=0;i<N_cores;i++)
	cen[i] = (CENTER2 *)malloc((div+1)*sizeof(CENTER2 ));

/*Put the centers in a convenient way for the paralization*/	
for(i=0;i<N_cores;i++)
	for(j=0;j<div;j++){
		for(k=0;k<3;k++)
			cen[i][j].pos[k] = centerh[(nh - 1) - (j*N_cores + i)].pos[k];
		cen[i][j].label = centerh[(nh - 1) - (j*N_cores + i)].label;
	}
		

for(i=0;i<rest;i++){
	for(k=0;k<3;k++)
		cen[i][div].pos[k] = centerh[(nh - 1) - (div*N_cores + i)].pos[k];
	cen[i][div].label = centerh[(nh - 1) - (div*N_cores + i)].label;
}

free(centerh);

/*Array with the final halos*/
hv = (HV *)malloc(nh*sizeof(HV));
np2 = 0;

printf("We have %d peaks to look!\n", nh);

omp_set_num_threads(N_cores);
#pragma omp parallel for private(i, j, k, l, cont, percent)
for(i=0;i<N_cores;i++){
percent = 0;
for(l=0;l<div+1;l++){

	/*Jump if end in this core*/
	if(l == div && i>=rest) continue;

	/*Declarates the main variables*/
	float *x_c, *dx, den, dist;
	int x, y, z, a, b, c, *xt, n;
	RAD *r;

	r = (RAD *)malloc(np*sizeof(RAD));
	x_c = (float *)malloc(3*sizeof(float));
	dx = (float *)malloc(3*sizeof(float));
	xt = (int *)malloc(3*sizeof(int));

	/*The center of the halo*/
	for(k=0;k<3;k++)	
		x_c[k] = cen[i][l].pos[k];

	cont = 0;
	indices(x_c, xt, Lb, nd);

	/*Evaluate the distance between all particles close to the center*/
	for(a=-nbh;a<=nbh;a++)
		for(b=-nbh;b<=nbh;b++)
			for(c=-nbh;c<=nbh;c++){
				x = xt[0] + a;
				y = xt[1] + b;
				z = xt[2] + c;
				if(x<0)	x += nd; if(x>=nd) x -= nd;
				if(y<0)	y += nd; if(y>=nd) y -= nd;
				if(z<0)	z += nd; if(z>=nd) z -= nd;
				
				for(j=0;j<cont1[x][y][z];j++){
					for (k=0;k<3;k++){
						dx[k] = x_c[k] - m[x][y][z][j].pos[k];	
						if(dx[k] < -L/2.0)
							dx[k] = dx[k] + L;
						if(dx[k] > L/2.0)
							dx[k] = dx[k] - L;
					}
					dist = 0.0;
					for (k=0;k<3;k++)
						dist += dx[k]*dx[k];
					dist = sqrt(dist);
					if(dist < Lch){	
						r[cont].label = m[x][y][z][j].label;
						r[cont].dist = dist;
						cont ++;
					}				
				}
			}	

	/*Ordening by distance to the core particle*/
	quickSort( r, 0, cont-1);

	den = 0.0;	/*Number density of the halo*/
	n = 1;		/*Number of particles in the halo*/
	a = 1;		/*Information about the three other central particles*/

	if(flag[r[0].label] != 0) a = 0;
		
	if(a == 1) flag[r[0].label] = cen[i][l].label + 1;
	

	/*Find the size and the number of particles of the halo*/
	while(a == 1){
		if(flag[r[n].label] != 0) break;
		den = 3.0*(n+1)/(4.0*PI*r[n].dist*r[n].dist*r[n].dist);
		if(den/n_m < Dh) break;
		flag[r[n].label] = cen[i][l].label + 1;
		n++;		
	}

	if(a == 1){
		#pragma omp critical
		{
		hv[np2].label = cen[i][l].label + 1;
		hv[np2].pos[0] = x_c[0];
		hv[np2].pos[1] = x_c[1];
		hv[np2].pos[2] = x_c[2];
		hv[np2].R = (r[n].dist + r[n-1].dist)/2;
		hv[np2].errR = (r[n].dist - r[n-1].dist)/2;
		hv[np2].np = n;
		hv[np2].den = den/n_m;
		np2 ++;
		}
	}

	free(x_c); free(dx); free(xt); free(r);

	if(percent < l*100/div){
		percent = l*100/div;
		printf("%d%% of the halos found in the porcessor %d\n", percent, i);
	}
}
}

/*Open the output halo's file*/
outh = fopen(outhfile,"w");
if (outh == NULL) {
	printf("Unable to open halo's output file %s.\n\n", outhfile);
	exit(0);
}

/*Print the main halo's information*/
fprintf(outh,"%d\n", np2);
for(i=0;i<np2;i++)
	fprintf(outh,"%d %f %f %f %f %f %d %f\n", hv[i].label, hv[i].pos[0], hv[i].pos[1], hv[i].pos[2], hv[i].R, hv[i].errR, hv[i].np, hv[i].den);

/*Free the center's positions*/
for(i=0;i<N_cores;i++)
	free(cen[i]);

free(cen);

fclose(outh);
free(hv);
}

/*****************************/
/*Start the grow of the voids*/
/*****************************/
if(DOhv == 1 || DOhv == 2){

/*Sort the troughts using its density*/
quickSort_center(centerv, 0, nv-1);

/*Number of voids for each core*/
div = nv/N_cores;
rest = nv%N_cores;

/*Alloc the array with the final centers*/
cen = (CENTER2 **)malloc(N_cores*sizeof(CENTER2 *));

for(i=0;i<N_cores;i++)
	cen[i] = (CENTER2 *)malloc((div+1)*sizeof(CENTER2 ));


/*Put the centers in a convenient way for the paralization*/	
for(i=0;i<N_cores;i++)
	for(j=0;j<div;j++){
		for(k=0;k<3;k++)
			cen[i][j].pos[k] = centerv[j*N_cores + i].pos[k];
		cen[i][j].label = centerv[j*N_cores + i].label;
	}

for(i=0;i<rest;i++){
	for(k=0;k<3;k++)
		cen[i][div].pos[k] = centerv[div*N_cores + i].pos[k];
	cen[i][div].label = centerv[div*N_cores + i].label;	
}

free(centerv);

/*Array with the final voids*/
hv = (HV *)malloc(nv*sizeof(HV));
np2 = 0;

printf("\nWe have %d troughs to look!\n", nv);

omp_set_num_threads(N_cores);
#pragma omp parallel for private(i, j, k, l, cont, percent)
for(i=0;i<N_cores;i++){
percent = 0;
for(l=0;l<div+1;l++){

	/*Jump if end in this core*/
	if(l == div && i>=rest) continue;

	/*Declarates the main variables*/
	float *x_c, *dx, den, dist;
	int x, y, z, a, b, c, *xt, n;
	RAD *r;

	r = (RAD *)malloc(np*sizeof(RAD));
	x_c = (float *)malloc(3*sizeof(float));
	dx = (float *)malloc(3*sizeof(float));
	xt = (int *)malloc(3*sizeof(int));

	/*The center of the void*/
	for(k=0;k<3;k++)	
		x_c[k] = cen[i][l].pos[k];

	cont = 0;
	indices(x_c, xt, Lb, nd);

	/*Evaluate the distance between all particles close to the center*/
	for(a=-nbv;a<=nbv;a++)
		for(b=-nbv;b<=nbv;b++)
			for(c=-nbv;c<=nbv;c++){
				x = xt[0] + a;
				y = xt[1] + b;
				z = xt[2] + c;
				if(x<0)	x += nd; if(x>=nd) x -= nd;
				if(y<0)	y += nd; if(y>=nd) y -= nd;
				if(z<0)	z += nd; if(z>=nd) z -= nd;
				
				for(j=0;j<cont1[x][y][z];j++){
					for (k=0;k<3;k++){
						dx[k] = x_c[k] - m[x][y][z][j].pos[k];	
						if(dx[k] < -L/2.0)
							dx[k] = dx[k] + L;
						if(dx[k] > L/2.0)
							dx[k] = dx[k] - L;
					}
					dist = 0.0;
					for (k=0;k<3;k++)
						dist += dx[k]*dx[k];
					dist = sqrt(dist);
					if(dist < Lcv){	
						r[cont].label = m[x][y][z][j].label;
						r[cont].dist = dist;
						cont ++;
					}				
				}
			}	

	/*Ordening by distance to the core particle*/
	quickSort( r, 0, cont-1);

	den = 0.0;	/*Number density of the halo*/
	n = 3;		/*Number of particles in the halo*/
	a = 1;		/*Information about the three other central particles*/

	for(k=0;k<3;k++)
		if(flag[r[k].label] != 0) a = 0;
		
	for(k=0;k<3 && a == 1;k++)
		flag[r[k].label] = -cen[i][l].label - 1;

	/*Find the size and the number of particles of the halo*/
	while(a == 1 && n < cont){
		if(flag[r[n].label] != 0) break;
		den = 3.0*(n+1)/(4.0*PI*r[n].dist*r[n].dist*r[n].dist);
		if(den/n_m > Dv) break;
		flag[r[n].label] = -cen[i][l].label - 1;
		n++;		
	}

	if(n == cont)
		printf("n = %d and cont = %d. r[cont].dist = %f\n", n, cont, r[cont-1].dist);

	if(a == 1){
		#pragma omp critical
		{	
		hv[np2].label = cen[i][l].label + 1;
		hv[np2].pos[0] = x_c[0];
		hv[np2].pos[1] = x_c[1];
		hv[np2].pos[2] = x_c[2];
		hv[np2].R = (r[n].dist + r[n-1].dist)/2;
		hv[np2].errR = (r[n].dist - r[n-1].dist)/2;
		hv[np2].np = n;
		hv[np2].den = den/n_m;
		np2 ++;
		}
	}

	free(x_c); free(dx); free(xt); free(r);

	if(percent < l*100/div){
		percent = l*100/div;
		printf("%d%% of the voids found in the porcessor %d\n", percent, i);
	}
}
}

/*Open the output void's file*/
outv = fopen(outvfile,"w");
if (outv == NULL) {
	printf("Unable to open the output file %s.\n\n", outvfile);
	exit(0);
}

/*Print the main void's information*/
fprintf(outv,"%d\n", np2);
for(i=0;i<np2;i++)
	fprintf(outv,"%d %f %f %f %f %f %d %f\n", hv[i].label, hv[i].pos[0], hv[i].pos[1], hv[i].pos[2], hv[i].R, hv[i].errR, hv[i].np, hv[i].den);

/*Free the center's positions*/
for(i=0;i<N_cores;i++)
	free(cen[i]);

free(cen);

fclose(outv);
free(hv);
}

/*Free the count1*/
for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		free(cont1[i][j]);

for(i=0;i<nd;i++)
	free(cont1[i]);

free(cont1);

/*Free the m structure*/
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

/*Open the particles properties file*/
outp = fopen(propfile, "w");
if (outp == NULL) {
	printf("Unable to open %s\n", propfile);
	exit(0);
}

/*Print the structure of the particle*/
fprintf(outp, "%d\n", np);
for(i=0;i<np;i++)
	fprintf(outp, "%d\n", flag[i]);

fclose(outp);

return 0;
}

