#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Define the cyclic sum*/
int cysum(int i, int j, int nd){
	int resp;

	resp = i+j;
	if(resp>=nd)	resp -= nd;
	if(resp<0)	resp += nd;

	return resp;
}

/*Define the cyclic sum for floats*/
float cysumf(float x, float y, float L){
	float resp;

	resp = x + y;
	if(resp>=L)	resp -= L;
	if(resp<0)	resp += L;

	return resp;
}

/*Give a indice for the partricle*/
void indices(float x[], int xt[], float Lb, int nd){
	xt[0] = (int) x[0]/Lb;
	xt[1] = (int) x[1]/Lb;
	xt[2] = (int) x[2]/Lb;
	if(xt[0]==nd)	xt[0] -=1;
	if(xt[1]==nd)	xt[1] -=1;
	if(xt[2]==nd)	xt[2] -=1;
}

int main(int argc,char *argv[])
{
FILE *part, ****out;
char partfile[100], outfile[100];
int i, j, k, l, nd, np, nptot, ix, iy, iz;
int ***cont, xt[3];
float L, Lb, Lc, xmin, xmax, ymin, ymax, zmin, zmax;
float p[3];

/*Check all inputs*/
if (argc != 4){
	printf("Wrong number of arguments.\n");
	printf("arg1: The preffix of the input file\n");
	printf("arg2: Number of divisions (at least 2)\n");
	printf("arg3: Extra size in the splits (in Mpc/h)\n\n");
	exit(0);
}

/*Read the inputs*/
sprintf(partfile, "%s", argv[1]);
nd = atoi(argv[2]);
Lc = atof(argv[3]);
if(nd < 2) nd = 2;

/*Open the particles file*/
part = fopen(partfile,"rb");
if (part == NULL) {
	printf("Unable to open %s\n", partfile);
	exit(0);
}

/*Read the total number of particles and the size of the box*/
fread(&np, sizeof(int), 1, part);	
fread(&L, sizeof(int), 1, part);
Lb = L/nd;

/*Alloc the outputs*/
out = (FILE ****)malloc(nd*sizeof(FILE ***));

for(i=0;i<nd;i++)
	out[i] = (FILE ***)malloc(nd*sizeof(FILE **));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		out[i][j] = (FILE **)malloc(nd*sizeof(FILE *));

/*Open the output files*/
nptot = 0;
for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++){
			sprintf(outfile, "%s_%d_%d_%d", argv[1], i, j, k);

			out[i][j][k] = fopen(outfile, "wb");
			if (out[i][j][k] == NULL) {
				printf("Unable to open %s\n", outfile);
				exit(0);
			}

			fwrite(&nptot, sizeof(int), 1, out[i][j][k]);
			fwrite(&np, sizeof(int), 1, out[i][j][k]);
			fwrite(&L, sizeof(float), 1, out[i][j][k]);

			xmin = i*Lb; xmax = i*Lb+Lb;
			ymin = j*Lb; ymax = j*Lb+Lb;
			zmin = k*Lb; zmax = k*Lb+Lb;

			fwrite(&xmin, sizeof(float), 1, out[i][j][k]);
			fwrite(&xmax, sizeof(float), 1, out[i][j][k]);
			fwrite(&ymin, sizeof(float), 1, out[i][j][k]);
			fwrite(&ymax, sizeof(float), 1, out[i][j][k]);
			fwrite(&zmin, sizeof(float), 1, out[i][j][k]);
			fwrite(&zmax, sizeof(float), 1, out[i][j][k]);
		}


/*Allocating the counts*/
cont = (int ***)malloc(nd*sizeof(int **));

for(i=0;i<nd;i++)
	cont[i] = (int **)malloc(nd*sizeof(int *));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		cont[i][j] = (int *)malloc(nd*sizeof(int));

for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++)
			cont[i][j][k] = 0;

/*Read the position of all particles, determines your sub-box(es) and save them in the output files*/
for (i=0;i<np;i++){
	for(j=0;j<3;j++)	
		fread(&p[j], sizeof(float), 1, part);

	indices(p, xt, Lb, nd);	

	for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[xt[0]][xt[1]][xt[2]]);
	cont[xt[0]][xt[1]][xt[2]] += 1;
	nptot ++;

	ix = 0;
	iy = 0;
	iz = 0;

	if(p[0] - xt[0]*Lb < Lc) ix = -1;
	if(p[0] - xt[0]*Lb > Lb-Lc) ix = 1;

	if(p[1] - xt[1]*Lb < Lc) iy = -1;
	if(p[1] - xt[1]*Lb > Lb-Lc) iy = 1;

	if(p[2] - xt[2]*Lb < Lc) iz = -1;
	if(p[2] - xt[2]*Lb > Lb-Lc) iz = 1;

	if(ix != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[cysum(xt[0], ix, nd)][xt[1]][xt[2]]);
		cont[cysum(xt[0], ix, nd)][xt[1]][xt[2]] += 1;
		nptot ++;
	}
	if(iy != 0){
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[xt[0]][cysum(xt[1], iy, nd)][xt[2]]); 
		cont[xt[0]][cysum(xt[1], iy, nd)][xt[2]] += 1; 
		nptot ++;
	}
	if(iz != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[xt[0]][xt[1]][cysum(xt[2], iz, nd)]);
		cont[xt[0]][xt[1]][cysum(xt[2], iz, nd)] += 1; 
		nptot ++;
	}
	if(ix != 0 && iy != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[cysum(xt[0], ix, nd)][cysum(xt[1], iy, nd)][xt[2]]);
		cont[cysum(xt[0], ix, nd)][cysum(xt[1], iy, nd)][xt[2]] += 1; 
		nptot ++;
	}
	if(ix != 0 && iz != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[cysum(xt[0], ix, nd)][xt[1]][cysum(xt[2], iz, nd)]);
		cont[cysum(xt[0], ix, nd)][xt[1]][cysum(xt[2], iz, nd)] += 1; 
		nptot ++;
	}
	if(iy != 0 && iz != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[xt[0]][cysum(xt[1], iy, nd)][cysum(xt[2], iz, nd)]);
		cont[xt[0]][cysum(xt[1], iy, nd)][cysum(xt[2], iz, nd)] += 1; 
		nptot ++;
	}
	if(ix != 0 && iy != 0 && iz != 0){ 
		for(j=0;j<3;j++) fwrite(&p[j], sizeof(float), 1, out[cysum(xt[0], ix, nd)][cysum(xt[1], iy, nd)][cysum(xt[2], iz, nd)]);
		cont[cysum(xt[0], ix, nd)][cysum(xt[1], iy, nd)][cysum(xt[2], iz, nd)] += 1; 
		nptot ++;
	}
}
fclose(part);

printf("Original number of particles: %d and total number of particles: %d\n", np, nptot);

/*Save the total numbar of particles in each output file*/
for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		for(k=0;k<nd;k++){
			rewind(out[i][j][k]);
			fwrite(&cont[i][j][k], sizeof(float), 1, out[i][j][k]);
		}

/*Free count*/
for(i=0;i<nd;i++)
	for(j=0;j<nd;j++)
		free(cont[i][j]);

for(i=0;i<nd;i++)
	free(cont[i]);

free(cont);

/*Close the outputs*/
for(i=0;i<nd;i++){
	for(j=0;j<nd;j++){
		for(k=0;k<nd;k++)
			fclose(out[i][j][k]);
		free(out[i][j]);
	}
	free(out[i]);
}
free(out);

return 0;
}
