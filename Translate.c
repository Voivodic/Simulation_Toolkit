#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
FILE *original, *output;
char posfile[50], outfile[50];
int i, j, np, h_vel;
float **p, v[3], L;

if (argc != 5) {
	printf("Wrong number of arguments.\n");
	printf("arg1: Position file.\n");
	printf("arg2: Output binary file.\n");
	printf("arg3: The box size.\n");	
	printf("arg4: The particle catalog also have velocities? Yes (1) or No (0).\n\n");
	exit(0);
}
sprintf(posfile,"%s", argv[1]);
sprintf(outfile,"%s", argv[2]);
h_vel = atoi(argv[4]);
L = atof(argv[3]);

/*Check the velocity option*/
if(h_vel != 1 && h_vel != 0){
	printf("You need to say if your catalog has velocity. Yes (1) or No (0)?\n");
	exit(0);
}

if((original = fopen(posfile,"r")) == NULL){
	printf("\nErro ao abrir o arquivo original.\n\n");
        return 1;
}             

if((output = fopen(outfile, "wb")) == NULL){
       	printf("\nErro ao abrir o arquivo de saida.\n\n");
      	return 1;
}

/*Read the total number of particles and alocate they*/
fscanf(original,"%d", &np);
p = (float **)malloc(np*sizeof(float *));
for(i=0;i<np;i++)
	p[i] = (float *)malloc(3*sizeof(float));

/*Read the position of all particles*/
if(h_vel == 1){
	for(i=0;i<np;i++){
		for (j=0;j<3;j++)	
			fscanf(original,"%f", &p[i][j]); 
		for(j=0;j<3;j++)
			fscanf(original,"%f", &v[j]);
	}
}

else{
	for(i=0;i<np;i++){
		for (j=0;j<3;j++)	
			fscanf(original,"%f", &p[i][j]); 
	}
}
	
fclose(original);

/*Print the particle information in a binary file*/
//fwrite(&np, sizeof(int), 1, output);
//fwrite(&L, sizeof(float), 1, output);
for(i=0;i<np;i++)
	for(j=0;j<3;j++)
		fwrite(&p[i][j], sizeof(float), 1, output);

fclose(output);
return 0;
}
