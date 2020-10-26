#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#define PI 3.141592653589793

/*Evaluate the ciclic sum of x and y*/
int mod(int x, int y, int nd){
	int resp;

	resp = x + y;

	if(resp<0) resp += nd;
	else if(resp>=nd) resp -= nd;

	return resp;
}

/*Define the cyclic sum for floats*/
double cysumf(double x, double y, double L){
	double resp;

	resp = x + y;
	if(resp>=L)	resp -= L;
	if(resp<0)	resp += L;

	return resp;
}

/*Give a indice for the partricle*/
void ind(double x[], int xt[], double Ld, int nd){
	xt[0] = floor(x[0]/Ld);
	xt[1] = floor(x[1]/Ld);
	xt[2] = floor(x[2]/Ld);
	if(xt[0]==nd)	xt[0] -=1;
	if(xt[1]==nd)	xt[1] -=1;
	if(xt[2]==nd)	xt[2] -=1;
}

/*Give the density to each grid using the NGP density assignment*/
void NGP(double *grid, double *pos, int nd, double L){
	double Ld =  L/nd;
	int xt[3];
	
	ind(pos, xt, Ld, nd);
	
	grid[xt[0]*nd*nd + xt[1]*nd + xt[2]] += 1.0;
}

/*Give the density to each grid using the CIC density assignment*/
void CIC(double *grid, double *pos, int nd, double L, double npart){
	double Ld =  L/nd;
	double dx[3], t[3];
	int xt[3], sign[3], i, j, k;
	
	ind(pos, xt, Ld, nd);
	
	for(i=0;i<3;i++){
		dx[i] = pos[i]/Ld - (xt[i] + 1.0/2.0);	

		if(dx[i]>=0.0){ 
			sign[i] = 1;
			t[i] = 1.0 - dx[i];
		}
		else{
			sign[i] = -1;
			dx[i] = -dx[i];
			t[i] = 1.0 - dx[i];
		}
	}

	grid[xt[0]*nd*nd + xt[1]*nd + xt[2]] += npart*t[0]*t[1]*t[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + xt[1]*nd + xt[2]] += npart*dx[0]*t[1]*t[2];
	grid[xt[0]*nd*nd + mod(xt[1],sign[1],nd)*nd + xt[2]] += npart*t[0]*dx[1]*t[2];
	grid[xt[0]*nd*nd + xt[1]*nd + mod(xt[2],sign[2],nd)] += npart*t[0]*t[1]*dx[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + mod(xt[1],sign[1],nd)*nd + xt[2]] += npart*dx[0]*dx[1]*t[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + xt[1]*nd + mod(xt[2],sign[2],nd)] += npart*dx[0]*t[1]*dx[2];
	grid[xt[0]*nd*nd + mod(xt[1],sign[1],nd)*nd + mod(xt[2],sign[2],nd)] += npart*t[0]*dx[1]*dx[2];
	grid[mod(xt[0],sign[0],nd)*nd*nd + mod(xt[1],sign[1],nd)*nd + mod(xt[2],sign[2],nd)] += npart*dx[0]*dx[1]*dx[2];
}

/*Define de sinc function*/
double sinc(double x){
	double resp;
	if(x!=0.0) resp = sin(x)/x;
	else resp = 1.0;	

	return resp;
}

/*Define window function*/
double W(double k1, double k2, double k3, double Lb){
	double resp;
	resp = sinc(k1*Lb/2.0)*sinc(k2*Lb/2.0)*sinc(k3*Lb/2.0);

	return pow(resp, 2);
}

/*Define the bin for the mode*/
int Indice(double k, double kmin, double dk, int LOG){
	int resp;
	double tmp;

	if(LOG == 0){
		tmp = (k - kmin)/dk;

		resp = floor(tmp);
		if(tmp - (double) resp >= 0.5)
			resp += 1;
	}

	else{
		tmp = (log10(k) - log10(kmin))/dk;

		resp = floor(tmp);
		if(tmp - (double) resp >= 0.5)
			resp += 1;
	}

	return resp;
}

int main(int argc,char *argv[])
{
FILE *first, *second, *out, *grid;
char firstfile[100], secondfile[100], outfile[100], gridfile[100];
int i, j, k, tmp;
int kind_PS, np1, np2, trashi, cont, nd, ng, N, Nk, *contk, ind, LOG, npart;
double L, X[3], trashf, x, y, z, R, errR, den, R_min, R_max, *grid1, *grid2, *Phh, *Phm_r, *Phm_i, *K, *Kmean, n_bar1, n_bar2, kn, Norm, kmin, kmax, dlogk, dk, kx, ky, kz, kmod, Lb;
float tmpf;
fftw_complex *out1, *out2, *gridk;
fftw_plan p1, p2;

if (argc != 10){
	printf("Wrong number of arguments.\n");
	printf("arg1: Name of the first file.\n");
	printf("arg2: Name of the second file.\n");
	printf("arg3: Name of the output file.\n");
	printf("arg4: Box size.\n");
	printf("arg5: What is the minimum radius?\n");	
	printf("arg6: What is the maximum radius?\n");	
	printf("arg7: You want to compute the matter (0), matter-halo/void (1) or just the auto halo/void (2)  spectrum?\n");
	printf("arg8: Number of division by dimension in the grid.\n");
	printf("arg9: Do a linear (0) or log (1) binning?\n\n");
	exit(0);
}

/*Get the name of all files*/
sprintf(firstfile,"%s", argv[1]);
sprintf(secondfile,"%s", argv[2]);
sprintf(outfile,"%s", argv[3]);

/*Some parameters of simulation and grid*/
R_min = atof(argv[5]);		/*Minimum radius of halos/voids*/
R_max = atof(argv[6]);		/*Maximum radius of halos/voids*/
kind_PS = atoi(argv[7]);	/*Kind of measurement*/
nd = atoi(argv[8]);		/*Number of division by dimensions in the grid*/
L = atof(argv[4]);		/*Needed in the auto power spectrum calculaitons*/
LOG = atoi(argv[9]);		/*Kind of binning in the final power spectrum*/
ng = nd*nd*nd;			/*Total number of grid cells*/
Lb = L/nd;			/*Size of each cell*/
Norm = pow(L/(nd*nd), 3.0/2.0);	/*Correct normalization of FFTWs*/
kn = 2.0*PI/L;
Nk = (int) nd/2;

/*Check the kind of power spectrum*/
if(kind_PS != 0 && kind_PS !=1 && kind_PS !=2){
	printf("You need to say what is the kind of power spectrums that you want to calculate! Matter (0),  Halo/void (1) or auto (2).?\n");
	exit(0);
}

/*Check the kind of binning*/
if(LOG != 0 && LOG != 1){
	printf("You need to put a valid binning! Linear (0) or log (1).?\n");
	exit(0);
}

/*Alloc the grids*/
grid1 = (double *)malloc(ng*sizeof(double));
grid2 = (double *)malloc(ng*sizeof(double));

/*Alloc teh FFTW stuff*/
out1 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
out2 = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));
gridk = (fftw_complex*) fftw_malloc(nd*nd*(nd/2+1)*sizeof(fftw_complex));

p1 = fftw_plan_dft_r2c_3d(nd, nd, nd, grid1, out1, FFTW_ESTIMATE); 
p2 = fftw_plan_dft_r2c_3d(nd, nd, nd, grid2, out2, FFTW_ESTIMATE); 

/*Compute the density grid of teh first (matter) file*/
if(kind_PS == 0 || kind_PS == 1){

/*Initialize the grid quantities*/
for(i=0;i<ng;i++){
	grid1[i] = 0.0;
	grid2[i] = 0.0;
}

/*Open the first file*/
first = fopen(firstfile,"rb");
if (first == NULL) {
	printf("Unable to open %s\n",firstfile);
	exit(0);
}

/*Read the total number of particles and the box size in the matter calogue*/
fscanf(first, "%d", &np1);

/*Save the information of the first file*/
for(i=0;i<np1;i++){	
	for(j=0;j<3;j++){	
		fread(&tmpf, sizeof(float), 1, first);
		X[j] = (double) tmpf;
	}

	CIC(grid1, X, nd, L, 1);
	for(j=0;j<3;j++)
		X[j] = cysumf(X[j], -Lb/2.0, L);
	CIC(grid2, X, nd, L, 1);
}
fclose(first);

/*Compute the overdensity (delta)*/
n_bar1 = ((double) np1)/(nd*nd*nd);
for(i=0;i<ng;i++){
	grid1[i] =  grid1[i]/n_bar1 - 1.0;
	grid2[i] =  grid2[i]/n_bar1 - 1.0;
}

kx = 0.0;
ky = 0.0;
for(i=0;i<ng;i++){
	kx += grid1[i];
	ky += grid1[i]*grid1[i];
}
kx = kx/ng;
ky = sqrt(ky)/ng - kx*kx;
printf("Mean = %lf and std = %lf\n", kx, ky);

/*Compute the FFT of the grids*/
fftw_execute(p1);
fftw_execute(p2);

/*Take the mean over the two grids*/
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			tmp = nd*(nd/2 + 1)*i + (nd/2 + 1)*j + k;
			gridk[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0/W(kx, ky, kz, Lb)*Norm;
			gridk[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0/W(kx, ky, kz, Lb)*Norm;
		}
	}
}
}

/*Compute the density grid of the second (halo/void) file*/
if(kind_PS == 1 || kind_PS == 2){

/*Open the second file*/
second = fopen(secondfile,"r");
if (second == NULL) {
	printf("Unable to open %s\n", secondfile);
	exit(0);
}

/*Read the total number of particles of the second catalog*/
fscanf(second,"%d", &np2);

/*Initialize the grid quantities*/
for(i=0;i<ng;i++){
	grid1[i] = 0.0;
	grid2[i] = 0.0;
}

/*Read the halo/void catalog*/
cont = 0;
for(i=0;i<np2;i++){
	fscanf(second, "%d %lf %lf %lf %lf %lf %d %lf", &trashi, &x, &y, &z, &R, &errR, &trashi, &den);
	//fscanf(second, "%lf %lf %lf %lf %lf %lf %lf %d", &x, &y, &z, &trashf, &trashf, &trashf, &R, &npart);
	//fscanf(second, "%lf %lf %lf %lf %lf %lf", &X[0], &X[1], &X[2], &trashf, &trashf, &trashf);

	if(R>=R_min && R<R_max){
		npart = 1;

		CIC(grid1, X, nd, L, 1);
		for(j=0;j<3;j++)
			X[j] = cysumf(X[j], -Lb/2.0, L);
		CIC(grid2, X, nd, L, 1);

		cont += npart;
	}
}
fclose(second);

n_bar2 = ((double) cont)/(nd*nd*nd);
printf("n_bar2 = %e\n", n_bar2);
/*Compute the overdensity (delta)*/
for(i=0;i<ng;i++){
	grid1[i] =  grid1[i]/n_bar2 - 1.0;
	grid2[i] =  grid2[i]/n_bar2 - 1.0;
}

kx = 0.0;
ky = 0.0;
for(i=0;i<ng;i++){
	kx += grid1[i];
	ky += grid1[i]*grid1[i];
}
kx = kx/ng;
ky = sqrt(ky)/ng - kx*kx;
printf("Mean = %lf and std = %lf\n", kx, ky);

/*Compute the FFT of the grids*/
fftw_execute(p1);
fftw_execute(p2);

/*Take the mean over the two grids*/
for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd/2+1;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			tmp = nd*(nd/2 + 1)*i + (nd/2 + 1)*j + k;
			out1[tmp][0] = (out1[tmp][0] + out2[tmp][0]*cos((kx + ky + kz)*Lb/2.0) + out2[tmp][1]*sin((kx + ky + kz)*Lb/2.0))/2.0/W(kx, ky, kz, Lb)*Norm;
			out1[tmp][1] = (out1[tmp][1] + out2[tmp][1]*cos((kx + ky + kz)*Lb/2.0) - out2[tmp][0]*sin((kx + ky + kz)*Lb/2.0))/2.0/W(kx, ky, kz, Lb)*Norm;
		}
	}
}
}

/*Free the density grids*/
free(grid1);
free(grid2);
fftw_free(out2);
fftw_free(p1);
fftw_free(p2);

/*Set the main values for the power spectrums*/
K = (double *)malloc(Nk*sizeof(double));
Kmean = (double *)malloc(Nk*sizeof(double));
contk = (int *)malloc(Nk*sizeof(int));
Phh = (double *)malloc(Nk*sizeof(double));

/*Set the final vectors*/
kmin = kn;
kmax = kn*nd/2.0;

if(LOG == 0){
	dk = (kmax - kmin)/(Nk-1);
	for(i=0;i<Nk;i++){
		Phh[i] = 0.0;
		contk[i] = 0;
		K[i] = kmin + i*dk;
		Kmean[i] = 0.0;
	}


}
else{
	dlogk = (log10(kmax) - log10(kmin))/(Nk-1);
	for(i=0;i<Nk;i++){
		Phh[i] = 0.0;
		contk[i] = 0;
		K[i] = pow(10.0, log10(kmin) + i*dlogk);
		Kmean[i] = 0.0;
	}
}

/*Sum the delta square in each k-bin for the Pmm*/
if(kind_PS == 0){

for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			kmod = sqrt(kx*kx + ky*ky + kz*kz);
			if(LOG == 0)	ind = Indice(kmod, kmin, dk, LOG);
			else	ind = Indice(kmod, kmin, dlogk, LOG);

			if(ind < Nk && ind >= 0){
				if(k < nd/2+1) tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				else tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + nd - k;

				Phh[ind] += (gridk[tmp][0]*gridk[tmp][0] + gridk[tmp][1]*gridk[tmp][1]);
				Kmean[ind] += kmod;
				contk[ind] += 1;
			}
		}
	}
}

/*Take the mean*/
for(i=0;i<Nk;i++)
	if(contk[i]>0){
		Phh[i] = Phh[i]/contk[i];
		Kmean[i] = Kmean[i]/contk[i];	
	}

/*Open the output file*/
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/*Save the Pmm in the out file*/
for(i=0;i<Nk;i++)
	fprintf(out, "%lf %lf %lf %d\n", K[i], Kmean[i], Phh[i], contk[i]);
fclose(out);
}

/*Sum the delta square in each k-bin for the Phm and Phh (or Pvm and Pvv)*/
if(kind_PS == 1){

/*Alloc and initiate the array for the cross power spectrum*/
Phm_r = (double *)malloc(Nk*sizeof(double));
Phm_i = (double *)malloc(Nk*sizeof(double));

for(i=0;i<Nk;i++){
	Phm_r[i] = 0.0;
	Phm_i[i] = 0.0;
}

for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			kmod = sqrt(kx*kx + ky*ky + kz*kz);
			if(LOG == 0)	ind = Indice(kmod, kmin, dk, LOG);
			else	ind = Indice(kmod, kmin, dlogk, LOG);

			if(ind < Nk && ind >= 0){
				if(k < nd/2+1) tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				else tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + nd - k;

				Phh[ind] += out1[tmp][0]*out1[tmp][0] + out1[tmp][1]*out1[tmp][1];
				Phm_r[ind] += gridk[tmp][0]*out1[tmp][0] + gridk[tmp][1]*out1[tmp][1];
				Phm_i[ind] += gridk[tmp][1]*out1[tmp][0] - gridk[tmp][0]*out1[tmp][1];
				Kmean[ind] += kmod;
				contk[ind] += 1;
			}
		}
	}
}

/*Take the mean*/
for(i=0;i<Nk;i++)
	if(contk[i]>0){
		Phh[i] = Phh[i]/contk[i];
		Phm_r[i] = Phm_r[i]/contk[i];
		Phm_i[i] = Phm_i[i]/contk[i];	
		Kmean[i] = Kmean[i]/contk[i];	
	}

/*Open the output file*/
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/*Save the Phm and Phh in the out file*/
for(i=0;i<Nk;i++)
	fprintf(out, "%lf %lf %lf %lf %lf %d\n", K[i], Kmean[i], Phm_r[i], Phm_i[i], Phh[i], contk[i]);
fclose(out);
free(Phm_r); free(Phm_i);
}

/*Sum the delta square in each k-bin for the Phh (or Pvv)*/
if(kind_PS == 2){

for(i=0;i<nd;i++){
	if(2*i<nd) kx = i*kn;
	else kx = (i-nd)*kn;

	for(j=0;j<nd;j++){
		if(2*j<nd) ky = j*kn;
		else ky = (j-nd)*kn;

		for(k=0;k<nd;k++){
			if(2*k<nd) kz = k*kn;
			else kz = (k-nd)*kn;

			kmod = sqrt(kx*kx + ky*ky + kz*kz);
			if(LOG == 0)	ind = Indice(kmod, kmin, dk, LOG);
			else	ind = Indice(kmod, kmin, dlogk, LOG);

			if(ind < Nk && ind >= 0){
				if(k < nd/2+1) tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + k;
				else tmp = i*nd*(nd/2 + 1) + j*(nd/2 + 1) + nd - k;
				Phh[ind] += out1[tmp][0]*out1[tmp][0] + out1[tmp][1]*out1[tmp][1];
				Kmean[ind] += kmod;
				contk[ind] += 1;
			}
		}
	}
}

/*Take the mean*/
for(i=0;i<Nk;i++)
	if(contk[i]>0){
		Phh[i] = Phh[i]/contk[i];
		Kmean[i] = Kmean[i]/contk[i];	
	}

/*Open the output file*/
out = fopen(outfile,"w");
if (out == NULL){
	printf("Unable to open %s\n", outfile);
	exit(0);
}

/*Save the Phh in the out file*/
for(i=0;i<Nk;i++)
	fprintf(out, "%lf %lf %lf %d\n", K[i], Kmean[i], Phh[i], contk[i]);
fclose(out);
}

/*Free the memory*/
free(K);
free(Phh);
fftw_free(out1);
if(kind_PS < 2)
	fftw_free(gridk);

return 0;
}
