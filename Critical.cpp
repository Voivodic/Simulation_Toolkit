#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <vector>
#include <cassert>

#define PI 3.141592653589793

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;

typedef Delaunay::Cell_handle                            Cell_handle;
typedef Delaunay::Vertex_handle                          Vertex_handle;

float dist2(Point x, Point y){
	float resp;

	resp = (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]);
	
	return sqrt(resp);
}

int Is_Inside(Point x, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){
	int resp;

	if(x[0]>xmin && x[0]<xmax && x[1]>ymin && x[1]<ymax && x[2]>zmin && x[2]<zmax)
		resp = 1;
	else
		resp = 0;

	return resp;
}

int main(int argc,char *argv[])
{
FILE *part, *crit;
char partfile[100], critfile[100];
int i, j, k, cont, N_cores, np, nh, nv, nptot, flag, ix, iy, iz;
float L, n_m, Dh, Dv, x[3], dist, distt, max, min, Den, xmin, xmax, ymin, ymax, zmin, zmax;
std::vector<Point> P;

/*Check all inputs*/
if (argc != 7){
	printf("Wrong number of arguments.\n");
	printf("arg1: The preffix of the input file\n");
	printf("arg2: Slice's number in the x direction\n");
	printf("arg3: Slice's number in the y direction\n");
	printf("arg4: Slice's number in the z direction\n");
	printf("arg5: Halo's overdensity.\n");
	printf("arg6: Void's overdensity.\n");
	exit(0);
}

/*Read the parameters given by the user*/
ix = atoi(argv[2]);		/*Slice number in the x direction*/
iy = atoi(argv[3]);		/*Slice number in the y direction*/
iz = atoi(argv[4]);		/*Slice number in the z direction*/
Dh = atof(argv[5]);		/*Halo's contrast density*/
Dv = atof(argv[6]);		/*Void's contrast density*/

/*Get the name of all files*/
sprintf(partfile,"%s_%d_%d_%d", argv[1], ix, iy, iz);
sprintf(critfile,"%s_%d_%d_%d.crit", argv[1], ix, iy, iz);

/*Open the particles file*/
part = fopen(partfile,"rb");
if (part == NULL) {
	printf("Unable to open %s\n",partfile);
	exit(0);
}

/*Read the informations in the particle file*/
fread(&nptot, sizeof(int), 1, part);	
fread(&np, sizeof(int), 1, part);	
fread(&L, sizeof(float), 1, part);

fread(&xmin, sizeof(float), 1, part);
fread(&xmax, sizeof(float), 1, part);
fread(&ymin, sizeof(float), 1, part);
fread(&ymax, sizeof(float), 1, part);
fread(&zmin, sizeof(float), 1, part);
fread(&zmax, sizeof(float), 1, part);

/*Some parameters of simulation and grid*/
n_m = np/(L*L*L);			/*The mean density in (Mpc/h)^-3*/

printf("\nDoing slice (%d, %d, %d) of catalogue %s\n", ix, iy, iz, argv[1]);
printf("np = %d, L = %f, n_m = %f and nptot = %d\n", np, L, n_m, nptot);

/*Read the position of all particles*/
printf("Reading the particle catalogue\n");
for (i=0;i<nptot;i++){
	for(j=0;j<3;j++)	
		fread(&x[j], sizeof(float), 1, part);

	P.push_back(Point(x[0], x[1], x[2]));
}
fclose(part);

/*Building the Delaunay triangulation*/
printf("Computting the Delaunay triangulation\n");
Delaunay T(P.begin(), P.end());

/*Open the output file*/
crit = fopen(critfile, "wb");
if (crit == NULL) {
	printf("Unable to open %s\n", critfile);
	exit(0);
}

/*Put the number of halos and voids in the file*/
nh = 0;
nv = 0;
fwrite(&nh, sizeof(int), 1, crit);
fwrite(&nv, sizeof(int), 1, crit);

/*Iterate over finite cells*/
Delaunay::Finite_cells_iterator cit;
Point Q, Q2, Z;
Vertex_handle vh;
Cell_handle ch;

cont = 0;
max = 0.0;
printf("Finding the void centers\n");
for(cit = T.finite_cells_begin(); cit != T.finite_cells_end(); cit++){
	/*Finding the dual point to this cell*/
	Q = T.dual(cit);
	cont++;
	if(Is_Inside(Q, xmin, xmax, ymin, ymax, zmin, zmax) == 0)
		continue;
	
	/*Computting the distance to some particle in the cell*/
	vh = cit->vertex(0);
	Z = vh->point();
	dist = dist2(Q, Z);

	/*Computting the density*/
	Den = 3.0/(n_m*PI*dist*dist*dist);

	/*Save the void centers*/	
	if(Den < Dv){
		if(dist > max) max = dist;

		for(i=0;i<3;i++){
			x[i] = (float) Q[i];
			fwrite(&x[i], sizeof(float), 1, crit);
		}
		fwrite(&Den, sizeof(float), 1, crit);
			
		nv ++;
	}
} 

/*Iterate over finite vertices*/
Delaunay::Finite_vertices_iterator verit;
//Delaunay::Finite_vertices_iterator it;

int num;
min = 1.0;
printf("Finding the halo centers\n");
for(verit = T.finite_vertices_begin(); verit != T.finite_vertices_end(); verit++){
	
	std::vector<Vertex_handle> verh;

	/*Check if this particle is inside the box*/
	Q = verit->point();
	if(Is_Inside(Q, xmin, xmax, ymin, ymax, zmin, zmax) == 0)
		continue;

	T.finite_adjacent_vertices (verit, std::back_inserter(verh));

	/*Find the far neighbour particle*/
	num = 0;
	distt = 0.0;
	for(i=0;i<verh.size();i++){
  		Q2 = verh[i]->point();
   		dist = dist2(Q, Q2);
		if(dist > distt)
			distt = dist;
		num ++;
	}

	Den = num/(4.0/3.0*n_m*PI*distt*distt*distt);

	//printf("Size = %d, num = %d, dist = %f and den = %f\n", verh.size(), num, distt, Den);

	/*Save the halo centers*/
	if(Den > Dh){
		if(distt < min) min = distt;

		for(i=0;i<3;i++){
			x[i] = (float) Q[i];
			fwrite(&x[i], sizeof(float), 1, crit);
		}
		fwrite(&Den, sizeof(float), 1, crit);
		
		nh ++;
	}
	verh.clear();
}

/*Update the number of halos and voids*/
rewind(crit);
fwrite(&nh, sizeof(int), 1, crit);
fwrite(&nv, sizeof(int), 1, crit);
fclose(crit);

if(cont == (int) T.number_of_finite_cells())
	printf("Finite cells = %d\n", cont);
else
	printf("Some problem, cont = %d and finite cells = %d\n", cont, (int) T.number_of_finite_cells());

printf("Maximun distance: %f\n", max);
printf("Minimum distance: %f\n", min);
printf("Number of halos: %d\nNumber of Voids: %d\n", nh, nv);

return 0;
}

