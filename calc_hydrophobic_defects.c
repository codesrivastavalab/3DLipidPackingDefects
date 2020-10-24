/*=================================================================================================================
 ********** Extracting the hydrophobic defects or lipid packing defects on a membrane monolayer surface ***********

 ********************************************* 23rd October 2020 **************************************************

 * Author: Madhusmita Tripathy 
 * 		Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
 *		Indian Institute of Science, Bangalore - 560012
 * 	email: madhu.cfl@gmail.com
 *
 *
 * Note: This code identifies the lattice grid points that constitute the hydrophobic defects in a lipid membrane  *
 * As input, one has to povide the box dimensions, the selection of lipid atoms from one leaflet and the reference *
 * atoms over the trajectory. In simulations, where the membrane CoM is not fixed and evolves over the trajectory, *
 * one has to provide the box centers over the trajectory as well. 
 * 
 * The output is the defect grid points over the trajectory. For checks, one can also output the lipid selection   *
 * and referece atom data.
 *
 * Here, we analyze a lipid membrane with 5 kinds of atoms, more can be added as needed.
 *
 * Cite: M. Tripathy, Subashini T., A. Srivastava, JCTC 2020
 *================================================================================================================*/
 # include <stdio.h>
 # include <math.h>
 # include <stdlib.h>
 # include <string.h>

 # define dim			3		/* Dimension */	

 /* all distances in Angstroms */
 # define boxz			50.0		/* Simulation box size along z-direction: should be big enough to contain all the lipids */ 	

 # define resolution		1.0		/* Lattice resolution */
 # define size_probe		1.4		/* Probe radius: radius of a water molecule */

 # define max_lat_points	500		/* maximum no of lattice voids */
 
 # define max_lines		50000		/* maximum entries in the atom selection file */

 # define size_C		1.7		/* vdw radius of carbon */
 # define size_H		1.2		/* hydrogen */
 # define size_O		1.52		/* oxygen */
 # define size_P		1.8		/* phosphorus */
 # define size_N		1.55		/* nitrogen */

 FILE *fp, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp11, *fp12, *fp_chk;
//================================================================================================================//
//***************************************************** main programme *******************************************//
//===============================================================================================================//
 int main()
 {
 int max_pixel, lat_pixel, pix_i, pix_j, pix_k, ***lat_pix_value, lo[4], hi[4];
 double lat_px[max_lat_points+1], lat_py[max_lat_points+1], lat_pz[max_lat_points+1];
 int lat_pt_counter;
 int dum, pol, bead, i, j, k, l, m, n, o, p, t, count;
 double x, y, z, x1, y1, z1, distance_x, distance_y, distance_z, distance;

 char ch, file[300], file1[300], file2[300], file3[300], file11[300], file12[300], str[10];
 int num_atoms, tp, type[max_lines], bx_count, num_ref_atoms, type_ref[max_lines];
 double bx, by, bz, x_lo, x_hi, y_lo, y_hi, boxx[max_lines], boxy[max_lines];
 double cx, cy, cz, cenx[max_lines], ceny[max_lines];
 double xu[max_lines], yu[max_lines], zu[max_lines];  
 double x_ref[max_lines], y_ref[max_lines], z_ref[max_lines];  
 int lat_pixel_x, lat_pixel_y, lat_pixel_z, z_max[max_lat_points+1][max_lat_points+1], ***defect_lat_pix_value;
 double radius[6], rad, distz;
 int ii, jj, kk;
 int run, start_run=1, max_run=3;

 double z_cutoff = 40.0;		/* The cutoff values used to select lipid atoms from one monolayer */

 lat_pixel_z = (int)(boxz/resolution);			

 /* vdw + probe radius for all atoms */
 radius[1] = size_C + size_probe;	/* Later, twice this value will be the extension of the local grid around each atom */
 radius[2] = size_H + size_probe;
 radius[3] = size_O + size_probe;
 radius[4] = size_P + size_probe;
 radius[5] = size_N + size_probe;

 
 /* Read in box dimensions */
 fp = fopen("input_files/box_dimension.dat", "r");
 i = 1;
 while((!feof(fp)))			
 {
	fscanf(fp,"%lf %lf %lf", &bx, &by, &bz);
	boxx[i] = bx;
	boxy[i] = by;
	i++;
 }
 fclose(fp);

 /* Read in box centers */
 fp = fopen("input_files/box_center.dat", "r");		/* Important when the bilayer CoM is not fixed */
 i = 1;							/* Parrticularly, in NPT simulation, where the system volume fluctuates */
 while((!feof(fp)))			
 {
	fscanf(fp,"%lf %lf %lf", &cx, &cy, &cz);
	cenx[i] = cx;
	ceny[i] = cy;
	i++;
 }
 fclose(fp);


 fp_chk = fopen("output_files/ref_atom_chk.dat", "w");

 system("mkdir all_coord");

 for(run=start_run; run<=max_run; run++)
 {
	x_hi = cenx[run] + (0.5*boxx[run]);		/* Estimate the simulation box extensions in X and Y directions */
	x_lo = cenx[run] - (0.5*boxx[run]);
	y_hi = ceny[run] + (0.5*boxy[run]);
	y_lo = ceny[run] - (0.5*boxy[run]);

	sprintf(file, "input_files/monolayer_atoms_%d.dat", run);	/* Read in lipid selection data */
 	fp = fopen(file, "r");	
 	if(fp == NULL)
 	{
		printf("cannot open file\n");
 	}
 	printf("Reading file\n");
  
 	count = 0;
 	while((!feof(fp)))			
 	{
		count++;

		fscanf(fp," %c %lf %lf %lf", &ch, &x, &y, &z); 
		if (ch == 'C')				/* Identify the atoms, to fix their vdW radius */
		{
			tp = 1;
		}
		else if (ch == 'H')
		{
			tp = 2;
		}
		else if (ch == 'O')
		{
			tp = 3;
		}
		else if (ch == 'P')
		{
			tp = 4;
		}
		else if (ch == 'N')
		{
			tp = 5;
		}
		else
		{
			printf("\nUnidentified atom type: %c\n", ch);
			exit(1);
		}

		type[count] = tp;

		if(x < x_lo)				/* Apply PBC along X and Y directions */
		{
			x = x + boxx[run];
		}
		else if (x > x_hi)
		{
			x = x - boxx[run];
		}
		if(y < y_lo)
		{
			y = y + boxy[run];
		}
		else if (y > y_hi)
		{
			y = y - boxy[run];
		}

		xu[count] = x - x_lo;			/* Move the box to the positive quadrant */
		yu[count] = y - y_lo;			/* with origin at one corner of the box */
		zu[count] = z - z_cutoff;		
 	}
 	fclose(fp);
 
 	num_atoms = count - 1;
 	printf("Number of atoms: %d\n", num_atoms);

	if(run%100 == 0)				/* Coordinates of the shifted lipid selection */
	{						/* Useful as a check, or to generate images */
		sprintf(file11, "all_coord/all_coord_%d.xyz", run);
 		fp11 = fopen(file11, "w");	

 		for(i=1; i<=num_atoms ; i++)
 		{
			fprintf(fp11, "%d\t %lf\t %lf\t %lf\n", type[i], xu[i], yu[i], zu[i]);
 		} 
 		fclose(fp11);  
	} 

 	/* Lattice pixels */
 	lat_pixel_x = (int)(boxx[run]/resolution);	/* Number of lattice points along the X and Y directions */		
 	lat_pixel_y = (int)(boxy[run]/resolution);			
 
 	printf("Grid: %d\t %d\t %d\n", lat_pixel_x, lat_pixel_y, lat_pixel_z);

 	lat_pix_value = (int ***)malloc((lat_pixel_x + 1)*sizeof(int**));	/* Allocate memory */
 	for (i=0; i<=lat_pixel_x; i++) 
 	{
		lat_pix_value[i] = (int **) malloc((lat_pixel_y + 1)*sizeof(int *));
		for (j=0; j<=lat_pixel_y; j++) 
		{
			lat_pix_value[i][j] = (int *)malloc((lat_pixel_z + 1)*sizeof(int));
		      	if(lat_pix_value[i][j] == NULL) 
		     	{
		     	 	free(lat_pix_value[i][j]);
		            	printf("cannot allocate memory for pix_value[%d][%d] \n", i, j);
		            	exit(1);
		      	} 
		}
 	}	 

	defect_lat_pix_value = (int ***)malloc((lat_pixel_x + 1)*sizeof(int**));	/* Allocate memory for the defect grids */
 	for (i=0; i<=lat_pixel_x; i++) 
 	{
		defect_lat_pix_value[i] = (int **) malloc((lat_pixel_y + 1)*sizeof(int *));
		for (j=0; j<=lat_pixel_y; j++) 
		{
			defect_lat_pix_value[i][j] = (int *)malloc((lat_pixel_z + 1)*sizeof(int));
		      	if(defect_lat_pix_value[i][j] == NULL) 
		     	{
		     	 	free(defect_lat_pix_value[i][j]);
		            	printf("cannot allocate memory for pix_value[%d][%d] \n", i, j);
		            	exit(1);
		      	} 
		}
 	} 

 	/* Initialize the lattice pixels: all lattice grids to 1 and defect grids to 0 */
 	for(i=0; i<=lat_pixel_x; i++)
 	{
		for(j=0; j<=lat_pixel_y; j++)
		{
			for(k=0; k<=lat_pixel_z; k++)
			{
				lat_pix_value[i][j][k] = 1;	
				defect_lat_pix_value[i][j][k] = 0;	
			}
		}
 	}
 
	/* Estimate the free volume of the lipid selection: grids in the box which are not occupied by lipid atoms */
	for(n=1; n<=num_atoms; n++)
 	{
		t = type[n];				/* For each atom in the selection, assign a local grid around it */
							/* with dimesion 2X(vdW+probe radius) in all 3 directions */
		lo[1] = (int)((xu[n]-radius[t])/resolution);	hi[1] = (int)((xu[n]+radius[t])/resolution);
		lo[2] = (int)((yu[n]-radius[t])/resolution);	hi[2] = (int)((yu[n]+radius[t])/resolution);
		lo[3] = (int)((zu[n]-radius[t])/resolution);	hi[3] = (int)((zu[n]+radius[t])/resolution);	

		if(lo[1] < 0)				/* These local grid ponints lie between 0 and maximum value of the large grid */
		{
			lo[1] = 0;
		}
		if(hi[1] > lat_pixel_x)
		{
			hi[1] = lat_pixel_x;
		}
		if(lo[2] < 0)
		{
			lo[2] = 0;
		}
		if(hi[2] > lat_pixel_y)
		{
			hi[2] = lat_pixel_y;
		}
		if(lo[3] < 0)
		{
			lo[3] = 0;
		}
		if(hi[3] > lat_pixel_z)
		{
			hi[3] = lat_pixel_z;
		}

		for(i=lo[1]; i<=hi[1]; i++)	/* Check for overlap between the atom and its surrounding local grid points */
		{		
			x1 = i*resolution;
	
			for(j=lo[2]; j<=hi[2]; j++)
			{
				y1 = j*resolution;

				for(k=lo[3]; k<=hi[3]; k++)
				{
					z1 = k*resolution;

					distance_x = x1 - xu[n];	
					distance_y = y1 - yu[n];	
					distance_z = z1 - zu[n];

					distance = sqrt(distance_x*distance_x + distance_y*distance_y + distance_z*distance_z);

					if(distance <= radius[t])  /* Modify the lattice pixels if overlapping */	
					{
						pix_i = i;	pix_j = j;	pix_k = k;
					
						lat_pix_value[pix_i][pix_j][pix_k] = 0; 
					} 		
				}
			}
		}
 	} 


	/* Read in reference atom data */
	sprintf(file1, "input_files/ref_atoms_%d.dat", run);
 	fp1 = fopen(file1, "r");	
	if(fp1 == NULL)
 	{
		printf("cannot open file\n");
 	}

 	printf("Reading file\n");
  
 	count = 0;
 	while((!feof(fp1)))			
 	{
		count++;

		fscanf(fp1," %s %lf %lf %lf", str, &x, &y, &z); 

		if(x < x_lo)					/* Apply PBC along X and Y directions */			
		{
			x = x + boxx[run];
		}
		else if (x > x_hi)
		{
			x = x - boxx[run];
		}
		if(y < y_lo)
		{
			y = y + boxy[run];
		}
		else if (y > y_hi)
		{
			y = y - boxy[run];
		}

		x_ref[count] = x - x_lo;			/* Move the coordinates to the positive quadrant */
		y_ref[count] = y - y_lo;
		z_ref[count] = z - z_cutoff;		
 	}
 	fclose(fp1);
 
	num_ref_atoms = count - 1;
 	printf("Number of reference atoms at setp %d: %d\n", run, num_ref_atoms);

	fprintf(fp_chk, "%d\t %d\n", run, num_ref_atoms);	/* Number of reference atoms over the trajectory */
								/* Useful in case where there are lipid flippings */
 	
	rad = radius[1] + (3.0 * size_probe);			/* Extension of a defect pocket around a reference atom */

	/* Identify the defect pockets: grids in the free volume which are near the lipid tails */
 	for(n=1; n<=num_ref_atoms; n++)
 	{
		/* Fix the local grid again */
		lo[1] = (int)((x_ref[n]-rad)/resolution);	hi[1] = (int)((x_ref[n]+rad)/resolution);
		lo[2] = (int)((y_ref[n]-rad)/resolution);	hi[2] = (int)((y_ref[n]+rad)/resolution);
		lo[3] = (int)((z_ref[n]-rad)/resolution);	hi[3] = (int)((z_ref[n]+rad)/resolution);	

		if(lo[1] < 0)			/* These local grid ponints lie between 0 and maximum value of the large grid */
		{
			lo[1] = 0;
		}
		if(hi[1] > lat_pixel_x)
		{
			hi[1] = lat_pixel_x;
		}
		if(lo[2] < 0)
		{
			lo[2] = 0;
		}
		if(hi[2] > lat_pixel_y)
		{
			hi[2] = lat_pixel_y;
		}
		if(lo[3] < 0)
		{
			lo[3] = 0;
		}
		if(hi[3] > lat_pixel_z)
		{
			hi[3] = lat_pixel_z;
		}
	
		for(i=lo[1]; i<=hi[1]; i++)		/* Check for overlap between the refernce atom and the free volume grid points */
		{		
			x1 = i*resolution;

			for(j=lo[2]; j<=hi[2]; j++)
			{
				y1 = j*resolution;
	
				for(k=lo[3]; k<=hi[3]; k++)
				{
					z1 = k*resolution;

					distance_x = x1 - x_ref[n];	
					distance_y = y1 - y_ref[n];	

					distance = sqrt(distance_x*distance_x + distance_y*distance_y);		
				
				/* If within a certain radius around the reference atom in XY plane and below the Z-coordinate  */
					if(distance < rad && z1 < z_ref[n] && lat_pix_value[i][j][k] == 1) 	
					{
						defect_lat_pix_value[i][j][k] = 1; 
					} 		
				}
			}
		}
 	}	

	sprintf(file3, "output_files/defects_%d.xyz", run);
 	fp3 = fopen(file3, "w");	
 	lat_pt_counter = 0;
 	for(i=0; i<=lat_pixel_x; i++)
 	{
		for(j=0; j<=lat_pixel_y; j++)
		{
			for(k=0; k<=lat_pixel_z; k++)
			{
				if(defect_lat_pix_value[i][j][k]==1)
				{
					lat_pt_counter++;				
					fprintf(fp3, "1\t %.2lf\t %.2lf\t %.2lf\n", i*resolution, j*resolution, k*resolution); 
				}
			}
		}
 	}
 	fclose(fp3); 
 	printf("Number of defect grid points: %d\n", lat_pt_counter); 
 }
 fclose(fp_chk);
 
 for(i=0; i<=lat_pixel_x; i++)							/* Free pointers */
 {
	for(j=0; j<=lat_pixel_y; j++)			
 	{
		free(lat_pix_value[i][j]);
		free(defect_lat_pix_value[i][j]);
	}
 }


}


