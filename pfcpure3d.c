#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <fftw3-mpi.h>

//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#define Pi (2.*acos(0.))

//const gsl_rng *gBaseRand;       /* global rand number generator */

//domain parameters
char run[3];			//label of current run
double atomsx,atomsy,atomsz;	//number of atoms in x,y and z
ptrdiff_t Nx,Ny,Nz;		//number of grid points in x,y and z
ptrdiff_t Nx2;			//padded number of grid points in x
double spacing;			//lattice spacing
double dx,dy,dz;		//numerical grid spacings
double n0;			//avg density value
char restartrun[3];		//run containing IC state
int itemp;			//thermal noise switch
double ampc,epdt;		//thermal noise amplitude unscaled, scaled

//time parameters
int totalTime, printFreq, eneFreq;	//simulation time, print freq, ene freq
double dt,ene1,ene2,ene3;				//time step, sys avg free energy
int restartFlag,restartTime;		//flag for restarting (set to 1 or 2) and time

//MPI parameters
int argc;
char **argv;
int myid,numprocs;
ptrdiff_t alloc_local, local_n0, local_0_start;

/**********************correlation variables for standard PFC**********************/
//double epsstart,omcut,rk0;	//epsstart=-r, omcut=min C2 value, rk0=C2(k=0)

//**********************correlation variables for XPFC**********************/
double sigmaT,omcut,rk0;	// DebyeWaller Temperature, omcut=min C2 value, rk0=C2(k=0)
double alpha1;			// Widths of 1st C2 peak
double rho1;			// Density of 1st family of planes
double beta1; 			// Number of planes within 1st family of planes

/**********************equation of motion variables****************************/
int icons,iwave;		//1=conserved dynamics,1=wave dynamics
double Mn;			//mobility constant for density (iwave=0)
double alphaw,betaw;		//wave speed, wave damping rate

/**********************BCC IC parameters***************************************/
double amp0,icpower;	//amp0=amplitude, icpower=exponent

double w,u;		//coefficients for cubic and quartic terms in energy
double gamma12, gamma13, gamma23, nc, sigmac;         //Coupling coefficients between layers
double facx,facy,facz;	//fourier factors for scaling lengths in k-space

//*******************Arrays for fields
//real space arrays
double *n1;			//density
double *n1fNL;			//non-linear term for density (gfr)
double *n2;			//density
double *n2fNL;			//non-linear term for density (gfr)
double *n3;			//density
double *n3fNL;			//non-linear term for density (gfr)
double *g1;			//wave term (g=dn/dt)
double *g2;			//wave term (g=dn/dt)
double *g3;			//wave term (g=dn/dt)
double *ran2;			//noise array2
double *ran3;			//noise array3

//complex arrays
fftw_complex *kn1;			//k-space density
fftw_complex *kn1fNL;			//k-space non-linear term
fftw_complex *kn2;			//k-space density
fftw_complex *kn2fNL;			//k-space non-linear term
fftw_complex *kn3;			//k-space density
fftw_complex *kn3fNL;			//k-space non-linear term
fftw_complex *omegak1,*k1arr;		//k-space C2 and k2 arrays
fftw_complex *omegak2,*k2arr;		//k-space C2 and k2 arrays
fftw_complex *omegak3,*k3arr;		//k-space C2 and k2 arrays

//fftw plans
fftw_plan planF_n1, planB_n1, planF_n2, planB_n2, planF_n3, planB_n3;		//forward and backward plan for density
fftw_plan planF_NL_n1, planF_NL_n2, planF_NL_n3;			//non-linear forward transform
fftw_plan planB_NL_n1, planB_NL_n2, planB_NL_n3;			//backward plan for energy
fftw_plan planF_g1, planF_g2, planF_g3;			//forward plan for g



void freeMemory()
{	
	free(n1);
	free(n2);
	free(n3);
	free(n1fNL);
	free(n2fNL);
	free(n3fNL);
	if (iwave==1)
	{
		free(g1);
		free(g2);
		free(g3);
	}
	if (itemp > 0) free(ran2);
	if (itemp > 0) free(ran3);

	fftw_free(kn1);
	fftw_free(kn2);
	fftw_free(kn3);
	fftw_free(kn1fNL);
	fftw_free(kn2fNL);
	fftw_free(kn3fNL);
	fftw_free(omegak1);
	fftw_free(omegak2);
	fftw_free(omegak3);
	fftw_free(k1arr);
	fftw_free(k2arr);
	fftw_free(k3arr);

	fftw_destroy_plan(planF_n1);
	fftw_destroy_plan(planF_n2);
	fftw_destroy_plan(planF_n3);
	fftw_destroy_plan(planB_n1);
	fftw_destroy_plan(planB_n2);
	fftw_destroy_plan(planB_n3);
	fftw_destroy_plan(planF_NL_n1);
	fftw_destroy_plan(planF_NL_n2);
	fftw_destroy_plan(planF_NL_n3);
	fftw_destroy_plan(planB_NL_n1);
	fftw_destroy_plan(planB_NL_n2);
	fftw_destroy_plan(planB_NL_n3);
	if (iwave==1) 
	{
		fftw_destroy_plan(planF_g1);
		fftw_destroy_plan(planF_g2);
		fftw_destroy_plan(planF_g3);
	}
}




void restart(int time,double *Array, char *var, char *filepre)
{
	double *tmp;
	double dum;
	char filename[BUFSIZ];
	MPI_Status status;
    int tag = 0;
    	FILE *fp;

    	ptrdiff_t i,j,k;

    	double * buffer_data = (double*) malloc((local_n0*Ny*Nx)*sizeof(double));

    	if(myid == 0)
    	{
    		printf("var=%s, filepre=%s, time=%d\n",var,filepre,time);
        	fflush(stdout);
    		sprintf(filename,"%s_%s_%d.dat",var,filepre,time);

    		// Read from file
    		fp = fopen(filename,"r");
    		double * temp_data = (double*) malloc((Nx*Ny*Nz)*sizeof(double));
        	float number = 0.;
    		if (fp != NULL)
        	{
            	for(k=0;k<Nz;k++)
            	for(j=0;j<Ny;j++)
            	for(i=0;i<Nx;i++)
            	{
            	    fscanf(fp, "%f ", &number);
            	    temp_data[i+j*Nx+k*Nx*Ny] = number;
            	}
        	}
        	else
        	{
        	    printf("error! file does not exist. \n");
        	}
        	fclose(fp);

        	// Transfer to self
        	for(k=0;k<local_n0;k++)
			for(j=0;j<Ny;j++)
			for(i=0;i<Nx;i++)
				buffer_data[i+j*Nx+k*Nx*Ny] = temp_data[i+j*Nx+k*Nx*Ny];

            // Transfer to slaves
    		for (int ii=1; ii<numprocs; ii++)
        	{
            	int offset, count;
            	// master asks for local information from slave
            	MPI_Recv(&offset, 1, MPI_INT, ii, tag, MPI_COMM_WORLD, &status);
            	MPI_Recv(&count, 1, MPI_INT, ii, tag, MPI_COMM_WORLD, &status);

            	// master sends data
            	MPI_Send(temp_data+offset, count, MPI_DOUBLE, ii, tag, MPI_COMM_WORLD);
        	}
        	free(temp_data);
    	}
    	else
    	{
        	// slaves send local information to master
        	int offset = local_0_start*Ny*Nx;
        	int count = local_n0*Ny*Nx;
        	MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        	MPI_Send(&count, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

        	// slave recv data
        	MPI_Recv(buffer_data, count, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    	}

    	
	for(k=0;k<local_n0;k++)
		for(j=0;j<Ny;j++)
			for(i=0;i<Nx;i++)
					Array[i+j*Nx2+k*Nx2*Ny] = buffer_data[i+j*Nx+k*Nx*Ny];

	free(buffer_data);
	MPI_Barrier(MPI_COMM_WORLD);
}





void output(char *var, int time,double *Array,char *filepre)
{
	char filename[BUFSIZ];
    	FILE *fp;
    int tag = 0;
    MPI_Status status;

	ptrdiff_t i,j,k;
	int kk; 

	double * resize = (double*) malloc((local_n0*Ny*Nx)*sizeof(double));
	for(k=0;k<local_n0;k++)
		for(j=0;j<Ny;j++)
			for(i=0;i<Nx;i++)
				resize[i+Nx*(j+k*Ny)] = Array[i+Nx2*(j+k*Ny)];

	if ( myid == 0 )
	{
		double * buffer = (double*) malloc((Nx*Ny*Nz)*sizeof(double));
		memcpy(buffer, resize, local_n0*Ny*Nx*sizeof(double));
		for (int i=1; i<numprocs; i++)
        {
            int count, offset;
            MPI_Recv(&offset, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&count, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(buffer + offset, count, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }
	

	// Output double precision text files with n and, if iwave=1, g
    	sprintf(filename,"%s_%s_%d.dat",var,filepre,time);
	   		fp = fopen(filename,"w");
		if(fp==NULL)
		{
			printf("Unable to open file for writing\n");
			exit(1);
		}
		for(k=0;k<Nz;k++)
			for(j=0;j<Ny;j++)
				for(i=0;i<Nx;i++)
					fprintf(fp,"%.5g\n",buffer[i+Nx*(j+k*Ny)]);
		free(buffer);
		fclose(fp);
	}
	else
    {
        int offset = local_0_start*Ny*Nx;
        int count = local_n0*Ny*Nx;
        MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(&count, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(resize, count, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
    free(resize);
}




void normalize(double *Array)
{
	ptrdiff_t i,j,k;
	ptrdiff_t index,index1,index2;

	//normalize the transform
	for(k=0;k<local_n0;k++)
	{
		index1 = Nx2*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + Nx2*j;
			for(i=0;i<Nx;i++)
			{
				index = i + index2;
				Array[index] /= (Nx*Ny*Nz);
			}
		}
	}
}




void timeStepn()
{
	ptrdiff_t i,j,k;
	ptrdiff_t index,index1,index2;

	for(k=0;k<local_n0;k++)
	{
		index1 = (Nx/2+1)*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + (Nx/2+1)*j;
			for(i=0;i<(Nx/2+1);i++)
			{
				index = i + index2;
				//STRONG SI TIME STEP
				kn1[index][0] = omegak1[index][0] * ( kn1[index][0] + k1arr[index][0]*(kn1fNL[index][0]) );
				kn1[index][1] = omegak1[index][1] * ( kn1[index][1] + k1arr[index][0]*(kn1fNL[index][1]) );
				kn2[index][0] = omegak2[index][0] * ( kn2[index][0] + k2arr[index][0]*(kn2fNL[index][0]) );
				kn2[index][1] = omegak2[index][1] * ( kn2[index][1] + k2arr[index][0]*(kn2fNL[index][1]) );
				kn3[index][0] = omegak3[index][0] * ( kn3[index][0] + k3arr[index][0]*(kn3fNL[index][0]) );
				kn3[index][1] = omegak3[index][1] * ( kn3[index][1] + k3arr[index][0]*(kn3fNL[index][1]) );
			}
		}
	}
}





void timeStepnWave()
{
	ptrdiff_t i,j,k;
	ptrdiff_t index,index1,index2;
	double rmult;

	// First part of algorithm
	//SI1
        rmult = dt*(1.-dt*betaw);
	//SI2,EXPLICIT
        //rmult = dt/(1.+dt*betaw);
	for(k=0;k<local_n0;k++)
	{
		index1 = (Nx/2+1)*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + (Nx/2+1)*j;
			for(i=0;i<(Nx/2+1);i++)
			{
				index = i + index2;
				//SI1,SI2
				kn1[index][0] = omegak1[index][0] * ( kn1[index][0] + dt*k1arr[index][0]*kn1fNL[index][0] );
				kn1[index][1] = omegak1[index][1] * ( kn1[index][1] + dt*k1arr[index][1]*kn1fNL[index][1] );
				kn2[index][0] = omegak2[index][0] * ( kn2[index][0] + dt*k2arr[index][0]*kn2fNL[index][0] );
				kn2[index][1] = omegak2[index][1] * ( kn2[index][1] + dt*k2arr[index][1]*kn2fNL[index][1] );
				kn3[index][0] = omegak3[index][0] * ( kn3[index][0] + dt*k3arr[index][0]*kn3fNL[index][0] );
				kn3[index][1] = omegak3[index][1] * ( kn3[index][1] + dt*k3arr[index][1]*kn3fNL[index][1] );
				//EXPLICIT
				//kn[index][0] = (omegak[index][0] + 1.) * kn[index][0] + rmult*karr[index][0]*knfNL[index][0];
				//kn[index][1] = (omegak[index][1] + 1.) * kn[index][1] + rmult*karr[index][1]*knfNL[index][1];
			}
		}
	}
	// Second part of algorithm
	fftw_execute(planF_g1);
	fftw_execute(planF_g2);
	fftw_execute(planF_g3);
	for(k=0;k<local_n0;k++)
	{
		index1 = (Nx/2+1)*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + (Nx/2+1)*j;
			for(i=0;i<(Nx/2+1);i++)
			{
				index = i + index2;
				//SI1,SI2
				kn1[index][0] = kn1[index][0] + rmult*omegak1[index][0]*kn1fNL[index][0];
				kn1[index][1] = kn1[index][1] + rmult*omegak1[index][1]*kn1fNL[index][1];
				kn2[index][0] = kn2[index][0] + rmult*omegak2[index][0]*kn2fNL[index][0];
				kn2[index][1] = kn2[index][1] + rmult*omegak2[index][1]*kn2fNL[index][1];
				kn3[index][0] = kn3[index][0] + rmult*omegak3[index][0]*kn3fNL[index][0];
				kn3[index][1] = kn3[index][1] + rmult*omegak3[index][1]*kn3fNL[index][1];
				//EXPLICIT
				//kn[index][0] = kn[index][0] + rmult*knfNL[index][0];
				//kn[index][1] = kn[index][1] + rmult*knfNL[index][1];
			}
		}
	}
	// Put old n array into old g array for updating g
	for(k=0;k<local_n0;k++)
	{
		index1 = Nx2*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + Nx2*j;
			for(i=0;i<Nx;i++)
			{
				index = i + index2;
				g1[index]=n1[index];
				g2[index]=n2[index];
				g3[index]=n3[index];
			}
		}
	}
}





void calcNL(int time)
{
	ptrdiff_t i,j,k,h;
	ptrdiff_t index,index1,index2;

        double rarrtop[2*2], sigma = 1.0, rnum, epdt2;
	int iup,idn,index1p,index2p,index2pp,indexip;
        MPI_Status status;

	// Noise term
    /*    if (itemp == 1)
	{

        //CONSERVED ADDITIVE NOISE
		iup = myid+1;
		idn = myid-1;
		if (myid == numprocs-1) iup = 0;
		if (myid == 0) idn = numprocs-1;
		for(k=0;k<local_n0;k++)
		{
			index1 = Nx2*k*Ny;
			for(j=0;j<Ny;j++)
			{
				index2 = index1 + Nx2*j;
				for(i=0;i<Nx;i++)
				{
					index=i + index2;

//      					rnum =gsl_ran_gaussian_ziggurat(gBaseRand,sigma);
//             	       			nfNL[index] = rnum;
//      					rnum =gsl_ran_gaussian_ziggurat(gBaseRand,sigma);
//             	       			ran2[index] = rnum;
//      					rnum =gsl_ran_gaussian_ziggurat(gBaseRand,sigma);
//             	       			ran3[index] = rnum;
				}
			}
		}
            	//SHARE BOTTOM PLANE, 3RD SET OF RANDOM NUMBERS
		for(j=0;j<Ny;j++)
		{
			index2 = Nx2*j;
			for(i=0;i<Nx;i++)
			{
				index = i + index2;
        			MPI_Send(&ran3[index],1,MPI_DOUBLE,idn,i+j*Nx,MPI_COMM_WORLD);
        			MPI_Recv(&rarrtop[i+j*Nx],1,MPI_DOUBLE,iup,i+j*Nx,MPI_COMM_WORLD,&status);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//ADD NOISE
              	epdt2 = epdt/dx;
		for(k=0;k<local_n0-1;k++)
		{
			index1 = Nx2*k*Ny;
			index1p = Nx2*(k+1)*Ny;
			for(j=0;j<Ny-1;j++)
			{
				index2 = index1 + Nx2*j;
				index2p = index1 + Nx2*(j+1);
				index2pp = index1p + Nx2*j;
				for(i=0;i<Nx-1;i++)
				{
				      index = i + index2;
             	       		      n[index]=n[index]+epdt2*(
                               			nfNL[index+1]-nfNL[index]+
                               			ran2[i+index2p]-ran2[index]+
                               			ran3[i+index2pp]-ran3[index]);
				}
			}
		}
		for(k=0;k<local_n0-1;k++)
		{
			index1 = Nx2*k*Ny;
			index1p = Nx2*(k+1)*Ny;
			for(j=Ny-1;j<Ny;j++)
			{
				index2 = index1 + Nx2*j;
				index2p = index1 + Nx2*0;
				index2pp = index1p + Nx2*j;
				for(i=0;i<Nx-1;i++)
				{
				      index = i + index2;
             	       		      n[index]=n[index]+epdt2*(
                               			nfNL[index+1]-nfNL[index]+
                               			ran2[i+index2p]-ran2[index]+
                               			ran3[i+index2pp]-ran3[index]);
				}
			}
		}
		for(k=0;k<local_n0-1;k++)
		{
			index1 = Nx2*k*Ny;
			index1p = Nx2*(k+1)*Ny;
			for(j=0;j<Ny-1;j++)
			{
				index2 = index1 + Nx2*j;
				index2p = index1 + Nx2*(j+1);
				index2pp = index1p + Nx2*j;
				for(i=Nx-1;i<Nx;i++)
				{
				      index = i + index2;
             	       		      n[index]=n[index]+epdt2*(
                              			nfNL[index2]-nfNL[index]+
                              			ran2[i+index2p]-ran2[index]+
                               			ran3[i+index2pp]-ran3[index]);
				}
			}
		}
		for(k=0;k<local_n0-1;k++)
		{
			index1 = Nx2*k*Ny;
			index1p = Nx2*(k+1)*Ny;
			for(j=Ny-1;j<Ny;j++)
			{
				index2 = index1 + Nx2*j;
				index2p = index1 + Nx2*0;
				index2pp = index1p + Nx2*j;
				for(i=Nx-1;i<Nx;i++)
				{
				      index = i + index2;
             	       		      n[index]=n[index]+epdt2*(
                               			nfNL[index2]-nfNL[index]+
                               			ran2[i+index2p]-ran2[index]+
                               			ran3[i+index2pp]-ran3[index]);
				}
			}
		}
		for(k=local_n0-1;k<local_n0;k++)
		{
			index1 = Nx2*k*Ny;
			for(j=0;j<Ny;j++)
			{
				index2 = index1 + Nx2*j;
				index2p = index1 + Nx2*(j+1);
                    		if (j == Ny-1) index2p = index1 + Nx2*0;
				index2pp = Nx*j;
				for(i=0;i<Nx;i++)
				{
				      index = i + index2;
                       		      indexip = 1+index;
                       		      if (i == Nx-1) indexip = index2;
             	       		      n[index]=n[index]+epdt2*(
                               			nfNL[indexip]-nfNL[index]+
                               			ran2[i+index2p]-ran2[index]+
                               			rarrtop[i+index2pp]-ran3[index]);
				}
			}
		}
	}*/

	// Bulk nonlinear terms
	for(k=0;k<local_n0;k++)
	{
		index1 = Nx2*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + Nx2*j;
			for(i=0;i<Nx;i++)
			{
				gamma13 = -0.06*sin(2*Pi*j/Ny);
				index = i + index2;
				double couple1,couple2,couple3,dcouple1,dcouple2,dcouple3;
				couple1 = (1+tanh((n1[index]-nc)/sigmac))/2.;
				dcouple1 = 1./cosh((n1[index]-nc)/sigmac)/cosh((n1[index]-nc)/sigmac)/2./sigmac;
				couple2 = (1+tanh((n2[index]-nc)/sigmac))/2.;
				dcouple2 = 1./cosh((n2[index]-nc)/sigmac)/cosh((n2[index]-nc)/sigmac)/2./sigmac;
				couple3 = (1+tanh((n3[index]-nc)/sigmac))/2.;
				dcouple3 = 1./cosh((n3[index]-nc)/sigmac)/cosh((n3[index]-nc)/sigmac)/2./sigmac;
				n1fNL[index] = n1[index]*n1[index]*( -0.5*w + u*n1[index]/3. ) + dcouple1*(gamma12*couple2+gamma13*couple3); 
				n2fNL[index] = n2[index]*n2[index]*( -0.5*w + u*n2[index]/3. ) + dcouple2*(gamma12*couple1+gamma23*couple3);
				n3fNL[index] = n3[index]*n3[index]*( -0.5*w + u*n3[index]/3. ) + dcouple3*(gamma13*couple1+gamma12*couple2);
			}
		}
	}
	fftw_execute(planF_NL_n1);
	fftw_execute(planF_NL_n2);
	fftw_execute(planF_NL_n3);
}




void calcCorrelations(int time)
{
	ptrdiff_t i,j,k,kk;
	ptrdiff_t index,index1,index2;
	ptrdiff_t zk;
	double kx,ky,kz,k2,rk;

	int Nx2p1 = Nx/2+1;
      	double dtMn = dt*Mn;
      	double dtMn2 = -1./dt/Mn;
      	double dtalph2 = -1./dt/alphaw/alphaw;
      	double dtalphbet = -dt*dt*alphaw*alphaw/(1.+betaw*dt);
      	double Pi32 = 32.*Pi*Pi;

	double omval;
	double qC1,qC2,PreC1,PreC2,DW1,DW2,kC1,kC2;
	double qC3,PreC3,DW3,kC3;

	char filename[BUFSIZ];
    	FILE *fp;

	// XPFC kernel for BCC - 1 peak
	qC1 = 2.*Pi*sqrt(3);		//first mode
	PreC1 = -0.5/(alpha1*alpha1);
	DW1 = exp(-0.5*sigmaT*sigmaT*qC1*qC1/(rho1*beta1) );
 
	for(k=0;k<local_n0;k++)
	{
		zk = k + local_0_start;

		if ( zk < Nz/2 )
			kz = zk*facz;
		else
			kz = (zk-Nz)*facz;

		for(j=0;j<Ny;j++)
		{
			if ( j < Ny/2 )
				ky = j*facy;
			else
				ky = (j-Ny)*facy;

			for(i=0;i<Nx2p1;i++)
			{

			   index = i + Nx2p1*(j+k*Ny);
			   kx = i*facx;
		  	   k2 = kx*kx + ky*ky + kz*kz;
			   if (k2 >= 0.0) rk = sqrt( k2 );
			   else rk = 0.0;

			   //XPFC kernel: omval=1-C2(k)
			   kC1 = DW1*exp( PreC1 * (rk-qC1)*(rk-qC1) );
			   omval = 1.-kC1;
		   
			   // HIGH k CUTTOFF: OTHERWISE XPFC FREE ENERGY VARIES WILDLY WITH dx, DOESNT CONVERGE
		   	   if (k2 > Pi32) omval=(1.-k2*2./Pi32)*(1.-k2*2./Pi32);

	                   if (omval > omcut) omval = omcut;
	                   if (k2 == 0.) rk0=omval;    //For energy calc
	                   if (icons == 1) omval=k2*omval;
	                   if (iwave == 0)
	                   {
				   			k1arr[index][0] = -dtMn*k2;
				   			k1arr[index][1] = -dtMn*k2;
				   			omegak1[index][0] = 1./(1. + dtMn*omval);
				   			omegak1[index][1] = 1./(1. + dtMn*omval);
				   			k2arr[index][0] = -dtMn*k2;
				   			k2arr[index][1] = -dtMn*k2;
				   			omegak2[index][0] = 1./(1. + dtMn*omval);
				   			omegak2[index][1] = 1./(1. + dtMn*omval);
				   			k3arr[index][0] = -dtMn*k2;
				   			k3arr[index][1] = -dtMn*k2;
				   			omegak3[index][0] = 1./(1. + dtMn*omval);
				   			omegak3[index][1] = 1./(1. + dtMn*omval);
	                   }
	                   else
	                   {
				   k1arr[index][0] = k2/dtalph2;
				   k1arr[index][1] = k2/dtalph2;
				   k2arr[index][0] = k2/dtalph2;
				   k2arr[index][1] = k2/dtalph2;
				   k3arr[index][0] = k2/dtalph2;
				   k3arr[index][1] = k2/dtalph2;
				   //SI1
				   omegak1[index][0] = 1./(1. - dt*omval/dtalph2);
				   omegak1[index][1] = 1./(1. - dt*omval/dtalph2);
				   omegak2[index][0] = 1./(1. - dt*omval/dtalph2);
				   omegak2[index][1] = 1./(1. - dt*omval/dtalph2);
				   omegak3[index][0] = 1./(1. - dt*omval/dtalph2);
				   omegak3[index][1] = 1./(1. - dt*omval/dtalph2);
				   //SI2
				   //omegak[index][0] = 1./(1. - dtalphbet*omval);
				   //omegak[index][1] = 1./(1. - dtalphbet*omval);
				   //EXPLICIT
				   //omegak[index][0] = dtalphbet*omval;
				   //omegak[index][1] = dtalphbet*omval;
	                   }

			}
		}
	}	

}




void initialize(int type)
{
	ptrdiff_t i,j,k;
	ptrdiff_t index;
	
	double qx,qy,qz;
    double x,y,z;
    double x1,y1,z1;
  	double sigma = 1.0, rnum;
  	double theta = 0.;

        MPI_Status status;

        qx = 2.*Pi*dx*sqrt(3);
       	qy = 2.*Pi*dy*sqrt(3);
        qz = 2.*Pi*dz*sqrt(3);

	// Input BCC crystal
	for(k=0;k<local_n0;k++)
	{
		z = (double)(k + local_0_start);
                for(j=0;j<Ny;j++)
                {
			y = (double)(j);
			for(i=0;i<Nx;i++)
			{
				x = (double)(i);
				index = i+Nx2*(j+k*Ny);

				if(y>=10 &&y <= Ny/2-10)
				{
					n1[index] = amp0*(cos(qy*(y-4*Pi/3./qy))+cos(-qy*(y-4*Pi/3./qy)*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*(y-4*Pi/3./qy)*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
					n2[index] = amp0*(cos(qy*y)+cos(-qy*y*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*y*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
					n3[index] = amp0*(cos(qy*(y-4*Pi/3./qy))+cos(-qy*(y-4*Pi/3./qy)*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*(y-4*Pi/3./qy)*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
				}
				else if(y >= Ny/2+10 && y<Ny-10)
				{
					y1 = y*cos(theta)+z*sin(theta);
					z1 = -y*sin(theta)+z*cos(theta);
					n1[index] = amp0*(cos(qy*(y1-2*Pi/3./qy))+cos(-qy*(y1-2*Pi/3./qy)*sin(-Pi/6.)+qz*(z1-2*sqrt(3)*Pi/3./qy)*cos(-Pi/6.))+cos(-qy*(y1-2*Pi/3./qy)*sin(Pi/6.)+qz*(z1-2*sqrt(3)*Pi/3./qy)*cos(Pi/6.)))+n0;
					n2[index] = amp0*(cos(qy*y1)+cos(qy*(-y1*sin(-Pi/6.)+z1*cos(-Pi/6.)))+cos(qy*(-y1*sin(Pi/6.)+z1*cos(Pi/6.))))+n0;
					n3[index] = amp0*(cos(qy*(y1-4*Pi/3./qy))+cos(-qy*(y1-4*Pi/3./qy)*sin(-Pi/6.)+qz*z1*cos(-Pi/6.))+cos(-qy*(y1-4*Pi/3./qy)*sin(Pi/6.)+qz*z1*cos(Pi/6.)))+n0;
				}
				else
				{
					n1[index] = n0;
					n2[index] = n0;
					n3[index] = n0;
				}

				n1[index] = amp0*(cos(qy*(y-4*Pi/3./qy))+cos(-qy*(y-4*Pi/3./qy)*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*(y-4*Pi/3./qy)*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
				n2[index] = amp0*(cos(qy*y)+cos(-qy*y*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*y*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
				n3[index] = amp0*(cos(qy*(y-4*Pi/3./qy))+cos(-qy*(y-4*Pi/3./qy)*sin(-Pi/6.)+qz*z*cos(-Pi/6.))+cos(-qy*(y-4*Pi/3./qy)*sin(Pi/6.)+qz*z*cos(Pi/6.)))+n0;
				
				if (iwave==1) 
				{
					g1[index] = 0.0;
					g2[index] = 0.0;
					g3[index] = 0.0;
//      					rnum =gsl_ran_gaussian_ziggurat(gBaseRand,sigma);
//					g[index] = rnum*0.0000001;
				}
			}
		}
	}
}



void setfftwPlans()
{
	/************** n **************/
	planF_n1 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n1, kn1, MPI_COMM_WORLD, FFTW_MEASURE);
	planF_n2 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n2, kn2, MPI_COMM_WORLD, FFTW_MEASURE);
	planF_n3 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n3, kn3, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_n1 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn1, n1, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_n2 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn2, n2, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_n3 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn3, n3, MPI_COMM_WORLD, FFTW_MEASURE);
	planF_NL_n1 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n1fNL, kn1fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	planF_NL_n2 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n2fNL, kn2fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	planF_NL_n3 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, n3fNL, kn3fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_NL_n1 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn1, n1fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_NL_n2 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn2, n2fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	planB_NL_n3 = fftw_mpi_plan_dft_c2r_3d(Nz, Ny, Nx, kn3, n3fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	if (iwave==1)
	{
		planF_g1 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, g1, kn1fNL, MPI_COMM_WORLD, FFTW_MEASURE);
		planF_g2 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, g2, kn2fNL, MPI_COMM_WORLD, FFTW_MEASURE);
		planF_g3 = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, g3, kn3fNL, MPI_COMM_WORLD, FFTW_MEASURE);
	}
}



void allocateArrays()
{
	/****************** n ******************/
	n1 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	n2 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	n3 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	n1fNL = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	n2fNL = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	n3fNL = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	if (iwave==1)
	{
		g1 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
		g2 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
		g3 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	}
	if (itemp > 0) 
	{
		ran2 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
		ran3 = (double *) fftw_malloc( sizeof (double)*(alloc_local+4*Ny*local_n0) );
	}
       
	kn1 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	kn2 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	kn3 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	kn1fNL = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	kn2fNL = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	kn3fNL = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	omegak1 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	omegak2 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	omegak3 = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	k1arr = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	k2arr = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	k3arr = (fftw_complex *) fftw_malloc( sizeof (fftw_complex)*alloc_local );
	printf("%d local_sizeof_real= %zu \n",myid,sizeof (double)*2*alloc_local);
	printf("%d local_sizeof_complex= %zu \n",myid,sizeof (fftw_complex)*alloc_local);
}



void fftwMPIsetup()
{
	//starting mpi daemons
	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	printf("\n myid: %d of numprocs: %d processors has been activated \n",myid,numprocs);

	fftw_mpi_init();
	
	//get the local (for current cpu) array sizes
	alloc_local = fftw_mpi_local_size_3d(Nz, Ny, Nx, MPI_COMM_WORLD,&local_n0, &local_0_start);

	printf(" %d local_n0=%zu local_0_start=%zu\n",myid,local_n0,local_0_start);
}




void domainParams()
{
	//from the lattice spacing and grid spacing, compute the number of grid points
	if (myid==0) printf("dx=%f dy=%f dz=%f numAtomsx=%f numAtomsy=%f numAtomsz=%f\n",dx, dy, dz, atomsx,atomsy,atomsz);

	//BCC
	Nx = (int)(floor(spacing/dx*atomsx+.5));
	Ny = (int)(floor(spacing/dy*atomsy+.5));
	Nz = (int)(floor(spacing/dz*atomsz+.5));

	if (myid==0) printf("number of grid points in x %zu\n",Nx);
	if (myid==0) printf("number of grid points in y %zu\n",Ny);
	if (myid==0) printf("number of grid points in z %zu\n",Nz);
	
	Nx2 = 2*(Nx/2+1);

	//set up fourier scaling factors
	facx = 2.*Pi/(Nx*dx);
	facy = 2.*Pi/(Ny*dy);
	facz = 2.*Pi/(Nz*dz);

	//thermal noise amplitude, dt scaling
	//epdt = ampc/sqrt(dt);	//nonconserved dynamics
	epdt = ampc*sqrt(dt);	//conserved dynamics
}




void inputVariables()
{
	char *line = malloc(BUFSIZ);
	FILE *in;
	in = fopen("pfcpure3d.in","rt");	
   
	if (myid==0) printf("reading input file\n");
	if (in == NULL)
	{
		printf("Either wrong file name or file does not exist\n");
		printf("Exiting simulation!\n");
	}
	
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%s", run);	
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf %lf", &atomsx,&atomsy,&atomsz);	
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf %lf %lf %lf",&spacing,&dx,&dy,&dz,&dt);	
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf",&amp0,&icpower);	
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf",&n0);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%d %d %d",&totalTime,&printFreq,&eneFreq);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%d %lf",&icons,&Mn);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%d %lf %lf",&iwave,&alphaw,&betaw);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%d %lf",&itemp,&ampc);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf %lf",&alpha1,&rho1,&beta1);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf",&sigmaT,&omcut);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf",&w,&u);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%lf %lf %lf %lf %lf",&gamma12,&gamma13,&gamma23,&nc,&sigmac);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%d %d",&restartFlag,&restartTime);
	line = fgets(line, BUFSIZ, in);
	sscanf (line, "%s", restartrun);
	fclose(in);

    if (myid==0) 
    {
		printf("done reading input file\n");
	
		printf("run=%s\n",run);
		printf("numAtomsx=%lf	numAtomsy=%lf	numAtomsz=%lf\n",atomsx,atomsy,atomsz);
		printf("spacing=%lf	dx=%lf dt=%lf\n",spacing,dx,dt);
		printf("amp0=%lf	icpower=%lf\n",amp0,icpower);
		printf("n0=%lf\n",n0);
		printf("totalTime=%d	printFreq=%d	eneFreq=%d\n",totalTime,printFreq,eneFreq);
		printf("icons=%d	Mn=%lf\n",icons,Mn);
		printf("iwave=%d	alphaw=%lf	betaw=%lf\n",iwave,alphaw,betaw);
		printf("itemp=%d	ampc=%lf\n",itemp,ampc);
		printf("alpha1=%lf	rho1=%lf	beta1=%lf\n",alpha1,rho1,beta1);
		printf("sigmaT=%lf	omcut=%lf\n",sigmaT,omcut);
		printf("w=%lf		u=%lf\n",w,u);
		printf("gamma12=%lf		gamma13=%lf    gamma23=%lf    nc=%lf    sigmac=%lf\n",gamma12,gamma13,gamma23,nc,sigmac);
		printf("rFlag=%d	rTime=%d\n",restartFlag,restartTime);
		printf("restartrun=%s\n",restartrun);
    }
	
}





void energy(int time,char *filepre)
{
	ptrdiff_t i,j,k,h;
	ptrdiff_t index,index1,index2;

        double mu1,avgdens1;
        double mu2,avgdens2;
        double mu3,avgdens3;
        double ened1,mud1,avgdensd1;
        double ened2,mud2,avgdensd2;
        double ened3,mud3,avgdensd3;
        int isource,k2;
        MPI_Status status;

	char filename[BUFSIZ];
    	FILE *fpene;

        if (myid==0)
        {
    		sprintf(filename,"%s.ene",filepre);
		fpene = fopen(filename,"a");
        }

	fftw_execute(planF_n1);		//forward transform density
	fftw_execute(planF_n2);		//forward transform density
	fftw_execute(planF_n3);		//forward transform density
	MPI_Barrier(MPI_COMM_WORLD);

	// Correlation part of energy
	for(k=0;k<local_n0;k++)
	{
		index1 = (Nx/2+1)*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + (Nx/2+1)*j;
			for(i=0;i<(Nx/2+1);i++)
			{
				index = i + index2;
                                if (k1arr[index][0] == 0.0)
                                {
                                    kn1[index][0] = kn1[index][0]*rk0;
                                    kn1[index][1] = kn1[index][1]*rk0;
                                }
		                else
                                {
                                    if (iwave==0)
                                    {
                 		        kn1[index][0] = -kn1[index][0]*(1./omegak1[index][0]-1.)/k1arr[index][0];
                 		        kn1[index][1] = -kn1[index][1]*(1./omegak1[index][1]-1.)/k1arr[index][1];
                                    }
                                    else if (iwave==1)
                                    {
			   		//SI1
                 		        kn1[index][0] = -kn1[index][0]*(1./omegak1[index][0]-1.)/k1arr[index][0]/dt;
                 		        kn1[index][1] = -kn1[index][1]*(1./omegak1[index][1]-1.)/k1arr[index][1]/dt;
			   		//SI2
                 		        //kn[index][0] = -kn[index][0]*(1./omegak[index][0]-1.)/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = -kn[index][1]*(1./omegak[index][1]-1.)/karr[index][1]/dt*(1.+betaw*dt);
			   		//EXPLICIT
                 		        //kn[index][0] = kn[index][0]*omegak[index][0]/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = kn[index][1]*omegak[index][1]/karr[index][1]/dt*(1.+betaw*dt);
                                    }
                                }
                                if (k2arr[index][0] == 0.0)
                                {
                                    kn2[index][0] = kn2[index][0]*rk0;
                                    kn2[index][1] = kn2[index][1]*rk0;
                                }
		                else
                                {
                                    if (iwave==0)
                                    {
                 		        kn2[index][0] = -kn2[index][0]*(1./omegak2[index][0]-1.)/k2arr[index][0];
                 		        kn2[index][1] = -kn2[index][1]*(1./omegak2[index][1]-1.)/k2arr[index][1];
                                    }
                                    else if (iwave==1)
                                    {
			   		//SI1
                 		        kn2[index][0] = -kn2[index][0]*(1./omegak2[index][0]-1.)/k2arr[index][0]/dt;
                 		        kn2[index][1] = -kn2[index][1]*(1./omegak2[index][1]-1.)/k2arr[index][1]/dt;
			   		//SI2
                 		        //kn[index][0] = -kn[index][0]*(1./omegak[index][0]-1.)/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = -kn[index][1]*(1./omegak[index][1]-1.)/karr[index][1]/dt*(1.+betaw*dt);
			   		//EXPLICIT
                 		        //kn[index][0] = kn[index][0]*omegak[index][0]/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = kn[index][1]*omegak[index][1]/karr[index][1]/dt*(1.+betaw*dt);
                                    }
                                }
                                if (k3arr[index][0] == 0.0)
                                {
                                    kn3[index][0] = kn3[index][0]*rk0;
                                    kn3[index][1] = kn3[index][1]*rk0;
                                }
		                else
                                {
                                    if (iwave==0)
                                    {
                 		        kn3[index][0] = -kn3[index][0]*(1./omegak3[index][0]-1.)/k3arr[index][0];
                 		        kn3[index][1] = -kn3[index][1]*(1./omegak3[index][1]-1.)/k3arr[index][1];
                                    }
                                    else if (iwave==1)
                                    {
			   		//SI1
                 		        kn3[index][0] = -kn3[index][0]*(1./omegak3[index][0]-1.)/k3arr[index][0]/dt;
                 		        kn3[index][1] = -kn3[index][1]*(1./omegak3[index][1]-1.)/k3arr[index][1]/dt;
			   		//SI2
                 		        //kn[index][0] = -kn[index][0]*(1./omegak[index][0]-1.)/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = -kn[index][1]*(1./omegak[index][1]-1.)/karr[index][1]/dt*(1.+betaw*dt);
			   		//EXPLICIT
                 		        //kn[index][0] = kn[index][0]*omegak[index][0]/karr[index][0]/dt*(1.+betaw*dt);
                 		        //kn[index][1] = kn[index][1]*omegak[index][1]/karr[index][1]/dt*(1.+betaw*dt);
                                    }
                                }
			}
		}
	}
	fftw_execute(planB_NL_n1);	//backward transform energy term
	normalize(n1fNL);		//normalize density after inverse transform
	fftw_execute(planB_NL_n2);	//backward transform energy term
	normalize(n2fNL);		//normalize density after inverse transform
	fftw_execute(planB_NL_n3);	//backward transform energy term
	normalize(n3fNL);		//normalize density after inverse transform
	MPI_Barrier(MPI_COMM_WORLD);

	// Add bulk terms to correlation part of energy
        ene1=0.;
        mu1=0.;
        avgdens1=0.;
        ene2=0.;
        mu2=0.;
        avgdens2=0.;
        ene3=0.;
        mu3=0.;
        avgdens3=0.;
	for(k=0;k<local_n0;k++)
	{
		index1 = Nx2*k*Ny;
           	k2 = k+local_0_start;
		    for(j=0;j<Ny;j++)
		    {
			index2 = index1 + Nx2*j;
			for(i=0;i<Nx;i++)
			{
				index = i + index2;
                               	ene1 = ene1 + .5*n1[index]*( n1fNL[index] + n1[index]*n1[index]*(-w/3. + u*n1[index]/6.) );
                               	mu1 = mu1 + n1fNL[index] + n1[index]*n1[index]*(-w/6. + u*n1[index]/12.);
                               	avgdens1 = avgdens1 + n1[index];
                               	ene2 = ene2 + .5*n2[index]*( n2fNL[index] + n2[index]*n2[index]*(-w/3. + u*n2[index]/6.) );
                               	mu2 = mu2 + n2fNL[index] + n2[index]*n2[index]*(-w/6. + u*n2[index]/12.);
                               	avgdens2 = avgdens2 + n2[index];
                               	ene3 = ene3 + .5*n3[index]*( n3fNL[index] + n3[index]*n3[index]*(-w/3. + u*n3[index]/6.) );
                               	mu3 = mu3 + n3fNL[index] + n3[index]*n3[index]*(-w/6. + u*n3[index]/12.);
                               	avgdens3 = avgdens3 + n3[index];
			}
		    }
	}

	// Tabulate and write avg energy, chemical potential, and density
        if (myid == 0)
	{
		for (isource=1; isource<numprocs; isource++)
		{
                     	MPI_Recv(&ened1,1,MPI_DOUBLE,isource,1,MPI_COMM_WORLD,&status);
                     	ene1 = ene1+ened1;
                     	MPI_Recv(&mud1,1,MPI_DOUBLE,isource,2,MPI_COMM_WORLD,&status);
                     	mu1 = mu1+mud1;
                     	MPI_Recv(&avgdensd1,1,MPI_DOUBLE,isource,3,MPI_COMM_WORLD,&status);
                     	avgdens1 = avgdens1+avgdensd1;
                     	MPI_Recv(&ened2,1,MPI_DOUBLE,isource,1,MPI_COMM_WORLD,&status);
                     	ene2 = ene2+ened2;
                     	MPI_Recv(&mud2,1,MPI_DOUBLE,isource,2,MPI_COMM_WORLD,&status);
                     	mu2 = mu2+mud2;
                     	MPI_Recv(&avgdensd2,1,MPI_DOUBLE,isource,3,MPI_COMM_WORLD,&status);
                     	avgdens2 = avgdens2+avgdensd2;
                     	MPI_Recv(&ened3,1,MPI_DOUBLE,isource,1,MPI_COMM_WORLD,&status);
                     	ene3 = ene3+ened3;
                     	MPI_Recv(&mud3,1,MPI_DOUBLE,isource,2,MPI_COMM_WORLD,&status);
                     	mu3 = mu3+mud3;
                     	MPI_Recv(&avgdensd3,1,MPI_DOUBLE,isource,3,MPI_COMM_WORLD,&status);
                     	avgdens3 = avgdens3+avgdensd3;
		}
                ene1 = ene1/Nx/Ny/Nz;
               	mu1  = mu1/Nx/Ny/Nz;
               	avgdens1 = avgdens1/Nx/Ny/Nz;
               	ene2 = ene2/Nx/Ny/Nz;
               	mu2  = mu2/Nx/Ny/Nz;
               	avgdens2 = avgdens2/Nx/Ny/Nz;
               	ene3 = ene3/Nx/Ny/Nz;
               	mu3  = mu3/Nx/Ny/Nz;
               	avgdens3 = avgdens3/Nx/Ny/Nz;

		fprintf(fpene,"%9f %.17g %.17g %.17g \n",time*dt,ene1,mu1,avgdens1);
	       	fclose(fpene);
	}
        else
	{
                MPI_Send(&ene1,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
                MPI_Send(&mu1,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
                MPI_Send(&avgdens1,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
                MPI_Send(&ene2,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
                MPI_Send(&mu2,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
                MPI_Send(&avgdens2,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
                MPI_Send(&ene3,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
                MPI_Send(&mu3,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
                MPI_Send(&avgdens3,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

}

void add_noise()
{
	ptrdiff_t i,j,k;
	ptrdiff_t index1,index2,index;
	double rnum;
	time_t t;
	srand((unsigned) time(&t)+myid);
	for(k=0;k<local_n0;k++)
	{
		index1 = Nx2*k*Ny;
		for(j=0;j<Ny;j++)
		{
			index2 = index1 + Nx2*j;
			for(i=0;i<Nx;i++)
			{
				index=i + index2;
				rnum = 2.*((double)rand() / RAND_MAX)-1.;
				n1[index] += epdt/dx*rnum;
				rnum = 2.*((double)rand() / RAND_MAX)-1.;
				n2[index] += epdt/dx*rnum;
				rnum = 2.*((double)rand() / RAND_MAX)-1.;
				n3[index] += epdt/dx*rnum;
			}
		}
	}
}

double rand_normal(double mean, double stddev)
{
    double n2 = 0.0;
    int n2_cached = 0;

    time_t t;
	srand((unsigned) time(&t)+myid);

    if (!n2_cached) {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x + y*y;
        } while (r==0.0 || r>1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    } else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}


int main(int argc, char* argv[])
{
	int t;	
	int i,j,k,index,index1,index2;	
  	unsigned long randSeed, nn = 100000000;
  	double sigma = 1.0, rnum, sum, stdsum;
	char filename[BUFSIZ];

	//read in the parameter file
	inputVariables();

	//set up domain parameters
	domainParams();

	//setup the fftw mpi environment
	fftwMPIsetup();

	//allocate real and k-space arrays
	allocateArrays();

	MPI_Barrier(MPI_COMM_WORLD);

	//set up fftw plans
	setfftwPlans();				

	MPI_Barrier(MPI_COMM_WORLD);

  	//specifying to use Mersenne twister MT-19937 as the uniform PRNG
//  	gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
//  	srand(time(NULL)+myid);          /* initialization for rand() */
//  	randSeed = rand();               /* returns a non-negative integer */
//      	printf ("myid=%d, seed=%lu\n", myid, randSeed);
//        fflush(stdout);
//  	gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */
  	//test generator
  	/*sum = 0.;
  	stdsum = 0.;
  	for (i = 0; i < nn; i++) {
      		rnum = gsl_ran_gaussian_ziggurat (gBaseRand, sigma);
      		sum = sum + rnum;
      		if (rnum > sigma || rnum < -sigma) stdsum = stdsum + 1.;
    	}
  	printf ("avg=%.17g\n",sum/nn);
  	printf ("stddev=%.17g\n",(nn-stdsum)/nn);
        fflush(stdout);*/

	//determine if this is a restart or initialization
	if (restartFlag == 0)
	{
		initialize(0);
		if (myid==0) printf("output 0\n");
        //energy(0,run);
		output("1",0,n1,run);
		output("2",0,n2,run);
		output("3",0,n3,run);
		if (iwave==1)
		{
			output("g1",0,g1,run);
			output("g2",0,g2,run);
			output("g3",0,g3,run);
		}
	}
	else
	{
		if (myid==0) 
		{
			printf("restart flag has been triggered\n");
			printf("restart from time iteration %d\n",restartTime);
		}
		restart(restartTime,n1,"1",restartrun);
		restart(restartTime,n2,"2",restartrun);
		restart(restartTime,n3,"3",restartrun);
		if (iwave==1)
		{
			restart(restartTime,g1,"g1",restartrun);
			restart(restartTime,g2,"g2",restartrun);
			restart(restartTime,g3,"g3",restartrun);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (myid==0) printf("beginning iterations\n");

	for(t=restartTime+1;t<totalTime;t++)
	{
		dy = (1+0.05*t/totalTime)*dy;
		//caluclate the SI correlation kernel and Laplacian arrays in kspace
		calcCorrelations(0);

		if(t%1000==1)
			//rand_normal(0.0,ampc);
			add_noise();

		fftw_execute(planF_n1);		//forward transform density
		fftw_execute(planF_n2);		//forward transform density
		fftw_execute(planF_n3);		//forward transform density
		calcNL(t);			//calculate non-linear terms
		MPI_Barrier(MPI_COMM_WORLD);
		if (iwave==0) timeStepn();
		else if (iwave==1) timeStepnWave();
		MPI_Barrier(MPI_COMM_WORLD);
		fftw_execute(planB_n1);		//inverse transform density
		fftw_execute(planB_n2);		//inverse transform density
		fftw_execute(planB_n3);		//inverse transform density
		normalize(n1);	//normalize density after inverse transform
		normalize(n2);	//normalize density after inverse transform
		normalize(n3);	//normalize density after inverse transform
		// Complete last part of wave algorithm
		if (iwave==1) 
        {
			for(k=0;k<local_n0;k++)
			{
				index1 = Nx2*k*Ny;
				for(j=0;j<Ny;j++)
				{
					index2 = index1 + Nx2*j;
					for(i=0;i<Nx;i++)
					{
						index = i + index2;
						g1[index]=(n1[index]-g1[index])/dt;
						g2[index]=(n2[index]-g2[index])/dt;
						g3[index]=(n3[index]-g3[index])/dt;
					}
				}
			}
        }
        
		MPI_Barrier(MPI_COMM_WORLD);

		// Output density and, if iwave=1, g
		if( t%printFreq == 0 )
		{
			if (myid==0) {printf("time = %d\n",t); fflush(stdout);}
			output("1",t,n1,run);
			output("2",t,n2,run);
			output("3",t,n3,run);
			if (iwave==1)
			{
				output("g1",t,g1,run);
				output("g2",t,g2,run);
				output("g3",t,g3,run);
			}
		}

		// Output energy
		//if( t%eneFreq == 0 ) energy(t,run);

	}

	MPI_Barrier(MPI_COMM_WORLD);
//  	gsl_rng_free(gBaseRand);	//free rand number generator
	freeMemory();  //free memory and destroy fftw plans and arrays

	//exit mpi environment
	MPI_Finalize();

	return 0;
}
