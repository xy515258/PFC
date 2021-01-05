#include<fstream>
#include<iostream>
#include<string>
#include <stdlib.h>
#include <math.h>

using namespace std;
#define PI 3.1415926

int main(int argc, char** argv)
{
    ifstream inf1,inf2,inf3;
    ofstream outf1;
    
    int ym, yp, zm, zp;
    int step, index,indexmy,indexpy,indexmz,indexpz,indexmymz,indexmypz,indexpymz,indexpypz;
    int atom_num = 0;
    int totalstep = 500000;
    double spacing = 0.81675;
    double dx = 0.0680625;
    double dy = 0.0680625;
    double dz = 0.0680625;
    double atomsx = 0.1;
    double atomsy = 80.;
    double atomsz = 80.;
    
    const int NX = (int)(floor(spacing/dx*atomsx+.5));
    const int NY = (int)(floor(spacing/dy*atomsy+.5));
    const int NZ = (int)(floor(spacing/dz*atomsz+.5));
    
    string number(argv[1]);
    double* n1 = new double[NX*NY*NZ];
    double* n2 = new double[NX*NY*NZ];
    double* n3 = new double[NX*NY*NZ];
    double* X = new double[NX*NY*NZ];
    double* Y = new double[NX*NY*NZ];
    double* Z = new double[NX*NY*NZ];
    int* type = new int[NX*NY*NZ];
    
    string filename1("1_a01_9000.dat");
    string filename2("2_a01_9000.dat");
    string filename3("3_a01_9000.dat");

    step = stoi(number);
    dy = (1+0.0*step/totalstep)*dy;
    char filenameout[sizeof "dump.00000000"];
    sprintf(filenameout, "dump.%08d", step);

    double xlo = -1/2.*1.633;
    double xhi = 1/2.*1.633;
    double ylo = -NY*dy/2.;
    double yhi = NY*dy/2.;
    double zlo = -NZ*dz/2.;
    double zhi = NZ*dz/2.;

    filename1.replace(6,4,number);
    filename2.replace(6,4,number);
    filename3.replace(6,4,number);
    
    inf1.open(filename1);
    if(inf1 == NULL)
    {
        perror ("The following error occurred");
        exit(0);
    }
    inf2.open(filename2);
    inf3.open(filename3);
    outf1.open(filenameout);
    
    for(int i=0; i<NX*NY*NZ; i++)
    {
        inf1>>n1[i];
        inf2>>n2[i];
        inf3>>n3[i];
    }
    
    for(int z=0; z<NZ; z++)
    for(int y=0; y<NY; y++)
    {
        index = y+NY*z;
        ym = y-1;
        yp = y+1;
        zm = z-1;
        zp = z+1;
        if(y==0)
            ym = NY-1;
        else if(y==NY-1)
            yp = 0;
        if(z==0)
            zm = NZ-1;
        else if(z==NZ-1)
            zp = 0;
        indexmy = ym+NY*z;
        indexpy = yp+NY*z;
        indexmz = y+NY*zm;
        indexpz = y+NY*zp;
        indexmymz = ym+NY*zm;
        indexmypz = ym+NY*zp;
        indexpymz = yp+NY*zm;
        indexpypz = yp+NY*zp;
        if((n1[index]>=n1[indexmy])&&(n1[index]>n1[indexpy])&&(n1[index]>=n1[indexmz])&&(n1[index]>n1[indexpz])&&(n1[index]>=n1[indexmymz])&&(n1[index]>n1[indexmypz])&&(n1[index]>=n1[indexpymz])&&(n1[index]>n1[indexpypz]))
        {
            atom_num += 1;
            type[atom_num] = 1;
            X[atom_num] = 0.16667;
            Y[atom_num] = (double)y/NY;
            Z[atom_num] = (double)z/NZ;
        }
        if((n2[index]>=n2[indexmy])&&(n2[index]>n2[indexpy])&&(n2[index]>=n2[indexmz])&&(n2[index]>n2[indexpz])&&(n2[index]>=n2[indexmymz])&&(n2[index]>n2[indexmypz])&&(n2[index]>=n2[indexpymz])&&(n2[index]>n2[indexpypz]))
        {
            atom_num += 1;
            type[atom_num] = 2;
            X[atom_num] = 0.5;
            Y[atom_num] = (double)y/NY;
            Z[atom_num] = (double)z/NZ;
        }
        
        if((n3[index]>=n3[indexmy])&&(n3[index]>n3[indexpy])&&(n3[index]>=n3[indexmz])&&(n3[index]>n3[indexpz])&&(n3[index]>=n3[indexmymz])&&(n3[index]>n3[indexmypz])&&(n3[index]>=n3[indexpymz])&&(n3[index]>n3[indexpypz]))
        {
            atom_num += 1;
            type[atom_num] = 3;
            X[atom_num] = 0.83333;
            Y[atom_num] = (double)y/NY;
            Z[atom_num] = (double)z/NZ;
        }
            
    }
    
    cout<<"n= "<<atom_num<<endl;
    
    outf1<<"ITEM: TIMESTEP"<<endl;
    outf1<<step<<endl;
    outf1<<"ITEM: NUMBER OF ATOMS"<<endl;
    outf1<<atom_num<<endl;
    outf1<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
    outf1<<xlo<<"  "<<xhi<<endl;
    outf1<<ylo<<"  "<<yhi<<endl;
    outf1<<zlo<<"  "<<zhi<<endl;
    outf1<<"ITEM: ATOMS id type xs ys zs"<<endl;
    
    for(int i=1;i<=atom_num;i++)
    {
        outf1<<i<<" "<<type[i]<<" "<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<endl;
    }
    
    
    inf1.close();
    inf2.close();
    inf3.close();
    outf1.close();
    return 0;
}

