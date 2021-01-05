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
    
    int step;
    
    double spacing = 0.81675;
    double dx = 0.0680625;
    double atomsx = 0.1;
    double atomsy = 80.;
    double atomsz = 80.;
    double dy = dx;
    double dz = dx;
    
    const int NX = (int)(floor(spacing/dx*atomsx+.5));
    const int NY = (int)(floor(spacing/dy*atomsy+.5));
    const int NZ = (int)(floor(spacing/dz*atomsz+.5));
    
    string number(argv[1]);
    double* n1 = new double[NX*NY*NZ];
    double* n2 = new double[NX*NY*NZ];
    double* n3 = new double[NX*NY*NZ];
    double* n = new double[NX*NY*NZ];
    string filename1("1_a01_9000.dat");
    string filename2("2_a01_9000.dat");
    string filename3("3_a01_9000.dat");

    step = stoi(number);
    char filenameout[sizeof "data_00000000.vtk"];
    sprintf(filenameout, "data_%08d.vtk", step);

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
    
    outf1 << "# vtk DataFile Version 2.0" << endl;
    outf1 << "n" << endl << "ASCII" << endl << "DATASET STRUCTURED_POINTS" << endl;
    outf1 << "DIMENSIONS " << NX << " " << NY << " " << NZ << endl;
    outf1 << "ASPECT_RATIO 1 1 1" << endl << "ORIGIN 0 0 0" << endl << "POINT_DATA " << NX*NY*NZ << endl;
    outf1 << "SCALARS n1 double" << endl << "LOOKUP_TABLE default" <<endl;
    for(int i=0; i<NX*NY*NZ; i++)
    {
        inf1>>n1[i];
        inf2>>n2[i];
        inf3>>n3[i];
        outf1<<n1[i] << " ";
    }
    outf1 << endl;
    
    outf1 << "SCALARS n2 double" << endl << "LOOKUP_TABLE default" <<endl;
    for(int i=0; i<NX*NY*NZ; i++)
    {
        outf1<<n2[i] << " ";
    }
    outf1 << endl;
    
    outf1 << "SCALARS n3 double" << endl << "LOOKUP_TABLE default" <<endl;
    for(int i=0; i<NX*NY*NZ; i++)
    {
        outf1<<n3[i] << " ";
    }
    outf1 << endl;
    
    outf1 << "SCALARS n double" << endl << "LOOKUP_TABLE default" <<endl;
    for(int i=0; i<NX*NY*NZ; i++)
    {
        n[i] = 0;
        if(n1[i]>0.5)
            n[i] = 1;
        if(n2[i]>0.5)
            n[i] = 2;
        if(n3[i]>0.5)
            n[i] = 3;
        outf1<<n[i] << " ";
    }
        
    inf1.close();
    inf2.close();
    inf3.close();
    outf1.close();
    
    
    return 0;
}

