#include<fstream>
#include<iostream>
#include<string>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;
#define PI 3.1415926

// Global variables
int grain_size;
int Ny;
int Nz;

//Check if all the atoms in this two grid point are in the same orientation
int check(vector<int>** location_list, double* orientation, int y0, int z0, int y1, int z1)
{
    int same_orientation = 0;
    int interface = 0;
    double avg_ori_0, avg_ori_1;
    avg_ori_0 = 0;
    avg_ori_1 = 0;
    if(location_list[y0][z0].size()&&location_list[y1][z1].size())
    {
        for(int p=0; p < location_list[y0][z0].size(); p++)
        {
            if(orientation[location_list[y0][z0][p]]==180)
                interface = 1;
            avg_ori_0 += orientation[location_list[y0][z0][p]]/location_list[y0][z0].size();
        }
        for(int q=0; q < location_list[y1][z1].size(); q++)
        {
            if(orientation[location_list[y1][z1][q]]==180)
                interface = 1;
            avg_ori_1 += orientation[location_list[y1][z1][q]]/location_list[y1][z1].size();
        }
        if(interface==0&&abs(avg_ori_0-avg_ori_1)<3)
            same_orientation = 1;
    }
    return same_orientation;
}

//DFS algorithm to visit all the grid points and label the points with same orientation
void DFS(int** grid, vector<int>** location_list, double* orientation, int y, int z)
{
    grid[y][z] = 1;
    grain_size += 1;
    if(y>0)
        if(check(location_list,orientation,y,z,y-1,z)&&grid[y-1][z]==0)
            DFS(grid,location_list,orientation,y-1,z);
    if(y==0)
        if(check(location_list,orientation,y,z,Ny-1,z)&&grid[Ny-1][z]==0)
            DFS(grid,location_list,orientation,Ny-1,z);
    if(y<Ny-1)
        if(check(location_list,orientation,y,z,y+1,z)&&grid[y+1][z]==0)
            DFS(grid,location_list,orientation,y+1,z);
    if(y==Ny-1)
        if(check(location_list,orientation,y,z,0,z)&&grid[0][z]==0)
            DFS(grid,location_list,orientation,0,z);
    if(z>0)
        if(check(location_list,orientation,y,z,y,z-1)&&grid[y][z-1]==0)
            DFS(grid,location_list,orientation,y,z-1);
    if(z==0)
        if(check(location_list,orientation,y,z,y,Nz-1)&&grid[y][Nz-1]==0)
            DFS(grid,location_list,orientation,y,Nz-1);
    if(z<Nz-1)
        if(check(location_list,orientation,y,z,y,z+1)&&grid[y][z+1]==0)
            DFS(grid,location_list,orientation,y,z+1);
    if(z==Nz-1)
        if(check(location_list,orientation,y,z,y,0)&&grid[y][0]==0)
            DFS(grid,location_list,orientation,y,0);
}

int main(int argc, char** argv)
{
    ifstream inf1;
    ofstream outf1,outf2;

    string number(argv[1]);
    int step = stoi(number);
    char filenamein[sizeof "dump.00000000"];
    sprintf(filenamein, "dump.%08d", step);
    char filenameout1[sizeof "angle.00000000"];
    sprintf(filenameout1, "angle.%08d", step);
    char filenameout2[sizeof "grain_number.txt"];
    sprintf(filenameout2, "grain_number.txt");
    
    string block;

    int timestep,totalatoms;
    double xlo,xhi,ylo,yhi,zlo,zhi;
    
    inf1.open(filenamein);
    if(inf1 == NULL)
    {
        perror ("The following error occurred");
        exit(0);
    }
    outf1.open(filenameout1);
    outf2.open(filenameout2,ios::app);
    
    inf1 >> block >> block;
    inf1 >> timestep;
    inf1 >> block >> block >> block >> block;
    inf1 >> totalatoms;
    inf1 >> block >> block >> block >> block >> block >> block;
    inf1 >> xlo >> xhi;
    inf1 >> ylo >> yhi;
    inf1 >> zlo >> zhi;
    inf1 >> block >> block >> block >> block >> block >> block >> block;
    
    int* id = new int[totalatoms];
    double* type = new double[totalatoms];
    double* x = new double[totalatoms];
    double* y = new double[totalatoms];
    double* z = new double[totalatoms];
    int* location_y = new int[totalatoms];
    int* location_z = new int[totalatoms];
    double* orientation = new double[totalatoms];
    
    //Make rough grid points
    double cut_off_max = 0.5;
    double cut_off_min = 0.3;
    Ny = (int)((yhi-ylo)/cut_off_max);  //The size of grid depends on cut-off
    Nz = (int)((zhi-zlo)/cut_off_max);
    //Used for the count of grain number, initially 0 but 1 once visited
    int** grid = new int*[Ny];
    //List of atom id in every grid point
    vector<int>** location_list = new vector<int>*[Ny];
    for(int i=0; i<Ny; i++)
    {
        grid[i] = new int[Nz];
        for(int j=0; j<Nz; j++)
            grid[i][j] = 0;  //Give initial value
        location_list[i] = new vector<int>[Nz];
    }
    
    for(int i=0; i<totalatoms; i++)
    {
        inf1>> id[i] >> type[i] >> x[i] >> y[i] >> z[i];
        location_y[i] = (int)(y[i]*Ny);
        location_z[i] = (int)(z[i]*Nz);
        //Put the atoms into the list of their grid point
        location_list[location_y[i]][location_z[i]].push_back(i);
    }
    
    outf1 << "ITEM: TIMESTEP" << endl;
    outf1 << timestep << endl;
    outf1 << "ITEM: NUMBER OF ATOMS" << endl;
    outf1 << totalatoms << endl;
    outf1 << "ITEM: BOX BOUNDS pp pp pp" << endl;
    outf1 << xlo << " " << xhi << endl;
    outf1 << ylo << " " << yhi << endl;
    outf1 << zlo << " " << zhi << endl;
    outf1 << "ITEM: ATOMS id type xs ys zs" << endl;
    
    double dist;
    
    for(int i=0;i<totalatoms;i++)
    {
        orientation[i] = 180;
    }
    
    // Now check every atom
    for(int i=0;i<totalatoms;i++)
    {
        // Check the mid layer atoms
        if(type[i]==2)
        {
            vector<int> neighbor_list;
            int close_atom_num = 0;
            // Check its surrounding grid points
            for(int ii=-1;ii<2;ii++)
            for(int jj=-1;jj<2;jj++)
            {
                int y0 = location_y[i] + ii;
                int z0 = location_z[i] + jj;
                if(y0<0)
                    y0 += Ny;
                else if(y0>=Ny)
                    y0 += -Ny;
                if(z0<0)
                    z0 += Nz;
                else if(z0>=Nz)
                    z0 += -Nz;
                // Check the other type atoms in the surrounding grid points
                for(int p=0; p < location_list[y0][z0].size(); p++)
                {
                    int j = location_list[y0][z0][p];
                    if(type[j]!=2)
                    {
                        double y_dist, z_dist;
                        y_dist = y[j]-y[i];
                        z_dist = z[j]-z[i];
                        if(y_dist>0.5)
                            y_dist += -1;
                        else if(y_dist<-0.5)
                            y_dist += 1;
                        if(z_dist>0.5)
                            z_dist += -1;
                        else if(z_dist<-0.5)
                            z_dist += 1;
                        dist = sqrt(y_dist*y_dist*(yhi-ylo)*(yhi-ylo)+z_dist*z_dist*(zhi-zlo)*(zhi-zlo));
                        // If the distance is within the cut-off range, then list as neighbor atom
                        if((dist-cut_off_min)*(dist-cut_off_max)<0)
                        {
                            neighbor_list.push_back(j);
                            close_atom_num += 1;
                        }
                    }
                }
            }
            // For the atoms with 6 neighbors, calculate it's orientation
            if(close_atom_num == 6)
            {
                orientation[i] = 0;
                // Calculate the orientation by the average angle of all type 1 atoms
                for(int k=0; k<neighbor_list.size(); k++)
                {
                    int p = neighbor_list[k];
                    if(type[p]==1)
                    {
                        double y_dist, z_dist;
                        y_dist = y[p]-y[i];
                        z_dist = z[p]-z[i];
                        if(y_dist>0.5)
                            y_dist += -1;
                        else if(y_dist<-0.5)
                            y_dist += 1;
                        if(z_dist>0.5)
                            z_dist += -1;
                        else if(z_dist<-0.5)
                            z_dist += 1;
                        orientation[i] += ((atan2((z_dist),(y_dist))+PI)*180/PI)/3.;
                    }
                }
                
                // The orientation should range from 0 to 120 degree
                int ori = (int)orientation[i];
                ori = ori/120;
                orientation[i] = orientation[i]-ori*120;
                
                // Set this orientation to all the neighbor atoms
                for(int k=0; k<6; k++)
                {
                    int p = neighbor_list[k];
                    orientation[p] = orientation[i];
                }
            }
        }
    }
    
    // Output the information
    for(int i=0; i<totalatoms; i++)
        outf1<<id[i]<<" "<<orientation[i]<<" "<<x[i]<<" "<<y[i]<<" "<<z[i]<<endl;
    
    //Calculate the average number of atoms on one grid point
    double avg_num = 0;
    for(int i=0; i<Ny; i++)
    for(int j=0; j<Nz; j++)
        avg_num += (double)location_list[i][j].size()/Ny/Nz;
    
    // Now count for the number of grains by DFS
    int count = 0;
    for(int i=0; i<Ny; i++)
    for(int j=0; j<Nz; j++)
        if(grid[i][j]==0)
        {
            grain_size = 0;
            DFS(grid,location_list,orientation,i,j);
            if(grain_size>=20)
                count += 1;
        }
    
    cout << "average number of atoms on grid point: " << avg_num << endl;
    cout << "number of grains: " << count << endl;
    
    outf2 << step << " " << count << endl;
    
    inf1.close();
    outf1.close();
    outf2.close();
    return 0;
}

