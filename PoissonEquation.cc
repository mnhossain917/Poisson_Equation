#include <iostream>
#include <math.h>
#include "DoubleInt.cc"
#include "GaussElmnFunction.cc"
// #include "lapacke.h"

using namespace std;

// Function Declaration 
void Kmat_LocalGenerate(double dx, double dy, double **Kmat_Local);
void Kmat_GlobalGenerate(double Lx, double Ly, int nnx, int nny,  double **Kmat_Global);
void Fvec_GlobalGenerate(double Lx, double Ly, int nnx, int nny,  double *Fvec_Global);
void KmatFvec_GlobalBcImp_Form(double Lx, double Ly, int nnx, int nny, double **Kmat_Global, double *Fvec_Global, int *NodeAtBoundry, int *NodeInBoundry, double **Kmat_Global_BcImpd, double *Fvec_Global_BcImpd);

// Forcing Function
double givenForceFunction(double x, double y) 
{ return (-(-(y*y)+y)*sin(3.14159*x));} 

// Shape Function
double ShapeFun_1(double a_x, double b_y, double x, double y) 
{ return ((1-(x/a_x))*(1-(y/b_y)));} 

double ShapeFun_2(double a_x, double b_y, double x, double y) 
{ return ((1-(x/a_x))*(y/b_y));} 

double ShapeFun_3(double a_x, double b_y, double x, double y) 
{ return ((x/a_x)*(y/b_y));} 

double ShapeFun_4(double a_x, double b_y, double x, double y) 
{ return ((x/a_x)*(1-(y/b_y)));} 


// int Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[])
{   
    double Lx=1, Ly=1; // length in X and Y direction
    int nnx=5, nny=5;  // number of DOF in X and Y direction
    int nNode=nnx*nny; // Total number of nodes
    
    // Global Stiffness Matrix 
    double **Kmat_Global;  
    Kmat_Global=new double*[nNode]; // Dynamically create array of pointers 
    for(int i=0;i<nNode;i++)
        {Kmat_Global[i]=new double[nNode];} // Dynamic allocation of memory
    
    Kmat_GlobalGenerate(Lx, Ly, nnx, nny, Kmat_Global);

    // Global Load vector
    double *Fvec_Global;
    Fvec_Global=new double[nNode];
    Fvec_GlobalGenerate(Lx, Ly, nnx, nny, Fvec_Global);
    
    // Array for storing Boundary nodes
    int nNodeBC=2*(nnx+nny-2);
    int *NodeAtBoundry;
    NodeAtBoundry=new int[nNodeBC];
    
    // Array for storing nodes inside boundary
    int nNodeInBC=nNode-(2*(nnx+nny-2));
    int *NodeInBoundry;
    NodeInBoundry=new int[nNodeInBC];
    
    // Memory allocation for global stiffness matrix 
    // after imposing boundary condition
    double **Kmat_Global_BcImpd;  
    Kmat_Global_BcImpd=new double*[nNodeInBC]; 
    for(int i=0;i<nNodeInBC;i++)
        {Kmat_Global_BcImpd[i]=new double[nNodeInBC];} 
    
    // Memory allocation for global Force vector
    // after imposing boundary condition
    double *Fvec_Global_BcImpd;
    Fvec_Global_BcImpd=new double[nNodeInBC];
  
    // Global stiffness matrix And Force vector
    // after imposing boundary condition
    KmatFvec_GlobalBcImp_Form( Lx,  Ly,  nnx,  nny, Kmat_Global, Fvec_Global, NodeAtBoundry, NodeInBoundry, Kmat_Global_BcImpd, Fvec_Global_BcImpd);
    
    // Solving Linear system of equation
    // Solution is stored in "Fvec_Global_BcImpd"
    GaussElmn_Seq(Kmat_Global_BcImpd, Fvec_Global_BcImpd, nNodeInBC, nNodeInBC);
    
    // Solution with boundary condition
    double *SolVecGlobal;
    SolVecGlobal=new double[nNode];
    
    for(int i=0;i<nNode;i++)
       {  SolVecGlobal[i]=0;}

    int i_bc;
    for(int i=0;i<nNodeInBC;i++)
    {  i_bc =NodeInBoundry[i]-1;
        SolVecGlobal[i_bc]=Fvec_Global_BcImpd[i]; }

    // Solution in Matrix form  
    double **SolMatGlobal;  
    SolMatGlobal=new double*[nny]; // Dynamically create array of pointers 
    for(int i=0;i<nny;i++)
        {SolMatGlobal[i]=new double[nnx];} // Dynamic allocation of memory
    
    int NodeNumCount=0;
    for(int j=0;j< nnx;j++)
    {  for(int i=nny-1;i> -1;i--)
        {SolMatGlobal[i][j]=SolVecGlobal[NodeNumCount];
        NodeNumCount=NodeNumCount+1;}
    } 

    // Display Solution
    for(int i=0;i<nny;i++)
    {  for(int j=0;j<nnx;j++)
         {cout<<SolMatGlobal[i][j]<<"\t";}
        cout<<"\n";
    }
    
}

// %%%%% Functions %%%%%%%

// Global Force Matrix Function
void KmatFvec_GlobalBcImp_Form(double Lx, double Ly, int nnx, int nny, double **Kmat_Global, double *Fvec_Global, int *NodeAtBoundry, int *NodeInBoundry, double **Kmat_Global_BcImpd, double *Fvec_Global_BcImpd)
{   int nNode=nnx*nny; // Total number of nodes
    double dx=Lx/(nnx-1), dy=Ly/(nny-1); // length and width of each elements
    int nxElm=(nnx-1), nyElm=(nny-1); 
    int nElm=nxElm*nyElm; //Total number of elements

    int nNodeLocal=4; // Number of Local nodes

    // GLobal Coordinate for 2D
    double **global_Coord;  
    global_Coord=new double*[2]; // Dynamically create array of pointers 
    for(int i=0;i<2;i++)
        {global_Coord[i]=new double[nNode];} // Dynamic allocation of memory
    
    int n_cord=0;
    for(int i=0; i<nnx; i++)
    {  for(int j=0; j<nny; j++)
       {  global_Coord[0][n_cord]=i*dx;
          global_Coord[1][n_cord]=j*dy;
          n_cord=n_cord+1;
       }
    }

    // Boundry nodes
    int nNodeBC=2*(nnx+nny-2);
    int count=0;
    for(int i=0;i<nNode;i++)
    { if(global_Coord[0][i]==0 || global_Coord[1][i]==0 
         || global_Coord[0][i]==Lx  || global_Coord[1][i]==Ly )
        { NodeAtBoundry[count]=i+1 ;
          count=count+1; }
    }
    
    // Node Inside BC
    int nNodeInBC=nNode-(2*(nnx+nny-2));       
    int count1=0;
    for(int i=0;i<nNode;i++)
    { if( global_Coord[0][i]>0 && global_Coord[0][i]<Lx
          && global_Coord[1][i]>0  && global_Coord[1][i]<Ly )
        { NodeInBoundry[count1]=i+1 ;
          count1=count1+1; }
    }
    
    // Global Matrix after imposing Boundry Condition
    int i_bc, j_bc;
    for(int i=0;i<nNodeInBC;i++)
      { for(int j=0;j<nNodeInBC;j++)
        {   i_bc =NodeInBoundry[i]-1;
            j_bc=NodeInBoundry[j]-1;
            Kmat_Global_BcImpd[i][j]=Kmat_Global[i_bc][j_bc];
        }
      }

    for(int i=0;i<nNodeInBC;i++)
       {  i_bc =NodeInBoundry[i]-1;
          {Fvec_Global_BcImpd[i]=Fvec_Global[i_bc];}
       }

}

// Global Force Matrix Function
void Fvec_GlobalGenerate(double Lx, double Ly, int nnx, int nny,  double *Fvec_Global)
{   int nNode=nnx*nny; 
    double dx=Lx/(nnx-1), dy=Ly/(nny-1);
    
    int nxElm=(nnx-1), nyElm=(nny-1);
    int nElm=nxElm*nyElm;

    int nNodeLocal=4;

    // Assigning memory for  element matrix
    double **Fvec_Local_Matform;  
    Fvec_Local_Matform=new double*[nNodeLocal]; // Dynamically create array of pointers 
    for(int i=0;i<nNodeLocal;i++)
        {Fvec_Local_Matform[i]=new double[nElm];} // Dynamic allocation of memory
    

    // GLobal Coordinate for 2D
    double **global_Coord;  
    global_Coord=new double*[2]; // Dynamically create array of pointers 
    for(int i=0;i<2;i++)
        {global_Coord[i]=new double[nNode];} // Dynamic allocation of memory
    
    int n_cord=0;
    for(int i=0; i<nnx; i++)
    {  for(int j=0; j<nny; j++)
       {  global_Coord[0][n_cord]=i*dx;
          global_Coord[1][n_cord]=j*dy;
          n_cord=n_cord+1;
       }
    }
    
    // Connecting Global and local node number 
    double **globalNumbers;
    globalNumbers= new double*[nny];
    for(int i=0;i<nny;i++)
       {globalNumbers[i]=new double[nnx];} 
    
    int NodeNum=1;
    for(int j=0; j<nnx;j++)
    { for(int i=0; i<nny;i++)
        { globalNumbers[i][j]=NodeNum; 
        NodeNum=NodeNum+1;}
    } 


    double **global_Num; 
    global_Num= new double*[nNodeLocal];
    for(int i=0;i<nNodeLocal;i++)
       {global_Num[i]=new double[nElm];} 
    
    int iElm=0;
    for(int i=0; i<nxElm;i++)
    {  for(int j=0; j<nyElm;j++)
       { global_Num[0][iElm]=globalNumbers[j][i];
         global_Num[1][iElm]=globalNumbers[j+1][i];
         global_Num[2][iElm]=globalNumbers[j+1][i+1];
         global_Num[3][iElm]=globalNumbers[j][i+1];
          
         iElm =iElm+1;

       }
    }
    
    int ii, jj;
    for(int k=0;k<nElm;k++)
    {   ii=global_Num[0][k]-1;
        double lx_k=global_Coord[0][ii];
        double ly_k=global_Coord[1][ii];
        double ux_k=lx_k+dx;
        double uy_k=ly_k+dy;
  
        Fvec_Local_Matform[0][k]=doubleIntegral( lx_k,  ux_k,  ly_k,  uy_k, ShapeFun_1,  givenForceFunction);
        Fvec_Local_Matform[1][k]=doubleIntegral( lx_k,  ux_k,  ly_k,  uy_k, ShapeFun_2,  givenForceFunction);
        Fvec_Local_Matform[2][k]=doubleIntegral( lx_k,  ux_k,  ly_k,  uy_k, ShapeFun_3,  givenForceFunction);
        Fvec_Local_Matform[3][k]=doubleIntegral( lx_k,  ux_k,  ly_k,  uy_k, ShapeFun_4,  givenForceFunction);
    }
    
    for(int i=0;i<nNode;i++)
    {Fvec_Global[i]=0;}

    for(int k=0;k<nElm;k++)
    {   for(int i=0;i<nNodeLocal;i++)
        { ii=global_Num[i][k]-1;
          Fvec_Global[ii]=Fvec_Global[ii]+Fvec_Local_Matform[i][k];
        }
    }

}


// Function for Global Stiffness Matrix 
void Kmat_GlobalGenerate(double Lx, double Ly, int nnx, int nny,  double **Kmat_Global)
{   int nNode=nnx*nny;
    double dx=Lx/(nnx-1), dy=Ly/(nny-1);
    int nxElm=(nnx-1), nyElm=(nny-1);
    int nElm=nxElm*nyElm;
    int nNodeLocal=4;

    // Assigning memory for  element matrix
    double **Kmat_Local;  
    Kmat_Local=new double*[nNodeLocal]; // Dynamically create array of pointers 
    for(int i=0;i<nNodeLocal;i++)
        {Kmat_Local[i]=new double[nNodeLocal];} // Dynamic allocation of memory
    
    // Formulation of Local Matrix 
    Kmat_LocalGenerate(dx, dy, Kmat_Local);
    
    // GLobal Coordinate for 2D
    double **global_Coord;  
    global_Coord=new double*[2]; // Dynamically create array of pointers 
    for(int i=0;i<2;i++)
        {global_Coord[i]=new double[nNode];} // Dynamic allocation of memory
    
    int n_cord=0;
    for(int i=0; i<nnx; i++)
    {  for(int j=0; j<nny; j++)
       {  global_Coord[0][n_cord]=i*dx;
          global_Coord[1][n_cord]=j*dy;
          n_cord=n_cord+1;
       }
    }
    
    // Connecting Global and local node number 
    double **globalNumbers;
    globalNumbers= new double*[nny];
    for(int i=0;i<nny;i++)
       {globalNumbers[i]=new double[nnx];} // Dynamic allocation of memory
    
    int NodeNum=1;
    for(int j=0; j<nnx;j++)
    { for(int i=0; i<nny;i++)
        { globalNumbers[i][j]=NodeNum; 
        NodeNum=NodeNum+1;}
    } 


    double **global_Num; 
    global_Num= new double*[nNodeLocal];
    for(int i=0;i<nNodeLocal;i++)
       {global_Num[i]=new double[nElm];} // Dynamic allocation of memory
    
    int iElm=0;
    for(int i=0; i<nxElm;i++)
    {  for(int j=0; j<nyElm;j++)
       { global_Num[0][iElm]=globalNumbers[j][i];
         global_Num[1][iElm]=globalNumbers[j+1][i];
         global_Num[2][iElm]=globalNumbers[j+1][i+1];
         global_Num[3][iElm]=globalNumbers[j][i+1];
          
         iElm =iElm+1;

       }
    }
    
    // Formulation of Global Matrix
    for(int i=0;i<nNode;i++)
    {  for(int j=0;j<nNode;j++)
        Kmat_Global[i][j]=0;
    }
    
    int ii, jj;
    for(int k=0;k<nElm;k++)
    {   for(int i=0;i<nNodeLocal;i++)
          {ii=global_Num[i][k]-1;
           for(int j=0;j<nNodeLocal;j++)
              {jj=global_Num[j][k]-1;
               Kmat_Global[ii][jj]=Kmat_Global[ii][jj]+Kmat_Local[i][j];
               }
          }
    }
}


//  Function for Local Stiffness Matrix 
void Kmat_LocalGenerate(double dx, double dy, double **Kmat_Local)
{   // Diagonal elements
    for (int i=0;i<4;i++)
        Kmat_Local[i][i]=(2*(dx*dx+dy*dy))/(6*dx*dy);

    // Upper traingular part
    Kmat_Local[0][1]= (dy*dy-2*dx*dx)/(6*dx*dy);
    Kmat_Local[0][2]= -(dx*dx+dy*dy)/(6*dx*dy);
    Kmat_Local[0][3]= (dx*dx-2*dy*dy)/(6*dx*dy);

    Kmat_Local[1][2]= (dx*dx-2*dy*dy)/(6*dx*dy);
    Kmat_Local[1][3]= -(dx*dx+dy*dy)/(6*dx*dy);

    Kmat_Local[2][3]= (dy*dy-2*dx*dx)/(6*dx*dy);

    // Lower traingular part
    Kmat_Local[1][0]= (dy*dy-2*dx*dx)/(6*dx*dy);
    Kmat_Local[2][0]= -(dx*dx+dy*dy)/(6*dx*dy);
    Kmat_Local[3][0]= (dx*dx-2*dy*dy)/(6*dx*dy);

    Kmat_Local[2][1]= (dx*dx-2*dy*dy)/(6*dx*dy);
    Kmat_Local[3][1]= -(dx*dx+dy*dy)/(6*dx*dy);

    Kmat_Local[3][2]= (dy*dy-2*dx*dx)/(6*dx*dy);
}