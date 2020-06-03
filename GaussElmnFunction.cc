// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function : Gauss Elemination Sequencial
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GaussElmn_Seq(double **A, double *b, int n_row, int n_col)
{ // Solving x=inv(A)*b using LU/ Gauss Elimination
  // Constructing Upper trangular matrix
  {for(int k=0; k<n_col-1; k++)
     {for(int i=k+1; i<n_row; i++)
        {double l_ik=A[i][k]/A[k][k];
         for(int j=k;j<n_col;j++)
         {A[i][j]=A[i][j]-l_ik*A[k][j];
          b[i]=b[i]-l_ik*b[k];}
        }
     }

  // Solving Ux=y 
  b[n_row-1]=b[n_row-1]/A[n_row-1][n_col-1];
  for(int k=n_row-2;k>=0;k--)
    {for(int j=k+1;j<n_col;j++)
      {b[k]-=A[k][j]*b[j];
        b[k]=b[k]/A[k][k];}
    }
  }}



