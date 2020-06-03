// Function to find the double integral value 
double doubleIntegral( double lx, double ux, 
					double ly, double uy, 
					double (*ShapeFun)(double,double,double,double), 
					double (*ForceFunction)(double, double) ) 
{  
    // Numner of points 
	int nx=100, ny=100; 

    double dx=(ux - lx);
    double dy=(uy - ly);
    double h=(ux - lx)/(nx-1);
    double k=(uy - ly)/(ny-1);
 
    double answer; 
    double **z;  
    z=new double*[nx]; // Dynamically create array of pointers 
    for(int i=0;i<nx;i++)
        {z[i]=new double[ny];} // Dynamic allocation of memory
    
	double *ax ;
	ax=new double[nx];

	// Values of function in Matrix format
	for (int i = 0; i < nx; i++) 
    { 	for (int j = 0; j < ny; j++) 
        { 	z[i][j] = ShapeFun(dx, dy, i * h, j * k) * ForceFunction(lx + i * h, ly + j * k); 
		} 
	} 

	// 
	for (int i = 0; i < nx; ++i) { 
		ax[i] = 0; 
		for (int j = 0; j < ny; ++j) { 
			if (j == 0 || j == ny - 1) 
				ax[i] += z[i][j]; 
			else if (j % 2 == 0) 
				ax[i] += 2 * z[i][j]; 
			else
				ax[i] += 4 * z[i][j]; 
		} 
		ax[i] *= (k / 3); 
	} 

	answer = 0; 

	// Final integral value 
	for (int i = 0; i < nx; ++i) { 
		if (i == 0 || i == nx - 1) 
			answer += ax[i]; 
		else if (i % 2 == 0) 
			answer += 2 * ax[i]; 
		else
			answer += 4 * ax[i]; 
	} 
	answer *= (h / 3); 

	return answer; 
} 