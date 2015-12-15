// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* y0);
void RKstep(double* const yn, const double* const y0,double* k1, double* k2, double* k3, double* k4, const double x, const double dx);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
        const int dim = 2;
	double dx = 0.1,x=0;
	const double L = 35;
	
	for(double p0 = 0.1; p0<=5.; p0+=0.2){
	
        double y0[dim] = {p0, 0.};
	double yn[dim];
	double theta, thetal, thetar,a,b1,b2,b4,k1[2],k2[2],k3[2],k4[2];
	x=0;

        out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0,k1,k2,k3,k4,x,dx);
		if(y0[1]>0 && yn[1]<0) break;
		
		
		
                for(int i=0; i<dim; i++) y0[i] = yn[i];
		out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << y0[2] << endl;
	}
              a=1.;
	      
	      thetal=0.;
	      thetar=1.;
	
	      while (abs(a)>1e-6){
	        theta=(thetal+thetar)/2.0;
	      
		b1=theta - 3.*theta*theta/2. + 2.*theta*theta*theta/3.;
	        b2=theta*theta - 2.*theta*theta*theta/3. ;
	        b4=-theta*theta/2. + 2.*theta*theta*theta/3. ;
	      
	        a = y0[1] + dx*(b1*k1[1]+b2*k2[1]+b2*k3[1]+b4*k4[1]);
		
		if (a>0.) thetal=theta;
		else thetar=theta;
		
	
	       }
	       
	       cout << p0 << "\t" << x + theta*dx -dx << endl;
	}
	
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0, double* k1, double* k2, double* k3, double* k4,
            const double x, const double dx)
{
	const int dim =2;
	

        for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
        f(k2);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3);

        for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
// Lorenz model
void f(double* y0)
{
	double y[2] = { y0[0], y0[1] };

        y0[0] = y[1];
	y0[1] = -y[0] /sqrt(1+y[0]*y[0]);

}
