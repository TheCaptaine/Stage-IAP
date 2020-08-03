#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cstdio>

using namespace std;

void rk4(int n, double x, double *y, double dx, int f(int, double, double*, double*, double), double q);
int luminosite(int n, double t, double *y, double *dy, double q);
double temperature(double L, double t);
double flux(double t, double T);
double magnitude(double f);

const double h = 6.62607004e-34; // constante de planck, unite : SI
const double Kb = 1.38064852e-23; // constante de boltzman, unite : SI
const double Cstefan = 5.670374e-8; // constante de stefan, unite : SI
const double Conversion_Parsec_en_m = 3.086e16; // unite : m
const double Conversion_Mo_en_kg = 2e30; // unite : kg
const double Conversion_erg_en_J = 1e-7; // unite J
const double c = 3e8; // unite : m/s
const double Eps_star = 1e10; // unite : w
const double t_star = 86.4; // unite : s
const float alpha = 1.3;

double lambda = 0; // unite : m 
double Mej = 0; // unite : kg
double Vexp = 0; // unite : m/s
double kappa = 0; // unite : m^2/kg
double q1 = 0.5;
double q2 = 0.5;

int main() {
	double q = 0;
	char namefile1[256];
	char namefile2[256];
	char namefile3[256];
	char namefile4[256];
	int n = 1;
	double t0 = 1000.; // unite : s
	double l0 = 1e18; // unite : w
	double lambda_valeurs[11] = {2.2e-6, 1.6e-6, 1.2e-6, 1e-6, 0.9e-6, 0.8e-6, 0.7e-6, 0.55e-6, 0.45e-6, 0.365e-6, 0.27e-6};
	double y[1] = {l0};
	int Tmax = 30 * (3600*24); // unite : s
	double dt = 0.01 * (3600*24); // unite : s
	ofstream files5("namefiles.res");
	for (int k = 0; k < 2; k++) {
		if (k == 0) { // equatorial
			Mej = 0.05 * Conversion_Mo_en_kg;
			Vexp = 0.149*c;
			kappa = 0.365;
			q = q1;
		}
		else { // polaire
			Mej = 2.3e-2 * Conversion_Mo_en_kg;
			Vexp = 0.256*c;
			kappa = 0.05;
			q = q2;
		}
		for (int j = 0; j < 11; j++) {
			lambda = lambda_valeurs[j];
			y[0] = l0;
			sprintf(namefile1, "luminosite_lambda%g_m%g_v%g_k%g.res", lambda, Mej, Vexp, kappa);
			sprintf(namefile2, "temperature_lambda%g_m%g_v%g_k%g.res", lambda, Mej, Vexp, kappa);
			sprintf(namefile3, "flux_lambda%g_m%g_v%g_k%g.res", lambda, Mej, Vexp, kappa);
			sprintf(namefile4, "magnitude_lambda%g_m%g_v%g_k%g.res", lambda, Mej, Vexp, kappa);
			ofstream files1(namefile1);
			ofstream files2(namefile2);
			ofstream files3(namefile3);
			ofstream files4(namefile4);
			for (double t = t0; t < Tmax; t += dt) {
				rk4(n, t, y, dt, luminosite, q);
				files1 << t / (3600*24) << " " << y[0] << endl;
				files2 << t / (3600*24) << " " << temperature(y[0], t) << endl;
				files3 << t / (3600*24) << " " << flux(t, temperature(y[0], t)) << endl;
				files4 << t / (3600*24) << " " << magnitude(flux(t, temperature(y[0], t))) << endl;
			}
			files1.close();
			files2.close();
			files3.close();
			files4.close();
			files5 << namefile1 << "\n" << namefile2 << "\n" << namefile3 << "\n" << namefile4 << endl;
		}
	}
	files5.close();
	FILE *f = popen("python3 codeaveccorrection.py", "w");
	fflush(f);
	fclose(f);
	return 0;
}

int luminosite(int n, double t, double *y, double *dy, double q) {
	double L_star = Eps_star*Mej*pow(t/t_star, -alpha);
	dy[0] = (-y[0] + L_star)/(3*q*kappa*Mej/(4*M_PI*c*Vexp*t));
	return 0;
}

double temperature(double L, double t) {
	return pow(L/(4*M_PI*Vexp*Vexp*t*t*Cstefan), 0.25);
}

double flux(double t, double T) {
	double nu = c/lambda; // unite : hz
	double D = 40e6 * Conversion_Parsec_en_m; // unite : m
	return (Vexp*t/D)*(Vexp*t/D)*(2*M_PI*h*nu*nu*nu/(c*c)*1/(exp(h*nu/(Kb*T))-1));
}

double magnitude(double f) {
	return -2.5*log10(f)-56.1; 
}

void rk4(int n, double x, double *y, double dx, int f(int, double, double*, double*, double), double q) {
	int i;
	double ddx;
	double d1[n], d2[n], d3[n], d4[n], yp[n];
	ddx = dx/2;
	f(n, x, y, d1, q);
	for (i = 0; i < n; i++)
		yp[i] = y[i] + d1[i] * ddx;
	f(n, x+ddx, yp, d2, q);
	for (i = 0; i < n; i++)
		yp[i] = y[i] + d2[i] * ddx;
	f(n, x+ddx, yp, d3, q);
	for (i = 0; i < n; i++)
		yp[i] = y[i] + d3[i] * dx;
	f(n, x+dx, yp, d4, q);
	for (i = 0; i < n; i++)
		y[i] = y[i] + dx * (d1[i] + 2*d2[i] + 2*d3[i] + d4[i])/6;
}
