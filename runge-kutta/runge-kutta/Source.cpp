#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>


#define eps 1e-6
#define Cx 0.3
//#define P 1000000
#define K 2.5
#define alf 0
#define phi 0
#define d 0.5
#define ro0 1.2250
#define g 9.8066
#define I 1e3
#define EarthR 6371000
#define fuel0 8e3
#define m0 1e4

int P = 1e5;
double s = M_PI*(d*d / 4.);

using namespace std;

double velocityFunc(double v, double m, double teta, double y);
double tetaFunc(double v, double m, double teta, double y);
double lengthFunc(double v, double teta, double y);
double hightFunc(double v, double teta);
double massFunc(double m);
double Gc(double m);
double X(double ro, double v);
double Y(double X);
double ro(double h);

void main()
{
	double x0 = 0, y0 = 0, v0 = 4.7, teta0 = 80;
	
	double x, y, m, v, teta;
	double t0 = 0, h = 0.0251, t;
	double maxH = -1;
	int count = 0;
	vector<double> k1, k2, k3, k4;		// {v, teta, x, y, m}

	ofstream fx("x.txt"), fy("y.txt"), fteta("teta.txt"), ft("t.txt"), fv("v.txt");
	fx << x0 << endl;
	fy << y0 << endl;
	fteta << teta0 << endl;
	ft << t0 << endl;
	fv << v0 << endl;


	// k1
	k1.emplace_back(h*velocityFunc(v0, m0, teta0, y0));
	k1.emplace_back(h*tetaFunc(v0, m0, teta0, y0));
	k1.emplace_back(h*lengthFunc(v0, teta0, y0));
	k1.emplace_back(h*hightFunc(v0, teta0));
	k1.emplace_back(h*massFunc(m0));
	
	// k2
	k2.emplace_back(h*velocityFunc(v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2.,  y0 + k1.at(4) / 2.));
	k2.emplace_back(h*tetaFunc(v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*lengthFunc(v0 + k1.at(0) / 2., teta0 + k1.at(2) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*hightFunc(v0 + k1.at(0) / 2., teta0 + k1.at(2) / 2.));
	k2.emplace_back(h*massFunc(m0 + k1.at(4) / 2.));

	// k3
	k3.emplace_back(h*velocityFunc(v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*tetaFunc(v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*lengthFunc(v0 + k2.at(0) / 2., teta0 + k2.at(2) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*hightFunc(v0 + k2.at(0) / 2., teta0 + k2.at(2) / 2.));
	k3.emplace_back(h*massFunc(m0 + k2.at(4) / 2.));
	
	// k4
	k4.emplace_back(h*velocityFunc(v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), y0 + k3.at(4)));
	k4.emplace_back(h*tetaFunc(v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), y0 + k3.at(4)));
	k4.emplace_back(h*lengthFunc(v0 + k3.at(0), teta0 + k3.at(2), y0 + k3.at(4)));
	k4.emplace_back(h*hightFunc(v0 + k3.at(0), teta0 + k3.at(2)));
	k4.emplace_back(h*massFunc(m0 + k3.at(4)));

	v = v0 + (k1.at(0) + 2 * k2.at(0) + 2 * k3.at(0) + k4.at(0)) / 6.;
	teta = teta0 + (k1.at(1) + 2 * k2.at(1) + 2 * k3.at(1) + k4.at(1)) / 6.;
	x = x0 + (k1.at(2) + 2 * k2.at(2) + 2 * k3.at(2) + k4.at(2)) / 6.;
	y = y0 + (k1.at(3) + 2 * k2.at(3) + 2 * k3.at(3) + k4.at(3)) / 6.;
	m = m0 + (k1.at(4) + 2 * k2.at(4) + 2 * k3.at(4) + k4.at(4)) / 6.;
	t = t0 + h;

	count++;

	fx << x << endl;
	fy << y << endl;

	while ((y) > eps)
	{
		maxH = y;

		// k1
		k1[0] = h*velocityFunc(v, m, teta, y);
		k1[1] = h*tetaFunc(v, m, teta, y);
		k1[2] = h*lengthFunc(v, teta, y);
		k1[3] = h*hightFunc(v, teta);
		k1[4] = h*massFunc(m);

		// k2
		k2[0] = h*velocityFunc(v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., y + k1.at(4) / 2.);
		k2[1] = h*tetaFunc(v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., y + k1.at(4) / 2.);
		k2[2] = h*lengthFunc(v + k1.at(0) / 2., teta + k1.at(2) / 2., y + k1.at(4) / 2.);
		k2[3] = h*hightFunc(v + k1.at(0) / 2., teta + k1.at(2) / 2.);
		k2[4] = h*massFunc(m + k1.at(4) / 2.);

		// k3
		k3[0] = h*velocityFunc(v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., y + k2.at(4) / 2.);
		k3[1] = h*tetaFunc(v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., y + k2.at(4) / 2.);
		k3[2] = h*lengthFunc(v + k2.at(0) / 2., teta + k2.at(2) / 2., y + k2.at(4) / 2.);
		k3[3] = h*hightFunc(v + k2.at(0) / 2., teta + k2.at(2) / 2.);
		k3[4] = h*massFunc(m + k2.at(4) / 2.);

		// k4
		k4[0] = h*velocityFunc(v + k3.at(0), m + k3.at(1), teta + k3.at(2), y + k3.at(4));
		k4[1] = h*tetaFunc(v + k3.at(0), m + k3.at(1), teta + k3.at(2), y + k3.at(4));
		k4[2] = h*lengthFunc(v + k3.at(0), teta + k3.at(2), y + k3.at(4));
		k4[3] = h*hightFunc(v + k3.at(0), teta + k3.at(2));
		k4[4] = h*massFunc(m + k3.at(4));


		v += (k1.at(0) + 2 * k2.at(0) + 2 * k3.at(0) + k4.at(0)) / 6.;
		teta += (k1.at(1) + 2 * k2.at(1) + 2 * k3.at(1) + k4.at(1)) / 6.;
		x += (k1.at(2) + 2 * k2.at(2) + 2 * k3.at(2) + k4.at(2)) / 6.;
		y += (k1.at(3) + 2 * k2.at(3) + 2 * k3.at(3) + k4.at(3)) / 6.;
		m += (k1.at(4) + 2 * k2.at(4) + 2 * k3.at(4) + k4.at(4)) / 6.;
		t += h;

		count++;

		fx << x << endl;
		fy << y << endl;
		fteta << teta << endl;
		ft << t << endl;
		fv << v << endl;
	}
	std::cout << "Hight = " << y / 1000. << " km" << endl
		<< "Teta = " << teta << endl
		<< "Mass = " << m << " kg" << endl
		<< "Time = " << t << " sec" << endl
		<< "Iteration count = " << count << endl;
	getchar();
}


double ro(double H)
{
	return (ro0*pow(M_E, -H / 7170.));
}
double X(double ro, double v)
{
	return (Cx*ro*v*v*s / 2.);
}
double Y(double X)
{
	return K*X;
}
double Gc(double m)
{
	if ((m - (m0 - fuel0)) > eps)
		return (P / I);
	else
	{
		P = 0;
		return 0;
	}
}

double velocityFunc(double v, double m, double teta, double y)
{
	return (P*cos(M_PI*(alf + phi) / 180.) - X(ro(y/1000.), v) - m*g*sin(M_PI*teta / 180.)) / m;
}
double tetaFunc(double v, double m, double teta, double y)
{
	return (P*sin(M_PI*(alf + phi) / 180.) + Y(X(ro(y/1000.), v)) - m*g*cos(M_PI*teta / 180.) + m*v*v*cos(M_PI*teta / 180.) / (EarthR + y)) / (m*v);
}
double lengthFunc(double v, double teta, double y)
{
	return abs(v)*cos(M_PI*teta / 180.)*(EarthR / (EarthR + y));
}
double hightFunc(double v, double teta)
{
	return v*sin(M_PI*teta / 180.);
}
double massFunc(double m)
{
	return -Gc(m);
}