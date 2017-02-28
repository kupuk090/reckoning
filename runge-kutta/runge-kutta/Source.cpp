#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>


#define Cx 0.3
//#define m0 100000
#define P 923400
#define K 2.5
#define alf 0
#define phi 0
#define d 0.5
#define ro0 1.2250
#define R 8.3145
#define m_aero 28.97
#define g 9.81
#define I 10000
#define EarthR 6371000

using namespace std;

double f1(double t, double v, double m, double teta, double x, double y);
double f2(double t, double v, double m, double teta, double x, double y);
double f3(double t, double v, double m, double teta, double x, double y);
double f4(double t, double v, double m, double teta, double x, double y);
double f5(double t, double v, double m, double teta, double x, double y);
double X(double ro, double v);
double Y(double X);
double temp(double h);
double ro(double h);
double s();

void main()
{
	double x0 = 0, y0 = 0, m0 = 100000, v0 = 20, teta0 = 90;
	double x, y, m, v, teta;
	double t0 = 0, h = 2, t;
	double maxH = -1;
	int count = 0;
	vector<double> k1, k2, k3, k4;

	k1.emplace_back(h*f1(t0, v0, m0, teta0, x0, y0));
	k1.emplace_back(h*f2(t0, v0, m0, teta0, x0, y0));
	k1.emplace_back(h*f3(t0, v0, m0, teta0, x0, y0));
	k1.emplace_back(h*f4(t0, v0, m0, teta0, x0, y0));
	k1.emplace_back(h*f5(t0, v0, m0, teta0, x0, y0));

	k2.emplace_back(h*f1(t0 + h / 2., v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., x0 + k1.at(3) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*f2(t0 + h / 2., v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., x0 + k1.at(3) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*f3(t0 + h / 2., v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., x0 + k1.at(3) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*f4(t0 + h / 2., v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., x0 + k1.at(3) / 2., y0 + k1.at(4) / 2.));
	k2.emplace_back(h*f5(t0 + h / 2., v0 + k1.at(0) / 2., m0 + k1.at(1) / 2., teta0 + k1.at(2) / 2., x0 + k1.at(3) / 2., y0 + k1.at(4) / 2.));

	k3.emplace_back(h*f1(t0 + h / 2., v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., x0 + k2.at(3) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*f2(t0 + h / 2., v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., x0 + k2.at(3) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*f3(t0 + h / 2., v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., x0 + k2.at(3) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*f4(t0 + h / 2., v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., x0 + k2.at(3) / 2., y0 + k2.at(4) / 2.));
	k3.emplace_back(h*f5(t0 + h / 2., v0 + k2.at(0) / 2., m0 + k2.at(1) / 2., teta0 + k2.at(2) / 2., x0 + k2.at(3) / 2., y0 + k2.at(4) / 2.));

	k4.emplace_back(h*f1(t0 + h, v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), x0 + k3.at(3), y0 + k3.at(4)));
	k4.emplace_back(h*f2(t0 + h, v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), x0 + k3.at(3), y0 + k3.at(4)));
	k4.emplace_back(h*f3(t0 + h, v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), x0 + k3.at(3), y0 + k3.at(4)));
	k4.emplace_back(h*f4(t0 + h, v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), x0 + k3.at(3), y0 + k3.at(4)));
	k4.emplace_back(h*f5(t0 + h, v0 + k3.at(0), m0 + k3.at(1), teta0 + k3.at(2), x0 + k3.at(3), y0 + k3.at(4)));

	v = v0 + (k1.at(0) + 2 * k2.at(0) + 2 * k3.at(0) + k4.at(0)) / 6.;
	teta = teta0 + (k1.at(1) + 2 * k2.at(1) + 2 * k3.at(1) + k4.at(1)) / 6.;
	x = x0 + (k1.at(2) + 2 * k2.at(2) + 2 * k3.at(2) + k4.at(2)) / 6.;
	y = y0 + (k1.at(3) + 2 * k2.at(3) + 2 * k3.at(3) + k4.at(3)) / 6.;
	m = m0 + (k1.at(4) + 2 * k2.at(4) + 2 * k3.at(4) + k4.at(4)) / 6.;
	t = t0 + h;

	count++;

	while (y > maxH)
	{
		maxH = y;

		k1[0] = h*f1(t, v, m, teta, x, y);
		k1[1] = h*f2(t, v, m, teta, x, y);
		k1[2] = h*f3(t, v, m, teta, x, y);
		k1[3] = h*f4(t, v, m, teta, x, y);
		k1[4] = h*f5(t, v, m, teta, x, y);

		k2[0] = h*f1(t + h / 2., v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., x + k1.at(3) / 2., y + k1.at(4) / 2.);
		k2[1] = h*f2(t + h / 2., v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., x + k1.at(3) / 2., y + k1.at(4) / 2.);
		k2[2] = h*f3(t + h / 2., v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., x + k1.at(3) / 2., y + k1.at(4) / 2.);
		k2[3] = h*f4(t + h / 2., v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., x + k1.at(3) / 2., y + k1.at(4) / 2.);
		k2[4] = h*f5(t + h / 2., v + k1.at(0) / 2., m + k1.at(1) / 2., teta + k1.at(2) / 2., x + k1.at(3) / 2., y + k1.at(4) / 2.);

		k3[0] = h*f1(t + h / 2., v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., x + k2.at(3) / 2., y + k2.at(4) / 2.);
		k3[1] = h*f2(t + h / 2., v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., x + k2.at(3) / 2., y + k2.at(4) / 2.);
		k3[2] = h*f3(t + h / 2., v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., x + k2.at(3) / 2., y + k2.at(4) / 2.);
		k3[3] = h*f4(t + h / 2., v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., x + k2.at(3) / 2., y + k2.at(4) / 2.);
		k3[4] = h*f5(t + h / 2., v + k2.at(0) / 2., m + k2.at(1) / 2., teta + k2.at(2) / 2., x + k2.at(3) / 2., y + k2.at(4) / 2.);

		k4[0] = h*f1(t + h, v + k3.at(0), m + k3.at(1), teta + k3.at(2), x + k3.at(3), y + k3.at(4));
		k4[1] = h*f2(t + h, v + k3.at(0), m + k3.at(1), teta + k3.at(2), x + k3.at(3), y + k3.at(4));
		k4[2] = h*f3(t + h, v + k3.at(0), m + k3.at(1), teta + k3.at(2), x + k3.at(3), y + k3.at(4));
		k4[3] = h*f4(t + h, v + k3.at(0), m + k3.at(1), teta + k3.at(2), x + k3.at(3), y + k3.at(4));
		k4[4] = h*f5(t + h, v + k3.at(0), m + k3.at(1), teta + k3.at(2), x + k3.at(3), y + k3.at(4));


		v += (k1.at(0) + 2 * k2.at(0) + 2 * k3.at(0) + k4.at(0)) / 6.;
		teta += (k1.at(1) + 2 * k2.at(1) + 2 * k3.at(1) + k4.at(1)) / 6.;
		x += (k1.at(2) + 2 * k2.at(2) + 2 * k3.at(2) + k4.at(2)) / 6.;
		y += (k1.at(3) + 2 * k2.at(3) + 2 * k3.at(3) + k4.at(3)) / 6.;
		m += (k1.at(4) + 2 * k2.at(4) + 2 * k3.at(4) + k4.at(4)) / 6.;
		t += h;

		count++;
	}
	cout << y / 1000. << " km" << endl << count;
	getchar();
}

double s()
{
	return M_PI*(d*d / 4.);
}
double temp(double H)
{
	return (0.0222*H*H - 2.6638*H + 276.4005);
}
double ro(double H)
{
	return (ro0*pow(M_E, -m_aero*g*H / (R*temp(H))));
}
double X(double ro, double v)
{
	return (Cx*ro*v*v*s() / 2.);
}
double Y(double X)
{
	return K*X;
}
double Gc()
{
	return (P / I);
}

double f1(double t, double v, double m, double teta, double x, double y)
{
	return (P*cos(M_PI*(alf + phi) / 180.) - X(ro(y/1000.), v) - m*g*sin(M_PI*teta / 180.)) / m;
}
double f2(double t, double v, double m, double teta, double x, double y)
{
	return (P*sin(M_PI*(alf + phi) / 180.) + Y(X(ro(y/1000.), v)) - m*g*cos(M_PI*teta / 180.) + m*v*v*cos(M_PI*teta / 180.) / (EarthR + y)) / (m*v);
}
double f3(double t, double v, double m, double teta, double x, double y)
{
	return v*cos(M_PI*teta / 180.)*(EarthR / (EarthR + y));
}
double f4(double t, double v, double m, double teta, double x, double y)
{
	return v*sin(M_PI*teta / 180.);
}
double f5(double t, double v, double m, double teta, double x, double y)
{
	return -Gc();
}