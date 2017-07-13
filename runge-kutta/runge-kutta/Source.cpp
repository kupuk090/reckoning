#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>



#define eps 1
#define eps_ 1e-3
#define Cx 0.3
#define K 3.5

//#define alf 0

#define phi 0
#define d 0.5
#define ro0 1.2250
#define g 9.8066
#define I1 2.2e3
#define I2 9e3
#define EarthR 6371000

//#define fuel0 8.2e3
//#define fuel01 5e3
#define fuel01 0
//#define fuel02 4.6e3
#define fuel02 9.6e3
#define m0 1.27e4
#define massVeh1 1e3
#define massVeh2 8e2

#define horizontalH 5.4e3
#define horizontalErr 400
#define horizontalFlying false

#define finalPoint 1.e6
#define finalErr 1e2

double alf = 0;
//double m0 = 1.27e4;

double P = 3e5;
double s = 2;
bool firstVehicle = true;
bool secondVehicle = false;

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

map<string, double> RungeKutta(const map<string, double>* data, const double step);
void ballisticTrajectory(map<string, double> itData, double shutoffH, double h, ofstream &fxls, map<string, double>* lastItter);
void constParamTrajectory(map<string, double> itData, double shutoffH, double h, ofstream & fxls, double* shutoffTime, map<string, double>* lastItter);

void main()
{
	double x0 = 0, y0 = 0, v0 = 4.7, teta0 = 80;

	double t0 = 0, h = 0.0025;
	double maxH = -1;
	int count = 0;
	int writeCount = 0;

	double shutoffH = -1;

	bool stopWriting = false;

	map<string, double> itData;
	map<string, double> lastItter;
	vector<map<string, double>> result;

	ofstream fxls("output.xls", ios::out);

	teta0 = teta0*M_PI/180.;

	fxls << "x\t" << "y\t" << "teta\t" << "t\t" << "v\t" << "X\t" << "Y\t" << "m\t"<< "ro\t" << endl
		<< x0/1000. << "\t"
		<< y0/1000. << "\t"
		<< teta0*180/M_PI << "\t"
		<< t0 << "\t"
		<< v0 << "\t"
		<< X(ro(y0), v0) << "\t"
		<< Y(X(ro(y0), v0)) << "\t" 
		<< m0 << "\t" 
		<< ro(y0) << "\t" << endl;
	

//==========================================================================================================================	
	//double left = 0, right = 2e5;
	//bool maxHightAchieved = false;
	//bool trigger1 = false, trigger2 = false;
	//double deltaY;
	//while (true)
	//{
	//	itData["velocity"] = v0;
	//	itData["teta"] = teta0;
	//	itData["length"] = x0;
	//	itData["hight"] = y0;
	//	itData["mass"] = m0;
	//	itData["time"] = t0;

	//	maxH = itData["hight"] - 1;
	//	P = 3e5;
	//	firstVehicle = true;
	//	secondVehicle = false;
	//	trigger1 = false;
	//	trigger2 = false;
	//	shutoffH = (left + right) / 2.;
	//	deltaY = -1;

	//	while (itData["hight"] > -1)
	//	{
	//		deltaY = Y(X(ro(itData["hight"]), itData["velocity"]))*cos(itData["teta"]) - X(ro(itData["hight"]), itData["velocity"])*sin(itData["teta"]) - itData["mass"] * g;

	//		if (maxHightAchieved && (fabs(itData["teta"]*180/M_PI) < eps) && (deltaY > 0))
	//			break;

	//		if (itData["hight"] < maxH)
	//			maxHightAchieved = true;
	//		else
	//			maxH = itData["hight"];

	//		// отстреливание ступеней
	//		if (!firstVehicle && !trigger1)
	//		{
	//			itData["mass"] = itData["mass"] - massVeh1;
	//			trigger1 = true;
	//		}
	//		if (!firstVehicle && !secondVehicle && trigger2)
	//		{
	//			itData["mass"] = itData["mass"] - massVeh2;
	//			trigger2 = true;
	//		}

	//		// условие накладываемое на ПВРД высотой
	//		if (secondVehicle && (itData["hight"] > 7.5e4))
	//		{
	//			P = 0;
	//			secondVehicle = false;
	//		}

	//		// отключение двигателя на высоте shuoffH
	//		if (fabs(itData["hight"] - shutoffH) < eps)
	//		{
	//			P = 0;
	//			firstVehicle = false;
	//			secondVehicle = false;
	//		}

	//		//runge-kutta
	//		itData = RungeKutta(&itData, h);
	//	}
	//	if (fabs(itData["hight"] - horizontalH) < horizontalErr)
	//		break;
	//	else
	//	{
	//		if (itData["hight"] > horizontalH)
	//			right = shutoffH;
	//		else
	//			left = shutoffH;
	//	}
	//}
//==========================================================================================================================

//==========================================================================================================================
	
	//itData["velocity"] = 1700;			// на 4 махах адекватные значения угла
	//itData["teta"] = x0;
	//itData["length"] = y0;
	//itData["hight"] = 2e4;
	////itData["mass"] = m0 - fuel01;
	//itData["mass"] = m0;
	//itData["time"] = t0;

	//firstVehicle = false;

	//constParamTrajectory(itData, shutoffH, h, fxls, &lastItter);

//==========================================================================================================================

	//itData["velocity"] = v0;
	//itData["teta"] = teta0;
	//itData["length"] = x0;
	//itData["hight"] = y0;
	//itData["mass"] = m0;
	//itData["time"] = t0;

	//P = 3e5;
	//firstVehicle = true;
	//secondVehicle = false;
	//itData = RungeKutta(&itData, h);
	//count++;

	//ballisticTrajectory(itData, shutoffH, h, fxls, &lastItter);

//==========================================================================================================================
	
	//double left = 45, right = 85;
	//double startingTeta;
	//int dihotomyCount = 0;
	//lastItter.clear();
	//while (true)
	//{
	//	startingTeta = (left + right) / 2.;

	//	itData["velocity"] = v0;
	//	itData["teta"] = startingTeta*M_PI / 180.;
	//	itData["length"] = x0;
	//	itData["hight"] = y0;
	//	itData["mass"] = m0;
	//	itData["time"] = t0;

	//	P = 3e5;
	//	firstVehicle = true;
	//	secondVehicle = false;

	//	fxls.open("output.xls", ios::out);
	//	fxls << "x\t" << "y\t" << "teta\t" << "t\t" << "v\t" << "X\t" << "Y\t" << "m\t" << "ro\t" << endl
	//		<< x0 / 1000. << "\t"
	//		<< y0 / 1000. << "\t"
	//		<< teta0 * 180 / M_PI << "\t"
	//		<< t0 << "\t"
	//		<< v0 << "\t"
	//		<< X(ro(y0), v0) << "\t"
	//		<< Y(X(ro(y0), v0)) << "\t"
	//		<< m0 << "\t"
	//		<< ro(y0) << "\t" << endl;

	//	ballisticTrajectory(itData, -1, h, fxls, &lastItter);
	//	dihotomyCount++;

	//	if (dihotomyCount >= 100)
	//		break;

	//	if (fabs(lastItter["length"] - finalPoint) < finalErr)
	//		break;
	//	else
	//	{
	//		if (lastItter["length"] > finalPoint)
	//			right = startingTeta;
	//		else
	//			left = startingTeta;
	//	}
	//}

	//std::cout << "StartingTeta = " << startingTeta << endl << endl;

//==========================================================================================================================

	double constV = 1700;
	double constH = 2e4;


	lastItter.clear();
	double maxshuoffTime;

	itData["velocity"] = constV;
	itData["teta"] = 0;
	itData["length"] = 0;
	itData["hight"] = constH;
	itData["mass"] = m0;
	itData["time"] = t0;

	firstVehicle = false;

	constParamTrajectory(itData, shutoffH, h, fxls, &maxshuoffTime, &lastItter);


	double left = 0, right = maxshuoffTime;
	double shutoffT;
	int dihotomyCount = 0;
	lastItter.clear();
	while (true)
	{
		shutoffT = (left + right) / 2.;

		itData["velocity"] = constV;
		itData["teta"] = 0;
		itData["length"] = 0;
		itData["hight"] = constH;
		itData["mass"] = m0;
		itData["time"] = t0;

		P = 3e5;
		alf = 0;
		firstVehicle = false;
		secondVehicle = true;


		double b1, b2;
		int writeCount = 1;
		int writingPeriod = 1 / h - 1;
		while (P > 0)
		{
			b1 = X(ro(itData["hight"]), itData["velocity"]);
			b2 = -Y(X(ro(itData["hight"]), itData["velocity"])) + itData["mass"] * g - itData["mass"] * itData["velocity"] * itData["velocity"] * cos(itData["teta"]) / (EarthR + itData["hight"]);
			alf = atan2(b2, b1) * 180 / M_PI;
			P = b1 / cos(alf*M_PI / 180.);

			itData = RungeKutta(&itData, h);

			if (writeCount == writingPeriod)
			{
				fxls << itData["length"] / 1000. << "\t"
					<< itData["hight"] / 1000. << "\t"
					<< itData["teta"] * 180 / M_PI << "\t"
					<< itData["time"] << "\t"
					<< itData["velocity"] << "\t"
					<< X(ro(itData["hight"]), itData["velocity"]) << "\t"
					<< Y(X(ro(itData["hight"]), itData["velocity"])) << "\t"
					<< itData["mass"] << "\t"
					<< ro(itData["hight"]) << "\t" << endl;

				writeCount = 0;
			}
			else
				writeCount++;

			if (fabs(itData["time"] - shutoffT) < h)
				break;
		}
		secondVehicle = false;
		//std::cout << "Alfa = " << alf << endl;
		alf = 0;
		P = 0;
		ballisticTrajectory(itData, shutoffH, h, fxls, &lastItter);


		dihotomyCount++;
		if (dihotomyCount >= 100)
			break;

		if (fabs(lastItter["length"] - finalPoint) < finalErr)
			break;
		else
		{
			if (lastItter["length"] > finalPoint)
				right = shutoffT;
			else
				left = shutoffT;
		}
	}

//==========================================================================================================================
	std::cout << "Hight = " << lastItter["hight"] / 1000. << " km" << endl
		<< "Length = " << lastItter["length"] / 1000. << " km" << endl
		<< "Teta = " << lastItter["teta"] *180/M_PI << endl
		<< "Mass = " << lastItter["mass"] << " kg" << endl
		<< "Time = " << lastItter["time"] << " sec" << endl
		<< "Iteration count = " << dihotomyCount << endl;

	getchar();

	fxls.close();
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
	if (((m > (m0 - (fuel01 + fuel02)))) && (P != 0))
	{
		if ((m > (m0 - fuel01)) && firstVehicle)
		{
			return (P / I1);
		}
		else
		{
			firstVehicle = false;
			secondVehicle = true;
			return (P / I2);
		}
	}
	//if ((m - (m0 - (fuel0))) > eps_)
	//	return (P / I1);
	else
	{
		P = 0;
		secondVehicle = false;
		return 0;
	}
}


double velocityFunc(double v, double m, double teta, double y)
{
	return (P*cos(M_PI*(alf + phi) / 180.) - X(ro(y), v) - m*g*sin(teta)) / m;
}
double tetaFunc(double v, double m, double teta, double y)
{
	return (P*sin(M_PI*(alf + phi) / 180.) + Y(X(ro(y), v)) - m*g*cos(teta) + (m*v*v*cos(teta) / (EarthR + y))) / (m*v);
}
double lengthFunc(double v, double teta, double y)
{
	return v*cos(teta)*(EarthR / (EarthR + y));
}
double hightFunc(double v, double teta)
{
	return v*sin(teta);
}
double massFunc(double m)
{
	return -Gc(m);
}


map<string, double> RungeKutta(const map<string, double>* data, const double step)
{
	double v, teta, x, y, m, t;
	double h = step;
	vector<double> k1(5), k2(5), k3(5), k4(5);		// {v, teta, x, y, m}
	map<string, double> result;

	v = data->at("velocity");
	teta = data->at("teta");
	x = data->at("length");
	y = data->at("hight");
	m = data->at("mass");
	t = data->at("time");

	// k1
	k1[0] = h*velocityFunc(v, m, teta, y);
	k1[1] = h*tetaFunc(v, m, teta, y);
	k1[2] = h*lengthFunc(v, teta, y);
	k1[3] = h*hightFunc(v, teta);
	k1[4] = h*massFunc(m);

	// k2
	k2[0] = h*velocityFunc(v + k1.at(0) / 2., m + k1.at(4) / 2., teta + k1.at(1) / 2., y + k1.at(3) / 2.);
	k2[1] = h*tetaFunc(v + k1.at(0) / 2., m + k1.at(4) / 2., teta + k1.at(1) / 2., y + k1.at(3) / 2.);
	k2[2] = h*lengthFunc(v + k1.at(0) / 2., teta + k1.at(1) / 2., y + k1.at(3) / 2.);
	k2[3] = h*hightFunc(v + k1.at(0) / 2., teta + k1.at(1) / 2.);
	k2[4] = h*massFunc(m + k1.at(4) / 2.);

	// k3
	k3[0] = h*velocityFunc(v + k2.at(0) / 2., m + k2.at(4) / 2., teta + k2.at(1) / 2., y + k2.at(3) / 2.);
	k3[1] = h*tetaFunc(v + k2.at(0) / 2., m + k2.at(4) / 2., teta + k2.at(1) / 2., y + k2.at(3) / 2.);
	k3[2] = h*lengthFunc(v + k2.at(0) / 2., teta + k2.at(1) / 2., y + k2.at(3) / 2.);
	k3[3] = h*hightFunc(v + k2.at(0) / 2., teta + k2.at(1) / 2.);
	k3[4] = h*massFunc(m + k2.at(4) / 2.);

	// k4
	k4[0] = h*velocityFunc(v + k3.at(0), m + k3.at(4), teta + k3.at(1), y + k3.at(3));
	k4[1] = h*tetaFunc(v + k3.at(0), m + k3.at(4), teta + k3.at(1), y + k3.at(3));
	k4[2] = h*lengthFunc(v + k3.at(0), teta + k3.at(1), y + k3.at(3));
	k4[3] = h*hightFunc(v + k3.at(0), teta + k3.at(1));
	k4[4] = h*massFunc(m + k3.at(4));


	v += (k1.at(0) + 2 * k2.at(0) + 2 * k3.at(0) + k4.at(0)) / 6.;
	teta += (k1.at(1) + 2 * k2.at(1) + 2 * k3.at(1) + k4.at(1)) / 6.;
	x += (k1.at(2) + 2 * k2.at(2) + 2 * k3.at(2) + k4.at(2)) / 6.;
	y += (k1.at(3) + 2 * k2.at(3) + 2 * k3.at(3) + k4.at(3)) / 6.;
	m += (k1.at(4) + 2 * k2.at(4) + 2 * k3.at(4) + k4.at(4)) / 6.;
	t += h;

	result.insert(make_pair("velocity", v));
	result.insert(make_pair("teta", teta));
	result.insert(make_pair("length", x));
	result.insert(make_pair("hight", y));
	result.insert(make_pair("mass", m));
	result.insert(make_pair("time", t));

	return result;
}
void ballisticTrajectory(map<string, double> itData, double shutoffH, double h, ofstream & fxls, map<string, double>* lastItter)
{
	int writeCount = 1;
	int writingPeriod = 1/h - 1;
	bool trigger1 = false, trigger2 = false;

	//runge-kutta
	itData = RungeKutta(&itData, h);

	while ((itData["hight"] + 1)> eps)
	{
		//// отстреливание ступеней
		//if (!firstVehicle && !trigger1)
		//{
		//	itData["mass"] = itData["mass"] - massVeh1;
		//	m0 -= massVeh1;
		//	trigger1 = true;
		//}
		//if (!firstVehicle && !secondVehicle && !trigger2)
		//{
		//	itData["mass"] = itData["mass"] - massVeh2;
		//	m0 -= massVeh2;
		//	trigger2 = true;
		//}

		// условие накладываемое на ПВРД высотой
		if (secondVehicle && (itData["hight"] > 7.5e4))
		{
			P = 0;
			secondVehicle = false;
		}

		// отключение двигателя на высоте shuoffH
		if ((fabs(itData["hight"] - shutoffH) < eps) && (shutoffH > 0))
		{
			P = 0;
			firstVehicle = false;
			secondVehicle = false;
		}
		
		//runge-kutta
		itData = RungeKutta(&itData, h);

		if (writeCount == writingPeriod)
		{
			fxls << itData["length"] / 1000. << "\t"
				<< itData["hight"] / 1000. << "\t"
				<< itData["teta"] * 180 / M_PI << "\t"
				<< itData["time"] << "\t"
				<< itData["velocity"] << "\t"
				<< X(ro(itData["hight"]), itData["velocity"]) << "\t"
				<< Y(X(ro(itData["hight"]), itData["velocity"])) << "\t"
				<< itData["mass"] << "\t"
				<< ro(itData["hight"]) << "\t" << endl;

			writeCount = 0;
		}
		else
			writeCount++;
	}

	fxls.close();

	lastItter->clear();
	lastItter->insert(itData.begin(), itData.end());
}
void constParamTrajectory(map<string, double> itData, double shutoffH, double h, ofstream & fxls, double* shutoffTime, map<string, double>* lastItter)
{
	double b1, b2;
	int writeCount = 1;
	int writingPeriod = 1 / h - 1;
	while (P > 0)
	{
		b1 = X(ro(itData["hight"]), itData["velocity"]);
		b2 = -Y(X(ro(itData["hight"]), itData["velocity"])) + itData["mass"] * g - itData["mass"] * itData["velocity"] * itData["velocity"] * cos(itData["teta"]) / (EarthR + itData["hight"]);
		alf = atan2(b2, b1) * 180 / M_PI;
		P = b1 / cos(alf*M_PI / 180.);

		itData = RungeKutta(&itData, h);

		if (writeCount == writingPeriod)
		{
			fxls << itData["length"] / 1000. << "\t"
				<< itData["hight"] / 1000. << "\t"
				<< itData["teta"] * 180 / M_PI << "\t"
				<< itData["time"] << "\t"
				<< itData["velocity"] << "\t"
				<< X(ro(itData["hight"]), itData["velocity"]) << "\t"
				<< Y(X(ro(itData["hight"]), itData["velocity"])) << "\t"
				<< itData["mass"] << "\t"
				<< ro(itData["hight"]) << "\t" << endl;

			writeCount = 0;
		}
		else
			writeCount++;
	}
	secondVehicle = false;
	*shutoffTime = itData["time"];
	std::cout << "Alfa = " << alf << endl;
	alf = 0;
	//itData["mass"] += massVeh1;

	ballisticTrajectory(itData, shutoffH, h, fxls, lastItter);
}