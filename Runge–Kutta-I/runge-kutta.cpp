#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include "matplotpp.h"
#include "../packages/nupengl.core.0.1.0.1/build/native/include/GL/glut.h"
#include <iomanip>

long double function(long double x, long double t)
{
	//return x + exp(t);//x(t) = t exp(t)
	return x / (t + 1);
}

long double comparisonFunction(long double t)
{
	return t * exp(t);
}

typedef long double (*RealFunction)(long double);

struct Point
{
	Point(long double xt, long double t)
	{
		this->xt = xt;
		this->t = t;
	}

	long double xt;
	long double t;
};



std::vector<Point*> rungeKuttaI(long double t0, long double tn, long double x0, long double h)
{
	long double k1;

	if (h < 1e-5)
		h = 1e-5;

	if (h > 0.01)
		h = 0.01;

	int n = fabs(tn - t0) / h;

	std::vector<Point*> Xt;
	Xt.push_back(new Point(x0, t0));//x(t0)=x0

	for (int i = 1; i <= n; i++)
	{
		k1 = h * function(Xt[i - 1]->xt, Xt[i - 1]->t);
		Xt.push_back(new Point(Xt[i - 1]->xt + k1, t0 + i * h));
	}
	return Xt;
}


std::vector<Point*> rungeKuttaI(long double t0, long double tn, long double x0, long double epsilon, long double LipschitzConstant, long double dfdt_t0_x0, long double dfdx_t0_x0)
{
	//http://isites.harvard.edu/fs/docs/icb.topic250025.files/NumericalTheor.pdf
	//|GTE| = epsilon <= (hM/2L) (exp(L(t-t0))-1) 
	//epsilon <= (hM/2L) (exp(L(B-A))-1)
	// h >= (2L epsilon )/( M (exp(L(B-A))-1) )
	long double M = dfdt_t0_x0+ dfdx_t0_x0*function(x0,t0);//upper bound on the second derivative of x(t) on the given interval
	long double L = LipschitzConstant;//the Lipschitz constant of f
	long double h = (2L*epsilon) / (M*(exp(L*(tn - t0)) - 1));
	return rungeKuttaI(t0, tn, x0, h);
}


std::vector<Point*> rungeKuttaI(long double t0, long double tn, long double x0, long double epsilon, long double LipschitzConstant, long double upperBoundOnTheSecondDerivativeOfX)
{
	//http://isites.harvard.edu/fs/docs/icb.topic250025.files/NumericalTheor.pdf
	//|GTE| = epsilon <= (hM/2L) (exp(L(t-t0))-1) 
	//epsilon <= (hM/2L) (exp(L(B-A))-1)
	// h >= (2L epsilon )/( M (exp(L(B-A))-1) )
	long double M = upperBoundOnTheSecondDerivativeOfX;//upper bound on the second derivative of x(t) on the given interval
	long double L = LipschitzConstant;//the Lipschitz constant of f
	long double h = (2L*epsilon) / (M*(exp(L*(tn - t0)) - 1));
	return rungeKuttaI(t0, tn, x0, h);
}

std::string toComputatorNET(std::vector<Point*> xt)
{
	std::string retXt = "var xt = array(";

	std::string retT = "var t = array(";

	for (int i = 0; i < xt.size(); i++)
	{
		retT += std::to_string(xt[i]->t) + ",";
		retXt += std::to_string(xt[i]->xt) + ",";
	}
	retT[retT.length() - 1] = ')';
	retXt[retXt.length() - 1] = ')';

	std::string ret = retT + ";\n" + retXt + ";\n" + "plot(t,xt);\nplot((real x) => x·exp(x),0,1,0,3);";
	return ret;
}


class ComparisonPlot : public MatPlot
{
public:
	ComparisonPlot(): comparisonFunction(nullptr)
	{
	}

	ComparisonPlot(std::vector<Point*> points, RealFunction comparisonFunction)
	{
		this->points = points;
		this->comparisonFunction = comparisonFunction;
	}


private:
	RealFunction comparisonFunction;
	std::vector<Point*> points;

	void DISPLAY() override
	{
		std::vector<double> xNumerical, xAnalytic, yNumerical, yAnalytic;
		double ymax=-1e32, ymin=1e32;
		for (int i = 0; i < points.size(); ++i)
		{
			if (points[i]->xt < ymin)
				ymin = points[i]->xt;
			else if(points[i]->xt > ymax)
				ymax = points[i]->xt;
			xNumerical.push_back(points[i]->t);
			yNumerical.push_back(points[i]->xt);
		}

		long double A = points[0]->t;
		long double B = points[points.size() - 1]->t;
		int N = 100000;
		for (int i = 0; i <= N; i++)
		{
			xAnalytic.push_back(A + i * (fabs(B - A) / N));
			yAnalytic.push_back(this->comparisonFunction(xAnalytic[i]));
		}

		title("Porownanie: czerwony - rozwiazanie numeryczne, zielony - rozwiazanie analityczne");
		axis(floor(A), ceil(B), floor(ymin)-1, ceil(ymax)+1);
		
		plot(xNumerical, yNumerical);
		set("r");
		//plot(xAnalytic, yAnalytic);
		//set("g");
	}
};

ComparisonPlot* plot;

void showChart(int argc, char* argv[], std::vector<Point*> xt)
{
	glutInit(&argc, argv);
	glutCreateWindow(100, 100, 800, 700, "Wykres - porownanie rozwiazania analitycznego z numerycznym");

	plot = new ComparisonPlot(xt, comparisonFunction);

	glutDisplayFunc([]() -> void
		{
			plot->display();
		});
	glutReshapeFunc([](int w, int h) -> void
		{
			plot->reshape(w, h);
		});
	glutMainLoop();
}

void saveToFile(std::vector<Point*> xt)
{
	std::setprecision(16);
	std::ofstream out("output.txt");
	out.precision(16);
	out << toComputatorNET(xt);
	out.close();
}

int main(int argc, char* argv[])
{
	long double epsilon, x0;
	long double LipschitzConstant, dfdt_t0_x0, dfdx_t0_x0,M;


	std::cin >> x0 >> epsilon >> LipschitzConstant >> dfdt_t0_x0 >> dfdx_t0_x0;
	std::vector<Point*> xt = rungeKuttaI(0, 1, x0, epsilon, LipschitzConstant, dfdt_t0_x0, dfdx_t0_x0);

	//std::cin >> x0 >> epsilon >> LipschitzConstant >> M;
	//std::vector<Point*> xt = rungeKuttaI(0, 1, x0, epsilon, LipschitzConstant, M);

	//std::cin >> x0 >> epsilon;
	//std::vector<Point*> xt = rungeKuttaI(0, 1, x0, epsilon);

	saveToFile(xt);

	showChart(argc, argv, xt);

	return 0;
}

