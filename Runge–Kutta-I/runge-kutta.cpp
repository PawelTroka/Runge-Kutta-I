#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include "matplotpp.h"
#include "../packages/nupengl.core.0.1.0.1/build/native/include/GL/glut.h"
#include <iomanip>

long double function(long double x, long double t)
{
	return x + exp(t);//x(t) = t exp(t)
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

std::vector<Point*> rungeKuttaI(long double A, long double B, long double x0, long double epsilon)
{
	long double h = epsilon, t0 = A;
	long double k1;
	int n = fabs(B - A) / h;

	std::vector<Point*> Xt;
	Xt.push_back(new Point(x0, t0));//x(t0)=x0

	for (int i = 1; i <= n; i++)
	{
		k1 = h * function(Xt[i - 1]->xt, Xt[i - 1]->t);
		Xt.push_back(new Point(Xt[i - 1]->xt + k1, t0 + i * h));
	}
	return Xt;
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

		for (int i = 0; i < points.size(); ++i)
		{
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

		plot(xNumerical, yNumerical);
		set("r");
		plot(xAnalytic, yAnalytic);
		set("g");
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

	std::cin >> x0 >> epsilon;

	std::vector<Point*> xt = rungeKuttaI(0, 1, x0, epsilon);

	saveToFile(xt);

	showChart(argc, argv, xt);

	return 0;
}

