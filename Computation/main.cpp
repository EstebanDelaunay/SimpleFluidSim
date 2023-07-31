#include <iostream>

using std::cout;
using std::endl;
using std::begin;
using std::end;

int main() 
{
	// Paramètres physiques
	constexpr double radius(0.5);
	constexpr double Re(200);
	constexpr double pi(3.14159);

	// Maillage espace
	constexpr int sizeR(100);
	constexpr int sizeTheta(150);

	// Paramètres de calcul
	constexpr int kMax(1e3);
	constexpr double eps(1e-5);

	// Maillage en coordonnée polaire
	constexpr double lengthR = 19.0 * radius;
	constexpr double stepSizeR = lengthR / (sizeR - 1.0);
	double r[sizeR]{0.0};
	for (double* p = r+1; p != end(r); p++) { *p = *(p - 1) + stepSizeR; }

	constexpr double lengthTheta = 2.0 * pi;
	constexpr double stepSizeTheta = lengthTheta / (sizeTheta - 1.0);
	double theta[sizeTheta]{ 0.0 };
	for (double* p = theta + 1; p != end(theta); p++) { *p = *(p - 1) + stepSizeTheta; }


	cout << "Done" << endl;

	return 0;
}