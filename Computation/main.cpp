#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::begin;
using std::end;
using std::sin;

int main()
{
	// Paramètres physiques
	constexpr double radius = 0.5, Re = 200.0, uInf = 1.0;
	constexpr double pi(3.14159);

	// Maillage espace
	constexpr size_t sizeR(40), sizeTheta(40);

	// Paramètres de calcul
	constexpr int kMax(100000);
	constexpr double eps(1e-5);

	// Maillage en coordonnée polaire
	constexpr double lengthR = 19.0 * radius;
	constexpr double stepSizeR = lengthR / (sizeR - 1.0);
	double r[sizeR]{ radius };
	for (double* p = r + 1; p != end(r); p++) { *p = *(p - 1) + stepSizeR; }

	constexpr double lengthTheta = 2.0 * pi;
	constexpr double stepSizeTheta = lengthTheta / (sizeTheta - 1.0);
	double theta[sizeTheta]{ 0.0 };
	for (double* p = theta + 1; p != end(theta); p++) { *p = *(p - 1) + stepSizeTheta; }


	double psi[sizeR][sizeTheta]{ 0.0 };
	double psiOld[sizeR][sizeTheta]{ 0.0 };

	// Conditions aux bords
	for (size_t i = 0; i != sizeR; ++i) {
		psi[i][0] = 20.0;
		psi[i][sizeTheta - 1] = 20.0;
	}

	for (size_t j = 0; j != sizeTheta; ++j) {
		psi[0][j] = 20;
		psi[sizeR - 1][j] = uInf * 19.0 * radius * sin(theta[j]) + 20.0;
	}

	int k = 0;
	double error = 1.0;
	while (k < kMax && error > eps)
	{
		for (size_t i = 0; i != sizeR; ++i) {
			for (size_t j = 0; j != sizeTheta; ++j) {
				psiOld[i][j] = psi[i][j];
			}
		}

		for (size_t i = 1; i != sizeR-1; ++i) {
			for (size_t j = 1; j != sizeTheta-1; ++j) {
				psi[i][j] = ((1.0 / (stepSizeR * stepSizeR) + 1.0 / (2.0 * r[i] * stepSizeR)) * psiOld[i + 1][j] + (1.0 / (stepSizeR * stepSizeR) - 1.0 / (2.0 * r[i] * stepSizeR)) * psiOld[i - 1][j] + (psiOld[i][j + 1] + psiOld[i][j - 1]) / ((r[i] * r[i]) * (stepSizeTheta * stepSizeTheta))) * 1.0 / (2.0 / (stepSizeR * stepSizeR) + 2.0 / ((r[i] * r[i]) * (stepSizeTheta * stepSizeTheta)));
			}
			psi[i][0] =     ((1.0 / (stepSizeR * stepSizeR) + 1.0 / (2.0 * r[i] * stepSizeR)) * psiOld[i + 1][0] + (1.0 / (stepSizeR * stepSizeR) - 1.0 / (2.0 * r[i] * stepSizeR)) * psiOld[i - 1][0] + (psiOld[i][1] + psiOld[i][sizeR - 2]) / ((r[i] * r[i]) * (stepSizeTheta * stepSizeTheta))) * 1.0 / (2.0 / (stepSizeR * stepSizeR) + 2.0 / ((r[i] * r[i]) * (stepSizeTheta * stepSizeTheta)));
			psi[i][sizeTheta-1] = psi[i][0];
		}

		error = 0.0;
		for (size_t i = 0; i != sizeR; ++i) {
			for (size_t j = 0; j != sizeTheta; ++j) {
				error += psi[i][j] - psiOld[i][j];
			}
		}
		//cout << error << " Iterations " << k++ << endl;
	}

	cout << "Done" << endl;

	return 0;
}