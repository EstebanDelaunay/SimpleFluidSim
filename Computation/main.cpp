#include <cmath>
#include <fstream>

using std::pow;

int main()
{
	// Ouverture du fichier pour export des données
	std::ofstream expFile("data", std::ofstream::out);

	// Paramètres physiques
	constexpr double radius = 0.5, Re = 200.0, uInf = 1.0;
	constexpr double pi(3.141592653589793);

	// Paramètres temps
	int timeCurrent(0);
	constexpr int timeMax(10000);
	const double stepSizeTime(1e-2), & dt = stepSizeTime;

	// Maillage espace
	constexpr size_t sizeR(50), sizeTheta(75);

	// Paramètres de calcul
	constexpr int kMax(10000);
	constexpr double eps(1e-5), om_psi(1.9);

	// Maillage en coordonnée polaire
	constexpr double lengthR = 19.0 * radius;
	const double stepSizeR = lengthR / (sizeR - 1.0), & dr = stepSizeR;
	double r[sizeR]{ radius };
	for (double* p = r + 1; p != std::end(r); p++) { *p = *(p - 1) + stepSizeR; }

	constexpr double lengthTheta = 2.0 * pi;
	const double stepSizeTheta = lengthTheta / (sizeTheta - 1.0), & dth = stepSizeTheta;
	double theta[sizeTheta]{ 0.0 };
	for (double* p = theta + 1; p != std::end(theta); p++) { *p = *(p - 1) + stepSizeTheta; }

	// Initialisation des matrices
	double psi[sizeR][sizeTheta]{ 0.0 };
	double psiOld[sizeR][sizeTheta]{ 0.0 };

	double w[sizeR][sizeTheta]{ 0.0 };
	double wOld[sizeR][sizeTheta]{ 0.0 };

	// Conditions aux bords
	for (size_t j = 0; j != sizeTheta; ++j) {
		psi[sizeR-1][j] = r[sizeR-1] * std::sin(theta[j]);
	}

	// Initialisation des constantes de calcul
	const double dr2 = pow(dr, 2), dth2 = pow(dth, 2);
	double g(0.0), f(0.0);

	double a[sizeR]{}, b[sizeR]{}, c[sizeR]{}, d[sizeR]{};
	for (size_t i = 0; i != sizeR; i++) {
		a[i] = 2.0 * (1.0 / dr2 + 1.0 / (pow(r[i],2.0) * dth2));
		b[i] = 1.0 / dr2 + 1.0 / (2.0 * r[i] * dr);
		c[i] = 1.0 / dr2 - 1.0 / (2.0 * r[i] * dr);
		d[i] = 1.0 / (pow(r[i], 2.0) * dth2);
	}

	// Calcul
	while (timeCurrent != timeMax) {
		timeCurrent++;
		
		for (size_t i = 0; i != sizeR; ++i) {
			for (size_t j = 0; j != sizeTheta; ++j) {
				wOld[i][j] = w[i][j];
			}
		}

		int k = 0;
		double error = 1.0;
		while (k != kMax && error > eps)
		{
			for (size_t i = 0; i != sizeR; ++i) {
				for (size_t j = 0; j != sizeTheta; ++j) {
					psiOld[i][j] = psi[i][j];
				}
			}

			for (size_t i = 1; i != sizeR - 1; ++i) {
				for (size_t j = 1; j != sizeTheta - 1; ++j) {
					psi[i][j] = (1.0 - om_psi) * psiOld[i][j] + (om_psi / a[i]) * (b[i] * psi[i + 1][j] + c[i] * psi[i - 1][j] + d[i] * (psi[i][j + 1] + psi[i][j - 1]        ) + w[i][j]);
				}
				psi[i][0] =     (1.0 - om_psi) * psiOld[i][0] + (om_psi / a[i]) * (b[i] * psi[i + 1][0] + c[i] * psi[i - 1][0] + d[i] * (psi[i][1]     + psi[i][sizeTheta - 2]) + w[i][0]);
				psi[i][sizeTheta - 1] = psi[i][0];

			}

			error = 0.0;
			for (size_t i = 0; i != sizeR; ++i) {
				for (size_t j = 0; j != sizeTheta; ++j) {
					error += psi[i][j] - psiOld[i][j];
				}
			}
		}

		
				
		for (size_t j = 0; j != sizeTheta; ++j) {
			w[0][j] = (psi[2][j] - 8.0 * psi[1][j]) / (2.0 * dr2);
		}

		for (size_t i = 1; i != sizeR - 1; ++i) {
			for (size_t j = 1; j != sizeTheta - 1; ++j) {
				g = (psi[i + 1][j] - psi[i - 1][j]) * (w[i][j + 1] - w[i][j - 1]) - (psi[i][j + 1] - psi[i][j - 1]) * (w[i + 1][j] - w[i - 1][j]);
				f = g / (4.0 * r[i] * dr * dth) + (2.0 / Re) * (b[i] * w[i + 1][j] + c[i] * w[i - 1][j] + d[i] * (w[i][j + 1] + w[i][j - 1]) - a[i] * w[i][j]);
				w[i][j] = wOld[i][j] + dt * f;
			}

			g = (psi[i + 1][0] - psi[i - 1][0]) * (w[i][1] - w[i][sizeTheta - 2]) - (psi[i][1] - psi[i][sizeTheta - 2]) * (w[i + 1][0] - w[i - 1][0]);
			f = g / (4.0 * r[i] * dr * dth) + (2.0 / Re) * (b[i] * w[i + 1][0] + c[i] * w[i - 1][0] + d[i] * (w[i][1] + w[i][sizeTheta - 2]) - a[i] * w[i][0]);
			w[i][0] = wOld[i][0] + dt * f;
			w[i][sizeTheta - 1] = w[i][0];

		}
		
	}

	// Exportation des données
	for (size_t i = 0; i != sizeR; ++i) {
		for (size_t j = 0; j != sizeTheta; ++j) {
			expFile << psi[i][j] << " ";
		}
		expFile << std::endl;
	}

	return 0;
}