#include "eigen3/Eigen/Eigenvalues"
class ExpoTrigFunction {
public:
	double E1, C1, B1, A1, w1, cos1, sin1, E2, C2, B2, A2, w2, cos2, sin2;
	
	ExpoTrigFunction(double E1, double C1, double B1, double A1, double w1, double cos1, double sin1, double E2, double C2, double B2, double A2, double w2, double cos2, double sin2) : E1(E1), C1(C1), B1(B1), A1(A1), w1(w1), cos1(cos1), sin1(sin1), E2(E2), C2(C2), B2(B2), A2(A2), w2(w2), cos2(cos2), sin2(sin2) {}
	
	ExpoTrigFunction() : E1{}, C1{}, B1{}, A1{}, w1{}, cos1{}, sin1{}, E2{}, C2{}, B2{}, A2{}, w2{}, cos2{}, sin2{} {}

	double operator()(double t) {
		return std::exp(E1*t)*(C1+B1*t+A1*t*t+cos1*std::cos(w1*t)+sin1*std::sin(w1*t))
				+ std::exp(E2*t)*(C2+B2*t+A2*t*t+cos2*std::cos(w2*t)+sin2*std::sin(w2*t));
	}
};

class Solver {
public:
	int n;
	Eigen::MatrixXd coeff;
	double damp;
	Eigen::VectorXd constfac;
	std::vector<ExpoTrigFunction> vz;
	std::vector<ExpoTrigFunction> z;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigens;
	Eigen::VectorXd xC;

	Solver() {}
	Solver(int n, const Eigen::MatrixXd &coeff, double damp, const Eigen::VectorXd &constfac) : n(n), coeff(coeff), damp(damp), constfac(constfac), vz(n), z(n), eigens(coeff), xC(n) {}

	void solve(const Eigen::VectorXd &x0, const Eigen::VectorXd &vx0, double t0=0.0) {
		double a = -damp; double s = a/2;
		Eigen::VectorXd ax0 = coeff * x0 + a * vx0 + constfac;
		
		Eigen::MatrixXd coeffnew = coeff + Eigen::MatrixXd::Identity(n, n)*(a*a)/4;

		double exp0 = std::exp(t0*a/2);
		Eigen::VectorXd vz0 = eigens.eigenvectors().transpose() * vx0 / exp0;
		Eigen::VectorXd az0 = eigens.eigenvectors().transpose() * (ax0 - vx0*a/2) / exp0;

		for (int i = 0; i < n; i++) {
			double K = eigens.eigenvalues()(i) + (a*a)/4;

			if (K < 0) {
				double omega = std::sqrt(-K);

				double cosomegat0 = std::cos(omega*t0);
				double sinomegat0 = std::sin(omega*t0);

				double C2 = sinomegat0 * vz0(i) + cosomegat0 * az0(i) / omega;
				double C1 = cosomegat0 * vz0(i) - sinomegat0 * az0(i) / omega;

				vz[i] = ExpoTrigFunction(0, 0, 0, 0, omega, C1, C2, 0, 0, 0, 0, 0, 0, 0);

				if (a == 0 && omega == 0)
					z[i] = ExpoTrigFunction(0, 0, C1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
				else
					z[i] = ExpoTrigFunction(s, 0, 0, 0, omega, (C1*s-C2*omega)/(s*s+omega*omega), (C1*omega+C2*s)/(s*s+omega*omega), 0, 0, 0, 0, 0, 0, 0);
			} else if (K > 0) {
				double omega = std::sqrt(K);
				
				double e1 = std::exp(omega*t0);
				double e2 = std::exp(-omega*t0);

				double C2 = (omega * vz0(i) - az0(i)) / (2 * omega * e2);
				double C1 = (omega * vz0(i) + az0(i)) / (2 * omega * e1);

				vz[i] = ExpoTrigFunction(omega, C1, 0, 0, 0, 0, 0, -omega, C2, 0, 0, 0, 0, 0);
				
				if ((a/2 + omega) == 0 && (a/2 - omega) == 0)
					z[i] = ExpoTrigFunction(0, 0, C1+C2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
				else if ((a/2 + omega) == 0)
					z[i] = ExpoTrigFunction(0, 0, C1, 0, 0, 0, 0, a/2 - omega, C2/(a/2 - omega), 0, 0, 0, 0, 0);
				else if ((a/2 - omega) == 0)
					z[i] = ExpoTrigFunction(a/2 + omega, C1/(a/2 + omega), 0, 0, 0, 0, 0, 0, 0, C2, 0, 0, 0, 0);
				else
					z[i] = ExpoTrigFunction(a/2 + omega, C1/(a/2 + omega), 0, 0, 0, 0, 0, a/2 - omega, C2/(a/2 - omega), 0, 0, 0, 0, 0);
			} else {
				vz[i] = ExpoTrigFunction(0, vz0(i), az0(i), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
				
				if (a == 0)
					z[i] = ExpoTrigFunction(0, 0, vz0(i), az0(i)/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
				else
					z[i] = ExpoTrigFunction(s, vz0(i)*2/a - az0(i)*4/(a*a), az0(i)*2/a, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			}
		}

		Eigen::VectorXd zt0withoutconst(n);
		for (int i = 0; i < n; i++)
			zt0withoutconst(i) = z[i](t0);

		xC = x0 - eigens.eigenvectors() * zt0withoutconst;
	}

	std::pair<Eigen::VectorXd, Eigen::VectorXd> soln(double t) {
		Eigen::VectorXd zt(n);
		Eigen::VectorXd vzt(n);
		for (int i = 0; i < n; i++) {
			zt(i) = z[i](t);
			vzt(i) = vz[i](t);
		}

		return std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(eigens.eigenvectors() * zt + xC, std::exp(-damp*t/2) * (eigens.eigenvectors() * vzt));
	}
};
