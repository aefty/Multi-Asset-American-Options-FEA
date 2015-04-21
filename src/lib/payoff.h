/**
 * initalCondition.h
 * =============
 * Inital Conditions term
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Author: Aryan Eftekhari & Edoardo Vecchi
 */


inline double Payoff_average (std::vector<double> &S, double strike_price, double time, const double rfr) {
	double discount  = std::exp(-1. * rfr * time);
	double sum = 0.;

	for (int i = 0; i < S.size(); ++i) {
		sum += S[i];
	};

	double pay = (sum / 2.0 - strike_price) * discount;

	std::cout << pay;

	return std::max(pay , 0.0);
}
