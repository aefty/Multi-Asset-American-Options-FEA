/**
 * initalCondition.h
 * =============
 * Initial Conditions term
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Author: Aryan Eftekhari
 */

template<int dim>
class Payoff {

  public:
	Payoff() {};

	/**
	 * Average Payoff function
	 * @param  S            [description]
	 * @param  strike_price [description]
	 * @param  time         [description]
	 * @param  rfr          [description]
	 * @return              [description]
	 */
	virtual double payoff_avg (std::vector<double> S, double strike_price, double time, const double rfr) const {
		double discount  = std::exp(-1. * rfr * time);
		double sum = 0.;

		for (int i = 0; i < S.size(); ++i) {
			sum += S[i];
		};

		double pay = (sum / 2.0 - strike_price) * discount;
		return std::max(pay , 0.0);
	};

	virtual double payoff_test() const {
		return 123.0;
	};

  private:

};
