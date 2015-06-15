
inline function PayOff(const Point<dim>  & p, double r, double t) {
	double temp = 0.0;

	for (int i = 0; i < dim; ++i) {
		temp += std::exp(p[i]) / dim;
	};

	return (std::max(1.0 - temp, 0.0) * std::exp(-1.0 * r * t));
};