/* ---------------------------------------------------------------------
 * Edited By: Aryan Eftekhari & Edoardo Vecchi 2015 (Based on modification done by Patrick Sanan, May 2015 )
 * Università della Svizzera italiana
 *
 * Copyright (C) 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2013
 */

template <int dim>
class IntrinsicValue : public Function<dim> {
  public:
	IntrinsicValue () : Function<dim>() {}
	virtual double value (const Point<dim>  & x , const unsigned int component) const {
		Assert (component == 0, ExcInternalError());

		double temp = 0.0;

		for (int i = 0; i < dim; ++i) {
			temp += std::exp(x[i]) / dim;
		};

		return std::max(1.0 - temp, 0.0);
	};
};
