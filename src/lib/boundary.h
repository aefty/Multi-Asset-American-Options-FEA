/**
 * boundary.h
 * =============
 * boundary conditions
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Edited by: Aryan Eftekhari & Edoardo Vecchi
 *
 * ---------------------------------------------------------------------
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
#include "problem/euro2d.h"

template<int dim, template <int> class PayoffClass>
class Boundary : public Function<dim> , public PayoffClass<dim> {
  public:
	Boundary() {};
	virtual double value (const Point<dim>  &S, const unsigned int component = 0) const {
		Assert(component == 0, ExcInternalError());
		double time = this->get_time();
		std::vector<double> S_temp ;

		if (dim == 3) {
			S_temp = {S[0], S[1], S[2]};
		}  else if (dim == 2) {
			S_temp = {S[0], S[1]};
		} else {
			S_temp = {S[0]};
		}

		return this->payoff_avg (S_temp, STRIKE_PRICE, time,  RFR);
	};
};


