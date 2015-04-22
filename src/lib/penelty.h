/**
 * source.h
 * =============
 * Penalty term
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Edited by: Aryan Eftekhari
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

template<int dim>
class Penelty : public Function<dim> {
  public:
	Penelty(): Function<dim>() {};

	virtual double value (const Point<dim> &p, const unsigned int component = 0) const {
		Assert (component == 0, ExcInternalError());
		Assert (dim == 2, ExcNotImplemented());

		const double time = this->get_time();

		//std::cout << "AMERICAN";

		return 10;
	};

  private:

};


// Dynamic Source
/*
template<int dim>
double RightHandSide<dim>::value (const Point<dim> &p, const unsigned int component) const {

	Assert (component == 0, ExcInternalError());
	Assert (dim == 2, ExcNotImplemented());

	const double time = this->get_time();
	const double point_within_period = (time / period - std::floor(time / period));

	if ((point_within_period >= 0.0) && (point_within_period <= 0.2)) {

		if ((p[0] > 0.5) && (p[1] > -0.5)){
			return 1;
		}else{
			return 0;
		}

	} else if ((point_within_period >= 0.5) && (point_within_period <= 0.7)) {

		if ((p[0] > -0.5) && (p[1] > 0.5)){
			return 1;
		}else{
			return 0;
		}

	} else{
		return 0;
	}
};
*/