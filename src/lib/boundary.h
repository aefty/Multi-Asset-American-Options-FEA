/**
 * boundary.h
 * =============
 * boundary conditions
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Edite by: Aryan Eftekhari & Edoardo Vecchi
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
class BoundaryValues : public Function<dim> {
  public:
	virtual double value (const Point<dim>  &p, const unsigned int component = 10) const;
};


template<int dim>
double BoundaryValues<dim>::value (const Point<dim> &p, const unsigned int component) const {
	Assert(component == 0, ExcInternalError());

	double time = this->get_time();


	if (time > .25) {
		return  1;
	} else {
		return -1;
	}
}

