/**
 * Transform.h
 * =============
 * Penalty term
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
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
class Transform : public Function<dim> {
	/**
	 * Variable Definition
	 */
	std::vector<double> tempVec;

  public:
	Transform() {};

	virtual double S2X (std::vector<double> S) const {

		tempVec.clear();

		for (int i = 0; i < dim; ++i) {
			tempVec.push_back(log(S[i]));
		}

		return tempVec;
	};

	virtual double X2S (std::vector<double> X) const {

		tempVec.clear();

		for (int i = 0; i < dim; ++i) {
			tempVec.push_back(exp (X[i]));
		}

		return tempVec;
	};

  private:
};