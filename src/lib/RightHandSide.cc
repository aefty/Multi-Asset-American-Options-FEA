template<int dim>
class RightHandSide : public Function<dim> {
 public:
  RightHandSide ()
    :
    Function<dim>(),
    period (0.1)
  {}

  virtual double value (const Point<dim> &p, const unsigned int component = 0) const;

 private:
  const double period;
};





template<int dim>
double RightHandSide<dim>::value (const Point<dim> &p, const unsigned int component) const {
  Assert (component == 0, ExcInternalError());
  Assert (dim == 2, ExcNotImplemented());

  const double time = this->get_time();
  const double point_within_period = (time / period - std::floor(time / period));

  if ((point_within_period >= 0.0) && (point_within_period <= 0.2)) {
    if ((p[0] > 0.5) && (p[1] > -0.5))
    { return 1; }
    else
    { return 0; }
  } else if ((point_within_period >= 0.5) && (point_within_period <= 0.7)) {
    if ((p[0] > -0.5) && (p[1] > 0.5))
    { return 1; }
    else
    { return 0; }
  } else
  { return 0; }
};
