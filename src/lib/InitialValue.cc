template <int dim>
class InitialValue : public Function<dim> {
 public:
   InitialValue () : Function<dim>() {}
   virtual double value (const Point<dim>  & x , const unsigned int component) const {
      Assert (component == 0, ExcInternalError());

      double temp = 0.0;

      for (int i = 0; i < dim; ++i) {
         temp += std::exp(x[i]) / dim;
      };

      return std::max(1.0 - temp, 0.0);
   };
};
