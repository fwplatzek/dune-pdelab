// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEPPARAMETER_HH
#define DUNE_PDELAB_MULTISTEPPARAMETER_HH

#include <string>

namespace Dune {
  namespace PDELab {

    /** \brief Base parameter class for multistep time stepping parameters.
     *
     * Abstract base class for the implementation of multistep method
     * \f[
     * \sum_{r=0}^R \alpha_{R-r} u^{n+1-r} = \Delta t \sum_{r=0}^R \beta_{R-r} f(t^{n+1-r}, u^{n+1-r})
     * \f]
     *
     * \tparam value_type C++ type of the floating point parameters
     * \tparam order_ Order of the ODE this scheme is appropriate for.
     */
    template<typename value_type_, unsigned order_>
    class MultiStepParameterInterface
    {
    public:
      /** \brief Export type of the parameters
       */
      typedef value_type_ RealType;

      /** \brief Return true if method is implicit.
       */
      virtual bool implicit() const = 0;

      /** \brief Return number of steps of the method.
       */
      virtual unsigned steps() const = 0;

      /** \brief Temporal order of the problems this methods is appropriate for
       */
      static const unsigned order = order_;

      /** \brief Return \f$\alpha\f$ coefficients.
       * \note that r \f$\in 0,...,R\f$
       */
      virtual RealType alpha(int r) const = 0;

      /** \brief Return \f$\beta\f$ coefficients.
       * \note that r \f$\in 0,...,R\f$
       */
      virtual RealType beta(int r) const = 0;

      /** \brief Return name of the scheme.
       */
      virtual std::string name() const = 0;

      /** \brief every abstract base class has a virtual destructor.
       */
      virtual ~MultiStepParameterInterface() {}

    }; // end class MultiStepParameterInterface

  } // end namespace PDELab
} // end namespace Dune
#endif
