// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEPPARAMETER_HH
#define DUNE_PDELAB_MULTISTEPPARAMETER_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Multi-step base class and multi-step parameters.
 */

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
     */
    template<typename value_type_>
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

      /** \brief Return maximal temporal order of the problems this method is appropriate for.
       */
      virtual unsigned order() const = 0;

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

    /**
     * \brief Parameter class for the BDF2 method
     * of order 1 ODE's.
     *
     * Implementation of the multi step method
     * \f[
     * u^{n+1} - \frac43 u^n + \frac13 u^{n-1} = \frac23 \Delta t f(t^{n+1}, u^{n+1})
     * \f]
     *
     * \tparam value_type C++ type of the floating point parameters
     *
     * A-stable implicit linear multistep method
     */

    template<typename value_type>
    class BDF2Scheme : public MultiStepParameterInterface<value_type>
    {
    public :

      BDF2Scheme()
      {
        A[0] = 1./3; A[1] = -4./3; A[2] = 1.0;
        B[0] = 0.0; B[1] = 0.0; B[2] = 2./3;
      }

      /** \brief Return true if method is implicit.
       */
      virtual bool implicit() const {
        return true;
      }

      /** \brief Return number of steps of the method.
       */
      virtual unsigned steps() const {
        return 2;
      }

      /** \brief Return order of the problem this method is appropriate for.
       *
       * Solves \f$\dot{u} = f(t,u)\f$ numerically
       */
      virtual unsigned order() const {
        return 1;
      }

      /** \brief Return \f$\alpha\f$ coefficients.
       *
       * Returns \f$\alpha_r\f$:
       * \f{align*}{
       * \alpha_0 &= \frac13 & \alpha_1 &= -\frac43 & \alpha_2 &= 1
       * \f}
       *
       * \note that r \f$\in 0,...,R\f$
       */
      virtual value_type alpha(int r) const {
        return A[r];
      }

      /** \brief Return \f$\beta\f$ coefficients.
       *
       * Returns \f$\beta_r\f$:
       * \f{align*}{
       * \beta_0 &= 0 & \beta_1 &= 0 & \beta_2 &= \frac23
       * \f}
       *
       * \note that r \f$\in 0,...,R\f$
       */
      virtual value_type beta(int r) const {
        return B[r];
      }

      /** \brief Return name of the scheme.
       */
      virtual std::string name() const {
        return std::string("BDF2 (Order 1)");
      }

    private :
      Dune::FieldVector<value_type,3> A;
      Dune::FieldVector<value_type,3> B;

    };

    /**
     * \brief Parameter class for the BDF3 method
     * of order 1 ODE's.
     *
     * Implementation of the multi step method
     * \f[
     * u^{n+1} - \frac{18}{11} u^n + \frac{9}{11} u^{n-1} - \frac{2}{11} u^{n-2} =
     * \frac{6}{11} \Delta t f(t^{n+1}, u^{n+1})
     * \f]
     *
     * \tparam value_type C++ type of the floating point parameters
     *
     * A(88°)-stable implicit linear multistep method
     */

    template<typename value_type>
    class BDF3Scheme : public MultiStepParameterInterface<value_type>
    {
    public :

      BDF3Scheme()
      {
        A[0] = -2./11; A[1] = 9./11; A[2] = -18./11; A[3] = 1.0;
        B[0] = 0.0; B[1] = 0.0; B[2] = 0.0; B[3] = 6./11;
      }

      /** \brief Return true if method is implicit.
       */
      virtual bool implicit() const {
        return true;
      }

      /** \brief Return number of steps of the method.
       */
      virtual unsigned steps() const {
        return 3;
      }

      /** \brief Return order of the problem this method is appropriate for.
       *
       * Solves \f$\dot{u} = f(t,u)\f$ numerically
       */
      virtual unsigned order() const {
        return 1;
      }

      /** \brief Return \f$\alpha\f$ coefficients.
       *
       * Returns \f$\alpha_r\f$:
       * \f{align*}{
       * \alpha_0 &= -\frac{2}{11} & \alpha_1 &= \frac{9}{11} & \alpha_2 &= -\frac{18}{11} & \alpha_3 &= 1
       * \f}
       *
       * \note that r \f$\in 0,...,R\f$
       */
      virtual value_type alpha(int r) const {
        return A[r];
      }

      /** \brief Return \f$\beta\f$ coefficients.
       *
       * Returns \f$\beta_r\f$:
       * \f{align*}{
       * \beta_0 &= 0 & \beta_1 &= 0 & \beta_2 &= 0 & \beta_3 &= \frac{6}{11}
       * \f}
       *
       * \note that r \f$\in 0,...,R\f$
       */
      virtual value_type beta(int r) const {
        return B[r];
      }

      /** \brief Return name of the scheme.
       */
      virtual std::string name() const {
        return std::string("BDF3 (Order 1)");
      }

    private :
      Dune::FieldVector<value_type,4> A;
      Dune::FieldVector<value_type,4> B;

    };

  } // end namespace PDELab
} // end namespace Dune
#endif
