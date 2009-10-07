// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_ELECTRODYNAMIC_HH
#define DUNE_PDELAB_ELECTRODYNAMIC_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/typetraits.hh>
#include<dune/grid/common/quadraturerules.hh>

#include"../common/geometrywrapper.hh"
#include"../gridoperatorspace/gridoperatorspace.hh"
#include"../gridoperatorspace/gridoperatorspaceutilities.hh"
#include"pattern.hh"
#include"flags.hh"

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    //! A local operator for a simple electrodynamic wave equation
    /**
     * Solve the equation
     * \f{align*}{
     *    \nabla\times\mu^{-1}\nabla\times\mathbf E
     *    +\epsilon\partial_t^2\mathbf E &= -\partial_t\mathbf J
     *    \qquad\text{in }\Omega \\
     *    \mathbf{\hat n}\times\mathbf E &= 0
     *    \qquad\text{on }\partial\Omega
     * \f}
     * (i.e. \f$\sigma=0\f$, impressed current \f$\mathbf J\f$ and PEC
     * boundary conditions
     * everywhere).  This leads to the following weak formulation:
     *
     * \f[
     *    \int_\Omega\{
     *       \mu^{-1}(\nabla\times\mathbf v)\cdot(\nabla\times\mathbf E)
     *       +\mathbf v\cdot\epsilon\partial_t^2\mathbf E
     *       +\mathbf v\cdot\partial_t\mathbf J
     *    \}dV = 0\qquad \forall \mathbf v
     * \f]
     * The boundary integral from the greens theorem vanishes because we will
     * impress dirichlet boundarycondition everywhere.  Inserting a base
     * \f$\{\mathbf N_i\}\f$ into \f$\mathbf v\f$ and \f$\mathbf E\f$:
     * \f[
     *    \sum_jT_{ij}(\partial_t^2u_j)+\sum_jS_{ij}u_j+f_i=0 \qquad\forall j
     * \f]
     * with
     * \f{align*}{
     *    T_{ij}&=\int_\Omega\epsilon\mathbf N_i\cdot\mathbf N_jdV \\
     *    S_{ij}&=\int_\Omega\mu^{-1}(\nabla\times\mathbf N_i)\cdot
     *                               (\nabla\times\mathbf N_j)dV   \\
     *    f_i&=\int_\Omega\mathbf N_i\cdot\partial_t\mathbf JdV    \\
     *    \mathbf E &= \sum_iu_i\mathbf N_i
     * \f}
     *
     * The time scheme is central differences from Jin (12.29).  With the
     * simplifications we have done above it looks like
     * \f[
     *    Tu^{n+1}=2Tu^n-Tu^{n-1}-(\Delta t)^2Su^n-(\Delta t)^2f^n = 0
     * \f]
     * Bringing this into the residual formulation we get
     * \f[
     *    r=T(u^{n+1}-2u^n+u^{n-1})+(\Delta t)^2Su^n+(\Delta t)^2f^n
     * \f]
     *
     * \note Currently \f$\epsilon\f$ is fixed to 1.
     *
     * \tparam Eps    Type of function to evaluate \f$\epsilon\f$
     * \tparam Mu     Type of function to evaluate \f$\mu\f$
     * \tparam DtJ    Type of function to evaluate \f$\partial_t\mathbf J\f$
     * \tparam GCV    Type of the global coefficient vector, used for storing
     *                pointers to \f$u^n\f$ and \f$u^{n-1}\f$
     */
    template<typename Eps, typename Mu, typename DtJ, typename GCV>
	class Electrodynamic
      : public NumericalJacobianApplyVolume<Electrodynamic<Eps, Mu, DtJ, GCV> >
      , public NumericalJacobianVolume<Electrodynamic<Eps, Mu, DtJ, GCV> >
      , public NumericalJacobianApplyBoundary<Electrodynamic<Eps, Mu, DtJ,
                                                             GCV> >
      , public NumericalJacobianBoundary<Electrodynamic<Eps, Mu, DtJ, GCV> >
      , public FullVolumePattern
      , public LocalOperatorDefaultFlags
	{
      //! size of FieldVector for holding the curl
      template <unsigned d>
      struct CurlTraits {
        static const unsigned dim =
          d == 1 ? 2 :
          d == 2 ? 1 :
          /*else*/ 3;
      };

      template<typename RF>
      static void jacobianToCurl(FieldVector<RF, 1> &curl,
                                 const FieldVector<FieldVector<RF, 2>, 2> &jacobian)
      {
        curl[0] = jacobian[1][0] - jacobian[0][1];
      }
      template<typename RF>
      static void jacobianToCurl(FieldVector<RF, 3> &curl,
                                 const FieldVector<FieldVector<RF, 3>, 3> &jacobian)
      {
        for(unsigned i = 0; i < 3; ++i)
          curl[i] = jacobian[(i+2)%3][(i+1)%3] - jacobian[(i+1)%3][(i+2)%3];
      }
	public:

      // pattern assembly flags
      enum { doPatternVolume = true };

	  // residual assembly flags
      enum { doAlphaVolume = true };

      //! Construct an Electrodynamic localoperator
      /**
       * \param eps_     Reference to function object to evaluate
       *                 \f$\epsilon\f$.
       * \param mu_      Reference to function object to evaluate \f$\mu\f$.
       * \param dtJ_     Reference to function object to evaluate
       *                 \f$\partial_t\mathbf J\f$, the time derivative of the
       *                 impressed current.
       * \param Delta_t_ Size of time step.
       * \param qorder_  Quadrature order to use.
       *
       * \note The references the the function objects should be valid for as
       *       long as this localoperators residual() method is used.
       */
      Electrodynamic(const Eps &eps_, const Mu& mu_, const DtJ& dtJ_,
                     double Delta_t_ = 0, int qorder_ = 2)
        : eps(eps_)
        , mu(mu_)
        , dtJ(dtJ_)
        , Ecur(0)
        , Eprev(0)
        , Delta_t(Delta_t_)
        , qorder(qorder_)
      {}

	  //! volume integral depending on test and ansatz functions
      /**
       * We support only Galerkin method lfsu==lfsv
       */
	  template<typename EG, typename LFS, typename X, typename R>
	  void alpha_volume (const EG& eg, const LFS& lfsu, const X& x, const LFS& lfsv, R& r) const
	  {
		dune_static_assert(LFS::Traits::LocalFiniteElementType::
                           Traits::LocalBasisType::Traits::dimRange == 3 ||
                           LFS::Traits::LocalFiniteElementType::
                           Traits::LocalBasisType::Traits::dimRange == 2,
                           "Works only in 2D or 3D");

		// domain and range field type
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainFieldType DF;
        dune_static_assert((Dune::is_same<typename EG::Geometry::ctype, DF>::value),
                           "Grids ctype and Finite Elements DomainFieldType must match");
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::DomainType D;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::RangeType RangeType;
		static const unsigned dimRange = LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::dimRange;
		typedef typename LFS::Traits::LocalFiniteElementType::
		  Traits::LocalBasisType::Traits::JacobianType JacobianType;
        typedef FieldVector<RF, CurlTraits<dimRange>::dim> CurlType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        std::vector<std::vector<DF> >
          T(lfsu.size(), std::vector<DF>(lfsu.size(), 0));
        std::vector<std::vector<DF> >
          S(lfsu.size(), std::vector<DF>(lfsu.size(), 0));
        std::vector<DF>
          f(lfsu.size(), 0);

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>&
          rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);

        // loop over quadrature points
        for(typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin();
            it!=rule.end(); ++it) {
          // values of basefunctions
          std::vector<RangeType> phi(lfsu.size());
          lfsu.localFiniteElement().localBasis()
            .evaluateFunctionGlobal(it->position(),phi,eg.geometry());

          // curl of the basefunctions
          std::vector<JacobianType> J(lfsu.size());
          lfsu.localFiniteElement().localBasis()
            .evaluateJacobianGlobal(it->position(),J,eg.geometry());

          std::vector<CurlType> rotphi(lfsu.size());
          for(unsigned i = 0; i < lfsu.size(); ++i)
            jacobianToCurl(rotphi[i], J[i]);

          // calculate T
          Dune::FieldVector<RF,1> epsval;
          eps.evaluate(eg.entity(), it->position(), epsval);

          RF factor = it->weight()
            * eg.geometry().integrationElement(it->position()) * epsval;

          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              T[i][j] += factor * (phi[i] * phi[j]);

          // calculate S
          Dune::FieldVector<RF,1> muval;
          mu.evaluate(eg.entity(), it->position(), muval);

          factor = it->weight()
            * eg.geometry().integrationElement(it->position()) / muval;
            
          for(unsigned i = 0; i < lfsu.size(); ++i)
            for(unsigned j = 0; j < lfsu.size(); ++j)
              S[i][j] += factor * (rotphi[i] * rotphi[j]);

          // calculate f
          typename DtJ::Traits::RangeType dtJval;
          dtJ.evaluate(eg.entity(), it->position(), epsval);

          factor = it->weight()
            * eg.geometry().integrationElement(it->position());

          for(unsigned i = 0; i < lfsu.size(); ++i)
            f[i] += factor * (phi[i] * dtJval);
        }

        // get coefficients from Ecur and Eprev
        X xprev;
        lfsu.vread(*Eprev, xprev);
        
        X xcur;
        lfsu.vread(*Ecur, xcur);

        // v1 = x - 2*xcur + xprev, just what is needed to right-multiply to
        //     the matrix T in the first term
        X v1(lfsu.size());
        for(unsigned i = 0; i < lfsu.size(); ++i)
          v1[i] = x[i] - 2*xcur[i] + xprev[i];

        double Delta_t2 = Delta_t*Delta_t;
        // !!! modify S -> Delta_t^2*S
        for(unsigned i = 0; i < lfsu.size(); ++i)
          for(unsigned j = 0; j < lfsu.size(); ++j)
            S[i][j] *= Delta_t2;

        // !!! modify f -> Delta_t^2*f
        for(unsigned i = 0; i < lfsu.size(); ++i)
          f[i] *= Delta_t2;

        // calculate residual
        for(unsigned i = 0; i < lfsu.size(); ++i)
          for(unsigned j = 0; j < lfsu.size(); ++j)
            r[i] += T[i][j] * v1[j] + S[i][j] * xcur[j] + f[i];
	  }

      //! set Eprev
      /**
       * \param Eprev_   The coefficients of the electric field for time step n-1
       */
      void setEprev(const GCV &Eprev_)
      {
        Eprev = &Eprev_;
      }

      //! set Ecur
      /**
       * \param Ecur_    The coefficients of the electric field for time step n
       */
      void setEcur(const GCV &Ecur_)
      {
        Ecur = &Ecur_;
      }

      //! set Delta_t
      /**
       * \param Delta_t_ The time step
       */
      void setDelta_t(double Delta_t_)
      {
        Delta_t = Delta_t_;
      }

    private:
      const Eps &eps;
      const Mu &mu;
      const DtJ &dtJ;
      const GCV *Ecur;
      const GCV *Eprev;
      double Delta_t;
      const int qorder;
	};

    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ELECTRODYNAMIC_HH
