// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_GRIDOPERATOR_HH
#define DUNE_PDELAB_MULTISTEP_GRIDOPERATOR_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepGridOperator
 */

//============================================
// TODO
// Implement member functions for explicit multi-step methods.
//============================================

#include <dune/pdelab/instationary/multistep.hh>
#include <dune/pdelab/gridoperator/multistep/localassembler.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/constraints/common/constraints.hh>

namespace Dune {
  namespace PDELab {

    /** \brief Multi-step grid-operator to be used in MultiStepMethod
     * \tparam GO0 Stationary grid-operator.
     * \tparam GO1 Instationary grid-operator of temporal order one (typically mass matrix involved).
     * \tparam implicit Determines whether method is implicit or not.
     */

    // TODO higher temporal order problems with more than two grid-operators
    template<typename GO0, typename GO1, bool implicit = true>
    class MultiStepGridOperator
    {
    public :

      //! The sparsity pattern container for the jacobian matrix
      typedef typename GO0::Pattern Pattern;

      //! The global UDG assembler type
      typedef typename GO0::Traits::Assembler Assembler;

      //! The local assembler types of the subordinate grid operators
      //! @{
      // TODO higher temporal order problems with more than two grid-operators
      typedef typename GO0::Traits::LocalAssembler LocalAssemblerDT0;
      typedef typename GO1::Traits::LocalAssembler LocalAssemblerDT1;
      //! @}

      //! The local assembler type
      typedef MultiStepLocalAssembler<MultiStepGridOperator,LocalAssemblerDT0,LocalAssemblerDT1> LocalAssembler;

      //! The BorderDOFExchanger
      typedef typename GO0::BorderDOFExchanger BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <typename GO0::Traits::TrialGridFunctionSpace,
       typename GO0::Traits::TestGridFunctionSpace,
       typename GO0::Traits::MatrixBackend,
       typename GO0::Traits::DomainField,
       typename GO0::Traits::RangeField,
       typename GO0::Traits::JacobianField,
       typename GO0::Traits::TrialGridFunctionSpaceConstraints,
       typename GO0::Traits::TestGridFunctionSpaceConstraints,
       Assembler,
       LocalAssembler> Traits;

      //! The io types of the operator
      //! @{
      typedef typename Traits::Domain Domain;
      typedef typename Traits::Range Range;
      typedef typename Traits::Jacobian Jacobian;
      //! @}

      template <typename MFT>
      struct MatrixContainer{
        typedef Jacobian Type;
      };

      //! The type for real number e.g. time
      typedef typename LocalAssembler::Real Real;

      //! The type of the one step method parameters
      typedef typename LocalAssembler::MultiStepParameters MultiStepParameters;

      /** constructor for non trivial constraints
       * \param method_ Parameter object. This chooses the actual method used.
       * \param go0_ stationary gridoperator object
       * \param go1_ instationary gridoperator object (typically mass matrix involved)
       */
      MultiStepGridOperator(const MultiStepParameters& method_, GO0 & go0_, GO1 & go1_) :
        method(&method_) ,
        global_assembler(go0_.assembler()) ,
        go0(go0_), go1(go1_) ,
        la0(go0_.localAssembler()), la1(go1_.localAssembler()) ,
        const_residual( go0_.testGridFunctionSpace() ) ,
        local_assembler(la0,la1, const_residual)
      {
        GO0::setupGridOperators(Dune::tie(go0_,go1_));
        local_assembler.setMethod(method_);
        if(!implicit)
          local_assembler.setDTAssemblingMode(LocalAssembler::DoNotAssembleDT);
      }

      /** \brief divides instationary term, i.e. mass term, by the time step size
       */
      void divideMassTermByDeltaT()
      {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::DivideOperator1ByDT);
      }

      /** \brief multiply stationary term by the time step size
       * \note this is the default behaviour
       */
      void multiplySpatialTermByDeltaT()
      {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");
        local_assembler.setDTAssemblingMode(LocalAssembler::MultiplyOperator0ByDT);
      }

      /** \brief get the trial grid function space
       */
      const typename Traits::TrialGridFunctionSpace& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      /** \brief get the test grid function space
       */
      const typename Traits::TestGridFunctionSpace& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

      /** \brief get dimension of trial function space
       */
      typename Traits::TrialGridFunctionSpace::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      /** \brief get dimension of test function space
       */
      typename Traits::TestGridFunctionSpace::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }

      /** \brief fill pattern of jacobian matrix
       */
      void fill_pattern(Pattern & p) const {
        if(implicit){
          typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
          PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
          global_assembler.assemble(pattern_engine);
        }
        //============================================
        // TODO
        // pattern in explicit case in else-condition
        //============================================
      }

      /** \brief Prepare for doing a step.
       */
      template<typename OldValues>
      void preStep(Real time_, Real dt_, const OldValues& oldValues)
      {
        // use existing methods and set the number of stages to 1
        local_assembler.preStep(time_, dt_, 1);
        // assemble constant part of residual
        typedef typename LocalAssembler::LocalPreStepAssemblerEngine PreStepEngine;
        PreStepEngine & prestep_engine = local_assembler.localPreStepAssemblerEngine(oldValues);
        global_assembler.assemble(prestep_engine);
      }

      /** \brief Assemble residual.
       */
      void residual(const Domain& x, Range& r) const {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");

        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      /** \brief Assemble jacobian.
       */
      void jacobian(const Domain& x, Jacobian& a) const {
        if(!implicit)
          DUNE_THROW(Dune::Exception,"This function should not be called in explicit mode");

        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
      }

      /** \brief Interpolate constrained values from a given function.
       */
      template<typename X, typename F>
      void interpolate(const X& xold, F& f, X& x) const
      {
        // set time in boundary value function
        f.setTime(local_assembler.timeAfterStep());

        go0.localAssembler().setTime(local_assembler.timeAfterStep());

        // interpolate
        go0.interpolate(xold,f,x);

        // copy non-constrained dofs from old time step
        Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(), xold, x);
      }

      /** \brief Set time stepping method.
       */
      void setMethod(const MultiStepParameters& method_)
      {
        method = &method_;
        local_assembler.setMethod(method_);
      }

      /** \brief To be called after step is completed.
       */
      void postStep()
      {
        la0.postStage();
        la1.postStage();
        la0.postStep();
        la1.postStep();
      }

      void update()
      {
        go0.update();
        go1.update();
        const_residual = Range(go0.testGridFunctionSpace());
      }

      const typename Traits::MatrixBackend& matrixBackend() const
      {
        return go0.matrixBackend();
      }

    private :
      Assembler & global_assembler;
      GO0 & go0;
      GO1 & go1;
      const MultiStepParameters *method;
      LocalAssemblerDT0 & la0;
      LocalAssemblerDT1 & la1;
      Range const_residual;
      mutable LocalAssembler local_assembler;

    }; // end class MultiStepGridOperator

  } // end namespace PDELab
} // end namespace Dune

#endif
