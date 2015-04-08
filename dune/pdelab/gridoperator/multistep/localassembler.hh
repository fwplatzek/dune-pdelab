// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_LOCAL_ASSEMBLER_HH
#define DUNE_PDELAB_MULTISTEP_LOCAL_ASSEMBLER_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepLocalAssembler
 */

#include <dune/typetree/typetree.hh>
#include <dune/pdelab/gridoperator/multistep/patternengine.hh>
#include <dune/pdelab/gridoperator/multistep/prestepengine.hh>
#include <dune/pdelab/gridoperator/multistep/residualengine.hh>
#include <dune/pdelab/gridoperator/multistep/jacobianengine.hh>
#include <dune/pdelab/instationary/multistepparameter.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>

namespace Dune {
  namespace PDELab {

    /** \brief The local assembler for multi-step methods.

        \tparam LA0 The local assembler for the temporal derivative term of order zero.
        \tparam La1 The local assembler for the temporal derivative term of order one.
    */
    template<typename GO, typename LA0, typename LA1>
    class MultiStepLocalAssembler :
      public Dune::PDELab::LocalAssemblerBase<
      typename GO::Traits::MatrixBackend,
      typename GO::Traits::TrialGridFunctionSpaceConstraints,
      typename GO::Traits::TestGridFunctionSpaceConstraints>
    {
    public :

      //! The types of the local assemblers of order one and zero
      typedef LA0 LocalAssemblerDT0;
      typedef LA1 LocalAssemblerDT1;

      typedef Dune::PDELab::LocalAssemblerTraits<GO> Traits;

      //! The base class
      typedef Dune::PDELab::LocalAssemblerBase
      <typename GO::Traits::MatrixBackend,
       typename GO::Traits::TrialGridFunctionSpaceConstraints,
       typename GO::Traits::TestGridFunctionSpaceConstraints> Base;

      //! The local assembler engines
      typedef MultiStepLocalPatternAssemblerEngine<MultiStepLocalAssembler> LocalPatternAssemblerEngine;
      typedef MultiStepLocalPreStepAssemblerEngine<MultiStepLocalAssembler> LocalPreStepAssemblerEngine;
      typedef MultiStepLocalResidualAssemblerEngine<MultiStepLocalAssembler> LocalResidualAssemblerEngine;
      typedef MultiStepLocalJacobianAssemblerEngine<MultiStepLocalAssembler> LocalJacobianAssemblerEngine;

      //======================
      // friend classes
      //======================
      friend class MultiStepLocalPatternAssemblerEngine<MultiStepLocalAssembler>;
      friend class MultiStepLocalPreStepAssemblerEngine<MultiStepLocalAssembler>;
      friend class MultiStepLocalResidualAssemblerEngine<MultiStepLocalAssembler>;
      friend class MultiStepLocalJacobianAssemblerEngine<MultiStepLocalAssembler>;

      void static_checks() {
        static_assert((is_same<typename LA0::Traits::Jacobian::Pattern,
                       typename LA1::Traits::Jacobian::Pattern>::value),
                      "Received two local assemblers which are non-compatible "
                      "due to different matrix pattern types");
        static_assert((is_same<typename LA0::Traits::Jacobian,
                       typename LA1::Traits::Jacobian>::value),
                      "Received two local assemblers which are non-compatible "
                      "due to different jacobian types");
        static_assert((is_same<typename LA0::Traits::Solution,
                       typename LA1::Traits::Solution>::value),
                      "Received two local assemblers which are non-compatible "
                      "due to different solution vector types");
        static_assert((is_same<typename LA0::Traits::Residual,
                            typename LA1::Traits::Residual>::value),
                           "Received two local assemblers which are non-compatible "
                           "due to different residual vector types");
      }

      //! The local operators type for real numbers e.g. time
      typedef typename Traits::RangeField Real;

      typedef Dune::PDELab::MultiStepParameterInterface<Real> MultiStepParameters;

      //! Constructor with empty constraints
      MultiStepLocalAssembler(LA0 & la0_, LA1 & la1_, typename Traits::Residual & const_residual_) :
        Base(la0_.trialConstraints(), la0_.testConstraints()) ,
        la0(la0_), la1(la1_) ,
        const_residual(const_residual_) ,
        time(0.0), dt_mode(MultiplyOperator0ByDT), steps(0)
        // TODO add engines here in constructor
      {
        static_checks();
      }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void preStep(Real time_, Real dt_, int stages_) {
        time = time_;
        dt = dt_;

        // This switch decides which term will be multiplied with dt
        if(dt_mode == DivideOperator1ByDT){
          dt_factor0 = 1.0;
          dt_factor1 = 1.0 / dt;
        }
        else if(dt_mode == MultiplyOperator0ByDT){
          dt_factor0 = dt;
          dt_factor1 = 1.0;
        }
        else{
          DUNE_THROW(Dune::Exception,"Unknown mode for assembling of time step size!");
        }

        la0.preStep(time_,dt_, stages_);
        la1.preStep(time_,dt_, stages_);
      }

      //! Set the multi-step method parameters
      void setMethod(const MultiStepParameters & method_) {
        msp_method = &method_;
        // NOTE Set number of steps only here, maybe we need to do that somewhere else, too.
        steps = msp_method->steps();
      }

      enum DTAssemblingMode { DivideOperator1ByDT, MultiplyOperator0ByDT };

      //! Determines whether the time step size is multiplied to the
      //! to the stationary term or divided in the instationary term.
      void setDTAssemblingMode(DTAssemblingMode dt_mode_) {
        dt_mode = dt_mode_;
      }

      //! Access time after step
      Real timeAfterStep() {
        return time+dt;
      }

      void setWeight(const Real weight) {
        la0.setWeight(weight);
        la1.setWeight(weight);
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPatternAssemblerEngine & localPatternAssemblerEngine
      (typename Traits::MatrixPattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine.
      //! This engine is completely configured and ready to use.
      LocalPreStepAssemblerEngine & localPreStepAssemblerEngine
      (const std::vector<typename Traits::Solution*>& x)
      {
        // TODO add here code
      }

      //! Returns a reference to the requested engine.
      //! This engine is completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (typename Traits::Residual & r, const typename Traits::Solution & x)
      {
        // TODO add here code
      }

      //! Returns a reference to the requested engine.
      //! This engine is completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
      {
        // TODO add here code
      }

      //============================================
      // TODO
      // Add engines for the explicit case.
      //============================================

    private :
      LA0 & la0;
      LA1 & la1;
      const MultiStepParameters* msp_method;
      typename Traits::Residual & const_residual;
      Real time;
      Real dt;
      Real dt_factor0, dt_factor1;
      DTAssemblingMode dt_mode;
      int steps;
      //============================================
      // TODO
      // Add engines here as private members.
      //============================================
    }; // end class MultiStepLocalAssembler

  } // end namespace PDELab
} // end namespace Dune

#endif
