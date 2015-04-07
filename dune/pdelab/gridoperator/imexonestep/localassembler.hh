// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_IMEX_ONESTEP_LOCAL_ASSEMBLER_HH
#define DUNE_PDELAB_IMEX_ONESTEP_LOCAL_ASSEMBLER_HH

/**
 * \author Pavel Hron, Marian Piatkowski
 * \file
 * \brief Implementation of the IMEXOneStepLocalAssembler
 */

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/imexonestep/patternengine.hh>
#include <dune/pdelab/gridoperator/imexonestep/prestageengine.hh>
#include <dune/pdelab/gridoperator/imexonestep/residualengine.hh>
#include <dune/pdelab/gridoperator/imexonestep/jacobianengine.hh>
#include <dune/pdelab/instationary/imexonestepparameter.hh>

namespace Dune {
  namespace PDELab {

    /** \brief The local assembler for the IMEX one-step methods.
     *
     * \tparam LA0 The local assembler for the implicitly treated temporal derivative term of order zero.
     * \tparam LA1 The local assembler for the explicitly treated temporal derivative term of order zero.
     * \tparam La2 The local assembler of the temporal derivative term of order one.
     */
    template<typename GO, typename LA0, typename LA1, typename LA2>
    class IMEXOneStepLocalAssembler
      : public Dune::PDELab::LocalAssemblerBase<
      typename GO::Traits::MatrixBackend,
      typename GO::Traits::TrialGridFunctionSpaceConstraints,
      typename GO::Traits::TestGridFunctionSpaceConstraints>
    {
    public:

      //! The types of the local assemblers of order one and zero
      typedef LA0 LocalAssemblerDT0;
      typedef LA1 LocalAssemblerDT1;
      typedef LA2 LocalAssemblerDT2;

      typedef Dune::PDELab::LocalAssemblerTraits<GO> Traits;

      //! The base class
      typedef Dune::PDELab::LocalAssemblerBase
      <typename GO::Traits::MatrixBackend,
       typename GO::Traits::TrialGridFunctionSpaceConstraints,
       typename GO::Traits::TestGridFunctionSpaceConstraints> Base;

      //! The local assembler engines
      //! @{
      typedef IMEXOneStepLocalPatternAssemblerEngine<IMEXOneStepLocalAssembler> LocalPatternAssemblerEngine;
      typedef IMEXOneStepLocalPreStageAssemblerEngine<IMEXOneStepLocalAssembler> LocalPreStageAssemblerEngine;
      typedef IMEXOneStepLocalResidualAssemblerEngine<IMEXOneStepLocalAssembler> LocalResidualAssemblerEngine;
      typedef IMEXOneStepLocalJacobianAssemblerEngine<IMEXOneStepLocalAssembler> LocalJacobianAssemblerEngine;

      friend class IMEXOneStepLocalPatternAssemblerEngine<IMEXOneStepLocalAssembler>;
      friend class IMEXOneStepLocalPreStageAssemblerEngine<IMEXOneStepLocalAssembler>;
      friend class IMEXOneStepLocalResidualAssemblerEngine<IMEXOneStepLocalAssembler>;
      friend class IMEXOneStepLocalJacobianAssemblerEngine<IMEXOneStepLocalAssembler>;
      //! @}

      void static_checks(){
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

      //! The type of the one step parameter object
      typedef Dune::PDELab::IMEXTimeSteppingParameterInterface<Real> OneStepParameters;

      //! Constructor with empty constraints
      IMEXOneStepLocalAssembler (LA0 & la0_, LA1 & la1_, LA2 & la2_, typename Traits::Residual & const_residual_)
        : Base(la0_.trialConstraints(),la0_.testConstraints()),
          la0(la0_), la1(la1_),la2(la2_),
          const_residual(const_residual_),
          time(0.0), dt_mode(MultiplyOperator0ByDT), stage(0),
          pattern_engine(*this), prestage_engine(*this), residual_engine(*this), jacobian_engine(*this)
      { static_checks(); }

      //! Notifies the local assembler about the current time of
      //! assembling. Should be called before assembling if the local
      //! operator has time dependencies.
      void preStep(Real time_, Real dt_, int stages_){
        time = time_;
        dt = dt_;

        // This switch decides which term will be multiplied with dt
        if(dt_mode == DivideOperator1ByDT){
          dt_factor0 = 1.0;
          dt_factor1 = 1.0;
          dt_factor2 = 1.0 / dt;
        }
        else if(dt_mode == MultiplyOperator0ByDT){
          dt_factor0 = dt;
          dt_factor1 = dt;
          dt_factor2 = 1.0;
        }
        else{
          DUNE_THROW(Dune::Exception,"Unknown mode for assembling of time step size!");
        }

        la0.preStep(time_,dt_, stages_);
        la1.preStep(time_,dt_, stages_);
        la2.preStep(time_,dt_, stages_);
      }

      //! Set the one step method parameters
      void setMethod(const OneStepParameters & method_){
        osp_method = & method_;
      }

      //! Set the current stage of the one step scheme
      void setStage(int stage_){
        stage = stage_;
      }

      enum DTAssemblingMode { DivideOperator1ByDT, MultiplyOperator0ByDT};

      //! Determines whether the time step size is multiplied to the
      //! mass term (first order time derivative) or the elliptic term
      //! (zero-th order time derivative).
      void setDTAssemblingMode(DTAssemblingMode dt_mode_){
        dt_mode = dt_mode_;
      }

      //! Access time at given stage
      Real timeAtStage(int stage_){
        return time+osp_method->d(stage_)*dt;
      }

      //! Access time at given stage
      Real timeAtStage(){
        return time+osp_method->d(stage)*dt;
      }

      void setWeight(const Real weight){
        la0.setWeight(weight);
        la1.setWeight(weight);
        la2.setWeight(weight);
      }

      //! Access methods which provid "ready to use" engines
      //! @{

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPatternAssemblerEngine & localPatternAssemblerEngine
      (typename Traits::MatrixPattern & p)
      {
        pattern_engine.setPattern(p);
        return pattern_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalPreStageAssemblerEngine & localPreStageAssemblerEngine
      (const std::vector<typename Traits::Solution*> & x)
      {
        prestage_engine.setSolutions(x);
        prestage_engine.setConstResidual(const_residual);
        return prestage_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalResidualAssemblerEngine & localResidualAssemblerEngine
      (typename Traits::Residual & r, const typename Traits::Solution & x)
      {
        residual_engine.setSolution(x);
        residual_engine.setConstResidual(const_residual);
        residual_engine.setResidual(r);
        return residual_engine;
      }

      //! Returns a reference to the requested engine. This engine is
      //! completely configured and ready to use.
      LocalJacobianAssemblerEngine & localJacobianAssemblerEngine
      (typename Traits::Jacobian & a, const typename Traits::Solution & x)
      {
        jacobian_engine.setSolution(x);
        jacobian_engine.setJacobian(a);
        return jacobian_engine;
      }


      //! @}

    private:

      //! The local assemblers for the temporal derivative of order
      //! one and zero
      //! @{
      LA0 & la0;
      LA1 & la1;
      LA2 & la2;
      //! @}

      //! The one step parameter object containing the generalized
      //! butcher tableau parameters
      const OneStepParameters * osp_method;

      //! The constant part of the residual
      typename Traits::Residual & const_residual;

      //! The current time of assembling
      Real time;

      //! The time step size
      Real dt;

      /** The time step factors for assembling. Depending on the value
          of \a dt_mode, it will hold:

          dt_factor0 = dt and dt_factor1 = 1.0

          or

          dt_factor0 = 1.0 and dt_factor1 = 1.0 / dt

          or

          dt_factor0 = 1.0 and dt_factor1 = 1.0 .

      */
      Real dt_factor0, dt_factor1, dt_factor2;

      //! Determines, whether the time step size will be multiplied
      //! with the time derivative term of first of zero-th order.
      DTAssemblingMode dt_mode;

      //! The current stage of the one step scheme
      int stage;

      //! The engine member objects
      //! @{
      LocalPatternAssemblerEngine  pattern_engine;
      LocalPreStageAssemblerEngine prestage_engine;
      LocalResidualAssemblerEngine residual_engine;
      LocalJacobianAssemblerEngine jacobian_engine;
      //! @}
    }; // end class IMEXOneStepLocalAssembler

  } // end namespace PDELab
} // end namespace Dune
#endif
