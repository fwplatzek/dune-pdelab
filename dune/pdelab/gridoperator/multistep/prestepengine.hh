// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_PRESTEPENGINE_HH
#define DUNE_PDELAB_MULTISTEP_PRESTEPENGINE_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepLocalPreStepAssemblerEngine
 */

#include <dune/pdelab/gridoperator/multistep/enginebase.hh>
#include <vector>
#include <cmath>

namespace Dune {
  namespace PDELab {

    /**
       \brief The local assembler engine for multi-step methods which
       assembles the constant part of the residual vector

       \tparam MSLA The local one step assembler

    */
    template<typename MSLA>
    class MultiStepLocalPreStepAssemblerEngine :
      public MultiStepLocalAssemblerEngineBase<
      MSLA,
      // TODO higher order problems with more than two grid-operators
      typename MSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
      typename MSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
      >
    {

      typedef MultiStepLocalAssemblerEngineBase<
        MSLA,
        typename MSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
        typename MSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine
        > BaseT;

      // TODO higher order problems with more than two grid-operators
      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;

    public:
      //! The type of the wrapping local assembler
      typedef MSLA LocalAssembler;

      //! Types of the subordinate assemblers and engines
      //! @{
      typedef typename MSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename MSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      typedef typename LocalAssemblerDT0::LocalResidualAssemblerEngine ResidualEngineDT0;
      typedef typename LocalAssemblerDT1::LocalResidualAssemblerEngine ResidualEngineDT1;
      //! @}

      //! The type of the residual vector
      typedef typename MSLA::Traits::Residual Residual;
      typedef typename Residual::ElementType ResidualElement;

      //! The type of the solution vector
      typedef typename MSLA::Traits::Solution Solution;
      typedef typename Solution::ElementType SolutionElement;

      //! The type for real numbers
      typedef typename MSLA::Real Real;

      //! The type of the solution container
      typedef std::vector<Solution*> Solutions;

      /**
         \brief Constructor

         \param [in] la_ The local assembler object which
         creates this engine
      */
      MultiStepLocalPreStepAssemblerEngine(LocalAssembler & la_)
        : BaseT(la_),
          invalid_residual(static_cast<Residual*>(0)),
          invalid_solutions(static_cast<Solutions*>(0)),
          const_residual_0(invalid_residual),
          const_residual_1(invalid_residual),
          solutions(invalid_solutions)
      {}

      // TODO higher order problems with more than two grid-operators
      //! Query methods for the global grid assembler
      //! @{
      bool requireSkeleton() const
      { return lae0->requireSkeleton() || lae1->requireSkeleton(); }
      //! @}

      //! Set current solution vector. Must be called before
      //! setConstResidual()! Should be called prior to assembling.
      void setSolutions(const Solutions & solutions_){
        solutions = &solutions_;
      }

      // TODO higher order problems with more than two grid-operators
      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setConstResiduals(Residual & const_residual_0_, Residual & const_residual_1_){
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*const_residual_0,*((*solutions)[0])));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*const_residual_1,*((*solutions)[0])));
      }

      // TODO higher order problems with more than two grid-operators
      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setConstResidual(Residual & const_residual_){
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solutions != invalid_solutions);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*const_residual_0,*((*solutions)[0])));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*const_residual_1,*((*solutions)[0])));
      }

      //! Methods for loading of the local function's
      //! coefficients. These methods are blocked. The loading of the
      //! coefficients is done in each assemble call.
      //!@{
      template<typename LFSU>
      void loadCoefficientsLFSUInside(const LFSU & lfsu_s){}
      template<typename LFSU>
      void loadCoefficientsLFSUOutside(const LFSU & lfsu_n){}
      template<typename LFSU>
      void loadCoefficientsLFSUCoupling(const LFSU & lfsu_c){}
      //! @}

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();

        *const_residual_0 = 0.0;
        *const_residual_1 = 0.0;

        // extract coefficients of the time step scheme
        // for the constant part of the residual
        // therefore the size is equal to la.steps
        alpha.resize(la.steps); // resize actually does nothing when the new size is equal to the old one
        beta.resize(la.steps);
        do0.resize(la.steps);
        do1.resize(la.steps);
        for(int i=0; i<la.steps; ++i) {
          alpha[i] = la.msp_method->alpha(i);
          beta[i] = la.msp_method->beta(i);
          do0[i] = (std::abs(beta[i]) > 1e-6);
          do1[i] = (std::abs(alpha[i]) > 1e-6);
        }

        // prepare local operators, set the stage to 1
        la.la0.preStage(la.time+la.dt,1);
        la.la1.preStage(la.time+la.dt,1);
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

      //! @ Assembling methods
      //! @{
      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu);
          lae1->loadCoefficientsLFSUInside(lfsu);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVVolume(eg,lfsu,lfsv);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVVolume(eg,lfsu,lfsv);
          }
        }
      }

      template<typename EG, typename LFSV>
      void assembleVVolume(const EG & eg, const LFSV & lfsv)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVVolume(eg,lfsv);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVVolume(eg,lfsv);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N>
      void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                              const LFSU_N & lfsu_n, const LFSV_N & lfsv_n)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUInside(lfsu_s);
          lae0->loadCoefficientsLFSUOutside(lfsu_n);
          lae1->loadCoefficientsLFSUOutside(lfsu_n);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVSkeleton(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n);
          }
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N>
      void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVSkeleton(ig,lfsv_s,lfsv_n);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVSkeleton(ig,lfsv_s,lfsv_n);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUInside(lfsu_s);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVBoundary(ig,lfsu_s,lfsv_s);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVBoundary(ig,lfsu_s,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSV_S>
      void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVBoundary(ig,lfsv_s);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVBoundary(ig,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S>
      void assembleUVProcessor(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUInside(lfsu_s);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVProcessor(ig,lfsu_s,lfsv_s);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVProcessor(ig,lfsu_s,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSV_S>
      void assembleVProcessor(const IG & ig, const LFSV_S & lfsv_s)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVProcessor(ig,lfsv_s);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVProcessor(ig,lfsv_s);
          }
        }
      }

      template<typename IG, typename LFSU_S, typename LFSV_S, typename LFSU_N, typename LFSV_N,
               typename LFSU_C, typename LFSV_C>
      void assembleUVEnrichedCoupling(const IG & ig,
                                             const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                             const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                             const LFSU_C & lfsu_c, const LFSV_C & lfsv_c)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu_s);
          lae1->loadCoefficientsLFSUInside(lfsu_s);

          lae0->loadCoefficientsLFSUOutside(lfsu_n);
          lae1->loadCoefficientsLFSUOutside(lfsu_n);

          lae0->loadCoefficientsLFSUCoupling(lfsu_c);
          lae1->loadCoefficientsLFSUCoupling(lfsu_c);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVEnrichedCoupling(ig,lfsu_s,lfsv_s,lfsu_n,lfsv_n,lfsu_c,lfsv_c);
          }
        }
      }

      template<typename IG, typename LFSV_S, typename LFSV_N, typename LFSV_C>
      void assembleVEnrichedCoupling(const IG & ig,
                                            const LFSV_S & lfsv_s,
                                            const LFSV_N & lfsv_n,
                                            const LFSV_C & lfsv_c)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVEnrichedCoupling(ig,lfsv_s,lfsv_n,lfsv_c);
          }
        }
      }

      template<typename EG, typename LFSU, typename LFSV>
      void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          lae0->setSolution(*((*solutions)[r]));
          lae1->setSolution(*((*solutions)[r]));

          lae0->loadCoefficientsLFSUInside(lfsu);
          lae1->loadCoefficientsLFSUInside(lfsu);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleUVVolumePostSkeleton(eg,lfsu,lfsv);
          }
        }
      }

      template<typename EG, typename LFSV>
      void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv)
      {
        for(int r=0; r<la.steps; ++r) {
          // reset the time in the local assembler
          // TODO pay attention for non-uniform time steps
          la.la0.setTime(la.time+(1-la.steps+r)*la.dt);
          la.la1.setTime(la.time+(1-la.steps+r)*la.dt);

          // check if coefficient in stationary part is zero ...
          // if so, skip the assembling
          if(do0[r]) {
            // TODO higher order problems needs power of dt_factor0
            la.la0.setWeight(beta[r]*la.dt_factor0);
            lae0->assembleVVolumePostSkeleton(eg,lfsv);
          }

          // check if coefficient in instationary part is zero ...
          // if so, skip the assembling
          if(do1[r]) {
            // TODO higher order problems needs power of dt_factor1
            la.la1.setWeight(alpha[r]*la.dt_factor1);
            lae1->assembleVVolumePostSkeleton(eg,lfsv);
          }
        }
      }
      //! @}

    private :
      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solutions * const invalid_solutions;

      //! Pointer to the current constant part residual vector in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      Residual * const_residual_0;
      Residual * const_residual_1;
      //! @}

      //! Pointer to the current residual vector in which to assemble
      const Solutions * solutions;

      //! Coefficients of time stepping scheme
      std::vector<Real> alpha;
      std::vector<Real> beta;
      std::vector<bool> do0;
      std::vector<bool> do1;
    }; // end class MultiStepLocalPreStepAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
