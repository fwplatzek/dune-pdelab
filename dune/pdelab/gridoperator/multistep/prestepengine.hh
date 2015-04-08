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
      typename MSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
      typename MSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine
      >
    {

      typedef MultiStepLocalAssemblerEngineBase<
        MSLA,
        typename MSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
        typename MSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine
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
        //============================================
        // TODO add code here
        //============================================
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv)
      {
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

      //============================================
      // TODO add assembling methods
      // Add all types of functions from the file
      // "gridoperator/onestep/prestageengine.hh", line 176
      //============================================
      //! @ Assembling methods
      //! @{

      //============================================
      // TODO add private members of this class
      //============================================
    private :
    }; // end class MultiStepLocalPreStepAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
