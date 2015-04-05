// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_PATTERNENGINE_HH
#define DUNE_PDELAB_MULTISTEP_PATTERNENGINE_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepLocalPatternAssemblerEngine
 */

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/multistep/enginebase.hh>

namespace Dune {
  namespace PDELab {

    /**
       \brief The local assembler engine for multi-step sub triangulations which
       creates the matrix pattern.

       \tparam MSLA The local udg assembler.

    */
    template<typename MSLA>
    class MultiStepLocalPatternAssemblerEngine :
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

      typedef typename MSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename MSLA::LocalAssemblerDT1 LocalAssemblerDT1;

      //! The type of the matrix pattern container
      typedef typename LocalAssembler::Traits::MatrixPattern Pattern;
      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      MultiStepLocalPatternAssemblerEngine(const LocalAssembler & la_)
        : BaseT(la_),
          invalid_pattern(static_cast<Pattern*>(0)), pattern(invalid_pattern)
      {}

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){

        // Set pointer to global pattern
        pattern = &pattern_;

        // Initialize the engines of the two wrapped local assemblers
        // TODO higher order problems with more than two grid-operators
        setLocalAssemblerEngineDT0(la.la0.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineDT1(la.la1.localPatternAssemblerEngine(pattern_));
      }


      //! @name Notification functions
      //! @{
      void preAssembly(){
        implicit = la.osp_method->implicit();

        // TODO higher order problems with more than two grid-operators
        lae0->preAssembly();
        lae1->preAssembly();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        // TODO higher order problems with more than two grid-operators
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;
    }; // end class MultiStepLocalPatternAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
