// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_IMEX_ONESTEP_PATTERNENGINE_HH
#define DUNE_PDELAB_IMEX_ONESTEP_PATTERNENGINE_HH

/**
 * \author Pavel Hron, Marian Piatkowski
 * \file
 * \brief Implementation of the local pattern assembler engine
 */

#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/imexonestep/enginebase.hh>

namespace Dune {
  namespace PDELab {

    /**
       \brief The local assembler engine for IMEXOneStep sub triangulations which
       creates the matrix pattern

       \tparam LA The local udg assembler

    */
    template<typename OSLA>
    class IMEXOneStepLocalPatternAssemblerEngine
      : public IMEXOneStepLocalAssemblerEngineBase<
      OSLA,
      typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
      typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine,
      typename OSLA::LocalAssemblerDT2::LocalPatternAssemblerEngine>
    {

      typedef IMEXOneStepLocalAssemblerEngineBase<
        OSLA,
        typename OSLA::LocalAssemblerDT0::LocalPatternAssemblerEngine,
        typename OSLA::LocalAssemblerDT1::LocalPatternAssemblerEngine,
        typename OSLA::LocalAssemblerDT2::LocalPatternAssemblerEngine
        > BaseT;

      using BaseT::la;
      using BaseT::lae0;
      using BaseT::lae1;
      using BaseT::lae2;
      using BaseT::implicit;
      using BaseT::setLocalAssemblerEngineDT0;
      using BaseT::setLocalAssemblerEngineDT1;
      using BaseT::setLocalAssemblerEngineDT2;

    public:
      //! The type of the wrapping local assembler
      typedef OSLA LocalAssembler;

      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;
      typedef typename OSLA::LocalAssemblerDT2 LocalAssemblerDT2;

      //! The type of the matrix pattern container
      typedef typename LocalAssembler::Traits::MatrixPattern Pattern;
      typedef Dune::PDELab::LocalSparsityPattern LocalPattern;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      IMEXOneStepLocalPatternAssemblerEngine(const LocalAssembler & la_)
        : BaseT(la_),
          invalid_pattern(static_cast<Pattern*>(0)), pattern(invalid_pattern)
      {}

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setPattern(Pattern & pattern_){

        // Set pointer to global pattern
        pattern = &pattern_;

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineDT1(la.la1.localPatternAssemblerEngine(pattern_));
        setLocalAssemblerEngineDT2(la.la2.localPatternAssemblerEngine(pattern_));
      }


      //! @name Notification functions
      //! @{
      void preAssembly(){
        implicit = la.osp_method->implicit();

        lae0->preAssembly();
        lae1->preAssembly();
        lae2->preAssembly();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
        lae2->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid solution pointer
      Pattern * const invalid_pattern;

      //! Pointer to the current matrix pattern container
      Pattern * pattern;

    }; // end class IMEXOneStepLocalPatternAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
