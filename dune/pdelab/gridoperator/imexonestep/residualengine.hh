// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_IMEX_ONESTEP_RESIDUALENGINE_HH
#define DUNE_PDELAB_IMEX_ONESTEP_RESIDUALENGINE_HH

/**
 * \author Pavel Hron, Marian Piatkowski
 * \file
 * \brief Implementation of the IMEXOneStepLocalResidualAssemblerEngine
 */

#include <dune/pdelab/gridoperator/imexonestep/enginebase.hh>
#include <cmath>

namespace Dune {
  namespace PDELab {

    /**
       \brief The local assembler engine for one step methods which
       assembles the residual vector

       \tparam LA The local one step assembler

    */
    template<typename OSLA>
    class IMEXOneStepLocalResidualAssemblerEngine
      : public IMEXOneStepLocalAssemblerEngineBase<
      OSLA,
      typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
      typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine,
      typename OSLA::LocalAssemblerDT2::LocalResidualAssemblerEngine>
    {

      typedef IMEXOneStepLocalAssemblerEngineBase<
        OSLA,
        typename OSLA::LocalAssemblerDT0::LocalResidualAssemblerEngine,
        typename OSLA::LocalAssemblerDT1::LocalResidualAssemblerEngine,
        typename OSLA::LocalAssemblerDT2::LocalResidualAssemblerEngine
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
      typedef OSLA IMEXOneStepLocalAssembler;

      //! Types of the subordinate assemblers and engines
      //! @{
      typedef typename OSLA::LocalAssemblerDT0 LocalAssemblerDT0;
      typedef typename OSLA::LocalAssemblerDT1 LocalAssemblerDT1;
      typedef typename OSLA::LocalAssemblerDT2 LocalAssemblerDT2;

      typedef typename LocalAssemblerDT0::LocalResidualAssemblerEngine ResidualEngineDT0;
      typedef typename LocalAssemblerDT1::LocalResidualAssemblerEngine ResidualEngineDT1;
      typedef typename LocalAssemblerDT2::LocalResidualAssemblerEngine ResidualEngineDT2;
      //! @}

      //! The type of the residual vector
      typedef typename OSLA::Traits::Residual Residual;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      typedef OSLA LocalAssembler;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      IMEXOneStepLocalResidualAssemblerEngine(const LocalAssembler & local_assembler_)
        : BaseT(local_assembler_),
          invalid_residual(static_cast<Residual*>(0)),
          invalid_solution(static_cast<Solution*>(0)),
          residual_0(invalid_residual),
          residual_1(invalid_residual),
          residual_2(invalid_residual),
          const_residual_0(invalid_residual),
          const_residual_1(invalid_residual),
          const_residual_2(invalid_residual),
          solution(invalid_solution)
      {}

      //! Set current solution vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setSolution(const Solution & solution_){
        solution = &solution_;
      }

      //! Set current const residual vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setConstResidual(const Residual &const_residual_){
        const_residual_0 = &const_residual_;
        const_residual_1 = &const_residual_;
        const_residual_2 = &const_residual_;
      }

      //! Set current const residual vector. Should be called prior to
      //! assembling.
      void setResidual(Residual & residual_){
        residual_0 = &residual_;
        residual_1 = &residual_;
        residual_2 = &residual_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*residual_0,*solution));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*residual_1,*solution));
        setLocalAssemblerEngineDT2(la.la2.localResidualAssemblerEngine(*residual_2,*solution));
      }

      //! Set current const residual vectors. Must be called before
      //! setResidual(). Should be called prior to assembling. Here,
      //! separate vectors are used for the operators corresponding to
      //! the time dervatives of order one and zero.
      void setConstResiduals(const Residual &const_residual_0_, const Residual &const_residual_1_, const Residual &const_residual_2_){
        const_residual_0 = &const_residual_0_;
        const_residual_1 = &const_residual_1_;
        const_residual_2 = &const_residual_2_;
      }

      //! Set current const residual vectors. Should be called prior
      //! to assembling. Here, separate vectors are used for the
      //! operators corresponding to the time dervatives of order one
      //! and zero.
      void setResiduals(Residual & residual_0_, Residual & residual_1_, Residual & residual_2_){
        residual_0 = &residual_0_;
        residual_1 = &residual_1_;
        residual_2 = &residual_2_;

        // Initialize the engines of the two wrapped local assemblers
        assert(solution != invalid_solution);
        setLocalAssemblerEngineDT0(la.la0.localResidualAssemblerEngine(*residual_0,*solution));
        setLocalAssemblerEngineDT1(la.la1.localResidualAssemblerEngine(*residual_1,*solution));
        setLocalAssemblerEngineDT2(la.la2.localResidualAssemblerEngine(*residual_2,*solution));
      }

      //! When multiple engines are combined in one assembling
      //! procedure, this method allows to reset the weights which may
      //! have been changed by the other engines.
      void setWeights(){
        la.la0.setWeight(b_rr * la.dt_factor0);
        la.la1.setWeight(eb_rr * la.dt_factor1);
        la.la2.setWeight(la.dt_factor2);
      }

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();
        lae2->preAssembly();

        // Extract the coefficients of the time step scheme
        b_rr = la.osp_method->b(la.stage,la.stage);
        d_r = la.osp_method->d(la.stage);
        eb_rr = la.osp_method->eb(la.stage,la.stage);
        ed_r = la.osp_method->ed(la.stage);
        implicit = std::abs(b_rr) > 1e-6;

        // prepare local operators for stage
        la.la0.setTime(la.time + d_r * la.dt);
        la.la1.setTime(la.time + ed_r * la.dt);
        la.la2.setTime(la.time + d_r * la.dt);

        setWeights();
      }

      //todo need revision
      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){

        // Update residual vectors with constant part
        assert(const_residual_0 != invalid_residual);
        assert(const_residual_1 != invalid_residual);
        assert(const_residual_2 != invalid_residual);
        *residual_0 += *const_residual_0;
        if(residual_0 != residual_1){
          assert(const_residual_0 != const_residual_1);
          *residual_1 += *const_residual_1;
          *residual_2 += *const_residual_2;
        }

        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
        lae2->postAssembly(gfsu,gfsv);
      }
      //! @}


    private:

      //! Default value indicating an invalid residual pointer
      Residual * const invalid_residual;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current constant part residual vector in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      Residual * residual_0;
      Residual * residual_1;
      Residual * residual_2;
      //! @}

      //! Pointer to the current constant part residual vectors in
      //! which to assemble the residual corresponding to the operator
      //! representing the time derivative of order zero and one.
      //! @{
      const Residual * const_residual_0;
      const Residual * const_residual_1;
      const Residual * const_residual_2;
      //! @}

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r, eb_rr, ed_r;

    }; // end class IMEXOneStepLocalResidualAssemblerEngine
  } // end namespace PDELab
} // end namespace Dune
#endif
