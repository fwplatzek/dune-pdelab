// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_IMEX_ONESTEP_JACOBIANENGINE_HH
#define DUNE_PDELAB_IMEX_ONESTEP_JACOBIANENGINE_HH

/**
 * \author Pavel Hron, Marian Piatkowski
 * \file
 * \brief Implementation of the IMEXOneStepLocalJacobianAssemblerEngine
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
    class IMEXOneStepLocalJacobianAssemblerEngine
      : public IMEXOneStepLocalAssemblerEngineBase<
      OSLA,
      typename OSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
      typename OSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine,
      typename OSLA::LocalAssemblerDT2::LocalJacobianAssemblerEngine>
    {

      typedef IMEXOneStepLocalAssemblerEngineBase<
        OSLA,
        typename OSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
        typename OSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine,
        typename OSLA::LocalAssemblerDT2::LocalJacobianAssemblerEngine
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

      typedef typename LocalAssemblerDT0::LocalJacobianAssemblerEngine JacobianEngineDT0;
      typedef typename LocalAssemblerDT1::LocalJacobianAssemblerEngine JacobianEngineDT1;
      typedef typename LocalAssemblerDT2::LocalJacobianAssemblerEngine JacobianEngineDT2;

      //! The type of the residual vector
      typedef typename OSLA::Traits::Jacobian Jacobian;

      //! The type of the solution vector
      typedef typename OSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename OSLA::Real Real;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      IMEXOneStepLocalJacobianAssemblerEngine(const LocalAssembler & local_assembler_)
        : BaseT(local_assembler_),
          invalid_jacobian(static_cast<Jacobian*>(0)),
          invalid_solution(static_cast<Solution*>(0)),
          jacobian(invalid_jacobian), solution(invalid_solution)
      {}


      //! Set current solution vector. Must be called before
      //! setResidual(). Should be called prior to assembling.
      void setSolution(const Solution & solution_){
        solution = &solution_;
      }

      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setJacobian(Jacobian & jacobian_){
        jacobian = &jacobian_;

        assert(solution != invalid_solution);

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localJacobianAssemblerEngine(*jacobian,*solution));
        setLocalAssemblerEngineDT1(la.la1.localJacobianAssemblerEngine(*jacobian,*solution));
        setLocalAssemblerEngineDT2(la.la2.localJacobianAssemblerEngine(*jacobian,*solution));
      }

      //! When multiple engines are combined in one assembling
      //! procedure, this method allows to reset the weights which may
      //! have been changed by the other engines.
      void setWeights(){
        la.la0.setWeight(b_rr * la.dt_factor0);
        la.la1.setWeight(eb_rr * la.dt_factor1);
        //la.la1.setWeight(0.0);
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

        // Here we only want to know whether this stage is implicit
        implicit = std::abs(b_rr) > 1e-6;

        // prepare local operators for stage
        la.la0.setTime(la.time + d_r * la.dt);
        la.la1.setTime(la.time + ed_r * la.dt);
        la.la2.setTime(la.time + d_r * la.dt);

        setWeights();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
        lae2->postAssembly(gfsu,gfsv);
      }
      //! @}

    private:

      //! Default value indicating an invalid residual pointer
      Jacobian * const invalid_jacobian;

      //! Default value indicating an invalid solution pointer
      Solution * const invalid_solution;

      //! Pointer to the current constant part residual vector in
      //! which to assemble
      Jacobian * jacobian;

      //! Pointer to the current residual vector in which to assemble
      const Solution * solution;

      //! Coefficients of time stepping scheme
      Real b_rr, d_r, eb_rr, ed_r;

    }; // end class IMEXOneStepLocalJacobianAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
