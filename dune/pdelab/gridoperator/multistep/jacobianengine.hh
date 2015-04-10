// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_JACOBIANENGINE_HH
#define DUNE_PDELAB_MULTISTEP_JACOBIANENGINE_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepLocalJacobianAssemblerEngine
 */

#include <dune/pdelab/gridoperator/multistep/enginebase.hh>
#include <cmath>

namespace Dune {
  namespace PDELab {

    /**
       \brief The local assembler engine for one step methods which
       assembles the residual vector

       \tparam LA The local one step assembler

    */
    template<typename MSLA>
    class MultiStepLocalJacobianAssemblerEngine :
      public MultiStepLocalAssemblerEngineBase<
      MSLA,
      // TODO higher order problems with more than two grid-operators
      typename MSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
      typename MSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine
      >
    {

      typedef MultiStepLocalAssemblerEngineBase<
        MSLA,
        typename MSLA::LocalAssemblerDT0::LocalJacobianAssemblerEngine,
        typename MSLA::LocalAssemblerDT1::LocalJacobianAssemblerEngine
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

      typedef typename LocalAssemblerDT0::LocalJacobianAssemblerEngine JacobianEngineDT0;
      typedef typename LocalAssemblerDT1::LocalJacobianAssemblerEngine JacobianEngineDT1;

      //! The type of the residual vector
      typedef typename MSLA::Traits::Jacobian Jacobian;

      //! The type of the solution vector
      typedef typename MSLA::Traits::Solution Solution;

      //! The type for real numbers
      typedef typename MSLA::Real Real;

      /**
         \brief Constructor

         \param [in] local_assembler_ The local assembler object which
         creates this engine
      */
      MultiStepLocalJacobianAssemblerEngine(const LocalAssembler & local_assembler_)
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

      // TODO higher order problems with more than two grid-operators
      //! Set current residual vector. Should be called prior to
      //! assembling.
      void setJacobian(Jacobian & jacobian_){
        jacobian = &jacobian_;

        assert(solution != invalid_solution);

        // Initialize the engines of the two wrapped local assemblers
        setLocalAssemblerEngineDT0(la.la0.localJacobianAssemblerEngine(*jacobian,*solution));
        setLocalAssemblerEngineDT1(la.la1.localJacobianAssemblerEngine(*jacobian,*solution));
      }

      //! When multiple engines are combined in one assembling
      //! procedure, this method allows to reset the weights which may
      //! have been changed by the other engines.
      void setWeights(){
        // TODO higher order problems needs power of dt_factor0
        la.la0.setWeight(beta_R * la.dt_factor0);
        // TODO think about normalization of the leading coeff alpha_R
        // TODO higher order problems needs power of dt_factor1
        la.la1.setWeight(alpha_R * la.dt_factor1);
      }

      //! Notifier functions, called immediately before and after assembling
      //! @{
      void preAssembly()
      {
        lae0->preAssembly();
        lae1->preAssembly();

        // extract coefficients of the time step scheme
        alpha_R = la.msp_method->alpha(la.steps);
        beta_R = la.msp_method->beta(la.steps);
        implicit = std::abs(beta_R) > 1e-6;

        // set time in local operators
        la.la0.setTime(la.time + la.dt);
        la.la0.setTime(la.time + la.dt);

        setWeights();
      }

      template<typename GFSU, typename GFSV>
      void postAssembly(const GFSU& gfsu, const GFSV& gfsv){
        lae0->postAssembly(gfsu,gfsv);
        lae1->postAssembly(gfsu,gfsv);
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
      Real alpha_R, beta_R;
    }; // end class MultiStepLocalJacobianAssemblerEngine

  } // end namespace PDELab
} // end namespace Dune
#endif
