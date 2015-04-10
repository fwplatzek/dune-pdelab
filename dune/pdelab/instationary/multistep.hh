// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_MULTISTEP_HH
#define DUNE_PDELAB_MULTISTEP_HH

/**
 * \author Marian Piatkowski
 * \file
 * \brief Implementation of the MultiStepMethod.
 */

//============================================
// TODO
// Implement ExplicitMultiStepMethod for explicit multi-step methods.
//============================================

#include <iostream>
#include <iomanip>

#include <dune/common/ios_state.hh>
#include <dune/pdelab/instationary/multistepparameter.hh>

namespace Dune {
  namespace PDELab {

    // status information of the solver (typically Newton)
    struct MultiStepMethodPartialResult
    {
      unsigned int timesteps;
      double assembler_time;
      double linear_solver_time;
      int linear_solver_iterations;
      int nonlinear_solver_iterations;

      MultiStepMethodPartialResult() :
        timesteps(0),
        assembler_time(0.0),
        linear_solver_time(0.0),
        linear_solver_iterations(0),
        nonlinear_solver_iterations(0)
      {}
    };

    struct MultiStepMethodResult
    {
      MultiStepMethodPartialResult total;
      MultiStepMethodPartialResult successful;
      MultiStepMethodResult() : total(), successful()
      {}
    };

    /** \brief Do one step of a time-stepping scheme.
     *
     * \tparam T          type to represent time values.
     * \tparam IGOS       multi-step gridoperator.
     * \tparam PDESOLVER  solver problem in each step (typically Newton).
     * \tparam TrlV       vector type to represent coefficients of solutions.
     * \tparam TstV       vector type to represent residuals.
     *
     */
    template<class T, class IGOS, class PDESOLVER, class TrlV, class TstV = TrlV>
    class MultiStepMethod
    {
      typedef typename PDESOLVER::Result PDESolverResult;

    public :
      typedef MultiStepMethodResult Result;

      /** construct a new multi-step scheme
       *
       * \param method_ Parameter object. This chooses the actual method used.
       * \param igos_ Multistep GridOperator.
       * \param pdesolver_ Solver object (typically Newton).
       *
       * The contructed method object stores references to the object it is
       * constructed with, so these objects should be valid for as long as the
       * constructed object is used (or until setMethod() is called, see
       * there).
       */
      MultiStepMethod(const MultiStepParameterInterface<T>& method_,
                      IGOS& igos_, PDESOLVER& pdesolver_)
        : method(&method_), igos(igos_), pdesolver(pdesolver_), verbosityLevel(1), timestep(1), res()
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
      }

      /** \brief change verbosity level; 0 means completely quiet
       */
      void setVerbosityLevel (int level)
      {
        if (igos.trialGridFunctionSpace().gridView().comm().rank()>0)
          verbosityLevel = 0;
        else
          verbosityLevel = level;
      }

      /** \brief change number of current step
       */
      void setStepNumber(int newstep) { timestep = newstep; }

      /** \brief access to the (non) linear solver
       */
      const PDESOLVER & getPDESolver() const
      {
        return pdesolver;
      }

      /** \brief Access to the (non) linear solver
       */
      PDESOLVER & getPDESolver()
      {
        return pdesolver;
      }

      const Result& result() const
      {
        return res;
      }

      /**
       * \brief redefine the method to be used; can be done before every step.
       *
       * \param method_ Parameter object.
       *
       * The MultiStepMethod object stores a reference to the method_ object.
       * The old method object is no longer referenced after this member
       * function returns.
       */
      void setMethod (const MultiStepParameterInterface<T>& method_)
      {
        method = &method_;
      }

      /**
       * \brief do on step.
       *
       * \param[in]  time start of time step.
       * \param[in]  dt time step size.
       * \param[in]  oldValues Vector of pointers to the old values. Must support
                     the expression *oldValues[i] such that *oldValues[0]
                     yields to a reference of the old value which lies in the past the most
                     and *oldValues[sizeof(oldValues)-1] yields to a reference
                     of the current solution.
       * \param[in,out] xnew new value at the end of the time step.
       */
      template<class OldValues>
      T apply(T time, T dt, const OldValues& oldValues, TrlV& xnew)
      {
        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if (verbosityLevel>=1) {
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << timestep
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(time, dt, oldValues);

        // solve stage
        try {
          pdesolver.apply(xnew);
        }
        catch(...)
          {
            // catch statistics if time step failed
            PDESolverResult pderes = pdesolver.result();
            res.total.assembler_time += pderes.assembler_time;
            res.total.linear_solver_time += pderes.linear_solver_time;
            res.total.linear_solver_iterations += pderes.linear_solver_iterations;
            res.total.nonlinear_solver_iterations += pderes.iterations;
            res.total.timesteps += 1;
            throw;
          }

        // step cleanup
        igos.postStep();

        // update statistics
        PDESolverResult pderes = pdesolver.result();
        res.total.assembler_time += pderes.assembler_time;
        res.total.linear_solver_time += pderes.linear_solver_time;
        res.total.linear_solver_iterations += pderes.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += pderes.iterations;
        res.total.timesteps += 1;

        res.successful.assembler_time += pderes.assembler_time;
        res.successful.linear_solver_time += pderes.linear_solver_time;
        res.successful.linear_solver_iterations += pderes.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += pderes.iterations;
        res.successful.timesteps += 1;

        if (verbosityLevel>=1) {
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        timestep++;
        return dt;

      } // end of method apply

      /**
       * \brief do one step.
       *
       * \param[in]  time start of time step.
       * \param[in]  dt suggested time step size.
       * \param[in]  oldValues Vector of pointers to the old values. Must support
                     the expression *oldValues[i] such that *oldValues[0]
                     yields to a reference of the old value which lies in the past the most
                     and *oldValues[sizeof(oldValues)-1] yields to a reference
                     of the current solution.
       * \param[in]  f function to interpolate boundary conditions from.
       * \param[in,out] xnew new value at the end of the time step.
       *
       * This is a version which interpolates constraints at the beginning of each time step.
       *
       */
      template<typename OldValues, typename F>
      T apply(T time, T dt, const OldValues& oldValues, F& f, TrlV& xnew)
      {
        // save formatting attributes
        ios_base_all_saver format_attribute_saver(std::cout);

        if (verbosityLevel>=1) {
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "TIME STEP [" << method->name() << "] "
                    << std::setw(6) << timestep
                    << " time (from): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time
                    << " dt: "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << dt
                    << " time (to): "
                    << std::setw(12) << std::setprecision(4) << std::scientific
                    << time+dt
                    << std::endl;
          std::cout.flags(oldflags);
        }

        // prepare assembler
        igos.preStep(time, dt, oldValues);

        // set boundary conditions and initial value
        //======================
        // NOTE Which one is more correct?
        //======================
        // igos.interpolate(*oldValues.back(),f,xnew);
        igos.interpolate(*(oldValues.back()),f,xnew);

        // solve stage
        try {
          pdesolver.apply(xnew);
        }
        catch(...)
          {
            // catch statistics if time step failed
            PDESolverResult pderes = pdesolver.result();
            res.total.assembler_time += pderes.assembler_time;
            res.total.linear_solver_time += pderes.linear_solver_time;
            res.total.linear_solver_iterations += pderes.linear_solver_iterations;
            res.total.nonlinear_solver_iterations += pderes.iterations;
            res.total.timesteps += 1;
            throw;
          }

        // step cleanup
        igos.postStep();

        // update statistics
        PDESolverResult pderes = pdesolver.result();
        res.total.assembler_time += pderes.assembler_time;
        res.total.linear_solver_time += pderes.linear_solver_time;
        res.total.linear_solver_iterations += pderes.linear_solver_iterations;
        res.total.nonlinear_solver_iterations += pderes.iterations;
        res.total.timesteps += 1;

        res.successful.assembler_time += pderes.assembler_time;
        res.successful.linear_solver_time += pderes.linear_solver_time;
        res.successful.linear_solver_iterations += pderes.linear_solver_iterations;
        res.successful.nonlinear_solver_iterations += pderes.iterations;
        res.successful.timesteps += 1;

        if (verbosityLevel>=1) {
          std::ios_base::fmtflags oldflags = std::cout.flags();
          std::cout << "::: timesteps      " << std::setw(6) << res.successful.timesteps
                    << " (" << res.total.timesteps << ")" << std::endl;
          std::cout << "::: nl iterations  " << std::setw(6) << res.successful.nonlinear_solver_iterations
                    << " (" << res.total.nonlinear_solver_iterations << ")" << std::endl;
          std::cout << "::: lin iterations " << std::setw(6) << res.successful.linear_solver_iterations
                    << " (" << res.total.linear_solver_iterations << ")" << std::endl;
          std::cout << "::: assemble time  " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.assembler_time << " (" << res.total.assembler_time << ")" << std::endl;
          std::cout << "::: lin solve time " << std::setw(12) << std::setprecision(4) << std::scientific
                    << res.successful.linear_solver_time << " (" << res.total.linear_solver_time << ")" << std::endl;
          std::cout.flags(oldflags);
        }

        timestep++;
        return dt;

      } // end of method apply which interpolates

    private :
      const MultiStepParameterInterface<T> *method;
      IGOS& igos;
      PDESOLVER& pdesolver;
      int verbosityLevel;
      int timestep;
      Result res;

    }; // end class MultiStepMethod

  } // end namespace PDELab
} // end namespace Dune
#endif
