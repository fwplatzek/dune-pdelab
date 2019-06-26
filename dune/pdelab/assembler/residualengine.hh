// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH
#define DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/assembler/context.hh>
#include <dune/pdelab/assembler/elementdata.hh>
#include <dune/pdelab/assembler/functionspacedata.hh>
#include <dune/pdelab/assembler/vectordata.hh>
#include <dune/pdelab/localoperator/guardedcalls.hh>
#include <dune/pdelab/gridfunctionspace/flavor.hh>
#include <dune/pdelab/assembler/enginebase.hh>

namespace Dune::PDELab::Experimental {

  template<
    typename TrialVector_,
    typename TestVector_,
    typename LOP,
    typename TrialConstraints_ = EmptyTransformation,
    typename TestConstraints_ = EmptyTransformation,
    typename EngineParameters = DefaultResidualEngineParameters<false,Galerkin::automatic>
    >
  class ResidualEngine
    : public InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>
    , public FunctionSpaceProvider<typename TrialVector_::GridFunctionSpace,
                                   typename TestVector_::GridFunctionSpace,
                                   TrialConstraints_,
                                   TestConstraints_,
                                   not LocalOperator::disableFunctionSpaceFlavors<LOP>()
                                   >
  {

    static constexpr bool enable_flavors = not LocalOperator::disableFunctionSpaceFlavors<LOP>();

    using IEB = InstationaryEngineBase<typename TrialVector_::value_type,EngineParameters::instationary>;

    using FSP = FunctionSpaceProvider<
      typename TrialVector_::GridFunctionSpace,
      typename TestVector_::GridFunctionSpace,
      TrialConstraints_,
      TestConstraints_,
      enable_flavors
      >;

    using Types = typename FSP::Types;

  public:

    using IEB::instationary;
    using IEB::stage;
    using IEB::setStage;
    using IEB::oneStepMethod;
    using IEB::setOneStepMethod;
    using IEB::timestepFactor;
    using IEB::timestepTimeFactor;

    using FSP::unconstrained;
    using FSP::trialConstraints;
    using FSP::testConstraints;

    using size_type        = std::size_t;

    using TrialVector      = TrialVector_;
    using TestVector       = TestVector_;

    using TestSpace        = typename Types::TestSpace;
    using TestConstraints  = typename Types::TestConstraints;

    using TrialSpace       = typename Types::TrialSpace;
    using TrialConstraints = typename Types::TrialConstraints;

    using EntitySet        = typename TestSpace::Traits::EntitySet;

    using TimeReal         = typename TrialVector::value_type;

    static constexpr
    std::bool_constant<EngineParameters::template galerkin<TrialSpace,TestSpace>>
    isGalerkin()
    {
      return {};
    }

  private:

    LOP* _lop;

    bool _stage_accept_mode = false;

    std::vector<std::shared_ptr<TestVector>> _residuals;
    std::vector<std::shared_ptr<TestVector>> _time_residuals;

    const TrialVector* _coefficient;
    TestVector* _residual      = nullptr;;
    TestVector* _time_residual = nullptr;

  public:

    template<typename ElementFlavor>
    struct Data
      : public Context::RootContext
    {

      using Flavor = ElementFlavor;
      using Engine = ResidualEngine;
      using EntitySet = typename Engine::EntitySet;

      static constexpr bool assembleVariablePart()
      {
        return EngineParameters::assembleVariablePart;
      }

      static constexpr bool assembleConstantPart()
      {
        return EngineParameters::assembleConstantPart;
      }

      static constexpr bool assembleOffDiagonalSkeletonPart()
      {
        return EngineParameters::assembleOffDiagonalSkeletonPart;
      }

      static constexpr bool assembleDiagonalSkeletonPart()
      {
        return EngineParameters::assembleDiagonalSkeletonPart;
      }

      static constexpr auto isGalerkin()
      {
        return ResidualEngine::isGalerkin();
      }

      static constexpr std::bool_constant<EngineParameters::fastDG> fastDG()
      {
        return {};
      }

      Data(Engine& engine)
        : _engine(engine)
      {}

      Engine& engine()
      {
        return _engine;
      }

      const Engine& engine() const
      {
        return _engine;
      }

    private:

      Engine& _engine;

    };


    static constexpr bool intersectionsTwoSided()
    {
      return LocalOperator::intersectionsTwoSided<LOP>();
    }

    static constexpr bool visitBoundaryIntersections()
    {
      return std::is_invocable_v<decltype(LocalOperator::boundaryIntegral()),LOP,PlaceHolder&>;
    }

    static constexpr bool visitSkeletonIntersections()
    {
      return std::is_invocable_v<decltype(LocalOperator::skeletonIntegral()),LOP,PlaceHolder&>;
    }

    static constexpr bool visitPeriodicIntersections()
    {
      return visitSkeletonIntersections();
    }

    static constexpr bool visitProcessorIntersections()
    {
      return visitSkeletonIntersections();
    }

    template<typename LOP_>
    ResidualEngine(
      const TrialVector& trial_vector,
      TestVector& test_vector,
      LOP_& lop,
      std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
      )
      : _lop(&lop)
      , _coefficient(&trial_vector)
      , _residual(&test_vector)
    {}

    ResidualEngine(
      const TrialVector& trial_vector,
      TestVector& test_vector,
      LOP& lop,
      const TrialConstraints& trial_constraints,
      const TestConstraints& test_constraints,
      EngineParameters = {}
      )
      : FSP(&trial_constraints,&test_constraints)
      , _lop(&lop)
      , _coefficient(&trial_vector)
      , _residual(&test_vector)
    {}

    template<typename LOP_>
    ResidualEngine(
      LOP_& lop,
      std::enable_if_t<unconstrained() and std::is_same_v<LOP_,LOP>,EngineParameters> = {}
      )
      : _lop(&lop)
      , _coefficient(nullptr)
      , _residual(nullptr)
    {}

    ResidualEngine(
      LOP& lop,
      const TrialConstraints& trial_constraints,
      const TestConstraints& test_constraints,
      EngineParameters = {}
      )
      : FSP(&trial_constraints,&test_constraints)
      , _lop(&lop)
      , _coefficient(nullptr)
      , _residual(nullptr)
    {}

    void setCoefficient(const TrialVector& coefficient)
    {
      _coefficient = &coefficient;
    }

    void setResidual(TestVector& residual)
    {
      _residual = &residual;
    }

    TestVector& residual()
    {
      if (_stage_accept_mode)
      {
        assert(_residuals.size() > stage());
        assert(_residuals[stage()]);
        return *_residuals[stage()];
      }
      else
      {
        assert(_residual);
        return *_residual;
      }
    }

    template<typename h = int>
    TestVector& timeResidual()
    {
      static_assert(Std::to_true_type_v<h> and instationary(),"Calling timeResidual() is only allowed in instationary mode");
      if (_stage_accept_mode)
      {
        assert(_time_residuals.size() > stage());
        assert(_time_residuals[stage()]);
        return *_time_residuals[stage()];
      }
      else
      {
        assert(_time_residual);
        return *_time_residual;
      }
    }

    const TrialVector& coefficient() const
    {
      return *_coefficient;
    }

    const TestSpace& testSpace() const
    {
      return _residual->gridFunctionSpace();
    }

    const TrialSpace& trialSpace() const
    {
      return _coefficient->gridFunctionSpace();
    }

    LOP& localOperator()
    {
      return *_lop;
    }

    const LOP& localOperator() const
    {
      return *_lop;
    }

    void setOneStepMethod(shared_ptr<const OneStep::Method<TimeReal>> one_step_method)
    {
      IEB::setOneStepMethod(one_step_method);
      _residuals.resize(oneStepMethod().stages());
      _time_residuals.resize(oneStepMethod().stages());
      for (auto& r : _residuals)
        r = std::make_shared<TestVector>(testSpace());
      for (auto& r : _time_residuals)
        r = std::make_shared<TestVector>(testSpace());
    }

    template<typename Assembler>
    bool acceptStage(int new_stage, Assembler& assembler, const TrialVector& solution)
    {
      if (new_stage == stage())
        return false;
      _stage_accept_mode = true;
      updateWeights();
      auto coefficient = _coefficient;
      _coefficient = &solution;
      auto residual = _residual;
      auto time_residual = _time_residual;
      assembler.assemble(*this);
      _residual = residual;
      _time_residual = time_residual;
      _stage_accept_mode = false;
      IEB::acceptStage(new_stage,assembler,solution);
      if (stage() == oneStepMethod().stages())
        _coefficient = coefficient;
      return true;
    }

    void updateWeights()
    {
      IEB::updateWeights();
      if (instationary() and _stage_accept_mode)
      {
        IEB::setWeight(IEB::timestepFactor());
        IEB::setTimeWeight(IEB::timestepTimeFactor());
      }
    }

    template<typename Assembler>
    auto insideElementContext(const Assembler& assembler)
    {
      return
        insideElement(
          extractElementContext(
            *_lop,
            elementTimeResidualData<TestVector,Flavor::Test>(
              std::bool_constant<instationary()>{},
              elementResidualData(
                vectorData<TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                  elementCoefficientData(
                    vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                      trialSpaceData(
                        testSpaceData(
                          elementGridData(
                            timeData(
                              Data<ElementFlavor::Inside<enable_flavors>>(*this)
                              )))))))))));
    }

    template<typename Assembler, typename InsideElementContext>
    auto elementContext(const Assembler& assembler, InsideElementContext inside_element_context)
    {
      return
        elementDomainData(
          outsideElement(
            extractElementContext(
              *_lop,
              elementTimeResidualData<TestVector,Flavor::Test>(
                std::bool_constant<instationary()>{},
                elementResidualData(
                  vectorData<TestVector,Flavor::Test,LocalViewDataMode::accumulate>(
                    elementCoefficientData(
                      vectorData<TrialVector,Flavor::Trial,LocalViewDataMode::read>(
                        trialSpaceData(
                          testSpaceData(
                            elementGridData(
                              timeData(
                                Data<ElementFlavor::Outside<enable_flavors>>(*this)
                                )))))))))),
            std::forward<InsideElementContext>(inside_element_context)
            ));
    }

    template<typename Assembler, typename ElementContext>
    auto intersectionContext(const Assembler& assembler, ElementContext element_context)
    {
      return
        extractContext(
          *_lop,
          intersectionDomainData(
            std::forward<ElementContext>(element_context)
            ));
    }

    template<typename Assembler>
    auto context(const Assembler& assembler)
    {
      return intersectionContext(assembler,elementContext(assembler,insideElementContext(assembler)));
    }

    template<typename Context>
    void start(Context& ctx)
    {
      invoke_if_possible(LocalOperator::start(),*_lop,ctx);
    }

    template<typename Context, typename Element, typename Index>
    bool skipElement(Context& ctx, const Element& element, Index index) const
    {
      return invoke_or(LocalOperator::skipElement(),false,*_lop,ctx,element,index);
    }

    template<typename Context, typename Element, typename Index>
    void startElement(Context& ctx, const Element& element, Index index) const
    {
      invoke_if_possible(LocalOperator::startElement(),*_lop,ctx,element,index);
    }

    template<typename Context>
    void volume(Context& ctx)
    {
      invoke_if_possible(LocalOperator::volumeIntegral(),*_lop,ctx.elementContext());
    }

    template<typename Context>
    void startIntersections(Context& ctx)
    {
      invoke_if_possible(LocalOperator::startIntersections(),*_lop,ctx.elementContext());
    }

    template<typename Context>
    void skeleton(Context& ctx)
    {
      invoke_if_possible(LocalOperator::skeletonIntegral(),*_lop,ctx.intersectionContext());
    }

    template<typename Context>
    void periodic(Context& ctx)
    {
      invoke_if_possible(LocalOperator::skeletonIntegral(),*_lop,ctx.intersectionContext());
    }

    template<typename Context>
    void boundary(Context& ctx)
    {
      invoke_if_possible(LocalOperator::boundaryIntegral(),*_lop,ctx.intersectionContext());
    }

    template<typename Context>
    void processor(Context& ctx)
    {
      invoke_if_possible(LocalOperator::boundaryIntegral(),*_lop,ctx.intersectionContext());
    }

    template<typename Context>
    void finishIntersections(Context& ctx)
    {
      invoke_if_possible(LocalOperator::finishIntersections(),*_lop,ctx.elementContext());
    }

    template<typename Context>
    void volumePostIntersections(Context& ctx)
    {
      invoke_if_possible(LocalOperator::volumeIntegralPostIntersections(),*_lop,ctx.elementContext());
    }

    template<typename Context, typename Element, typename Index>
    void finishElement(Context& ctx, const Element& element, Index index) const
    {
      invoke_if_possible(LocalOperator::finishElement(),*_lop,ctx,element,index);
    }

    template<typename Context>
    void finish(Context& ctx)
    {
      invoke_if_possible(LocalOperator::finish(),*_lop,ctx);
      if (instationary() and not _stage_accept_mode)
      {
        const auto& osm = oneStepMethod();
        for (int r = 0; r < stage(); ++r) {
          if (osm.timeDerivativeActive(stage(),r))
            _time_residual->axpy(osm.timeDerivativeWeight(stage(),r)*timestepTimeFactor(),*_time_residuals[r]);
          if (osm.active(stage(),r))
            _residual->axpy(osm.weight(stage(),r)*timestepFactor(),*_residuals[r]);
        }
      }
      constrain_residual(testConstraints(),residual());
      if constexpr(instationary())
        if (_residual != _time_residual)
          constrain_residual(testConstraints(),timeResidual());
    }

    template<typename Context>
    void result(Context& ctx)
    {}

  };

  template<typename Coefficient, typename Residual, typename LOP>
  ResidualEngine(
      const Coefficient&,
      Residual&,
      LOP&
    )
    -> ResidualEngine<
      Coefficient,
      Residual,
      LOP,
      EmptyTransformation,
      EmptyTransformation,
      DefaultResidualEngineParameters<false,Galerkin::automatic>
      >;


  template<typename Coefficient, typename Residual, typename LOP, typename EngineParameters>
  ResidualEngine(
      const Coefficient&,
      Residual&,
      LOP&,
      EngineParameters
    )
    -> ResidualEngine<
      Coefficient,
      Residual,
      LOP,
      EmptyTransformation,
      EmptyTransformation,
      EngineParameters
      >;

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_RESIDUALENGINE_HH