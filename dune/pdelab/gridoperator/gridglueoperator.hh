#ifndef DUNE_PDELAB_GRIDGLUEOPERATOR_HH
#define DUNE_PDELAB_GRIDGLUEOPERATOR_HH

#include "default/assembler.hh"
#include <dune/pdelab/gridfunctionspace/subspace.hh>

// #include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
// #include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
// #include <dune/pdelab/gridfunctionspace/interpolate.hh>
// #include <dune/common/tupleutility.hh>

namespace Dune{
  namespace PDELab{

    template<typename GFSU, typename GFSV,
             typename CU, typename CV>
    class GridGlueAssembler
    {
      GFSU gfsu_;
      GFSV gfsv_;
      const CU & cu_;
      const CV & cv_;

      typedef TypeTree::TreePath<0> Path0;
      typedef TypeTree::TreePath<1> Path1;
      typedef GridFunctionSubSpace<GFSU,Path0> GFSU0;
      typedef GridFunctionSubSpace<GFSV,Path0> GFSV0;
      typedef GridFunctionSubSpace<GFSU,Path1> GFSU1;
      typedef GridFunctionSubSpace<GFSV,Path1> GFSV1;

      GFSU0 gfsu0_;
      GFSV0 gfsv0_;
      GFSU1 gfsu1_;
      GFSV1 gfsv1_;

      typedef DefaultAssembler<GFSU0,GFSV0,CU,CV> Patch0Asm;
      typedef DefaultAssembler<GFSU1,GFSV1,CU,CV> Patch1Asm;

      Patch0Asm patch0asm_;
      Patch1Asm patch1asm_;

    public:
      GridGlueAssembler (const GFSU& gfsu, const GFSV& gfsv, const CU& cu, const CV& cv)
        : gfsu_(gfsu)
        , gfsv_(gfsv)
        , cu_(cu)
        , cv_(cv)
        , gfsu0_(gfsu)
        , gfsv0_(gfsv)
        , gfsu1_(gfsu)
        , gfsv1_(gfsv)
        , patch0asm_(gfsu0_,gfsv0_,cu_,cv_)
        , patch1asm_(gfsu1_,gfsv1_,cu_,cv_)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return gfsu_;
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return gfsv_;
      }

      template<class LocalAssemblerEngine>
      void assemble(LocalAssemblerEngine & assembler_engine) const
      {
        LocalFunctionSpace<GFSU> lfsu(gfsu_);
        LocalFunctionSpace<GFSV> lfsv(gfsv_);
        LocalFunctionSpace<GFSU> rlfsu(gfsu_);
        LocalFunctionSpace<GFSV> rlfsv(gfsv_);

        patch0asm_.assemble(assembler_engine);
        patch1asm_.assemble(assembler_engine);

        assert(& gfsu_.gridGlue() == & gfsv_.gridGlue());

        for (auto iit = gfsu_.gridGlue().template ibegin<0>();
             iit != gfsu_.gridGlue().template iend<0>();
             ++iit)
        {
          typedef typename GFSU::Traits::GridGlue::Grid0Patch::GridView::template Codim<0>::Entity Grid0Element;
          typedef typename GFSU::Traits::GridGlue::Grid1Patch::GridView::template Codim<0>::Entity Grid1Element;
          typedef GridGlueContext<Grid0Element,GFS_DOM0> Ctx0;
          typedef GridGlueContext<Grid1Element,GFS_DOM1> Ctx1;
          Ctx0 ctx0(*iit->inside());
          Ctx1 ctx1(*iit->outside());
          lfsu.bind(ctx0);
          rlfsu.bind(*iit);
          lfsv.bind(ctx0);
          rlfsv.bind(*iit);

          // coupling_0_to_1(lfsu.template child<0>(),lfsv.template child<0>(),
          //     rlfsu.template child<1>(),rlfsv.template child<1>());

          lfsu.bind(ctx1);
          // coupling_1_to_0(lfsu.template child<1>(),lfsv.template child<1>(),
          //     rlfsu.template child<0>(),rlfsv.template child<0>());
        }

        // communicate();

        // for (auto iit = gfsu.gridGlue().template ibegin<1>();
        //      iit != gfsu.gridGlue().template iend<1>();
        //      ++iit)
        //     coupling_0_from_1(lfsu.template child<0>(),lfsv.template child<0>(),
        //         rlfsu.template child<1>(),rlfsv.template child<1>());
        //     coupling_1_from_0(lfsu.template child<1>(),lfsv.template child<1>(),
        //         rlfsu.template child<0>(),rlfsv.template child<0>());
        // }
      }
    };

    /**
       \brief Standard grid operator implementation

       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam MB The matrix backend to be used for representation of the jacobian
       \tparam DF The domain field type of the operator
       \tparam RF The range field type of the operator
       \tparam JF The jacobian field type
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)
       \tparam nonoverlapping_mode Switch for nonoverlapping grids

    */
    template<typename GFSU, typename GFSV, typename LOP,
             typename MB, typename DF, typename RF, typename JF,
             typename CU=Dune::PDELab::EmptyTransformation,
             typename CV=Dune::PDELab::EmptyTransformation,
             bool nonoverlapping_mode = false>
    class GridGlueOperator
    {
    public:

      //! The global assembler type
      typedef GridGlueAssembler<GFSU,GFSV,CU,CV> Assembler;

      //! The type of the domain (solution).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSU,DF>::Type Domain;
      //! The type of the range (residual).
      typedef typename Dune::PDELab::BackendVectorSelector<GFSV,RF>::Type Range;
      //! The type of the jacobian.
      typedef typename Dune::PDELab::BackendMatrixSelector<MB,Domain,Range,JF>::Type Jacobian;

      //! The sparsity pattern container for the jacobian matrix
      typedef typename Jacobian::Pattern Pattern;

      //! The local assembler type
      typedef DefaultLocalAssembler<GridGlueOperator,LOP,nonoverlapping_mode>
      LocalAssembler;

      typedef typename SelectType<
        nonoverlapping_mode,
        NonOverlappingBorderDOFExchanger<GridGlueOperator>,
        OverlappingBorderDOFExchanger<GridGlueOperator>
        >::Type BorderDOFExchanger;

      //! The grid operator traits
      typedef Dune::PDELab::GridOperatorTraits
      <GFSU,GFSV,MB,DF,RF,JF,CU,CV,Assembler,LocalAssembler> Traits;

      template <typename MFT>
      struct MatrixContainer{
        typedef typename Traits::Jacobian Type;
      };

      //! Constructor for non trivial constraints
      GridGlueOperator(const GFSU & gfsu_, const CU & cu_, const GFSV & gfsv_, const CV & cv_, LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_,cu_,cv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_, cu_, cv_,dof_exchanger)
        , backend(mb_)
      {}

      //! Constructor for empty constraints
      GridGlueOperator(const GFSU & gfsu_, const GFSV & gfsv_, LOP & lop_, const MB& mb_ = MB())
        : global_assembler(gfsu_,gfsv_)
        , dof_exchanger(make_shared<BorderDOFExchanger>(*this))
        , local_assembler(lop_,dof_exchanger)
        , backend(mb_)
      {}

      //! Get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return global_assembler.trialGridFunctionSpace();
      }

      //! Get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return global_assembler.testGridFunctionSpace();
      }

      //! Get dimension of space u
      typename GFSU::Traits::SizeType globalSizeU () const
      {
        return trialGridFunctionSpace().globalSize();
      }

      //! Get dimension of space v
      typename GFSV::Traits::SizeType globalSizeV () const
      {
        return testGridFunctionSpace().globalSize();
      }

      Assembler & assembler() { return global_assembler; }

      const Assembler & assembler() const { return global_assembler; }

      LocalAssembler & localAssembler() const { return local_assembler; }


      //! Visitor which is called in the method setupGridOperators for
      //! each tuple element.
      template <typename GridOperatorTuple>
      struct SetupGridOperator {
        SetupGridOperator()
          : index(0), size(Dune::tuple_size<GridOperatorTuple>::value) {}

        template <typename T>
        void visit(T& elem) {
          elem.localAssembler().doPreProcessing = index == 0;
          elem.localAssembler().doPostProcessing = index == size-1;
          ++index;
        }

        int index;
        const int size;
      };

      //! Method to set up a number of grid operators which are used
      //! in a joint assembling. It is assumed that all operators are
      //! specializations of the same template type
      template<typename GridOperatorTuple>
      static void setupGridOperators(GridOperatorTuple tuple)
      {
        Dune::ForEachValue<GridOperatorTuple> forEach(tuple);
        SetupGridOperator<GridOperatorTuple> setup_visitor;
        forEach.apply(setup_visitor);
      }

      //! Interpolate the constrained dofs from given function
      template<typename F, typename X>
      void interpolate (const X& xold, F& f, X& x) const
      {
        DUNE_THROW(NotImplemented, "interpolate doesn't work for GridGlueOperator");
        // // Interpolate f into grid function space and set corresponding coefficients
        // Dune::PDELab::interpolate(f,global_assembler.trialGridFunctionSpace(),x);

        // // Copy non-constrained dofs from old time step
        // Dune::PDELab::copy_nonconstrained_dofs(local_assembler.trialConstraints(),xold,x);
      }

      //! Fill pattern of jacobian matrix
      void fill_pattern(Pattern & p) const {
        typedef typename LocalAssembler::LocalPatternAssemblerEngine PatternEngine;
        PatternEngine & pattern_engine = local_assembler.localPatternAssemblerEngine(p);
        global_assembler.assemble(pattern_engine);
      }

      //! Assemble residual
      void residual(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalResidualAssemblerEngine ResidualEngine;
        ResidualEngine & residual_engine = local_assembler.localResidualAssemblerEngine(r,x);
        global_assembler.assemble(residual_engine);
      }

      //! Assembler jacobian
      void jacobian(const Domain & x, Jacobian & a) const {
        typedef typename LocalAssembler::LocalJacobianAssemblerEngine JacobianEngine;
        JacobianEngine & jacobian_engine = local_assembler.localJacobianAssemblerEngine(a,x);
        global_assembler.assemble(jacobian_engine);
      }

      //! Apply jacobian matrix without explicitly assembling it
      void jacobian_apply(const Domain & x, Range & r) const {
        typedef typename LocalAssembler::LocalJacobianApplyAssemblerEngine JacobianApplyEngine;
        JacobianApplyEngine & jacobian_apply_engine = local_assembler.localJacobianApplyAssemblerEngine(r,x);
        global_assembler.assemble(jacobian_apply_engine);
      }

      void make_consistent(Jacobian& a) const {
        // we assume to work on consistent meshes
      }

    private:
      Assembler global_assembler;
      shared_ptr<BorderDOFExchanger> dof_exchanger;

      mutable LocalAssembler local_assembler;
      MB backend;

    };

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEOPERATOR_HH
