// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONER_HH
#define DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONER_HH

#include <dune/common/shared_ptr.hh>
#include <dune/istl/preconditioners.hh>

namespace Dune {
  namespace PDELab {
    namespace istl {

      template<typename X, typename Y>
      class DiagonalBlockPreconditioner
        : public Preconditioner<X,Y>
      {

      public:

        static const SolverCategory::Category category = SolverCategory::sequential;

        typedef std::size_t size_type;

        typedef Preconditioner<
          typename X::block_type,
          typename Y::block_type
          > Block;

        typedef shared_ptr<Block> BlockStorage;

        DiagonalBlockPreconditioner(std::initializer_list<BlockStorage> blocks)
          : _block_count(blocks.size())
          , _blocks(blocks)
        {}

        virtual void pre(X& x, Y& y)
        {
          for (size_type i = 0; i < _block_count; ++i)
            _blocks[i]->pre(x[i],y[i]);
        }

        virtual void apply(X& v, const Y& d)
        {
          for (size_type i = 0; i < _block_count; ++i)
            _blocks[i]->apply(v[i],d[i]);
        }

        virtual void post(X& x)
        {
          for (size_type i = 0; i < _block_count; ++i)
            _blocks[i]->post(x[i]);
        }

        virtual ~DiagonalBlockPreconditioner()
        {}

      private:

        const size_type _block_count;
        std::vector<shared_ptr<Block> > _blocks;

      };


      template<typename X, typename Y>
      struct build_heterogeneous_preconditioner_blocks;

      template<typename... XB, typename... YB>
      struct build_heterogeneous_preconditioner_blocks<
        HeterogeneousVector<XB...>,
        HeterogeneousVector<YB...>
        >
      {

        typedef tuple<Preconditioner<XB,YB>...> type;
        typedef tuple<shared_ptr<Preconditioner<XB,YB> >...> storage_type;

      };


      template<typename X, typename Y>
      class HeterogeneousDiagonalBlockPreconditioner
        : public Preconditioner<X,Y>
      {

      public:

        static const SolverCategory::Category category = SolverCategory::sequential;

        typedef std::size_t size_type;

        typedef typename build_heterogeneous_preconditioner_blocks<X,Y>::storage_type BlockStorage;

        HeterogeneousDiagonalBlockPreconditioner(BlockStorage&& block_storage)
          : _block_storage(std::forward<BlockStorage>(block_storage))
        {}


        virtual void pre(X& x, Y& y)
        {
          pre_blocked(x,y,TypeTree::tuple_indices(_block_storage));
        }

        template<std::size_t... i>
        void pre_blocked(X& x, Y& y, TypeTree::index_pack<i...>)
        {
          TypeTree::discard((get<i>(_block_storage)->pre(x.template block<i>(),y.template block<i>()),0)...);
        }


        virtual void apply(X& v, const Y& d)
        {
          apply_blocked(v,d,TypeTree::tuple_indices(_block_storage));
        }

        template<std::size_t... i>
        void apply_blocked(X& v, const Y& d, TypeTree::index_pack<i...>)
        {
          TypeTree::discard((get<i>(_block_storage)->apply(v.template block<i>(),d.template block<i>()),0)...);
        }


        virtual void post(X& x)
        {
          post_blocked(x,TypeTree::tuple_indices(_block_storage));
        }

        template<std::size_t... i>
        void post_blocked(X& x, TypeTree::index_pack<i...>)
        {
          TypeTree::discard((get<i>(_block_storage)->post(x.template block<i>()),0)...);
        }


        virtual ~HeterogeneousDiagonalBlockPreconditioner()
        {}

      private:

        BlockStorage _block_storage;

      };


      template<typename X, typename Y>
      class NoOpPreconditioner
        : public Preconditioner<X,Y>
      {

      public:

        static const SolverCategory::Category category = SolverCategory::sequential;

        virtual ~NoOpPreconditioner()
        {}

        virtual void pre(X& x, Y& y)
        {}

        virtual void post(X& x)
        {}

        virtual void apply(X& x, const Y& y)
        {
          x = y;
        }

      };


      template<typename M, typename X, typename Y>
      class SuperLUPreconditioner
        : public Preconditioner<X,Y>
      {

      public:

        static const SolverCategory::Category category = SolverCategory::sequential;

        SuperLUPreconditioner(const M& m)
          : _solver(m)
        {}

        virtual ~SuperLUPreconditioner()
        {}

        virtual void pre(X& x, Y& y)
        {
          _y.resize(y.size(),false);
        }

        virtual void post(X& x)
        {}

        virtual void apply(X& x, const Y& y)
        {
          _y = y;
          Dune::InverseOperatorResult res;
          std::cout << "update in: " << x.two_norm() << " defect in: " << y.two_norm();
          _solver.apply(x,_y,res);
          std::cout << " update out: " << x.two_norm() << " defect out: " << _y.two_norm() << std::endl;
        }

      private:

        SuperLU<M> _solver;
        Y _y;

      };

    } // namespace istl
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_PRECONDITIONER_HH
