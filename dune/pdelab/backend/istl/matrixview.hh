// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXVIEW_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/nullptr.hh>

namespace Dune {
  namespace PDELab {

    namespace istl {

      template<typename M_, typename RowCache, typename ColCache>
      class ConstMatrixView
      {

      public:

        typedef typename remove_const<M_>::type Container;

        dune_static_assert(
          (is_same<
             typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,
             typename Container::TestGridFunctionSpace
             >::value),
          "The RowCache passed to LocalView must belong to the underlying GFSV"
          );

        dune_static_assert(
          (is_same<
             typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,
             typename Container::TrialGridFunctionSpace
             >::value),
          "The ColCache passed to LocalView must belong to the underlying GFSU"
          );

      public:

        typedef typename Container::field_type E;
        typedef typename Container::size_type size_type;

        typedef E ElementType;

        typedef RowCache RowIndexCache;
        typedef ColCache ColIndexCache;

        typedef typename RowCache::LocalFunctionSpace LFSV;
        typedef typename ColCache::LocalFunctionSpace LFSU;

        typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
        typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

        typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
        typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

        ConstMatrixView()
          : _container(nullptr)
          , _row_cache(nullptr)
          , _col_cache(nullptr)
        {}

        ConstMatrixView(M_& container)
          : _container(&container)
          , _row_cache(nullptr)
          , _col_cache(nullptr)
        {}

        const RowIndexCache& rowIndexCache() const
        {
          assert(_row_cache);
          return *_row_cache;
        }

        const ColIndexCache& colIndexCache() const
        {
          assert(_col_cache);
          return *_col_cache;
        }

        void attach(M_& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const RowCache& row_cache, const ColCache& col_cache)
        {
          _row_cache = &row_cache;
          _col_cache = &col_cache;
        }

        void unbind()
        {}

        size_type N() const
        {
          return rowIndexCache().size();
        }

        size_type M() const
        {
          return colIndexCache().size();
        }

        template<typename LC, typename Accessor>
        void read(LC& local_container, Accessor access) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              local_container.getEntry(i,j) = container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          read(local_container,container().accessor(typename Container::template State<Container::built>()));
        }


        const ElementType& operator()(size_type i, size_type j) const
        {
          return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        const ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j) const
        {
          return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        const ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j) const
        {
          return container()(i,j);
        }

        const ElementType& operator()(const RowContainerIndex& i, size_type j) const
        {
          return container()(i,colIndexCache().containerIndex(j));
        }

        const ElementType& operator()(size_type i, const ColContainerIndex& j) const
        {
          return container()(rowIndexCache().containerIndex(i),j);
        }

        const Container& container() const
        {
          return *_container;
        }


        const Container& global_container() const DUNE_DEPRECATED_MSG("global_container() is deprecated, use container() instead.")
        {
          return *_container;
        }

      protected:

        M_* _container;
        const RowCache* _row_cache;
        const ColCache* _col_cache;

      };


      template<typename M_, typename RowCache, typename ColCache>
      class MatrixView
        : public ConstMatrixView<M_,RowCache,ColCache>
      {

        typedef ConstMatrixView<M_,RowCache,ColCache> BaseT;

      public:

        typedef M_ Container;
        typedef typename Container::ElementType ElementType;
        typedef typename Container::size_type size_type;

        typedef RowCache RowIndexCache;
        typedef ColCache ColIndexCache;

        typedef typename RowCache::LocalFunctionSpace LFSV;
        typedef typename ColCache::LocalFunctionSpace LFSU;

        typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
        typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

        typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
        typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

        using BaseT::rowIndexCache;
        using BaseT::colIndexCache;
        using BaseT::N;
        using BaseT::M;

        // Explicitly pull in operator() from the base class to work around a problem
        // with clang not finding the const overloads of the operator from the base class.
        using BaseT::operator();

        MatrixView()
        {}

        MatrixView(Container& container)
          : BaseT(container)
        {}

        void commit()
        {}


        template<typename LC, typename Accessor>
        void read(LC& local_container, Accessor access) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              local_container.getEntry(i,j) = access(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        template<typename LC, typename Accessor>
        void read(LC& local_container, Accessor access)
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              local_container.getEntry(i,j) = accessor(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        template<typename LC>
        void read(LC& local_container)
        {
          if (container().state() == Container::built)
            read(local_container,container().accessor(typename Container::template State<Container::built>()));
          else
            read(local_container,container().accessor(typename Container::template State<Container::building>()));
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          if (container().state() == Container::built)
            read(local_container,container().accessor(typename Container::template State<Container::built>()));
          else
            DUNE_THROW(Exception,"not supported");
        }


        template<typename LC, typename Accessor>
        void write(const LC& local_container, Accessor access)
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              access(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) = local_container.getEntry(i,j);
        }

        template<typename LC>
        void write(const LC& local_container)
        {
          if (container().state() == Container::built)
            write(local_container,container().accessor(typename Container::template State<Container::built>()));
          else
            write(local_container,container().accessor(typename Container::template State<Container::building>()));
        }


        template<typename LC, typename Accessor>
        void add(const LC& local_container, Accessor access)
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += local_container.getEntry(i,j);
        }

        template<typename LC>
        void add(const LC& local_container)
        {
          if (container().state() == Container::built)
            add(local_container,container().accessor(typename Container::template State<Container::built>()));
          else
            add(local_container,container().accessor(typename Container::template State<Container::building>()));
        }


        ElementType& operator()(size_type i, size_type j)
        {
          return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j)
        {
          return container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j));
        }

        ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j)
        {
          return container()(i,j);
        }

        ElementType& operator()(const RowContainerIndex& i, size_type j)
        {
          return container()(i,colIndexCache().containerIndex(j));
        }

        ElementType& operator()(size_type i, const ColContainerIndex& j)
        {
          return container()(rowIndexCache().containerIndex(i),j);
        }

        void add(size_type i, size_type j, const ElementType& v)
        {
          container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += v;
        }

        void add(const RowDOFIndex& i, const ColDOFIndex& j, const ElementType& v)
        {
          container()(rowIndexCache().containerIndex(i),colIndexCache().containerIndex(j)) += v;
        }

        void add(const RowContainerIndex& i, const ColContainerIndex& j, const ElementType& v)
        {
          container()(i,j) += v;
        }

        void add(const RowContainerIndex& i, size_type j, const ElementType& v)
        {
          container()(i,colIndexCache().containerIndex(j)) += v;
        }

        void add(size_type i, const ColContainerIndex& j, const ElementType& v)
        {
          container()(rowIndexCache().containerIndex(i),j) += v;
        }

        Container& container()
        {
          return *(this->_container);
        }

        Container& global_container() DUNE_DEPRECATED_MSG("global_container() is deprecated, use container() instead.")
        {
          return *(this->_container);
        }

      };

    } // namespace istl

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_MATRIXVIEW_HH
