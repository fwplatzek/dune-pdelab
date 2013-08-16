// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_HETEROGENEOUSMATRIXBACKEND_HH
#define DUNE_PDELAB_HETEROGENEOUSMATRIXBACKEND_HH

#include <dune/pdelab/backend/istl/heterogeneousvector.hh>
#include <dune/pdelab/backend/istl/heterogeneousmatrix.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>

namespace Dune {
  namespace PDELab {

    template<typename E, typename RV, typename CV, std::size_t j, std::size_t m, typename... Columns>
    struct build_matrix_row
    {

      typedef typename istl::build_matrix_type<E,RV,typename CV::template Block<j>::type>::type block_type;

      typedef typename build_matrix_row<E,RV,CV,j+1,m,Columns...,block_type>::type type;

    };

    template<typename E, typename RV, typename CV, std::size_t m, typename... Columns>
    struct build_matrix_row<E,RV,CV,m,m,Columns...>
    {
      typedef HeterogeneousVector<Columns...> type;
    };

    template<typename E, typename RV, typename CV, std::size_t i, std::size_t n, typename... Rows>
    struct build_matrix_rows
    {
      typedef typename build_matrix_row<E,typename RV::template Block<i>::type,CV,0,CV::block_count>::type row_type;

      typedef typename build_matrix_rows<E,RV,CV,i+1,n,Rows...,row_type>::type type;
    };

    template<typename E, typename RV, typename CV, std::size_t n, typename... Rows>
    struct build_matrix_rows<E,RV,CV,n,n,Rows...>
    {
      typedef HeterogeneousMatrix<Rows...> type;
    };


    template<typename E, typename RV, typename CV>
    struct build_heterogeneous_matrix
    {

      typedef typename build_matrix_rows<E,RV,CV,0,RV::block_count>::type type;

    };


    template<typename RowOrdering, typename ColOrdering, typename SubPattern_>
    class DensePattern
    {

    public:

      typedef SubPattern_ SubPattern;
      typedef std::size_t size_type;

      template<typename RI, typename CI>
      void add_link(const RI& ri, const CI& ci)
      {
        recursive_add_entry(ri.view(),ci.view());
      }

      template<typename RI, typename CI>
      void recursive_add_entry(const RI& ri, const CI& ci)
      {
        block(ri.back(),ci.back()).recursive_add_entry(ri.back_popped(),ci.back_popped());
      }

      DensePattern(const RowOrdering& row_ordering, const ColOrdering& col_ordering)
        : _row_ordering(row_ordering)
        , _col_ordering(col_ordering)
        , _col_count(col_ordering.blockCount())
      {
        for(std::size_t i = 0; i < row_ordering.blockCount(); ++i)
          for(std::size_t j = 0; j < _col_count; ++j)
            {
              _blocks.emplace_back(row_ordering.childOrdering(i),col_ordering.childOrdering(j));
              _blocks.back().resize(row_ordering.childOrdering(i).blockCount());
            }
      }

      SubPattern& block(size_type i, size_type j)
      {
        return _blocks[i * _col_count + j];
      }

      const SubPattern& block(size_type i, size_type j) const
      {
        return _blocks[i * _col_count + j];
      }

    private:

      const RowOrdering& _row_ordering;
      const ColOrdering& _col_ordering;
      const size_type _col_count;
      std::vector<SubPattern> _blocks;

    };


    template<typename GFSV, typename GFSU, typename C>
    class HeterogeneousMatrixContainer
    {

    public:

      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::field_type field_type;
      typedef typename C::size_type size_type;

      static const size_type row_count = C::row_block_count;
      static const size_type col_count = C::col_block_count;

      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;

      typedef typename GFSV::Ordering::Traits::ContainerIndex RowIndex;
      typedef typename GFSU::Ordering::Traits::ContainerIndex ColIndex;

      typedef DensePattern<
        OrderingBase<
          typename GFSV::Ordering::Traits::DOFIndex,
          typename GFSV::Ordering::Traits::GlobalDOFIndex,
          typename GFSV::Ordering::Traits::ContainerIndex
          >,
        OrderingBase<
          typename GFSU::Ordering::Traits::DOFIndex,
          typename GFSU::Ordering::Traits::GlobalDOFIndex,
          typename GFSU::Ordering::Traits::ContainerIndex
          >,
        typename istl::build_pattern_type<
          typename C::template Block<0,0>::type,
          GFSV,
          GFSU,
          typename GFSV::Ordering::ContainerAllocationTag
          >::type
        > Pattern;

      template<typename RowCache, typename ColCache>
      class LocalView
      {

        dune_static_assert((is_same<typename RowCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSV>::value),
                           "The RowCache passed to LocalView must belong to the underlying GFSV");

        dune_static_assert((is_same<typename ColCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFSU>::value),
                           "The ColCache passed to LocalView must belong to the underlying GFSU");

      public:

        typedef typename C::field_type E;
        typedef E ElementType;

        typedef RowCache RowIndexCache;
        typedef ColCache ColIndexCache;

        typedef typename RowCache::LocalFunctionSpace LFSV;
        typedef typename ColCache::LocalFunctionSpace LFSU;

        typedef typename LFSV::Traits::DOFIndex RowDOFIndex;
        typedef typename LFSV::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex RowContainerIndex;

        typedef typename LFSU::Traits::DOFIndex ColDOFIndex;
        typedef typename LFSU::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ColContainerIndex;

        LocalView()
          : _container(nullptr)
          , _row_cache(nullptr)
          , _col_cache(nullptr)
        {}

        LocalView(HeterogeneousMatrixContainer& container)
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

        void attach(HeterogeneousMatrixContainer& container)
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
        {
        }

        size_type N() const
        {
          return _row_cache->size();
        }

        size_type M() const
        {
          return _col_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              local_container.getEntry(i,j) = (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        template<typename LC>
        void write(const LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) = local_container.getEntry(i,j);
        }

        template<typename LC>
        void add(const LC& local_container) const
        {
          for (size_type i = 0; i < N(); ++i)
            for (size_type j = 0; j < M(); ++j)
              (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += local_container.getEntry(i,j);
        }

        void commit()
        {
        }


        ElementType& operator()(size_type i, size_type j)
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        const ElementType& operator()(size_type i, size_type j) const
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j)
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        const ElementType& operator()(const RowDOFIndex& i, const ColDOFIndex& j) const
        {
          return (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j));
        }

        ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j)
        {
          return (*_container)(i,j);
        }

        const ElementType& operator()(const RowContainerIndex& i, const ColContainerIndex& j) const
        {
          return (*_container)(i,j);
        }

        ElementType& operator()(const RowContainerIndex& i, size_type j)
        {
          return (*_container)(i,_col_cache->containerIndex(j));
        }

        const ElementType& operator()(const RowContainerIndex& i, size_type j) const
        {
          return (*_container)(i,_col_cache->containerIndex(j));
        }

        ElementType& operator()(size_type i, const ColContainerIndex& j)
        {
          return (*_container)(_row_cache->containerIndex(i),j);
        }

        const ElementType& operator()(size_type i, const ColContainerIndex& j) const
        {
          return (*_container)(_row_cache->containerIndex(i),j);
        }


        void add(size_type i, size_type j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
        }

        void add(const RowDOFIndex& i, const ColDOFIndex& j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),_col_cache->containerIndex(j)) += v;
        }

        void add(const RowContainerIndex& i, const ColContainerIndex& j, const ElementType& v)
        {
          (*_container)(i,j) += v;
        }

        void add(const RowContainerIndex& i, size_type j, const ElementType& v)
        {
          (*_container)(i,_col_cache->containerIndex(j)) += v;
        }

        void add(size_type i, const ColContainerIndex& j, const ElementType& v)
        {
          (*_container)(_row_cache->containerIndex(i),j) += v;
        }

        HeterogeneousMatrixContainer& container()
        {
          return *_container;
        }

        const HeterogeneousMatrixContainer& container() const
        {
          return *_container;
        }

      private:

        HeterogeneousMatrixContainer* _container;
        const RowCache* _row_cache;
        const ColCache* _col_cache;

      };

      template<typename GO>
      HeterogeneousMatrixContainer (const GO& go)
      {
        init_accessors();
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix_blocks(go.testGridFunctionSpace().ordering(),
                               go.trialGridFunctionSpace().ordering(),
                               pattern);
      }

      template<typename GO>
      HeterogeneousMatrixContainer (const GO& go, const E& e)
      {
        init_accessors();
        Pattern pattern(go.testGridFunctionSpace().ordering(),go.trialGridFunctionSpace().ordering());
        go.fill_pattern(pattern);
        allocate_matrix_blocks(go.testGridFunctionSpace().ordering(),
                               go.trialGridFunctionSpace().ordering(),
                               pattern);
        _container = e;
      }

      size_type N() const
      {
        return _container.N();
      }

      size_type M() const
      {
        return _container.M();
      }

      HeterogeneousMatrixContainer& operator= (const E& e)
      {
        _container = e;
        return *this;
      }

      HeterogeneousMatrixContainer& operator*= (const E& e)
      {
        _container *= e;
        return *this;
      }

      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        return block(ri.back(),ci.back()).access_element(_container,ri,ci);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        return block(ri.back(),ci.back()).access_element(_container,ri,ci);
      }

      const Container& base() const
      {
        return _container;
      }

      Container& base()
      {
        return _container;
      }

      void flush()
      {}

      void finalize()
      {}

      void clear_row(const RowIndex& ri, const E& diagonal_entry)
      {
        for(size_type j = 0; j < col_count; ++j)
          block(ri.back(),j).clear_row(_container,ri);
        (*this)(ri,ri) = diagonal_entry;
      }

    private:

      struct DynamicBlockAccessorBase
      {

      public:

        virtual E& access_element(Container& c, const RowIndex& ri, const ColIndex& ci) const = 0;
        virtual const E& access_element(const Container& c, const RowIndex& ri, const ColIndex& ci) const = 0;
        virtual void clear_row(Container& c, const RowIndex& ri) const = 0;
        virtual void allocate(Container& c,
                              const typename GFSV::Ordering& row_ordering,
                              const typename GFSU::Ordering& col_ordering,
                              const Pattern& pattern
                              ) const = 0;

        virtual ~DynamicBlockAccessorBase()
        {}

      };

      template<std::size_t i, std::size_t j>
      class DynamicBlockAccessor
        : public DynamicBlockAccessorBase
      {

      public:

        virtual E& access_element(Container& c, const RowIndex& ri, const ColIndex& ci) const
        {
          return istl::access_matrix_element(istl::container_tag(c.template block<i,j>()),
                                             c.template block<i,j>(),
                                             ri,ci,ri.size()-2,ci.size()-2);
        }

        virtual const E& access_element(const Container& c, const RowIndex& ri, const ColIndex& ci) const
        {
          return istl::access_matrix_element(istl::container_tag(c.template block<i,j>()),
                                             c.template block<i,j>(),
                                             ri,ci,ri.size()-2,ci.size()-2);
        }

        virtual void clear_row(Container& c, const RowIndex& ri) const
        {
          istl::clear_matrix_row(istl::container_tag(c.template block<i,j>()),c.template block<i,j>(),ri,ri.size()-2);
        }

        virtual void allocate(Container& c,
                              const typename GFSV::Ordering& row_ordering,
                              const typename GFSU::Ordering& col_ordering,
                              const Pattern& pattern
                              ) const
        {
          istl::allocate_matrix(row_ordering.childOrdering(i),
                                col_ordering.childOrdering(j),
                                pattern.block(i,j),
                                c.template block<i,j>());
        }

      };


      void allocate_matrix_blocks(const typename GFSV::Ordering& row_ordering,
                                  const typename GFSU::Ordering& col_ordering,
                                  const Pattern& pattern
                                  )
      {
        for (size_type i = 0; i < row_count; ++i)
          for (size_type j = 0; j < row_count; ++j)
            block(i,j).allocate(_container,row_ordering,col_ordering,pattern);
      }

      static const DynamicBlockAccessorBase& block(size_type i, size_type j)
      {
        return *_blocks[i * col_count + j];
      }

      static void init_accessors()
      {
        if (!_blocks[0])
          create_accessor_rows(TypeTree::index_range<row_count>());
      }

      template<std::size_t... i>
      static void create_accessor_rows(TypeTree::index_pack<i...>)
      {
        TypeTree::discard((create_accessor_row<i>(TypeTree::index_range<col_count>()),0)...);
      }

      template<std::size_t i, std::size_t... j>
      static void create_accessor_row(TypeTree::index_pack<j...>)
      {
        TypeTree::discard((create_accessor<i,j>(),0)...);
      }

      template<std::size_t i, std::size_t j>
      static void create_accessor()
      {
        static DynamicBlockAccessor<i,j> accessor;
        _blocks[i * col_count + j] = &accessor;
      }

      Container _container;

      static const DynamicBlockAccessorBase* _blocks[row_count * col_count];

    };

    template<typename GFSV, typename GFSU, typename C>
    const typename HeterogeneousMatrixContainer<GFSV,GFSU,C>::DynamicBlockAccessorBase*
    HeterogeneousMatrixContainer<GFSV,GFSU,C>::_blocks[HeterogeneousMatrixContainer<GFSV,GFSU,C>::row_count * HeterogeneousMatrixContainer<GFSV,GFSU,C>::col_count] = { nullptr };

    struct HeterogeneousMatrixBackend
    {

      typedef std::size_t size_type;

      template<typename VV, typename VU, typename E>
      struct MatrixHelper
      {
        typedef HeterogeneousMatrixContainer<typename VV::GridFunctionSpace,typename VU::GridFunctionSpace,typename build_heterogeneous_matrix<E,typename VV::Container,typename VU::Container>::type > type;
      };
    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_HETEROGENEOUSMATRIXBACKEND_HH
