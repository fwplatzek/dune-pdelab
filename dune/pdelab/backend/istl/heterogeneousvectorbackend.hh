// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_HETEROGENEOUSVECTORBACKEND_HH
#define DUNE_PDELAB_HETEROGENEOUSVECTORBACKEND_HH

#include <dune/pdelab/backend/istl/heterogeneousvector.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>


namespace Dune {
  namespace PDELab {

    struct heterogeneous_vector_backend_tag {};

    struct HeterogeneousVectorBackend
    {

      typedef heterogeneous_vector_backend_tag tag;

      typedef std::size_t size_type;

      static const size_type blockSize = 1;

      struct Traits
      {
        static const size_type block_size = 1;
        static const bool blocked = true;
        static const size_type max_blocking_depth = 1;
      };

      bool blocked() const
      {
        return true;
      }
    };


    template<typename LocalContainer, typename Container, typename LFSCache, typename Functor>
    struct transfer_data
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::StaticTraversal
    {

      template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(LFS&& lfs, Child&& child, TreePath tp, ChildIndex child_index)
      {
        auto& child_container = _container.template block<ChildIndex::value>();
        for (std::size_t i = 0, end = child.size(); i < end; ++i)
          {
            auto li = child.localIndex(i);
            auto& ci = _lfs_cache.containerIndex(li);
            // Here the local container spans the complete root LocalFunctionSpace,
            // so we need to translate the local indices for the local container.
            _f(accessBaseContainer(_local_container)[li],istl::access_vector_element(istl::container_tag(child_container),child_container,ci,ci.size()-2));
          }
      }

      transfer_data(LocalContainer& local_container,
                    Container& container,
                    const LFSCache& lfs_cache,
                    Functor f = Functor())
        : _local_container(local_container)
        , _container(container)
        , _lfs_cache(lfs_cache)
        , _f(f)
      {}

      LocalContainer& _local_container;
      Container& _container;
      const LFSCache& _lfs_cache;
      Functor _f;

    };


    template<typename LocalContainer, typename Container, typename LFSCache, typename TargetLFS, typename Functor>
    struct transfer_data_filtered
      : public TypeTree::DirectChildrenVisitor
      , public TypeTree::StaticTraversal
    {

      template<typename LFS, typename Child, typename TreePath, typename ChildIndex>
      void afterChild(LFS&& lfs, Child&& child, TreePath tp, ChildIndex ci)
      {
        if (_selected_child != ChildIndex::value)
          return;

        auto& child_container = _container.template block<ChildIndex::value>();
        for (std::size_t i = 0, end = _target_lfs.size(); i < end; ++i)
          {
            auto li = _target_lfs.localIndex(i);
            auto& ci = _lfs_cache.container_index(li);
            // Here the local container only spans the subspace of the target  LocalFunctionSpace,
            // so we must not translate the local indices for the local container.
            _f(accessBaseContainer(_local_container)[i],access_istl_vector_element(child_container,ci,ci.size()-2));
          }
      }

      transfer_data_filtered(LocalContainer& local_container,
                             Container& container,
                             const LFSCache& lfs_cache,
                             const TargetLFS& target_lfs,
                             std::size_t selected_child,
                             Functor f = Functor())
        : _local_container(local_container)
        , _container(container)
        , _lfs_cache(lfs_cache)
        , _target_lfs(target_lfs)
        , _selected_child(selected_child)
        , _f(f)
      {}

      LocalContainer& _local_container;
      Container& _container;
      const LFSCache& _lfs_cache;
      const TargetLFS& _target_lfs;
      const std::size_t _selected_child;
      Functor _f;

    };


    struct read_entry
    {
      template<typename LE, typename GE>
      void operator()(LE& local_entry, const GE& global_entry) const
      {
        local_entry = global_entry;
      }
    };

    struct write_entry
    {
      template<typename LE, typename GE>
      void operator()(const LE& local_entry, GE& global_entry) const
      {
        global_entry = local_entry;
      }
    };

    struct add_entry
    {
      template<typename LE, typename GE>
      void operator()(const LE& local_entry, GE& global_entry) const
      {
        global_entry += local_entry;
      }
    };

    template<typename LocalContainer, typename Container, typename LFSCache>
    using read_data = transfer_data<LocalContainer,Container,LFSCache,read_entry>;

    template<typename LocalContainer, typename Container, typename LFSCache>
    using write_data = transfer_data<LocalContainer,Container,LFSCache,write_entry>;

    template<typename LocalContainer, typename Container, typename LFSCache>
    using add_data = transfer_data<LocalContainer,Container,LFSCache,add_entry>;

    template<typename LocalContainer, typename Container, typename LFSCache, typename ChildLFS>
    using read_data_filtered = transfer_data_filtered<LocalContainer,Container,LFSCache,ChildLFS,read_entry>;

    template<typename LocalContainer, typename Container, typename LFSCache, typename ChildLFS>
    using write_data_filtered = transfer_data_filtered<LocalContainer,Container,LFSCache,ChildLFS,write_entry>;

    template<typename LocalContainer, typename Container, typename LFSCache, typename ChildLFS>
    using add_data_filtered = transfer_data_filtered<LocalContainer,Container,LFSCache,ChildLFS,add_entry>;



    template<typename GFS, typename C>
    class DynamicBlockAccessorBase
    {

    public:

      typedef C Container;
      typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;
      typedef typename C::field_type E;

      virtual E& access_element(Container& c, const ContainerIndex& ci) const = 0;
      virtual const E& access_element(const Container& c, const ContainerIndex& ci) const = 0;
      virtual void allocate(Container& c, const typename GFS::Ordering& ordering) const = 0;

      virtual ~DynamicBlockAccessorBase()
      {}

    };

    template<typename GFS, typename C, std::size_t i>
    class DynamicBlockAccessor
      : public DynamicBlockAccessorBase<GFS,C>
    {

    public:

      typedef typename C::field_type E;
      typedef typename GFS::Ordering::Traits::ContainerIndex CI;

      virtual E& access_element(C& container, const CI& ci) const
      {
        return istl::access_vector_element(istl::container_tag(container.template block<i>()),container.template block<i>(),ci,ci.size()-2);
      }

      virtual const E& access_element(const C& container, const CI& ci) const
      {
        return istl::access_vector_element(istl::container_tag(container.template block<i>()),container.template block<i>(),ci,ci.size()-2);
      }

      virtual void allocate(C& container, const typename GFS::Ordering& ordering) const
      {
        istl::resize_vector(istl::container_tag(container.template block<i>()),container.template block<i>(),ordering.childOrdering(i).blockCount(),false);
        istl::dispatch_vector_allocation(ordering.childOrdering(i),
                                         container.template block<i>(),
                                         typename GFS::Ordering::ContainerAllocationTag()
                                         );
      }

    };


    template<typename GFS, typename C>
    class HeterogeneousVectorContainer
    {

    public:
      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef GFS GridFunctionSpace;
      typedef Container BaseT;
      typedef typename Container::field_type field_type;
      typedef typename Container::size_type size_type;

      typedef typename GFS::Ordering::Traits::ContainerIndex ContainerIndex;

      template<size_type i>
      using Block = typename Container::template Block<i>;

      template<typename LFSCache>
      struct LocalView
      {

        //dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
        //                   "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

        typedef E ElementType;
        //typedef typename LFSCache::LocalFunctionSpace LFS;
        typedef typename LFSCache::DOFIndex DOFIndex;
        typedef typename LFSCache::ContainerIndex ContainerIndex;

        LocalView()
          : _container(nullptr)
          , _lfs_cache(nullptr)
        {}

        LocalView(HeterogeneousVectorContainer& container)
          : _container(&container)
          , _lfs_cache(nullptr)
        {}

        void attach(HeterogeneousVectorContainer& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const LFSCache& lfs_cache)
        {
          _lfs_cache = &lfs_cache;
        }

        void unbind()
        {
        }

        size_type size() const
        {
          return _lfs_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          read_data<
            LC,
            const Container,
            LFSCache
            > visitor(local_container,*_container,*_lfs_cache);

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename LC>
        void write(const LC& local_container)
        {
          write_data<
            const LC,
            Container,
            LFSCache
            > visitor(local_container,*_container,*_lfs_cache);

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename LC>
        void add(const LC& local_container)
        {
          add_data<
            const LC,
            Container,
            LFSCache
            > visitor(local_container,*_container,*_lfs_cache);

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename ChildLFS, typename LC>
        void read(const ChildLFS& child_lfs, LC& local_container) const
        {
          if (child_lfs.size() == 0)
            return;

          read_data_filtered<
            LC,
            const Container,
            LFSCache,
            ChildLFS
            > visitor(local_container,
                      *_container,
                      *_lfs_cache,
                      child_lfs,
                      _lfs_cache->container_index(child_lfs.localIndex(0)).back()
                      );

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename ChildLFS, typename LC>
        void write(const ChildLFS& child_lfs, const LC& local_container)
        {
          if (child_lfs.size() == 0)
            return;

          write_data_filtered<
            const LC,
            Container,
            LFSCache,
            ChildLFS
            > visitor(local_container,
                      *_container,
                      *_lfs_cache,
                      child_lfs,
                      _lfs_cache->container_index(child_lfs.localIndex(0)).back()
                      );

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename ChildLFS, typename LC>
        void add(const ChildLFS& child_lfs, const LC& local_container)
        {
          if (child_lfs.size() == 0)
            return;

          add_data_filtered<
            const LC,
            Container,
            LFSCache,
            ChildLFS
            > visitor(local_container,
                      *_container,
                      *_lfs_cache,
                      child_lfs,
                      _lfs_cache->container_index(child_lfs.localIndex(0)).back()
                      );

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }


        template<typename ChildLFS, typename LC>
        void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
        {
          read(child_lfs,local_container);
        }

        template<typename ChildLFS, typename LC>
        void write_sub_container(const ChildLFS& child_lfs, const LC& local_container)
        {
          write(child_lfs,local_container);
        }

        template<typename ChildLFS, typename LC>
        void add_sub_container(const ChildLFS& child_lfs, const LC& local_container)
        {
          add(child_lfs,local_container);
        }

        void commit()
        {
        }

        ElementType& operator[](size_type i)
        {
          return (*_container)[_lfs_cache->container_index(i)];
        }

        const ElementType& operator[](size_type i) const
        {
          return (*_container)[_lfs_cache->container_index(i)];
        }

        ElementType& operator[](const DOFIndex& di)
        {
          return (*_container)[_lfs_cache->container_index(di)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return (*_container)[_lfs_cache->container_index(di)];
        }

        ElementType& operator[](const ContainerIndex& ci)
        {
          return (*_container)[ci];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci];
        }


        /*
        ElementType& operator[](size_type i)
        {
          auto& ci = _lfs_cache->container_index(i);
          return child(ci.back()).entry(i - _lfs_cache->offset(ci.back()));
        }

        const ElementType& operator[](size_type i) const
        {
          auto& ci = _lfs_cache->container_index(i);
          return child(ci.back()).entry(i - _lfs_cache->offset(ci.back()));
        }

        ElementType& operator[](const DOFIndex& di)
        {
          auto& ci = _lfs_cache->container_index(i);
          return (*this)[ci];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          auto& ci = _lfs_cache->container_index(i);
          return (*this)[ci];
        }

        ElementType& operator[](const ContainerIndex& ci)
        {
          size_type i;
          bool extended_block;
          std::tie(i,extended_block) = _lfs_cache->local_index(ci);
          if (extended_block)
            return child(ci.back()).extended_entry(i - _lfs_cache->extended_offset(ci.back()));
          else
            return child(ci.back()).entry(i - _lfs_cache->offset(ci.back()));
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          size_type i;
          bool extended_block;
          std::tie(i,extended_block) = _lfs_cache->local_index(ci);
          if (extended_block)
            return child(ci.back()).extended_entry(i - _lfs_cache->extended_offset(ci.back()));
          else
            return child(ci.back()).entry(i - _lfs_cache->offset(ci.back()));
        }
        */

        HeterogeneousVectorContainer& container()
        {
          return *_container;
        }

        const HeterogeneousVectorContainer& container() const
        {
          return *_container;
        }

        const LFSCache& cache() const
        {
          return *_lfs_cache;
        }

      private:

        HeterogeneousVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      template<typename LFSCache>
      struct ConstLocalView
      {

        dune_static_assert((is_same<typename LFSCache::LocalFunctionSpace::Traits::GridFunctionSpace,GFS>::value),
                           "The LocalFunctionSpace passed to LocalView must belong to the underlying GridFunctionSpace.");

        typedef E ElementType;
        typedef typename LFSCache::LocalFunctionSpace LFS;
        typedef LFS LocalFunctionSpace;
        typedef typename LFS::Traits::DOFIndex DOFIndex;
        typedef typename LFS::Traits::GridFunctionSpace::Ordering::Traits::ContainerIndex ContainerIndex;

        ConstLocalView()
          : _container(nullptr)
          , _lfs_cache(nullptr)
        {}

        ConstLocalView(const HeterogeneousVectorContainer& container)
          : _container(&container)
          , _lfs_cache(nullptr)
        {}

        void attach(const HeterogeneousVectorContainer& container)
        {
          _container = &container;
        }

        void detach()
        {
          _container = nullptr;
        }

        void bind(const LFSCache& lfs_cache)
        {
          _lfs_cache = &lfs_cache;
        }

        void unbind()
        {
        }

        size_type size() const
        {
          return _lfs_cache->size();
        }

        template<typename LC>
        void read(LC& local_container) const
        {
          read_data<
            LC,
            const Container,
            LFSCache
            > visitor(local_container,*_container,*_lfs_cache);

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename ChildLFS, typename LC>
        void read(const ChildLFS& child_lfs, LC& local_container) const
        {
          if (child_lfs.size() == 0)
            return;

          read_data_filtered<
            LC,
            const Container,
            LFSCache,
            ChildLFS
            > visitor(local_container,
                      *_container,
                      *_lfs_cache,
                      child_lfs,
                      _lfs_cache->container_index(child_lfs.local_index(0)).back()
                      );

          TypeTree::applyToTree(_lfs_cache->localFunctionSpace(),visitor);
        }

        template<typename ChildLFS, typename LC>
        void read_sub_container(const ChildLFS& child_lfs, LC& local_container) const
        {
          read(child_lfs,local_container);
        }


        const ElementType& operator[](size_type i) const
        {
          return (*_container)[_lfs_cache->container_index(i)];
        }

        const ElementType& operator[](const DOFIndex& di) const
        {
          return (*_container)[_lfs_cache->container_index(di)];
        }

        const ElementType& operator[](const ContainerIndex& ci) const
        {
          return (*_container)[ci];
        }

        const HeterogeneousVectorContainer& container() const
        {
          return *_container;
        }


      private:

        const HeterogeneousVectorContainer* _container;
        const LFSCache* _lfs_cache;

      };


      HeterogeneousVectorContainer(const HeterogeneousVectorContainer& rhs)
        : _gfs(rhs._gfs)
        , _container(make_shared<Container>())
      {
        init_accessors();
        allocate_vector_blocks();
        (*_container) = rhs.base();
      }

      HeterogeneousVectorContainer(const GFS& gfs)
        : _gfs(gfs)
        , _container(make_shared<Container>())
      {
        init_accessors();
        allocate_vector_blocks();
      }

      HeterogeneousVectorContainer(const GFS& gfs, const E& e)
        : _gfs(gfs)
        , _container(make_shared<Container>())
      {
        init_accessors();
        allocate_vector_blocks();
        (*_container)=e;
      }

      void detach()
      {
        _container.reset();
      }

      void attach(shared_ptr<Container> container)
      {
        _container = container;
      }

      bool attached() const
      {
        return bool(_container);
      }

      shared_ptr<Container> storage()
      {
        return _container;
      }

      size_type N() const
      {
        return _container->N();
      }

      HeterogeneousVectorContainer& operator= (const HeterogeneousVectorContainer& r)
      {
        (*_container) = r.base();
        return *this;
      }

      HeterogeneousVectorContainer& operator= (const E& e)
      {
        (*_container)=e;
        return *this;
      }

      HeterogeneousVectorContainer& operator*= (const E& e)
      {
        (*_container)*=e;
        return *this;
      }


      HeterogeneousVectorContainer& operator+= (const E& e)
      {
        (*_container)+=e;
        return *this;
      }

      HeterogeneousVectorContainer& operator+= (const HeterogeneousVectorContainer& e)
      {
        (*_container) += e.base();
        return *this;
      }

      HeterogeneousVectorContainer& operator-= (const HeterogeneousVectorContainer& e)
      {
        (*_container) -= e.base();
        return *this;
      }

      /*
      block_type& block(std::size_t i)
      {
        return (*_container)[i];
      }

      const block_type& block(std::size_t i) const
      {
        return (*_container)[i];
      }
      */

      E& operator[](const ContainerIndex& ci)
      {
        return block(ci.back()).access_element(*_container,ci);
      }

      const E& operator[](const ContainerIndex& ci) const
      {
        return block(ci.back()).access_element(*_container,ci);
      }

      typename Dune::template FieldTraits<E>::real_type two_norm() const
      {
        return _container->two_norm();
      }

      /*
      typename Dune::template FieldTraits<E>::real_type one_norm() const
      {
        return _container->one_norm();
      }

      typename Dune::template FieldTraits<E>::real_type infinity_norm() const
      {
        return _container->infinity_norm();
      }
      */

      E operator*(const HeterogeneousVectorContainer& y) const
      {
        return (*_container)*y.base();
      }

      E dot(const HeterogeneousVectorContainer& y) const
      {
        return _container->dot(y.base());
      }


      HeterogeneousVectorContainer& axpy(const E& a, const HeterogeneousVectorContainer& y)
      {
        _container->axpy(a, y.base());
        return *this;
      }

      // for debugging and AMG access
      Container& base ()
      {
        return *_container;
      }

      const Container& base () const
      {
        return *_container;
      }

      operator Container&()
      {
        return *_container;
      }

      operator const Container&() const
      {
        return *_container;
      }

      /*
      iterator begin()
      {
        return _container->begin();
      }


      const_iterator begin() const
      {
        return _container->begin();
      }

      iterator end()
      {
        return _container->end();
      }


      const_iterator end() const
      {
        return _container->end();
      }

      size_t flatsize() const
      {
        return _container->dim();
      }
      */

      /*
      template<typename X>
      void std_copy_to (std::vector<X>& x) const
      {
        // FIXME: do this hierachically
        size_t n = flatsize();
        x.resize(n);
        for (size_t i=0; i<n; i++)
          x[i] = (*_container)[i][i];
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        // FIXME: do this hierachically
        //test if x has the same size as the container
        assert (x.size() == flatsize());
        for (size_t i=0; i<flatsize(); i++)
          (*_container)[i][i] = x[i];
      }
      */

    private:


      void allocate_vector_blocks()
      {
        for (size_type i = 0; i < Container::block_count; ++i)
          block(i).allocate(*_container,_gfs.ordering());
      }

      static const DynamicBlockAccessorBase<GFS,Container>& block(size_type i)
      {
        return *_blocks[i];
      }

      static void init_accessors()
      {
        if (!_blocks[0])
          create_accessors(TypeTree::index_range<Container::block_count>());
      }

      template<std::size_t... i>
      static void create_accessors(TypeTree::index_pack<i...>)
      {
        TypeTree::discard((create_accessor<i>(),0)...);
      }

      template<std::size_t i>
      static void create_accessor()
      {
        static DynamicBlockAccessor<GFS,Container,i> accessor;
        _blocks[i] = &accessor;
      }

      const GFS& _gfs;
      shared_ptr<Container> _container;
      static const DynamicBlockAccessorBase<GFS,Container>* _blocks[Container::block_count];
    };

    template<typename GFS, typename Container>
    const DynamicBlockAccessorBase<GFS,Container>*
    HeterogeneousVectorContainer<GFS,Container>::_blocks[Container::block_count] = { nullptr };


    namespace istl {

      // TMPs for deducing heterogeneous vector structure from GFS backends

      template<typename E>
      struct extract_vector
      {

        template<typename Node, typename TreePath>
        struct doVisit
        {
          // visit only root and direct children (to extract ISTL vector types)
          static const bool value = TypeTree::TreePathSize<TreePath>::value < 2;
        };

        template<typename Node>
        struct extract_istl_vector
        {

          typedef typename TypeTree::AccumulateType<
            Node,
            vector_creation_policy<E>
            >::type vector_descriptor;

          // calculate the ISTL vector type for the current child
          typedef typename vector_descriptor::vector_type type;

        };

        template<typename Node>
        struct extract_heterogeneous_vector
        {

          // export void as special marker type for the root node
          // - required in parent-child reduction
          typedef void type;
        };

        template<typename Node, typename TreePath>
        struct visit
          : public std::conditional<TypeTree::TreePathSize<TreePath>::value == 0,
                                    extract_heterogeneous_vector<Node>,
                                    extract_istl_vector<Node>
                                    >::type
        {};

      };

      // sibling reduction functor
      struct combine_istl_blocks
      {

        // We have to partially specialize this thing, as the first argument is a
        // tuple and we have to get at its arguments
        template<typename D1, typename D2>
        struct reduce;

      };

      template<typename D2, typename... D1>
      struct combine_istl_blocks::reduce<tuple<D1...>,D2>
      {
        typedef tuple<D1...,D2> type;
      };

      template<typename D2, typename... D1>
      struct combine_istl_blocks::reduce<tuple<D1...>,tuple<D2> >
      {
        typedef tuple<D1...,D2> type;
      };


      // sibling reduction functor
      struct build_heterogeneous_vector
      {

        // We have to partially specialize this thing, as the first argument is a
        // tuple and we have to get at its arguments
        template<typename Blocks, typename Root>
        struct reduce;

      };

      template<typename Root, typename... Blocks>
      struct build_heterogeneous_vector::reduce<tuple<Blocks...>,Root>
      {

        typedef typename std::conditional<
          is_same<Root,void>::value,
          HeterogeneousVector<Blocks...>,
          tuple<Root>
          >::type type;

      };

      // policy describing the GFS tree -> ISTL vector reduction
      template<typename E>
      struct heterogeneous_vector_creation_policy
        : public TypeTree::TypeAccumulationPolicy<extract_vector<E>,
                                                  combine_istl_blocks,
                                                  tuple<>,
                                                  build_heterogeneous_vector,
                                                  TypeTree::bottom_up_reduction>
      {};

    } // namespace istl

    // helper struct invoking the GFS tree -> ISTL vector reduction
    template<typename GFS, typename E>
    struct HeterogeneousVectorSelectorHelper
    {

      typedef typename TypeTree::AccumulateType<
        GFS,
        istl::heterogeneous_vector_creation_policy<E>
        >::type vector_type;

      typedef HeterogeneousVectorContainer<GFS,vector_type> Type;

    };

    template<typename GFS, typename E>
    struct BackendVectorSelectorHelper<HeterogeneousVectorBackend,GFS,E>
      : public HeterogeneousVectorSelectorHelper<GFS,E>
    {};

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_HETEROGENEOUSVECTORBACKEND_HH
