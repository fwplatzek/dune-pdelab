// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
#define DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH

#include <dune/pdelab/backend/tags.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/matrixhelpers2.hh>
#include <dune/pdelab/backend/istl/matrixview.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>

namespace Dune {
  namespace PDELab {

    namespace istl {

      struct MatrixState
      {
        enum type {
          building,
          built
        };

#if HAVE_TEMPLATE_ALIASES

        template<type v>
        using State = integral_constant<type,v>;

#else

        template<type v>
        struct State
          : public integral_constant<type,v>
        {};

#endif

      };

    } // namespace istl

    template<typename GFSV, typename GFSU, typename C>
    class ISTLMatrixContainer
      : public istl::MatrixState
    {

    public:

      typedef typename C::field_type ElementType;
      typedef ElementType E;
      typedef C Container;
      typedef C BaseT;
      typedef typename C::field_type field_type;
      typedef typename C::block_type block_type;
      typedef typename C::size_type size_type;

      typedef GFSU TrialGridFunctionSpace;
      typedef GFSV TestGridFunctionSpace;

      typedef typename GFSV::Ordering RowOrdering;
      typedef typename GFSU::Ordering ColOrdering;

      typedef typename RowOrdering::Traits::ContainerIndex RowIndex;
      typedef typename ColOrdering::Traits::ContainerIndex ColIndex;

      struct Pattern {};

#if HAVE_TEMPLATE_ALIASES

      template<typename RowCache, typename ColCache>
      using LocalView = istl::MatrixView<ISTLMatrixContainer,RowCache,ColCache>;

      template<typename RowCache, typename ColCache>
      using ConstLocalView = istl::ConstMatrixView<const ISTLMatrixContainer,RowCache,ColCache>;

#else

      template<typename RowCache, typename ColCache>
      struct LocalView
        : public UncachedMatrixView<ISTLMatrixContainer,RowCache,ColCache>
      {

        LocalView()
        {}

        LocalView(ISTLMatrixContainer& mc)
          : UncachedMatrixView<ISTLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

      template<typename RowCache, typename ColCache>
      struct ConstLocalView
        : public ConstUncachedMatrixView<const ISTLMatrixContainer,RowCache,ColCache>
      {

        ConstLocalView()
        {}

        ConstLocalView(const ISTLMatrixContainer& mc)
          : ConstUncachedMatrixView<const ISTLMatrixContainer,RowCache,ColCache>(mc)
        {}

      };

#endif // HAVE_TEMPLATE_ALIASES


      struct BuildModeAccessor
      {

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return istl::access_matrix_element2(istl::container_tag(_wrapper.base()),
                                              _wrapper.base(),
                                              _wrapper,
                                              ri,
                                              ci,
                                              ri.size()-1,
                                              ci.size()-1,
                                              0);
        }

        explicit BuildModeAccessor(ISTLMatrixContainer& wrapper)
          : _wrapper(wrapper)
        {}

        ISTLMatrixContainer& _wrapper;

      };


      struct ReadOnlyAccessor
      {

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return istl::access_matrix_element(istl::container_tag(_wrapper.base()),_wrapper.base(),ri,ci,ri.size()-1,ci.size()-1);
        }

        explicit ReadOnlyAccessor(const ISTLMatrixContainer& wrapper)
          : _wrapper(wrapper)
        {}

        const ISTLMatrixContainer& _wrapper;

      };


      struct ReadWriteAccessor
      {

        const E& operator()(const RowIndex& ri, const ColIndex& ci) const
        {
          return istl::access_matrix_element(istl::container_tag(_wrapper.base()),_wrapper.base(),ri,ci,ri.size()-1,ci.size()-1);
        }

        E& operator()(const RowIndex& ri, const ColIndex& ci)
        {
          return istl::access_matrix_element(istl::container_tag(_wrapper.base()),_wrapper.base(),ri,ci,ri.size()-1,ci.size()-1);
        }

        explicit ReadWriteAccessor(ISTLMatrixContainer& wrapper)
          : _wrapper(wrapper)
        {}

        ISTLMatrixContainer& _wrapper;

      };



      template<typename GO>
      explicit ISTLMatrixContainer (const GO& go, typename Container::size_type avg_row_size, double overflow)
        : _state(building)
        , _container(make_shared<Container>(
                       go.testGridFunctionSpace().ordering().blockCount(),
                       go.trialGridFunctionSpace().ordering().blockCount(),
                       avg_row_size,
                       overflow,
                       Container::mymode
                       )
                     )
      {
        std::cout << _container->N() << " x " << _container->M() << std::endl;
      }

      /*
      template<typename GO>
      ISTLMatrixContainer (const GO& go)
        : _state(building)
        , _container(make_shared<Container>())
      {
        allocate_matrix(go.testGridFunctionSpace().ordering(),
                        go.trialGridFunctionSpace().ordering(),
                        pattern,
                        *_container);
        _container = e;
      }
      */

      //! Creates an ISTLMatrixContainer without allocating an underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::unattached_container = tags::unattached_container())
        : _state(building)
      {}

      //! Creates an ISTLMatrixContainer with an empty underlying ISTL matrix.
      explicit ISTLMatrixContainer (tags::attached_container)
        : _state(building)
        , _container(make_shared<Container>())
      {}

      ISTLMatrixContainer(const ISTLMatrixContainer& rhs)
        : _state(rhs._state)
        , _container(make_shared<Container>(*(rhs._container)))
        , _row_ordering(rhs._row_ordering)
        , _col_ordering(rhs._col_ordering)
        , _avg_row_size(rhs._avg_row_size)
        , _overflow(rhs._overflow)
      {}

      ISTLMatrixContainer& operator=(const ISTLMatrixContainer& rhs)
      {
        if (this == &rhs)
          return *this;
        _state = rhs._state;
        if (attached())
          {
            (*_container) = (*(rhs._container));
          }
        else
          {
            _container = make_shared<Container>(*(rhs._container));
          }
        return *this;
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

      const shared_ptr<Container>& storage() const
      {
        return _container;
      }

      size_type N() const
      {
        return _container->N();
      }

      size_type M() const
      {
        return _container->M();
      }

      ISTLMatrixContainer& operator= (const E& e)
      {
        (*_container) = e;
        return *this;
      }

      ISTLMatrixContainer& operator*= (const E& e)
      {
        (*_container) *= e;
        return *this;
      }

      BuildModeAccessor accessor(State<building>)
      {
        return BuildModeAccessor(*this);
      }

      ReadOnlyAccessor accessor(State<built>) const
      {
        return ReadOnlyAccessor(*this);
      }

      ReadWriteAccessor accessor(State<built>)
      {
        return ReadWriteAccessor(*this);
      }


      E& operator()(const RowIndex& ri, const ColIndex& ci)
      {
        if (_state == built)
          return accessor(State<built>())(ri,ci);
        else
          return accessor(State<building>())(ri,ci);
      }

      const E& operator()(const RowIndex& ri, const ColIndex& ci) const
      {
        if (_state == built)
          return accessor(State<built>())(ri,ci);
        else
          DUNE_THROW(Exception,"Read-only access not supported before finalizing pattern creation");
      }


      const Container& base() const
      {
        return *_container;
      }

      Container& base()
      {
        return *_container;
      }

      void flush()
      {}

      void finalize()
      {
        if (_state == built)
          return;
        istl::finalize(istl::container_tag(*_container),true_type(),*_container);
        _state = built;
      }

      void clear_row(const RowIndex& ri, const E& diagonal_entry)
      {
        istl::clear_matrix_row(istl::container_tag(*_container),*_container,ri,ri.size()-1);
        (*this)(ri,ri) = diagonal_entry;
      }

      typename Container::size_type averageRowSize() const
      {
        return _avg_row_size;
      }

      void setAverageRowSize(typename Container::size_type avg_row_size)
      {
        _avg_row_size = avg_row_size;
      }

      double overflow() const
      {
        return _overflow;
      }

      void setOverflow(double overflow)
      {
        _overflow = overflow;
      }

    private:

      istl::MatrixState::type _state;
      shared_ptr<Container> _container;
      shared_ptr<const RowOrdering> _row_ordering;
      shared_ptr<const ColOrdering> _col_ordering;
      typename Container::size_type _avg_row_size;
      double _overflow;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTLMATRIXBACKEND_HH
