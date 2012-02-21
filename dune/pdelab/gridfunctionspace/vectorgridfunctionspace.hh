// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH
#define DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH

#include <cstddef>

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/powernode.hh>
#include <dune/pdelab/gridfunctionspace/lexicographicordering.hh>
#include <dune/pdelab/gridfunctionspace/orderingbase.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // power grid function space
    //=======================================

    /** \brief base class for tuples of grid function spaces
        product of identical grid function spaces
        base class that holds implementation of the methods

        PGFS(T,k) = {T}^k

        \tparam T the underlying are all grid function spaces
        \tparam k power factor
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
    */
    template<typename T, std::size_t k,
             typename Backend,
             typename OrderingTag = LexicographicOrderingTag>
    class VectorGridFunctionSpace :
      public TypeTree::PowerNode<T,k>,
      public PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace<T, k, Backend, OrderingTag>,
        typename T::Traits::GridViewType,
        Backend,
        OrderingTag,
        k
      >
    {

    public:

      typedef VectorGridFunctionSpaceTag ImplementationTag;

      typedef TypeTree::PowerNode<T,k> BaseT;

      typedef PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        typename T::Traits::GridViewType,
        Backend,
        OrderingTag,
        k
        > ImplementationBase;

      friend class PowerCompositeGridFunctionSpaceBase<
        VectorGridFunctionSpace,
        typename T::Traits::GridViewType,
        Backend,
        OrderingTag,
        k>;

      typedef TypeTree::TransformTree<VectorGridFunctionSpace,
                                      gfs_to_ordering<VectorGridFunctionSpace>
                                      > ordering_transformation;

    public:

      typedef typename ordering_transformation::Type Ordering;

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      VectorGridFunctionSpace(T& c, const Backend& backend = Backend())
        : BaseT(c)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               const Backend& backend = Backend())
        : BaseT(c0,c1)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               T& c5,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4,c5)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               T& c5,
                               T& c6,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               T& c5,
                               T& c6,
                               T& c7,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               T& c5,
                               T& c6,
                               T& c7,
                               T& c8,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
        , ImplementationBase(backend)
      {
      }

      VectorGridFunctionSpace (T& c0,
                               T& c1,
                               T& c2,
                               T& c3,
                               T& c4,
                               T& c5,
                               T& c6,
                               T& c7,
                               T& c8,
                               T& c9,
                               const Backend& backend = Backend())
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
        , ImplementationBase(backend)
      {
      }

      shared_ptr<Ordering> ordering() const
      {
        if (!_ordering)
          _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
        return _ordering;
      }

    private:

      mutable shared_ptr<Ordering> _ordering;

    };

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_VECTORGRIDFUNCTIONSPACE_HH
