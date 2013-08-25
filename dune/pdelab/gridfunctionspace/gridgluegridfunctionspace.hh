// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDGLUEGRIDFUNTIONSPACE_HH
#define DUNE_PDELAB_GRIDGLUEGRIDFUNTIONSPACE_HH

#include "gridglueremotelfs.hh"

#include <dune/common/shared_ptr.hh>

#include <dune/pdelab/common/typetree/compositenodemacros.hh>
#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>

namespace Dune {
  namespace PDELab {

    //=======================================
    // composite grid function space
    //=======================================

    //! \addtogroup GridFunctionSpace grid function space
    //! \ingroup PDELab
    //! \{

    struct GridGlueGridFunctionSpaceTag : public CompositeGridFunctionSpaceTag {};

    //! Trait class for the GridGlue grid function space
    template<typename GG, typename GFS1, typename GFS2, typename B, typename O>
    struct GridGlueGridFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = 1,
        //! \brief number of child spaces
        noChilds = 2
      };

      const static std::size_t CHILDREN = 2;

      //! \brief the GridGlue where the coupling takes places
      typedef GG GridGlue;

      //! \brief the grid view where grid function is defined upon
      typedef GridGlue GridViewType;

      typedef GridGlue GridView;

      //! \brief vector backend
      typedef B BackendType;

      typedef B Backend;

      //! \brief mapper
      typedef O MapperType;

      typedef O OrderingTag;

      //! \brief short cut for size type exported by Backend
      typedef typename B::size_type SizeType;
    };

    /** \brief base class for tuples of grid function spaces
        base class that holds implementation of the methods
        this is the default version with lexicographic ordering
        \tparam Mapper is the ordering parameter. Use e.g.
        \link GridFunctionSpaceLexicographicMapper GridFunctionSpaceLexicographicMapper \endlink
        or \link  GridFunctionSpaceComponentBlockwiseMapper  GridFunctionSpaceComponentBlockwiseMapper \endlink
        or \link  GridFunctionSpaceBlockwiseMapper  GridFunctionSpaceBlockwiseMapper \endlink
        or \link  GridFunctionSpaceDynamicBlockwiseMapper  GridFunctionSpaceDynamicBlockwiseMapper \endlink
        \tparam Ti are all grid function spaces
    */
    template<typename GG,
             typename GFS1,
             typename GFS2,
             typename Backend,
             typename OT = LexicographicOrderingTag>
    class GridGlueGridFunctionSpace
      : public TypeTree::VariadicCompositeNode<GFS1,GFS2>
      , public GridFunctionSpaceBase<
                 GridGlueGridFunctionSpace<GG,GFS1,GFS2,Backend,OT>,
                 GridGlueGridFunctionSpaceTraits<GG,GFS1,GFS2,Backend,OT>
                 >
//      , public DataHandleProvider< GridGlueGridFunctionSpace<GFS1,GFS2,Backend,OrderingTag> >
    {
    public:

      typedef GridGlueGridFunctionSpaceTag ImplementationTag;
      typedef OT OrderingTag;

      typedef TypeTree::VariadicCompositeNode<GFS1,GFS2> BaseT;

    private:

      typedef GridFunctionSpaceBase<
                 GridGlueGridFunctionSpace<GG,GFS1,GFS2,Backend,OrderingTag>,
                 GridGlueGridFunctionSpaceTraits<GG,GFS1,GFS2,Backend,OrderingTag>
                 > ImplementationBase;

      template<typename,typename>
      friend class GridFunctionSpaceBase;

    public:

      //! export traits class
      typedef typename ImplementationBase::Traits Traits;

      GridGlueGridFunctionSpace (GG& glue,
                                 GFS1& c0,
                                 GFS2& c1,
                                 const Backend& backend = Backend(),
                                 const OrderingTag ordering_tag = OrderingTag())
        : BaseT(c0,c1)
        , ImplementationBase(backend,ordering_tag)
        , _glue(glue)
      {}

#if 1
    private:
      typedef TypeTree::TransformTree<GridGlueGridFunctionSpace,
                                      gfs_to_ordering<GridGlueGridFunctionSpace>
                                      > ordering_transformation;
    public:

      typedef typename ordering_transformation::Type Ordering;

      //! Direct access to the DOF ordering.
      const Ordering &ordering() const
      {
        return *orderingStorage();
      }

      //! Direct access to the DOF ordering.
      Ordering &ordering()
      {
        return *orderingStorage();
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<const Ordering> orderingStorage() const
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

      //! Direct access to the storage of the DOF ordering.
      shared_ptr<Ordering> orderingStorage()
      {
        if (!this->isRootSpace())
          {
            DUNE_THROW(GridFunctionSpaceHierarchyError,
                       "Ordering can only be obtained for root space in GridFunctionSpace tree.");
          }
        if (!_ordering)
          {
            create_ordering();
            this->update(*_ordering);
          }
        return _ordering;
      }

    private:

      // This method here is to avoid a double update of the Ordering when the user calls
      // GFS::update() before GFS::ordering().
      void create_ordering() const
      {
        _ordering = make_shared<Ordering>(ordering_transformation::transform(*this));
      }

      mutable shared_ptr<Ordering> _ordering;

#endif

    private:
      GG& _glue;

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEGRIDFUNTIONSPACE_HH
