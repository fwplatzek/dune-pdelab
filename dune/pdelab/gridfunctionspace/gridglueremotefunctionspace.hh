#ifndef DUNE_PDELAB_GRIDGLUEGRIDREMOTEFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDGLUEGRIDREMOTEFUNCTIONSPACE_HH

#include <dune/common/nullptr.hh>
#include <dune/localfunctions/monom.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/utility.hh>
#include <dune/pdelab/ordering/transformations.hh>

#include "gridgluetags.hh"

namespace Dune {
  namespace PDELab {

    //! Trait class for the GridGlue grid function space
    template<typename GFS>
    struct RemoteLeafFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = false,
      };

      typedef GFS RemoteGridFunctionSpace;

      //! \brief the grid view where grid function is defined upon
      typedef typename GFS::Traits::GridView GridViewType;
      typedef typename GFS::Traits::GridView GridView;

      //! vector backend
      typedef typename GFS::Traits::Backend BackendType;
      typedef typename GFS::Traits::Backend Backend;

      typedef typename GFS::Traits::SizeType SizeType;
      typedef /* typename GFS::Traits::OrderingTag */
      LexicographicOrderingTag OrderingTag;

    };

    template<typename GFS
#ifdef STORE_GRIDGLUE
             , typename GG
#endif
             >
    class RemoteLeafFunctionSpace
      : public TypeTree::LeafNode
      , public GridFunctionSpaceBase<
                 RemoteLeafFunctionSpace<GFS>,
                 RemoteLeafFunctionSpaceTraits<GFS>
                 >
    {
      typedef TypeTree::TransformTree<RemoteLeafFunctionSpace,gfs_to_ordering<RemoteLeafFunctionSpace> > ordering_transformation;

      typedef GridFunctionSpaceBase<
        RemoteLeafFunctionSpace<GFS>,
        RemoteLeafFunctionSpaceTraits<GFS>
        > BaseT;
    public:

//       RemoteLeafFunctionSpace(const RemoteLeafFunctionSpace& other)
//         : BaseT(other)
// #ifdef STORE_GRIDGLUE
//         , gridglue_(other.gridGlue())
// #endif
//       {
//       }

      template<typename Transformation>
      RemoteLeafFunctionSpace(const GFS& gfs, const Transformation& t)
        : BaseT(gfs.backend(),
          LexicographicOrderingTag())
          // gfs.orderingTag())
        , remoteGfs_(stackobject_to_shared_ptr(gfs))
#ifdef STORE_GRIDGLUE
        , gridglue_(t.gridglue)
#endif
      {}

      template<typename Transformation>
      RemoteLeafFunctionSpace(shared_ptr<const GFS> pgfs, const Transformation& t)
        : BaseT(pgfs->backend(),
          LexicographicOrderingTag())
          // pgfs->orderingTag())
        , remoteGfs_(pgfs)
#ifdef STORE_GRIDGLUE
        , gridglue_(t.gridglue)
#endif
      {}

      //! export Traits class
      typedef RemoteLeafFunctionSpaceTraits<GFS> Traits;
      typedef typename GFS::SizeTag SizeTag;
      typedef typename GFS::OrderingTag OrderingTag;
      typedef RemoteLeafFunctionSpaceTag ImplementationTag;

      typedef typename ordering_transformation::Type Ordering;

      shared_ptr<const GFS> remoteGfs() const { return remoteGfs_; }

    private:

      shared_ptr<const GFS> remoteGfs_;

#ifdef STORE_GRIDGLUE
      typedef GG GridGlue;

      const GridGlue & gridGlue() const { return gridglue_; }

    private:
      const GridGlue & gridglue_;
#endif
    };

    // Register LeafGFS -> RemoteLeafGFS transformation
    template<typename GridFunctionSpace, typename Params>
    Dune::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_remote_gfs<Params>,
      RemoteLeafFunctionSpace<GridFunctionSpace>
      >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_remote_gfs<Params>* t, LeafGridFunctionSpaceTag* tag);

    // register PowerGFS -> RemoteLocalFunctionSpace transformation
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_remote_gfs_template
    {
      // typedef typename Transformation::GridGlue Params;
      template<typename TC>
      struct result
      {
        typedef PowerGridFunctionSpace<
          typename Dune::TypeTree::TransformTree<TC,
                                                 Dune::PDELab::gfs_to_remote_gfs<void> >::Type,
          SourceNode::CHILDREN,
          typename SourceNode::Backend,
          typename SourceNode::OrderingTag> type;
      };
    };
    template<typename PowerGridFunctionSpace, typename Params>
    Dune::TypeTree::TemplatizedGenericPowerNodeTransformation<
      PowerGridFunctionSpace,
      gfs_to_remote_gfs<Params>,
      power_gfs_to_remote_gfs_template<PowerGridFunctionSpace,gfs_to_remote_gfs<Params> >::template result
      >
    registerNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_remote_gfs<Params>* t, PowerGridFunctionSpaceTag* tag);

    // register CompositeGFS -> RemoteLocalFunctionSpace transformation (variadic version)
    template<typename SourceNode, typename Transformation>
    struct variadic_composite_gfs_to_remote_gfs_template
    {
      // typedef typename Transformation::GridGlue Params;
      template<typename... TC>
      struct result
      {
        typedef CompositeGridFunctionSpace<typename SourceNode::Backend,
                                           typename SourceNode::OrderingTag,
                                           typename Dune::TypeTree::TransformTree<TC,
                                                                                  Dune::PDELab::gfs_to_remote_gfs<void> >::Type...> type;
      };
    };
    template<typename CompositeGridFunctionSpace, typename Params>
    Dune::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_remote_gfs<Params>,
      variadic_composite_gfs_to_remote_gfs_template<CompositeGridFunctionSpace,gfs_to_remote_gfs<Params> >::template result
      >
    registerNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_remote_gfs<Params>* t, CompositeGridFunctionSpaceTag* tag);

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEGRIDREMOTEFUNCTIONSPACE_HH
