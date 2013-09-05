#ifndef DUNE_PDELAB_GRIDGLUEGRIDREMOTEFUNCTIONSPACE_HH
#define DUNE_PDELAB_GRIDGLUEGRIDREMOTEFUNCTIONSPACE_HH

#include <dune/common/nullptr.hh>
#include <dune/localfunctions/monom.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/utility.hh>

#include "gridgluetags.hh"

namespace Dune {
  namespace PDELab {

    //! Trait class for the GridGlue grid function space
    template<typename GFS>
    struct RemoteFunctionSpaceTraits
    {
      enum{
        //! \brief True if this grid function space is composed of others.
        isComposite = false,
      };

      //! \brief the grid view where grid function is defined upon
      typedef void GridViewType;
      typedef void GridView;

      //! vector backend
      typedef typename GFS::Traits::Backend BackendType;
      typedef typename GFS::Traits::Backend Backend;

    };

    template<typename GFS>
    class RemoteLeafFunctionSpace
      : public TypeTree::LeafNode
      // , public GridFunctionSpaceBase<
      //            RemoteLeafFunctionSpace<GFS>,
      //            RemoteLeafFunctionSpaceTraits<GFS>
      //            >
    {
      typedef TypeTree::TransformTree<RemoteLeafFunctionSpace,gfs_to_ordering<RemoteLeafFunctionSpace> > ordering_transformation;
    public:

      //! export Traits class
      typedef RemoteFunctionSpaceTraits<GFS> Traits;
      typedef typename GFS::SizeTag SizeTag;
      typedef typename GFS::OrderingTag OrderingTag;
      typedef LeafGridFunctionSpaceTag ImplementationTag;

      typedef typename ordering_transformation::Type Ordering;
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
      typedef typename Transformation::GridFunctionSpace Params;
      template<typename TC>
      struct result
      {
        typedef PowerGridFunctionSpace<
          typename Dune::TypeTree::TransformTree<TC,
                                                 Dune::PDELab::gfs_to_remote_gfs<TC> >::Type,
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
      typedef typename Transformation::GridFunctionSpace Params;
      template<typename... TC>
      struct result
      {
        typedef CompositeGridFunctionSpace<typename SourceNode::Backend,
                                           typename SourceNode::OrderingTag,
                                           typename Dune::TypeTree::TransformTree<TC,
                                                                                  Dune::PDELab::gfs_to_remote_gfs<TC> >::Type...> type;
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
