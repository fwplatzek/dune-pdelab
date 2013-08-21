#ifndef DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH
#define DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH

namespace Dune {
  namespace PDELab {

    struct NoGFS {};
    struct NoIndex {};
    struct gfs_to_remote_lfs {};

    template<>
    struct gfs_to_lfs<NoGFS> {
      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef NoIndex DOFIndex;
    };

    template<typename Element>
    struct RemoteLFSComputeSizeVisitor :
      public ComputeSizeVisitor<Element>
    {
      typedef typename ComputeSizeVisitor<Element>::Entity Entity;

      template<typename Node, typename TreePath>
      void leaf(Node& node, TreePath treePath)
      {
        node.offset = this->offset;
        // Node::FESwitch::setStore(node.pfe, node.finiteElement());
        // node.n = Node::FESwitch::basis(*node.pfe).size();
        this->offset += node.n;
      }

      RemoteLFSComputeSizeVisitor(const Entity& entity, std::size_t offset = 0) :
        ComputeSizeVisitor<Element>(entity,offset)
      {}
    };

    template<typename LFS /*, typename GG */>
    class RemoteLFS : public LFS
    {
    public:
      template <typename... T>
      RemoteLFS(T&&... t) : LFS(std::forward<T>(t)...) {}

      size_t maxLocalSize () { return this->localSize(); };

      void bind() const;

      /*
      DGFEM finiteElement() const { ... };

    private:
      template<>
      void bind (NodeType& node,
        const typename GG:RemoteIntersection& e)
      {
        typedef typename LocalFunctionSpaceBaseNode<GFS,DOFIndex>::Traits::Element Element;
        assert(&node == this);

        // compute sizes
        RemoteLFSComputeSizeVisitor<Element> csv(e);
        TypeTree::applyToTree(node,csv);
      }
      */
    };

    // Register LeafGFS -> RemoteLFS transformation
    template<typename GridFunctionSpace>
    Dune::PDELab::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_remote_lfs,
      RemoteLFS<typename TypeTree::TransformTree<GridFunctionSpace, gfs_to_lfs<NoGFS> >::Type>
      >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_remote_lfs* t, LeafGridFunctionSpaceTag* tag);

    // transformation template, we need a custom template in order to inject the DOFIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_remote_lfs_template
    {
      template<typename TC>
      struct result
      {
        typedef RemoteLFS< PowerLocalFunctionSpaceNode<SourceNode,NoIndex,TC,SourceNode::CHILDREN> > type;
      };
    };

    template<typename PowerGridFunctionSpace>
    Dune::PDELab::TypeTree::TemplatizedGenericPowerNodeTransformation<
      PowerGridFunctionSpace,
      gfs_to_remote_lfs,
      power_gfs_to_remote_lfs_template<PowerGridFunctionSpace,gfs_to_remote_lfs>::template result
      >
    registerNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_remote_lfs* t, PowerGridFunctionSpaceTag* tag);

    // transformation template, we need a custom template in order to inject the MultiIndex type into the LocalFunctionSpace
    template<typename SourceNode, typename Transformation>
    struct variadic_composite_gfs_to_remote_lfs_template
    {
      template<typename... TC>
      struct result
      {
        typedef RemoteLFS< CompositeLocalFunctionSpaceNode<SourceNode,NoIndex,TC...> > type;
      };
    };

    // register CompositeGFS -> LocalFunctionSpace transformation (variadic version)
    template<typename CompositeGridFunctionSpace>
    Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_remote_lfs,
      variadic_composite_gfs_to_remote_lfs_template<CompositeGridFunctionSpace,gfs_to_remote_lfs>::template result
      >
    registerNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_remote_lfs* t, CompositeGridFunctionSpaceTag* tag);

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH
