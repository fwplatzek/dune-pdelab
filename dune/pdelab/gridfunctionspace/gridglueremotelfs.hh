#ifndef DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH
#define DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH

#include <dune/common/nullptr.hh>
#include <dune/localfunctions/monom.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/utility.hh>
#include <dune/pdelab/gridfunctionspace/gridgluelocalfunctionspace.hh>

#include "gridgluetags.hh"

namespace Dune {
  namespace PDELab {

    struct NoObject {};

    template<>
    struct gfs_to_lfs<NoObject> {
      typedef NoObject Params;
      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef NoObject DOFIndex;
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
        node.pfe = nullptr;
        node.n = Node::FESwitch::basis(node.finiteElement()).size();
        this->offset += node.n;
      }

      RemoteLFSComputeSizeVisitor(const Element& e,std::size_t offset = 0) :
        ComputeSizeVisitor<Element>(e,offset)
      {}
    };

    template<typename LFS>
    class RemoteLFS : public LFS
    {
    public:
      template <typename... T>
      RemoteLFS(T&&... t) :
        LFS(std::forward<T>(t)...)
      {}

      size_t maxLocalSize () { return this->localSize(); };

      template<typename>
      friend struct ComputeSizeVisitor;

    //   /** \todo extract RIT from GG
    //    */
    //   template<typename RIT>
    //   void bind(const RIT & rit) const
    //   {
    //     this->bind(*this, rit);
    //   }
    // private:
    //   template<typename NodeType, typename RIT>
    //   void bind (NodeType& node,
    //     const RIT& rit)
    //   {
    //     assert(&node == this);

    //     // compute sizes
    //     RemoteLFSComputeSizeVisitor<RIT> csv(rit);
    //     TypeTree::applyToTree(node,csv);
    //   }
    };

    //! traits for single component local function space
    template<typename LFS /*, typename GG */>
    struct RemoteLeafLocalFunctionSpaceTraits :
      public LFS::Traits
    {
      typedef typename LFS::Traits::GridFunctionSpace RemoteGridFunctionSpace;
      typedef typename LFS::Traits::FiniteElementType RemoteFiniteElementType;
      typedef typename LFS::FESwitch RemoteFESwitch;
      typedef BasisInterfaceSwitch<RemoteFESwitch> RemoteBasisSwitch;
      typedef typename RemoteBasisSwitch::Domain DomainType;
      typedef typename RemoteBasisSwitch::Range RangeType;
      enum { dim = RemoteGridFunctionSpace::GridView::dimension-1 };
      enum { order = 3 };
      enum { diffOrder = 1 };

      //! Type of local finite element
      typedef MonomLocalFiniteElement<DomainType, RangeType, dim, order, diffOrder> FiniteElementType;
      typedef FiniteElementType FiniteElement;

      //! \brief Type of constraints engine
      typedef NoConstraints ConstraintsType;
      typedef ConstraintsType Constraints;
    };

    // Register LeafGFS -> RemoteLFS transformation
    template<typename GridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      RemoteLFS<typename TypeTree::TransformTree<GridFunctionSpace, gfs_to_lfs<Params> >::Type>
      >
    registerNodeTransformation(GridFunctionSpace* gfs, gfs_to_remote_lfs<Params>* t, LeafGridFunctionSpaceTag* tag);

    // register PowerGFS -> RemoteLocalFunctionSpace transformation
    template<typename SourceNode, typename Transformation>
    struct power_gfs_to_remote_lfs_template
    {
      typedef typename Transformation::GridFunctionSpace Params;
      template<typename TC>
      struct result
      {
        typedef RemoteLFS< PowerLocalFunctionSpaceNode<SourceNode,typename gfs_to_lfs<Params>::DOFIndex,TC,SourceNode::CHILDREN> > type;
      };
    };
    template<typename PowerGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericPowerNodeTransformation<
      PowerGridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      power_gfs_to_remote_lfs_template<PowerGridFunctionSpace,gfs_to_remote_lfs<Params> >::template result
      >
    registerNodeTransformation(PowerGridFunctionSpace* pgfs, gfs_to_remote_lfs<Params>* t, PowerGridFunctionSpaceTag* tag);

    // register CompositeGFS -> RemoteLocalFunctionSpace transformation (variadic version)
    template<typename SourceNode, typename Transformation>
    struct variadic_composite_gfs_to_remote_lfs_template
    {
      typedef typename Transformation::GridFunctionSpace Params;
      template<typename... TC>
      struct result
      {
        typedef RemoteLFS< CompositeLocalFunctionSpaceNode<SourceNode,typename gfs_to_lfs<Params>::DOFIndex,TC...> > type;
      };
    };
    template<typename CompositeGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      CompositeGridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      variadic_composite_gfs_to_remote_lfs_template<CompositeGridFunctionSpace,gfs_to_remote_lfs<Params> >::template result
      >
    registerNodeTransformation(CompositeGridFunctionSpace* cgfs, gfs_to_remote_lfs<Params>* t, CompositeGridFunctionSpaceTag* tag);

    // register GridGlueGFS -> RemoteLocalFunctionSpace transformation (variadic version)
    template<typename SourceNode, typename Transformation>
    struct gridglue_gfs_to_remote_lfs_template
    {
      typedef typename Transformation::GridFunctionSpace Params;
      template<typename TC0, typename TC1, typename... DUMMY>
      struct result
      {
        typedef RemoteLFS< GridGlueLocalFunctionSpaceNode<SourceNode,Params,TC0,TC1> > type;
      };
    };
    template<typename GridGlueGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      GridGlueGridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      gridglue_gfs_to_remote_lfs_template<GridGlueGridFunctionSpace,gfs_to_remote_lfs<Params> >::template result
      >
    registerNodeTransformation(GridGlueGridFunctionSpace* cgfs, gfs_to_remote_lfs<Params>* t, GridGlueGridFunctionSpaceTag* tag);

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH
