#ifndef DUNE_PDELAB_GRIDGLUELFS_HH
#define DUNE_PDELAB_GRIDGLUELFS_HH

#include "localfunctionspace.hh"
#include "gridgluegridfunctionspace.hh"
#include <dune/grid-glue/adapter/gridglue.hh>

namespace Dune {
  namespace PDELab {

    // template<typename GG>
    // struct gfs_to_remote_lfs {
    //   typedef GG GridGlue;
    //   GridGlue _gridGlue;
    // };

    template<typename GFS,
             typename DI>
    struct GridGlueLocalFunctionSpaceBaseTraits
    {
      //! \brief Type of the underlying grid function space
      typedef GFS GridFunctionSpaceType;

      //! \brief Type of the underlying grid function space
      typedef GFS GridFunctionSpace;

      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridViewType;

      //! \brief Type of the grid view that the underlying grid function space is defined on.
      typedef typename GFS::Traits::GridViewType GridView;

      //! \brief Type to store indices from Backend
      typedef typename GFS::Traits::SizeType SizeType;

      //! \brief Type of container to store indices
      typedef typename std::vector<SizeType> IndexContainer;

      //! \brief Type of MultiIndex associated with this LocalFunctionSpace.
      typedef DI DOFIndex;

      //! \brief Type of container to store multiindices.
      typedef typename std::vector<DI> DOFIndexContainer;
    };

    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename LFS0, typename LFS1>
    class GridGlueLocalFunctionSpaceNode :
      // TODO: Mixin RemoteLFSes
      public TypeTree::VariadicCompositeNode<LFS0,LFS1>
      , public LocalFunctionSpaceBaseNode<GFS,DOFIndex>
    {
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef TypeTree::VariadicCompositeNode<LFS0,LFS1> TreeNode;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

    public:
      typedef GridGlueLocalFunctionSpaceBaseTraits<GFS,DOFIndex> Traits;

      typedef PowerLocalFunctionSpaceTag ImplementationTag;

      //! \brief initialize with grid function space
      template<typename Transformation>
      GridGlueLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
                                      const Transformation& t,
                                      Dune::shared_ptr<LFS0> lfs0,
                                      Dune::shared_ptr<LFS1> lfs1)
        : TreeNode(lfs0,lfs1)
        , BaseT(gfs)
      {}

      template<typename Transformation>
      GridGlueLocalFunctionSpaceNode (const GFS& gfs,
                                      const Transformation& t,
                                      Dune::shared_ptr<LFS0> lfs0,
                                      Dune::shared_ptr<LFS1> lfs1)
        : TreeNode(lfs0,lfs1)
        , BaseT(stackobject_to_shared_ptr(gfs))
      {}

      // //! \brief initialize with grid function space
      // template<typename Transformation, typename... DUMMY>
      // GridGlueLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
      //                                 const Transformation& t,
      //                                 Dune::shared_ptr<LFS0> lfs0,
      //                                 Dune::shared_ptr<LFS1> lfs1,
      //                                 DUMMY&... other)
      //   : TreeNode(lfs0,lfs1)
      //   , BaseT(gfs)
      // {}

      // template<typename Transformation, typename... DUMMY>
      // GridGlueLocalFunctionSpaceNode (const GFS& gfs,
      //                                 const Transformation& t,
      //                                 Dune::shared_ptr<LFS0> lfs0,
      //                                 Dune::shared_ptr<LFS1> lfs1,
      //                                 DUMMY&... other)
      //   : TreeNode(lfs0,lfs1)
      //   , BaseT(stackobject_to_shared_ptr(gfs))
      // {}

    public:
      //! \brief bind local function space to one of the GridGlue contextes (sub-domain cell or remote intersection)
      template<typename Context>
      void bind(const Context & c) const
      {
        this->bind(*this, this->template child<0>(), c);
        this->bind(*this, this->template child<1>(), c);
        // These are the RemoteLocalFunctionSpace sub-trees
        // this->bind(*this, this->template child<2>(), c);
        // this->bind(*this, this->template child<3>(), c);
      }

    private:
      typedef typename GFS::Traits::GridGlue GridGlue;
      typedef integral_constant<int,0> Patch0Tag;
      typedef integral_constant<int,1> Patch1Tag;
      typedef typename GridGlue::Intersection CouplingIntersection;
      template<typename NodeType, typename ChildNodeType>
      void bind (NodeType& node, ChildNodeType& child,
        const CouplingIntersection& rit)
      {
        assert(&node == this);
        // compute sizes
        RemoteLFSComputeSizeVisitor<CouplingIntersection> csv(rit);
        TypeTree::applyToTree(node,csv);
      }
      template<typename NodeType, typename ChildNodeType>
      void bind (NodeType& node, ChildNodeType& child,
        const typename ChildNodeType::Traits::Element& e)
      {
#if 0
        assert(&node == this);

        typedef typename Traits::Element Elementy;

        // compute sizes
        ComputeSizeVisitor<Element> csv(e);
        TypeTree::applyToTree(node,csv);

        // initialize iterators and fill indices
        FillIndicesVisitor<Element> fiv(e);
        TypeTree::applyToTree(node,fiv);
#endif
      }
    };

    // register GridGlueGFS -> LocalFunctionSpace transformation (variadic version)
    template<typename GridGlueGridFunctionSpace, typename Params>
    struct GridGlueGridFunctionSpaceToLocalFunctionSpaceNodeTransformation
    {
      typedef gfs_to_lfs<Params> Transformation;

      static const bool recursive = true;

      template<typename TC0, typename TC1, typename... DUMMY>
      struct result
      {
        typedef GridGlueLocalFunctionSpaceNode<GridGlueGridFunctionSpace,
                                               typename Transformation::DOFIndex,
                                               TC0,TC1> type;
        typedef shared_ptr<type> storage_type;
      };

      template<typename TC0, typename TC1, typename... DUMMY>
      static typename result<TC0,TC1>::type transform(const GridGlueGridFunctionSpace& s,
        const Transformation& t, shared_ptr<TC0> c0, shared_ptr<TC1> c1, DUMMY&&...)
      {
        return typename result<TC0,TC1>::type(s,t,c0,c1);
      }

      template<typename TC0, typename TC1, typename... DUMMY>
      static typename result<TC0,TC1>::storage_type transform_storage(shared_ptr<const GridGlueGridFunctionSpace> s,
        const Transformation& t, shared_ptr<TC0> c0, shared_ptr<TC1> c1, DUMMY&&...)
      {
        return make_shared<typename result<TC0,TC1>::type>(s,t,c0,c1);
      }

    };
    template<typename GridGlueGridFunctionSpace, typename Params>
    GridGlueGridFunctionSpaceToLocalFunctionSpaceNodeTransformation<GridGlueGridFunctionSpace, Params>
    registerNodeTransformation(GridGlueGridFunctionSpace* cgfs, gfs_to_lfs<Params>* t, GridGlueGridFunctionSpaceTag* tag);

  } // end PDELab
} // end Dune

#endif // DUNE_PDELAB_GRIDGLUELFS_HH
