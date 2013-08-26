#ifndef DUNE_PDELAB_GRIDGLUELFS_HH
#define DUNE_PDELAB_GRIDGLUELFS_HH

#include "localfunctionspace.hh"
#include "gridgluetags.hh"
#include "gridgluegridfunctionspace.hh"
#include <dune/grid-glue/adapter/gridglue.hh>

namespace Dune {
  namespace PDELab {

    template<typename GFS,
             typename DI>
    struct GridGlueLocalFunctionSpaceBaseTraits
      : public LocalFunctionSpaceBaseTraits<GFS,DI>
    {
      typedef typename GFS::Traits::GridGlue GridGlue;
    };

// #define GRIDGLUELFSMIXIN
    // local function space for a power grid function space
    template<typename GFS, typename DOFIndex, typename LFS0, typename LFS1>
    class GridGlueLocalFunctionSpaceNode :
      public TypeTree::VariadicCompositeNode<LFS0,
                                             LFS1
#ifdef GRIDGLUELFSMIXIN
                                             ,
                                             // Mixin RemoteLFSes
                                             typename Dune::PDELab::TypeTree::TransformTree<typename GFS::template Child<0>::Type,
                                                                                            Dune::PDELab::gfs_to_remote_lfs>::Type,
                                             typename Dune::PDELab::TypeTree::TransformTree<typename GFS::template Child<1>::Type,
                                                                                            Dune::PDELab::gfs_to_remote_lfs>::Type
#endif
                                             >
      , public LocalFunctionSpaceBaseNode<GFS,DOFIndex>
    {
      typedef GridGlueLocalFunctionSpaceNode<GFS, DOFIndex, LFS0, LFS1> This;
      typedef LocalFunctionSpaceBaseNode<GFS,DOFIndex> BaseT;
      typedef TypeTree::VariadicCompositeNode<LFS0,
                                              LFS1
#ifdef GRIDGLUELFSMIXIN
                                              ,
                                              typename Dune::PDELab::TypeTree::TransformTree<typename GFS::template Child<0>::Type,
                                                                                             Dune::PDELab::gfs_to_remote_lfs>::Type,
                                              typename Dune::PDELab::TypeTree::TransformTree<typename GFS::template Child<1>::Type,
                                                                                             Dune::PDELab::gfs_to_remote_lfs>::Type
#endif
                                              > TreeNode;

      template<typename>
      friend struct PropagateGlobalStorageVisitor;

      template<typename>
      friend struct ClearSizeVisitor;

      template<typename>
      friend struct ComputeSizeVisitor;

      template<typename>
      friend struct FillIndicesVisitor;

      template<typename ChildGFS>
      static
      Dune::shared_ptr<typename Dune::PDELab::TypeTree::TransformTree<ChildGFS,Dune::PDELab::gfs_to_remote_lfs>::Type>
      createRemoteLocalFunctionSpace(Dune::shared_ptr<const ChildGFS> cgfs)
      {
        Dune::PDELab::gfs_to_remote_lfs trafo;
        return Dune::PDELab::TypeTree::TransformTree<ChildGFS,Dune::PDELab::gfs_to_remote_lfs>::transform_storage(cgfs, trafo);
      }

    public:
      typedef GridGlueLocalFunctionSpaceBaseTraits<GFS,DOFIndex> Traits;

      typedef PowerLocalFunctionSpaceTag ImplementationTag;

      //! \brief initialize with grid function space
      template<typename Transformation>
      GridGlueLocalFunctionSpaceNode (shared_ptr<const GFS> gfs,
                                      const Transformation& t,
                                      Dune::shared_ptr<LFS0> lfs0,
                                      Dune::shared_ptr<LFS1> lfs1)
        : TreeNode(lfs0,lfs1
#ifdef GRIDGLUELFSMIXIN
          ,
                   createRemoteLocalFunctionSpace(
                     gfs->template childStorage<0>()),
                   createRemoteLocalFunctionSpace(
                     gfs->template childStorage<1>())
#endif
          )
        , BaseT(gfs)
      {}

      template<typename Transformation>
      GridGlueLocalFunctionSpaceNode (const GFS& gfs,
                                      const Transformation& t,
                                      Dune::shared_ptr<LFS0> lfs0,
                                      Dune::shared_ptr<LFS1> lfs1)
        : TreeNode(lfs0,lfs1
#ifdef GRIDGLUELFSMIXIN
          ,
                   createRemoteLocalFunctionSpace(
                     gfs.template childStorage<0>()),
                   createRemoteLocalFunctionSpace(
                     gfs.template childStorage<1>())
#endif
          )
        , BaseT(stackobject_to_shared_ptr(gfs))
      {}

    public:
      //! \brief bind local function space to one of the GridGlue contextes (sub-domain cell or remote intersection)
      template<typename Context>
      void bind(const Context & c) const
      {
        this->bind(*this, this->template child<0>(), c);
        this->bind(*this, this->template child<1>(), c);
#ifdef GRIDGLUELFSMIXIN
        // These are the RemoteLocalFunctionSpace sub-trees
        this->bind(*this, this->template child<2>(), c);
        this->bind(*this, this->template child<3>(), c);
#endif
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
#if 0
        assert(&node == this);
        // compute sizes
        RemoteLFSComputeSizeVisitor<CouplingIntersection> csv(rit);
        TypeTree::applyToTree(node,csv);
#endif
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
    template<typename SourceNode, typename Transformation>
    struct gridglue_gfs_to_lfs_template
    {
      template<typename TC0, typename TC1, typename... DUMMY>
      struct result
      {
        typedef GridGlueLocalFunctionSpaceNode<SourceNode,typename Transformation::DOFIndex,TC0,TC1> type;
      };
    };
    template<typename GridGlueGridFunctionSpace, typename Params>
    Dune::PDELab::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      GridGlueGridFunctionSpace,
      gfs_to_lfs<Params>,
      gridglue_gfs_to_lfs_template< GridGlueGridFunctionSpace,gfs_to_lfs<Params> >::template result
      >
    registerNodeTransformation(GridGlueGridFunctionSpace* cgfs, gfs_to_lfs<Params>* t, GridGlueGridFunctionSpaceTag* tag);

  } // end PDELab
} // end Dune

#endif // DUNE_PDELAB_GRIDGLUELFS_HH
