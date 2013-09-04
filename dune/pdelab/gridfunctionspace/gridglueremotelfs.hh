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
      friend struct RemoteLFSComputeSizeVisitor;
    };

    //! traits for single component local function space
    template<typename LFS /*, typename GG */>
    struct RemoteLeafLocalFunctionSpaceTraits :
      public LFS::Traits
    {
      typedef typename LFS::Traits::GridFunctionSpace RemoteGridFunctionSpace;
      typedef typename LFS::Traits::FiniteElementType RemoteFiniteElementType;
      typedef FiniteElementInterfaceSwitch<
        RemoteFiniteElementType> RemoteFESwitch;
      typedef BasisInterfaceSwitch<
        typename RemoteFESwitch::Basis> RemoteBasisSwitch;
      // TODO extract DF from GridGlue
      typedef typename RemoteBasisSwitch::DomainField DomainField;
      typedef typename RemoteBasisSwitch::RangeField RangeField;
      // TODO extract dimension from GridGlue
      enum { dim = RemoteGridFunctionSpace::Traits::GridView::dimension-1 };
      enum { order = 1 };
      enum { diffOrder = 1 };

      //! Type of local finite element
      typedef MonomLocalFiniteElement<DomainField, RangeField, dim, order, diffOrder> FiniteElementType;
      typedef FiniteElementType FiniteElement;

      //! \brief Type of constraints engine
      typedef NoConstraints ConstraintsType;
      typedef ConstraintsType Constraints;
    };

    template<typename LFS /*, typename GG */>
    struct RemoteLeafLocalFunctionSpace
      : public LFS
      // : public LocalFunctionSpaceBaseNode<GG,DOFIndex>
      // , public TypeTree::LeafNode
    {
    public:
      template<typename>
      friend struct RemoteLFSComputeSizeVisitor;

      typedef RemoteLeafLocalFunctionSpaceTraits<LFS /*, GG*/> Traits;
      typedef FiniteElementInterfaceSwitch<
        typename Traits::FiniteElementType
        > FESwitch;

      template <typename... T>
      RemoteLeafLocalFunctionSpace(T&&... t)
        : LFS(std::forward<T>(t)...)
        , gt(GeometryType::simplex, Traits::dim)
        , fem(gt)
      {}

      //! get finite element
      const typename Traits::FiniteElementType& finiteElement () const
      {
        return fem;
      }

      //! \brief get local finite element
      const typename Traits::FiniteElementType& localFiniteElement () const
        DUNE_DEPRECATED
      {
        return fem;
      }

      //! \brief get constraints engine
      const typename Traits::ConstraintsType& constraints () const
      {
        // TODO return empty CC
//        return this->pgfs->constraints();
      }

      //! Calculates the multiindices associated with the given entity.
      template<typename Entity, typename DOFIndexIterator>
      void dofIndices(const Entity& e, DOFIndexIterator it, DOFIndexIterator endit)
      {
        // get layout of entity
        const typename FESwitch::Coefficients &coeffs =
          FESwitch::coefficients(fem);

        // evaluate consecutive index of remote intersection
        // typename GG::IndexType index = e.index();
        std::size_t index = e.index();

        for (std::size_t i = 0; i < std::size_t(coeffs.size()); ++i, ++it)
        {
          // store data
          typedef typename Traits::GridFunctionSpace GFS;
          GFS::Ordering::Traits::DOFIndexAccessor::store(*it,gt,index,i);

          // make sure we don't write past the end of the iterator range
          assert(it != endit);
        }
      }

      template<typename GC, typename LC>
      void insert_constraints (const LC& lc, GC& gc) const
      {
        // // LC and GC are maps of maps
        // typedef typename LC::const_iterator local_col_iterator;
        // typedef typename LC::value_type::second_type::const_iterator local_row_iterator;
        // typedef typename GC::iterator global_col_iterator;
        // typedef typename GC::value_type::second_type global_row_type;

        // for (local_col_iterator cit=lc.begin(); cit!=lc.end(); ++cit)
        // {

        //   // look up entry in global map, if not found, insert an empty one.
        //   global_col_iterator gcit = gc.insert(std::make_pair(std::ref(this->dofIndex(cit->first)),global_row_type())).first;

        //   // copy row to global container with transformed indices
        //   for (local_row_iterator rit=(cit->second).begin(); rit!=(cit->second).end(); ++rit)
        //     gcit->second[this->dofIndex(rit->first)] = rit->second;
        // }
      }

    protected:
      GeometryType gt;
      typename Traits::FiniteElementType fem;
    };

    // Register LeafGFS -> RemoteLFS transformation
    template<typename GridFunctionSpace, typename Params>
    Dune::TypeTree::GenericLeafNodeTransformation<
      GridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      RemoteLeafLocalFunctionSpace<typename TypeTree::TransformTree<GridFunctionSpace, gfs_to_lfs<Params> >::Type>
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
    Dune::TypeTree::TemplatizedGenericPowerNodeTransformation<
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
    Dune::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
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
    Dune::TypeTree::TemplatizedGenericVariadicCompositeNodeTransformation<
      GridGlueGridFunctionSpace,
      gfs_to_remote_lfs<Params>,
      gridglue_gfs_to_remote_lfs_template<GridGlueGridFunctionSpace,gfs_to_remote_lfs<Params> >::template result
      >
    registerNodeTransformation(GridGlueGridFunctionSpace* cgfs, gfs_to_remote_lfs<Params>* t, GridGlueGridFunctionSpaceTag* tag);

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_GRIDGLUEGRIDREMOTELFS_HH
