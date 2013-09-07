// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_REMOTEORDERING_HH
#define DUNE_PDELAB_ORDERING_REMOTEORDERING_HH

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! Gridview ordering for leaf spaces
    template<typename GFS, typename Transformation>
    class RemoteOrdering
      : public TypeTree::LeafNode
      , public VirtualOrderingBase<typename Transformation::DOFIndex,
                                   typename Transformation::ContainerIndex>
      , public OrderingBase<typename Transformation::DOFIndex,
                            typename Transformation::ContainerIndex>
    {
    public:
      typedef OrderingTraits<typename Transformation::DOFIndex,
                             typename Transformation::ContainerIndex> Traits;

      static const bool has_dynamic_ordering_children = false;

      static const bool consume_tree_index = false;

    protected:

      typedef TypeTree::LeafNode NodeT;

      typedef OrderingBase<typename Transformation::DOFIndex,
                           typename Transformation::ContainerIndex> BaseT;

    public:

      RemoteOrdering(const GFS& gfs, const Transformation & t)
        : BaseT(*this,true,const_cast<GFS*>(&gfs),this)
      {}

      RemoteOrdering(shared_ptr<const GFS>& gfs, const Transformation & t)
        : BaseT(*this,true,const_cast<GFS*>(gfs.get()),this)
      {}

#ifndef DOXYGEN

// we need to override the default copy / move ctor to fix the delegate pointer, but that is
// hardly interesting to our users...

      RemoteOrdering(const RemoteOrdering& r)
        : BaseT(r)
      {
        this->setDelegate(this);
      }

#if HAVE_RVALUE_REFERENCES

      RemoteOrdering(RemoteOrdering&& r)
        : BaseT(std::move(r))
      {
        this->setDelegate(this);
      }

#endif // HAVE_RVALUE_REFERENCES

#endif // DOXYGEN

      virtual void map_index_dynamic(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {
        mapIndex(di,ci);
      }

      typename Traits::ContainerIndex mapIndex(const typename Traits::DOFIndex& di) const
      {
        typename Traits::ContainerIndex ci;
        mapIndex(di.view(),ci);
        return ci;
      }

      void mapIndex(typename Traits::DOFIndexView di, typename Traits::ContainerIndex& ci) const
      {

        const typename Traits::SizeType geometry_type_index = Traits::DOFIndexAccessor::geometryType(di);
        const typename Traits::SizeType entity_index = Traits::DOFIndexAccessor::entityIndex(di);
        assert (di.treeIndex().size() == 1);
        ci.push_back(di.treeIndex().back());
        ci.push_back(entity_index);
      }

#if 0
      template<typename ItIn, typename ItOut>
      void map_lfs_indices(const ItIn begin, const ItIn end, ItOut out) const
      {
        typedef typename Traits::SizeType size_type;

        for (ItIn in = begin; in != end; ++in, ++out)
        {
          const size_type entity_index = Traits::DOFIndexAccessor::entityIndex(*in);
          assert(in->treeIndex().size() == 1);
          out->push_back(in->treeIndex().back());
          out->push_back(entity_index);
        }
      }

      template<typename CIOutIterator>
      typename Traits::SizeType
      extract_entity_indices(const typename Traits::DOFIndex::EntityIndex& ei,
                             typename Traits::SizeType child_index,
                             CIOutIterator ci_out, const CIOutIterator ci_end) const
      {
        typedef typename Traits::SizeType size_type;

        const size_type entity_index = Traits::DOFIndexAccessor::GeometryIndex::entityIndex(ei);

        for (size_type i = 0; i < size; ++i, ++ci_out)
        {
          ci_out->push_back(i);
          ci_out->push_back(entity_index);
        }
        return 0;
      }
#endif

      using BaseT::fixedSize;
    };


   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFORDERING_HH
