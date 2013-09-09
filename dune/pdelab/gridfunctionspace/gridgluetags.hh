#ifndef DUNE_PDELAB_GRIDGLUETAGS_HH
#define DUNE_PDELAB_GRIDGLUETAGS_HH

#include <dune/common/typetraits.hh>
#include "tags.hh"

#define GRIDGLUEGFSMIXIN
#define GRIDGLUELFSMIXIN
// #define STORE_GRIDGLUE

namespace Dune
{
  namespace PDELab
  {

    template<>
    struct gfs_to_lfs<Empty> {
      typedef Empty Params;
      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef Empty DOFIndex;
    };

    // template<typename GG>
    // struct gfs_to_remote_lfs {
    //   typedef GG GridGlue;
    //   GridGlue _gridGlue;
    // };
    template<typename GFS>
    struct gfs_to_remote_lfs {
      typedef GFS GridFunctionSpace;
    };

    template<typename GG>
    struct gfs_to_remote_gfs {
      // typedef GG GridGlue;
      // gfs_to_remote_gfs(const GG & gg) : gridglue(gg) {}
      // const GridGlue & gridglue;
    };

    struct RemoteLeafFunctionSpaceTag : public LeafGridFunctionSpaceTag {};
    struct GridGlueGridFunctionSpaceTag : public CompositeGridFunctionSpaceTag {};

  } // end PDELab
} // end Dune

#endif // DUNE_PDELAB_GRIDGLUETAGS_HH
