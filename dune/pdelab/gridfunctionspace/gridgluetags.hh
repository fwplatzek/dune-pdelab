#ifndef DUNE_PDELAB_GRIDGLUETAGS_HH
#define DUNE_PDELAB_GRIDGLUETAGS_HH

#include "tags.hh"

#define GRIDGLUELFSMIXIN

namespace Dune
{
  namespace PDELab
  {

    struct NoObject {};

    template<>
    struct gfs_to_lfs<NoObject> {
      typedef NoObject Params;
      //! The MultiIndex type that will be used in the resulting LocalFunctionSpace tree.
      //typedef Dune::PDELab::MultiIndex<std::size_t,TypeTree::TreeInfo<GFS>::depth> MultiIndex;
      typedef NoObject DOFIndex;
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

    template<typename GFS>
    struct gfs_to_remote_gfs {
      typedef GFS GridFunctionSpace;
    };

    struct GridGlueGridFunctionSpaceTag : public CompositeGridFunctionSpaceTag {};

  } // end PDELab
} // end Dune

#endif // DUNE_PDELAB_GRIDGLUETAGS_HH
