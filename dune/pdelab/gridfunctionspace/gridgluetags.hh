#ifndef DUNE_PDELAB_GRIDGLUETAGS_HH
#define DUNE_PDELAB_GRIDGLUETAGS_HH

#include "tags.hh"

namespace Dune
{
  namespace PDELab
  {

    // template<typename GG>
    // struct gfs_to_remote_lfs {
    //   typedef GG GridGlue;
    //   GridGlue _gridGlue;
    // };
    struct gfs_to_remote_lfs {};

    struct GridGlueGridFunctionSpaceTag : public CompositeGridFunctionSpaceTag {};

  } // end PDELab
} // end Dune

#endif // DUNE_PDELAB_GRIDGLUETAGS_HH
