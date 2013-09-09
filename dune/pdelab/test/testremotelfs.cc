// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/p1fem.hh>
#include <dune/pdelab/finiteelementmap/q12dfem.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/q22dfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridglueremotelfs.hh>
#include <dune/pdelab/gridfunctionspace/gridgluegridfunctionspace.hh>

// test function trees
template<int dim>
struct test;

template<>
struct test<2> {

  template<typename Params, class GFS>
  static void _testremotelfs(GFS & gfs);

  template<class GFS>
  static void testremotelfs(GFS & gfs);

  template<class GV>
  static void testleafgridfunction(const GV& gv)
  {
    // instantiate finite element maps
    Dune::GeometryType gt;
    gt.makeCube(2);
    typedef Dune::PDELab::P0LocalFiniteElementMap<float,double,GV::dimension> P0FEM;
    P0FEM p0fem(gt);
    typedef Dune::PDELab::Q12DLocalFiniteElementMap<float,double> Q12DFEM;
    Q12DFEM q12dfem;
    typedef Dune::PDELab::Q22DLocalFiniteElementMap<float,double> Q22DFEM;
    Q22DFEM q22dfem;

    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM> P0GFS;
    P0GFS p0gfs(gv,p0fem);
    testremotelfs(p0gfs);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q12DFEM> GFS1;
    GFS1 gfs1(gv,q12dfem);
    testremotelfs(gfs1);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM> GFS2;
    GFS2 gfs2(gv,q22dfem);
    testremotelfs(gfs2);
    typedef Dune::PDELab::GridFunctionSpace<GV,Q22DFEM,Dune::PDELab::NoConstraints,
      Dune::PDELab::ISTLVectorBackend<> > GFS3;
    GFS3 gfs3(gv,q22dfem);
    testremotelfs(gfs3);

    // test power
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,2,Dune::PDELab::ISTLVectorBackend<> > PGFS2;
    PGFS2 pgfs2(gfs2,gfs2);
    testremotelfs(pgfs2);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,3,Dune::PDELab::ISTLVectorBackend<> > PGFS3;
    PGFS3 pgfs3(gfs2,gfs2,gfs2);
    testremotelfs(pgfs3);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,4,Dune::PDELab::ISTLVectorBackend<> > PGFS4;
    PGFS4 pgfs4(gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs4);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,5,Dune::PDELab::ISTLVectorBackend<> > PGFS5;
    PGFS5 pgfs5(gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs5);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,6,Dune::PDELab::ISTLVectorBackend<> > PGFS6;
    PGFS6 pgfs6(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs6);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,7,Dune::PDELab::ISTLVectorBackend<> > PGFS7;
    PGFS7 pgfs7(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs7);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,8,Dune::PDELab::ISTLVectorBackend<> > PGFS8;
    PGFS8 pgfs8(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs8);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,9,Dune::PDELab::ISTLVectorBackend<> > PGFS9;
    PGFS9 pgfs9(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs9);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,10,Dune::PDELab::ISTLVectorBackend<> > PGFS10;
    PGFS10 pgfs10(gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2,gfs2);
    testremotelfs(pgfs10);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::ISTLVectorBackend<> > PGFS17;
    PGFS17 pgfs17(gfs2);
    testremotelfs(pgfs17);
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS2,17,Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab:: EntityBlockedOrderingTag> PGFS17B;
    PGFS17B pgfs17b(gfs2);
    testremotelfs(pgfs17b);

    // test composite
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2> CGFS2;
    CGFS2 cgfs2(gfs1,pgfs2);
    testremotelfs(cgfs2);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2> CGFS3;
    CGFS3 cgfs3(gfs1,pgfs2,cgfs2);
    testremotelfs(cgfs3);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3> CGFS4;
    CGFS4 cgfs4(gfs1,pgfs2,cgfs2,cgfs3);
    testremotelfs(cgfs4);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4> CGFS5;
    CGFS5 cgfs5(gfs1,pgfs2,cgfs2,cgfs3,cgfs4);
    testremotelfs(cgfs5);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5> CGFS6;
    CGFS6 cgfs6(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5);
    testremotelfs(cgfs6);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6> CGFS7;
    CGFS7 cgfs7(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6);
    testremotelfs(cgfs7);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7> CGFS8;
    CGFS8 cgfs8(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7);
    testremotelfs(cgfs8);
    typedef Dune::PDELab::CompositeGridFunctionSpace<Dune::PDELab::ISTLVectorBackend<>,
        Dune::PDELab::LexicographicOrderingTag,GFS1,PGFS2,CGFS2,CGFS3,CGFS4,CGFS5,CGFS6,CGFS7,CGFS8> CGFS9;
    CGFS9 cgfs9(gfs1,pgfs2,cgfs2,cgfs3,cgfs4,cgfs5,cgfs6,cgfs7,cgfs8);
    testremotelfs(cgfs9);
  }
};

template<typename Params, class GFS>
void test<2>::_testremotelfs(GFS & gfs)
{
  typedef typename Dune::TypeTree::TransformTree<GFS,Dune::PDELab::gfs_to_remote_lfs<Params> > Trafo;
  Dune::PDELab::gfs_to_remote_lfs<Params> trafo;
  typedef typename Trafo::Type RemoteLFS;
  RemoteLFS rlfs = Trafo::transform(gfs, trafo);
  Trafo::transform_storage(Dune::stackobject_to_shared_ptr(gfs), trafo);
}

template<class GFS>
void test<2>::testremotelfs(GFS & gfs)
{
  _testremotelfs<GFS>(gfs);
  _testremotelfs<Dune::Empty>(gfs);
}

template<class GV>
void testleafgridfunction(const GV& gv)
{
  test<GV::dimension>::testleafgridfunction(gv);
}


int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // 2D
    {
      std::cout << "2D tests" << std::endl;
      // need a grid in order to test grid functions
      Dune::FieldVector<double,2> L(1.0);
      Dune::FieldVector<int,2> N(1);
      Dune::FieldVector<bool,2> B(false);
      Dune::YaspGrid<2> grid(L,N,B,0);
      grid.globalRefine(1);

      testleafgridfunction(grid.leafView());
    }

    // test passed
    return 0;

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
