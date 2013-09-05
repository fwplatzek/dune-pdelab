// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#if ! HAVE_PSURFACE
#error we need psurface
#endif

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/classname.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/test/couplingtest.hh>

#include "../finiteelementmap/q1fem.hh"
#include "../gridfunctionspace/gridfunctionspace.hh"
#include "../gridfunctionspace/gridfunctionspaceutilities.hh"
#include "../gridfunctionspace/interpolate.hh"
#include "../constraints/conforming.hh"
#include "../common/function.hh"
#include "../common/vtkexport.hh"
#include "../backend/istlvectorbackend.hh"
#include "../backend/istlmatrixbackend.hh"
#include "../gridoperator/gridoperator.hh"
#include "../backend/istlsolverbackend.hh"
#include "../localoperator/laplacedirichletp12d.hh"
#include "../localoperator/poisson.hh"
#include "../gridfunctionspace/vtk.hh"
#include <dune/pdelab/gridfunctionspace/gridgluegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridgluelocalfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include "../gridoperator/gridglueoperator.hh"

using namespace Dune;

template <class GridView>
class VerticalFaceDescriptor
  : public ExtractorPredicate<GridView,1>
{
public:
  VerticalFaceDescriptor(double sliceCoord)
    : sliceCoord_(sliceCoord)
  {}

  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& eptr,
                        unsigned int face) const
  {
    const int dim = GridView::dimension;
    const Dune::GenericReferenceElement<double,dim>& refElement = Dune::GenericReferenceElements<double, dim>::general(eptr->type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(eptr->geometry().corner(refElement.subEntity(face,1,i,dim))[0] - sliceCoord_) > 1e-6 )
        return false;

    return true;
  }

private:
  double sliceCoord_;
};

struct LFSCheck {
  template<Dune::PDELab::GridGlueContextTag TAG, typename LFS, typename CTX>
  static void check(LFS & lfs, const CTX & c, std::string info)
  {
    lfs.bind(Dune::PDELab::GridGlueContext<CTX,TAG>(c));
    assert((TAG == Dune::PDELab::GFS_DOM0   && lfs.template child<0>().size() != 0) || lfs.template child<0>().size() == 0);
    assert((TAG == Dune::PDELab::GFS_DOM1   && lfs.template child<1>().size() != 0) || lfs.template child<1>().size() == 0);
    assert((TAG == Dune::PDELab::TRACE_DOM0 && lfs.template child<2>().size() != 0) || lfs.template child<2>().size() == 0);
    assert((TAG == Dune::PDELab::TRACE_DOM1 && lfs.template child<3>().size() != 0) || lfs.template child<3>().size() == 0);
    std::cout << info << "\t"
              << lfs.template child<0>().size() << "\t"
              << lfs.template child<1>().size() << "\t"
              << lfs.template child<2>().size() << "\t"
              << lfs.template child<3>().size() << "\n";
    // if (TAG == Dune::PDELab::TRACE_DOM0 || TAG == Dune::PDELab::TRACE_DOM1)
      for (unsigned int i=0; i<lfs.size(); i++)
        std::cout << lfs.dofIndex(i) << "\n";
  }
};

template<Dune::PDELab::GridGlueContextTag MAINTAG>
struct LFSSubCheck {
  template<Dune::PDELab::GridGlueContextTag TAG, typename LFS, typename CTX>
  static typename enable_if<TAG!=MAINTAG>::type
  check(LFS & lfs, const CTX & c, std::string info)
  {
    std::cout << info << "\t"
              << "skipped" << "\n";
  }
  template<Dune::PDELab::GridGlueContextTag TAG, typename LFS, typename CTX>
  static typename enable_if<TAG==MAINTAG>::type
  check(LFS & lfs, const CTX & c, std::string info)
  {
    lfs.bind(c);
    std::cout << info << "\t"
              << lfs.size() << "\n";
    for (unsigned int i=0; i<lfs.size(); i++)
      std::cout << lfs.dofIndex(i) << "\n";
  }
};

template<typename ContextOperator, typename LFS, typename GG>
void testlfs(LFS & lfs, const GG & gg)
{
  typedef typename GG::Grid0Patch Patch0;
  typedef typename GG::Grid1Patch Patch1;

  {
    typedef typename Patch0::GridView GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIt;
    for (ElementIt eit = gg.template gridView<0>().template begin<0>();
         eit != gg.template gridView<0>().template end<0>(); ++eit)
    {
      ContextOperator::template check<Dune::PDELab::GFS_DOM0>(lfs,*eit,"GFS0");
    }
  }

  {
    typedef typename Patch1::GridView GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIt;
    for (ElementIt eit = gg.template gridView<1>().template begin<0>();
         eit != gg.template gridView<1>().template end<0>(); ++eit)
    {
      ContextOperator::template check<Dune::PDELab::GFS_DOM1>(lfs,*eit,"GFS1");
    }
  }

  typedef GridGlueView<Patch0,Patch1,0> GGView0; // from 0 to 1
  typedef GridGlueView<Patch0,Patch1,1> GGView1; // from 1 to 0
  {
    typedef typename GGView0::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    for (IntersectionIterator iit = gg.template ibegin<0>();
         iit != gg.template iend<0>(); ++iit)
    {
      ContextOperator::template check<Dune::PDELab::GFS_DOM0>(lfs,*iit->inside(),"GFS0");
      ContextOperator::template check<Dune::PDELab::TRACE_DOM0>(lfs,*iit,"TRACE0");
      if (iit->neighbor())
        ContextOperator::template check<Dune::PDELab::GFS_DOM1>(lfs,*iit->outside(),"GFS1");
    }
  }

  {
    typedef typename GGView1::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    for (IntersectionIterator iit = gg.template ibegin<1>();
         iit != gg.template iend<1>(); ++iit)
    {
      ContextOperator::template check<Dune::PDELab::GFS_DOM1>(lfs,*iit->inside(),"GFS1");
      ContextOperator::template check<Dune::PDELab::TRACE_DOM1>(lfs,*iit,"TRACE1");
      if (iit->neighbor())
        ContextOperator::template check<Dune::PDELab::GFS_DOM0>(lfs,*iit->outside(),"GFS0");
    }
  }

}

namespace Dune {
  namespace PDELab {

    template<typename A0, typename A1, typename C0, typename C1, typename D0, typename D1>
    class GridGlueLocalProblem :
      public NumericalJacobianApplyVolume< GridGlueLocalProblem<A0,A1,C0,C1,D0,D1> >,
      public NumericalJacobianVolume< GridGlueLocalProblem<A0,A1,C0,C1,D0,D1> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doLambdaVolume = true };
      enum { doLambdaBoundary = true };

      GridGlueLocalProblem ()
      {}

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
      }

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
      }

      // boundary integral independen of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
      {
      }
    };
  } // end namespace Dune
} // end namespace Dune

// boundary grid function selecting boundary conditions
class ConstraintsParameters
  : public Dune::PDELab::DirichletConstraintsParameters
{

public:

  template<typename I>
  bool isDirichlet(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return true;
  }

  template<typename I>
  bool isNeumann(const I & ig, const Dune::FieldVector<typename I::ctype, I::dimension-1> & x) const
  {
    return false;
  }

};

template <int dim>
void testNonMatchingCubeGrids()
{

  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(2);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType cubeGrid0(elements, lower, upper);

  elements = 3;
  lower[0] += 1;
  upper[0] += 1;

  GridType cubeGrid1(elements, lower, upper);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType::LevelGridView DomGridView;
  typedef typename GridType::LevelGridView TarGridView;

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

  // ///////////////////////////////////////////
  //   Setup sub domain GFS
  // ///////////////////////////////////////////

  // constants and types
  typedef typename GlueType::ctype DF;
  typedef double RF;

  // make finite element map
  typedef Dune::PDELab::Q1LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem;

  // backends
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTLVectorBackend<> Backend;
  typedef Dune::PDELab::ISTLVectorBackend<Dune::PDELab::ISTLParameters::dynamic_blocking> BlockingBackend;

  // Dom GFS
  typedef Dune::PDELab::GridFunctionSpace<DomGridView,FEM,CON,Backend> GFSDOM;
  GFSDOM gfsdom(glue.template gridView<0>(),fem);
  gfsdom.name("dom");

  // Tar GFS
  typedef Dune::PDELab::GridFunctionSpace<TarGridView,FEM,CON,Backend> GFSTAR;
  GFSTAR gfstar(glue.template gridView<1>(),fem);
  gfstar.name("tar");

  typedef typename Dune::TypeTree::TransformTree<GFSDOM,
                                                 Dune::PDELab::gfs_to_remote_gfs<GFSDOM> > RemoteGFSTrafo;
  Dune::PDELab::gfs_to_remote_gfs<GFSDOM> trafo;
  typedef typename RemoteGFSTrafo::Type RemoteGFSDOM;
  RemoteGFSDOM remotegfsdom = RemoteGFSTrafo::transform(gfsdom, trafo);

  // GridGlue GFS
  typedef Dune::PDELab::GridGlueGridFunctionSpace<GlueType,GFSDOM,GFSTAR,BlockingBackend> GlueGFS;
  GlueGFS gluegfs(glue,gfsdom,gfstar);

  typedef Dune::PDELab::LocalFunctionSpace<GlueGFS> GlueLFS;
  GlueLFS gluelfs(gluegfs);

  // try creating remoteLFS
  {
    typedef typename Dune::TypeTree::TransformTree<GFSTAR,Dune::PDELab::gfs_to_remote_lfs<GlueGFS> > Trafo;
    Dune::PDELab::gfs_to_remote_lfs<GlueGFS> trafo;
    typedef typename Trafo::Type RemoteLFS;
    RemoteLFS remotelfs = Trafo::transform(gfstar, trafo);
  }
  {
    typedef typename Dune::TypeTree::TransformTree<GlueGFS,Dune::PDELab::gfs_to_remote_lfs<GlueGFS> > Trafo;
    Dune::PDELab::gfs_to_remote_lfs<GlueGFS> trafo;
    typedef typename Trafo::Type RemoteLFS;
    RemoteLFS remotelfs = Trafo::transform(gluegfs, trafo);
  }
  {
    typedef typename Dune::TypeTree::TransformTree<GlueGFS,Dune::PDELab::gfs_to_remote_lfs<Dune::PDELab::NoObject> > Trafo;
    Dune::PDELab::gfs_to_remote_lfs<Dune::PDELab::NoObject> trafo;
    typedef typename Trafo::Type RemoteLFS;
    RemoteLFS remotelfs = Trafo::transform(gluegfs, trafo);
  }

  // try getting subspaces
  typedef Dune::TypeTree::TreePath<0> Path0;
  typedef Dune::TypeTree::TreePath<1> Path1;
  typedef Dune::PDELab::GridFunctionSubSpace<GlueGFS, Path0> GlueGFS0;
  typedef Dune::PDELab::GridFunctionSubSpace<GlueGFS, Path1> GlueGFS1;
  typedef Dune::PDELab::LocalFunctionSpace<GlueGFS0> GlueLFS0;
  typedef Dune::PDELab::LocalFunctionSpace<GlueGFS1> GlueLFS1;
  GlueGFS0 gluegfs0(gluegfs);
  GlueGFS1 gluegfs1(gluegfs);
  GlueLFS0 gluelfs0(gluegfs0);
  GlueLFS1 gluelfs1(gluegfs1);

  // test the features of the local function space
  std::cout << "====================== GlueLFS =================\n";
  testlfs<LFSCheck>(gluelfs, glue);
  std::cout << "====================== done ====================\n";
  std::cout << "====================== GlueLFS0 ================\n";
  testlfs<LFSSubCheck<Dune::PDELab::GFS_DOM0> >(gluelfs0, glue);
  std::cout << "====================== done ====================\n";
  std::cout << "====================== GlueLFS1 ================\n";
  testlfs<LFSSubCheck<Dune::PDELab::GFS_DOM1> >(gluelfs1, glue);
  std::cout << "====================== done ====================\n";

  Dune::PDELab::EmptyTransformation constraints;
  Dune::PDELab::GridGlueAssembler<GlueGFS,GlueGFS,Dune::PDELab::EmptyTransformation,Dune::PDELab::EmptyTransformation>
    assembler(gluegfs,gluegfs,constraints,constraints);

  // make constraints map and initialize it from a function
  typedef double RF;
  typedef typename GlueGFS::template ConstraintsContainer<RF>::Type C;
  C cc;
  cc.clear();
  ConstraintsParameters constraintsparameters;
  Dune::PDELab::constraints(constraintsparameters,gluegfs0,cc);
  Dune::PDELab::constraints(constraintsparameters,gluegfs1,cc);

  // make grid operator
  typedef Dune::PDELab::GridGlueLocalProblem<int,int,int,int,int,int> LOP;
  LOP lop;
  typedef Dune::PDELab::GridGlueOperator<GlueGFS,GlueGFS,LOP,
                                     Dune::PDELab::ISTLMatrixBackend,
                                     double,double,double,
                                     C,C> GridOperator;
  GridOperator gridoperator(gluegfs,cc,gluegfs,cc,lop);

  typedef typename GridOperator::Traits::Jacobian M;
  std::cout << "Allocate Siffness Matrix " << Dune::className<typename M::BaseT>() << std::endl;
  M mat(gridoperator, 0.0);
  Dune::printmatrix(std::cout, mat.base(), "full_matrix", "r");

  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gluegfs, 0.0);
  Dune::printvector(std::cout, x0.base(), "vec", "r");
  gridoperator.jacobian(x0,mat);

//  Dune::printmatrix(std::cout, mat.base(), "full_matrix", "r");
}

int main(int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  std::cout << "============================================================\n";
  testNonMatchingCubeGrids<2>();
  std::cout << "============================================================\n";

  return 0;
}
catch (Exception e) {
  int i = 0; char** c = 0;
  std::cout << Dune::MPIHelper::instance(i,c).rank() << ": " << e << std::endl;
  return 1;
}
