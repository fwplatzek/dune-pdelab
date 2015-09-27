// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/laplacedirichletp12d.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/function/helpers.hh>

#include"gridexamples.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON>
void poisson (const GV& gv, const FEM& fem, std::string filename, int q)
{
  using Dune::PDELab::Backend::native;

  // constants and types
  typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

  // make function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::istl::VectorBackend<> > GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");

  // make model problem
  namespace Param = Dune::PDELab::Parameters;
  namespace ModelParam = Dune::PDELab::ConvectionDiffusionParameters;
  using namespace Dune::TypeTree::Indices;
  using Dune::PDELab::replaceEntityByBind;

  // grid typedefs
  using BCType = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;
  using BasicParam = Dune::PDELab::ConvectionDiffusionModelProblem<GV,R>;
  using Element = typename BasicParam::Traits::ElementType;
  using Domain = Dune::FieldVector<typename GV::ctype,GV::dimensionworld>;
  using Intersection = typename BasicParam::Traits::IntersectionType;
  using IntersectionDomain = typename BasicParam::Traits::IntersectionDomainType;

  // Define parameter functions f,g,j and \partial\Omega_D/N
  auto f = Dune::Functions::makeGridViewFunction (
    [] (const Domain& x)
    {
      return (x[0]>0.25 && x[0]<0.375 && x[1]>0.25 && x[1]<0.375) ? 50.0 : 0.0;
    },
    gv);
  auto g = Dune::Functions::makeGridViewFunction (
    [] (const Domain& x) { using std::exp; return exp(-x.two_norm2()); },
    gv);
  auto bctype =
    [] (const Intersection& is, const IntersectionDomain& x) -> BCType
    {
      auto xglobal = is.geometry().global(x);

      if (xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
      if (xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6)
      {
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      }
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    };
  auto j = [] (const Intersection& is, const IntersectionDomain& x) -> R
    {
      auto xglobal = is.geometry().global(x);
      if (xglobal[0] > 1.0 - 1E-6 && xglobal[1] > 0.5 + 1E-6) {
        return -5.0;
      } else {
        return 0.0;
      }
    };

  auto problem =
    Param::overwrite(
      BasicParam(),
      ModelParam::defineSourceTerm(localFunction(f)),
      ModelParam::defineBoundaryCondition(bctype),
      ModelParam::defineDirichletBoundaryValue(localFunction(g)),
      ModelParam::defineNeumannBoundaryValue(j)
      );

  // make constraints map and initialize it from a function
  typedef typename GFS::template ConstraintsContainer<R>::Type C;
  C cg;
  Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<decltype(problem)> bcond (gv,problem);
  Dune::PDELab::constraints(bcond,gfs,cg);

  // make local operator
  typedef Dune::PDELab::ConvectionDiffusionFEM<decltype(problem),FEM> LOP;
  LOP lop(problem);

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(27);

  // make grid operator
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
                                     MBE,
                                     double,double,double,
                                     C,C> GridOperator;
  GridOperator gridoperator(gfs,cg,gfs,cg,lop,mbe);

  // make coefficent Vector and initialize it from a function
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Domain DV;
  DV x0(gfs);

  // initialize DOFs from Dirichlet extension
  Dune::PDELab::interpolate(g,gfs,x0);
  Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);

  // represent operator as a matrix
  // There is some weird shuffling around here - please leave it in,
  // it's there to test the copy constructor and assignment operator of the
  // matrix wrapper
  typedef typename GridOperator::Traits::Jacobian M;
  M m(gridoperator);
  gridoperator.jacobian(x0,m);

  // evaluate residual w.r.t initial guess
  typedef typename GridOperator::Traits::Range RV;
  RV r(gfs);
  r = 0.0;
  gridoperator.residual(x0,r);

  // solve the jacobian system
  DV x(gfs,0.0);
  {
    using CG = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
    using Solver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator,CG,DV>;
    CG cg(1000,2);
    Solver solver(gridoperator, cg, x, 1e-10, 1e-99, false);
    solver.apply();
  }

  // output grid function with VTKWriter
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
  vtkwriter.write(filename,Dune::VTK::ascii);
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    // YaspGrid Q1 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q1_2d",2);
    }

    // YaspGrid Q2 2D test
    {
      // make grid
      Dune::FieldVector<double,2> L(1.0);
      Dune::array<int,2> N(Dune::fill_array<int,2>(1));
      Dune::YaspGrid<2> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q2_2d",2);
    }

    // YaspGrid Q1 3D test
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q1_3d",2);
    }

    // YaspGrid Q2 3D test
    {
      // make grid
      Dune::FieldVector<double,3> L(1.0);
      Dune::array<int,3> N(Dune::fill_array<int,3>(1));
      Dune::YaspGrid<3> grid(L,N);
      grid.globalRefine(3);

      // get view
      typedef Dune::YaspGrid<3>::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q2_3d",2);
    }

    // UG Pk 2D test
#if HAVE_UG
    {
      // make grid
      std::shared_ptr<Dune::UGGrid<2> > grid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
      grid->globalRefine(4);

      // get view
      typedef Dune::UGGrid<2>::LeafGridView GV;
      const GV& gv=grid->leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_UG_Pk_2d",q);
    }
#endif

#if HAVE_ALBERTA
    {
      // make grid
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_Alberta_Pk_2d",q);
    }
#endif

#if HAVE_DUNE_ALUGRID
    {
      // make grid
      ALUUnitSquare grid;
      grid.globalRefine(4);

      // get view
      typedef ALUUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_ALU_Pk_2d",q);
    }
#endif

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
