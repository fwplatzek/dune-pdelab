#include <config.h>

#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/common/parameters.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/function/helpers.hh>

template <typename GV, typename Param>
void assemble(const GV & gv, Param & p)
{
    for (const auto & e : elements(gv))
    {
        Dune::FieldVector<double,2> x(0.5);
        p.bind(e);
        std::cout << x << "\t"
                  << p.f(x) << "\t"
                  << p.bctype(e,x) << std::endl;
    }
    std::cout << "---------------------------" << std::endl;
}

template<typename... Args>
Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type
freeFunctionBCType (const Args & ... )
{
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
};

int main(int argc, char** argv)
    try{
        //Maybe initialize Mpi
        Dune::MPIHelper::instance(argc, argv);

        // need a grid in order to test grid functions
        Dune::FieldVector<double,2> L(1.0);
        Dune::array<int,2> N(Dune::fill_array<int,2>(1));
        Dune::YaspGrid<2> grid(L,N);
        grid.globalRefine(1);

        // setup parameters
        using GV = Dune::YaspGrid<2>::LeafGridView;
        using RF = float;
        using BasicParam = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;

        using Element = typename BasicParam::Traits::ElementType;
        using Domain = typename BasicParam::Traits::DomainType;

        auto f = [](Element e, Domain x) { return e.geometry().global(x).two_norm(); };

        using namespace Dune::TypeTree::Indices;
        using Dune::PDELab::replaceEntityByBind;

        namespace Param = Dune::PDELab::Parameters;
        namespace ModelParam = Dune::PDELab::ConvectionDiffusionParameters;
        auto param =
            Param::merge(
                ModelParam::defineSourceTerm(
                    replaceEntityByBind(f,_1)
                    ),
                ModelParam::defineBoundaryCondition(freeFunctionBCType<Element,Domain>)
                );
        assemble(grid.leafGridView(),param);

        // refine / modify existing parameter class
        auto param2 =
            Param::overwrite(param, // BasicParam(),
                ModelParam::defineSourceTerm(
                    replaceEntityByBind(
                        [](Domain,const Element &) { return 1.0; } ,
                        _2)
                    )
                );
        assemble(grid.leafGridView(),param2);
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
