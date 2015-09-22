#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/pdelab/common/parameters.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>

template <typename Param>
void assemble(const Param & p)
{
    for (int i = 0; i <= 10; i++)
    {
        double x = i * 0.1;

        std::cout << x << "\t" << p.A(x) << "\t" << p.bc(x) << std::endl;
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
        grid.globalRefine(6);

        // setup parameters
        using GV = Dune::YaspGrid<2>::LeafGridView;
        using RF = float;
        using BasicParam = Dune::PDELab::ConvectionDiffusionModelProblem<GV,RF>;

        using Element = typename BasicParam::Traits::ElementType;
        using Domain = typename BasicParam::Traits::DomainType;

        auto f = [](Element, Domain x) { return x*x; };

        namespace Param = Dune::PDELab::Parameters;
        namespace ModelParam = Dune::PDELab::ConvectionDiffusionParameters;
        auto param =
            Param::overwrite(
                Param::merge(
                    ModelParam::defineSinkTerm([]() { return 1.0; } ),
                    ModelParam::defineBoundaryCondition(freeFunctionBCType<Element,Domain>)
                    ),
                ModelParam::defineSinkTerm(f)
                );
        // assemble(param);

        // refine / modify existing parameter class
        auto param2 =
            Param::overwrite(BasicParam(),
                ModelParam::defineBoundaryCondition(freeFunctionBCType<Element,Domain>)
                );
        // assemble(param2);

    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
