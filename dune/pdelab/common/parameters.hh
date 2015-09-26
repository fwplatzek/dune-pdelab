#ifndef DUNE_PDELAB_COMMON_PARAMETERS_HH
#define DUNE_PDELAB_COMMON_PARAMETERS_HH

/**
   \file parameters.hh Add support for functions style definition of LocalOperator parameter classes

   \todo implement bind and unbind for elements
   \todo implement bind and unbind for intersection
   \todo implement setTime
   \todo share storage for time, entity, intersection, where possible
 */

#include <utility>
#include <type_traits>

namespace Dune {
namespace PDELab {
    namespace Parameters {
        namespace Imp {

            struct NoParameter {};

#undef DefinePDELabParameterName
#define DefinePDELabParameterName(Name,Variable)                        \
            template<typename T, typename Base = ::Dune::PDELab::Parameters::Imp::NoParameter> \
            struct Name##Parameter : public Base                        \
            {                                                           \
                template <typename B>                                   \
                    using bindBaseClass = Name##Parameter<T,B>;         \
                Name##Parameter (const Name##Parameter &) = default;    \
                Name##Parameter (T && t) :                              \
                    Variable(std::forward<T>(t)) {}                     \
                Name##Parameter (const Base & b, const Name##Parameter<T,::Dune::PDELab::Parameters::Imp::NoParameter> & other) : \
                    Base(b),                                            \
                    Variable(other.Variable) {}                         \
                typename std::decay<T>::type Variable;                  \
            };                                                          \
            template<typename T>                                        \
            Name##Parameter<T> define##Name(T && t)                     \
            {                                                           \
                return { std::forward<T>(t) };                          \
            }

            template<typename Head>
            Head
            merge(const Head & head)
            {
                return head;
            };

            template<typename Head, typename... Pack>
            auto
            merge(const Head & head, Pack && ... pack)
                -> typename Head::template bindBaseClass<decltype(merge(pack...))>
            {
                auto remainder = merge(pack...);
                return typename Head::template bindBaseClass<decltype(remainder)>(remainder,head);
            };

        }

        /**
           \brief Combine individual parameters in to one parameter class

           \code{.cpp}
           auto param =
               Parameters::merge(
                   ConvectionDiffusionParameters::definePermeability(K),
                   ConvectionDiffusionParameters::defineBoundaryCondition(bc)
               );
           \endcode
         */
        template<typename... Pack>
        // Imp::MergedParameters<Pack...>
        auto
        merge(Pack && ... pack)
            -> decltype(Imp::merge(pack...))
        {
            return Imp::merge(pack...);
        };

        /**
           \brief Overwrite te implementation of a particular parameter a given parameter class

           \code{.cpp}
           auto param_new =
               Parameters::overwrite(
                   param,
                   ConvectionDiffusionParameters::defineBoundaryCondition(bc)
               );
           \endcode
         */
        template<typename Base, typename P>
        typename P::template bindBaseClass<Base>
        overwrite(const Base & b, const P & p)
        {
            return typename P::template bindBaseClass<Base>(b,p);
        };
    }

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_PARAMETERS_HH
