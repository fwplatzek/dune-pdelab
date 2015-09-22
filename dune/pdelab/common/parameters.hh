#ifndef DUNE_PDELAB_COMMON_PARAMETERS_HH
#define DUNE_PDELAB_COMMON_PARAMETERS_HH

#include <utility>
#include <tuple>
#include <type_traits>

namespace Dune {
namespace PDELab {
    namespace Parameters {
        namespace Imp {

            struct NoParameter {};

#undef DefinePDELabParameterName
#define DefinePDELabParameterName(Name,Variable)                        \
            template<typename T, typename Base = Imp::NoParameter>      \
            struct Name##Parameter : public Base                        \
            {                                                           \
                template <typename B>                                   \
                    using bindBaseClass = Name##Parameter<T,B>;         \
                Name##Parameter (const Name##Parameter &) = default;    \
                Name##Parameter (T && t) :                              \
                    Variable(std::forward<T>(t)) {}                     \
                Name##Parameter (const Base & b, const Name##Parameter<T,Imp::NoParameter> & other) : \
                    Base(b),                                            \
                    Variable(other.Variable) {}                         \
                typename std::decay<T>::type Variable;                  \
            };                                                          \
            template<typename T>                                        \
            Name##Parameter<T> define##Name(T && t)                     \
            {                                                           \
                return { std::forward<T>(t) };                          \
            }

            // we use multiple inheritance to create a struct containing all classes in the Parameter Pack
            template <typename... Pack> struct Unroll;

            template <typename Head, typename ... Next>
            struct Unroll<Head, Next...>:
                public std::decay<Head>::type,
                public Unroll<Next...>
            {
                Unroll(Head && h, Next && ... n) :
                    std::decay<Head>::type(std::forward<Head>(h)),
                    Unroll<Next...>(std::forward<Next>(n)...)
                {}
            };

            template <>
            struct Unroll<>
            {};

            template<typename... Pack>
            struct MergedParameters :
                public Unroll<Pack...>
            {
                MergedParameters(Pack && ... pack) :
                    Unroll<Pack...>(std::forward<Pack>(pack)...)
                {}
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
        Imp::MergedParameters<Pack...>
        merge(Pack && ... pack)
        {
            return Imp::MergedParameters<Pack...>(std::forward<Pack>(pack)...);
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
