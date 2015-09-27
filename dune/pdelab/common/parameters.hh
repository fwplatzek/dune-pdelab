#ifndef DUNE_PDELAB_COMMON_PARAMETERS_HH
#define DUNE_PDELAB_COMMON_PARAMETERS_HH

/**
   \file parameters.hh Add support for functions style definition of LocalOperator parameter classes

   \todo implement bind and unbind for elements
   \todo implement bind and unbind for intersection
   \todo implement setTime
   \todo share storage for time, entity, intersection, where possible
 */

#include <dune/functions/common/concept.hh>
#include <utility>
#include <type_traits>

namespace Dune {
namespace PDELab {
    namespace Parameters {
        namespace Imp {

            struct NoParameter {};

            // concepts for local functions

            struct HasBindMethod
            {
                template<class F, class C>
                auto require(F&& f, C&& c) -> decltype(
                    f.bind(c)
                    );
            };

            struct HasUnbindMethod
            {
                template<class F>
                auto require(F&& f) -> decltype(
                    f.unbind()
                    );
            };

            // forward calls to the local function _iff_ the appropriate method is available

            template<typename F, typename C,
                     typename std::enable_if< Functions::Concept::models<HasBindMethod, F,C>(), int>::type = 0>
            void forwardBind(F && f, C && c)
            {
                f.bind(std::forward<C>(c));
            }

            template<typename F, typename C,
                     typename std::enable_if< not Functions::Concept::models<HasBindMethod, F,C>(), int>::type = 0>
            void forwardBind(F && f, C && c) {}

            template<typename F,
                     typename std::enable_if< Functions::Concept::models<HasUnbindMethod, F>(), int>::type = 0>
            void forwardUnbind(F && f)
            {
                f.unbind();
            }

            template<typename F,
                     typename std::enable_if< not Functions::Concept::models<HasUnbindMethod, F>(), int>::type = 0>
            void forwardUnbind(F && f) {}

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
                template<typename C>                                    \
                    void bind(C && c) {                                 \
                    ::Dune::PDELab::Parameters::Imp::forwardBind(Variable, \
                        std::forward<C>(c));                            \
                    ::Dune::PDELab::Parameters::Imp::forwardBind(static_cast<Base&>(*this), \
                        std::forward<C>(c));                            \
                }                                                       \
                void unbind() {                                         \
                    ::Dune::PDELab::Parameters::Imp::forwardUnbind(Variable); \
                    ::Dune::PDELab::Parameters::Imp::forwardUnbind(static_cast<Base&>(*this)); \
                }                                                       \
            };                                                          \
            template<typename T>                                        \
            Name##Parameter<T> define##Name(T && t)                     \
            {                                                           \
                return { std::forward<T>(t) };                          \
            }

            template<typename Base, typename Head>
            auto
            merge(const Base & b, const Head & head)
                -> typename Head::template bindBaseClass<Base>
            {
                return typename Head::template bindBaseClass<Base>(b,head);
            };

            template<typename Base, typename Head, typename... Pack>
            auto
            merge(const Base & b, const Head & head, Pack && ... pack)
                -> typename Head::template bindBaseClass<decltype(merge(b,pack...))>
            {
                // forward the base to the end of our chain
                auto remainder = merge(b,std::forward<Pack>(pack)...);
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
        auto
        merge(Pack && ... pack)
            -> decltype(Imp::merge(Imp::NoParameter(),std::forward<Pack>(pack)...))
        {
            return Imp::merge(Imp::NoParameter(),std::forward<Pack>(pack)...);
            // return Imp::merge(std::forward<Pack> pack...);
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
        template<typename Base, typename... Pack>
        auto
        overwrite(const Base & b, Pack && ... pack)
            -> decltype(Imp::merge(b,std::forward<Pack>(pack)...))
        {
            return Imp::merge(b,std::forward<Pack>(pack)...);
        };
    }

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_COMMON_PARAMETERS_HH
