#ifndef DUNE_PDELAB_FUNCTION_HELPERS
#define DUNE_PDELAB_FUNCTION_HELPERS

#include <tuple>
#include <utility>
#include <dune/common/std/utility.hh>
#include <dune/functions/common/signature.hh>

namespace Dune {
namespace PDELab {
    namespace Imp {

        template<typename F>
        struct IsCallable;

#ifndef DOXYGEN
        template<typename F>
        struct IsCallable
        {
            struct yes { std::size_t dummy[2]; };
            struct no  { std::size_t dummy[1]; };

            template<typename C>
            static yes test(const decltype(&C::operator()) *);
            template<typename C>
            static no  test(...);

            enum { value = (sizeof(test<F>(0)) == sizeof(yes)) };
        };

        template<typename R, typename... D>
        struct IsCallable<R(D...)>
        {
            enum { value = true };
        };

        template<typename R, typename... D>
        struct IsCallable<R(*)(D...)>
        {
            enum { value = true };
        };
#endif

        /**
         * \brief Helper class to deduce the signature of a callable
         */
        template<class Signature, bool isCallable = IsCallable<Signature>::value >
        struct ExtendedSignatureTraits {};

#ifndef DOXYGEN
        template<class Signature>
        struct ExtendedSignatureTraits<Signature, false>
        {
            static_assert(Dune::AlwaysFalse<Signature>::value, "Trying to get the signature of something that is not callable");
        };

        /** \brief deduce the signature of the operator() of a class T */
        template<class T>
        struct ExtendedSignatureTraits<T, true>
            : public ExtendedSignatureTraits<decltype(&T::operator()), true>
        {};

        /** \brief deduce the signature of an arbitrary const member function of class C */
        template <typename C, typename R, typename... Parameters>
        struct ExtendedSignatureTraits<R(C::*)(Parameters...) const, true>
            : public ExtendedSignatureTraits<R(Parameters...), true>
        {};

        /** \brief deduce the signature of an arbitrary member function of class C */
        template <typename C, typename R, typename... Parameters>
        struct ExtendedSignatureTraits<R(C::*)(Parameters...), true>
            : public ExtendedSignatureTraits<R(Parameters...), true>
        {};

        /** \brief extract domain and range from a free functions pointer */
        template <typename R, typename... Parameters>
        struct ExtendedSignatureTraits<R(*)(Parameters...), true>
            : public ExtendedSignatureTraits<R(Parameters...), true>
        {};

        /** \brief extract domain and range from a signature (works only for free functions) */
        template<class R, typename... Params>
        struct ExtendedSignatureTraits<R(Params...), true>
        {
            using Range = R;
            using Parameters = std::tuple<Params...>;

            using RawRange = typename std::decay<Range>::type;
            using RawParameters = std::tuple<typename std::decay<Params>::type...>;
        };
#endif


        template<std::size_t O, std::size_t N, std::size_t ...S>
        struct make_index_range_impl : make_index_range_impl<O, N-1, O+N-1, S...> { };

        template<std::size_t O, std::size_t ...S>
        struct make_index_range_impl<O, 0, S...> {
            typedef Std::index_sequence<S...> type;
        };

        template<std::size_t O, std::size_t E>
        static constexpr inline typename make_index_range_impl<O, E-O>::type make_index_range ()
        {
            return typename make_index_range_impl<O, E-O>::type();
        }

        /**
           \todo actually derive the list of parameters for operator() from F
        */
        template <typename F, std::size_t N, class Imp>
        class ReplaceParameterByMember
        {
            using Fnkt = typename std::decay<F>::type;
        public:
            using Parameters = typename ExtendedSignatureTraits<Fnkt>::Parameters;
            using Range = typename ExtendedSignatureTraits<Fnkt>::RawRange;

            ReplaceParameterByMember(F && f)
                : function_(f)
            {}

            //! interleave arguments with the entity
            template<typename ... Args>
            Range
            operator()(Args && ... args) const
            {
                static_assert(sizeof...(Args) == std::tuple_size<Parameters>::value-1, "Wrong number of arguments");
                // split parameter at position N
                auto pre = make_index_range<0,N>();
                auto post = make_index_range<N,sizeof...(Args)>();
                return callFunc(pre,post,std::tuple<Args...>(std::forward<Args>(args)...));
            }

        private:
            const Imp& asImp () const {return static_cast<const Imp &>(*this);}

            template<std::size_t ... Pre,
                     std::size_t ... Post,
                     typename ... Args>
            Range callFunc(Std::index_sequence<Pre...>,
                Std::index_sequence<Post...>,
                std::tuple<Args...> && args) const
            {
                return function_(std::get<Pre>(args) ...,
                    asImp().getParameter(),
                    std::get<Post>(args) ...
                    );
            }
            template<std::size_t ... Post,
                     typename ... Args>
            Range callFunc(Std::index_sequence<>,
                Std::index_sequence<Post...>,
                std::tuple<Args...> && args) const
            {
                return function_(
                    asImp().getParameter(),
                    std::get<Post>(args) ...
                    );
            }
            template<std::size_t ... Pre,
                     typename ... Args>
            Range callFunc(Std::index_sequence<Pre...>,
                Std::index_sequence<>,
                std::tuple<Args...> && args) const
            {
                return function_(std::get<Pre>(args) ...,
                    asImp().getParameter()
                    );
            }
            template<typename ... Args>
            Range callFunc(Std::index_sequence<>,
                Std::index_sequence<>,
                const std::tuple<Args...> && args) const
            {
                return asImp().getParameter();
            }

            // modified function
            Fnkt function_;
        };

        template <typename F, std::size_t N>
        class ReplaceEntityByBind :
            public ReplaceParameterByMember<F,N,ReplaceEntityByBind<F,N>>
        {
            using Base = ReplaceParameterByMember<F,N,ReplaceEntityByBind<F,N>>;
            friend Base;
        public:
            using Base::Base;
            using Entity = typename std::tuple_element<N, typename Base::Parameters>::type;
            void bind (const Entity & entity)
            {
                entity_ = &entity;
            }
        private:
            const Entity & getParameter() const
            {
                return *entity_;
            }
            // entity type
            const Entity * entity_;
        };

        template <typename F, std::size_t N>
        class ReplaceTime :
            public ReplaceParameterByMember<F,N,ReplaceTime<F,N>>
        {
            using Base = ReplaceParameterByMember<F,N,ReplaceTime<F,N>>;
            friend Base;
        public:
            using Base::Base;
            void setTime (double t)
            {
                t_ = t;
            }
        private:
            double getParameter() const
            {
                return t_;
            }
            // entity type
            double t_;
        };

    } // end namespace Imp


    /**
       \brief construct a wrapper function to bind the value of bind(entity) to the N'th paramter of the callable

       the local function interface expects a function to be bound to
       an entity by calling f.bind(entity). Given a function directly
       takes entity as a parameter, we create new function, which
       forwards all parameters and interleaves these with the entity
       specified in bind.

       \code
       using namespace Dune::TypeTree::Indices;
       auto f = [](Element e, Domain x) { return e.geometry().global(x).two_norm(); };
       auto lf = replaceEntityByBind(f,_1);
       \endcode

       \param f a callable
       \tparam N which parameter is the entity
     */
    template<typename F, std::size_t N = 1>
    Imp::ReplaceEntityByBind<F,N-1>
    replaceEntityByBind(F && f, TypeTree::index_constant<N>)
    {
        // we follow the bind syntax and number parameters starting from 1!
        // -> thus we use N-1 :-)
        return Imp::ReplaceEntityByBind<F,N-1>(std::forward<F>(f));
    }

    template<typename F, std::size_t N = 1>
    Imp::ReplaceTime<F,N-1>
    replaceTime(F && f, TypeTree::index_constant<N>)
    {
        // we follow the bind syntax and number parameters starting from 1!
        // -> thus we use N-1 :-)
        return Imp::ReplaceTime<F,N-1>(std::forward<F>(f));
    }

} // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FUNCTION_HELPERS
