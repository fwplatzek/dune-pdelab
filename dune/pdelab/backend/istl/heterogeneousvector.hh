#ifndef DUNE_HETEROGENEOUSVECTOR_HH
#define DUNE_HETEROGENEOUSVECTOR_HH

#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <dune/common/tuples.hh>

#include <dune/pdelab/common/typetree/utility.hh>

// forward decl
namespace Dune {
  template<typename... T>
  class HeterogeneousVector;

  template<typename... T>
  class HeterogeneousMatrix;
}

#include <dune/istl/gsetc.hh>

namespace Dune {

  namespace impl {

    using Dune::PDELab::TypeTree::discard;
    using Dune::PDELab::TypeTree::index_pack;

    template<typename V, typename F, std::size_t... i>
    void unary_operator(V&& v, index_pack<i...> indices, F&& f)
    {
      discard(f(get<i>(std::forward<V>(v)))...);
    }

    template<typename V1, typename V2, typename F, std::size_t... i>
    void binary_operator(V1&& a, V2&& b, index_pack<i...> indices, F&& f)
    {
      discard(f(get<i>(std::forward<V1>(a)),get<i>(std::forward<V2>(b)))...);
    }

    template<typename V, typename F, std::size_t... i>
    void indexed_unary_operator(V&& v, index_pack<i...> indices, F&& f)
    {
      discard(f(get<i>(std::forward<V>(v)),i)...);
    }

    template<typename V1, typename V2, typename F, std::size_t... i>
    void indexed_binary_operator(V1&& a, V2&& b, index_pack<i...> indices, F&& f)
    {
      discard(f(get<i>(std::forward<V1>(a)),get<i>(std::forward<V2>(b)),i)...);
    }


    template<typename E>
    struct assignment
    {

      template<typename V>
      int operator()(V& v) const
      {
        v = _e;
        return 0;
      }

      assignment(const E& e)
        : _e(e)
      {}

      const E& _e;

    };


    struct addition
    {
      template<typename V>
      int operator()(V& a, const V& b)
      {
        a += b;
        return 0;
      }
    };

    struct subtraction
    {
      template<typename V>
      int operator()(V& a, const V& b)
      {
        a -= b;
        return 0;
      }
    };


    template<typename E>
    struct axpy
    {

      template<typename V>
      int operator()(V& x, const V& y)
      {
        x.axpy(_a,y);
        return 0;
      }

      axpy(const E& a)
        : _a(a)
      {}

      const E& _a;

    };

    template<typename E>
    struct scaling
    {

      template<typename V>
      int operator()(V& v) const
      {
        v *= _e;
        return 0;
      }

      scaling(const E& e)
        : _e(e)
      {}

      const E& _e;

    };


    template<typename E, std::size_t N>
    struct scalar_product
    {

      template<typename V>
      int operator()(const V& x, const V& y, std::size_t i)
      {
        _result[i] = x * y;
        return 0;
      }

      scalar_product()
      {}

      E result() const
      {
        return std::accumulate(_result,_result + N,E(0));
      }

      E _result[N];

    };

    template<typename E, std::size_t N>
    struct two_norm_squared
    {

      template<typename V>
      int operator()(const V& x, std::size_t i)
      {
        _result[i] = x.two_norm2();
        return 0;
      }

      two_norm_squared()
      {}

      E result() const
      {
        return std::accumulate(_result,_result + N,E(0));
      }

      E _result[N];

    };

  } // namespace impl



  template<typename... T>
  class HeterogeneousVector
  {

    template<typename...>
    friend class HeterogeneousMatrix;

    typedef typename Dune::PDELab::TypeTree::index_pack_builder<sizeof...(T)>::type indices;

  public:

    static const std::size_t block_count = sizeof...(T);

    typedef tuple<T...> Storage;

    typedef typename tuple_element<0,Storage>::type::field_type field_type;
    typedef typename std::size_t size_type;

    template<size_type i>
    struct Block
    {
      typedef typename tuple_element<i,Storage>::type type;
    };

    /**
     * number of elements
     */
    const int count()
    {
      return sizeof...(T);
    }

    /**
     * assignment operator
     */
    template<typename E>
    void operator=(const E& newval)
    {
      impl::unary_operator(_storage,indices(),impl::assignment<E>(newval));
    }

    /**
     * operator for vector addition
     */
    void operator+=(const HeterogeneousVector& y)
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::addition());
    }

    /**
     * operator for vector subtraction
     */
    void operator-=(const HeterogeneousVector& y)
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::subtraction());
    }

    template<typename E>
    void operator*=(const E& w)
    {
      impl::unary_operator(_storage,indices(),impl::scaling<E>(w));
    }


    // TODO: Fix operator* / dot product!
    field_type operator*(const HeterogeneousVector& y) const
    {
      impl::scalar_product<field_type,sizeof...(T)> sp;
      impl::indexed_binary_operator(_storage,y._storage,indices(),sp);
      return sp.result();
    }

    field_type dot(const HeterogeneousVector& y) const
    {
      impl::scalar_product<field_type,sizeof...(T)> sp;
      impl::indexed_binary_operator(_storage,y._storage,indices(),sp);
      return sp.result();
    }

    /**
     * two-norm^2
     */
    double two_norm2() const
    {
      impl::two_norm_squared<double,sizeof...(T)> norm;
      impl::indexed_unary_operator(_storage,indices(),norm);
      return norm.result();
    }

    /**
     * the real two-norm
     */
    double two_norm() const
    {
      return sqrt(two_norm2());
    }

    /**
     * axpy operation on this vector (*this += a * y)
     */
    template<typename E>
    void axpy(const E& a, const HeterogeneousVector& y)
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::axpy<E>(a));
    }

    Storage& storage()
    {
      return _storage;
    }

    const Storage& storage() const
    {
      return _storage;
    }

    template<size_type i>
    typename tuple_element<i,Storage>::type& block()
    {
      return get<i>(_storage);
    }

    template<size_type i>
    const typename tuple_element<i,Storage>::type& block() const
    {
      return get<i>(_storage);
    }

  private:

    Storage _storage;

  };

} // end namespace

#endif // DUNE_HETEROGENEOUSVECTOR_HH
