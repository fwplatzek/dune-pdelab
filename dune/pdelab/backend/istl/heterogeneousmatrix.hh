#ifndef DUNE_HETEROGENEOUSMATRIX_HH
#define DUNE_HETEROGENEOUSMATRIX_HH

#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility>
#include <dune/common/tuples.hh>

#include <dune/pdelab/common/typetree/utility.hh>
#include <dune/pdelab/backend/istl/heterogeneousvector.hh>

// forward decl
namespace Dune {
  template<typename... T>
  class HeterogeneousMatrix;
}

#include <dune/istl/gsetc.hh>

namespace Dune {

  namespace impl {

    template<typename A, typename X, typename Y>
    void umv_row(const A& a, const X& x, Y& y, std::integral_constant<std::size_t,0>)
    {}

    template<typename A, typename X, typename Y, typename I>
    void umv_row(const A& a, const X& x, Y& y, I)
    {
      const std::size_t col = A::block_count - I::value;
      get<col>(a.storage()).umv(get<col>(x),y);
      umv_row(a,x,y,std::integral_constant<std::size_t,I::value-1>());
    }

    template<typename X>
    struct umv
    {

      template<typename ARow, typename YRow>
      int operator()(const ARow& a, YRow& y) const
      {
        umv_row(a,_x,y,std::integral_constant<std::size_t,ARow::block_count>());
        return 0;
      }

      umv(const X& x)
        : _x(x)
      {}

      const X& _x;
    };


    template<typename A, typename X, typename Y>
    void mmv_row(const A& a, const X& x, Y& y, std::integral_constant<std::size_t,0>)
    {}

    template<typename A, typename X, typename Y, typename I>
    void mmv_row(const A& a, const X& x, Y& y, I)
    {
      const std::size_t col = A::block_count - I::value;
      get<col>(a.storage()).mmv(get<col>(x),y);
      mmv_row(a,x,y,std::integral_constant<std::size_t,I::value-1>());
    }

    template<typename X>
    struct mmv
    {

      template<typename ARow, typename YRow>
      int operator()(const ARow& a, YRow& y) const
      {
        mmv_row(a,_x,y,std::integral_constant<std::size_t,ARow::block_count>());
        return 0;
      }

      mmv(const X& x)
        : _x(x)
      {}

      const X& _x;
    };


    template<typename A, typename W, typename X, typename Y>
    void usmv_row(const A& a, const W& alpha, const X& x, Y& y, std::integral_constant<std::size_t,0>)
    {}

    template<typename A, typename W, typename X, typename Y, typename I>
    void usmv_row(const A& a, const W& alpha, const X& x, Y& y, I)
    {
      const std::size_t col = A::block_count - I::value;
      get<col>(a.storage()).usmv(alpha,get<col>(x),y);
      usmv_row(a,alpha,x,y,std::integral_constant<std::size_t,I::value-1>());
    }

    template<typename X, typename W>
    struct usmv
    {

      template<typename ARow, typename YRow>
      int operator()(const ARow& a, YRow& y) const
      {
        usmv_row(a,_alpha,_x,y,std::integral_constant<std::size_t,ARow::block_count>());
        return 0;
      }

      usmv(const X& x, const W& alpha)
        : _x(x)
        , _alpha(alpha)
      {}

      const X& _x;
      const W& _alpha;

    };

  } // namespace impl



  template<typename... T>
  class HeterogeneousMatrix
  {

    typedef typename Dune::PDELab::TypeTree::index_pack_builder<sizeof...(T)>::type indices;

  public:

    typedef tuple<T...> Storage;

    typedef typename tuple_element<0,Storage>::type::field_type field_type;
    typedef typename tuple_element<0,Storage>::type::size_type size_type;

    template<size_type i>
    struct Row
    {
      typedef typename tuple_element<i,Storage>::type type;
    };

    template<size_type i, size_type j>
    struct Block
    {
      typedef typename tuple_element<i,Storage>::type::template Block<j>::type type;
    };

    static const size_type row_block_count = sizeof...(T);
    static const size_type col_block_count = Row<0>::type::block_count;

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
     * y = A x
     */
    template<typename X, typename Y>
    void mv (const X& x, Y& y) const
    {
      y = 0; //reset y (for mv uses umv)

      impl::binary_operator(_storage,y._storage,indices(),impl::umv<typename X::Storage>(x._storage));
    }

    /**
     * y += A x
     */
    template<typename X, typename Y>
    void umv (const X& x, Y& y) const
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::umv<typename X::Storage>(x._storage));
    }

    /**
     * y -= A x
     */
    template<typename X, typename Y>
    void mmv (const X& x, Y& y) const
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::mmv<typename X::Storage>(x._storage));
    }

    //! y += alpha A x
    template<typename W, typename X, typename Y>
    void usmv (const W& alpha, const X& x, Y& y) const
    {
      impl::binary_operator(_storage,y._storage,indices(),impl::usmv<typename X::Storage,W>(x._storage,alpha));
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
    typename Row<i>::type& row()
    {
      return get<i>(_storage);
    }

    template<size_type i>
    const typename Row<i>::type& row() const
    {
      return get<i>(_storage);
    }

    template<size_type i, size_type j>
    typename Block<i,j>::type& block()
    {
      return get<i>(_storage).template block<j>();
    }

    template<size_type i, size_type j>
    const typename Block<i,j>::type& block() const
    {
      return get<i>(_storage).template block<j>();
    }


  private:

    Storage _storage;

  };

} // end namespace

#endif // DUNE_HETEROGENEOUSMATRIX_HH
