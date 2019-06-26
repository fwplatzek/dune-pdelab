// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH
#define DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH

#include <cstddef>
#include <optional>

#include <dune/common/iteratorfacades.hh>
#include <dune/common/reservedvector.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/identitygeometry.hh>

namespace Dune::PDELab::Experimental {

  template<typename Geometry>
  class ElementEmbedding
  {

  public:

    using size_type                        = std::size_t;
    using Global                           = Geometry;
    using Field                            = typename Geometry::ctype;
    using Element                          = IdentityGeometry<Field,Geometry::mydimension>;
    using Local                            = Element;
    using Inside                           = Element;
    using LocalCoordinate                  = typename Geometry::LocalCoordinate;
    using ElementCoordinate                = LocalCoordinate;
    using GlobalCoordinate                 = typename Geometry::GlobalCoordinate;
    using JacobianTransposed               = typename Geometry::JacobianTransposed;
    using JacobianInverseTransposed        = typename Geometry::JacobianInverseTransposed;
    using InsideJacobianTransposed         = JacobianTransposed;
    using InsideJacobianInverseTransposed  = JacobianInverseTransposed;
    using OutsideJacobianTransposed        = int;
    using OutsideJacobianInverseTransposed = int;
    using EmbeddingDescriptor              = ReservedVector<
      LocalCoordinate,
      (1 << (Geometry::mydimension-1))
      >;

    static constexpr int dimLocal = Geometry::mydimension;
    static constexpr int dimWorld = Geometry::coorddimension;

    const Global& global() const
    {
      return *_global;
    }

    Element local() const
    {
      return Local(global().type());
    }

    Element inside() const
    {
      return local();
    }

    LocalCoordinate inside(const LocalCoordinate& local) const
    {
      return local;
    }

    template<typename P>
    InsideJacobianTransposed jacobianTransposed(const P& p) const
    {
      return global().jacobianTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianInverseTransposed jacobianInverseTransposed(const P& p) const
    {
      return global().jacobianInverseTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianTransposed insideJacobianTransposed(const P& p) const
    {
      return global().jacobianTransposed(p.inside());
    }

    template<typename P>
    InsideJacobianInverseTransposed insideJacobianInverseTransposed(const P& p) const
    {
      return global().jacobianInverseTransposed(p.inside());
    }

    template<typename P>
    Field integrationElement(const P& p) const
    {
      return global().integrationElement(p.local());
    }

    ElementEmbedding(const Geometry& geo)
      : _global(&geo)
    {}

    EmbeddingDescriptor insideDescriptor() const
    {
      return {};
    }

  private:

    const Global* _global;

  };



  template<typename Embedding_, typename Derived>
  class GeometryCachingPoint
  {

  public:

    using Embedding                        = Embedding_;
    using Field                            = typename Embedding::Field;
    using LocalCoordinate                  = typename Embedding::LocalCoordinate;
    using ElementCoordinate                = typename Embedding::ElementCoordinate;
    using GlobalCoordinate                 = typename Embedding::GlobalCoordinate;
    using JacobianTransposed               = typename Embedding::JacobianTransposed;
    using JacobianInverseTransposed        = typename Embedding::JacobianInverseTransposed;
    using InsideJacobianTransposed         = typename Embedding::InsideJacobianTransposed;
    using InsideJacobianInverseTransposed  = typename Embedding::InsideJacobianInverseTransposed;
    using OutsideJacobianTransposed        = typename Embedding::OutsideJacobianTransposed;
    using OutsideJacobianInverseTransposed = typename Embedding::OutsideJacobianInverseTransposed;

    const LocalCoordinate& local() const
    {
      return static_cast<const Derived&>(*this).local();
    }

    operator const LocalCoordinate&() const
    {
      return local();
    }

    Field integrationElement() const
    {
      if (!_integration_element)
        return _integration_element.emplace(embedding().integrationElement(*this));
      else
        return *_integration_element;
    }

    const ElementCoordinate& element() const
    {
      return inside();
    }

    const ElementCoordinate& inside() const
    {
      if (!_inside)
        return _inside.emplace(embedding().inside(local()));
      else
        return *_inside;
    }

    const ElementCoordinate& outside() const
    {
      if (!_outside)
        return _outside.emplace(embedding().outside(local()));
      else
        return *_outside;
    }

    const GlobalCoordinate& global() const
    {
      if (!_global)
        return _global.emplace(embedding().global(local()));
      else
        return *_global;
    }

    const JacobianTransposed& jacobianTransposed() const
    {
      if (!_jacobian_transposed)
        return _jacobian_transposed.emplace(embedding().jacobianTransposed(*this));
      else
        return *_jacobian_transposed;
    }

    const JacobianInverseTransposed& jacobianInverseTransposed() const
    {
      if (!_jacobian_inverse_transposed)
        return _jacobian_inverse_transposed.emplace(embedding().jacobianInverseTransposed(*this));
      else
        return *_jacobian_inverse_transposed;
    }

    const InsideJacobianTransposed& insideJacobianTransposed() const
    {
      if (!_inside_jacobian_transposed)
        return _inside_jacobian_transposed.emplace(embedding().insideJacobianTransposed(*this));
      else
        return *_inside_jacobian_transposed;
    }

    const InsideJacobianInverseTransposed& insideJacobianInverseTransposed() const
    {
      if (!_inside_jacobian_inverse_transposed)
        return _inside_jacobian_inverse_transposed.emplace(embedding().insideJacobianInverseTransposed(*this));
      else
        return *_inside_jacobian_inverse_transposed;
    }

    const OutsideJacobianTransposed& outsideJacobianTransposed() const
    {
      if (!_outside_jacobian_transposed)
        return _outside_jacobian_transposed.emplace(embedding().outsideJacobianTransposed(*this));
      else
        return *_outside_jacobian_transposed;
    }

    const OutsideJacobianInverseTransposed& outsideJacobianInverseTransposed() const
    {
      if (!_outside_jacobian_inverse_transposed)
        return _outside_jacobian_inverse_transposed.emplace(embedding().outsideJacobianInverseTransposed(*this));
      else
        return *_outside_jacobian_inverse_transposed;
    }

    const GlobalCoordinate& unitOuterNormal() const
    {
      if (!_unit_outer_normal)
        return _unit_outer_normal.emplace(embedding().unitOuterNormal(*this));
      else
        return *_unit_outer_normal;
    }

    Embedding embedding() const
    {
      return static_cast<const Derived&>(*this).embedding();
    }

    friend std::ostream& operator<<(std::ostream& os, const GeometryCachingPoint& p)
    {
      os << "EP(local=" << p.local()
         << ", inside=" << p.inside()
         << ", global=" << p.global()
         << ")";
      return os;
    }

  private:

    mutable std::optional<ElementCoordinate> _inside;
    mutable std::optional<ElementCoordinate> _outside;
    mutable std::optional<GlobalCoordinate> _global;
    mutable std::optional<Field> _integration_element;
    mutable std::optional<JacobianTransposed> _jacobian_transposed;
    mutable std::optional<JacobianInverseTransposed> _jacobian_inverse_transposed;
    mutable std::optional<InsideJacobianTransposed> _inside_jacobian_transposed;
    mutable std::optional<InsideJacobianInverseTransposed> _inside_jacobian_inverse_transposed;
    mutable std::optional<OutsideJacobianTransposed> _outside_jacobian_transposed;
    mutable std::optional<OutsideJacobianInverseTransposed> _outside_jacobian_inverse_transposed;
    mutable std::optional<GlobalCoordinate> _unit_outer_normal;

  };



  template<typename Embedding_>
  class EmbeddedPoint
    : public GeometryCachingPoint<Embedding_,EmbeddedPoint<Embedding_>>
  {

  public:

    using Embedding        = Embedding_;
    using LocalCoordinate  = typename Embedding::Global::LocalCoordinate;

    const LocalCoordinate& local() const
    {
      return _local;
    }

    const Embedding& embedding() const
    {
      return _embedding;
    }

    EmbeddedPoint(const LocalCoordinate& local, const Embedding& embedding)
      : _local(local)
      , _embedding(embedding)
    {}

  private:

    LocalCoordinate _local;
    Embedding _embedding;

  };


  template<typename Rule_>
  class QuadraturePoint
  {

  public:

    using Rule                      = Rule_;
    using Embedding                 = typename Rule::Embedding;
    using Native                    = typename Rule::Native::value_type;
    using LocalCoordinate           = typename Rule::LocalCoordinate;
    using ElementCoordinate         = typename Rule::ElementCoordinate;
    using GlobalCoordinate          = typename Rule::GlobalCoordinate;
    using Index                     = typename Rule::Index;
    using Field                     = typename Rule::Field;
    using JacobianTransposed        = typename Embedding::JacobianTransposed;
    using JacobianInverseTransposed = typename Embedding::JacobianInverseTransposed;
    using InsideJacobianTransposed         = typename Embedding::InsideJacobianTransposed;
    using InsideJacobianInverseTransposed  = typename Embedding::InsideJacobianInverseTransposed;
    using OutsideJacobianTransposed        = typename Embedding::OutsideJacobianTransposed;
    using OutsideJacobianInverseTransposed = typename Embedding::OutsideJacobianInverseTransposed;

    const Native& native() const
    {
      return _qp;
    }

    const LocalCoordinate& local() const
    {
      return _qp.position();
    }

    operator const LocalCoordinate&() const
    {
      return _qp.position();
    }

    Index index() const
    {
      return _index;
    }

    Field ruleWeight() const
    {
      return _qp.weight();
    }

    Field weight() const
    {
      return ruleWeight() * integrationElement();
    }

    Field integrationElement() const
    {
      if (!_integration_element)
        return _integration_element.emplace(_rule.integrationElement(*this));
      else
        return *_integration_element;
    }

    const ElementCoordinate& element() const
    {
      return inside();
    }

    const ElementCoordinate& inside() const
    {
      if (!_inside)
        return _inside.emplace(_rule.inside(local()));
      else
        return *_inside;
    }

    const ElementCoordinate& outside() const
    {
      if (!_outside)
        return _outside.emplace(_rule.outside(local()));
      else
        return *_outside;
    }

    const GlobalCoordinate& global() const
    {
      if (!_global)
        return _global.emplace(_rule.global(local()));
      else
        return *_global;
    }

    const JacobianTransposed& jacobianTransposed() const
    {
      if (!_jacobian_transposed)
        return _jacobian_transposed.emplace(_rule.jacobianTransposed(*this));
      else
        return *_jacobian_transposed;
    }

    const JacobianInverseTransposed& jacobianInverseTransposed() const
    {
      if (!_jacobian_inverse_transposed)
        return _jacobian_inverse_transposed.emplace(_rule.jacobianInverseTransposed(*this));
      else
        return *_jacobian_inverse_transposed;
    }

    const InsideJacobianTransposed& insideJacobianTransposed() const
    {
      if (!_inside_jacobian_transposed)
        return _inside_jacobian_transposed.emplace(_rule.insideJacobianTransposed(*this));
      else
        return *_inside_jacobian_transposed;
    }

    const InsideJacobianInverseTransposed& insideJacobianInverseTransposed() const
    {
      if (!_inside_jacobian_inverse_transposed)
        return _inside_jacobian_inverse_transposed.emplace(_rule.insideJacobianInverseTransposed(*this));
      else
        return *_inside_jacobian_inverse_transposed;
    }

    const OutsideJacobianTransposed& outsideJacobianTransposed() const
    {
      if (!_outside_jacobian_transposed)
        return _outside_jacobian_transposed.emplace(_rule.outsideJacobianTransposed(*this));
      else
        return *_outside_jacobian_transposed;
    }

    const OutsideJacobianInverseTransposed& outsideJacobianInverseTransposed() const
    {
      if (!_outside_jacobian_inverse_transposed)
        return _outside_jacobian_inverse_transposed.emplace(_rule.outsideJacobianInverseTransposed(*this));
      else
        return *_outside_jacobian_inverse_transposed;
    }

    const GlobalCoordinate& unitOuterNormal() const
    {
      if (!_unit_outer_normal)
        return _unit_outer_normal.emplace(_rule.unitOuterNormal(*this));
      else
        return *_unit_outer_normal;
    }

    const Rule& rule() const
    {
      return _rule;
    }

    Embedding embedding() const
    {
      return _rule.embedding();
    }

    friend std::ostream& operator<<(std::ostream& os, const QuadraturePoint& p)
    {
      os << "QP(local=" << p.local()
         << ", inside=" << p.inside()
         << ", global=" << p.global()
         << ", weight=" << p.weight()
         << ")";
      return os;
    }

    QuadraturePoint(const Native& qp, Index index, const Rule& rule)
      : _qp(qp)
      , _index(index)
      , _rule(rule)
    {}

  private:

    const Native& _qp;
    Index _index;
    mutable std::optional<ElementCoordinate> _inside;
    mutable std::optional<ElementCoordinate> _outside;
    mutable std::optional<GlobalCoordinate> _global;
    mutable std::optional<Field> _integration_element;
    mutable std::optional<JacobianTransposed> _jacobian_transposed;
    mutable std::optional<JacobianInverseTransposed> _jacobian_inverse_transposed;
    mutable std::optional<InsideJacobianTransposed> _inside_jacobian_transposed;
    mutable std::optional<InsideJacobianInverseTransposed> _inside_jacobian_inverse_transposed;
    mutable std::optional<OutsideJacobianTransposed> _outside_jacobian_transposed;
    mutable std::optional<OutsideJacobianInverseTransposed> _outside_jacobian_inverse_transposed;
    mutable std::optional<GlobalCoordinate> _unit_outer_normal;
    const Rule& _rule;

  };


  template<typename Context_, typename QR, typename Embedding_>
  class QuadratureRule
  {

    friend class QuadraturePoint<QuadratureRule>;

  public:

    using Context           = Context_;
    using Native            = QR;
    using Embedding         = Embedding_;
    using Field             = typename Native::value_type::Field;
    using LocalCoordinate   = typename Native::value_type::Vector;
    using ElementCoordinate = typename Embedding::Element::GlobalCoordinate;
    using GlobalCoordinate  = typename Embedding::Global::GlobalCoordinate;

    //! The size type used by the container.
    using size_type         = typename Native::size_type;
    using Index             = size_type;

    using QuadraturePoint   = PDELab::Experimental::QuadraturePoint<QuadratureRule>;
    using value_type        = QuadraturePoint;

    class iterator
      : public RandomAccessIteratorFacade<iterator,QuadraturePoint,QuadraturePoint>
    {

      friend class RandomAccessIteratorFacade<iterator,QuadraturePoint,QuadraturePoint>;

      using NativeIterator = typename Native::const_iterator;

    public:

      // Add support for returning non-references from iterator.
      // We need a little bit of magic to make operator->() work for this iterator
      // because we return a temporary object from dereference(), and the standard
      // implementation of operator->() in the facade tries to take the address of
      // that temporary, which the compiler will vehemently object to... ;-)
      //
      // So I borrowed the following neat little trick from Boost's iterator library:
      // The proxy object stores a copy of the temporary View object, and operator()->
      // returns the proxy object to the caller. As mandated by the standard, the compiler
      // will then attempt to repeat the operator->() on the returned object and get the
      // address of the copy stored in the (temporary) proxy object. That proxy object
      // is guaranteed to live until the next sequence point, and that is precisely as
      // long as we have to guarantee the validity of the pointer to our View object.
      // Problem solved - and another example of how difficult it is to get this low-level
      // stuff implemented on the same level as Boost...
      struct proxy
      {

        explicit proxy(QuadratureRule&& v)
          : _tmp(v)
        {}

        QuadratureRule* operator->()
        {
          return &_tmp;
        }

        QuadratureRule _tmp;
      };

      // The proxy object will stand in as a pointer
      using pointer = proxy;

      iterator()
        : _iterator()
        , _index(0)
        , _rule(nullptr)
      {}

      iterator(NativeIterator it, std::size_t index, const QuadratureRule& rule)
        : _iterator(it)
        , _index(index)
        , _rule(&rule)
      {}

      bool equals(const iterator& other) const
      {
        return _iterator == other._iterator;
      }

      void increment()
      {
        ++_iterator;
        ++_index;
      }

      void decrement()
      {
        --_iterator;
        --_index;
      }

      void advance(int n)
      {
        _iterator += n;
        _index += n;
      }

      std::ptrdiff_t distanceTo(iterator& other) const
      {
        return other._iterator - _iterator;
      }

      QuadraturePoint dereference() const
      {
        return QuadraturePoint(*_iterator,_index,*_rule);
      }

      pointer operator->() const
      {
        return pointer(dereference());
      }

    private:

      NativeIterator _iterator;
      std::size_t _index;
      const QuadratureRule* _rule;

    };

    //! A const iterator over the quadrature points.
    using const_iterator   = iterator;

    //! Returns the maximum polynomial order up to which this rule is exact.
    int order() const
    {
      return _quadrature_rule->order();
    }

    //! Returns the geometry type that this rule is valid for.
    GeometryType type() const
    {
      return _quadrature_rule->type();
    }

    //! Returns the number of quadrature points.
    size_type size() const
    {
      return _quadrature_rule->size();
    }

    //! Returns an iterator pointing to the first quadrature point.
    const_iterator begin()
    {
      assert(not _started);
      _started = true;
      _ctx->beginQuadrature(*this);
      auto it = _quadrature_rule->begin();
      return {it,0,*this};
    }

    const_iterator util_begin()
    {
      assert(_started);
      auto it = _quadrature_rule->begin();
      return {it,0,*this};
    }

    //! Returns an iterator pointing after the last quadrature point.
    const_iterator end()
    {
      return {_quadrature_rule->end(),size(),*this};
    }

#ifndef DOXYGEN

    QuadratureRule(Context& ctx, const QR& quadrature_rule, Embedding&& embedding)
      : _ctx(&ctx)
      , _started(false)
      , _quadrature_rule(&quadrature_rule)
      , _embedding(std::move(embedding))
    {}

    ~QuadratureRule()
    {
      if (_started)
        _ctx->endQuadrature(*this);
    }

#endif

    Embedding embedding() const
    {
      return _embedding;
    }

  private:

    ElementCoordinate inside(const LocalCoordinate& coord) const
    {
      return _embedding.inside().global(coord);
    }

    ElementCoordinate outside(const LocalCoordinate& coord) const
    {
      return _embedding.outside().global(coord);
    }

    GlobalCoordinate global(const LocalCoordinate& coord) const
    {
      return _embedding.global().global(coord);
    }

    Field integrationElement(const QuadraturePoint& qp) const
    {
      return _embedding.global().integrationElement(qp.local());
    }

    auto jacobianTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.global().jacobianTransposed(qp.local());
    }

    auto jacobianInverseTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.global().jacobianInverseTransposed(qp.local());
    }

    auto insideJacobianTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.insideJacobianTransposed(qp);
    }

    auto insideJacobianInverseTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.insideJacobianInverseTransposed(qp);
    }

    auto outsideJacobianTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.outsideJacobianTransposed(qp);
    }

    auto outsideJacobianInverseTransposed(const QuadraturePoint& qp) const
    {
      return _embedding.outsideJacobianInverseTransposed(qp);
    }

    auto unitOuterNormal(const QuadraturePoint& qp) const
    {
      return _embedding.unitOuterNormal(qp);
    }

    Context* _ctx;
    bool _started;
    const QR* _quadrature_rule;
    Embedding _embedding;

  };

} // namespace Dune::PDELab::Experimental

#endif // DUNE_PDELAB_ASSEMBLER_QUADRATURERULE_HH