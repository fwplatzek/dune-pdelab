// -*- tab-width: 2; indent-tabs-mode: nil -*-
// vi: set et ts=2 sw=2 sts=2:

#ifndef DUNE_PDELAB_IMEX_ONESTEP_PARAMETER_HH
#define DUNE_PDELAB_IMEX_ONESTEP_PARAMETER_HH

/**
 * \author Pavel Hron, Marian Piatkowski
 * \file
 * \brief Implementation of the IMEX time stepping parameter classes
 */

#include <string>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune {
  namespace PDELab {

    //! Base parameter class for time stepping scheme parameters
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXTimeSteppingParameterInterface
    {
    public:
      typedef R RealType;

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const = 0;

      /*! \brief Return number of stages of the method
       */
      virtual unsigned s () const = 0;

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const = 0;

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const = 0;

      /*! \brief Return entries of the E matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const = 0;

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,r
      */
      virtual R d (int r) const = 0;

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,r
      */
      virtual R ed (int r) const = 0;

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const = 0;

      //! every abstract base class has a virtual destructor
      virtual ~IMEXTimeSteppingParameterInterface () {}
    };


    //! Parameters specifying IMEX-Euler method
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXEulerParameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXEulerParameter ()
      {
        D[0] = 0.0;  D[1] = 1.0;
        ED[0] = 0.0;  ED[1] = 1.0;
        A[0][0] = -1.0; A[0][1] = 1.0;
        B[0][0] = 0.0;  B[0][1] = 1.0;
        EB[0][0] = 1.0;  EB[0][1] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 1;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Euler");
      }

    private:
      Dune::FieldVector<R,2> D;
      Dune::FieldVector<R,2> ED;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
      Dune::FieldMatrix<R,1,2> EB;
    };

    //! Parameters specifying IMEX-Trapez method
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXTrapezParameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXTrapezParameter ()
      {
        D[0] = 0.0;  D[1] = 1.0;  D[2] = 1.0;
        ED[0] = 0.0;  ED[1] = 1.0;  ED[2] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0;
        A[1][0] = -0.5; A[1][1] = -0.5; A[1][2] = 1.0;

        B[0][0] = 0.5;  B[0][1] = 0.5; B[0][2] = 0.0;
        B[1][0] = 0.25;  B[1][1] = 0.25; B[1][2] = 0.0;

        EB[0][0] = 1.0;  EB[0][1] = 0.0; EB[0][2] = 0.0;
        EB[1][0] = 0.0;  EB[1][1] = 0.5; EB[1][2] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 2;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Trapez");
      }

    private:
      Dune::FieldVector<R,3> D;
      Dune::FieldVector<R,3> ED;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
      Dune::FieldMatrix<R,2,3> EB;
    };


    //! Parameters specifying IMEX-theta method [Koto2006]
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXThetaParameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXThetaParameter (R theta_=0.5):
        theta(theta_)
      {
        D[0] = 0.0;  D[1] = 1.0;
        ED[0] = 0.0;  ED[1] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0;

        B[0][0] = 1-theta;  B[0][1] = theta;

        EB[0][0] = 1.0;  EB[0][1] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 1;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Theta");
      }

    private:
      R theta;
      Dune::FieldVector<R,2> D;
      Dune::FieldVector<R,2> ED;
      Dune::FieldMatrix<R,1,2> A;
      Dune::FieldMatrix<R,1,2> B;
      Dune::FieldMatrix<R,1,2> EB;
    };

    //! Parameters specifying IMEX-2-order, 3 stage, method [Pareschi2004]??
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXPareschi2Parameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXPareschi2Parameter ()
      {


        D[0] = 0.0;  D[1] = 1;  D[2] = 0.5; D[3] = 1.0;
        ED[0] = 0.0;  ED[1] = 1;  ED[2] = 0.5; ED[3] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0;  A[0][3] = 0.0;
        A[1][0] = -0.5; A[1][1] = -0.5; A[1][2] = 1.0; A[1][3] = 0.0;
        A[2][0] = -1./3.; A[2][1] = -1./3.; A[2][2] = -1./3.; A[2][3] = 1;

        B[0][0] = 0.;  B[0][1] =1; B[0][2] = 0.0; B[0][3] = 0.0;
        B[1][0] = 0.;  B[1][1] =-1; B[1][2] = 1.; B[1][3] = 0.0;
        B[2][0] = 0;  B[2][1] =-7./6.; B[2][2] = 2./3.; B[2][3] = 1;

        EB[0][0] = 1.;  EB[0][1] =0; EB[0][2] = 0.0; EB[0][3] = 0.0;
        EB[1][0] = 0.;  EB[1][1] =0.; EB[1][2] = 0.; EB[1][3] = 0.0;
        EB[2][0] = -0.5;  EB[2][1] =0; EB[2][2] = 1.; EB[2][3] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 3;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Pareschi2");
      }

    private:
      Dune::FieldVector<R,4> D;
      Dune::FieldVector<R,4> ED;
      Dune::FieldMatrix<R,3,4> A;
      Dune::FieldMatrix<R,3,4> B;
      Dune::FieldMatrix<R,3,4> EB;
    };

    //! Parameters specifying IMEX-Alexander2 method [Ascher97]
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXAlexander2Parameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXAlexander2Parameter ()
      {
        const R alpha = 2-std::sqrt(2.0);
        const R gamma = 0.5*alpha;
        const R delta = 1-1./alpha;

        D[0] = 0.0;  D[1] = gamma;  D[2] = 1.0;
        ED[0] = 0.0;  ED[1] = gamma;  ED[2] = 1.0;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0;
        A[1][0] = -0.5; A[1][1] = -0.5; A[1][2] = 1.0;

        B[0][0] = 0.;  B[0][1] = gamma; B[0][2] = 0.0;
        B[1][0] = 0.;  B[1][1] = 1-3.*gamma/2.; B[1][2] = gamma;

        EB[0][0] = gamma;  EB[0][1] = 0.0; EB[0][2] = 0.0;
        EB[1][0] = delta-gamma/2.;  EB[1][1] = 1-delta; EB[1][2] = 0.0;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 2;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Alexander2");
      }

    private:
      Dune::FieldVector<R,3> D;
      Dune::FieldVector<R,3> ED;
      Dune::FieldMatrix<R,2,3> A;
      Dune::FieldMatrix<R,2,3> B;
      Dune::FieldMatrix<R,2,3> EB;
    };


    //! Parameters specifying IMEX-RK, order 3, four-stage [Ascher97]
    /**
     * \tparam R C++ type of the floating point parameters
     */
    template<class R>
    class IMEXAscher3Parameter : public IMEXTimeSteppingParameterInterface<R>
    {
    public:

      IMEXAscher3Parameter ()
      {


        D[0] = 0.0;  D[1] = 1./2.;  D[2] = 2./3.; D[3] = 1./2.;  D[4] = 1.;
        ED[0] = 0.0;  ED[1] = 1./2.;  ED[2] = 2./3.; ED[3] = 1./2.;  ED[4] = 1.;

        A[0][0] = -1.0; A[0][1] = 1.0; A[0][2] = 0.0;  A[0][3] = 0.0; A[0][4] = 0.0;
        A[1][0] = -0.5; A[1][1] = -0.5; A[1][2] = 1.0; A[1][3] = 0.0; A[1][4] = 0.0;
        A[2][0] = -1./3.; A[2][1] = -1./3.; A[2][2] = -1./3.;  A[2][3] = 1.0; A[2][4] = 0.0;
        A[3][0] = -1./4.; A[3][1] = -1./4.; A[3][2] = -1./4.; A[3][3] = -1./4.; A[3][4] = 1.0;


        B[0][0] = 0.;  B[0][1] = 1./2.; B[0][2] = 0.; B[0][3] = 0.; B[0][4] = 0.;
        B[1][0] = 0.;  B[1][1] = -1./12.; B[1][2] = 1./2.; B[1][3] = 0.; B[1][4] = 0.;
        B[2][0] = 0.;  B[2][1] = -13./18.; B[2][2] = 1./3.; B[2][3] = 1./2.; B[2][4] = 0.;
        B[3][0] = 0.;  B[3][1] = 35./24.; B[3][2] = -7./4.; B[3][3] = 3./8.; B[3][4] = 1./2;

        EB[0][0] = 1./2.;  EB[0][1] = 0.; EB[0][2] = 0.; EB[0][3] = 0.; EB[0][4] = 0.;
        EB[1][0] = 13./36.;  EB[1][1] = 1./18.; EB[1][2] = 0.; EB[1][3] = 0.; EB[1][4] = 0.;
        EB[2][0] = 25./54.;  EB[2][1] = -23./27.; EB[2][2] = 1./2.; EB[2][3] = 0.; EB[2][4] = 0.;
        EB[3][0] = -17./72;  EB[3][1] = 35./18.; EB[3][2] = 5./8.; EB[3][3] = -7./4.; EB[3][4] = 0.;
      }

      /*! \brief Return true if method is implicit
       */
      virtual bool implicit () const
      {
        return true;
      }

      /*! \brief Return number of stages s of the method
       */
      virtual unsigned s () const
      {
        return 4;
      }

      /*! \brief Return entries of the A matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R a (int r, int i) const
      {
        return A[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R b (int r, int i) const
      {
        return B[r-1][i];
      }

      /*! \brief Return entries of the B matrix
        \note that r ∈ 1,...,s and i ∈ 0,...,r
      */
      virtual R eb (int r, int i) const
      {
        return EB[r-1][i];
      }

      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R d (int i) const
      {
        return D[i];
      }


      /*! \brief Return entries of the d Vector
        \note that i ∈ 0,...,s
      */
      virtual R ed (int i) const
      {
        return ED[i];
      }

      /*! \brief Return name of the scheme
       */
      virtual std::string name () const
      {
        return std::string("IMEX Ascher (4,4,3)");
      }

    private:
      Dune::FieldVector<R,5> D;
      Dune::FieldVector<R,5> ED;
      Dune::FieldMatrix<R,4,5> A;
      Dune::FieldMatrix<R,4,5> B;
      Dune::FieldMatrix<R,4,5> EB;
    };

  } // end namespace PDELab
} // end namespace Dune
#endif
