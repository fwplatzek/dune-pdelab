// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_ISTLMATRIXBACKEND_HH
#define DUNE_ISTLMATRIXBACKEND_HH

#include<utility>
#include<vector>
#include<set>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/bcrsmatrix.hh>

#include"../gridoperatorspace/localmatrix.hh"

namespace Dune {
  namespace PDELab {

	//! ISTL backend for LinearOperatorSpace
	template<int ROWBLOCKSIZE, int COLBLOCKSIZE>
	class ISTLBCRSMatrixBackend
	{
	public:
	  //! The size type
	  typedef typename Dune::BCRSMatrix< Dune::FieldMatrix<float,1,1> >::size_type size_type;

	  //! container construction
	  template<typename T, typename E>
	  class Matrix : public Dune::BCRSMatrix< Dune::FieldMatrix<E,ROWBLOCKSIZE,COLBLOCKSIZE> >
	  {
		typedef Dune::FieldMatrix<E,ROWBLOCKSIZE,COLBLOCKSIZE> M;

	  public:
		typedef typename Dune::BCRSMatrix<M>::size_type size_type;
		typedef E ElementType;
		typedef Dune::BCRSMatrix<M> BaseT;

        //! construct container
		Matrix (const T& t) 
		  : BaseT(t.globalSizeV()/ROWBLOCKSIZE,t.globalSizeU()/COLBLOCKSIZE,
				  Dune::BCRSMatrix<M>::random) 
		{
		  Pattern pattern(t.globalSizeV(),t.globalSizeU());
		  t.fill_pattern(pattern);

		  // first dummy code: build full matrix
		  for (size_t i=0; i<t.globalSizeV(); ++i)
			this->setrowsize(i,pattern[i].size());
		  this->endrowsizes();

		  for (size_t i=0; i<t.globalSizeV(); ++i)
			{
			  for (typename std::set<size_type>::iterator it=pattern[i].begin(); 
				   it!=pattern[i].end(); ++it)
				this->addindex(i,*it);
			}
		  this->endindices();

		  // code for full matrix
		  //       for (int i=0; i<t.globalSizeV(); ++i)
		  //         this->setrowsize(i,t.globalSizeU());
		  //       this->endrowsizes();

		  //       for (int i=0; i<t.globalSizeV(); ++i)
		  //         for (int j=0; j<t.globalSizeU(); ++j)
		  //           this->addindex(i,j);
		  //       this->endindices();

		}

        //! set from element
		Matrix& operator= (const E& x)
		{
		  BaseT::operator=(x);
		  return *this;
		}

		//! for debugging and AMG access
		BaseT& base ()
		{
		  return *this;
		}

		//! for debugging and AMG access
		const BaseT& base () const
		{
		  return *this;
		}
	  };

	  //! extract type of container element 
	  template<class C>
	  struct Value
	  {
		typedef typename C::field_type Type;
	  };

	  //! type to store sparsity pattern of block indices
	  class Pattern : public std::vector< std::set<size_type> >
	  {
		typedef std::vector< std::set<size_type> > BaseT;
	  public:
		Pattern (size_type m_, size_type n_)
		{
		  this->resize(m_);
		}

		void add_link (size_type i, size_type j)
		{
		  (*this)[i/ROWBLOCKSIZE].insert(j/COLBLOCKSIZE);
		}
	  };

	  //==========================
	  // individual element access
	  //==========================


	  //! get const_reference to container element
	  /**
           * \note this method does not depend on T!
           */
	  template<typename C>
	  static const typename C::field_type& const_access (const C& c, size_type i, size_type j)
	  {
		return c[i/ROWBLOCKSIZE][j/COLBLOCKSIZE][i%ROWBLOCKSIZE][j%COLBLOCKSIZE];
	  }

	  //! get non const_reference to container element 
	  /**
           * \note this method does not depend on T!
           */
	  template<typename C>
	  static typename C::field_type& access (C& c, size_type i, size_type j)
	  {
		return c[i/ROWBLOCKSIZE][j/COLBLOCKSIZE][i%ROWBLOCKSIZE][j%COLBLOCKSIZE];
	  }


	  //=================
	  // submatrix access
	  //=================

// 	  // read a submatrix given by global indices
// 	  // we can assume C to be std::vector
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void read (const C& c, 
// 						const RI& row_index, const CI& col_index, 
// 						LocalMatrix<T>& submatrix)
// 	  {
// 		submatrix.resize(row_index.size(),col_index.size());
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			{
// 			  int I = row_index[i];
// 			  int J = col_index[j];
// 			  submatrix(i,j) = c[I/ROWBLOCKSIZE][J/COLBLOCKSIZE][I%ROWBLOCKSIZE][J%COLBLOCKSIZE];
// 			}
// 	  }

// 	  // write a submatrix given by global indices
// 	  // we can assume C to be std::vector
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void write (const RI& row_index, const CI& col_index, 
// 						 const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			{
// 			  int I = row_index[i];
// 			  int J = col_index[j];
// 			  c[I/ROWBLOCKSIZE][J/COLBLOCKSIZE][I%ROWBLOCKSIZE][J%COLBLOCKSIZE] = submatrix(i,j);
// 			}
// 	  }

// 	  // write a submatrix given by global indices
// 	  // we can assume C to be std::vector
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void add (const RI& row_index, const CI& col_index, 
// 					   const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			{
// 			  int I = row_index[i];
// 			  int J = col_index[j];
// 			  c[I/ROWBLOCKSIZE][J/COLBLOCKSIZE][I%ROWBLOCKSIZE][J%COLBLOCKSIZE] += submatrix(i,j);
// 			}
// 	  }

	  //! clear one row of the matrix
	  template<typename C, typename RI>
	  static void clear_row (RI i, C& c)
	  {
		typedef typename C::ColIterator coliterator;
		for (coliterator j=c[i].begin(); j!=c[i].end(); ++j)
		  for (int jj=0; jj<COLBLOCKSIZE; jj++)
			(*j)[i%ROWBLOCKSIZE][jj] = 0;
	  }
	};

  } // namespace PDELab
} // namespace Dune

#endif
