// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS2_HH
#define DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS2_HH

#include<utility>
#include<vector>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/pdelab/common/unordered_map.hh>
#include <dune/pdelab/common/unordered_set.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/backend/istl/tags.hh>

namespace Dune {
  namespace PDELab {

#ifndef DOXYGEN

    namespace istl {

      template<typename RI, typename CI, typename Block, typename Container>
      typename Block::field_type&
      access_matrix_element2(tags::field_matrix_1_1, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        assert(i == -1);
        assert(j == -1);
        return b[0][0];
      }

      template<typename RI, typename CI, typename Block, typename Container>
      typename Block::field_type&
      access_matrix_element2(tags::field_matrix_n_m, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        assert(i == 0);
        assert(j == 0);
        return b[ri[0]][ci[0]];
      }

      template<typename RI, typename CI, typename Block, typename Container>
      typename Block::field_type&
      access_matrix_element2(tags::field_matrix_1_m, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        assert(i == -1);
        assert(j == 0);
        return b[0][ci[0]];
      }

      template<typename RI, typename CI, typename Block, typename Container>
      typename Block::field_type&
      access_matrix_element2(tags::field_matrix_n_1, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        assert(i == 0);
        assert(j == -1);
        return b[ri[0]][0];
      }

      template<typename RI, typename CI, typename Block, typename Container>
      void setup_matrix(tags::bcrs_matrix, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        if (b.N() > 0)
          return;
        b.setBuildMode(Block::mymode);
        b.setmymodeparameters(c.averageRowSize(),c.overflow());
        b.setSize(c.rowOrdering().childBlockCount(ri,i,depth),c.colOrdering().childBlockCount(ci,j,depth));
        std::cout << b.N() << " x " << b.M() << std::endl;
      }

      template<typename RI, typename CI, typename Block, typename Container>
      void setup_matrix(tags::field_matrix, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {}

      template<typename RI, typename CI, typename Block, typename Container>
      typename Block::field_type&
      access_matrix_element2(tags::bcrs_matrix, Block& b, const Container& c, const RI& ri, const CI& ci, int i, int j, int depth)
      {
        auto& e = b.entry(ri[i],ci[j]);
        setup_matrix(container_tag(e),e,c,ri,ci,i,j,depth);
        return access_matrix_element2(container_tag(e),e,c,ri,ci,i-1,j-1,depth+1);
      }


      template<typename Block>
      void finalize(tags::field_matrix,false_type,Block& b)
      {
      }

      template<typename Block>
      void finalize(tags::bcrs_matrix,true_type,Block& b)
      {
        for (typename Block::iterator row_it = b.begin(),
               row_end = b.end();
             row_it != row_end;
             ++row_it)
          for (typename Block::row_type::iterator col_it = row_it->begin(),
                 col_end = row_it->end();
               col_it != col_end;
               ++col_it)
            finalize(container_tag(*col_it),var_requires_pattern(*col_it),*col_it);

        if (b.buildStage() != Block::built)
          {
            auto stats = b.compress();
            std::cout << "avg = " << stats.avg << std::endl
                      << "max = " << stats.maximum << std::endl
                      << "ovt = " << stats.overflow_total << std::endl
                      << "ovr = " << stats.mem_ratio << std::endl;
          }
      }

      template<typename Block>
      void finalize(tags::bcrs_matrix,false_type,Block& b)
      {
        if (b.buildStage() != Block::built)
          {
            auto stats = b.compress();
            std::cout << "avg = " << stats.avg << std::endl
                      << "max = " << stats.maximum << std::endl
                      << "ovt = " << stats.overflow_total << std::endl
                      << "ovr = " << stats.mem_ratio << std::endl;
          }
      }



    } // namespace istl

#endif // DOXYGEN

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_BACKEND_ISTL_MATRIXHELPERS2_HH
