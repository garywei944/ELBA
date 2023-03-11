// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef ELBA_KMERINTERSECTSR_HPP
#define ELBA_KMERINTERSECTSR_HPP

#include "../ReadOverlap.hpp"
#include "../ParallelOps.hpp"
#include "../Defines.hpp"

#include <cstdlib>
#include <cmath>
#include <functional>

std::pair<int, int> distance(const std::pair<PosInRead, PosInRead>& l, const std::pair<PosInRead, PosInRead>& r) {
	return {std::abs(static_cast<int>(l.first - r.first)), std::abs(static_cast<int>(l.second - r.second))};
}

bool operator>(const std::pair<int, int>& l, const int& c) {
	if(l.first > c && l.second > c) return true;
	else return false;
}

// int begQ = (cnt == 0)? cks->first.first : cks->second.first;
// int begT = (cnt == 0)? cks->first.second : cks->second.second;

namespace elba {
  template<typename IN, typename OUT>
  struct KmerIntersect {
    static OUT id() {
      OUT a;
      return a;
    }

    static bool returnedSAID() { return false; }

    static OUT add(const OUT &arg1, const OUT &arg2)
    {
  #ifdef TWOSEED
      OUT res(arg1.count + arg2.count); /* ReadOverlap constructor will need to take count as int input */

      res.begQs[0] = arg1.begQs[0];
      res.begQs[1] = arg2.begQs[0];

      res.begTs[0] = arg1.begTs[0];
      res.begTs[1] = arg2.begTs[0];

      return res;
  #else
      #error "require TWOSEED"
  #endif
    }

    static OUT multiply(const IN &arg1, const IN &arg2)
    {
      OUT a;

  #ifdef TWOSEED
      a.begQs[0] = arg1;
      a.begTs[0] = arg2;
  #else
      #error "require TWOSEED"
  #endif

      return a;
    }

    static void axpy(IN a, const IN &x, OUT &y) {
      y = add(y, multiply(a, x));
    }

    static MPI_Op mpi_op() {
      static MPI_Op mpiop;
      static bool exists = false;
      if (exists)
        return mpiop;
      else {
        MPI_Op_create(MPI_func, true, &mpiop);
        exists = true;
        return mpiop;
      }
    }

    static void
    MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
      for (int i = 0; i < *len; ++i) {
        *((OUT) inoutvec + i) = add(*((OUT) invec + i), *((OUT) inoutvec + 1));
      }

    }
  };
}
#endif //ELBA_KMERINTERSECTSR_HPP
