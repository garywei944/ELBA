#ifndef ELBA_XDROPALIGNER_HPP
#define ELBA_XDROPALIGNER_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

class XDropAligner : public PairwiseFunction
{
public:

  XDropAligner(ScoringScheme scheme, ushort seedlen, int xdrop);

  void apply_batch(seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
                   seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
                   uint64_t *lids,
                   uint64_t col_offset,
                   uint64_t row_offset,
                   PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
                   std::ofstream &lfs,
                   const bool noAlign,
                   ushort k,
                   uint64_t nreads,
                   std::vector<int64_t>& ContainedSeqPerProc,
                   float ratioScoreOverlap = 0.99,
                   int debugThr = 50) override;
private:
  ScoringScheme scheme;
  ushort seedlen;
  int xdrop;
};

#endif
