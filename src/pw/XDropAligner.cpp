
#include "../../include/pw/XDropAligner.hpp"

XDropAligner::XDropAligner(ScoringScheme scheme, ushort seedlen, int xdrop) :
PairwiseFunction(), scheme(scheme), seedlen(seedlen), xdrop(xdrop) {}

void XDropAligner::apply_batch(seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
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
                               float ratioScoreOverlap,
                               int debugThr)
{
}
