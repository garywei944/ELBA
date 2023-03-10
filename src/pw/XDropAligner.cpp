
#include "../../include/pw/XDropAligner.hpp"
#include <vector>
#include <algorithm>

typedef enum
{
    BAD_ALIGNMENT, /* bad alignment */
    FIRST_CONTAINED,
    SECOND_CONTAINED,
    FIRST_TO_SECOND_OVERLAP,
    SECOND_TO_FIRST_OVERLAP
} overlap_class_t;

struct xseed_t
{
    int begQ, endQ, begT, endT, score;
    bool rc;

    xseed_t() : begQ(0), endQ(0), begT(0), endT(0), score(-1), rc(false) {}
};

void classify_alignment(const xseed_t& ai, int lenQ, int lenT, overlap_class_t& kind)
{
    if (ai.score <= 0)
    {
        kind = BAD_ALIGNMENT;
        return;
    }

    int begTr = ai.rc? lenT - ai.endT : ai.begT;
    int endTr = ai.rc? lenT - ai.begT : ai.endT;

    int maplen = ((ai.endT - ai.begT) + (ai.endQ - ai.begQ)) / 2;
    int overhang = std::min(ai.begQ, begTr) + std::min(lenQ - ai.endQ, lenT - endTr);
    int overlap = maplen + overhang;

    float my_thr = (1 - DELTACHERNOFF) * (0.99 * (overlap + 0.0));

    if (ai.begQ <= begTr && lenQ - ai.endQ <= lenT - endTr)
    {
        kind = FIRST_CONTAINED;
    }
    else if (ai.begQ >= begTr && lenQ - ai.endQ >= lenT - endTr)
    {
        kind = SECOND_CONTAINED;
    }
    else if (ai.score < my_thr || overlap < 500)
    {
        kind = BAD_ALIGNMENT;
    }
    else if (ai.begQ > begTr)
    {
        kind = FIRST_TO_SECOND_OVERLAP;
    }
    else
    {
        kind = SECOND_TO_FIRST_OVERLAP;
    }
}


int _extend_seed_one_direction(seqan::Dna5String& seqQ, seqan::Dna5String& seqT, bool extleft, xseed_t& xseed, int mat, int mis, int gap, int dropoff)
{
    int lenQ = length(seqQ);
    int lenT = length(seqT);

    seqan::Dna5StringReverseComplement seqTr(seqT);

    int lenQ_ext = extleft? xseed.begQ : lenQ - xseed.endQ;
    int lenT_ext = extleft? xseed.begT : lenT - xseed.endT;

    int cols = lenQ_ext + 1;
    int rows = lenT_ext + 1;

    if (rows == 1 || cols == 1) return 0;

    constexpr int int_min = std::numeric_limits<int>::min();

    int len = 2 * std::max(cols, rows);
    int min_err_score = int_min / len;
    gap = std::max(gap, min_err_score);
    mis = std::max(mis, min_err_score);
    int undef = int_min - gap;

    std::vector<int> ad1, ad2, ad3, tmp;

    int min_col = 1, max_col = 2;
    int offset1 = 0, offset2 = 0, offset3 = 0;

    ad2.push_back(0);
    int best_ext_col = 0, best_ext_row = 0, best_ext_score = 0;

    ad3.resize(2);
    ad3[0] = ad3[1] = (-gap > dropoff)? undef : gap;

    int ad_no = 1, best = 0;
    int offsetQ = xseed.endQ;
    int offsetT = xseed.endT;

    while (min_col < max_col)
    {
        ++ad_no;
        tmp = std::move(ad1);
        ad1 = std::move(ad2);
        ad2 = std::move(ad3);
        ad3 = std::move(tmp);

        offset1 = offset2;
        offset2 = offset3;
        offset3 = min_col - 1;

        ad3.resize(max_col + 1 - offset3);
        ad3[0] = ad3[max_col - offset3] = undef;

        if (ad_no * gap > best - dropoff)
        {
            if (offset3 == 0) ad3[0] = ad_no * gap;
            if (ad_no - max_col == 0) ad3[max_col - offset3] = ad_no * gap;
        }

        int ad_best = ad_no * gap;

        for (int col = min_col; col < max_col; ++col)
        {
            int i3 = col - offset3;
            int i2 = col - offset2;
            int i1 = col - offset1;

            int posQ, posT;

            posQ = extleft? cols - 1 - col : col - 1 + offsetQ;
            posT = extleft? rows - 1 + col - ad_no : ad_no - col - 1 + offsetT;

            int temp = std::max(ad2[i2-1], ad2[i2]) + gap;
            int temp2 = ad1[i1-1] + ((seqQ[posQ] == (xseed.rc? seqTr[posT] : seqT[posT]))? mat : mis);
            temp = std::max(temp, temp2);

            if (temp < best - dropoff)
            {
                ad3[i3] = undef;
            }
            else
            {
                ad3[i3] = temp;
                ad_best = std::max(ad_best, temp);
            }

            if (temp > best)
            {
                best_ext_col = col;
                best_ext_row = ad_no - best_ext_col;
                best_ext_score = ad3[best_ext_col - offset3];
                assert(best_ext_score == temp);
            }
        }

        best = std::max(best, ad_best);

        while (min_col - offset3 < ad3.size() && ad3[min_col - offset3] == undef &&
               min_col - offset2 - 1 < ad2.size() && ad2[min_col - offset2 - 1] == undef)
        {
            ++min_col;
        }

        while (max_col - offset3 > 0 && ad3[max_col - offset3 - 1] == undef && ad2[max_col - offset2 - 1] == undef)
            --max_col;

        ++max_col;

        min_col = std::max(min_col, ad_no + 2 - rows);
        max_col = std::min(max_col, cols);
    }

    int ext_col = ad3.size() + offset3 - 2;
    int ext_row = ad_no - ext_col;
    int ext_score = ad3[ext_col - offset3];

    if (ext_score == undef)
    {
        if (ad2[ad2.size()-2] != undef)
        {
            ext_col = ad2.size() + offset2 - 2;
            ext_row = ad_no - 1 - ext_col;
            ext_score = ad2[ext_col - offset2];
        }
        else if (ad2.size() > 2 && ad2[ad2.size()-3] != undef)
        {
            ext_col = ad2.size() + offset2 - 3;
            ext_row = ad_no - 1 - ext_col;
            ext_score = ad2[ext_col - offset3];
        }
    }

    if (ext_score == undef)
    {
        for (int i = 0; i < ad1.size(); ++i)
        {
            if (ad1[i] > ext_score)
            {
                ext_score = ad1[i];
                ext_col = i + offset1;
                ext_row = ad_no - 2 - ext_col;
            }
        }
    }

    if (best_ext_score != undef)
    {
        if (extleft)
        {
            xseed.begT -= best_ext_row;
            xseed.begQ -= best_ext_col;
        }
        else
        {
            xseed.endT += best_ext_row;
            xseed.endQ += best_ext_col;
        }
    }

    return best_ext_score;
}

static inline int _xdrop_seed_and_extend_l(seqan::Dna5String& seqQ, seqan::Dna5String& seqT, int mat, int mis, int gap, int dropoff, const xseed_t& xseed, int& begQ_ext, int& begT_ext)
{
    xseed_t result = xseed;

    int lscore = _extend_seed_one_direction(seqQ, seqT, true, result, mat, mis, gap, dropoff);

    begQ_ext = result.begQ;
    begT_ext = result.begT;

    return lscore;
}

static inline int _xdrop_seed_and_extend_r(seqan::Dna5String& seqQ, seqan::Dna5String& seqT, int mat, int mis, int gap, int dropoff, const xseed_t& xseed, int& endQ_ext, int& endT_ext)
{
    xseed_t result = xseed;

    int rscore = _extend_seed_one_direction(seqQ, seqT, false, result, mat, mis, gap, dropoff);

    endQ_ext = result.endQ;
    endT_ext = result.endT;

    return rscore;
}

int
_xdrop_aligner
(
    seqan::Dna5String& seqQ,
    seqan::Dna5String& seqT,
    int begQ,
    int begT,
    int seedlen,
    int mat,
    int mis,
    int gap,
    int dropoff,
    xseed_t& result
)
{
    seqan::Dna5StringReverseComplement seqTr(seqT);

    xseed_t xseed;

    int lenQ = length(seqQ);
    int lenT = length(seqT);

    if (!(seedlen&1))
        return -1;

    if (begQ < 0 || begQ + seedlen > lenQ)
        return -1;

    if (begT < 0 || begT + seedlen > lenT)
        return -1;

    if (begQ == 0 && begT == 0)
        return -1;

    bool rc = (seqQ[begQ + (seedlen>>1)] != seqT[begT + (seedlen>>1)]);

    for (int i = 0; i < seedlen; ++i)
    {
        if (seqQ[begQ + i] != (rc? seqTr[lenT - begT - seedlen + i] : seqT[begT + i]))
            return -1;
    }

    xseed.begQ = begQ;
    xseed.endQ = xseed.begQ + seedlen;

    xseed.begT = rc? lenT - begT - seedlen : begT;
    xseed.endT = xseed.begT + seedlen;

    xseed.rc = rc;

    int begQ_ext, begT_ext, lscore;
    int endQ_ext, endT_ext, rscore;

    lscore = _xdrop_seed_and_extend_l(seqQ, seqT, mat, mis, gap, dropoff, xseed, begQ_ext, begT_ext);
    rscore = _xdrop_seed_and_extend_r(seqQ, seqT, mat, mis, gap, dropoff, xseed, endQ_ext, endT_ext);

    int score = lscore + rscore + mat * seedlen;

    result.begQ = begQ_ext;
    result.endQ = endQ_ext;

    result.begT = rc? lenT - endT_ext : begT_ext;
    result.endT = rc? lenT - begT_ext : endT_ext;

    result.rc = rc;
    result.score = score;

    return score;
}


XDropAligner::XDropAligner(ScoringScheme scheme, ushort seedlen, int dropoff) :
nalignments(0), scheme(scheme), seedlen(seedlen), dropoff(dropoff) {}

void XDropAligner::apply(uint64_t l_col_idx, uint64_t g_col_idx,
                         uint64_t l_row_idx, uint64_t g_row_idx,
                         seqan::Dna5String *seq_h, seqan::Dna5String *seq_v, ushort k,
                         elba::CommonKmers &cks, std::stringstream& ss)
{}

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
                               std::vector<int64_t>& ContainedSeqPerBatch,
                               float ratioScoreOverlap,
                               int debugThr)
{
    uint64_t npairs = seqan::length(seqsh); /* the number of sequence pairs whose shared seeds we need to extend */

    std::vector<xseed_t> ai(npairs);

    int mat = seqan::scoreMatch(scheme);
    int mis = seqan::scoreMismatch(scheme);
    int gap = seqan::scoreGap(scheme);

    for (int cnt = 0; cnt < 2; ++cnt)
    {
        auto start_time = std::chrono::system_clock::now();

        for (uint64_t i = 0; i < npairs; ++i)
        {
            elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

            /* TODO: check the logic here */

            int begQ = (cnt == 0)? cks->first.first : cks->second.first;
            int begT = (cnt == 0)? cks->first.second : cks->second.second;

            xseed_t result;

            auto start_time = std::chrono::system_clock::now();
            int score = _xdrop_aligner(seqan::source(seqsv[i]), seqan::source(seqsh[i]), begQ, begT, seedlen, mat, mis, gap, dropoff, result);
            auto end_time = std::chrono::system_clock::now();

            add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

            bool rc = result.rc;
            int lenQ = length(seqan::source(seqsv[i]));
            int lenT = length(seqan::source(seqsh[i]));
            int64_t idxQ = std::get<0>(mattuples[lids[i]]) + row_offset;
            int64_t idxT = std::get<1>(mattuples[lids[i]]) + col_offset;
            int begQr = ai[i].begQ;
            int endQr = ai[i].endQ;
            int begTr = rc? lenT - ai[i].endT : ai[i].begT;
            int endTr = rc? lenT - ai[i].begT : ai[i].endT;

            if (ai[i].score < result.score)
                ai[i] = result;

            /* std::cout << idxQ << "\t" << lenQ << "\t" << begQr << "\t" << endQr << "\t" << static_cast<int>(rc) << "\t" << idxT << "\t" << lenT << "\t" << ai[i].begT << "\t" << ai[i].endT << "\t" << score << std::endl; */
        }

        auto end_time = std::chrono::system_clock::now();
        add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());
    }

    auto start_time = std::chrono::system_clock::now();

    std::vector<std::vector<int64_t>> ContainedSeqPerThread(1);

    for (uint64_t i = 0; i < npairs; ++i)
    {
        overlap_class_t kind;

        int lenQ = length(seqan::source(seqsv[i]));
        int lenT = length(seqan::source(seqsh[i]));

        int64_t idxQ = std::get<0>(mattuples[lids[i]]) + row_offset;
        int64_t idxT = std::get<1>(mattuples[lids[i]]) + col_offset;

        classify_alignment(ai[i], lenQ, lenT, kind);

        bool rc = ai[i].rc;
        int score = ai[i].score;
        int begQr = ai[i].begQ;
        int endQr = ai[i].endQ;
        int begTr = rc? lenT - ai[i].endT : ai[i].begT;
        int endTr = rc? lenT - ai[i].begT : ai[i].endT;

        elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

        cks->first.first = begQr;
        cks->first.second = endQr;
        cks->second.first = begTr;
        cks->second.second = endTr;

        cks->lenv = lenQ;
        cks->lenh = lenT;
        cks->score = score;
        cks->rc = rc;
        cks->passed = false;

        if (kind != BAD_ALIGNMENT)
        {
            if (kind == FIRST_CONTAINED)
            {
                ContainedSeqPerThread[0].push_back(idxQ);
            }
            else if (kind == SECOND_CONTAINED)
            {
                ContainedSeqPerThread[0].push_back(idxT);
            }
            else if (kind == FIRST_TO_SECOND_OVERLAP)
            {
                cks->dir  = rc? 0 : 1;
                cks->dirT = rc? 0 : 2;
                cks->sfx  = ((lenT - endTr) - (lenQ - endQr));
                cks->sfxT = begQr - begTr;
                cks->passed = true;
            }
            else
            {
                cks->dir  = rc? 3 : 2;
                cks->dirT = rc? 3 : 1;
                cks->sfx  = begTr - begQr;
                cks->sfxT = ((lenQ - endQr) - (lenT - endTr));
                cks->passed = true;
            }
        }
    }

    int readcount = 0;
    readcount += ContainedSeqPerThread[0].size();

    unsigned int readssofar = 0;
    ContainedSeqPerBatch.resize(readcount);

    // Concatenate per-thread result
    for(int t = 0; t < 1; ++t)
    {
        copy(ContainedSeqPerThread[t].begin(), ContainedSeqPerThread[t].end(), ContainedSeqPerBatch.begin() + readssofar);
        readssofar += ContainedSeqPerThread[t].size();
    }

    auto end_time = std::chrono::system_clock::now();
    add_time("XA:StringOp", (ms_t(end_time - start_time)).count());
}

void XDropAligner::add_time(std::string type, double duration) {
  // may be called threaded
  int tid = 0;
  #ifdef THREADED
  tid = omp_get_thread_num();
  #endif

  if (types[tid].find(type) != types[tid].end()){
    size_t idx = types[tid][type];
    ++counts[tid][idx];
    times[tid][idx] += duration;
  } else {
    types[tid][type] = times[tid].size();
    counts[tid].push_back(1);
    times[tid].push_back(duration);
  }
}

void XDropAligner::print_avg_times(std::shared_ptr<ParallelOps> parops, std::ofstream &lfs)
{
  // counts and times can be empty if there were non non-zeros in that
  // process grid cell. Therefore, we need to fix the collective call as follows.
  int counts_size = 0;
  for (auto el : counts)
	  counts_size = std::max(counts_size, static_cast<int>(el.size()));

  MPI_Allreduce(MPI_IN_PLACE, &counts_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // max over all threads
  std::vector<uint64_t> counts_tmp(counts_size, 0);
  for (auto &el : counts)
    for (int i = 0; i < el.size(); ++i)
      counts_tmp[i] = std::max(counts_tmp[i], el[i]);

  std::vector<double> times_tmp(counts_size, 0);
  for (auto &el : times)
    for (int i = 0; i < el.size(); ++i)
      times_tmp[i] = std::max(times_tmp[i], el[i]);

  std::unordered_map<std::string, size_t> *types_tmp = NULL;
  for (auto &el : types)
  {
    if (el.size() != 0)
    {
      types_tmp = &el;
      break;
    }
  }

  std::vector<double> maxtime(counts_size, 0);
  std::vector<double> mintime(counts_size, 0);
  std::vector<uint64_t> maxcount(counts_size, 0);
  std::vector<uint64_t> mincount(counts_size, 0);

  // min, max time
  MPI_Reduce(&(times_tmp[0]), &(maxtime[0]), times_tmp.size(), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&(times_tmp[0]), &(mintime[0]), times_tmp.size(), MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  // avg count, time
  MPI_Allreduce(MPI_IN_PLACE, &(counts_tmp[0]), counts_tmp.size(),
      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(times_tmp[0]), times_tmp.size(),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (parops->world_proc_rank == 0 && types_tmp != NULL)
  {
    for (auto &type : *types_tmp)
    {
      std::cout << "  PWF:" << type.first << " avg per proc: "
                << (times_tmp[type.second] / parops->world_procs_count)
                << " ms" << std::endl;

      std::cout << "  PWF:" << type.first << " max per proc: "
                << maxtime[type.second]
                << " ms" << std::endl;

      std::cout << "  PWF:" << type.first << " min per proc: "
                << mintime[type.second]
                << " ms" << std::endl;
    }
  }
}
