#ifndef PRUNE_CHIMERAS_H_
#define PRUNE_CHIMERAS_H_

#include <cassert>
#include <vector>
#include <tuple>

using namespace elba;

class PileupVector
{
public:
    std::vector<int> pileup;

    PileupVector(int read_length) : pileup(read_length, 0) {}
    PileupVector(const std::vector<int>& v, int offset, int size) : pileup(size, 0)
    {
        for (int i = 0; i < size; ++i)
            pileup[i] = v[offset + i];
    }

    int Length() const { return pileup.size(); }

    void AddInterval(int begin, int end)
    {
        assert((begin > 0 && end < Length()));
        for (int i = begin; i < end; ++i) pileup[i]++;
    }

    ~PileupVector() = default;
};

std::vector<int> PackPileupVector(const std::vector<PileupVector>& pvs, std::vector<int>& lens, int& size)
{
    size = 0;
    lens.clear();
    std::vector<int> packed;
    for (const PileupVector& pv : pvs)
    {
        int len = pv.Length();
        lens.push_back(len);
        size += len;
        for (int i = 0; i < len; ++i)
            packed.push_back(pv.pileup[i]);

    }

    return packed;
}

std::vector<PileupVector> UnpackPileupVector(const std::vector<int>& packed, const std::vector<int>& lens)
{
    std::vector<PileupVector> pvs;

    int offset = 0;
    for (int i = 0; i < lens.size(); ++i)
    {
        pvs.emplace_back(packed, offset, lens[i]);
        offset += lens[i];
    }

    return pvs;
}

void add_gaps(int begin, int end, int length, std::vector<std::tuple<int, int>>& middle, std::vector<std::tuple<int, int>>& extremity)
{
    if (begin != end)
    {
        if (!begin || end == length) extremity.push_back({begin, end});
        else middle.push_back({begin, end});
    }
}

void ReportChimeras(std::shared_ptr<DistributedFastaData> dfd, const std::vector<PileupVector>& pileups)
{
    static const int coverage_min = 0;

    uint64_t col_seq_start_idx = dfd->col_seq_start_idx;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for (int i = 0; i < pileups.size(); ++i)
    {
        std::cout << myrank << ":: " << i << std::endl;
        const PileupVector& pv = pileups[i];

        std::vector<std::tuple<int, int>> middle, extremity;

        bool in_gap = true;
        int begin = 0, end = 0;

        for (int j = 0; j < pv.Length(); ++j)
        {
            if (pv.pileup[j] <= coverage_min && !in_gap)
            {
                end = 0, begin = j;
                in_gap = true;
            }

            if (pv.pileup[j] > coverage_min && in_gap)
            {
                end = j;
                in_gap = false;
                add_gaps(begin, end, pv.Length(), middle, extremity);
            }
        }

        if (in_gap)
        {
            end = pv.Length();
            add_gaps(begin, end, pv.Length(), middle, extremity);
        }

        if (middle.size() > 0)
        {
            std::cout << myrank << ":: chimeric: " << (i+col_seq_start_idx+1) << std::endl;
        }
    }
}

std::vector<PileupVector> GetReadPileup(std::shared_ptr<DistributedFastaData> dfd, PSpMat<ReadOverlap>::MPI_DCCols*& Rmat, const std::shared_ptr<ParallelOps>& parops)
{
    auto commgrid = Rmat->getcommgrid();

    auto spSeq = Rmat->seq();

    std::vector<PileupVector> local_pileups;
    uint64_t local_ncols = spSeq.getncol();

    for (uint64_t i = 0; i < local_ncols; ++i)
    {
        int read_length = length(*(dfd->col_seq(i)));
        local_pileups.emplace_back(read_length);
    }

    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
    {
        int colid = colit.colid();

        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
        {
            int rowid = nzit.rowid();
            ReadOverlap o = nzit.value();
            if (o.rc)
            {
                std::cout << rowid+dfd->row_seq_start_idx+1 << "\t" << o.l[0] << "\t" << o.b[0] << "\t" << o.e[0] << "\t-\t" <<
                             colid+dfd->col_seq_start_idx+1 << "\t" << o.l[1] << "\t" << (o.l[1] - o.e[1]) << "\t" << (o.l[1] - o.b[1]) << std::endl;
            }
            else
            {
                std::cout << rowid+dfd->row_seq_start_idx+1 << "\t" << o.l[0] << "\t" << o.b[0] << "\t" << o.e[0] << "\t+\t" <<
                             colid+dfd->col_seq_start_idx+1 << "\t" << o.l[1] << "\t" << o.b[1] << "\t" << o.e[1] << std::endl;
            }

            if (!o.rc) local_pileups[colid].AddInterval(o.b[1], o.e[1]);
            else local_pileups[colid].AddInterval(o.l[1] - o.e[1], o.l[1] - o.b[1]);
        }
    }

    int size;
    std::vector<int> lens;
    std::vector<int> packed = PackPileupVector(local_pileups, lens, size);

    assert(static_cast<int>(packed.size()) == size);

    MPI_Allreduce(MPI_IN_PLACE, packed.data(), size, MPI_INT, MPI_SUM, commgrid->GetColWorld());

    return UnpackPileupVector(packed, lens);
}


#endif
