#include "../include/PruneChimeras.hpp"
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <mpi.h>

using namespace elba;

PileupVector::PileupVector(int read_length) : pileup(read_length, 0) {}
PileupVector::PileupVector(const std::vector<int>& v, int offset, int size) : pileup(size, 0)
{
    for (int i = 0; i < size; ++i)
        pileup[i] = v[offset + i];
}

int PileupVector::Length() const { return pileup.size(); }

void PileupVector::AddInterval(int begin, int end)
{
    assert((begin >= 0 && end < Length()));
    for (int i = begin; i < end; ++i) pileup[i]++;
}

int PackPileupVectors(const std::vector<PileupVector>& pvs, std::vector<int>& packed, std::vector<int>& lens)
{
    int pv_size = pvs.size();
    int packed_size = 0;
    packed.clear();
    lens.reserve(pv_size);

    for (int i = 0; i < pv_size; ++i)
    {
        int len = pvs[i].Length();
        lens.push_back(len);
        packed_size += len;
    }

    packed.reserve(packed_size);

    for (int i = 0; i < pv_size; ++i)
        for (int j = 0; j < lens[i]; ++j)
            packed.push_back(pvs[i].pileup[j]);

    return packed_size;
}

void UnpackPileupVectors(const std::vector<int>& packed, const std::vector<int>& lens, std::vector<PileupVector>& pvs)
{
    int size = lens.size();
    pvs.clear();
    pvs.reserve(size);
    int offset = 0;

    for (int i = 0; i < size; ++i)
    {
        pvs.emplace_back(packed, offset, lens[i]);
        offset += lens[i];
    }
}

void add_gaps(int begin, int end, int length, std::vector<std::tuple<int, int>>& middle, std::vector<std::tuple<int, int>>& extremity)
{
    if (begin != end)
    {
        if (!begin || end == length) extremity.push_back({begin, end});
        else middle.push_back({begin, end});
    }
}

void ReportChimeras(const std::shared_ptr<CommGrid>& commgrid, std::shared_ptr<DistributedFastaData> dfd, const std::vector<PileupVector>& pileups)
{
    static const int coverage_min = 0;

    uint64_t col_seq_start_idx = dfd->col_seq_start_idx;

    int myrank;
    MPI_Comm_rank(commgrid->GetWorld(), &myrank);

    for (int i = 0; i < pileups.size(); ++i)
    {
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
    int myrank = commgrid->GetRank();

    auto spSeq = Rmat->seq();

    /* Want to calculate pileup vectors for each read. Do this by first
     * computing the pileup vectors for reads in the column range of the local
     * processor, by going through each local column and computing the pileup
     * as a function of the overlapping intervals of the column read for each
     * nonzero in that locally stored column */
    std::vector<PileupVector> local_pileups;
    uint64_t local_ncols = spSeq.getncol();
    local_pileups.reserve(local_ncols);

    for (uint64_t i = 0; i < local_ncols; ++i)
    {
        int read_length = length(*(dfd->col_seq(i)));
        local_pileups.emplace_back(read_length);
    }

    std::stringstream stream_name;
    stream_name << "elba_rank_" << myrank << ".paf";
    std::ofstream paf_stream(stream_name.str());

    std::stringstream line;

    //if (myrank == 0)
    //{
    //    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
    //    {
    //        int colid = colit.colid();

    //        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
    //        {
    //            std::cout << "AA\t" << rowid << "\t" << colid << "|t" 
    //        }
    //    }
    //}

    /* iterate over every local column */
    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
    {
        int colid = colit.colid();

        /* iterate over every row that has an overlap with current column */
        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
        {
            ReadOverlap o = nzit.value();
            int rowid = nzit.rowid();
            std::cout << "myrank=" << Rmat->getcommgrid()->GetRank() << "; 1 :: local_pileups[" << colid+1 << "].AddInterval(" << o.b[1] << ", " << o.e[1] << ") (" << rowid+1 << ") " << std::endl;
            local_pileups[colid].AddInterval(o.b[1], o.e[1]);
            std::cout << "myrank=" << Rmat->getcommgrid()->GetRank() << "; 2 :: local_pileups[" << colid+1 << "].AddInterval(" << o.b[1] << ", " << o.e[1] << ") (" << rowid+1 << ") " << std::endl;
        }
    }


    paf_stream.close();

    std::vector<int> lens;
    std::vector<int> packed;
    int size = PackPileupVectors(local_pileups, packed, lens);

    MPI_Allreduce(MPI_IN_PLACE, packed.data(), size, MPI_INT, MPI_SUM, commgrid->GetColWorld());

    UnpackPileupVectors(packed, lens, local_pileups);

    return local_pileups;
}
