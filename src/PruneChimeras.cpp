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
    assert((begin > 0 && end < Length()));
    for (int i = begin; i < end; ++i) pileup[i]++;
}

int PackPileupVector(const std::vector<PileupVector>& pvs, std::vector<int>& packed, std::vector<int>& lens)
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

    return pv_size;
}

// std::vector<int> PackPileupVector(const std::vector<PileupVector>& pvs, std::vector<int>& lens, int& size)
// {
    // size = 0;
    // lens.clear();
    // std::vector<int> packed;
    // for (int i = 0; i < pvs.size(); ++i)
    // {
        // int len = pvs[i].Length();
        // lens.push_back(len);
        // size += len;
        // for (int j = 0; j < len; ++j)
            // packed.push_back(pvs[i].pileup[j]);
    // }

    // return packed;
// }

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

// std::vector<PileupVector> AddPileups(const std::vector<PileupVector>& pv1, const std::vector<PileupVector>& pv2)
// {
    // int len = pv1.size();
    // assert(len == static_cast<int>(pv2.size()));

    // std::vector<PileupVector> pv;

    // for (int k = 0; k < len; ++k)
    // {
        // int plen = pv1[k].Length();
        // assert(plen == pv2[k].Length());

        // pv.emplace_back(plen);

        // for (int i = 0; i < plen; ++i)
        // {
            // pv[k].pileup[i] = pv1[k].pileup[i] + pv2[k].pileup[i];
        // }
    // }

    // return pv;
// }


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

//std::vector<PileupVector> GetReadPileup2(std::shared_ptr<DistributedFastaData> dfd, PSpMat<ReadOverlap>::MPI_DCCols*& Rmat, const std::shared_ptr<ParallelOps>& parops, std::vector<PileupVector>& row_pvs, std::vector<PileupVector>& col_pvs)
//{
//    auto commgrid = Rmat->getcommgrid();
//    int myrank = commgrid->GetRank();
//    auto spSeq = Rmat->seq();
//
//    std::vector<PileupVector> row_pileups, col_pileups;
//    uint64_t ncols = spSeq.getncol();
//    uint64_t nrows = spSeq.getnrow();
//
//    for (uint64_t i = 0; i < nrows; ++i)
//    {
//        int read_length = length(*dfd->row_seq(i));
//        row_pileups.emplace_back(read_length);
//    }
//
//    for (uint64_t j = 0; j < ncols; ++j)
//    {
//        int read_length = length(*dfd->col_seq(j));
//        col_pileups.emplace_back(read_length);
//    }
//
//    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
//    {
//        int colid = colit.colid();
//
//        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
//        {
//
//            int rowid = nzit.rowid();
//            ReadOverlap o = nzit.value();
//
//            row_pileups[rowid].AddInterval(o.b[0], o.e[0]);
//
//            if (!o.rc) col_pileups[colid].AddInterval(o.b[1], o.e[1]);
//            else col_pileups[colid].AddInterval(o.l[1] - o.e[1], o.l[1] - o.b[1]);
//        }
//    }
//
//    int rowsize, colsize;
//    std::vector<int> rowlens, collens;
//    std::vector<int> rowpacked = PackPileupVector(row_pileups, rowlens, rowsize);
//    std::vector<int> colpacked = PackPileupVector(col_pileups, collens, colsize);
//
//    MPI_Allreduce(MPI_IN_PLACE, rowpacked.data(), rowsize, MPI_INT, MPI_SUM, commgrid->GetRowWorld());
//    MPI_Allreduce(MPI_IN_PLACE, colpacked.data(), colsize, MPI_INT, MPI_SUM, commgrid->GetColWorld());
//
//    // row_pvs = UnpackPileupVector(rowpacked, rowlens);
//    // col_pvs = UnpackPileupVector(colpacked, collens);
//
//    // if (commgrid->GetRankInProcRow() == commgrid->GetRankInProcCol())
//    // {
//        // pvs = row_pvs + col_pvs;
//        // broadcast()
//    // }
//
//    //int size;
//    //std::vector<int> lens;
//    //std::vector<int> packed = PackPileupVector(local_pileups, lens, size);
//
//    //assert(static_cast<int>(packed.size()) == size);
//
//    //MPI_Allreduce(MPI_IN_PLACE, packed.data(), size, MPI_INT, MPI_SUM, commgrid->GetColWorld());
//
//    //return UnpackPileupVector(packed, lens);
//
//}

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

    for (uint64_t i = 0; i < local_ncols; ++i)
    {
        int read_length = length(*(dfd->col_seq(i)));
        local_pileups.emplace_back(read_length);
    }

    std::stringstream stream_name;
    stream_name << "elba_rank_" << myrank << ".paf";
    std::ofstream paf_stream(stream_name.str());

    /* iterate over every local column */
    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
    {
        int colid = colit.colid();

        /* iterate over every row that has an overlap with current column */
        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
        {
            int rowid = nzit.rowid();
            ReadOverlap o = nzit.value();
            if (o.rc)
            {
                paf_stream << rowid+dfd->row_seq_start_idx+1 << "\t" << o.l[0] << "\t" << o.b[0] << "\t" << o.e[0] << "\t-\t" <<
                              colid+dfd->col_seq_start_idx+1 << "\t" << o.l[1] << "\t" << (o.l[1] - o.e[1]) << "\t" << (o.l[1] - o.b[1]) << std::endl;
            }
            else
            {
                paf_stream << rowid+dfd->row_seq_start_idx+1 << "\t" << o.l[0] << "\t" << o.b[0] << "\t" << o.e[0] << "\t+\t" <<
                              colid+dfd->col_seq_start_idx+1 << "\t" << o.l[1] << "\t" << o.b[1] << "\t" << o.e[1] << std::endl;
            }

            /* append overlap interval of column */
            if (!o.rc) local_pileups[colid].AddInterval(o.b[1], o.e[1]);
            else local_pileups[colid].AddInterval(o.l[1] - o.e[1], o.l[1] - o.b[1]);
        }
    }

    paf_stream.close();

    /* In order to get the correct overall pileup vector, we need to sum for each column the  */

    return local_pileups;
    //int size;
    //std::vector<int> lens;
    //std::vector<int> packed = PackPileupVector(local_pileups, lens, size);

    //assert(static_cast<int>(packed.size()) == size);

    //MPI_Allreduce(MPI_IN_PLACE, packed.data(), size, MPI_INT, MPI_SUM, commgrid->GetColWorld());

    //return UnpackPileupVector(packed, lens);
}
