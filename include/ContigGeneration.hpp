#ifndef CONTIG_GENERATION_HPP_
#define CONTIG_GENERATION_HPP_

#include <map>
#include <cassert>

#include "TraceUtils.hpp"
#include "ReadOverlap.hpp"
#include "Utils.hpp"
#include "FastaData.hpp"
#include "DistributedFastaData.hpp"

using namespace combblas;

typedef int64_t IType; /* index type used in this file */

struct DistReadInfo
{
    std::shared_ptr<CommGrid> commgrid;
    MPI_Comm world;
    int myrank, nprocs;

    FastaData *lfd;
    std::vector<IType> offsets, numreads;
    IType cached_left_val, cached_right_val;
    int cached_index;

    DistReadInfo(std::shared_ptr<CommGrid> commgrid, FastaData *lfd);
    int GetReadOwner(IType gidx);
};


IType GetRead2Contigs
(
    SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G,
    FullyDistVec<IType,IType>& Read2Contigs,
    DistReadInfo& di,
    TraceUtils& tu
);

FullyDistVec<IType,IType> GetContigSizes
(
    const FullyDistVec<IType,IType>& Read2Contigs,
    const IType NumContigs,
    DistReadInfo& d
);

/* @func GetLocalRead2Procs      determine the local read to processor assignments for load
 *                               balancing contig generation.
 *
 * @param Read2Contigs           distributed read to contig id assignments.
 * @param AllContigSizesSorted   global vector of (contig id, size) tuples, sorted by size.
 * @param
 *
 * @description
 * Using the @Read2Contigs assignments and globally shared @AllContigSizesSorted vector of
 * used contig ids (and sizes), we compute the partitioning of used contigs to processor ranks using
 * greedy multiway number partitioning optimization algorithm. We call 'small'-indices the
 * indices of used contigs within the smaller @AllContigSizesSorted vector, and we call
 * 'large'-indices the corresponding global contig ids originally calculated by connected
 * components. The small indices are used for the optimization algorithm since they enable
 * us a headache-free way to run this algorithm and to broadcast the partitioning to
 * every processor (note: the optimization algorithm is only run on the root rank).
 *
 * Once the optimization algorithm is finished on the root rank, the used contig-to-processor
 * assignments are broadcast to everyone, and the mapping of the small indices to the large
 * indices is broadcast as well.
 *
 * We then need to create a hashmap which inverts the small-to-large map and gives us a
 * large-to-small map, so that we know which global contig-id we are referring to when
 * we are looking at a small index.
 *
 * At this point, we want to finally create a local vector for each processor which tells
 * us where the reads that our processor owns are supposed to end up. The idea is that since
 * we know which reads belong to which contigs, and we now know which contigs are to be sent
 * to which processors, we automatically know which reads are to be sent to which processors.
 * So we simply go through our local read-to-contig assignments and, using the globally computed
 * contig partititioning, determine where our local reads are supposed to go.
 *
 * @return local read to processor assignments */
std::vector<IType> GetLocalRead2Procs
(
    FullyDistVec<IType,IType>& Read2Contigs,
    std::vector<std::tuple<IType,IType>>& AllContigSizesSorted,
    const IType NumUsedContigs,
    DistReadInfo& di,
    TraceUtils& t
);

/* @func ReadExchange          communicate reads that each processor doesn't have access to,
                               but which needs in order to generate contigs.
 *
 * @param LocalRead2Procs      local read-to-processor assignments.
 * @param charbuf_info [ref]   hashmap of global read-indices to (charbuffer-offset, read-length) tuples.
 *
 * @description
 * Plan is to send a packed char buffer of reads to each processor and to send a
 * map @charbuf_info which has the associated offset within the char buffer and the
 * length of the read, so that read sequences can be quickly looked up within the
 * char buffer. We do this by communicating two objects: the char buffer,
 * and an array of info-tuples ordered according to how the char buffer itself is ordered.
 *
 * First step is to loop through @LocalRead2Procs vector for reads which need
 * to be sent to a different processor. For each such read, we collect the information
 * regarding where it is to be sent, where its offset is within the @lfd_buffer, and
 * how long it is. We also calculate the char buffer sendcounts in the same loop.
 * We then pack up the tuple and counts information and send them around.
 *
 * Then we allocate the char send and receive buffers because we now know how long
 * they should each be, and then we fill the char send buffer. We then send
 * the char buffers around.
 *
 * Finally, we use the tuple buffers to construct the local @charbuf_info map on
 * each process, and then we are done.
 *
 * @return char buffer pointer
 */
const char * ReadExchange
(
    std::vector<IType>& LocalRead2Procs,
    std::unordered_map<IType, std::tuple<IType, ushort>>& charbuf_info,
    DistReadInfo& di,
    TraceUtils& tu
);

/* @func LocalAssembly           Assemble contigs from locally induced subgraph and read sequence information.
 *
 * @param ContigChains           locally induced subgraph of contig chains.
 * @param LocalContigReadIdxs    mapping of local graph indices to original global read indices.
 * @param charbuf                char buffer of received reads (not including reads already here).
 * @param charbuf_info           information for char buffer lookup/reading.
 *
 * @description
 *
 * @returns vector of completed contigs */
std::vector<std::string> LocalAssembly
(
    SpCCols<IType,ReadOverlap>& ContigChains,
    std::vector<IType>& LocalContigReadIdxs,
    const char* charbuf,
    std::unordered_map<IType, std::tuple<IType, ushort>> charbuf_info,
    DistReadInfo& di
);

void AppendContig
(
    std::string& contig,
    const char *buf,
    IType offset,
    ushort len,
    IType start,
    IType end
);

int MPI_Alltoallv_str
(
    const char *sendbuf,
    const std::vector<IType>& sendcounts,
    const std::vector<IType>& sdispls,
    char *recvbuf,
    const std::vector<IType>& recvcounts,
    const std::vector<IType>& rdispls,
    MPI_Comm comm
);

/* @func CreateContig   Assemble contigs from distributed string graph.
 *
 * @param G             combblas distributed string graph.
 * @param dfd           distributed fasta data object.
 *
 * @return vector<string> of contigs.
 */

std::vector<std::string> CreateContig
(
    SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G,
    std::shared_ptr<DistributedFastaData> dfd,
    std::string& myoutput,
    const std::shared_ptr<TimePod>& tp,
    TraceUtils& tu
);


template <class IT, class NT, class DER>
FullyDistVec<IT,IT> LastNzRowIdxPerCol(const SpParMat<IT,NT,DER>& A);

IType RemoveBridgeVertices(SpParMat<IType,IType,SpDCCols<IType,IType>>& A);

IType KTipsRemoval
(
    SpParMat<IType,IType,SpDCCols<IType,IType>>& A,
    const FullyDistVec<IType,IType>& degrees,
    const IType l,
    TraceUtils& tu
);

/* @func GetRead2Contigs       determines which reads belong to which which contigs.
 *
 * @param G                    distributed string graph.
 * @Param Read2Contigs   [ref] computed assignments vector.
 *
 * @description
 * Contig branching points are found by CombBLAS row reduction (to find vertex degrees),
 * and then selecting rows with more than 2 nonzeros. These are temporarily zeroed from
 * the matrix so that connected components can be run on the graph with maximum degree of 2,
 * which gives us a distributed vector assigning each vertex to a contig.
 *
 * @return number of contigs */
IType GetRead2Contigs
(
    SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G,
    FullyDistVec<IType,IType>& Read2Contigs,
    DistReadInfo& di,
    TraceUtils& tu
);

/* @func GetContigSizes     calculates the number of reads within each contig.
 *
 * @param Read2Contigs      read-to-contig assignments.
 * @param NumContigs        total number of contigs.
 *
 * @return distributed vector of contig sizes */
FullyDistVec<IType,IType> GetContigSizes
(
    const FullyDistVec<IType,IType>& Read2Contigs,
    const IType NumContigs,
    DistReadInfo& di
);

/* @func GetAllContigSizesSorted
 *
 * @param ContigSizes
 * @param NumUsedContigs
 * @param minsize
 *
 * @description
 * Using the distributed vector of contig sizes as input, find all contigs with
 * >= @minsize reads, and All-gather them so that each processor has a vector
 * of (contig-id (global), contig-size) tuples. Then each processor sorts their
 * vector by contig-size and returns the result.
 *
 * @return */
std::vector<std::tuple<IType,IType>> GetAllContigSizesSorted
(
    FullyDistVec<IType,IType>& ContigSizes,
    IType& NumUsedContigs,
    IType minsize,
    DistReadInfo& di,
    TraceUtils& tu
);

/* @func ImposeMyReadDistribution
 *
 * @param assignments combblas distributed assignments vector
 *
 * @description
 * Takes a combblas distributed assignments vector and returns the local section
 * of the vector redistributed according to the read distribution of ELBA.
 *
 * @return redistributed local section of assignments vector.
 */
std::vector<IType> ImposeMyReadDistribution
(
    FullyDistVec<IType,IType>& assignments,
    DistReadInfo& di
);

void PruneTips
(
    SpParMat<IType,IType,SpDCCols<IType,IType>>& A,
    int ktip,
    TraceUtils& tu
);

#endif
