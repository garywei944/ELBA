#include "../include/mpifa.h"
#include "../include/mstring.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h> /* isspace() */

/*
 * Make MPI_SIZE_T an alias to the correct MPI_Datatype for size_t.
 */
#ifndef MPI_SIZE_T
#if SIZE_MAX == ULONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#else
#error "size_t must be unsigned long"
#endif
#endif

#define mpifa_check_oom(ptr) do { \
    if ((ptr) == NULL) { \
        int _rank_id; \
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank_id); \
        if (_rank_id == 0) fprintf(stderr, "mpifa: out of memory\n"); \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    } \
} while (0)

/*
 * Get the next power-of-two greater than or equal to x.
 */
static inline size_t up_size_t(size_t x)
{
    if (x != 0)
        x--;

    x |= x>>1;
    x |= x>>2;
    x |= x>>4;
    x |= x>>8;
    x |= x>>16;
    x |= x>>32;

    return x+1;
}

/*
 * Function: `read_faidx_file`
 *
 * Parameters:
 *     -fai: faidx_t object that will be read into.
 *     -fname: FAIDX file name (extension is usually .fai)
 *     -comm: MPI communicator for all the processes which will
 *            be receiving a chunk of the FAIDX records.
 *
 * Returns: 0 on success, -1 if failure
 */
int read_faidx_file(faidx_t *fai, const char *fname, MPI_Comm comm)
{
    int nprocs;       /* number of processes in comm                           */
    int myrank;       /* my process id in comm                                 */
    int *sendcounts;  /* MPI_Scatterv sendcounts for FAIDX records (root only) */
    int *displs;      /* MPI_Scatterv displs for FAIDX records (root only)     */
    int recvcount;    /* MPI_Scatterv recvcount for FAIDX records              */
    int ierr;         /* MPI error code                                        */

    size_t pos;         /* byte position within faidx file buffer (root only)         */
    size_t remain;      /* number of bytes remaining in faidx file buffer (root only) */
    size_t avail_recs;  /* number of records available in memory (root only)          */
    size_t linesize;    /* current line size (root only)                              */

    char *buf;   /* faidx file buffer (root only)         */
    char *ptr;   /* faidx file buffer pointer (root only) */
    char *next;  /* address of next '\n' char (root only) */

    MPI_File faidx_fh;
    MPI_Offset filesize;

    faidx_rec_t *grecs, *myrecs, *rec; /* global and local faidx_rec_t arrays */
    size_t num_recs; /* number of local FAIDX records */
    size_t lwidth;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    int ndiffs;
    string_store_t names;

    if (myrank == 0)
    {
        names = STRING_STORE_INIT;
        /*
         * Root process slurps in the entire FAIDX file.
         */
        ierr = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &faidx_fh);

        if (ierr != MPI_SUCCESS)
        {
            fprintf(stderr, "error: MPI_File_open failed on FAIDX file named '%s'\n", fname);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        MPI_File_get_size(faidx_fh, &filesize);

        buf = malloc(filesize);
        mpifa_check_oom(buf);

        MPI_File_read(faidx_fh, buf, filesize, MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_close(&faidx_fh);

        /*
         * Parse each FAIDX record from the file buffer.
         */
        ptr = buf;
        remain = filesize;
        pos = num_recs = avail_recs = lwidth = 0;
        grecs = NULL;

        ndiffs = 0;
        size_t lwidth_last = 0;

        while (pos < filesize)
        {
            next = memchr(ptr, '\n', remain);

            if (next == NULL) /* TODO: this could be problematic if file does not have '\n' as its last character */
                break;

            *next = (char)0; /* TODO: otherwise sscanf might be a problem (not actually positive but there is no downside, test later */
            linesize = next - ptr;
            assert(linesize <= remain);

            if (num_recs + 1 > avail_recs)
            {
                avail_recs = up_size_t(num_recs + 1);
                avail_recs = avail_recs > 256? avail_recs : 256;
                grecs = realloc(grecs, avail_recs * sizeof(faidx_rec_t));
            }

            size_t lwidth2;
            rec = &grecs[num_recs++];
            sscanf(ptr, "%*s %zu %zu %zu %zu", &rec->readlen, &rec->fpos, &rec->lwidth, &lwidth2);
            assert(lwidth2-rec->lwidth==1);

            if (lwidth2 != lwidth_last)
                ndiffs += !!(ndiffs < 2);

            lwidth_last = lwidth2;

            size_t namelen;

            for (namelen = 0; namelen < linesize; ++namelen)
                if (isspace(ptr[namelen]))
                    break;

            sstore_push(&names, ptr, namelen);

            pos += linesize;
            ptr = &buf[++pos];
            remain = filesize - pos;
        }

        free(buf);
        grecs = realloc(grecs, num_recs * sizeof(faidx_rec_t));

        sendcounts = malloc(nprocs * sizeof(int));
        displs = malloc(nprocs * sizeof(int));
        displs[0] = 0;

        /*
         * Processes with ids 0..nprocs-2 each get floor(num_recs/nprocs) FAIDX
         * records each. The process with id nprocs-1 gets all the remaining
         * sequences at the end.
         */
        for (int i = 0; i < nprocs-1; ++i)
        {
            sendcounts[i] = num_recs / nprocs;
            displs[i+1] = displs[i] + sendcounts[i];
        }

        sendcounts[nprocs-1] = num_recs - (nprocs-1)*(num_recs/nprocs);
    }

    MPI_Datatype faidx_rec_mpi_t;
    MPI_Type_contiguous(3, MPI_SIZE_T, &faidx_rec_mpi_t);
    MPI_Type_commit(&faidx_rec_mpi_t);

    /*
     * Root process tells each process how many FAIDX records it will be sent.
     */
    MPI_Scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, comm);

    /*
     * Each process will receive its assigned portion of FAIDX records
     * at the location pointed to by myrecs.
     */
    myrecs = malloc(recvcount * sizeof(faidx_rec_t));
    MPI_Scatterv(grecs, sendcounts, displs, faidx_rec_mpi_t, myrecs, recvcount, faidx_rec_mpi_t, 0, comm);

    MPI_Type_free(&faidx_rec_mpi_t);

    if (myrank == 0)
    {
        free(grecs);
        free(sendcounts);
        free(displs);
    }

    MPI_Bcast(&ndiffs, 1, MPI_INT, 0, comm);

    fai->recs = myrecs;
    fai->numrecs = recvcount;
    fai->lwidth = ndiffs == 2? 0 : myrecs[0].lwidth;
    fai->comm = comm;

    sstore_mpi_scatter(&names, &fai->names, 0, comm);

    return 0;
}

size_t faidx_size(const faidx_t fai)
{
    size_t size = fai.numrecs;
    MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_SIZE_T, MPI_SUM, fai.comm);
    return size;
}

void log_faidx(const faidx_t fai, const char *log_fname)
{
    MPI_Comm comm = fai.comm;

    int nprocs, myrank;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    string_t buf = STRING_INIT;

    size_t offset;
    size_t numreads = fai.numrecs;
    MPI_Exscan(&numreads, &offset, 1, MPI_SIZE_T, MPI_SUM, comm);
    if (!myrank) offset = 0;

    for (size_t i = 0; i < numreads; ++i)
    {
        faidx_rec_t rec = fai.recs[i];
        string_catf(&buf, "%d.%ld\t%.*s\t%ld\t%ld\t%ld\n", myrank, i, sstore_get_string_length(fai.names, i), sstore_get_string(fai.names, i), rec.fpos, rec.readlen, rec.lwidth);
    }

    MPI_File f;
    MPI_File_open(comm, log_fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
    MPI_File_write_shared(f, buf.buf, (int)buf.len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&f);

    string_destroy(buf);
}
