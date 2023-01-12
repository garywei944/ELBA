#ifndef MPIFA_H_
#define MPIFA_H_

#include <mpi.h>
#include <limits.h> /* ULONG_MAX */
#include "mstring.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_faidx_rec faidx_rec_t;

/*
 * FAIDX structure that stores an array of FAIDX records as well
 * as the line width used across all the records. Two kinds of line
 * widths are supported: 0 indicates that the sequences are each stored
 * entirely on their own line. Any other value greater than 0 indicates
 * that the sequences have newline characters inserted every @lwidth
 * characters.
 *
 * Warning: fai index files technically support the possibility that each
 *          individual sequence can have its own line width definition. But
 *          this is unnecessary and wastes having to pass around an array
 *          of line widths for each sequence. Therefore that case is unsupported
 *          and any FASTA file like that will cause undefined behaviour (likely a crash).
 */
typedef struct s_faidx
{
    string_store_t names;  /* local storage of sequence names */
    faidx_rec_t *recs;       /* local number of array of faidx records */
    size_t numrecs;        /* local number of records */
    size_t lwidth;         /* global FASTA line width */
    MPI_Comm comm;         /* MPI communicator */
} faidx_t;

/*
 * FAIDX record structure that identifies a read's length and its
 * offset within the given FASTA file.
 */
typedef struct s_faidx_rec
{
    size_t readlen; /* length of indexed sequence */
    size_t fpos;    /* position of read within FASTA */
    size_t lwidth;  /* line width */
} faidx_rec_t;

/*
 * Read a FAIDX file and distribute to all processes in the
 * MPI communicator.
 */
int read_faidx_file(faidx_t *fai, const char *fname, MPI_Comm comm);

/*
 * Global number of FAIDX records across all processes.
 */
size_t faidx_size(const faidx_t fai);

#ifdef __cplusplus
}
#endif

#endif
