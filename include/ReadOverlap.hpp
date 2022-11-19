#ifndef READ_OVERLAPS_H_
#define READ_OVERLAPS_H_

#include "kmer/CommonKmers.hpp"
#include "DistributedFastaData.hpp"

using namespace elba;
using namespace combblas;

struct ReadOverlap
{
    int sfx, sfxT, dir, dirT;
    int b[2], e[2], l[2], coords[2], sfxpath[4];
    bool transpose, rc;
    int score;

    void SetPathInf();
    ReadOverlap();

    ReadOverlap(const ReadOverlap& rhs);
    ReadOverlap(const CommonKmers& cks);

    bool is_invalid() const;

    bool arrows(int& t, int& h) const;

    operator bool() const { return true; } /* For T = R in transitive reduction */

    friend bool operator==(const ReadOverlap& lhs, const ReadOverlap& rhs)
    {
        return (lhs.sfx == rhs.sfx && lhs.dir == rhs.dir);
    }

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap o = b;
        return o;
    }
};

int intplus(int a, int b);

struct Tupleize : std::unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& e)
    {
        ReadOverlap out = e;
        switch (e.dir) {
            case 0:
                out.coords[0] = e.b[0];
                out.coords[1] = e.l[1] - e.b[1];
                break;
            case 1:
                out.coords[0] = e.transpose? e.l[0] - e.e[0] : e.b[0];
                out.coords[1] = e.transpose? e.l[1] - e.e[1] : e.b[1];
                break;
            case 2:
                out.coords[0] = e.transpose? e.l[0] - e.b[0] : e.e[0];
                out.coords[1] = e.transpose? e.l[1] - e.b[1] : e.e[1];
                break;
            case 3:
                out.coords[0] = e.e[0];
                out.coords[1] = e.l[1] - e.e[1];
                break;
            default:
                break;
        }
        return out;
    }
};

struct ReadOverlapGraphHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.score << "\t" << e.l[0] << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << e.l[1] << "\t" << e.b[1] << "\t" << e.e[1] << "\t" << e.dir << "\t" << e.sfx;
    }
};


struct PafHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        char strand = e.rc? '-' : '+';
        os << (row+1) << "\t" << e.l[0] << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << strand << "\t" << (col+1) << "\t" << e.l[1] << "\t" << (e.rc? e.l[1]-e.e[1] : e.b[1]) << "\t" << (e.rc? e.l[1]-e.b[1] : e.e[1]);
    }
};

struct ReadOverlapDiskHandler
{
    /* The read overlap disk handler should only be called on an overlap graph
     * which has not been symmetricized yet. That is, the adjacency matrix is
     * upper triangular. */

    ReadOverlap getNoNum(int64_t row, int64_t col) { return ReadOverlap(); }

    template <typename c, typename t>
    ReadOverlap read(std::basic_istream<c,t>& is, int64_t row, int64_t col)
    {
        ReadOverlap e;
        e.transpose = false;

        int rc;
        is >> rc >> e.dir >> e.dirT >> e.sfx >> e.sfxT >> e.b[0] >> e.e[0] >> e.b[1] >> e.e[1] >> e.l[0] >> e.l[1];
        e.rc = static_cast<bool>(rc);

        e.SetPathInf();
        e.sfxpath[e.dir] = e.sfx;

        return e;
    }

    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << static_cast<int>(e.rc) << "\t" << e.dir << "\t" << e.dirT << "\t" << e.sfx << "\t" << e.sfxT << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << e.b[1] << "\t" << e.e[1] << "\t" << e.l[0] << "\t" << e.l[1];
    }

     void binaryfill(FILE *rfFile, int64_t& row, int64_t& col, ReadOverlap& e) {}
     size_t entrylength() { return sizeof(ReadOverlap); }
};


#endif
