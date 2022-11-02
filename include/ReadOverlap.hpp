#ifndef READ_OVERLAPS_H_
#define READ_OVERLAPS_H_

#include <limits>
#include <cassert>
#include <iostream>
#include "kmer/CommonKmers.hpp"

using namespace elba;

static constexpr int MAX_INT = std::numeric_limits<int>::max();
extern int xdrop;

struct ReadOverlap
{
    int sfx, sfxT, dir, dirT;
    int b[2], e[2], l[2], coords[2], sfxpath[4];
    bool transpose, rc;
    int score;

    void SetPathInf() { sfxpath[0] = sfxpath[1] = sfxpath[2] = sfxpath[3] = MAX_INT; }

    ReadOverlap() : sfx(0), dir(-1), score(-1), transpose(false)
    {
        SetPathInf();
    }

    ReadOverlap(const ReadOverlap& rhs)
        : sfx(rhs.sfx), sfxT(rhs.sfxT), dir(rhs.dir), dirT(rhs.dirT), transpose(rhs.transpose), rc(rhs.rc), score(rhs.score)
    {
        b[0] = rhs.b[0]; b[1] = rhs.b[1];
        e[0] = rhs.e[0]; e[1] = rhs.e[1];
        l[0] = rhs.l[0]; l[1] = rhs.l[1];

        coords[0] = rhs.coords[0];
        coords[1] = rhs.coords[1];

        for (int i = 0; i < 4; ++i)
            sfxpath[i] = rhs.sfxpath[i]; 
    }

    ReadOverlap(const CommonKmers& cks) : transpose(false), score(static_cast<int>(cks.score))
    {
        b[0] = cks.first.first;  b[1] = cks.second.first;
        e[0] = cks.first.second; e[1] = cks.second.second;
        l[0] = cks.lenv;         l[1] = cks.lenh;

        rc = cks.rc;

        sfx  = cks.sfx;
        sfxT = cks.sfxT;
        dir  = cks.dir;
        dirT = cks.dirT;

        SetPathInf();

        sfxpath[dir] = sfx;
    }

    bool is_invalid() const { return (dir == -1); }

    bool arrows(int& t, int& h) const
    {
        if (is_invalid())
            return false;

        t = (dir>>1)&1;
        h = dir&1;

        return true;
    }

    operator bool() const { return true; } /* for T = R in transitive reduction */

    friend bool operator==(const ReadOverlap& lhs, const ReadOverlap& rhs)
    {
        return (lhs.sfx == rhs.sfx && lhs.dir == rhs.dir);
    }

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap myobj = b;
        return myobj;
    }
};

int intplus(int a, int b)
{
    return (a == MAX_INT || b == MAX_INT) ? MAX_INT : a + b;
}

struct Tupleize : unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& e)
    {
        ReadOverlap out = e;
        switch (e.dir) {
            case 0:
                out.coords[0] = e.b[0] + xdrop;
                out.coords[1] = e.l[1] - e.b[1];
                break;
            case 1:
                out.coords[0] = (e.transpose)? (e.l[0] - e.e[0] + xdrop) : (e.b[0] + xdrop);
                out.coords[1] = (e.transpose)? (e.l[1] - e.e[1]) : (e.b[1]);
                break;
            case 2:
                out.coords[0] = (e.transpose)? (e.l[0] - e.b[0] - xdrop) : (e.e[0] - xdrop);
                out.coords[1] = (e.transpose)? (e.l[1] - e.b[1]) : (e.e[1]);
                break;
            case 3:
                out.coords[0] = e.e[0] - xdrop;
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
