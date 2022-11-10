#include "../include/ReadOverlap.hpp"
#include <limits>
#include <cassert>
#include <iostream>

static constexpr int MAX_INT = std::numeric_limits<int>::max();

int intplus(int a, int b)
{
    return (a == MAX_INT || b == MAX_INT) ? MAX_INT : a + b;
}

void ReadOverlap::SetPathInf() { sfxpath[0] = sfxpath[1] = sfxpath[2] = sfxpath[3] = MAX_INT; }

ReadOverlap::ReadOverlap() : sfx(0), dir(-1), score(-1), transpose(false) { SetPathInf(); }

ReadOverlap::ReadOverlap(const ReadOverlap& rhs) : sfx(rhs.sfx), sfxT(rhs.sfxT), dir(rhs.dir), dirT(rhs.dirT), transpose(rhs.transpose), rc(rhs.rc), score(rhs.score)
{
    b[0] = rhs.b[0]; b[1] = rhs.b[1];
    e[0] = rhs.e[0]; e[1] = rhs.e[1];
    l[0] = rhs.l[0]; l[1] = rhs.l[1];

    coords[0] = rhs.coords[0];
    coords[1] = rhs.coords[1];

    for (int i = 0; i < 4; ++i)
        sfxpath[i] = rhs.sfxpath[i];
}

ReadOverlap::ReadOverlap(const CommonKmers& cks) : transpose(false), score(static_cast<int>(cks.score))
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

bool ReadOverlap::is_invalid() const { return (dir == -1); }

bool ReadOverlap::arrows(int& t, int& h) const
{
    if (is_invalid())
        return false;

    t = (dir>>1)&1;
    h = dir&1;

    return true;
}
