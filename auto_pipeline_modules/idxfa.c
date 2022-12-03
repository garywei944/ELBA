#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>

char progname[128];
#define BUFLEN (16384)

int usage()
{
    fprintf(stderr, "\nUsage: %s <reads.fa>\n\n", progname);
    return -1;
}

struct string_t { char *s; size_t l, m; };
struct stream_t { unsigned char *buf; int beg, end, is_eof, fd; };

static inline size_t rndup32(size_t n)
{
    n--;
    for (int i = 1; i <= 16; i <<= 1) n |= (n<<1);
    return ++n;
}

#define unlikely(x) __builtin_expect((x),0)

static inline int stream_getc(struct stream_t *fs)
{
    if (unlikely(!!fs->is_eof && fs->beg >= fs->end)) return -1;
    if (unlikely(fs->beg >= fs->end))
    {
        fs->beg = 0, fs->end = read(fs->fd, fs->buf, BUFLEN);
        if (fs->end < BUFLEN) fs->is_eof = 1;
        if (!fs->end) return -1;
    }
    return fs->buf[fs->beg++];
}

static inline int stream_getline(struct string_t *s, struct stream_t *fs, int append)
{
    if (!append)
    {
        s->l = 0;
        s->s[0] = 0;
    }

    int c;

    while ((c = stream_getc(fs)) != -1)
    {
        if (unlikely(s->m < s->l+1))
        {
            s->m = rndup32(s->m);
            s->s = realloc(s->s, s->m);
        }

        if (unlikely(c == '\n'))
        {
            s->s[s->l] = 0;
            return 1;
        }
        else
        {
            s->s[s->l++] = c;
        }
    }

    s->s[s->l] = 0;
    return 0;
}

int main(int argc, char *argv[])
{
    strncpy(progname, argv[0], 128);

    if (argc < 1 || argc > 2) return usage();

    int fd = (argc == 1)? fileno(stdin) : open(argv[1], O_RDONLY);

    if (fd < 0)
    {
        fprintf(stderr,"error: unable to open '%s' [terminating]\n", (argc==1)? "stdin" : argv[1]);
        return usage();
    }

    struct stream_t fs;
    struct string_t line;

    fs.buf = malloc(BUFLEN);
    fs.beg = fs.end = fs.is_eof = 0;
    fs.fd = fd;

    line.l = 0;
    line.m = 16;
    line.s = malloc(line.m);
    line.s[0] = 0;

    char name[1024];
    int i, j, l;

    i = l = 0;
    name[0] = 0;

    while (stream_getline(&line, &fs, 0))
    {
        if (line.s[0] == '>')
        {
            if (l > 0)
            {
                printf("%d\t%d\t%s\n", ++i, l, name);
                l = 0;
            }


            for (j = 1; j < 1024 && j < line.l; ++j)
                if (!isspace(line.s[j]))
                    name[j-1] = line.s[j];
                else break;

            name[j-1] = 0;
        }
        else l += line.l;
    }

    if (l > 0) printf("%d\t%d\t%s\n", ++i, l, name);

    free(line.s);
    free(fs.buf);

    if (fs.fd != fileno(stdin)) close(fs.fd);

    return 0;
}
