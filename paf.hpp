#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct
{
    void *fp;
    kstring_t buf;
} paf_file_t;

typedef struct
{
    const char *qn, *tn;             // these point to the input string; NOT allocated
    uint32_t ql, qs, qe, tl, ts, te; // query len/start/end, target len/start/end
    uint32_t ml;                     // number of residue matches
    uint32_t rev;                    // strand
    uint32_t bl;                     // alignment block length
} paf_rec_t;

paf_file_t *paf_open(const char *fn)
{
    kstream_t *ks;
    gzFile fp;
    paf_file_t *pf;
    fp = gzopen(fn, "r");
    if (fp == 0)
        return 0;
    ks = ks_init(fp);
    pf = (paf_file_t *)calloc(1, sizeof(paf_file_t));
    pf->fp = ks;
    return pf;
}
int paf_close(paf_file_t *pf)
{
    kstream_t *ks;
    if (pf == 0)
        return 0;
    free(pf->buf.s);
    ks = (kstream_t *)pf->fp;
    gzclose(ks->f);
    ks_destroy(ks);
    free(pf);
    return 0;
}

int paf_parse(int l, char *s, paf_rec_t *pr) // s must be NULL terminated
{
    // on return: <0 for failure; 0 for success; >0 for filtered
    char *q, *r;
    int i, t;
    for (i = t = 0, q = s; i <= l; ++i)
    {
        if (i < l && s[i] != '\t')
            continue;
        s[i] = 0;
        if (t == 0)
            pr->qn = q;
        else if (t == 1)
            pr->ql = strtol(q, &r, 10);
        else if (t == 2)
            pr->qs = strtol(q, &r, 10);
        else if (t == 3)
            pr->qe = strtol(q, &r, 10);
        else if (t == 4)
            pr->rev = (*q == '-');
        else if (t == 5)
            pr->tn = q;
        else if (t == 6)
            pr->tl = strtol(q, &r, 10);
        else if (t == 7)
            pr->ts = strtol(q, &r, 10);
        else if (t == 8)
            pr->te = strtol(q, &r, 10);
        else if (t == 9)
            pr->ml = strtol(q, &r, 10);
        else if (t == 10)
            pr->bl = strtol(q, &r, 10);
        ++t, q = i < l ? &s[i + 1] : 0;
    }
    if (t < 10)
        return -1;
    return 0;
}

int paf_read(paf_file_t *pf, paf_rec_t *r)
{
    int ret;
file_read_more:
    ret = ks_getuntil((kstream_t *)pf->fp, KS_SEP_LINE, &pf->buf, 0);
    if (ret < 0)
        return ret;
    ret = paf_parse(pf->buf.l, pf->buf.s, r);
    if (ret < 0)
        goto file_read_more;
    return ret;
}