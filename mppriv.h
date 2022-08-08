#ifndef MPPRIV_H
#define MPPRIV_H

#include "miniprot.h"

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifdef __cplusplus
extern "C" {
#endif

double mp_realtime(void);
double mp_cputime(void);
long mp_peakrss(void);
char *mp_strdup(const char *src);

#ifdef __cplusplus
}
#endif

#endif
