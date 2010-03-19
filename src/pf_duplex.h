#ifndef __INC_PF_DUPLEX_H__
#define __INC_PF_DUPLEX_H__

extern double **pr_duplex;
extern double **pr_duplex2;
double pf_duplex(const char *s1, const char *s2);
void free_pf_duplex();

#endif  /*  __INC_PF_DUPLEX_H__ */
