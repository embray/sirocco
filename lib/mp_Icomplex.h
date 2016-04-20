#ifndef _mp_Icomplex_H
#define _mp_Icomplex_H
#include "mp_interval.h"

typedef struct {
	mp_interval_t real;		// INTREVAL REAL PART
	mp_interval_t imag;		// INTERVAL IMAG PART
} mp_Icomplex_t;


void mp_ICFadd (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2);		// ROP = OP1 + OP2
void mp_ICFsub (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2);		// ROP = OP1 - OP2
void mp_ICFmul (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2);		// ROP = OP1 x OP2
void mp_ICFdiv (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2);		// ROP = OP1 / OP2


void mp_ICFprintf (mp_Icomplex_t x, char *s);


#endif

