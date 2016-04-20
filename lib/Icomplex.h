#ifndef _Icomplex_H
#define _Icomplex_H
#include "interval.h"

typedef struct {
	interval_t real;		// INTREVAL REAL PART
	interval_t imag;		// INTERVAL IMAG PART
} Icomplex_t;


void ICFadd (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2);		// ROP = OP1 + OP2
void ICFsub (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2);		// ROP = OP1 - OP2
void ICFmul (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2);		// ROP = OP1 x OP2
void ICFinv (Icomplex_t *rop, Icomplex_t op);						// ROP = 1 / OP
void ICFdiv (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2);		// ROP = OP1 / OP2
void ICFdiv2 (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2);		// ROP = OP1 / OP2  // OLD

#endif

