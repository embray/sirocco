#ifndef _mp_complex_H
#define _mp_complex_H

typedef struct {
	mpfr_t real;				// REAL PART
	mpfr_t imag;				// IMAG PART
} mp_complex_t;


void mp_Cadd (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2);		// ROP = OP1 + OP2
void mp_Csub (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2);		// ROP = OP1 - OP2
void mp_Cmul (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2);		// ROP = OP1 x OP2
void mp_Cdiv (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2);		// ROP = OP1 / OP2

#endif

