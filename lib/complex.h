#ifndef _complex_H
#define _complex_H

typedef struct {
	double real;				// REAL PART
	double imag;				// IMAG PART
} complex_t;


// ADDITION OF TWO COMPLEX_T NUMBERS
void Cadd (complex_t *rop, complex_t op1, complex_t op2);		// ROP = OP1 + OP2

// SUBSTRACTION OF TWO COMPLEX_T NUMBERS
void Csub (complex_t *rop, complex_t op1, complex_t op2);		// ROP = OP1 - OP2

// PRODUCT OF TWO COMPLEX_T NUMBERS
void Cmul (complex_t *rop, complex_t op1, complex_t op2);		// ROP = OP1 x OP2

// DIVISION OF TWO COMPLEX_T NUMBERS
void Cdiv (complex_t *rop, complex_t op1, complex_t op2);		// ROP = OP1 / OP2

// SIMPLE PRINT FUNCTION
void Cprintf (complex_t op, char *s);

#endif



