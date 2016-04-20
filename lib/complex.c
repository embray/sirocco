#include <stdio.h>
#include "complex.h"


// ADDITION OF TWO COMPLEX_T NUMBERS
void Cadd (complex_t *rop, complex_t op1, complex_t op2) {
	rop->real = op1.real + op2.real;
	rop->imag = op1.imag + op2.imag;
}


// SUBSTRACTION OF TWO COMPLEX_T NUMBERS
void Csub (complex_t *rop, complex_t op1, complex_t op2) {
	rop->real = op1.real - op2.real;
	rop->imag = op1.imag - op2.imag;
}

// PRODUCT OF TWO COMPLEX_T NUMBERS
void Cmul (complex_t *rop, complex_t op1, complex_t op2) {
	rop->real = op1.real * op2.real - op1.imag * op2.imag;
	rop->imag = op1.real * op2.imag + op1.imag * op2.real;
}


// DIVISION OF TWO COMPLEX_T NUMBERS
void Cdiv (complex_t *rop, complex_t op1, complex_t op2) {
	/* TAKE INTO ACCOUNT THAT DENOMINATOR OF REAL AND IMAGINARY PART IS THE SAME*/
	double factor = 1.0 / (op2.real*op2.real + op2.imag*op2.imag);

    rop->real =(op1.real * op2.real + op1.imag * op2.imag) * factor;
	rop->imag = (-op1.real * op2.imag + op1.imag * op2.real) * factor;
}


// SIMPLE PRINT FUNCTION
void Cprintf (complex_t op, char *s) {
	printf ("%s = %.16le + %.16le i\n", s, op.real, op.imag);
}


/*
int main () {
	complex_t x, y, z;
	
	x.real = 4.0;
	x.imag = -3.0;
	
	y.real = -3.0;
	y.imag = 5.0;
		
	printf ("-------------------------------------------------------------\n");
	Cprintf (x, "x  ");printf ("-------------------------------------------------------------\n");
	Cprintf (y, "y  ");printf ("-------------------------------------------------------------\n");
	
	Cadd (&z, x, y);
	Cprintf (z, "x+y");printf ("-------------------------------------------------------------\n");
	Csub (&z, x, y);
	Cprintf (z, "x-y");printf ("-------------------------------------------------------------\n");
	Cmul (&z, x, y);
	Cprintf (z, "x*y");printf ("-------------------------------------------------------------\n");
	Cdiv (&z, x, y);
	Cprintf (z, "x/y");printf ("-------------------------------------------------------------\n");
}*/
