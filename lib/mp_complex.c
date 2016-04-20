#include <stdio.h>
#include <mpfr.h>
#include "mp_complex.h"





void mp_Cadd (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2) {		// ROP = OP1 + OP2
	//(*rop).real = op1.real + op2.real;
	mpfr_add (rop->real, op1.real, op2.real, MPFR_RNDN); 
	//(*rop).imag = op1.imag + op2.imag;
	mpfr_add (rop->imag, op1.imag, op2.imag, MPFR_RNDN); 
}

void mp_Csub (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2) {		// ROP = OP1 - OP2
	//(*rop).real = op1.real - op2.real;
	mpfr_sub (rop->real, op1.real, op2.real, MPFR_RNDN); 
	//(*rop).imag = op1.imag - op2.imag;
	mpfr_sub (rop->imag, op1.imag, op2.imag, MPFR_RNDN); 
}

void mp_Cmul (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2) {		// ROP = OP1 x OP2
	mpfr_t aux; mpfr_init (aux);	
	mpfr_t ropReal; mpfr_init (ropReal);
	//(*rop).real = op1.real * op2.real - op1.imag * op2.imag;
	mpfr_mul (ropReal, op1.real, op2.real, MPFR_RNDN);
	mpfr_mul (aux, op1.imag, op2.imag, MPFR_RNDN);
	mpfr_sub (ropReal, ropReal, aux, MPFR_RNDN);	
	
	//(*rop).imag = op1.real * op2.imag + op1.imag * op2.real;	
	mpfr_mul (rop->imag, op1.real, op2.imag, MPFR_RNDN);
	mpfr_mul (aux, op1.imag, op2.real, MPFR_RNDN);
	mpfr_add (rop->imag, rop->imag, aux, MPFR_RNDN);

	mpfr_set (rop->real, ropReal, MPFR_RNDN);

	mpfr_clear (aux);
	mpfr_clear (ropReal);
}


void mp_Cdiv (mp_complex_t *rop, mp_complex_t op1, mp_complex_t op2) {		// ROP = OP1 / OP2
	mpfr_t aux, factor; 
	mpfr_init (aux); mpfr_init (factor);	
	mpfr_t ropReal; mpfr_init (ropReal);
	/* TAKE INTO ACCOUNT THAT DENOMINATOR OF REAL AND IMAGINARY PART IS THE SAME*/
	//double factor = 1.0 / (op2.real*op2.real + op2.imag*op2.imag);
	mpfr_mul (factor, op2.real, op2.real, MPFR_RNDN);
	mpfr_mul (aux, op2.imag, op2.imag, MPFR_RNDN);
	mpfr_add (factor, factor, aux, MPFR_RNDN);
	mpfr_ui_div (factor, 1UL, factor, MPFR_RNDN);

    //printf ("fac = %.15le\n", factor);
    //(*rop).real = (op1.real * op2.real + op1.imag * op2.imag) * factor;
    mpfr_mul (ropReal, op1.real, op2.real, MPFR_RNDN);
    mpfr_mul (aux, op1.imag, op2.imag, MPFR_RNDN);
    mpfr_add (ropReal, ropReal, aux, MPFR_RNDN);
    mpfr_mul (ropReal, ropReal, factor, MPFR_RNDN);
    
	//(*rop).imag = (-op1.real * op2.imag + op1.imag * op2.real) * factor;
	mpfr_mul (rop->imag, op1.imag, op2.real, MPFR_RNDN);
	mpfr_mul (aux, op1.real, op2.imag, MPFR_RNDN);
	mpfr_sub (rop->imag, rop->imag, aux, MPFR_RNDN);
	mpfr_mul (rop->imag, rop->imag, factor, MPFR_RNDN);
	
	mpfr_set (rop->real, ropReal, MPFR_RNDN);
	
	mpfr_clear (aux);
	mpfr_clear (factor);
	mpfr_clear (ropReal);
}

void mp_Cinit (mp_complex_t *op) {
	mpfr_init (op->real);
	mpfr_init (op->imag);
}


void mp_Cprintf (mp_complex_t op, char *s) {
	printf ("%s = ", s);
	mpfr_out_str (stdout, 10, (size_t) (mpfr_get_default_prec () / 3) + 1, op.real, MPFR_RNDN);
	printf (" + \n");
	mpfr_out_str (stdout, 10, (size_t) (mpfr_get_default_prec () / 3) + 1, op.imag, MPFR_RNDN);
	printf (" i \n");
}


/*

int main () {
	mpfr_set_default_prec (100);	// 100 bits -> 30 decimal digits
	mp_complex_t x, y, z;
	mp_Cinit (&x);
	mp_Cinit (&y);
	mp_Cinit (&z);
	
	mpfr_set_si (x.real, 4L, MPFR_RNDN);
	mpfr_set_si (x.imag, -3L, MPFR_RNDN);
	
	mpfr_set_si (y.real, -3L, MPFR_RNDN);
	mpfr_set_si (y.imag, 5L, MPFR_RNDN);
		
	printf ("-------------------------------------------------------------\n");
	mp_Cprintf (x, "x  ");printf ("-------------------------------------------------------------\n");
	mp_Cprintf (y, "y  ");printf ("-------------------------------------------------------------\n");
	
	mp_Cadd (&z, x, y);
	mp_Cprintf (z, "x+y");printf ("-------------------------------------------------------------\n");
	mp_Csub (&z, x, y);
	mp_Cprintf (z, "x-y");printf ("-------------------------------------------------------------\n");
	mp_Cmul (&z, x, y);
	mp_Cprintf (z, "x*y");printf ("-------------------------------------------------------------\n");
	mp_Cdiv (&z, x, y);
	mp_Cprintf (z, "x/y");printf ("-------------------------------------------------------------\n");
}
*/
