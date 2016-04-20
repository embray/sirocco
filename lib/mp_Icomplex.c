#include <stdio.h>
#include <mpfr.h>
#include "mp_Icomplex.h"


void mp_ICFadd (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2) {
	mp_Iadd (&(rop->imag), op1.imag, op2.imag);
	mp_Iadd (&(rop->real), op1.real, op2.real);
}

void mp_ICFsub (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2) {
	mp_Isub (&(rop->imag), op1.imag, op2.imag);
	mp_Isub (&(rop->real), op1.real, op2.real);
}

void mp_ICFmul (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2) {
	mp_interval_t aux;
	mpfr_init (aux.a); mpfr_init (aux.b);
	mp_interval_t ropReal;
	mpfr_init (ropReal.a); mpfr_init (ropReal.b);
	
	
	mp_Imul (&aux, op1.imag, op2.imag);
	mp_Imul (&ropReal, op1.real, op2.real);
	mp_Isub (&ropReal, ropReal, aux);

	mp_Imul (&aux, op1.imag, op2.real);
	mp_Imul (&(rop->imag), op1.real, op2.imag);
	mp_Iadd (&(rop->imag), rop->imag, aux);	
	
	
	mpfr_set (rop->real.a, ropReal.a, MPFR_RNDN);
	mpfr_set (rop->real.b, ropReal.b, MPFR_RNDN);
	
	mpfr_clear (aux.a); mpfr_clear (aux.b);
	mpfr_clear (ropReal.a); mpfr_clear (ropReal.b);
}

void mp_ICFdiv (mp_Icomplex_t *rop, mp_Icomplex_t op1, mp_Icomplex_t op2) {
	/* TAKE INTO ACCOUNT THAT DENOMINATOR OF REAL AND IMAGINARY PART IS THE SAME*/
	
	
	mp_interval_t factor, aux;
	mpfr_init (aux.a); mpfr_init (aux.b);
	mpfr_init (factor.a); mpfr_init (factor.b);	
	mp_interval_t ropReal;
	mpfr_init (ropReal.a); mpfr_init (ropReal.b);
	
	mp_Imul (&aux, op2.imag, op2.imag);
	mp_Imul (&factor, op2.real, op2.real);
	mp_Iadd (&factor, factor, aux);
	mp_Iinv (&factor, factor);
	


	mp_Imul (&aux, op1.imag, op2.imag);
	mp_Imul (&ropReal, op1.real, op2.real);
	mp_Iadd (&ropReal, ropReal, aux);
	mp_Imul (&ropReal, ropReal, factor);
	



	mp_Imul (&aux, op1.real, op2.imag);
	mp_Imul (&(rop->imag), op1.imag, op2.real);
	mp_Isub (&(rop->imag), rop->imag, aux);
	mp_Imul (&(rop->imag), rop->imag, factor);
	
	
	mpfr_set (rop->real.a, ropReal.a, MPFR_RNDD);
	mpfr_set (rop->real.b, ropReal.b, MPFR_RNDU);
	
	mpfr_clear (aux.a); mpfr_clear (aux.b);
	mpfr_clear (factor.a); mpfr_clear (factor.b);
	mpfr_clear (ropReal.a); mpfr_clear (ropReal.b);
}


void mp_ICFprintf (mp_Icomplex_t x, char *s) {
	printf ("%s = [", s);
	mpfr_out_str (stdout, 10, 20, x.real.a, MPFR_RNDN);
	printf (", ");
	mpfr_out_str (stdout, 10, 20, x.real.b, MPFR_RNDN);
	printf ("] + \n");
	printf ("\t\t[");
	mpfr_out_str (stdout, 10, 20, x.imag.a, MPFR_RNDN);
	printf (", ");
	mpfr_out_str (stdout, 10, 20, x.imag.b, MPFR_RNDN);
	printf ("] i\n");

}
/*
int main () {
	mpfr_set_default_prec (60);
	mp_Icomplex_t x[3];
	int i;
	for (i=0; i<3; i++) {
		mpfr_init (x[i].real.a); mpfr_init (x[i].real.b);
		mpfr_init (x[i].imag.a); mpfr_init (x[i].imag.b);	
	}
	
	mpfr_set_str (x[0].real.a, "1", 10, MPFR_RNDN);
	mpfr_set_str (x[0].real.b, "2", 10, MPFR_RNDN);
	mpfr_set_str (x[0].imag.a, "3", 10, MPFR_RNDN);
	mpfr_set_str (x[0].imag.b, "4", 10, MPFR_RNDN);
	
	mpfr_set_str (x[1].real.a, "1", 10, MPFR_RNDN);
	mpfr_set_str (x[1].real.b, "2", 10, MPFR_RNDN);
	mpfr_set_str (x[1].imag.a, "3", 10, MPFR_RNDN);
	mpfr_set_str (x[1].imag.b, "4", 10, MPFR_RNDN);
	
	mp_ICFprintf (x[0], "\nx  "); printf ("-----------------\n");
	mp_ICFprintf (x[1], "y  "); printf ("-----------------\n");
	
	mp_ICFdiv (&x[2], x[0], x[1]);
	mp_ICFprintf (x[2], "x/x"); printf ("-----------------\n");
	

}
*/

