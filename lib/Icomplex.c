#include <stdio.h>
#include <fenv.h>
#include "Icomplex.h"
#include "interval.h"


#define MIN(a,b) (((a) < (b))? (a) : (b))
#define MAX(a,b) (((a) > (b))? (a) : (b))

void ICFadd (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2) {
	Iadd (&(rop->real), op1.real, op2.real);
	Iadd (&(rop->imag), op1.imag, op2.imag);
}

void ICFsub (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2) {
	Isub (&(rop->imag), op1.imag, op2.imag);
	Isub (&(rop->real), op1.real, op2.real);
}

void ICFmul (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2) {
	interval_t aux, ropReal;

	Imul (&aux, op1.imag, op2.imag);
	Imul (&ropReal, op1.real, op2.real);
	Isub (&ropReal, ropReal, aux);
	
	

	Imul (&aux, op1.imag, op2.real);
	Imul (&(rop->imag), op1.real, op2.imag);
	Iadd (&(rop->imag), rop->imag, aux);	
	rop->real = ropReal;
}


// inverse of an interval [a+b] * c*i
void ICFinv_hor(Icomplex_t *rop, double a, double b, double c) {
	int oldRoundingMode = fegetround ();
	// REAL INTERVAL
	if (c == 0.0) {
		rop->imag.a = rop->imag.b = 0.0;
		fesetround (FE_DOWNWARD);
		rop->real.a = 1.0 / b;
		fesetround (FE_UPWARD);
		rop->real.b = 1.0 / a;
		
		fesetround (oldRoundingMode);
		return;
	}
	// INTERVAL WITH POSITIVE IMAGINARY PART
	if (c > 0.0) {
		double  r= 0.5 / c;
		// case a<0<b
		if (a < 0.0 && b > 0.0) {
			fesetround (FE_DOWNWARD);
			rop->imag.a = -1.0 /c;
			fesetround (FE_UPWARD);
			rop->imag.b = MAX (-c / (a*a + c*c), -c / (b*b + c*c));
			
			if (0.5/c - c / (a*a+c*c) >= 0.0) {
				fesetround (FE_DOWNWARD);
				rop->real.a = -0.5/c;
			}
			else {
				fesetround (FE_DOWNWARD);
				rop->real.a = a / (a*a+c*c);
			}
			
			if (c / (b*b+c*c) - 0.5/c <= 0.0) {
				fesetround (FE_UPWARD);
				rop->real.b = 0.5/c;
			}
			else {
				fesetround (FE_UPWARD);
				rop->real.b = b / (b*b + c*c);
			}
			// RETURN
			fesetround (oldRoundingMode);
			return;
			
		}
		
		// case 0<=a<b
		if (0<=a) {
			fesetround (FE_UPWARD);
			rop->imag.b = -c / (b*b+c*c);
			fesetround (FE_DOWNWARD);
			rop->imag.a = -c / (a*a+c*c);
		
			rop->real.a = MIN (a/(a*a + c*c), b/(b*b + c*c));
			
			
			if (a<=c && c <= b) {
				fesetround (FE_UPWARD);
				rop->real.b = 0.5/c;
			} else {
				fesetround (FE_UPWARD);
				rop->real.b = MAX (a/(a*a + c*c), b/(b*b + c*c));
			}
			
		
			// RETURN
			fesetround (oldRoundingMode);
			return;
		
		}
		
		// case a<b <= 0
		if (0>=b) {
			fesetround (FE_DOWNWARD);
			rop->imag.a = -c / (b*b+c*c);
			fesetround (FE_UPWARD);
			rop->imag.b = -c / (a*a+c*c);
		
			rop->real.b = MAX (a/(a*a + c*c), b/(b*b + c*c));
			
			
			if (a<=-c && -c <= b) {
				fesetround (FE_DOWNWARD);
				rop->real.a = -0.5/c;
			} else {
				fesetround (FE_DOWNWARD);
				rop->real.a = MIN (a/(a*a + c*c), b/(b*b + c*c));
			}
			
		
			// RETURN
			fesetround (oldRoundingMode);
			return;
		
		}		
	}
	
	// UP TO HERE, C<0
	ICFinv_hor (rop, a, b, -c);
	double aux = rop->imag.b;
	rop->imag.b = -rop->imag.a;
	rop->imag.a = -aux;
	fesetround (oldRoundingMode);
	return;
	
	}


// inverse of an interval x = c + i*[a+b]
void ICFinv_ver(Icomplex_t *rop, double a, double b, double c) {
	// (x)^-1 = i (i x)^-1
	
	// ix has its  (a,b,c) as  (-b,-a, c) from x
	ICFinv_hor (rop, -b, -a, c);		// rop stores inverse of (ix)
	
	interval_t aux = rop->imag;
	rop->imag = rop->real;
	rop->real.a = -aux.b;
	rop->real.b = -aux.a;
	
	return;	

}


void ICFinv (Icomplex_t *rop, Icomplex_t op) {
	// compute inverse of four edges
	Icomplex_t inv_edge[4];
	
	ICFinv_hor (&inv_edge[0], op.real.a, op.real.b, op.imag.a);
	ICFinv_hor (&inv_edge[1], op.real.a, op.real.b, op.imag.b);
	ICFinv_ver (&inv_edge[2], op.imag.a, op.imag.b, op.real.a);
	ICFinv_ver (&inv_edge[3], op.imag.a, op.imag.b, op.real.b);
	
	rop->real.a = MIN (MIN (MIN(inv_edge[0].real.a, inv_edge[1].real.a), 
				inv_edge[2].real.a), inv_edge[3].real.a);
	rop->real.b = MAX (MAX (MAX(inv_edge[0].real.b, inv_edge[1].real.b), 
				inv_edge[2].real.b), inv_edge[3].real.b);
	
	rop->imag.a = MIN (MIN (MIN(inv_edge[0].imag.a, inv_edge[1].imag.a), 
				inv_edge[2].imag.a), inv_edge[3].imag.a);
	rop->imag.b = MAX (MAX (MAX(inv_edge[0].imag.b, inv_edge[1].imag.b), 
				inv_edge[2].imag.b), inv_edge[3].imag.b);


	return;
	
}

void ICFdiv (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2) {
	Icomplex_t inverse;
	ICFinv (&inverse, op2);
	ICFmul (rop, op1, inverse);

}

void ICFdiv2 (Icomplex_t *rop, Icomplex_t op1, Icomplex_t op2) {
	/* TAKE INTO ACCOUNT THAT DENOMINATOR OF REAL AND IMAGINARY PART IS THE SAME*/
	interval_t factor, aux, ropReal;
	Imul (&aux, op2.imag, op2.imag);
	Imul (&factor, op2.real, op2.real);
	Iadd (&factor, factor, aux);
	
	
	Iinv (&factor, factor);

	Imul (&aux, op1.imag, op2.imag);
	Imul (&ropReal, op1.real, op2.real);
	Iadd (&ropReal, ropReal, aux);
	Imul (&ropReal, ropReal, factor);

	Imul (&aux, op1.real, op2.imag);
	Imul (&(rop->imag), op1.imag, op2.real);
	Isub (&(rop->imag), rop->imag, aux);
	Imul (&(rop->imag), rop->imag, factor);	
	
	rop->real = ropReal;
}


void ICFprintf (Icomplex_t x, char *s) {
	printf ("%s = [%.16le, %.16le] + \n[%.16le, %.16le] i\n", 
			s, x.real.a, x.real.b, x.imag.a, x.imag.b);

}
/*
int main () {
	Icomplex_t x[3];
	int i;

		
	x[0].real.a = 1.0;
	x[0].real.b = 2.0;
	x[0].imag.a = 3.0;
	x[0].imag.b = 4.0;
	
	x[1].real.a = 1.0;
	x[1].real.b = 2.0;
	x[1].imag.a = 3.0;
	x[1].imag.b = 4.0;
	
	ICFprintf (x[0], "\nx  "); printf ("-----------------\n");
	ICFprintf (x[1], "y  "); printf ("-----------------\n");
	
	ICFdiv (&x[2], x[0], x[1]);
	ICFprintf (x[2], "x/x"); printf ("-----------------\n");


	x[1].real.a = -3.1;
	x[1].real.b = -2.5;
	x[1].imag.a = 0.15;
	x[1].imag.b = 5.3;

	x[2].real.a = -0.96;
	x[2].real.b = 2.15;
	x[2].imag.a = 9.1;
	x[2].imag.b = 9.7;

	

	ICFinv (&x[0], x[2]);
	ICFprintf (x[0], "inverse"); printf ("-----------------\n");

	ICFdiv (&x[0], x[1], x[2]);
	ICFprintf (x[0], "new"); printf ("-----------------\n");
	printf ("diametro = %e %e\n\n", x[0].real.b - x[0].real.a, x[0].imag.b - x[0].imag.a);
	
	ICFdiv2 (&x[0], x[1], x[2]);
	ICFprintf (x[0], "old"); printf ("-----------------\n");
	printf ("diametro = %e %e\n\n", x[0].real.b - x[0].real.a, x[0].imag.b - x[0].imag.a);

		

}

*/

