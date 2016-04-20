#include <stdio.h>
#include <stdlib.h>


#include "complexPol.h"

unsigned int *combinatorial;

// computes combinatorial number (m over n)
/*void startCombinatorialNumbers (int degree) {
	combinatorial = (unsigned int *) malloc (((degree*(degree+1))/2) * sizeof (unsigned int));
	unsigned int i;		// counter for m
	unsigned int j; 		// counter for n

	for (i=0; i<=degree; i++) {
		combinatorial[(i*(i+1))/2] = combinatorial[(i*(i+3))/2] = 1;
		for (j=1; j<i; j++)
			combinatorial[(i*(i+1))/2+j] = combinatorial[((i-1)*i)/2+j-1] + combinatorial[((i-1)*i)/2+j];
		combinatorial[((i+1)*(i+2))/2-1] = 1;
	}

}*/


void freePolynomial (polynomial_t polynomial) {
	// frees memory for coeficients
	free (polynomial->coef);
	// frees memory for the polynomial structure
	free (polynomial);
}


polynomial_t newZeroPolynomial (int degree) {
	polynomial_t rop;

	// rop points to address of a struct of _polynomial_t
	rop = (polynomial_t) malloc (sizeof (struct _polynomial_t));

	// set order of new polynomial
	rop->degree = degree;

	// allocate space for coefficients
	// number of coefficients = (order+1)*(order+2)/2
	int i;						// counter for coeficients
	int ncoef = ((degree+1)*(degree+2))/2;
	rop->coef = (complex_t *) malloc (ncoef * sizeof (complex_t));

	// set coeficients to zero.
	for (i=0; i<ncoef; i++) (rop->coef)[i].real = (rop->coef)[i].imag = 0.0;

	return rop;
}

polynomial_t newPolynomial (complex_t *coef, int degree) {
	polynomial_t rop;
	// rop points to address of a struct of _polynomial_t
	rop = (polynomial_t) malloc (sizeof (struct _polynomial_t));

	// set order of new polynomial
	rop->degree = degree;

	// allocate space for coefficients
	// number of coefficients = (order+1)*(order+2)/2
	int i;						// counter for coeficients
	int ncoef = ((degree+1)*(degree+2))/2;
	rop->coef = (complex_t *) malloc (ncoef * sizeof (complex_t));

	// copy coeficients to new polynomial.
	for (i=0; i<ncoef; i++) (rop->coef)[i] = coef[i];

	return rop;
}


void evalPol (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree

	int k;					// counter for exponent of variables of monomial

	complex_t aux;				// auxiliar variable to store monomials

    // compute all powers of x and y (from 0 to degree)
    complex_t xPower[polynomial->degree+1], yPower[polynomial->degree+1];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=0; i<=polynomial->degree; i++) for (j=0; j<=i; j++) {
		// starts monomial with coeficient
		aux = (polynomial->coef)[(i*(i+1))/2 + j];
		// multiplies by x^(i-j)
		Cmul (&aux, aux, xPower[i-j]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j]);
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}

void evalPolX (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree


	complex_t aux;				// auxiliar variable to store monomials

    
    // compute all powers of x and y (from 0 to degree-1)
    complex_t xPower[polynomial->degree], yPower[polynomial->degree];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-1; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
         
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=1; i<=polynomial->degree; i++) for (j=0; j<i; j++) {
		// starts monomial with (i-j)
        aux.real = (double) (i-j); aux.imag = 0.0;
        // multiplies by coefficient
		Cmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
		// multiplies by x^(i-j-1)
		Cmul (&aux, aux, xPower[i-j-1]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j]);
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}


void evalPolY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree


	complex_t aux;				// auxiliar variable to store monomials
   
    
    // compute all powers of x and y (from 0 to degree-1)
    complex_t xPower[polynomial->degree], yPower[polynomial->degree];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-1; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
         
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=1; i<=polynomial->degree; i++) for (j=1; j<=i; j++) {
		// starts monomial with (j)
        aux.real = (double) (j); aux.imag = 0.0;
        // multiplies by coefficient
		Cmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
		// multiplies by x^(i-j-1)
		Cmul (&aux, aux, xPower[i-j]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j-1]);        
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}


void evalPolYY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree


	complex_t aux;				// auxiliar variable to store monomials
   
    
    // compute all powers of x and y (from 0 to degree-1)
    complex_t xPower[polynomial->degree-1], yPower[polynomial->degree-1];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-2; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
         
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=2; i<=polynomial->degree; i++) for (j=2; j<=i; j++) {
		// starts monomial with (j)
        	aux.real = j*(j-1.0); aux.imag = 0.0;
        // multiplies by coefficient
		Cmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
		// multiplies by x^(i-j-1)
		Cmul (&aux, aux, xPower[i-j]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j-2]);        
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}


void evalPolXY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree


	complex_t aux;				// auxiliar variable to store monomials
   
    
    // compute all powers of x and y (from 0 to degree-1)
    complex_t xPower[polynomial->degree-1], yPower[polynomial->degree-1];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-2; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
         
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=2; i<=polynomial->degree; i++) for (j=1; j<i; j++) {
		// starts monomial with (j)
        	aux.real = (double) (j*(i-j)); aux.imag = 0.0;
        // multiplies by coefficient
		Cmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
		// multiplies by x^(i-j-1)
		Cmul (&aux, aux, xPower[i-j-1]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j-1]);        
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}


void evalPolXX (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y) {
	// coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
	// with degree of second variable 'j' -> coef of x^(i-j) y^j
	int i;					// counter along degree of polynomial
	int j;					// counter for monomials of equal degree


	complex_t aux;				// auxiliar variable to store monomials
   
    
    // compute all powers of x and y (from 0 to degree-1)
    complex_t xPower[polynomial->degree-1], yPower[polynomial->degree-1];
    xPower[0].real = 1.0; 
    xPower[0].imag = 0.0;
    yPower[0].real = 1.0;
    yPower[0].imag = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-2; j++) {
        Cmul (&xPower[j+1], xPower[j], x);
        Cmul (&yPower[j+1], yPower[j], y);
    }
         
	// sets output to 0+0i
	rop->real = rop->imag = 0.0;

	// loop for each monomial of each degree
	for (i=2; i<=polynomial->degree; i++) for (j=0; j<i-1; j++) {
		// starts monomial with (j)
        	aux.real = (double) (i-j)*(i-j-1.0); aux.imag = 0.0;
        // multiplies by coefficient
		Cmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
		// multiplies by x^(i-j-1)
		Cmul (&aux, aux, xPower[i-j-2]);
		// multiplies by y^(j)
		Cmul (&aux, aux, yPower[j]);        
		// accumulate to output
		Cadd (rop, *rop, aux);
	}
}


