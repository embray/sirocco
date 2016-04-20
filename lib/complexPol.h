#ifndef _complexPol_H
#define _complexPol_H

#include "complex.h"


//COEFICIENTS OF THE 2VAR-POLYNOMIAL pol (x,y) ARE STORED AS:
// [1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3, ...]


typedef struct _polynomial_t {
	complex_t *coef;			// COEFICIENTS OF THE POLYNOMIAL
	int degree;				// MAX DEGREE OF THE POLINOMIAL
}* polynomial_t;


/* ALLOCATES MEMORY FOR COMBINATORIAL NUMBERS AND INITIALIZES THEM */
//void startCombinatorialNumbers (int degree);

/* GIVEN A DEGREE, ALLOCATES SPACE MEMORY FOR AN EMPTY POLYNOMIAL OF THAT DEGREE */
polynomial_t newZeroPolynomial (int degree);

/* GIVEN AN ARRAY OF COMPLEX_T COEFFICIENTS, AND THE DEGREE OF THE POLINOMIAL, ALLOCATES
   SPACE IN MEMORY FOR THE POLYNOMIAL AND RETURNS A POINTER TO THIS ADDRESS */
polynomial_t newPolynomial (complex_t *coef, int degree);

/* GIVEN A POLYNOMIAL, FREES THE ALLOCATED MEMORY FOR ITS COEFICIENTS AND THE POLYNOMIAL */
void freePolynomial (polynomial_t polynomial);


/* EVALUATION OF A POLYNOMIAL AND ITS DERIVATIVES */

/* ROP = POLYNOMIAL (Z1, Z2) */
void evalPol (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);
/* ROP = pX (Z1, Z2) -> EVALUATES THE DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO x EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalPolX (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);
/* ROP = pY (Z1, Z2) -> EVALUATES THE DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO y EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalPolY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);
/* ROP = pXX(Z1, Z2) -> EVALUATES THE SECOND DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO x TWO TIMES EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalPolXX (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);
/* ROP = pYY(Z1, Z2) -> EVALUATES THE SECOND DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO x and y EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalPolXY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);
/* ROP = pYY(Z1, Z2) -> EVALUATES THE SECOND DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO y TWO TIMES EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalPolYY (complex_t *rop, polynomial_t polynomial, complex_t x, complex_t y);



#endif


