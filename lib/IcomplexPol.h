#ifndef _IcomplexPol_H
#define _IcomplexPol_H

#include "complexPol.h"
#include "Icomplex.h"


typedef struct _Ipolynomial_t {
	Icomplex_t *coef;			// COEFICIENTS OF THE POLYNOMIAL
	int degree;				// MAX DEGREE OF THE POLINOMIAL
}* Ipolynomial_t;




void startCombinatorialNumbers (int degree);
void freeCombinatorialNumbers ();


/* GIVEN A DEGREE, ALLOCATES SPACE MEMORY FOR AN EMPTY POLYNOMIAL OF THAT DEGREE */
Ipolynomial_t newZeroIPolynomial (int degree);

/* GIVEN AN ARRAY OF COMPLEX_T COEFFICIENTS, AND THE DEGREE OF THE POLINOMIAL, ALLOCATES
   SPACE IN MEMORY FOR THE POLYNOMIAL AND RETURNS A POINTER TO THIS ADDRESS */
Ipolynomial_t newIPolynomial (Icomplex_t *coef, int degree);

/* GIVEN A COMPLEX POLYNOMIAL WITH EXACT VALUES, CREATES THE INTERVAL COMPLEX POLYNOMIAL */
Ipolynomial_t convertIpolynomial (polynomial_t p);

/* GIVEN A POLYNOMIAL, FREES THE ALLOCATED MEMORY FOR ITS COEFICIENTS AND THE POLYNOMIAL */
void freeIPolynomial (Ipolynomial_t polynomial);

/* ROP = POLYNOMIAL (Z1, Z2) */
void evalIPolHorner (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);

/* ROP = POLYNOMIAL (Z1, Z2) */
void evalIPol (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);

/* ROP = pX (Z1, Z2) -> EVALUATES THE DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO x EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalIPolX (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);

/* ROP = pY (Z1, Z2) -> EVALUATES THE DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO y EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalIPolY (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);

/* ROP = pY (Z1, Z2) -> EVALUATES THE DERIVATIVE OF THE POLYNOMIAL WITH RESPECT TO y EVALUATED IN (x,y). NO NEED OF COMPUTING */
void evalIPolYHorner (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);
void evalIPolYHornerXY (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);
void evalIPolYHornerYX (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y);

/* ROP = q (Z1,Z1) -> EVALUATES THE DE LINEAR TRANSFORMATION OF THE ORIGINAL POLYNOMIAL q(x,y)=p(x+x0,y+a(x+x0)). */
void evalTranslatedIPol (Icomplex_t *rop, Ipolynomial_t p, Icomplex_t x, Icomplex_t y, Icomplex_t a);


/* computes q(x,y) = p(x,y+ax) */
void IlinearVarChange (Ipolynomial_t q, Ipolynomial_t p, Icomplex_t a);

/* computes h(x,y) = q(x+x0,y)*/
void IlinearVarChange2 (Ipolynomial_t h, Ipolynomial_t q, Icomplex_t x0);


#endif

