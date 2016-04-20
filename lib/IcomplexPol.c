#include <stdio.h>
#include <stdlib.h>
#include "IcomplexPol.h"

#define MIN(a,b) (((a) < (b))? (a) : (b))
#define MAX(a,b) (((a) > (b))? (a) : (b))

static unsigned int *combinatorial;


/* COMPUTES ALL COMBINATORIAL NUMBERS */
/* (i over j) is stored in combinatorial[(i*(i+1))/2 + j] */
void startCombinatorialNumbers (int degree) {
    combinatorial = (unsigned int *) malloc ((((degree+2)*(degree+1))/2) * sizeof (unsigned int));
    unsigned int i;		// counter for m
    unsigned int j; 		// counter for n
    
    for (i=0; i<=degree; i++) {
        combinatorial[(i*(i+1))/2] = combinatorial[(i*(i+3))/2] = 1;        // FIRST AND LAST = 1
        
        // THE REST IS COMPUTED FROM THE RECURRENCE OF COMBINATORIAL NUMBERS.
        for (j=1; j<i; j++)
            combinatorial[(i*(i+1))/2+j] = combinatorial[((i-1)*i)/2+j-1] + combinatorial[((i-1)*i)/2+j];
        combinatorial[((i+1)*(i+2))/2-1] = 1;
    }    
}

void freeCombinatorialNumbers () {
    free (combinatorial);
}

Ipolynomial_t newZeroIPolynomial (int degree) {
    Ipolynomial_t rop;
    
    // rop points to address of a struct of _polynomial_t
    rop = (Ipolynomial_t) malloc (sizeof (struct _Ipolynomial_t));
    
    // set order of new polynomial
    rop->degree = degree;
    
    // allocate space for coefficients
    // number of coefficients = (order+1)*(order+2)/2
    int i;						// counter for coeficients
    int ncoef = ((degree+1)*(degree+2))/2;
    rop->coef = (Icomplex_t *) malloc (ncoef * sizeof (Icomplex_t));
    
    
    // set coeficients to zero.
    for (i=0; i<ncoef; i++) {
        (rop->coef)[i].real.a = (rop->coef)[i].real.b = 0.0;
        (rop->coef)[i].imag.a = (rop->coef)[i].imag.b = 0.0;
    }
    
    return rop;
}

Ipolynomial_t newIPolynomial (Icomplex_t *coef, int degree) {
    Ipolynomial_t rop;
    // rop points to address of a struct of _polynomial_t
    rop = (Ipolynomial_t) malloc (sizeof (struct _Ipolynomial_t));
    
    // set order of new polynomial
    rop->degree = degree;
    
    // allocate space for coefficients
    // number of coefficients = (order+1)*(order+2)/2
    int i;						// counter for coeficients
    int ncoef = ((degree+1)*(degree+2))/2;
    rop->coef = (Icomplex_t *) malloc (ncoef * sizeof (Icomplex_t));
    
    // copy coeficients to new polynomial.
    for (i=0; i<ncoef; i++) (rop->coef)[i] = coef[i];
    
    return rop;
}

void freeIPolynomial (Ipolynomial_t polynomial) {
    // frees memory for coeficients
    free (polynomial->coef);
    // frees memory for the polynomial structure
    free (polynomial);
}



void evalIPolHornerXY (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    //evaluate polynomial as an element of K[x][y] in Horner form
    //that is, each coefficient of y is a polynomial on x
    int i; //i will represent the degree on y
    int j; //j will represent the on x
    int d = polynomial->degree;
    Icomplex_t coef;
    *rop=(polynomial->coef)[(d+2)*(d+1)/2-1];
    for (i=d-1; i>=0; --i){
        ICFmul(rop, *rop, y);
        coef = (polynomial->coef)[d*(d+1)/2+i];
        for (j=d-1-i; j>=0; j--){
            ICFmul(&coef, coef, x);
            ICFadd(&coef, coef,  (polynomial->coef)[(i+j)*(i+j+1)/2+i]);
        }
        ICFadd(rop, *rop, coef);
    }
    return;
}



void evalIPolHornerYX (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    //evaluate polynomial as an element of K[y][x] in Horner form
    //that is, each coefficient of x is a polynomial on y
    int i; //i will represent the degree on x
    int j; //j will represent the on y
    int d = polynomial->degree;
    Icomplex_t coef;
    *rop=(polynomial->coef)[(d)*(d+1)/2];
    for (i=d-1; i>=0; --i){
        ICFmul(rop, *rop, x);
        coef = (polynomial->coef)[d*(d+1)/2+d-i];
        for (j=d-i-1; j>=0; --j){
            ICFmul(&coef, coef, y);
            ICFadd(&coef, coef, (polynomial->coef)[(i+j)*(i+j+1)/2+j]);
        }
        ICFadd(rop, *rop, coef);
    }
}

void evalIPolHorner (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    Icomplex_t ropx, ropy, ropc;
    evalIPolHornerXY (&ropx, polynomial, x, y);
    evalIPolHornerYX (&ropy, polynomial, x, y);
    evalIPol (&ropc, polynomial, x, y);
    rop->real.a = MAX(MAX(ropx.real.a, ropy.real.a), ropc.real.a);
    rop->real.b = MIN(MIN(ropx.real.b, ropy.real.b), ropc.real.b);
    rop->imag.a = MAX(MAX(ropx.imag.a, ropy.imag.a), ropc.imag.a);
    rop->imag.b = MIN(MIN(ropx.imag.b, ropy.imag.b), ropc.imag.b);
    return;
}


void evalIPol (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    // coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
    // with degree of second variable 'j' -> coef of x^(i-j) y^j
    int i;					// counter along degree of polynomial
    int j;					// counter for monomials of equal degree
    
    int k;					// counter for exponent of variables of monomial
    
    Icomplex_t aux;				// auxiliar variable to store monomials
    
    // compute all powers of x and y (from 0 to degree)
    Icomplex_t xPower[polynomial->degree+1], yPower[polynomial->degree+1];
    xPower[0].real.a = 1.0; xPower[0].real.b = 1.0; 
    xPower[0].imag.a = 0.0; xPower[0].imag.b = 0.0;
    yPower[0].real.a = 1.0; yPower[0].real.b = 1.0;
    yPower[0].imag.a = 0.0; yPower[0].imag.b = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree; j++) {
        ICFmul (&xPower[j+1], xPower[j], x);
        ICFmul (&yPower[j+1], yPower[j], y);
    }
    // sets output to 0+0i
    rop->real.a = rop->real.b = rop->imag.a = rop->imag.b = 0.0;
    
    // loop for each monomial of each degree
    for (i=0; i<=polynomial->degree; i++) for (j=0; j<=i; j++) {
        // starts monomial with coeficient
        aux = (polynomial->coef)[(i*(i+1))/2 + j];
        // multiplies by x^(i-j)
        ICFmul (&aux, aux, xPower[i-j]);
        // multiplies by y^(j)
        ICFmul (&aux, aux, yPower[j]);
        // accumulate to output
        ICFadd (rop, *rop, aux);
    }
}

void evalIPolX (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    // coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
    // with degree of second variable 'j' -> coef of x^(i-j) y^j
    int i;					// counter along degree of polynomial
    int j;					// counter for monomials of equal degree
    
    
    Icomplex_t aux;				// auxiliar variable to store monomials
    
    
    // compute all powers of x and y (from 0 to degree-1)
    Icomplex_t xPower[polynomial->degree], yPower[polynomial->degree];
    xPower[0].real.a = 1.0; xPower[0].real.b = 1.0; 
    xPower[0].imag.a = 0.0; xPower[0].imag.b = 0.0;
    yPower[0].real.a = 1.0; yPower[0].real.b = 1.0;
    yPower[0].imag.a = 0.0; yPower[0].imag.b = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-1; j++) {
        ICFmul (&xPower[j+1], xPower[j], x);
        ICFmul (&yPower[j+1], yPower[j], y);
    }
    
    // sets output to 0+0i
    rop->real.a = rop->real.b = rop->imag.a = rop->imag.b = 0.0;
    
    // loop for each monomial of each degree
    for (i=1; i<=polynomial->degree; i++) for (j=0; j<i; j++) {
        // starts monomial with (i-j)
        aux.real.a = aux.real.b = (double) (i-j); aux.imag.a = aux.imag.b = 0.0;
        // multiplies by coefficient
        ICFmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
        // multiplies by x^(i-j-1)
        ICFmul (&aux, aux, xPower[i-j-1]);
        // multiplies by y^(j)
        ICFmul (&aux, aux, yPower[j]);
        // accumulate to output
        ICFadd (rop, *rop, aux);
    }
}


void evalIPolYHornerXY (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    //evaluate the derivative of polynomial w.r.t. y as an element of K[x][y] in Horner form
    //that is, each coefficient of y is a polynomial on x
    int i; //i will represent the degree on y
    int j; //j will represent the on x
    int d = polynomial->degree; // degree of original polynomial
    Icomplex_t coef, aux;
    aux.real.a = aux.real.b = d;
    aux.imag.a = aux.imag.b = 0.0;
    *rop=(polynomial->coef)[(d+2)*(d+1)/2-1];
    ICFmul(rop, *rop, aux);
    for (i=d-1; i>0; --i){
        ICFmul(rop, *rop, y);
        coef = (polynomial->coef)[d*(d+1)/2+i];
        for (j=d-1-i; j>=0; j--){
            ICFmul(&coef, coef, x);
            ICFadd(&coef, coef,  (polynomial->coef)[(i+j)*(i+j+1)/2+i]);
        }
        aux.real.a = aux.real.b = i;
        ICFmul(&coef, coef, aux);
        ICFadd(rop, *rop, coef);
    }
    return;
}


void evalIPolYHornerYX (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    //evaluate the derivative of  polynomial wrt y as an element of K[y][x] in Horner form
    //that is, each coefficient of x is a polynomial on y
    int i; //i will represent the degree on x on the orginal polynomial
    int j; //j will represent the on y on the original polynomial
    int d = polynomial->degree; //degree of the original polynomial
    Icomplex_t coef, aux;
    *rop=(polynomial->coef)[(d)*(d+1)/2+1];
    for (i=d-2; i>=0; --i){
        ICFmul(rop, *rop, x);
        coef = (polynomial->coef)[d*(d+1)/2+d-i];
        aux.real.a = aux.real.b = d-i;
        aux.imag.a = aux.imag.b = 0.0;
        ICFmul(&coef, coef, aux);
        for (j=d-i-1; j>0; --j){
            ICFmul(&coef, coef, y);
            aux.real.a = aux.real.b = j;
            aux.imag.a = aux.imag.b = 0.0;
            ICFmul(&aux, (polynomial->coef)[(i+j)*(i+j+1)/2+j], aux);
            ICFadd(&coef, coef, aux);
        }
        ICFadd(rop, *rop, coef);
    }
}


void evalIPolY (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    // coef[(i*(i+1))/2 + j] is coeficient of monomial of degree 'i',
    // with degree of second variable 'j' -> coef of x^(i-j) y^j
    int i;					// counter along degree of polynomial
    int j;					// counter for monomials of equal degree
    
    
    Icomplex_t aux;				// auxiliar variable to store monomials
    
    
    // compute all powers of x and y (from 0 to degree-1)
    Icomplex_t xPower[polynomial->degree], yPower[polynomial->degree];
    xPower[0].real.a = 1.0; xPower[0].real.b = 1.0; 
    xPower[0].imag.a = 0.0; xPower[0].imag.b = 0.0;
    yPower[0].real.a = 1.0; yPower[0].real.b = 1.0;
    yPower[0].imag.a = 0.0; yPower[0].imag.b = 0.0;
    xPower[1] = x;
    yPower[1] = y;
    for (j=1; j<polynomial->degree-1; j++) {
        ICFmul (&xPower[j+1], xPower[j], x);
        ICFmul (&yPower[j+1], yPower[j], y);
    }
    
    // sets output to 0+0i
    rop->real.a = rop->real.b = rop->imag.a = rop->imag.b = 0.0;
    
    // loop for each monomial of each degree
    for (i=1; i<=polynomial->degree; i++) for (j=1; j<=i; j++) {
        // starts monomial with (j)
        aux.real.a = aux.real.b = (double) (j); aux.imag.a = aux.imag.b = 0.0;
        // multiplies by coefficient
        ICFmul (&aux, aux, (polynomial->coef)[(i*(i+1))/2 + j]);
        // multiplies by x^(i-j-1)
        ICFmul (&aux, aux, xPower[i-j]);
        // multiplies by y^(j)
        ICFmul (&aux, aux, yPower[j-1]);        
        // accumulate to output
        ICFadd (rop, *rop, aux);
    }
}

void evalIPolYHorner (Icomplex_t *rop, Ipolynomial_t polynomial, Icomplex_t x, Icomplex_t y) {
    Icomplex_t ropx, ropy, ropc;
    evalIPolYHornerXY (&ropx, polynomial, x, y);
    evalIPolYHornerYX (&ropy, polynomial, x, y);
    evalIPolY (&ropc, polynomial, x, y);
    rop->real.a = MAX(MAX(ropx.real.a, ropy.real.a), ropc.real.a);
    rop->real.b = MIN(MIN(ropx.real.b, ropy.real.b), ropc.real.b);
    rop->imag.a = MAX(MAX(ropx.imag.a, ropy.imag.a), ropc.imag.a);
    rop->imag.b = MIN(MIN(ropx.imag.b, ropy.imag.b), ropc.imag.b);
    return;
}


/* computes q(x,y) = p(x,y+ax) */
void IlinearVarChange (Ipolynomial_t q, Ipolynomial_t p, Icomplex_t a) {
    
    Icomplex_t combNumber;			// combinatorial number as complex_t
    Icomplex_t aux;				// auxiliar variable for internal sum
    
    
    combNumber.imag.a = combNumber.imag.b = 0.0;
    
    int i;					// counter for degree
    int j;					// counter for exponent of y
    int k;					// counter for elements modifying x^(i-j) y^j
    
    // for each degree i, only coefficients of monomials with j in {0, ..., i-1} are changed
    // for each j, only monomials with degree i and k in {j+1, ..., i} modify it
    // for each k, coef[(i*(i+1))/2+k] x^(i-k) (y+ax)^k
    // from binomial expansion, the term with x^(i-j) y^j has multiplier:
    // coef[(i*(i+1))/2+k] a^(k-j)  (k over j)

    
    Icomplex_t aPower[p->degree+1];
    aPower[0].real.a = aPower[0].real.b = 1.0;
    aPower[0].imag.a = aPower[0].imag.b = 0.0;
    aPower[1] = a;
    for (j=1; j<p->degree; j++) 
        ICFmul (&aPower[j+1], aPower[j], a);    
    
    for (i=0; i<=p->degree; i++) {
        for (j=0; j<=i; j++) {
            // basis coeficient
            q->coef[(i*(i+1))/2+j] = p->coef[(i*(i+1))/2+j];
            // minimum power of a is one
            for (k=j+1; k<=i; k++) {
                // set combinatorial number in complex_t variable
                combNumber.real.a = combNumber.real.b = combinatorial[(k*(k+1))/2+j];
                // aux = a^(k-j) (k over j)
                ICFmul (&aux, aPower[k-j], combNumber);
                // aux = coef[(i*(i+1))/2+k] * a^(k-j) (k over j)
                ICFmul (&aux, aux, p->coef[(i*(i+1))/2+k]);
                // updates new coeficient
                ICFadd (&(q->coef[(i*(i+1))/2+j]), q->coef[(i*(i+1))/2+j], aux);
            }
        }
    }
    
}

/* computes h(x,y) = q(x+x0,y) */
void IlinearVarChange2 (Ipolynomial_t h, Ipolynomial_t q, Icomplex_t x0) {
    
    Icomplex_t combNumber;			// combinatorial number as complex_t
    Icomplex_t aux;				// auxiliar variable for internal sum
    
    
    combNumber.imag.a = combNumber.imag.b = 0.0;
    
    int i;					// counter for degree
    int j;					// counter for exponent of y
    int k;					// counter for elements modifying x^(i-j) y^j
    
    // for each degree i, only coefficients of monomials with j in {0, ..., i-1} are changed
    // for each j, only monomials with degree i and k in {j+1, ..., i} modify it
    // for each k, coef[(i*(i+1))/2+k] x^(i-k) (y+ax)^k
    // from binomial expansion, the term with x^(i-j) y^j has multiplier:
    // coef[(i*(i+1))/2+k] a^(k-j)  (k over j)
    
    Icomplex_t x0Power[q->degree+1];
    x0Power[0].real.a = x0Power[0].real.b = 1.0;
    x0Power[0].imag.a = x0Power[0].imag.b = 0.0;
    x0Power[1] = x0;
    for (j=1; j<q->degree; j++) 
        ICFmul (&x0Power[j+1], x0Power[j], x0);
    
    
    
    for (i=0; i<=q->degree; i++) {
        for (j=0; j<=i; j++) {
            // basis coeficient
            h->coef[(i*(i+1))/2+j] = q->coef[(i*(i+1))/2+j];
            // minimum power of a is one
            for (k=i+1; k<=q->degree; k++) {
                // set combinatorial number in complex_t variable
                combNumber.real.a = combNumber.real.b = combinatorial[((k-j)*(k-j+1))/2 + k-i];
                // aux = a^(k-j) (k over j)
                ICFmul (&aux, x0Power[k-i], combNumber);
                // aux = coef[(i*(i+1))/2+k] * a^(k-j) (k over j)
                ICFmul (&aux, aux, q->coef[(k*(k+1))/2+j]);
                // updates new coeficient
                ICFadd (&(h->coef[(i*(i+1))/2+j]), h->coef[(i*(i+1))/2+j], aux);
            }
        }
    }
    
}

/*
int main () {
        double _coef[] = {0.000000000000000, -0.000000000000000, 0.000000000000000, 
-0.000000000000000, 0.000000000000000, -0.000000000000000, 
0.000000000000000, -0.000000000000000, -128.000000000000, 
-128.000000000000, -96.0000000000000, -96.0000000000000, 
0.000000000000000, -0.000000000000000, 0.000000000000000, 
-0.000000000000000, 72.0000000000000, 72.0000000000000, 
304.000000000000, 304.000000000000, 12.0000000000000, 
12.0000000000000, 180.000000000000, 180.000000000000, 
0.000000000000000, -0.000000000000000, 0.000000000000000, 
-0.000000000000000, 90.0000000000000, 90.0000000000000, 
-182.500000000000, -182.500000000000, 133.500000000000, 
133.500000000000, -190.500000000000, -190.500000000000, 
42.0000000000000, 42.0000000000000, -44.0000000000000, 
-44.0000000000000, 0.000000000000000, -0.000000000000000, 
0.000000000000000, -0.000000000000000, -40.2500000000000, 
-40.2500000000000, 18.2500000000000, 18.2500000000000, 
-72.0000000000000, -72.0000000000000, 21.0000000000000, 
21.0000000000000, -38.5000000000000, -38.5000000000000, 
5.50000000000000, 5.50000000000000, -6.00000000000000, 
-6.00000000000000, 0.000000000000000, -0.000000000000000};
        int i;
        Icomplex_t Icoef[15];
        for (i=0; i<15; i++) {
            Icoef[i].real.a = _coef[4*i];
            Icoef[i].real.b = _coef[4*i+1];
            Icoef[i].imag.a = _coef[4*i+2];
            Icoef[i].imag.b = _coef[4*i+3];
        }
        Ipolynomial_t p = newIPolynomial (Icoef, 4);
        
        Icomplex_t x = {{2.001,2.001102}, {3,3}};
        Icomplex_t y = {{4,4}, {-4.0001,-4}};
        
        Icomplex_t rop0;
        
        
        evalIPolY (&rop0, p, x, y);
        printf ("rop = [%e, %e] + i[%e,%e\n", rop0.real.a, rop0.real.b, rop0.imag.a, rop0.imag.b);
        printf ("radio = %e  %e\n", rop0.real.b - rop0.real.a, rop0.imag.b - rop0.imag.a);
        printf ("------------------------------------\n");
        
        
        evalIPolYHornerYX (&rop0, p, x, y);
        printf ("rop = [%e, %e] + i[%e,%e\n", rop0.real.a, rop0.real.b, rop0.imag.a, rop0.imag.b);
        printf ("radio = %e  %e\n", rop0.real.b - rop0.real.a, rop0.imag.b - rop0.imag.a);
        printf ("------------------------------------\n");
        
        
        evalIPolYHornerXY (&rop0, p, x, y);
        printf ("rop = [%e, %e] + i[%e,%e\n", rop0.real.a, rop0.real.b, rop0.imag.a, rop0.imag.b);
        printf ("radio = %e  %e\n", rop0.real.b - rop0.real.a, rop0.imag.b - rop0.imag.a);
        printf ("------------------------------------\n");
        
}*/
