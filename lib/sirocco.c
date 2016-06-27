#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include "complex.h"
#include "complexPol.h"
#include "Icomplex.h"
#include "IcomplexPol.h"
#include "list.h"


/******************************************
// f(z1,z2) = z1^3+2*z1*z2-3*z2-z1;
// COEFS = {{1,0,0,0}, {0,2,0}, {-1,3}, {0}}
// z1 = x + iy
// z2 = t(1+0.01i)
******************************************/

#define MAX(a,b)	(((a) < (b))? (b) : (a))
#define MIN(a,b)	(((a) > (b))? (b) : (a))

// PERFORM A MINIMUM OF 5 ITERATIONS OF NEWTON METHOD TO CORRECT THE APPROXIMATION OF THE ROOT (SUPPOSED TO BE GOOD ENOUGH)
void correctRoot (polynomial_t f, complex_t *x0, complex_t *y0) {
	complex_t evalF, aux, aux1;
	int count = 0;
	double last, new;			// sizes of corrections
	new = 9999;
	do {
		last = new;
		evalPol (&evalF, f, *x0, *y0);
		evalPolY (&aux1, f, *x0, *y0);
		Cdiv (&aux, evalF, aux1);
			new = fabs(aux.real) + fabs(aux.imag);
		Csub (y0, *y0, aux);
	  	++count;
	} while ((fabs(last - new)  > 1.0e-14 || count < 6 ) && count < 50);
	return;
}


// a = -fx / fy
// computes the derivative at x0, of the implicit curve (f(x,y(x)) = 0)
void getA (Icomplex_t *Ia, polynomial_t f, complex_t x0, complex_t y0) {
	complex_t a, aux;
	evalPolX (&a, f, x0, y0);
	evalPolY (&aux, f, x0, y0);
	Cdiv (&a, a, aux);
	a.real = -a.real; a.imag = -a.imag;
	(*Ia).real.a = (*Ia).real.b = a.real;
	(*Ia).imag.a = (*Ia).imag.b = a.imag;
}


// PERFORMS NEWTON TEST
int newtonTest (Ipolynomial_t If, Icomplex_t Ix0, Icomplex_t Iy0, Icomplex_t IY0) {
	Icomplex_t N, Iaux;
	
	
	evalIPolHorner (&N, If, Ix0, Iy0);
	evalIPolYHorner (&Iaux, If, Ix0, IY0); 
	
	// check we do not perform division by an interval containing 0. 
	// in this case, return fail the test
	if (isZeroContained(Iaux.real) && isZeroContained(Iaux.imag)) return 0;

	ICFdiv (&N, N, Iaux);
	ICFsub (&N, Iy0, N);

	if (isSubsetEqual (N.real, IY0.real) && isSubsetEqual(N.imag, IY0.imag))
		return 1;
	else 
		return 0;
}

// RETURNS WHETHER NEWTON EVALUATION N IS CONTAINED IN THE INTERVAL ENCLOSURE IY0 OF Iy0. 
// BESIDES REPEATS WITH A ENCLOSURE THREE TIMES BIGER. RETURNS SUCCESS IF BOTH ARE SATISFIED
int validatePoint (Ipolynomial_t If, Icomplex_t Ix0, Icomplex_t Iy0, Icomplex_t IY0) {
	Icomplex_t IY1;

	if (!newtonTest (If, Ix0, Iy0, IY0)) return 0;

	// COMPUTES AN ENCLOSURE THREE TIMES BIGGER FOR N2
	int oldRoundingMode = fegetround ();

	// CREATE A 3-TIMES BIGGER ENCLOSURE (ROUNDING TO BIGGER)
	fesetround (FE_DOWNWARD);
	IY1.real.a = 2*IY0.real.a - IY0.real.b;
	IY1.imag.a = 2*IY0.imag.a - IY0.imag.b;
	fesetround (FE_UPWARD);
	IY1.imag.b = 2*IY0.imag.b - IY0.imag.a;
	IY1.real.b = 2*IY0.real.b - IY0.real.a;
	fesetround (oldRoundingMode);
	
	if (!newtonTest (If, Ix0, Iy0, IY1)) return 0;

	return 1;
	
}

#define validateStep	validatePoint



// CORE FUNCTION OF THE LIBRARY
double* homotopyPath (int degree, 				// degree of the input polynomial
		double *_coef,							// list of coefficients interval and complex
		double _y0R,							// real part of approximation of root. f(0,y0) = 0
		double _y0I								// imag part of approximation of the root
		) {
	startCombinatorialNumbers(degree);			// initialize combinatorial numbers (compute once)



	list_t outputList = {0, NULL};				// initialize dynamic list for output
	double dataToList[3] = {0.0, _y0R, _y0I};
	appendData (&outputList, dataToList);
	
	
	int i, j, k;								// counters
	double stepsize = 0.0001;					// step size
	double eps;									// radius of the box
	
	polynomial_t f;								// complex polynomials (with fp arithmetic)
												// to be used in Newton method correction
	
	Ipolynomial_t If, Ig, Ih;					// complex polynomials (with Interval arithmetic)
												// to be used in validated homotopy
	
	
	complex_t a;								// y' (x0) tangent vector to the zero-curve
	Icomplex_t Ia;								// y'(x0) as Icomplex_t
	
	
	complex_t x0;								// x0 of the root (x0.imag = 0)
	complex_t y0;								// approximated y0 of the root
	Icomplex_t Ix0;								// Interval representation for x0
	Icomplex_t Iy0;								// Interval representation for y0
	Icomplex_t IY0;								// Interval enclosure for y0 for Newton Method
	Icomplex_t Iy1, IY1;						// Interval representation of zeros for rotated polynomial
	Icomplex_t Ix1;
	
	
	complex_t aux, aux1, aux2, auxy0;			// auxiliar variable for complex_t computation
	Icomplex_t Iaux;							// auxiliar varialbe for Icomplex_t computations
	
	Icomplex_t IlastY0;							// to prove that the boxes are connected
		IlastY0.real.a = -INFINITY;				// REMOVE???
		IlastY0.real.b = INFINITY;
		IlastY0.imag.a = -INFINITY;
		IlastY0.imag.b = INFINITY;
	
	
	int nCoef = ((degree+1)*(degree+2)) / 2;	// number of monomials
	complex_t coef[nCoef];						// coefficient list for floating point representation of
												// the polynomial.
	Icomplex_t Icoef[nCoef]; 					// coefficient list for Interval representation of
												// the polynomial.



	/***********************/
	/* DATA INITIALIZATION */
	/***********************/
	
	// INITIALIZE COEFICIENT DATA FROM ARGUMENT INTO APPROPIATE DATA STRUCTURES
	for (i=0; i<nCoef; i++)
		for (j=0; j<4; j++) {
			Icoef[i].real.a = _coef[4*i];
			Icoef[i].real.b = _coef[4*i+1];
			Icoef[i].imag.a = _coef[4*i+2];
			Icoef[i].imag.b = _coef[4*i+3];
			coef[i].real = (_coef[4*i] + _coef[4*i+1]) / 2.0;
			coef[i].imag = (_coef[4*i+2] + _coef[4*i+3]) / 2.0;
		}
	

	// INITIALIZE FLOATING POINT COMPLEX APPROXIMATION OF THE ROOT
	x0.real = 0.0; x0.imag = 0.0;
	y0.real = _y0R; y0.imag = _y0I;	
		
	
	// INITIALIZE AND ALLOCATE SPACE FOR POLYNOMIALS
	f = newPolynomial (coef, degree);
	If = newIPolynomial (Icoef, degree);
	Ig = newZeroIPolynomial (degree);
	Ih = newZeroIPolynomial (degree);	


#ifdef PROFILING	
	int nsteps = 0;								// number of steps performed
	int nrejectedEps = 0;						// number of rejected validations of boxes
	int nrejectedSteps = 0;						// number of rejected validations of tubes
#endif



	/*******************/
	/* CORRECT ROOT	   */
	/*******************/
	//correctRoot (f, &x0, &y0);

	
	/*******************/
	/** MAIN LOOP	   */
	/*******************/
	while (x0.real < 1.0) {

#ifdef PROFILING
	nsteps++; // used for profiling only
#endif

		/***********************************************************/
		/*** BOX SIZE ESTIMATION. F RESTRICTED TO X0 IS INJECTIVE. */
		/***********************************************************/
		evalPolY (&aux1, f, x0, y0);
		evalPolYY (&aux2, f, x0, y0);
		if (fabs (aux2.real) + fabs (aux.imag) < 1e-10) {
			eps = 0.5;
		} else {
			Cdiv (&aux, aux1, aux2);
			eps = sqrt(aux.real*aux.real + aux.imag * aux.imag)/8.0;
		}

		/**************************/
		/*** VALIDATE y0 AS POINT.*/
		/**************************/
		Ix0.real.a = Ix0.real.b = x0.real;
		Ix0.imag.a = Ix0.imag.b = 0.0;
		Iy0.real.a = Iy0.real.b = y0.real;
		Iy0.imag.a = Iy0.imag.b = y0.imag;
		IY0 = Iy0;
		IY0.real.a -= eps;
		IY0.real.b += eps;
		IY0.imag.a -= eps;
		IY0.imag.b += eps;
		while (validatePoint (If, Ix0, Iy0, IY0) == 0) {
		
#ifdef PROFILING
	nrejectedEps++; // used for profiling
#endif
			eps *= 0.5;
			IY0 = Iy0;
			IY0.real.a -= eps;
			IY0.real.b += eps;
			IY0.imag.a -= eps;
			IY0.imag.b += eps;
		}
		
		
		// COMPUTE a = -fx(x0,y0) / fy(x0.y0)
		// both fp and Interval
		getA (&Ia, f, x0, y0);
		a.real = Ia.real.a; a.imag = Ia.imag.a;

		/************************/
		/* ESTIMATE STEP SIZE ***/
		/************************/
		// we have: aux1 = fy, aux2 = fyy
		// store all in aux2
		Cmul (&aux2, aux2, a);
		Cmul (&aux2, aux2, a);
		
		evalPolXY (&aux, f, x0, y0);
		Cmul (&aux, aux, a);
		aux.real *= 2.0; aux.imag *= 2.0;
		Cadd (&aux2, aux2, aux);

		evalPolXX (&aux, f, x0, y0);
		Cadd (&aux2, aux2, aux);
		
		Cdiv (&aux, aux2, aux1);
		stepsize = sqrt (aux.real*aux.real + aux.imag*aux.imag);
		// DETECT INFLEXION POINT
		if (stepsize < 1.0e-10)
			stepsize=1.0;
		else 
			stepsize = sqrt(eps / stepsize);
		if (stepsize + x0.real > 1.0) stepsize = 1.000001 - x0.real;
		
		
		/*******************/	
		/** VALIDATE STEP  */
		/*******************/
		Iy1 = Iy0;
		IY1 = Iy1;
		IY1.real.a -= eps;
		IY1.real.b += eps;
		IY1.imag.a -= eps;
		IY1.imag.b += eps;

		// PERFORM TRANSLATION AND ROTATION OF THE POLYNOMIAL
		IlinearVarChange2 (Ig, If, Ix0);
		IlinearVarChange (Ih, Ig, Ia);
		
		Ix0.real.b = Ix0.real.a + stepsize;

		Ix1.imag.a = Ix1.imag.b = 0.0;
		Ix1.real.a = 0.0;
		Ix1.real.b = stepsize;	  

		while (validateStep (Ih, Ix1, Iy1, IY1) == 0) {
#ifdef PROFILING
	nrejectedSteps++; // used for profiling
#endif
			stepsize *= 0.5;
			eps*=0.95;
			IY1 = Iy1;
			IY1.real.a -= eps;
			IY1.real.b += eps;
			IY1.imag.a -= eps;
			IY1.imag.b += eps;
			Ix1.real.b = stepsize;
			if (stepsize < 1.0e-13) {
				/*******************/
				/* CLEAN VARIABLES */
				/*******************/
				//deleteList (&outputList);
				freeCombinatorialNumbers ();
				free (f->coef);
				free (If->coef);
				free (Ig->coef);
				free (Ih->coef);
				free (f);
				free (If);
				free (Ig);
				free (Ih);
				return NULL;
		}
		
#ifdef DEVELOPER
	printf ("x0 = %e\tstep = %e\n", x0.real, stepsize);
#endif



		/**************************/	
		/** VALIDATE FINAL STEP   */
		/**************************/
		if (x0.real + stepsize > 1.0) {
			stepsize = 1.0 - x0.real;
			Ix1.real.b = stepsize;
			
#ifdef DEVELOPER
	printf ("step final = %e\n", stepsize);
#endif

			while (validateStep (Ih, Ix1, Iy1, IY1) == 0);  
		}		
			

		/********************/
		/* UPDATE VARIABLES */
		/********************/
		x0.real = Ix0.real.a + stepsize;
		y0.real = y0.real + Ia.real.a * stepsize;
		y0.imag = y0.imag + Ia.imag.a * stepsize;
		
		
		// STORE LAST VALIDATED BOX TO CHECK WE HAVE NOT JUMPED TO ANOTHER THREAD
		IlastY0 = Iy0;
		IlastY0.real.a += Ia.real.a * stepsize;
		IlastY0.real.b += Ia.real.b * stepsize;
		IlastY0.imag.a += Ia.imag.a * stepsize;
		IlastY0.imag.b += Ia.imag.b * stepsize;
		IlastY0.real.a -= eps; IlastY0.real.b += eps;
		IlastY0.imag.a -= eps; IlastY0.imag.b += eps;

		

		/*******************/
		/* CORRECT ROOT	   */
		/*******************/
		correctRoot (f, &x0, &y0);

		/**********************************************/
		/* CHECK WE HAVE NOT JUMPED TO ANOTHER ROOT   */
		/**********************************************/
		if (!isContained (y0.real, IlastY0.real) || 
					!isContained (y0.imag, IlastY0.imag)) {
			printf ("error! Jumped to other thread!\n");
			/*******************/
			/* CLEAN VARIABLES */
			/*******************/
			//deleteList (&outputList);
			freeCombinatorialNumbers ();
			free (f->coef);
			free (If->coef);
			free (Ig->coef);
			free (Ih->coef);
			free (f);
			free (If);
			free (Ig);
			free (Ih);
			return NULL;
		}
		
		/*********************************/
		/* APPEND VALIDATED DATA TO LIST */
		/*********************************/
		dataToList[0] = x0.real; dataToList[1] = y0.real; dataToList[2] = y0.imag;
		appendData (&outputList, dataToList);
		
		if (x0.real + stepsize > 1.0) stepsize = 0.0;
		
		
	}
#ifdef PROFILING
	printf ("x0 = %e, y0 = %.15le%+.15le i\n", x0.real, y0.real, y0.imag);
	printf ("\t Number of steps		  = %i\n", nsteps);
	printf ("\t Number of rejected Eps= %i\n", nrejectedEps);
	printf ("\t Number of rejected steps = %i\n", nrejectedSteps);
#endif
	
	
	/***********************************************/
	/* PREPARE OUTPUT AS A BINARY ARRAY OF DOUBLES */
	/***********************************************/
	
	double *rop = malloc ((3 * outputList.len + 1) * sizeof (double));
	rop[0] = outputList.len;
	node_t *index = outputList.first;
	for (i=0; i<outputList.len; i++) {
		rop[3*i+1] = index->x[0];
		rop[3*i+2] = index->x[1];
		rop[3*i+3] = index->x[2];
		index = index->next;
	}

	/*******************/
	/* CLEAN VARIABLES */
	/*******************/
	//deleteList (&outputList);
	freeCombinatorialNumbers ();
	free (f->coef);
	free (If->coef);
	free (Ig->coef);
	free (Ih->coef);
	free (f);
	free (If);
	free (Ig);
	free (Ih);

	/*****************************************************************/
	/* OUTPUT VALIDATED PIECEWISE LINEAR APPROXIMATION OF THE THREAD */
	/*****************************************************************/	
	return rop;
}


