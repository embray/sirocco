#include <stdio.h>
#include <fenv.h>
#include <math.h>
#include "interval.h"


#define MIN(a,b) (((a) < (b))? (a) : (b))
#define MAX(a,b) (((a) > (b))? (a) : (b))

int IntersectInterval (interval_t *rop, interval_t op1, interval_t op2) {
	rop->a = MAX (op1.a, op2.a);
	rop->b = MIN (op1.b, op2.b);
	if (rop->a < rop->b) return 1;
	return 0;
}

int isZeroContained (interval_t op) {
	if(op.a <= 0 && op.b >= 0) return 1;
	return 0;
}

int isContained (double op1, interval_t op2) {
	if (op1 >= op2.a && op1 <= op2.b) return 1;
	return 0;
}

int isEqual (interval_t op1, interval_t op2) {
	if (op1.a == op2.a && op1.b == op2.b) return 1;
	return 0;
}

int isSubset (interval_t op1, interval_t op2) {
	if (op1.a >= op2.a && op1.b <= op2.b && !isEqual (op1, op2)) return 1;
	return 0;
}

int isSubsetEqual (interval_t op1, interval_t op2) {
	if (op1.a >= op2.a && op1.b <= op2.b) return 1;
	return 0;
}


void Iadd (interval_t *rop, interval_t op1, interval_t op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = op1.a + op2.a;

	fesetround (FE_UPWARD);
	rop->b = op1.b + op2.b;

	fesetround (oldRoundingMode);
}

void IaddC (interval_t *rop, interval_t op1, double op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = op1.a + op2;

	fesetround (FE_UPWARD);
	rop->b = op1.b + op2;

	fesetround (oldRoundingMode);
}

void Isub (interval_t *rop, interval_t op1, interval_t op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = op1.a - op2.b;

	fesetround (FE_UPWARD);
	rop->b = op1.b - op2.a;

	fesetround (oldRoundingMode);
}

void IsubC (interval_t *rop, interval_t op1, double op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = op1.a - op2;

	fesetround (FE_UPWARD);
	rop->b = op1.b - op2;

	fesetround (oldRoundingMode);
}

void ICsub (interval_t *rop, double op1, interval_t op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = op1 - op2.b;

	fesetround (FE_UPWARD);
	rop->b = op1 - op2.a;

	fesetround (oldRoundingMode);
}

void Imul (interval_t *rop, interval_t op1, interval_t op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = MIN (MIN (MIN (op1.a * op2.a, op1.a * op2.b), op1.b * op2.a), op1.b * op2.b);

	fesetround (FE_UPWARD);
	rop->b = MAX (MAX (MAX (op1.a * op2.a, op1.a * op2.b), op1.b * op2.a), op1.b * op2.b);

	fesetround (oldRoundingMode);
}


void ImulC (interval_t *rop, interval_t op1, double op2) {
	int oldRoundingMode = fegetround ();

	fesetround (FE_DOWNWARD);
	rop->a = MIN (op1.a * op2, op1.b * op2);

	fesetround (FE_UPWARD);
	rop->b = MAX (op1.a * op2, op1.b * op2);

	fesetround (oldRoundingMode);
}

void Iinv (interval_t *rop, interval_t op) {
	int oldRoundingMode = fegetround ();

	if (isZeroContained (op)) {
		//printf ("Divided by zero. Bad Result!\n");
		if (op.a == 0) {
			fesetround (FE_DOWNWARD);
			rop->a = 1.0 / op.b;
			rop->b = +INFINITY;
			return;
		}
		if (op.b == 0) {
			fesetround (FE_UPWARD);
			rop->b = 1.0 / op.a;
			rop->a = -INFINITY;
			return;
		}
		rop->a = -INFINITY;
		rop->b = +INFINITY;
		return;
	}
	
	fesetround (FE_DOWNWARD);
	rop->a = 1.0 / op.b;

	fesetround (FE_UPWARD);
	rop->b = 1.0 / op.a;

	fesetround (oldRoundingMode);
}

void IinvC (interval_t *rop, double op) {
	if (op == 0) {
		rop->a = NAN;
		rop->b = NAN;
		return;
	}

	int oldRoundingMode = fegetround ();

	fesetround (FE_UPWARD);
	rop->a = 1.0 / op;
	fesetround (FE_DOWNWARD);
	rop->a = 1.0 / op;

	fesetround (oldRoundingMode);
}

void Idiv (interval_t *rop, interval_t op1, interval_t op2) {
	interval_t invOp2;
	Iinv (&invOp2, op2);
	Imul (rop, op1, invOp2);
}


void Isqr (interval_t *rop, interval_t op) {
	int oldRoundingMode = fegetround ();

	if (isContained (0, op)) {
		rop->a = 0.0;
		fesetround (FE_UPWARD);
		rop->b = MAX (op.a * op.a, op.b * op.b);
	}
	else {
		fesetround (FE_DOWNWARD);
		rop->a = MIN (op.a * op.a, op.b * op.b);

		fesetround (FE_UPWARD);
		rop->b = MAX (op.a * op.a, op.b * op.b);
	}

	fesetround (oldRoundingMode);
}


void Iprintf (interval_t x, char *s) {
	printf ("%s = [%.16le, %.16le]\n", s, x.a, x.b);
}



/*

int main () {
	interval_t x, y, z;
	double c;
	
	x.a = -1.999999999999999999999999999;
	x.b = 2.999999999999999999999999999;
	y.a = 2.33;
	y.b = 3.66;
	c = -1.35;
	
	printf ("---------------------\n");
	Iprintf (x, "x  "); printf ("---------------------\n");
	Iprintf (y, "y  "); printf ("---------------------\n");
	printf ("c   = -1.35\n"); printf ("---------------------\n");
	
	
	Idiv (&z, x, y);
	Iprintf (z, "x/y"); printf ("---------------------\n");
	Idiv (&z, y, x);
	Iprintf (z, "y/x"); printf ("---------------------\n");

}
*/

