#include <stdio.h>
#include <mpfr.h>
#include "mp_interval.h"


int mp_IntersectInterval (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2) {
	//rop->a = MAX (op1.a, op2.a);
	mpfr_max (rop->a, op1.a, op2.a, MPFR_RNDD);
	
	//rop->b = MIN (op1.b, op2.b);
	mpfr_min (rop->b, op1.b, op2.b, MPFR_RNDU);
	
	//if (rop->a < rop->b) return 1;
	if (mpfr_cmp (rop->a, rop->b) < 0) return 1;
	return 0;
}


int mp_isZeroContained (mp_interval_t op2) {
	/* CHECKS WHETHER ZERO IS CONTAINED IN INTERVAL op2*/
	
	// CHECKS IF ZERO IS LOWER THAN LEFT BOUND
	if (mpfr_cmp_si (op2.a, 0L) > 0) return 0;
	// CHECKS IF ZERO IS GREATER THAN RIGHT BOUND
	if (mpfr_cmp_si (op2.b, 0L) < 0) return 0;
	
	// OTHERWISE IS INSIDE
	return 1;
}

int mp_isContained (mpfr_t op1, mp_interval_t op2) {
	//if (op1 >= op2.a && op1 <= op2.b) return 1;
	if (mpfr_cmp (op1, op2.a) >= 0 && mpfr_cmp (op1, op2.b) <= 0) return 1;
	return 0;
}


int mp_isEqual (mp_interval_t op1, mp_interval_t op2) {
	// if (op1.a == op2.a && op1.b == op2.b) return 1;
	if(mpfr_cmp (op1.a, op2.a) == 0 && mpfr_cmp (op1.b, op2.b) == 0) return 1;
	return 0;
}


int mp_isSubset (mp_interval_t op1, mp_interval_t op2) {
	//if (op1.a >= op2.a && op1.b <= op2.b && !isEqual (op1, op2)) return 1;
	if (mpfr_cmp (op1.a, op2.a) >= 0 && mpfr_cmp (op1.b, op2.b) <= 0 && !mp_isEqual (op1, op2)) return 1;
	return 0;
}

int mp_isSubsetEqual (mp_interval_t op1, mp_interval_t op2) {
	// if (op1.a >= op2.a && op1.b <= op2.b) return 1;
	if (mpfr_cmp (op1.a, op2.a) >= 0 && mpfr_cmp (op1.b, op2.b) <= 0) return 1;
	return 0;
}



void mp_Iadd (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2) {
	// LEFT BOUNDARIES ADD, ROUNDING DOWNWARDS
	mpfr_add (rop->a, op1.a, op2.a, MPFR_RNDD);
	// RIGHT BOUNDARIES ADD, ROUNDING UPWARDS
	mpfr_add (rop->b, op1.b, op2.b, MPFR_RNDU);
}

void mp_IaddC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2) {
	// LEFT BOUNDARY, ROUNDING DOWNWARDS
	mpfr_add (rop->a, op1.a, op2, MPFR_RNDD);

	// RIGHT BOUNDARY, ROUNDING UPWARDS
	mpfr_add (rop->b, op1.b, op2, MPFR_RNDU);
}

void mp_Isub (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2) {
	// LEFT BOUNDARY, ROUNDING DOWNWARDS
	mpfr_sub (rop->a, op1.a, op2.b, MPFR_RNDD);
	// RIGHT BOUNDARY, ROUNDING UPWARDS
	mpfr_sub (rop->b, op1.b, op2.a, MPFR_RNDU);
}

void mp_IsubC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2) {
	// LEFT BOUNDARY, ROUNDING DOWNWARDS
	mpfr_sub (rop->a, op1.a, op2, MPFR_RNDD);
	// RIGHT BOUNDARY, ROUNDING UPWARDS
	mpfr_sub (rop->b, op1.b, op2, MPFR_RNDU);
}

void mp_ICsub (mp_interval_t *rop, mpfr_t op1, mp_interval_t op2) {
	// LEFT BOUNDARY, ROUNDING DOWNWARDS
	mpfr_sub (rop->a, op1, op2.b, MPFR_RNDD);
	// RIGHT BOUNDARY, ROUNDING UPWARDS
	mpfr_sub (rop->b, op1, op2.a, MPFR_RNDU);
}

void mp_Imul (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2) {
	// DECLARE AUX MPFR_T VARIABLE TO STORE CROSSED PRODUCTS
	mpfr_t aux; mpfr_init (aux);
	mpfr_t ropA; mpfr_init (ropA);
	mpfr_t ropB; mpfr_init (ropB);
	// FIRST CANDIDATE TO LEFT BOUNDARY
	mpfr_mul (ropA, op1.a, op2.a, MPFR_RNDD);
		mpfr_mul (aux, op1.a, op2.b, MPFR_RNDD); mpfr_min (ropA, ropA, aux, MPFR_RNDD);
		mpfr_mul (aux, op1.b, op2.a, MPFR_RNDD); mpfr_min (ropA, ropA, aux, MPFR_RNDD);
		mpfr_mul (aux, op1.b, op2.b, MPFR_RNDD); mpfr_min (ropA, ropA, aux, MPFR_RNDD);


	// FIRST CANDIDATE TO RIGHT BOUNDARY
	mpfr_mul (ropB, op1.a, op2.a, MPFR_RNDU);
		mpfr_mul (aux, op1.a, op2.b, MPFR_RNDU); mpfr_max (ropB, ropB, aux, MPFR_RNDU);
		mpfr_mul (aux, op1.b, op2.a, MPFR_RNDU); mpfr_max (ropB, ropB, aux, MPFR_RNDU);
		mpfr_mul (aux, op1.b, op2.b, MPFR_RNDU); mpfr_max (rop->b, ropB, aux, MPFR_RNDU);

	mpfr_set (rop->a, ropA, MPFR_RNDD);

	// CLEAN AUX VARIABLE
	mpfr_clear (aux);
	mpfr_clear (ropA);
	mpfr_clear (ropB);
}


void mp_ImulC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2) {
	// DECLARE AUX MPFR_T VARIABLE TO STORE CROSSED PRODUCTS
	mpfr_t aux; mpfr_init (aux);
	mpfr_t ropA; mpfr_init (ropA);
	mpfr_t ropB; mpfr_init (ropB);
	
	// FIRST CANDIDATE TO LEFT BOUNDARY
	mpfr_mul (ropA, op1.a, op2, MPFR_RNDD);
		mpfr_mul (aux, op1.b, op2, MPFR_RNDD); mpfr_min (ropA, ropA, aux, MPFR_RNDD);

	// FIRST CANDIDATE TO RIGHT BOUNDARY
	mpfr_mul (ropB, op1.a, op2, MPFR_RNDU);
		mpfr_mul (aux, op1.b, op2, MPFR_RNDU); mpfr_max (rop->b, ropB, aux, MPFR_RNDU);


	mpfr_set (rop->a, ropA, MPFR_RNDD);
	
	// CLEAN AUX VARIABLE
	mpfr_clear (aux);
	mpfr_clear (ropA);
	mpfr_clear (ropB);
}




void mp_Iinv (mp_interval_t *rop, mp_interval_t op) {
	// CHECK IF OP = [0,0]
	if (mpfr_cmp_si (op.a, 0) == 0 && mpfr_cmp_si (op.b, 0) == 0) {
		mpfr_set_str (rop->a, "nan", 10, MPFR_RNDD);
		mpfr_set_str (rop->b, "nan", 10, MPFR_RNDU);
		return;
	}

	// CHECKS WHETHER ZERO IS CONTAINED IN OP
	if (mp_isZeroContained (op)) {
		// IF ZERO IS LEFT BOUNDARY
		if (mpfr_cmp_si (op.a, 0) == 0) {
			mpfr_si_div (rop->a, 1, op.b, MPFR_RNDD);
			mpfr_set_str (rop->b, "+inf", 10, MPFR_RNDU);
			return;
		}
		// IF ZERO IS RIGHT BOUNDARY
		if (mpfr_cmp_si (op.b, 0) == 0) {
			mpfr_si_div (rop->b, 1, op.a, MPFR_RNDU);
			mpfr_set_str (rop->a, "-inf", 10, MPFR_RNDD);
			return;
		} 
		// ZERO IS EXTRICTLY CONTAINED IN op
		mpfr_set_str (rop->a, "-inf", 10, MPFR_RNDD);
		mpfr_set_str (rop->b, "+inf", 10, MPFR_RNDU);
		return;
		
	} 
	// ZERO IS NOT IN op
	mpfr_t ropA; mpfr_init (ropA);
	
	mpfr_si_div (ropA, 1, op.b, MPFR_RNDD);
	mpfr_si_div (rop->b, 1, op.a, MPFR_RNDU);	

	mpfr_set (rop->a, ropA, MPFR_RNDD);
	
	mpfr_clear (ropA);
}

void mp_IinvC (mp_interval_t *rop, mpfr_t op) {
	// CHECKS WHETHER op is 0
	if (mpfr_cmp_si (op, 0) == 0) {
		mpfr_set_str (rop->a, "nan", 10, MPFR_RNDD);
		mpfr_set_str (rop->b, "nan", 10, MPFR_RNDU);
		return;
	}
	// UP TO HERE, op is NOT ZERO
	mpfr_si_div (rop->a, 1, op, MPFR_RNDD);
	mpfr_si_div (rop->b, 1, op, MPFR_RNDU);	
}



void mp_Idiv (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2) {
	// AUXILIAR mp_interval_t VARIABLE TO STORE INVERSE OF op2
	mp_interval_t invOp2; mpfr_init (invOp2.a); mpfr_init (invOp2.b);
	
	mp_Iinv (&invOp2, op2);
	mp_Imul (rop, op1, invOp2);
	
	// CLEAN AUXILIAR VARIABLE
	mpfr_clear (invOp2.a); mpfr_clear (invOp2.b);
}


void mp_Iprintf (mp_interval_t x, char *s) {
	printf ("%s = [", s);
	mpfr_out_str (stdout, 10, 16, x.a, MPFR_RNDD);
	printf (", ");
	mpfr_out_str (stdout, 10, 16, x.b, MPFR_RNDU);
	printf ("]\n");
}

int main () {
	mpfr_set_default_prec (60);
	mp_interval_t x, y, z;
	mpfr_t c;
	mpfr_init (x.a);	mpfr_init (x.b);
	mpfr_init (y.a);	mpfr_init (y.b);
	mpfr_init (z.a);	mpfr_init (z.b);
	mpfr_init (c);
	
	mpfr_set_str (x.a, "-1.99999999999999999999999999999999999999999999999999999999999999", 10, MPFR_RNDD);
	mpfr_set_str (x.b, "2.99999999999999999999999999999999999999999999999999999999999999", 10, MPFR_RNDU);
	mpfr_set_str (y.a, "2.33", 10, MPFR_RNDD);
	mpfr_set_str (y.b, "3.66", 10, MPFR_RNDD);
	mpfr_set_str (c, "-1.35", 10, MPFR_RNDN);
	
	printf ("---------------------\n");
	mp_Iprintf (x, "x  "); printf ("---------------------\n");
	mp_Iprintf (y, "y  "); printf ("---------------------\n");
	printf ("c   = -1.35\n"); printf ("---------------------\n");
	
	
	mp_Idiv (&z, x, y);
	mp_Iprintf (z, "x/y"); printf ("---------------------\n");
	mp_Idiv (&z, y, x);
	mp_Iprintf (z, "y/x"); printf ("---------------------\n");

}





/*
void IstreamPrint (FILE *fout, interval_t op) {
	fprintf (fout, "[%.17le, %.17le]", op.a, op.b);
	return;
}

void IstreamPrintln (FILE *fout, interval_t op) {
	fprintf (fout, "[%.17le, %.17le]\n", op.a, op.b);
	return;
}

int isContained (double op1, interval_t op2) {
	return (op1 >= op2.a && op1 <= op2.b);
}

int isEqual (interval_t op1, interval_t op2) {
	return (op1.a == op2.a && op1.b == op2.b);
}

int isSubset (interval_t op1, interval_t op2) {
	return (op1.a >= op2.a && op1.b <= op2.b && !isEqual (op1, op2));
}

int isSubsetEqual (interval_t op1, interval_t op2) {
	return (op1.a >= op2.a && op1.b <= op2.b);
}*/



/*
void Isqr (mp_interval_t *rop, mp_interval_t op) {
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

void Ilog (mp_interval_t *rop, mp_interval_t op) {
	if (op.a < 0.0) {
		rop->a = rop->b = NAN;
	}
	else {
		int oldRoundingMode = fegetround ();
		if (op.a == 0.0) {
			rop->a = -INFINITY;
			fesetround (FE_UPWARD);
			rop->b = log (op.b);
		} else {
			fesetround (FE_DOWNWARD);
			rop->a = log (op.a);
			fesetround (FE_UPWARD);
			rop->b = log (op.b);
		}
		fesetround (oldRoundingMode);
	}
}

void IlogC (mp_interval_t *rop, double op) {
	if (op < 0.0) {
		rop->a = rop->b = NAN;
	}
	else if (op == 0.0) {
		rop->a = rop->b = -INFINITY;
	}
	else {
		int oldRoundingMode = fegetround ();
		fesetround (FE_DOWNWARD);
		rop->a = log (op);
		fesetround (FE_UPWARD);
		rop->b = log (op);

		fesetround (oldRoundingMode);
	}
}*/




