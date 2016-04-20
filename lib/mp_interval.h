#ifndef _mp_interval_t_H
#define _mp_interval_t_H


typedef struct {
	mpfr_t a;
	mpfr_t b;
} mp_interval_t;


/********************
 MISCELLANEOUS FUNCTIONS
 **********************/
 
// COMPUTES INTERSECTION. returns 1 IF NOT EMPTY
int mp_IntersectInterval (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2);

// PRINTS INTERVAL
void mp_Iprintf (mp_interval_t x, char *s);
 


/*************************
 BOOLEAN FUNCTIONS
 **********************/

// TRUE IF ZERO IS CONTAINED IN OP
int mp_isZeroContained (mp_interval_t op);

// TRUE IF VALUE OP1 IS CONTAINED IN interval_t OP2
int mp_isContained (mpfr_t op1, mp_interval_t op2);			
	
// TRUE IF interval_t OP1 IS EQUAL TO interval_t OP2
int mp_isEqual (mp_interval_t op1, mp_interval_t op2);			

// TRUE IF interval_t OP1 IS SUBSET NOT EQUAL TO interval_t OP2
int mp_isSubset (mp_interval_t op1, mp_interval_t op2);				

// TRUE IF interval_t OP1 IS SUBSET OR EQUAL TO interval_t OP2
int mp_isSubsetEqual (mp_interval_t op1, mp_interval_t op2);		


// ARITHMETIC FUNCTIONS
void mp_Iadd (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2);		// ROP = OP1 + OP2
void mp_IaddC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2);		// ROP = OP1 + OP2	OP2 IS mpfr_t
void mp_Isub (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2);		// ROP = OP1 - OP2
void mp_IsubC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2);		// ROP = OP1 - OP2	OP2 IS mpfr_t
void mp_ICsub (mp_interval_t *rop, mpfr_t op1, mp_interval_t op2);		// ROP = OP1 - OP2	OP1 IS mpfr_t
void mp_Imul (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2);		// ROP = OP1 * OP2
void mp_ImulC (mp_interval_t *rop, mp_interval_t op1, mpfr_t op2);		// ROP = OP1 * OP2	OP2 IS mpfr_t
void mp_Iinv (mp_interval_t *rop, mp_interval_t op);				// ROP = 1 / OP



// TO DO
/*void mp_Idiv (mp_interval_t *rop, mp_interval_t op1, mp_interval_t op2);		// ROP = OP1 / OP2
void mp_IinvC_d (mp_interval_t *rop, mpfr_t op);				// ROP = 1 / OP		OP IS mpfr_t
void mp_Isqr (mp_interval_t *rop, mp_interval_t op);				// ROP = OP^2
void Ilog (mp_interval_t *rop, mp_interval_t op);				// ROP = log (OP)
void IlogC_d (mp_interval_t *rop, mpfr_t op);				// ROP = log (OP)	OP IS mpfr_t*/

#endif

