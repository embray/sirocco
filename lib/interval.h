#ifndef _interval_t_H
#define _interval_t_H


typedef struct {
	volatile double a;
	volatile double b;
} interval_t;


/**********************
 MISCELLANEOUS FUNCTIONS
 ************************/

// COMPUTES INTERSECTION. returns 1 IF NOT EMPTY
int IntersectInterval (interval_t *rop, interval_t op1, interval_t op2);

// PRINTS interval_t [a,b]
void Iprintf (interval_t x, char *s);


/*************************
 BOOLEAN FUNCTIONS
 **********************/

// TRUE IF ZERO IS CONTAINED IN OP
int isZeroContained (interval_t op);

// TRUE IF VALUE OP1 IS CONTAINED IN interval_t OP2
int isContained (double op1, interval_t op2);			
	
// TRUE IF interval_t OP1 IS EQUAL TO interval_t OP2
int isEqual (interval_t op1, interval_t op2);			

// TRUE IF interval_t OP1 IS SUBSET NOT EQUAL TO interval_t OP2
int isSubset (interval_t op1, interval_t op2);				

// TRUE IF interval_t OP1 IS SUBSET OR EQUAL TO interval_t OP2
int isSubsetEqual (interval_t op1, interval_t op2);			


/*******************
 ARITHMETIC FUNCTIONS
 ********************/
// ROP = OP1 + OP2
void Iadd (interval_t *rop, interval_t op1, interval_t op2);	

// ROP = OP1 + OP2	OP2 IS DOUBLE	
void IaddC (interval_t *rop, interval_t op1, double op2);		

// ROP = OP1 - OP2
void Isub (interval_t *rop, interval_t op1, interval_t op2);	

// ROP = OP1 - OP2	OP2 IS DOUBLE
void IsubC (interval_t *rop, interval_t op1, double op2);	

// ROP = OP1 - OP2	OP1 IS DOUBLE
void ICsub (interval_t *rop, double op1, interval_t op2);	

// ROP = OP1 * OP2
void Imul (interval_t *rop, interval_t op1, interval_t op2);	

// ROP = OP1 * OP2	OP2 IS DOUBLE
void ImulC (interval_t *rop, interval_t op1, double op2);		

// ROP = 1 / OP
void Iinv (interval_t *rop, interval_t op);				

// ROP = 1 / OP		OP IS DOUBLE
void IinvC (interval_t *rop, double op);				

// ROP = OP1 / OP2
void Idiv (interval_t *rop, interval_t op1, interval_t op2);	

// ROP = OP^2
void Isqr (interval_t *rop, interval_t op);

#endif

