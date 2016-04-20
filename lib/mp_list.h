#ifndef _mp_list_h
#define _mp_list_h

typedef struct _mp_node_t {
	mpfr_t x;
	struct _mp_node_t *next;

} mp_node_t;

typedef struct _mp_list {
	unsigned len;
	mp_node_t *first;
} mp_list_t;


void mp_appendData (mp_list_t *list, mpfr_t x);

void mp_deleteFirstElement (mp_list_t *list);

void mp_deleteLastElement (mp_list_t *list);

void mp_deleteList (mp_list_t *list);

void mp_printfList (mp_list_t list, char *s);


#endif
