#ifndef _list_h
#define _list_h

typedef struct _node_t {
	double x[3];
	struct _node_t *next;

} node_t;

typedef struct _list {
	unsigned len;
	node_t *first;
} list_t;


void appendData (list_t *list, double x[3]);

void deleteFirstElement (list_t *list);

void deleteLastElement (list_t *list);

void deleteList (list_t *list);

void printfList (list_t list, char *s);


#endif
