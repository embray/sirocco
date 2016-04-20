#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "mp_list.h"


void mp_appendData (mp_list_t *list, mpfr_t x) {
    
	mp_node_t *newNode = (mp_node_t *) malloc (sizeof (mp_node_t));
	mpfr_init (newNode->x);
	mpfr_set (newNode->x, x, MPFR_RNDN);    
	newNode->next = NULL;
	
	if (list->len == 0) {
		list->len = 1;
		list->first = newNode;
	} else {
		mp_node_t *index = list->first;
		while (index->next != NULL) 
			index = index->next;
		index->next = newNode;
		++(list->len);
	}
	return;
}


void mp_deleteFirstElement (mp_list_t *list) {
	if (list->len == 0) {
		printf ("Error. Empty list\n");
		return;
	}
	if (list->len == 1) {
		mpfr_clear (list->first->x);
		free (list->first);
		list->first = NULL;
		list->len = 0;
		return;
	}
	
	mp_node_t *first = list->first->next;
	mpfr_clear (list->first->x);
	free (list->first);
	list->first = first;
	--(list->len);
	return;
}

void mp_deleteLastElement (mp_list_t *list) {
	if (list->len == 0) {
		printf ("Error. Lista vacia\n");
		return;
	}
	if (list->len == 1) {
		mpfr_clear (list->first->x);
		free (list->first);
		list->first = NULL;
		list->len = 0;
	}
	mp_node_t *index = list->first;
	while (index->next->next != NULL) 
		index = index->next;
	mpfr_clear (index->next->x);
	free (index->next);
	index->next = NULL;
	--(list->len);
	return;
	
}

void mp_deleteList (mp_list_t *list) {
	while (list->first != NULL)
		mp_deleteFirstElement (list);
	return;
}


void mp_printfList (mp_list_t list, char *s) {
	int i;
	if (list.len == 0) {
		printf ("%s = {}\n", s);
		return;
	}
	
	if (list.len < 6) {
		mp_node_t *index = list.first;
		printf ("%s = {", s);
		mpfr_out_str (stdout, 10, 8, index->x, MPFR_RNDN);
		while (index->next != NULL) {
			index = index->next;
			printf (", ");
			mpfr_out_str (stdout, 10, 8, index->x, MPFR_RNDN);
		}
		printf ("}\n");
		return;
	}
	
	mp_node_t *index = list.first;
	printf ("%s = {", s);
	mpfr_out_str (stdout, 10, 8, index->x, MPFR_RNDN);
	for (i=0; i<3; i++) {
		index = index->next;
		printf (", ");
		mpfr_out_str (stdout, 10, 8, index->x, MPFR_RNDN);
	}
	while (index->next != NULL) {
		index = index->next;
	}
	printf (", ... , ");
	mpfr_out_str (stdout, 10, 8, index->x, MPFR_RNDN);
	printf ("}\n");
	return;


}

void mp_appendData_i (mp_list_t *list, int x) {
	mpfr_t data; mpfr_init (data);
	mpfr_set_si (data, x, MPFR_RNDN);
	mp_appendData (list, data);	
	mpfr_clear (data);
}

int main () {
	mp_list_t list = {0, NULL};

	mp_appendData_i (&list, 2);  mp_appendData_i (&list, 4);	 mp_appendData_i (&list, 6); 
	mp_appendData_i (&list, 8);  mp_appendData_i (&list, 10); mp_appendData_i (&list, 12); 
	mp_appendData_i (&list, 14); mp_appendData_i (&list, 16);

	mp_printfList (list, "l");
	mp_deleteList (&list);
	mp_printfList (list, "l");
	mp_appendData_i (&list, 1); mp_appendData_i (&list, 3); mp_appendData_i (&list, 5); mp_appendData_i (&list, 7); 
	mp_printfList (list, "l");
	mp_deleteFirstElement (&list);
	mp_printfList (list, "l");
	mp_deleteLastElement (&list);
	mp_printfList (list, "l");
	mp_deleteList (&list);

}

