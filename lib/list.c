#include <stdio.h>
#include <stdlib.h>
#include "list.h"


void appendData (list_t *list, double x[3]) {
    int i;
	node_t *newNode = (node_t *) malloc (sizeof (node_t));
	for (i=0; i<3; ++i) newNode->x[i] = x[i];    
	newNode->next = NULL;
	
	if (list->len == 0) {
		list->len = 1;
		list->first = newNode;
	} else {
		node_t *index = list->first;
		while (index->next != NULL) 
			index = index->next;
		index->next = newNode;
		++(list->len);
	}
	return;
}


void deleteFirstElement (list_t *list) {
	if (list->len == 0) {
		printf ("Error. Empty list\n");
		return;
	}
	if (list->len == 1) {
		free (list->first);
		list->first = NULL;
		list->len = 0;
		return;
	}
	
	node_t *first = list->first->next;
	free (list->first);
	list->first = first;
	--(list->len);
	return;
}

void deleteLastElement (list_t *list) {
	if (list->len == 0) {
		printf ("Error. Lista vacia\n");
		return;
	}
	if (list->len == 1) {
		free (list->first);
		list->first = NULL;
		list->len = 0;
	}
	node_t *index = list->first;
	while (index->next->next != NULL) 
		index = index->next;
	free (index->next);
	index->next = NULL;
	--(list->len);
	return;
	
}

void deleteList (list_t *list) {
	while (list->first != NULL)
		deleteFirstElement (list);
	return;
}


void printfList (list_t list, char *s) {
	int i;
	if (list.len == 0) {
		printf ("%s = {}\n", s);
		return;
	}
	
	if (list.len < 3) {
		node_t *index = list.first;
		printf ("%s = {%.5e, %.5e, %.5e", s, index->x[0], index->x[1], index->x[2]);
		while (index->next != NULL) {
			index = index->next;
			printf (", %.5e, %.5e, %.5e", index->x[0], index->x[1], index->x[2]);
		}
		printf ("}\n");
		return;
	}
	
	node_t *index = list.first;
	printf ("%s = {%.5e, %.5e, %.5e", s, index->x[0], index->x[1], index->x[2]);
	for (i=0; i<2; i++) {
		index = index->next;
		printf (", %.5e, %.5e, %.5e", index->x[0], index->x[1], index->x[2]);
	}
	while (index->next != NULL) {
		index = index->next;
	}
	printf (", ... , %.5e, %.5e, %.5e", index->x[0], index->x[1], index->x[2]);
	return;


}


/*int main () {
	list_t list = {0, NULL};
	
	appendData (&list, (double[3]) {2,2.1,2.2}); appendData (&list, (double[3]) {4,4.1,4.2}); 
	appendData (&list, (double[3]) {6,6.1,6.2}); appendData (&list, (double[3]) {8,8.1,8.2}); 
	appendData (&list, (double[3]) {10,10.1,10.2}); appendData (&list, (double[3]) {12,12.1,12.2}); 
	appendData (&list, (double[3]) {14,14.1,14.2}); appendData (&list, (double[3]) {16,16.1,16.2}); 	
	printfList (list, "l");
	
	deleteList (&list);
	printfList (list, "l");
	appendData (&list, (double[3]) {1,1.1,1.2}); appendData (&list, (double[3]) {3,3.1,3.2});
	appendData (&list, (double[3]) {5,5.1,5.2}); appendData (&list, (double[3]) {7,7.1,7.2}); 
	printfList (list, "l");
	deleteFirstElement (&list);
	printfList (list, "l");
	deleteLastElement (&list);
	printfList (list, "l");

}*/

