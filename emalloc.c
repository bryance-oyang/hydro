/** @file
 * @brief error checking malloc
 */
#include "emalloc.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <signal.h>

/**
 * @brief error checking malloc
 *
 * also zeros out allocated memory
 *
 * @param b size of memory to be allocated
 * @return pointer to allocated memory
 */
void *emalloc(size_t b)
{
	void *p;

	if ((p = malloc(b)) == NULL) {
		fprintf(stderr, "%s\n", strerror(ENOMEM));
		fflush(stderr);
		raise(SIGSEGV);
	}
	memset(p, 0, b);

	return p;
}
