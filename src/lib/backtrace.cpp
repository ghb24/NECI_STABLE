// Copyright (c) 2013, Ali Alavi unless otherwise noted.
// This program is integrated in Molpro with the permission of George Booth and Ali Alavi
  
#ifdef _MOLCAS_
#include "molcas_wrapper.h"
#endif

#ifndef __CYGWIN__
#include <execinfo.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define lenof(x) ((sizeof((x)))/(sizeof((x)[0])))

#ifdef __cplusplus
extern "C"
#endif
void print_backtrace_neci ()
{
#ifndef __CYGWIN__
	void * buf[30];
	int n = backtrace (buf, lenof(buf));
	char ** strs = backtrace_symbols (buf, n);

	printf ("-----------------------------------\n");
	printf ("Writing Backtrace\n");
	printf ("-----------------------------------\n");
	// n.b. We don't include the function 'print_backtrace' in the backtrace.
	int i;
	for (i = 1; i < n; i++)
		printf("%d: %s\n", i-1, strs[i]);
	printf ("-----------------------------------\n");
	free (strs);
#endif
}


// Wrapper to make NAG happy
#ifdef __cplusplus
extern "C"
#endif
size_t strlen_wrap (const char * str )
{
	return strlen (str);
}


