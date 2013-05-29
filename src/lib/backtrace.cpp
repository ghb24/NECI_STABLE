// Copyright (c) 2013, Ali Alavi
// This program is integrated in Molpro with the permission of George Booth and Ali Alavi
  
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

#define lenof(x) ((sizeof((x)))/(sizeof((x)[0])))

extern "C" void print_backtrace_neci ()
{
	void * buf[30];
	int n = backtrace (buf, lenof(buf));
	char ** strs = backtrace_symbols (buf, n);

	printf ("-----------------------------------\n");
	printf ("Writing Backtrace\n");
	printf ("-----------------------------------\n");
	// n.b. We don't include the function 'print_backtrace' in the backtrace.
	for (int i = 1; i < n; i++)
		printf("%d: %s\n", i-1, strs[i]);
	printf ("-----------------------------------\n");
	free (strs);
}

