#include <signal.h>

#ifdef __cplusplus
extern "C"
#endif
void neci_sigint(int signo);

// Clear all set signals
#ifdef __cplusplus
extern "C"
#endif
void clear_signals() {
	signal(SIGINT, SIG_DFL);
}

// Initialise signals
// --> Catch SIGINT and refer it to the f90 signal handler
#ifdef __cplusplus
extern "C"
#endif
void init_signals_helper() {
	clear_signals();
	signal(SIGINT, neci_sigint);
}

