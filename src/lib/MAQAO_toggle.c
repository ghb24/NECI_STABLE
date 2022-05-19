#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <stdbool.h>

void toggle_lprof() {
    /* static bool profiling = false; */
    /* if (getenv("LPROF_UG_ON") != NULL) { */
    /*     kill(getppid(), SIGTSTP); */
    /*     printf("MAQAO profiling switched %s\n", profiling ? "OFF" : "ON"); */
    /*     profiling = ! profiling; */
    /* } */
}
