#define _GNU_SOURCE
#include "rt_util.h"
#include <errno.h>
#include <sched.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#ifdef __QNX__
#include <sys/neutrino.h>
#endif

void sleep_until_ns(int64_t t) {
    struct timespec ts = { .tv_sec = t / 1000000000LL, .tv_nsec = t % 1000000000LL };
    while (clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &ts, NULL) == EINTR) { }
}

struct tramp { void *(*fn)(void *); void *arg; int cpu; };
static void *trampoline(void *p) {
    struct tramp t = *(struct tramp *)p;
    ((struct tramp *)p)->fn = NULL;          /* signal: copied */
    if (t.cpu >= 0) {
#ifdef __QNX__
        /* QNX Neutrino: restrict this thread to one core (untested port path). */
        ThreadCtl(_NTO_TCTL_RUNMASK, (void *)(uintptr_t)(1u << t.cpu));
#else
        cpu_set_t set; CPU_ZERO(&set); CPU_SET(t.cpu, &set);
        pthread_setaffinity_np(pthread_self(), sizeof set, &set);
#endif
    }
    return t.fn(t.arg);
}

int rt_thread_create(pthread_t *th, int fifo, int prio, int cpu, void *(*fn)(void *), void *arg) {
    static struct tramp slots[16]; static int nslot;
    struct tramp *t = &slots[nslot++ % 16];
    t->fn = fn; t->arg = arg; t->cpu = cpu;
    pthread_attr_t at; pthread_attr_init(&at);
    pthread_attr_setstacksize(&at, 1 << 20);
    if (fifo) {
        struct sched_param sp = { .sched_priority = prio };
        pthread_attr_setinheritsched(&at, PTHREAD_EXPLICIT_SCHED);
        pthread_attr_setschedpolicy(&at, SCHED_FIFO);
        pthread_attr_setschedparam(&at, &sp);
    }
    int rc = pthread_create(th, &at, trampoline, t);
    pthread_attr_destroy(&at);
    if (rc) fprintf(stderr, "pthread_create(fifo=%d prio=%d): %s\n", fifo, prio, strerror(rc));
    return rc;
}

int rt_mutex_init(pthread_mutex_t *mu, int pi) {
    pthread_mutexattr_t ma; pthread_mutexattr_init(&ma);
    if (pi) pthread_mutexattr_setprotocol(&ma, PTHREAD_PRIO_INHERIT);
    int rc = pthread_mutex_init(mu, &ma);
    pthread_mutexattr_destroy(&ma);
    return rc;
}

int rt_cond_init(pthread_cond_t *cv) {
    pthread_condattr_t ca; pthread_condattr_init(&ca);
    pthread_condattr_setclock(&ca, CLOCK_MONOTONIC);
    int rc = pthread_cond_init(cv, &ca);
    pthread_condattr_destroy(&ca);
    return rc;
}

int rt_lock_memory(void) {
    if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) { perror("mlockall"); return -1; }
    volatile char stack[256 * 1024];
    for (size_t i = 0; i < sizeof stack; i += 4096) stack[i] = 0;
    return 0;
}
