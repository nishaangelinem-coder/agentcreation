#ifndef RT_UTIL_H
#define RT_UTIL_H
#include <pthread.h>
#include <stdint.h>
#include <time.h>

static inline int64_t now_ns(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}
static inline int64_t thread_cpu_ns(void) {
    struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}
void sleep_until_ns(int64_t t_abs_ns);
/* fifo!=0: SCHED_FIFO at 'prio'; fifo==0: default time-sharing policy. cpu<0: no pinning. */
int  rt_thread_create(pthread_t *th, int fifo, int prio, int cpu, void *(*fn)(void *), void *arg);
int  rt_mutex_init(pthread_mutex_t *mu, int prio_inherit);
int  rt_cond_init(pthread_cond_t *cv);       /* condvar on CLOCK_MONOTONIC */
int  rt_lock_memory(void);                   /* mlockall + stack prefault  */
#endif
