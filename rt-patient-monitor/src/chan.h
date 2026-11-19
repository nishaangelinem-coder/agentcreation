/* chan.h - bounded FIFO channel between threads (non-blocking send, blocking receive).
 * The mutex optionally uses the priority-inheritance protocol. */
#ifndef CHAN_H
#define CHAN_H
#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
typedef struct {
    pthread_mutex_t mu; pthread_cond_t cv;
    size_t esz, cap, head, cnt, high_water;
    uint64_t drops; int closed;
    unsigned char *buf;
} chan_t;
int  chan_init(chan_t *c, size_t esz, size_t cap, int prio_inherit);
int  chan_push(chan_t *c, const void *e);   /* 0 ok, -1 full (element dropped, counted) */
int  chan_pop(chan_t *c, void *e);          /* 0 ok, -1 closed and empty               */
void chan_close(chan_t *c);
#endif
