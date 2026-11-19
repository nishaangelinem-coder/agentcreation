#include "chan.h"
#include "rt_util.h"
#include <stdlib.h>
#include <string.h>

int chan_init(chan_t *c, size_t esz, size_t cap, int pi) {
    memset(c, 0, sizeof *c);
    c->esz = esz; c->cap = cap;
    c->buf = calloc(cap, esz);               /* allocated once, before the RT phase */
    if (!c->buf) return -1;
    rt_mutex_init(&c->mu, pi); rt_cond_init(&c->cv);
    return 0;
}
int chan_push(chan_t *c, const void *e) {
    int rc = 0;
    pthread_mutex_lock(&c->mu);
    if (c->cnt == c->cap) { c->drops++; rc = -1; }
    else {
        memcpy(c->buf + ((c->head + c->cnt) % c->cap) * c->esz, e, c->esz);
        if (++c->cnt > c->high_water) c->high_water = c->cnt;
        pthread_cond_signal(&c->cv);
    }
    pthread_mutex_unlock(&c->mu);
    return rc;
}
int chan_pop(chan_t *c, void *e) {
    pthread_mutex_lock(&c->mu);
    while (c->cnt == 0 && !c->closed) pthread_cond_wait(&c->cv, &c->mu);
    if (c->cnt == 0) { pthread_mutex_unlock(&c->mu); return -1; }
    memcpy(e, c->buf + c->head * c->esz, c->esz);
    c->head = (c->head + 1) % c->cap; c->cnt--;
    pthread_mutex_unlock(&c->mu);
    return 0;
}
void chan_close(chan_t *c) {
    pthread_mutex_lock(&c->mu); c->closed = 1; pthread_cond_broadcast(&c->cv); pthread_mutex_unlock(&c->mu);
}
