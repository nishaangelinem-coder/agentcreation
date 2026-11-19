/* loadgen.c - reproducible interference load: CPU-bound, memory-churn and file-I/O threads
 * in the default time-sharing class.  Usage: loadgen <cpu_threads> <mem_threads> <io_threads> <seconds> */
#define _GNU_SOURCE
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static atomic_int run = 1;
static void *cpu_hog(void *a) { (void)a; volatile double x = 1.0; while (atomic_load(&run)) for (int i = 0; i < 100000; i++) x = sqrt(x * 1.0000001 + 1.0); return NULL; }
static void *mem_hog(void *a) {
    (void)a; const size_t sz = 48u << 20;
    while (atomic_load(&run)) { char *p = malloc(sz); if (!p) continue; for (size_t i = 0; i < sz; i += 4096) p[i] = (char)i; free(p); }
    return NULL;
}
static void *io_hog(void *a) {
    char fn[64]; snprintf(fn, sizeof fn, "/tmp/loadgen_%d_%ld.tmp", (int)getpid(), (long)(size_t)a);
    static char blk[1 << 20]; memset(blk, 0xA5, sizeof blk);
    while (atomic_load(&run)) {
        int fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0600); if (fd < 0) break;
        for (int i = 0; i < 16 && atomic_load(&run); i++) if (write(fd, blk, sizeof blk) < 0) break;
        fsync(fd); close(fd);
    }
    unlink(fn); return NULL;
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: loadgen cpu mem io seconds\n"); return 2; }
    int nc = atoi(argv[1]), nm = atoi(argv[2]), ni = atoi(argv[3]); double sec = atof(argv[4]);
    pthread_t th[256]; int n = 0;
    for (int i = 0; i < nc && n < 256; i++) pthread_create(&th[n++], NULL, cpu_hog, NULL);
    for (int i = 0; i < nm && n < 256; i++) pthread_create(&th[n++], NULL, mem_hog, NULL);
    for (long i = 0; i < ni && n < 256; i++) pthread_create(&th[n++], NULL, io_hog, (void *)i);
    usleep((useconds_t)(sec * 1e6)); atomic_store(&run, 0);
    for (int i = 0; i < n; i++) pthread_join(th[i], NULL);
    return 0;
}
