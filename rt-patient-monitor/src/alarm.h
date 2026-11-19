/* alarm.h - alarm-system state machine (pure logic, no OS calls; unit-testable off-target).
 * Priority scheme and latching behaviour are modelled on IEC 60601-1-8 concepts. */
#ifndef ALARM_H
#define ALARM_H
#include "common.h"
typedef struct {
    uint8_t cond, latched, active;
    uint64_t k_set; int64_t t_rel_set, t_fire, t_clear;
    uint32_t n_fired;
} alarm_slot_t;
typedef struct { alarm_slot_t s[AL_N]; alarm_prio_t out; } alarm_mgr_t;

extern const char *const alarm_name[AL_N];
extern const alarm_prio_t alarm_prio[AL_N];
extern const uint8_t alarm_latching[AL_N];

void alarm_mgr_init(alarm_mgr_t *m);
/* Returns 1 when the event starts a new annunciation (alarm "fires"). */
int  alarm_mgr_event(alarm_mgr_t *m, alarm_id_t id, int set, uint64_t k, int64_t t_rel, int64_t now);
/* Operator acknowledge: releases a latched alarm only if its condition has ended. Returns 1 if released. */
int  alarm_mgr_ack(alarm_mgr_t *m, alarm_id_t id);
alarm_prio_t alarm_mgr_output(const alarm_mgr_t *m);   /* highest priority among active alarms */
#endif
