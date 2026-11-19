#include "alarm.h"
#include <string.h>

const char *const alarm_name[AL_N] = { "ASYSTOLE", "TACHY", "BRADY", "DESAT", "SPO2_LOW",
                                       "LEADS_OFF", "LOW_PERF", "SYS_STALL", "DATA_LOSS" };
const alarm_prio_t alarm_prio[AL_N] = { PR_HIGH, PR_HIGH, PR_HIGH, PR_HIGH, PR_MED,
                                        PR_MED, PR_LOW, PR_HIGH, PR_HIGH };
/* High-priority physiological alarms latch until acknowledged; technical alarms do not. */
const uint8_t alarm_latching[AL_N] = { 1, 1, 1, 1, 0, 0, 0, 0, 0 };

void alarm_mgr_init(alarm_mgr_t *m) { memset(m, 0, sizeof *m); }

alarm_prio_t alarm_mgr_output(const alarm_mgr_t *m) {
    alarm_prio_t p = PR_NONE;
    for (int i = 0; i < AL_N; i++) if (m->s[i].active && alarm_prio[i] > p) p = alarm_prio[i];
    return p;
}

int alarm_mgr_event(alarm_mgr_t *m, alarm_id_t id, int set, uint64_t k, int64_t t_rel, int64_t now) {
    alarm_slot_t *a = &m->s[id]; int fired = 0;
    if (set) {
        if (!a->cond) {
            a->cond = 1;
            if (!a->active) { a->active = 1; a->k_set = k; a->t_rel_set = t_rel; a->t_fire = now; a->n_fired++; fired = 1; }
            if (alarm_latching[id]) a->latched = 1;
        }
    } else if (a->cond) {
        a->cond = 0; a->t_clear = now;
        if (!a->latched) a->active = 0;
    }
    m->out = alarm_mgr_output(m);
    return fired;
}

int alarm_mgr_ack(alarm_mgr_t *m, alarm_id_t id) {
    alarm_slot_t *a = &m->s[id];
    if (a->latched && !a->cond) { a->latched = 0; a->active = 0; m->out = alarm_mgr_output(m); return 1; }
    return 0;
}
