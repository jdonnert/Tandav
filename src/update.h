#ifndef UPDATE_H
#define UPDATE_H

#include "includes.h"
#include "cosmology.h"
#include "domain.h"
#include "accel.h"
#include "IO/io.h"
#include "properties.h"
#include "periodic.h"
#include "vector.h"
#include "Gravity_FMM/fmm.h"
#include "Gravity_Tree/tree.h" // <-- add your module .h below



enum Update_Parameters {
	BEFORE_PRESTEP,
	BEFORE_MAIN_LOOP,
	BEFORE_STEP,
	BEFORE_FIRST_KICK,
	BEFORE_SNAPSHOT,
	BEFORE_DRIFT,
	AFTER_DRIFT,
	BEFORE_DOMAIN_UPDATE,
	BEFORE_FORCES,
	BEFORE_SECOND_KICK,
	AFTER_STEP,
	RESTART_CONTINUE
};

void Update(enum Update_Parameters stage);

#endif // UPDATE_H
