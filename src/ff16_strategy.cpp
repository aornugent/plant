// FF16_Strategy is now a class template (FF16_Strategy_<S>) with all method
// bodies defined inline in inst/include/plant/models/ff16_strategy.h, so the
// active (AD) instantiation can see them. Nothing remains to compile out-of-line
// here; this TU is kept so the build's file list is unchanged.
#include <plant/models/ff16_strategy.h>
