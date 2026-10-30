#include <limits.h>
#include "defs.h"
#include "blend.h"

void blend_pixel(pixel* s, pixel* d, pixel* b) {
    float a = ( (float)(s->alpha) ) / USHRT_MAX;
    
    d->red   = (a * s->red  ) + ( (1 - a) * b->red  );
    d->green = (a * s->green) + ( (1 - a) * b->green);
    d->blue  = (a * s->blue) +  ( (1 - a) * b->blue );
    d->alpha = USHRT_MAX; // opaque
}
