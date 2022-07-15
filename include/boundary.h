#ifndef BOUNDARY_H
#define BOUNDARY_H

void boundx(Grid*);
void boundy(Grid*, double t=0.0);
void boundz(Grid*);

void boundx_dust(Grid*);
void boundy_dust(Grid*, double t=0.0);
void boundz_dust(Grid*);

void boundx2(Grid*);
void boundy2(Grid*, double t=0.0);
void boundz2(Grid*);
#endif
