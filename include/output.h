#ifndef OUTPUT_H
#define OUTPUT_H

string create_label();
void write_check_point(ofstream&,double,Grid*);
double load_grid(Grid*,string fname);

#endif
