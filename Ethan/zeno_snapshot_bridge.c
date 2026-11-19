/* Original format adapter using the Zeno snapshot API.
 * Zeno: Joshua E. Barnes, https://github.com/joshuabarnes/zeno
 * Compile with -DDOUBLEPREC and link libNBody/libClib.
 */
#include "stdinc.h"
#include "getparam.h"
#include "filestruct.h"
#include "vectmath.h"
#include "phatbody.h"
#include "strset.h"
#include <stdint.h>

string defv[] = {";NPZ-runner snapshot adapter", "mode=import", "in=???", "out=???", "VERSION=1.0", NULL};
int main(int argc, string argv[]) {
    string tags[] = {MassTag, PosTag, VelTag, NULL}, got[MaxBodyFields];
    bodyptr bodies = NULL, p;
    int n, i, k;
    real t = 0;
    double row[7];
    stream input, output;
    initparam(argv, defv);
    layout_body(tags, Precision, NDIM);
    input = stropen(getparam("in"), "r");
    output = stropen(getparam("out"), "w");
    if (streq(getparam("mode"), "import")) {
        if (fscanf(input, "%d", &n) != 1 || n < 1) error("Invalid particle count\n");
        bodies = allocate(n * SizeofBody);
        for (i=0; i<n; i++) {
            p = NthBody(bodies, i);
            for (k=0; k<7; k++) if (fscanf(input, "%lf", row+k) != 1) error("Incomplete input\n");
            Mass(p)=row[0];
            for(k=0;k<3;k++) { Pos(p)[k]=row[k+1]; Vel(p)[k]=row[k+4]; }
        }
        put_history(output);
        put_snap(output, &bodies, &n, &t, tags);
    } else if (streq(getparam("mode"), "export")) {
        get_history(input);
        while (get_snap(input, &bodies, &n, &t, got, FALSE, NULL)) {
            int32_t count = n; double time = t;
            if(!set_member(got,MassTag) || !set_member(got,PosTag) || !set_member(got,VelTag)) error("Missing frame fields\n");
            if(fwrite(&count,4,1,output)!=1 || fwrite(&time,8,1,output)!=1) error("Write failed\n");
            for(i=0;i<n;i++) {
                p=NthBody(bodies,i); row[0]=Mass(p);
                for(k=0;k<3;k++) {row[k+1]=Pos(p)[k];row[k+4]=Vel(p)[k];}
                if(fwrite(row,8,7,output)!=7) error("Write failed\n");
            }
            skip_history(input);
        }
    } else error("Unknown mode\n");
    fclose(input); fclose(output); return 0;
}
