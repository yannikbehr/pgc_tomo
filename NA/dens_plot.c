#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int dplot(double *mdls, double *data, double* kosu1, double *kosu2, 
		int ne, int nd, int nlayer, int prmidx, double dmax, double dmin,
		double dreso, double smax, double smin, double sreso, double mf,
		double *smean, double *wtmean, int *idmax){

	int imod,i;
	int id1, id2, iv1, iv2, id, iv;
	double thick, vold, v, dh, dold;
	int rows, columns;
	rows = (int)((dmax-dmin)/dreso);
	columns = (int)((smax-smin)/sreso);
	/* include end points */
	rows++;
	columns++;
	/* iterate over inversion runs */
	for(imod=0;imod<ne;imod++){
		dh = 0.0;
		/* iterate over layers */
		for(i=0;i<nlayer;i++){
			thick = mdls[nd*imod+i];
			v = mdls[nd*imod+prmidx*nlayer+i];
			if(v > smax)v=smax;
			dold = dh;
			dh = dh + thick;
			if(thick < 0.001)dh = dmax;
			if(dh > dmax)dh = dmax;
			/* calculate grid indices; round to the closest grid point*/
			id1 = (int)(((dold-dmin)/dreso) + 0.5);
			id2 = (int)(((dh-dmin)/dreso) + 0.5);
			iv = (int)(((v-smin)/sreso)+0.5);
			/* if id1 and id2 are different increase
	 cell count for each cell between id1 and id2 by 1*/
			for(id=id1;id<=id2;id++){
				if(data[imod] <= mf){
					kosu1[id*columns+iv]++;
					smean[id] = smean[id]+(v/data[imod]);
					wtmean[id] = wtmean[id] + 1./data[imod];
					if(id>*idmax)*idmax = id;
				}
				kosu2[id*columns+iv]++;
			}
			if(i>0){
				iv1 = (int)(((vold-smin)/sreso)+0.5);
				iv2 = (int)(((v-smin)/sreso)+0.5);
				/* if iv1 and iv2 are different increase
	   cell count for each cell between iv1 and iv2 by 1*/
				if(iv1<iv2){
					for(iv=iv1;iv<=iv2;iv++){
						if(data[imod] <= mf){
							kosu1[id1*columns+iv]++;
							//if(id1>*idmax) *idmax = id1;
						}
						kosu2[id1*columns+iv]++;
					}
				}else if(iv1>iv2){
					for(iv=iv1;iv>=iv2;iv--){
						if(data[imod] <= mf){
							kosu1[id1*columns+iv]++;
							//if(id1>*idmax) *idmax = id1;
						}
						kosu2[id1*columns+iv]++;
					}
				}
			}
			vold = v;
		}
	}
	for(i=0;i<*idmax;i++){
		if(smean[i] > 0){
			smean[i] = smean[i]/wtmean[i];
		}
	}
	return 1;
}

