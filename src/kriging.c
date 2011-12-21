#include <R.h>
#include <math.h>


/* Point-in-polygon */
int pointinpoly(int n, double *X, double *Y, double x, double y) {
  int i, j, c;

  for(i=0,j=n-1,c=0; i<n; j=i++) {
    if(((Y[i]>y)!=(Y[j]>y)) && (x<((X[j]-X[i])*(y-Y[i])/(Y[j]-Y[i])+X[i]))) {
      c = !c;
    }
  }
  return c;
}

/* Mathematical matrix functions */
/* n x p, p x m */
void mmult(double *x, double *y, int *n, int *p, int *m, double *z) {
  /* x: n x p
     y: p x m
     0: n x m
  */
  int i, j, k;

  for(i=0;i<*n;i++) {
    for(j=0;j<*m;j++) {
      z[i+j*(*n)] = 0;
      for(k=0;k<*p;k++) {
        z[i+j*(*n)] += x[i+k*(*n)] * y[k+j*(*p)];
      }
    }
  }
}

/* Variogram model functions */
double sphermodel(double h, double nugget, double range, double sill) {
  if(h>=range) return(sill);
  return( (sill-nugget) * (1.5*(h/range) - 0.5*pow(h/range, 3)) + nugget );
}

double expmodel(double h, double nugget, double range, double sill, double a) {
  return( (sill-nugget) * (1 - exp(-h/(range*a))) + nugget );
}

double gaussmodel(double h, double nugget, double range, double sill, double a) {
  if(h>range) return(nugget);
  return( (sill-nugget) * (1-exp(-(pow(h, 2))/(pow(range, 2)*a))) * nugget );
}

/* n-dimensional euclidian distance matrices */
void onedimdist(double *x, int *n, double *v) {
  int i, j;
  for(i=0;i<*n;i++) {
    for(j=0;j<*n;j++) {
      v[j+(*n)*i] = pow(x[i]-x[j], 2);
    }
  } 
}

void twodimdist(double *x, double *y, int *n, double *d) {
  int i, j;
  for(i=0;i<*n;i++) {
    for(j=0;j<*n;j++) {
      d[j+(*n)*i] = pow(pow(x[i]-x[j], 2) + pow(y[i]-y[j], 2), 0.5);
    }
  }
}

/* Model fitting/gridding/prediction functions */

void krigfit(double *d, double *nugget, double *range, double *sill, int *model, int *n, double *a) {
  /* Models:
       0: Spherical (default)
       1: Exponential
       2: Gaussian
  */
  int i, j;

  a[(*n+1)*(*n+1) - 1] = 0;
  if(*model==0) {
    for(i=0;i<(*n);i++) {
      for(j=0;j<(*n);j++) {
        a[j+(*n+1)*i] = sphermodel(d[j+(*n)*i], *nugget, *range, *sill);
      }
    }
  }

  if(*model==1) {
    for(i=0;i<(*n);i++) {
      for(j=0;j<(*n);j++) {
        a[j+(*n+1)*i] = expmodel(d[j+(*n)*i], *nugget, *range, *sill, 0.33333);
      }
    }
  }

  if(*model==1) {
    for(i=0;i<(*n);i++) {
      for(j=0;j<(*n);j++) {
        a[j+(*n+1)*i] = gaussmodel(d[j+(*n)*i], *nugget, *range, *sill, 0.33333);
      }
    }
  }


}

void kriggrid(double *blx, double *bly, double *pixel, int *xpixels, int *ypixels, double *gx, double *gy) {
  int i, j;
  int k = 0;

  for(i=0;i<*xpixels;i++) {
    for(j=0;j<*ypixels;j++) {
      gx[k] = *blx + i*(*pixel);
      gy[k] = *bly + j*(*pixel);
      k++;
    }
  }
}

void krigpolygons(int *n, double *X, double *Y, int *nlist, int *npoly, double *polygonsx, double *polygonsy, int *Gn, double *Gx, double *Gy) {
  int i, j, k;
  for(i=0, k=0; i<*nlist;i++) {
    for(j=0; j<*n;j++) {
      if(pointinpoly(npoly[i+1]-npoly[i], polygonsx+npoly[i], polygonsy+npoly[i], X[j], Y[j])) {
        Gx[k] = X[j];
        Gy[k] = Y[j];
        k++;
        (*Gn)++;
      }
    }
  }
}

void krigpred(double *x, double *y, double *response, double *gx, double *gy, double *inva, double *nugget, double *range, double *sill, int *model, int *n, int *ng, int *nplus, int *one, double *gpred) {
  int i, j, k;
  double xdist;
  double R[*nplus];
  double pred[1];
  double invaXR[*nplus];

  if(*model==0) {
    for(i=0;i<(*ng);i++) {
      for(j=0;j<(*n);j++) {
        xdist = pow(pow(x[j]-gx[i], 2) + pow(y[j]-gy[i], 2), 0.5);
        R[j] = sphermodel(xdist, *nugget, *range, *sill);
      }
      R[*n] = 1;

      /* Y' * (A^-1 * R) */
      mmult(inva, R, nplus, nplus, one, invaXR);
      mmult(response, invaXR, one, n, one, pred);

      gpred[i] = pred[0];
    }
  }

  if(*model==1) {
    for(i=0;i<(*ng);i++) {
      for(j=0;j<(*n);j++) {
        xdist = pow(pow(x[j]-gx[i], 2) + pow(y[j]-gy[i], 2), 0.5);
        R[j] = expmodel(xdist, *nugget, *range, *sill, 0.33333);
      }
      R[*n] = 1;

      /* Y' * (A^-1 * R) */
      mmult(inva, R, nplus, nplus, one, invaXR);
      mmult(response, invaXR, one, n, one, pred);

      gpred[i] = pred[0];
    }
  }

  if(*model==2) {
    for(i=0;i<(*ng);i++) {
      for(j=0;j<(*n);j++) {
        xdist = pow(pow(x[j]-gx[i], 2) + pow(y[j]-gy[i], 2), 0.5);
        R[j] = gaussmodel(xdist, *nugget, *range, *sill, 0.33333);
      }
      R[*n] = 1;

      /* Y' * (A^-1 * R) */
      mmult(inva, R, nplus, nplus, one, invaXR);
      mmult(response, invaXR, one, n, one, pred);

      gpred[i] = pred[0];
    }
  }
}

/* Image method loop optimized */
void krigimage(int *x, int *y, double *z, int *n, int *a, int *b, int *na, int *nb, double *c, int *d) {
  int i, j, k;

  for(i=0;i<*na;i++) {
    for(j=0;j<*nb;j++) {
      for(k=0;k<*n;k++) {
        if(x[k]==a[i]&&y[k]==b[j]) {
          c[j+i*(*nb)] = z[k];
          d[j+i*(*nb)] = 1;
	  break;
        }
      }
    }
  }
}
