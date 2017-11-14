/* Relative orientation program.
   Written October 4, 1998 (modified November 6, 2000) by Bon A. Dewitt
   Modified November 14, 2017 by Arnaud Palha
   This source code is placed in the public domain and may be freely
   copied and modified.
   No warranties are expressed or implied. User accepts all risks
   associated with its use.
   See comments at the end of this source code file for brief
   description and instructions for use.
*/
#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include	<math.h>




/* Maximum number of object points. Change as needed */
#define		MAXPTS		100

//#define		M_PI		3.141592653589793238

/* Macro to access upper triangular matrix stored as 1D array */
#define		INDUT(i,j)	( (((j) * ((j)-1)) / 2) + i )




typedef struct {
	double	m11, m12, m13,
			m21, m22, m23,
			m31, m32, m33,
			so, sp, sk,
			co, cp, ck,
			omega, phi, kappa,
			xl, yl, zl, f;
} PhoParamType;


typedef struct {
	double	xleft, yleft, xleft_res, yleft_res,
			xright, yright, xright_res, yright_res,
			X, Y, Z;
	char	name[24];
} PointType;




void ReadData( char *rootname, PhoParamType *Pphoto, PointType *points, int *Pnumpts );
void AllocateNormals( double **Pnorm, double **Prhs, int numpts );
void solve( double *a, double *b, int n, int invflag );
void RotationMatrix(PhoParamType *Pphoto);
double FormNormals( double *norm, double *rhs, PhoParamType leftphoto,
		PhoParamType rightphoto, PointType *points, int numpts );
void AuxVars( PointType point, PhoParamType photo, double *Pdx, double *Pdy,
			 double *Pdz, double *Pr, double *Ps, double *Pq );
void PartXLYLZL( double part[3][7], PhoParamType photo, double r, double s, double q );
double EpsRes( double *epsilon, double x, double y, double *Pxres, double *Pyres,
			  double f, double r, double s, double q);
void AddCorrections( PhoParamType *Prightphoto, PointType *points, double *rhs,
					int numpts, int iter, double s0, FILE *itfile );
void OutputResults( char *rootname, PhoParamType leftphoto, PhoParamType rightphoto,
				   PointType *points, double *norm, double s0, int numpts );
void pause(void);




int main(void)
{
	char			rootname[50], filename[50];
	FILE			*itfile;
	PhoParamType	leftphoto, rightphoto;
	PointType		points[MAXPTS];
	int				numpts, iter=0, converge=0, diverge=0;
	double			*norm, *rhs, s0=1.0e35, s0old;

	ReadData( rootname, &rightphoto, points, &numpts );

	AllocateNormals( &norm, &rhs, numpts );

	/* set parameters for left photo */
	leftphoto.omega = 0;
	leftphoto.phi = 0;
	leftphoto.kappa = 0;
	leftphoto.xl = 0;
	leftphoto.yl = 0;
	leftphoto.zl = rightphoto.zl;
	leftphoto.f = rightphoto.f;
	RotationMatrix(&leftphoto);

	/* Open file for iteration results */
	strcpy(filename, rootname);
	strcat(filename, ".itr");
	itfile = fopen(filename, "w");

	do {
		iter++;			/* Increment iteration count */
		s0old = s0;		/* Save former value of s0 for convergence test */
		RotationMatrix( &rightphoto); /* Compute rotation matrix elements */
		s0 = FormNormals(norm, rhs, leftphoto, rightphoto, points, numpts);
		printf("ITERATION %d     S0 = %6.5lf\n", iter, s0);

		/* Check for convergence or divergence */
		if (fabs(s0old-s0)/s0 < 0.0001) converge=1;
		else if (s0 > s0old) diverge=1;

		/* Solve for corrections
		   If converged or diverged call for inverse also */
		solve(norm, rhs, 5+numpts*3, converge|diverge);
		AddCorrections( &rightphoto, points, rhs, numpts, iter, s0, itfile );
	} while (!converge && !diverge);

	OutputResults( rootname, leftphoto, rightphoto, points, norm, s0, numpts );
	fclose(itfile);
	pause();
  return 0;
}





void ReadData( char *rootname, PhoParamType *Pphoto, PointType *points, int *Pnumpts )
{
	char	filename[50], tempstr[50];
	FILE	*infile;
	int		NumRead, i;
	double	xl, yl, xr, yr, ParallaxSum;

	printf("Enter root name of relative orientation data file (.dat assumed) --> ");
	scanf("%s", rootname);
	strcpy(filename, rootname);
	strcat(filename, ".dat");
	printf("\n");

	infile = fopen(filename, "r");
	if (infile == NULL) {
		printf("ERROR. Could not open file %s\n", filename);
		pause();
		exit(1);
	}

	/* read focal length */
	fscanf(infile, "%lf", &Pphoto->f);

	/* read points until end of file is reached */
	*Pnumpts = 0;
	do {
		NumRead = fscanf(infile, "%s %lf %lf %lf %lf", tempstr, &xl, &yl, &xr, &yr);
		if (NumRead == 5) {
			/* check for array overflow */
			if (*Pnumpts == MAXPTS) {
				printf("ERROR. More than %d points in data file\n", MAXPTS);
				pause();
				exit(1);
			}
			/* store data for this point and increase count */
			strcpy( points[*Pnumpts].name, tempstr );
			points[*Pnumpts].xleft = xl;
			points[*Pnumpts].yleft = yl;
			points[*Pnumpts].xright = xr;
			points[*Pnumpts].yright = yr;
			/* set initial approximations for object coordinates */
			points[*Pnumpts].X = xl;
			points[*Pnumpts].Y = yl;
			points[*Pnumpts].Z = 0;
			(*Pnumpts)++;

		}
	} while (NumRead == 5);

	fclose(infile);

	/* Quit if not enough pass points */
	if (*Pnumpts < 5) {
		printf("ERROR. Fewer than 5 pass points in data file\n");
		pause();
		exit(1);
	}

	/* set approximations for right photo exterior orientation parameters */
	Pphoto->omega = 0;
	Pphoto->phi = 0;
	Pphoto->kappa = 0;
	Pphoto->yl = 0;
	Pphoto->zl = Pphoto->f;

	/* Set constraint for XL-right based on average parallax of object points.
	   (assumes flight-line axis coordinates) */
	ParallaxSum=0;
	for (i=0; i<*Pnumpts; i++)
		ParallaxSum += points[i].xleft - points[i].xright;
	Pphoto->xl = ParallaxSum / *Pnumpts;

}





void AllocateNormals( double **Pnorm, double **Prhs, int numpts )
{
	int		numunk = 3*numpts + 5;

	*Pnorm = calloc( numunk*(numunk+1)/2+1, sizeof(double) );
	*Prhs = calloc( numunk+1, sizeof(double) );
	if ( (*Pnorm == NULL) || (*Prhs == NULL) ) {
		printf("Dynamic memory allocation failure.\n");
		pause();
		exit(1);
	}
}





void solve( double *a, double *b, int n, int invflag )
/* Solution and inverse by recursive partitioning for upper triangular
   normal equation matrix. When complete, array 'b' is overwritten by
   solution. If 'invflag' is true (i.e. non-zero) then array 'a' is
   overwritten by the inverse also.
*/
{
	int		piv, j, i;
	double	r, *s;

	/* Allocate scratch array for inverse if needed */
	if (invflag) {
		s = calloc(n+1, sizeof(double));
		if (s == NULL) {
			printf("ERROR. Insufficient memory.\n");
			pause();
			exit(1);
		}
	}

	/* Forward elimination */
	for (piv=1; piv<=n; piv++) {
		for (i=piv+1; i<=n; i++) {
			for (j=i; j<=n; j++)
				a[INDUT(i,j)] -= a[INDUT(piv,i)] * a[INDUT(piv,j)]/a[INDUT(piv,piv)];
			b[i] -= a[INDUT(piv,i)] * b[piv]/a[INDUT(piv,piv)];
		}
	}
	/* Back substitution and inverse (if requested) */
	for (piv=n; piv>0; piv--) {
		for (j=piv+1; j<=n; j++) b[piv] -= a[INDUT(piv,j)]*b[j];
		b[piv] /= a[INDUT(piv,piv)];
		if (invflag) {
			for (j=piv+1; j<=n; j++) {
				s[j] = 0.0;
				for (i=piv+1; i<=j; i++) s[j] -= a[INDUT(piv,i)]*a[INDUT(i,j)];
				for (i=j+1; i<=n; i++) s[j] -= a[INDUT(piv,i)]*a[INDUT(j,i)];
				s[j] /= a[INDUT(piv,piv)];
			}
			r = 1.0;
			for (j=piv+1; j<=n; j++) {
				r -= a[INDUT(piv,j)]*s[j];
				a[INDUT(piv,j)] = s[j];
			}
			a[INDUT(piv,piv)] = r / a[INDUT(piv,piv)];
		}
	}
	if (invflag) free(s);
}





void RotationMatrix(PhoParamType *Pphoto)
{
	/* Compute trig functions */
	Pphoto->so = sin(Pphoto->omega);
	Pphoto->co = cos(Pphoto->omega);
	Pphoto->sp = sin(Pphoto->phi);
	Pphoto->cp = cos(Pphoto->phi);
	Pphoto->sk = sin(Pphoto->kappa);
	Pphoto->ck = cos(Pphoto->kappa);

	/* Compute rotation matrix elements */
	Pphoto->m11 = Pphoto->cp * Pphoto->ck;
	Pphoto->m12 = Pphoto->so * Pphoto->sp * Pphoto->ck + Pphoto->co * Pphoto->sk;
	Pphoto->m13 = -Pphoto->co * Pphoto->sp * Pphoto->ck + Pphoto->so * Pphoto->sk;
	Pphoto->m21 = -Pphoto->cp * Pphoto->sk;
	Pphoto->m22 = -Pphoto->so * Pphoto->sp * Pphoto->sk + Pphoto->co * Pphoto->ck;
	Pphoto->m23 = Pphoto->co * Pphoto->sp * Pphoto->sk + Pphoto->so * Pphoto->ck;
	Pphoto->m31 = Pphoto->sp;
	Pphoto->m32 = -Pphoto->so * Pphoto->cp;
	Pphoto->m33 = Pphoto->co * Pphoto->cp;
}





double FormNormals( double *norm, double *rhs, PhoParamType leftphoto,
		PhoParamType rightphoto, PointType *points, int numpts )
{
	int		i, j, pt;
	double	lpart[3][7], lepsilon[3], rpart[3][7], repsilon[3],
			bdot[3][6], bdotdot[3][4], r, s, q, dx, dy, dz, sumres2=0;

	/* Zero out normal equations */
	for (i=1; i<=5+numpts*3; i++) {
		rhs[i] = 0;
		for (j=i; j<=5+numpts*3; j++) norm[INDUT(i,j)] = 0;
	}


	for (pt=0; pt<numpts; pt++) {
		/* Compute auxiliary variables for left photo */
		AuxVars( points[pt], leftphoto, &dx, &dy, &dz, &r, &s, &q);

		/* Compute partials with respect to unknowns for left photo */
		PartXLYLZL( lpart, leftphoto, r, s, q );

		/* Form epsilon (measured - computed), residuals, accumulate residuals-squared */
		sumres2 += EpsRes( lepsilon, points[pt].xleft, points[pt].yleft,
				&(points[pt].xleft_res), &(points[pt].yleft_res), leftphoto.f, r, s, q );

		/* Compute auxiliary variables for right photo */
		AuxVars( points[pt], rightphoto, &dx, &dy, &dz, &r, &s, &q);

		/* Compute partials with respect to unknowns for right photo */
		PartXLYLZL( rpart, rightphoto, r, s, q );
		rpart[1][1] = (rightphoto.f/(q*q)) *
			(r*(-rightphoto.m33*dy + rightphoto.m32*dz) -
			 q*(-rightphoto.m13*dy + rightphoto.m12*dz));
		rpart[1][2] = (rightphoto.f/(q*q)) *
			(r*(rightphoto.cp*dx + rightphoto.so*rightphoto.sp*dy -
			    rightphoto.co*rightphoto.sp*dz) -
			 q*(-rightphoto.sp*rightphoto.ck*dx +
			    rightphoto.so*rightphoto.cp*rightphoto.ck*dy -
				rightphoto.co*rightphoto.cp*rightphoto.ck*dz));
		rpart[1][3] = -rightphoto.f * s / q;
		rpart[2][1] = (rightphoto.f/(q*q)) *
			(s*(-rightphoto.m33*dy + rightphoto.m32*dz) -
			 q*(-rightphoto.m23*dy + rightphoto.m22*dz));
		rpart[2][2] = (rightphoto.f/(q*q)) *
			(s*(rightphoto.cp*dx + rightphoto.so*rightphoto.sp*dy -
			    rightphoto.co*rightphoto.sp*dz) -
			 q*(rightphoto.sp*rightphoto.sk*dx -
			    rightphoto.so*rightphoto.cp*rightphoto.sk*dy +
				rightphoto.co*rightphoto.cp*rightphoto.sk*dz));
		rpart[2][3] = rightphoto.f * r/ q;

		/* Form epsilon (measured - computed), residuals, accumulate residuals-squared */
		sumres2 += EpsRes( repsilon, points[pt].xright, points[pt].yright,
				&(points[pt].xright_res), &(points[pt].yright_res),
				rightphoto.f, r, s, q );

		/* Add contribution of this point to normal equations */
		/* Left photo first */
		for (i=1; i<=2; i++)
			for (j=1; j<=3; j++) bdotdot[i][j] = -lpart[i][j+3];
		for (i=1; i<=3; i++) {
			for (j=i; j<=3; j++)
				norm[INDUT(5+pt*3+i, 5+pt*3+j)] +=  bdotdot[1][i]*bdotdot[1][j] +
													bdotdot[2][i]*bdotdot[2][j];
			rhs[5+pt*3+i] += bdotdot[1][i]*lepsilon[1] + bdotdot[2][i]*lepsilon[2];
		}
		/* Right photo second */
		for (i=1; i<=2; i++) {
			for (j=1; j<=3; j++) bdot[i][j] = rpart[i][j];
			for (j=4; j<=5; j++) bdot[i][j] = rpart[i][j+1];
			for (j=1; j<=3; j++) bdotdot[i][j] = -rpart[i][j+3];
		}
		for (i=1; i<=5; i++) {
			for (j=i; j<=5; j++)
				norm[INDUT(i,j)] += bdot[1][i]*bdot[1][j] + bdot[2][i]*bdot[2][j];
			rhs[i] += bdot[1][i]*repsilon[1] + bdot[2][i]*repsilon[2];
			for (j=1; j<=3; j++)
				norm[INDUT(i,5+pt*3+j)] +=  bdot[1][i]*bdotdot[1][j] +
											bdot[2][i]*bdotdot[2][j];
		}
		for (i=1; i<=3; i++) {
			for (j=i; j<=3; j++)
				norm[INDUT(5+pt*3+i, 5+pt*3+j)] +=  bdotdot[1][i]*bdotdot[1][j] +
													bdotdot[2][i]*bdotdot[2][j];
			rhs[5+pt*3+i] += bdotdot[1][i]*repsilon[1] + bdotdot[2][i]*repsilon[2];
		}
	}

	/* Compute and return standard error of unit weight */
	if (numpts>5) return sqrt(sumres2/(numpts-5));
	else return 1.0e30; /* Zero degrees of freedom, s0 = infinity */
}





void AuxVars( PointType point, PhoParamType photo, double *Pdx, double *Pdy,
			 double *Pdz, double *Pr, double *Ps, double *Pq )
{
	*Pdx = point.X - photo.xl;
	*Pdy = point.Y - photo.yl;
	*Pdz = point.Z - photo.zl;
	*Pr = photo.m11*(*Pdx) + photo.m12*(*Pdy) + photo.m13*(*Pdz);
	*Ps = photo.m21*(*Pdx) + photo.m22*(*Pdy) + photo.m23*(*Pdz);
	*Pq = photo.m31*(*Pdx) + photo.m32*(*Pdy) + photo.m33*(*Pdz);
}




void PartXLYLZL( double part[3][7], PhoParamType photo, double r, double s, double q )
{
	part[1][4] = -photo.f * (r*photo.m31 - q*photo.m11) / (q*q);
	part[1][5] = -photo.f * (r*photo.m32 - q*photo.m12) / (q*q);
	part[1][6] = -photo.f * (r*photo.m33 - q*photo.m13) / (q*q);
	part[2][4] = -photo.f * (s*photo.m31 - q*photo.m21) / (q*q);
	part[2][5] = -photo.f * (s*photo.m32 - q*photo.m22) / (q*q);
	part[2][6] = -photo.f * (s*photo.m33 - q*photo.m23) / (q*q);
}




double EpsRes( double *epsilon, double x, double y, double *Pxres, double *Pyres,
			  double f, double r, double s, double q)
{
	/* Form epsilon (measured - computed) */
	epsilon[1] = x + f * r / q;
	epsilon[2] = y + f * s / q;

	/* Save residuals */
	*Pxres = -epsilon[1];
	*Pyres = -epsilon[2];

	/* Return sum of residuals-squared */
	return (epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2]);
}




void AddCorrections( PhoParamType *Prightphoto, PointType *points, double *rhs,
					int numpts, int iter, double s0, FILE *itfile )
{
	int	pt;

	fprintf(itfile, "ITERATION: %d      S0 estimate: %6.5lf\n\n\n", iter, s0);
	Prightphoto->omega += rhs[1];
	Prightphoto->phi += rhs[2];
	Prightphoto->kappa += rhs[3];
	Prightphoto->yl += rhs[4];
	Prightphoto->zl += rhs[5];
	fprintf(itfile, "Right photo exterior orientation parameters:\n\n"
		"%5s  %10s  %10s\n", "Param", "Correction", "New Approx");
	fprintf(itfile, "%5s  %10.5lf  %10.5lf  (degrees)\n", "omega",
		rhs[1] * 180 / M_PI, Prightphoto->omega * 180 / M_PI );
	fprintf(itfile, "%5s  %10.5lf  %10.5lf  (degrees)\n", " phi ",
		rhs[2] * 180 / M_PI, Prightphoto->phi * 180 / M_PI );
	fprintf(itfile, "%5s  %10.5lf  %10.5lf  (degrees)\n", "kappa",
		rhs[3] * 180 / M_PI, Prightphoto->kappa * 180 / M_PI );
	fprintf(itfile, "%5s  %10s  %10.5lf\n", "XL  ", " ", Prightphoto->xl);
	fprintf(itfile, "%5s  %10.5lf  %10.5lf\n", "YL  ", rhs[4], Prightphoto->yl);
	fprintf(itfile, "%5s  %10.5lf  %10.5lf\n", "ZL  ", rhs[5], Prightphoto->zl);

	fprintf(itfile, "\n\nObject space coordinates:\n\n"
		"%8s %10s %10s %10s  %10s %10s %10s\n", "Point", "delta X ", "delta Y ",
		"delta Z ", "X new ", "Y new ", "Z new ");
	for (pt=0; pt<numpts; pt++) {
		points[pt].X += rhs[5+pt*3+1];
		points[pt].Y += rhs[5+pt*3+2];
		points[pt].Z += rhs[5+pt*3+3];
		fprintf(itfile, "%8s %10.5lf %10.5lf %10.5lf  %10.5lf %10.5lf %10.5lf\n",
			points[pt].name, rhs[5+pt*3+1], rhs[5+pt*3+2], rhs[5+pt*3+3],
			points[pt].X, points[pt].Y, points[pt].Z );
	}
	fprintf(itfile, "\n\n\n\n");
}





void OutputResults( char *rootname, PhoParamType leftphoto, PhoParamType rightphoto,
				   PointType *points, double *norm, double s0, int numpts )
{
	FILE	*outfile;
	char	filename[50];
	int		i;
	double	lsumx2=0, lsumy2=0, rsumx2=0, rsumy2=0;

	strcpy(filename, rootname);
	strcat(filename, ".out");
	outfile = fopen(filename, "w");
	printf("\nResults are in file %s\n", filename);

	fprintf(outfile, "Exterior orientation parameters:\n\n");
	fprintf(outfile, "%-10s  %10s  %10s  %8s\n", "Parameter", "Left pho",
		"Right pho", "SD right");
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf  %8.4lf\n", "Omega(deg)",
		leftphoto.omega*180/M_PI, rightphoto.omega*180/M_PI,
		s0*sqrt(norm[INDUT(1,1)])*180/M_PI);
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf  %8.4lf\n", "Phi(deg)",
		leftphoto.phi*180/M_PI, rightphoto.phi*180/M_PI,
		s0*sqrt(norm[INDUT(2,2)])*180/M_PI );
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf  %8.4lf\n", "Kappa(deg)",
		leftphoto.kappa*180/M_PI, rightphoto.kappa*180/M_PI,
		s0*sqrt(norm[INDUT(3,3)])*180/M_PI );
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf\n", "XL", leftphoto.xl, rightphoto.xl);
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf  %8.4lf\n", "YL",
		leftphoto.yl, rightphoto.yl, s0*sqrt(norm[INDUT(4,4)]) );
	fprintf(outfile, "%-10s  %10.4lf  %10.4lf  %8.4lf\n", "ZL",
		leftphoto.zl, rightphoto.zl, s0*sqrt(norm[INDUT(5,5)]) );

	fprintf(outfile, "\n\n\nObject space coordinates:\n\n%8s %9s %9s %9s %7s %7s %7s\n",
		"point", "X   ", "Y   ", "Z   ", "sdX ", "sdY ", "sdZ ");
	for (i=0; i<numpts; i++) {
		fprintf(outfile, "%8s %9.4lf %9.4lf %9.4lf %7.4lf %7.4lf %7.4lf\n",
				points[i].name, points[i].X, points[i].Y, points[i].Z,
				s0*sqrt(norm[INDUT(5+i*3+1,5+i*3+1)]),
				s0*sqrt(norm[INDUT(5+i*3+2,5+i*3+2)]),
				s0*sqrt(norm[INDUT(5+i*3+3,5+i*3+3)]) );
	}

	fprintf(outfile, "\n\n\nPhoto coordinate residuals:\n\n"
		"%8s %7s %7s %7s %7s\n", "point", "xl-res", "yl-res", "xr-res", "yr-res");
	for (i=0; i<numpts; i++) {
		fprintf(outfile, "%8s %7.4lf %7.4lf %7.4lf %7.4lf\n",
				points[i].name, points[i].xleft_res, points[i].yleft_res,
				points[i].xright_res, points[i].yright_res );
		lsumx2 += points[i].xleft_res * points[i].xleft_res;
		lsumy2 += points[i].yleft_res * points[i].yleft_res;
		rsumx2 += points[i].xright_res * points[i].xright_res;
		rsumy2 += points[i].yright_res * points[i].yright_res;
	}
	fprintf(outfile, "\n%8s %7.4lf %7.4lf %7.4lf %7.4lf\n",
		"RMS", sqrt(lsumx2/numpts), sqrt(lsumy2/numpts),
		sqrt(lsumx2/numpts), sqrt(lsumy2/numpts) );

	fprintf(outfile, "\n\n\nStandard error of unit weight: %7.4lf\n"
		"Degrees of freedom: %d\n", s0, numpts-5);
	fclose(outfile);
}





void pause(void)
{
  printf("Press the ENTER key to end program.\n");
  char c = 0;
  scanf("%c", &c);
  do {
    c = getchar();
  } while(c!='\n');
}






/*
This program performs a relative orientation solution for a stereopair of
near vertical photos by least squares. No weights are used for photo
coordinates or object space coordinates. All x and y photo coordinates
must be corrected for principal point offsets, lens distortions, and
atmospheric refraction as needed. Data file must be plain ASCII text and
have a .dat extension. The format of the file is as follows:

  Line 1: <focal length>
  Lines 2 through end of file:
  <point name> <left photo x> <left photo y> <right photo x> <right photo y>

Sample data file:


152.113
a -4.878   1.974  -97.920  -2.923
b 89.307   2.709   -1.507  -1.856
c  0.261  84.144  -90.917  78.970
d 90.334  83.843   -1.571  79.470
e -4.668 -86.821 -100.060 -95.748
f 88.599 -85.274   -0.965 -94.319


Output is sent to a file with the same root name and a .out extension.
Iteration information is sent to a file with a .itr extension.

*/
