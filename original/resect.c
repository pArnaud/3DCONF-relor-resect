/* Single photo space resection program.
   Written September 27, 1998 (Modified November 6, 2000) by Bon A. Dewitt
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
#include	<malloc.h>
#include	<math.h>
#include	<conio.h>



/* Maximum number of image points. Change as needed */
#define		MAXPTS		100

#define		M_PI		3.141592653589793238

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
	double	x, y, xres, yres,
			X, Y, Z;
	char	name[24];
} PointType;




void ReadData( char *rootname, PhoParamType *Pphoto, PointType *points, int *Pnumpts );
void ComputeApproximations( PhoParamType *Pphoto, PointType *points, int numpts );
void solve( double *a, double *b, int n, int invflag );
void RotationMatrix(PhoParamType *Pphoto);
double FormNormals( double *norm, double *rhs, PhoParamType photo, PointType *points,
				   int numpts );
void AddCorrections( PhoParamType *Pphoto, double *rhs );
void OutputResults( char *rootname, PhoParamType photo, PointType *points,
				   double *norm, double s0, int numpts );
void pause(void);



void main(void)
{
	char			rootname[50];
	PhoParamType	photo;
	PointType		points[MAXPTS];
	int				numpts, iter=0, converge=0, diverge=0;
	double			norm[22], rhs[7], s0=1.0e30, s0old;

	ReadData( rootname, &photo, points, &numpts );
	ComputeApproximations( &photo, points, numpts );

	do {
		iter++;			/* Increment iteration count */
		s0old = s0;		/* Save former value of s0 for convergence test */
		RotationMatrix(&photo);	/* Compute rotation matrix elements */
		s0 = FormNormals(norm, rhs, photo, points, numpts);
		printf("ITERATION %d     S0 = %6.5lf\n", iter, s0);

		/* Check for convergence or divergence */
		if (fabs(s0old-s0)/s0 < 0.0001) converge=1;
		else if (s0 > s0old) diverge=1;

		/* Solve for corrections
		   If converged or diverged call for inverse also */
		solve(norm, rhs, 6, converge|diverge);
		AddCorrections( &photo, rhs );
	} while (!converge && !diverge);

	OutputResults( rootname, photo, points, norm, s0, numpts );
	pause();
}





void ReadData( char *rootname, PhoParamType *Pphoto, PointType *points, int *Pnumpts )
{
	char	filename[50], tempstr[50];
	FILE	*infile;
	int		NumRead;
	double	x, y, X, Y, Z;

	printf("Enter root name of resection data file (.dat assumed) --> ");
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
		NumRead = fscanf(infile, "%s %lf %lf %lf %lf %lf", tempstr, &x, &y, &X, &Y, &Z);
		if (NumRead == 6) {
			/* check for array overflow */
			if (*Pnumpts == MAXPTS) {
				printf("ERROR. More than %d points in data file\n", MAXPTS);
				pause();
				exit(1);
			}
			/* store data for this point and increase count */
			strcpy( points[*Pnumpts].name, tempstr );
			points[*Pnumpts].x = x;
			points[*Pnumpts].y = y;
			points[*Pnumpts].X = X;
			points[*Pnumpts].Y = Y;
			points[*Pnumpts].Z = Z;
			(*Pnumpts)++;

		}
	} while (NumRead == 6);

	fclose(infile);

	/* Quit if not enough control points */
	if (*Pnumpts < 3) {
		printf("ERROR. Fewer than 3 control points in data file\n");
		pause();
		exit(1);
	}
}





void ComputeApproximations( PhoParamType *Pphoto, PointType *points, int numpts )
{
	int		i, j, count=0;
	double	H, sumH=0, dx, dy, cx, cy, a, b, c, dij2, sqrtterm,
			norm[11], rhs[5];

	/* Assume vertical photo */
	Pphoto->omega = 0;
	Pphoto->phi = 0;

	/* Compute height by method of Section 6-9 of text
	   Use each possible pair of points and take average */
	for (i=0; i<numpts-1; i++) {
		for (j=i+1; j<numpts; j++) {
			dx = points[j].x - points[i].x;
			dy = points[j].y - points[i].y;
			cx = points[i].x * points[i].Z  -  points[j].x * points[j].Z;
			cy = points[i].y * points[i].Z  -  points[j].y * points[j].Z;
			dij2 = (points[j].X-points[i].X)*(points[j].X-points[i].X) +
				   (points[j].Y-points[i].Y)*(points[j].Y-points[i].Y);
			/* Quadratic formula coefficients */
			a = dx*dx + dy*dy;
			b = 2.0 * (dx*cx + dy*cy);
			c = cx*cx + cy*cy - Pphoto->f * Pphoto->f * dij2;
			/* Solve quadratic */
			sqrtterm = b*b-4*a*c;
			if (sqrtterm < 0.0) {
				printf("ERROR. Flying height calculation failed due to square root\n"
						"of a negative number.\n");
				pause();
				exit(1);
			}
			H = (-b + sqrt(sqrtterm))/(2*a);
			/* Accumulate subtotal and increment count for average calculation */
			sumH += H;
			count++;
		}
	}
	/* ZL = average H */
	Pphoto->zl = sumH / count;

	/* Compute 2D conformal transformation using ground coordinates from
	   a vertical photo and actual ground coordinates.
	   Transformation coefficients Tx, Ty, and theta will give
	   xl, yl, and kappa, respectively */

	/* First, zero out normal equations */
	for (i=1; i<=4; i++) {
		rhs[i] = 0;
		for (j=i; j<=4; j++) norm[INDUT(i,j)] = 0;
	}

	/* Form normal equations directly, one point at a time */
	for (i=0; i<numpts; i++) {
		double	xvert, yvert;
		xvert = points[i].x * (Pphoto->zl - points[i].Z) / Pphoto->f;
		yvert = points[i].y * (Pphoto->zl - points[i].Z) / Pphoto->f;
		norm[INDUT(1,1)] += xvert*xvert + yvert*yvert;
		norm[INDUT(1,3)] += xvert;
		norm[INDUT(1,4)] += yvert;
		rhs[1] += xvert*points[i].X + yvert*points[i].Y;
		rhs[2] += xvert*points[i].Y - yvert*points[i].X;
		rhs[3] += points[i].X;
		rhs[4] += points[i].Y;
	}
	norm[INDUT(2,2)] = norm[INDUT(1,1)];
	norm[INDUT(2,3)] = -norm[INDUT(1,4)];
	norm[INDUT(2,4)] = norm[INDUT(1,3)];
	norm[INDUT(3,3)] = numpts;
	norm[INDUT(4,4)] = numpts;

	/* Normal equations are formed, now solve */
	solve( norm, rhs, 4, 0 );

	/* Obtain transformation parameters for initial approximations */
	Pphoto->kappa = atan2( rhs[2], rhs[1] );
	Pphoto->xl = rhs[3];
	Pphoto->yl = rhs[4];
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





double FormNormals( double *norm, double *rhs, PhoParamType photo, PointType *points,
				   int numpts )
{
	int		i, j, pt;
	double	bdot[3][7], epsilon[3], r, s, q, dx, dy, dz, sumres2=0;

	/* Zero out normal equations */
	for (i=1; i<=6; i++) {
		rhs[i] = 0;
		for (j=i; j<=6; j++) norm[INDUT(i,j)] = 0;
	}


	for (pt=0; pt<numpts; pt++) {
		/* Compute auxiliary variables */
		dx = points[pt].X - photo.xl;
		dy = points[pt].Y - photo.yl;
		dz = points[pt].Z - photo.zl;
		r = photo.m11*dx + photo.m12*dy + photo.m13*dz;
		s = photo.m21*dx + photo.m22*dy + photo.m23*dz;
		q = photo.m31*dx + photo.m32*dy + photo.m33*dz;

		/* Form b-dot (partials with respect to unknowns) */
		bdot[1][1] = (photo.f/(q*q)) *
			(r*(-photo.m33*dy + photo.m32*dz) - q*(-photo.m13*dy + photo.m12*dz));
		bdot[1][2] = (photo.f/(q*q)) *
			(r*(photo.cp*dx + photo.so*photo.sp*dy - photo.co*photo.sp*dz) -
			 q*(-photo.sp*photo.ck*dx + photo.so*photo.cp*photo.ck*dy -
				photo.co*photo.cp*photo.ck*dz));
		bdot[1][3] = -photo.f * s / q;
		bdot[1][4] = -photo.f * (r*photo.m31 - q*photo.m11) / (q*q);
		bdot[1][5] = -photo.f * (r*photo.m32 - q*photo.m12) / (q*q);
		bdot[1][6] = -photo.f * (r*photo.m33 - q*photo.m13) / (q*q);
		bdot[2][1] = (photo.f/(q*q)) *
			(s*(-photo.m33*dy + photo.m32*dz) - q*(-photo.m23*dy + photo.m22*dz));
		bdot[2][2] = (photo.f/(q*q)) *
			(s*(photo.cp*dx + photo.so*photo.sp*dy - photo.co*photo.sp*dz) -
			 q*(photo.sp*photo.sk*dx - photo.so*photo.cp*photo.sk*dy +
				photo.co*photo.cp*photo.sk*dz));
		bdot[2][3] = photo.f * r/ q;
		bdot[2][4] = -photo.f * (s*photo.m31 - q*photo.m21) / (q*q);
		bdot[2][5] = -photo.f * (s*photo.m32 - q*photo.m22) / (q*q);
		bdot[2][6] = -photo.f * (s*photo.m33 - q*photo.m23) / (q*q);

		/* Form epsilon (measured - computed) */
		epsilon[1] = points[pt].x + photo.f * r / q;
		epsilon[2] = points[pt].y + photo.f * s / q;

		/* Add contribution of this point to normal equations */
		for (i=1; i<=6; i++) {
			for (j=i; j<=6; j++)
				norm[INDUT(i,j)] += bdot[1][i]*bdot[1][j] + bdot[2][i]*bdot[2][j];
			rhs[i] += bdot[1][i]*epsilon[1] + bdot[2][i]*epsilon[2];
		}

		/* Accumulate subtotal of residuals-squared */
		sumres2 += epsilon[1]*epsilon[1] + epsilon[2]*epsilon[2];

		/* Save residuals */
		points[pt].xres = -epsilon[1];
		points[pt].yres = -epsilon[2];
	}

	/* Compute and return standard error of unit weight */
	if (numpts>3) return sqrt(sumres2/(2*numpts));
	else return 1.0e30; /* Zero degrees of freedom, s0 = infinity */
}





void AddCorrections( PhoParamType *Pphoto, double *rhs )
{
	Pphoto->omega += rhs[1];
	Pphoto->phi += rhs[2];
	Pphoto->kappa += rhs[3];
	Pphoto->xl += rhs[4];
	Pphoto->yl += rhs[5];
	Pphoto->zl += rhs[6];
}





void OutputResults( char *rootname, PhoParamType photo, PointType *points,
				   double *norm, double s0, int numpts )
{
	FILE	*outfile;
	char	filename[50];
	int		i;
	double	sumx=0, sumy=0, sumx2=0, sumy2=0;

	strcpy(filename, rootname);
	strcat(filename, ".out");
	outfile = fopen(filename, "w");
	printf("\nResults are in file %s\n", filename);

	fprintf(outfile, "Exterior orientation parameters:\n\n");
	fprintf(outfile, "Omega = %10.4lf  +- %6.4lf deg\n", photo.omega*180/M_PI,
		s0*sqrt(norm[INDUT(1,1)])*180/M_PI);
	fprintf(outfile, "Phi   = %10.4lf  +- %6.4lf deg\n", photo.phi*180/M_PI,
		s0*sqrt(norm[INDUT(2,2)])*180/M_PI );
	fprintf(outfile, "Kappa = %10.4lf  +- %6.4lf deg\n", photo.kappa*180/M_PI,
		s0*sqrt(norm[INDUT(3,3)])*180/M_PI );
	fprintf(outfile, "XL    = %10.4lf  +- %6.4lf\n", photo.xl,
		s0*sqrt(norm[INDUT(4,4)]) );
	fprintf(outfile, "YL    = %10.4lf  +- %6.4lf\n", photo.yl,
		s0*sqrt(norm[INDUT(5,5)]) );
	fprintf(outfile, "ZL    = %10.4lf  +- %6.4lf\n", photo.zl,
		s0*sqrt(norm[INDUT(6,6)]) );
	fprintf(outfile, "\n\n\nPhoto coordinate residuals:\n\n"
		"%8s %7s %7s\n", "point", "x-res", "y-res");
	for (i=0; i<numpts; i++) {
		fprintf(outfile, "%8s %7.4lf %7.4lf\n",
				points[i].name, points[i].xres, points[i].yres);
		sumx += points[i].xres;
		sumy += points[i].yres;
		sumx2 += points[i].xres * points[i].xres;
		sumy2 += points[i].yres * points[i].yres;
	}
	fprintf(outfile, "\n%8s %7.4lf %7.4lf\n%8s %7.4lf %7.4lf\n",
		"average", sumx/numpts, sumy/numpts,
		"RMS", sqrt(sumx2/numpts), sqrt(sumy2/numpts) );
	fprintf(outfile, "\n\n\nStandard error of unit weight: %7.4lf\n"
		"Degrees of freedom: %d\n", s0, 2*numpts-6);
	fclose(outfile);
}





void pause(void)
{
	printf("Press a key to end program.");
	getch();
}





/*
This program performs a space resection solution for a near vertical photo
by least squares. No weights are used for photo coordinates or object space
coordinates. All x and y photo coordinates must be corrected for principal
point offsets, lens distortions, and atmospheric refraction as needed.
Data file must be plain ASCII text and have a .dat extension. The format
of the file is as follows:

  Line 1: focal length
  Lines 2 through end of file:
		<point name> <photo x> <photo y> <ground X> <ground Y> <ground Z>

Sample data file:


152.916
tn08   86.421  -83.977  1268.1022  1455.0274  -14.3939
ts08 -100.916   92.582   732.1811   545.3437  -14.7009
re08  -98.322  -89.161  1454.5532   731.6659  -14.3509
rw08   78.812   98.123   545.2449  1268.2324  -14.6639
0000   -8.641    5.630  1000.0000  1000.0000  -14.4540


Output is sent to a file with the same root name and a .out extension.

*/
