/*  Three dimensional conformal coordinate transformation program
	Written October 25, 1998 (modified November 6, 2000) by Bon A. Dewitt
 	Modified November 14, 2017 by Arnaud Palha
	This source code is placed in the public domain and may be freely
	copied and modified.
	No warranties are expressed or implied. User accepts all risks
	associated with its use.
	See comments at the end of this source code file for brief
	description and instructions for use.
*/
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>




/* Maximum numbers of common and unknown points. Change as needed */
#define		MAXCOM		200
#define		MAXUNK		400


/* Macro to access upper triangular matrix stored as 1D array */
#define     INDUT(i,j)  ( (((j) * ((j)-1)) / 2) + i )

/* Macro to compute the square of a number */
#define		SQR(x)		((x)*(x))




typedef struct {
	char	name[16];
	double	xcon, ycon, zcon, xarb, yarb, zarb, xres, yres, zres;
} CommonPointType;

typedef struct {
	char	name[16];
	double	x, y, z;
} UnkPointType;

typedef struct {
    double  m11, m12, m13,
            m21, m22, m23,
            m31, m32, m33,
            so, sp, sk,
            co, cp, ck,
			st, ss, sa,
			ct, cs, ca;
} RotationMatrixType;




void ReadData( char *rootname, CommonPointType *ComPoints, UnkPointType *UnkPoints,
			  int *Pnumcom, int *Pnumunk );
void InitApprox( CommonPointType *ComPoints, int numcom, double *params, FILE *outf);
void solve( double *a, double *b, int n, int invflag );
void RotationMatrixOPK(double omega, double phi, double kappa,
					   RotationMatrixType *PRotMat);
void RotationMatrixTSA(double tilt, double swing, double azimuth,
					   RotationMatrixType *PRotMat);
double FormNormals(double *norm, double *rhs, RotationMatrixType RotMat,
				   CommonPointType *ComPoints, int numcom, double *params);
void FormAI( double AI[4][8], RotationMatrixType RotMat, double scale,
			double x, double y, double z );
void AddCorrections( double *rhs, double *params );
void PrintIter( double *rhs, double *params, CommonPointType *ComPoints, int numcom,
			   double s0, int iter, FILE *outf, char pnames[8][8] );
void PrintResults( FILE *outf, double *norm, double s0, double *params, int numunk,
				  UnkPointType *UnkPoints, char pnames[8][8] );
void pause(void);




int main(void)
{
	char				rootname[50], filename[50];
	int					numcom, numunk, iter=0, converge=0, diverge=0;
	double				s0=1.0e30, s0old, norm[29], rhs[8], params[8];
	CommonPointType		ComPoints[MAXCOM];
	UnkPointType		UnkPoints[MAXUNK];
	RotationMatrixType	RotMat;
	FILE				*outf;
	char				pnames[8][8] = {"", "scale", "omega", "phi", "kappa",
										"Tx", "Ty", "Tz"};

	ReadData( rootname, ComPoints, UnkPoints, &numcom, &numunk );

	/* Open output file (use .out extension) */
	strcpy(filename, rootname);
	strcat(filename, ".out");
	outf = fopen(filename, "w");

	/* Compute initial approximations for unknown parameters */
	InitApprox( ComPoints, numcom, params, outf );

	do {
		iter++;         /* Increment iteration count */
		s0old = s0;     /* Save former value of s0 for convergence test */
		/* Compute rotation matrix elements */
		RotationMatrixOPK( params[2], params[3], params[4], &RotMat);
		s0 = FormNormals(norm, rhs, RotMat, ComPoints, numcom, params);
		printf("ITERATION %d     S0 = %6.5lf\n", iter, s0);

		/* Check for convergence or divergence */
		if (fabs(s0old-s0)/s0 < 0.0001) converge=1;
		else if (s0 > s0old) diverge=1;

		/* Solve for corrections
		   If converged or diverged call for inverse also */
		solve(norm, rhs, 7, converge|diverge);
		AddCorrections( rhs, params );
		PrintIter( rhs, params, ComPoints, numcom, s0, iter, outf, pnames );
	} while (!converge && !diverge);

	PrintResults( outf, norm, s0, params, numunk, UnkPoints, pnames );

	fclose(outf);
	pause();
  return 0;
}




void ReadData( char *rootname, CommonPointType *ComPoints, UnkPointType *UnkPoints,
			  int *Pnumcom, int *Pnumunk )
{
    char    filename[50], tempstr[50];
    FILE    *infile;
    int     done;

	/* Get name of input file and open if it exists */
    printf("Enter root name of data file (.dat assumed) --> ");
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

	/* Read coordinates of common points until the end flag character, # */
	*Pnumcom = 0;
	done = 0;
	do {
		fscanf(infile, "%s", tempstr);
		if (strcmp(tempstr, "#") == 0) done = 1;
		else {
			if (*Pnumcom == MAXCOM) {
				printf("ERROR. More than %d common points in data file.\n", MAXCOM);
				pause();
				exit(1);
			}
			else {
				strcpy( ComPoints[*Pnumcom].name, tempstr );
				fscanf(infile, "%lf %lf %lf %lf %lf %lf", &(ComPoints[*Pnumcom].xarb),
					&(ComPoints[*Pnumcom].yarb), &(ComPoints[*Pnumcom].zarb),
					&(ComPoints[*Pnumcom].xcon), &(ComPoints[*Pnumcom].ycon),
					&(ComPoints[*Pnumcom].zcon) );
				(*Pnumcom)++;
			}
		}
	} while (!done);

	/* Be sure there are enough points */
	if (*Pnumcom < 3) {
		printf("ERROR. Less than 3 common points.\n");
		pause();
		exit(1);
	}

	/* Read coordinates of unknown points until the end flag character, # */
	*Pnumunk = 0;
	done = 0;
	do {
		fscanf(infile, "%s", tempstr);
		if (strcmp(tempstr, "#") == 0) done = 1;
		else {
			if (*Pnumunk == MAXUNK) {
				printf("ERROR. More than %d unknown points in data file.\n", MAXUNK);
				pause();
				exit(1);
			}
			else {
				strcpy( UnkPoints[*Pnumunk].name, tempstr );
				fscanf(infile, "%lf %lf %lf", &(UnkPoints[*Pnumunk].x),
					&(UnkPoints[*Pnumunk].y), &(UnkPoints[*Pnumunk].z) );
				(*Pnumunk)++;
			}
		}
	} while (!done);

	fclose(infile);

}




void InitApprox( CommonPointType *ComPoints, int numcom, double *params, FILE *outf )
/*	This function computes initial approximations for scale, omega, phi, kappa,
	Tx, Ty, and Tz. The method used is from:

	  Dewitt, B.A.: "Initial Approximations for the Three-Dimensional Conformal
			Coordinate Transformation, Photogrammetric Engineering and Remote
			Sensing, vol. 62, no. 1, p. 79, 1996.
*/
{
	int		numscale, i, j, pt1, pt2, pt3, ind1, ind2, ind3;
	double	sumscale, distcon, distarb, maxaltsq, dsq12, dsq23, dsq13, a2, b2, c2, h2,
			a_arb, b_arb, c_arb, a_con, b_con, c_con, arb_tilt, arb_az,
			con_tilt, con_az, x_arb1, y_arb1, x_arb2, y_arb2,
			x_con1, y_con1, x_con2, y_con2, azimuth_con, azimuth_arb, swing,
			txsum, tysum, tzsum;
	RotationMatrixType	ArbRotMat, ConRotMat, FullRotMat;

	/* Calculate approximation for scale using all pairs of points */
	numscale = 0;
	sumscale = 0;
	for (i=0; i<numcom-1; i++) {
		for (j=i+1; j<numcom; j++) {
			distcon = sqrt(	SQR(ComPoints[i].xcon - ComPoints[j].xcon) +
							SQR(ComPoints[i].ycon - ComPoints[j].ycon) +
							SQR(ComPoints[i].zcon - ComPoints[j].zcon)   );
			distarb = sqrt(	SQR(ComPoints[i].xarb - ComPoints[j].xarb) +
							SQR(ComPoints[i].yarb - ComPoints[j].yarb) +
							SQR(ComPoints[i].zarb - ComPoints[j].zarb)   );
			sumscale += distcon / distarb;
			numscale++;
		}
	}
	params[1] = sumscale / numscale;

	/*  Find geometrically strongest triangle of 3 points
		Strength is based on triangle with the largest altitude
		Altitude is defined as the perpendicular distance from the longest
		leg (or extension thereof) to the point not on the longest leg */
	maxaltsq = 0;
	for (ind1=0; ind1<numcom-2; ind1++) {
		for (ind2=ind1+1; ind2<numcom-1; ind2++) {
			dsq12 = SQR(ComPoints[ind1].xcon - ComPoints[ind2].xcon) +
					SQR(ComPoints[ind1].ycon - ComPoints[ind2].ycon) +
					SQR(ComPoints[ind1].zcon - ComPoints[ind2].zcon);
			for (ind3=ind2+1; ind3<numcom; ind3++) {
				dsq13 = SQR(ComPoints[ind1].xcon - ComPoints[ind3].xcon) +
						SQR(ComPoints[ind1].ycon - ComPoints[ind3].ycon) +
						SQR(ComPoints[ind1].zcon - ComPoints[ind3].zcon);
				dsq23 = SQR(ComPoints[ind2].xcon - ComPoints[ind3].xcon) +
						SQR(ComPoints[ind2].ycon - ComPoints[ind3].ycon) +
						SQR(ComPoints[ind2].zcon - ComPoints[ind3].zcon);

				if ((dsq12 >= dsq13) && (dsq12 >= dsq23)) {
					c2 = dsq12;
					a2 = dsq13;
					b2 = dsq23;
				}
				else {
					if ((dsq13 >= dsq12) && (dsq13 >= dsq23)) {
						c2 = dsq13;
						a2 = dsq12;
						b2 = dsq23;
					}
					else {
						c2 = dsq23;
						a2 = dsq12;
						b2 = dsq13;
					}
				}
				h2 = (2*(c2*(a2+b2)+a2*b2)-a2*a2-b2*b2-c2*c2)/(4*c2);
				if (h2 < 0) {
					printf("ERROR IN ALTITUDE COMPUTATION\n");
					pause();
					exit(1);
				}
				if (h2 > maxaltsq) {
					pt1 = ind1;
					pt2 = ind2;
					pt3 = ind3;
					maxaltsq = h2;
				}
			}
		}
	}

	/*  Compute coefficients of equation of plane through the 3 points
		in the arbitrary system */
	a_arb = ComPoints[pt1].yarb*ComPoints[pt2].zarb +
			ComPoints[pt2].yarb*ComPoints[pt3].zarb +
			ComPoints[pt3].yarb*ComPoints[pt1].zarb -
			ComPoints[pt1].yarb*ComPoints[pt3].zarb -
			ComPoints[pt2].yarb*ComPoints[pt1].zarb -
			ComPoints[pt3].yarb*ComPoints[pt2].zarb;
	b_arb = -ComPoints[pt1].xarb*ComPoints[pt2].zarb -
			ComPoints[pt2].xarb*ComPoints[pt3].zarb -
			ComPoints[pt3].xarb*ComPoints[pt1].zarb +
			ComPoints[pt1].xarb*ComPoints[pt3].zarb +
			ComPoints[pt2].xarb*ComPoints[pt1].zarb +
			ComPoints[pt3].xarb*ComPoints[pt2].zarb;
	c_arb = ComPoints[pt1].xarb*ComPoints[pt2].yarb +
			ComPoints[pt2].xarb*ComPoints[pt3].yarb +
			ComPoints[pt3].xarb*ComPoints[pt1].yarb -
			ComPoints[pt1].xarb*ComPoints[pt3].yarb -
			ComPoints[pt2].xarb*ComPoints[pt1].yarb -
			ComPoints[pt3].xarb*ComPoints[pt2].yarb;


	/*  Compute coefficients of equation of plane through the 3 points
		in the control system */
	a_con = ComPoints[pt1].ycon*ComPoints[pt2].zcon +
			ComPoints[pt2].ycon*ComPoints[pt3].zcon +
			ComPoints[pt3].ycon*ComPoints[pt1].zcon -
			ComPoints[pt1].ycon*ComPoints[pt3].zcon -
			ComPoints[pt2].ycon*ComPoints[pt1].zcon -
			ComPoints[pt3].ycon*ComPoints[pt2].zcon;
	b_con = -ComPoints[pt1].xcon*ComPoints[pt2].zcon -
			ComPoints[pt2].xcon*ComPoints[pt3].zcon -
			ComPoints[pt3].xcon*ComPoints[pt1].zcon +
			ComPoints[pt1].xcon*ComPoints[pt3].zcon +
			ComPoints[pt2].xcon*ComPoints[pt1].zcon +
			ComPoints[pt3].xcon*ComPoints[pt2].zcon;
	c_con = ComPoints[pt1].xcon*ComPoints[pt2].ycon +
			ComPoints[pt2].xcon*ComPoints[pt3].ycon +
			ComPoints[pt3].xcon*ComPoints[pt1].ycon -
			ComPoints[pt1].xcon*ComPoints[pt3].ycon -
			ComPoints[pt2].xcon*ComPoints[pt1].ycon -
			ComPoints[pt3].xcon*ComPoints[pt2].ycon;

	/* Compute tilt & azimuth of plane in arbitrary system */
	arb_tilt = atan2( c_arb, hypot(a_arb, b_arb) ) + M_PI/2;
	arb_az = atan2( a_arb, b_arb );
	/* Compute tilt & azimuth of plane in control system */
	con_tilt = atan2( c_con, hypot(a_con, b_con) ) + M_PI/2;
	con_az = atan2( a_con, b_con );
	/*  Compute the corresponding rotation matrices with swing = 0 */
	RotationMatrixTSA(arb_tilt, 0.0, arb_az, &ArbRotMat);
	RotationMatrixTSA(con_tilt, 0.0, con_az, &ConRotMat);
	/*  Rotate arbitrary and control coordinates for points 1 and 2
		to get a line in the arbitrary system and the corresponding
		line in the control system.  Don't need to rotate Z because
		with these rotations all of the z values must be the same. */
	x_arb1 =	ArbRotMat.m11*ComPoints[pt1].xarb +
				ArbRotMat.m12*ComPoints[pt1].yarb +
				ArbRotMat.m13*ComPoints[pt1].zarb;
	y_arb1 =	ArbRotMat.m21*ComPoints[pt1].xarb +
				ArbRotMat.m22*ComPoints[pt1].yarb +
				ArbRotMat.m23*ComPoints[pt1].zarb;
	x_arb2 =	ArbRotMat.m11*ComPoints[pt2].xarb +
				ArbRotMat.m12*ComPoints[pt2].yarb +
				ArbRotMat.m13*ComPoints[pt2].zarb;
	y_arb2 =	ArbRotMat.m21*ComPoints[pt2].xarb +
				ArbRotMat.m22*ComPoints[pt2].yarb +
				ArbRotMat.m23*ComPoints[pt2].zarb;
	x_con1 =	ConRotMat.m11*ComPoints[pt1].xcon +
				ConRotMat.m12*ComPoints[pt1].ycon +
				ConRotMat.m13*ComPoints[pt1].zcon;
	y_con1 =	ConRotMat.m21*ComPoints[pt1].xcon +
				ConRotMat.m22*ComPoints[pt1].ycon +
				ConRotMat.m23*ComPoints[pt1].zcon;
	x_con2 =	ConRotMat.m11*ComPoints[pt2].xcon +
				ConRotMat.m12*ComPoints[pt2].ycon +
				ConRotMat.m13*ComPoints[pt2].zcon;
	y_con2 =	ConRotMat.m21*ComPoints[pt2].xcon +
				ConRotMat.m22*ComPoints[pt2].ycon +
				ConRotMat.m23*ComPoints[pt2].zcon;
	/* Get swing by subtracting azimuths */
	azimuth_con = atan2( x_con2-x_con1, y_con2-y_con1 );
	azimuth_arb = atan2( x_arb2-x_arb1, y_arb2-y_arb1 );
	swing = azimuth_con - azimuth_arb;
	RotationMatrixTSA(arb_tilt, swing, arb_az, &ArbRotMat);
	/*  Now compute (ConRotMat:transpose * ArbRotMat):transpose
		This is overall rotation matrix */
	FullRotMat.m11 = ConRotMat.m11*ArbRotMat.m11 + ConRotMat.m21*ArbRotMat.m21 +
					 ConRotMat.m31*ArbRotMat.m31;
	FullRotMat.m21 = ConRotMat.m11*ArbRotMat.m12 + ConRotMat.m21*ArbRotMat.m22 +
					 ConRotMat.m31*ArbRotMat.m32;
	FullRotMat.m31 = ConRotMat.m11*ArbRotMat.m13 + ConRotMat.m21*ArbRotMat.m23 +
					 ConRotMat.m31*ArbRotMat.m33;
	FullRotMat.m12 = ConRotMat.m12*ArbRotMat.m11 + ConRotMat.m22*ArbRotMat.m21 +
					 ConRotMat.m32*ArbRotMat.m31;
	FullRotMat.m22 = ConRotMat.m12*ArbRotMat.m12 + ConRotMat.m22*ArbRotMat.m22 +
					 ConRotMat.m32*ArbRotMat.m32;
	FullRotMat.m32 = ConRotMat.m12*ArbRotMat.m13 + ConRotMat.m22*ArbRotMat.m23 +
					 ConRotMat.m32*ArbRotMat.m33;
	FullRotMat.m13 = ConRotMat.m13*ArbRotMat.m11 + ConRotMat.m23*ArbRotMat.m21 +
					 ConRotMat.m33*ArbRotMat.m31;
	FullRotMat.m23 = ConRotMat.m13*ArbRotMat.m12 + ConRotMat.m23*ArbRotMat.m22 +
					 ConRotMat.m33*ArbRotMat.m32;
	FullRotMat.m33 = ConRotMat.m13*ArbRotMat.m13 + ConRotMat.m23*ArbRotMat.m23 +
					 ConRotMat.m33*ArbRotMat.m33;

	/* Compute omega, phi, kappa from rotation matrix */
	params[3] = asin( FullRotMat.m31 ); /* phi */
	if (fabs(FullRotMat.m31) < 1.0) { /* typical case, compute omega and kappa */
		params[2] = atan2( -FullRotMat.m32, FullRotMat.m33 ); /* omega */
		params[4] = atan2( -FullRotMat.m21, FullRotMat.m11 ); /*kappa */
	}
	else { /* omega and kappa are undefined, so define them */
		params[2] = 0.0; /* omega */
		params[4] = atan2( FullRotMat.m12, FullRotMat.m22 ); /* kappa */
	}

	/* Compute average Tx, Ty, and Tz using all common points */
	txsum = 0.0;
	tysum = 0.0;
	tzsum = 0.0;
	for (i=0; i<numcom; i++) {
		txsum += ComPoints[i].xcon - params[1] * (FullRotMat.m11*ComPoints[i].xarb +
			FullRotMat.m21*ComPoints[i].yarb + FullRotMat.m31*ComPoints[i].zarb);
		tysum += ComPoints[i].ycon - params[1] * (FullRotMat.m12*ComPoints[i].xarb +
			FullRotMat.m22*ComPoints[i].yarb + FullRotMat.m32*ComPoints[i].zarb);
		tzsum += ComPoints[i].zcon - params[1] * (FullRotMat.m13*ComPoints[i].xarb +
			FullRotMat.m23*ComPoints[i].yarb + FullRotMat.m33*ComPoints[i].zarb);
	}
	params[5] = txsum/numcom;
	params[6] = tysum/numcom;
	params[7] = tzsum/numcom;

	/* Print initial approximations to output file */
	fprintf(outf, "Initial approximations:\n");
	fprintf(outf, "Scale = %9.5lf\n", params[1]);
	fprintf(outf, "Omega = %9.4lf degrees\n", params[2]*180/M_PI);
	fprintf(outf, "Phi   = %9.4lf degrees\n", params[3]*180/M_PI);
	fprintf(outf, "Kappa = %9.4lf degrees\n", params[4]*180/M_PI);
	fprintf(outf, "TX    = %9.3lf\n", params[5]);
	fprintf(outf, "TY    = %9.3lf\n", params[6]);
	fprintf(outf, "TZ    = %9.3lf\n\n\n\n", params[7]);
}




void solve( double *a, double *b, int n, int invflag )
/* Solution and inverse by recursive partitioning for upper triangular
   normal equation matrix. When complete, array 'b' is overwritten by
   solution. If 'invflag' is true (i.e. non-zero) then array 'a' is
   overwritten by the inverse also.
*/
{
    int     piv, j, i;
    double  r, *s;

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





void RotationMatrixOPK(double omega, double phi, double kappa,
					   RotationMatrixType *PRotMat)
/* Calculate rotation matrix based on omega, phi, and kappa */
{
    /* Compute trig functions */
    PRotMat->so = sin(omega);
    PRotMat->co = cos(omega);
    PRotMat->sp = sin(phi);
    PRotMat->cp = cos(phi);
    PRotMat->sk = sin(kappa);
    PRotMat->ck = cos(kappa);

    /* Compute rotation matrix elements */
    PRotMat->m11 = PRotMat->cp * PRotMat->ck;
    PRotMat->m12 = PRotMat->so * PRotMat->sp * PRotMat->ck + PRotMat->co * PRotMat->sk;
    PRotMat->m13 = -PRotMat->co * PRotMat->sp * PRotMat->ck + PRotMat->so * PRotMat->sk;
    PRotMat->m21 = -PRotMat->cp * PRotMat->sk;
    PRotMat->m22 = -PRotMat->so * PRotMat->sp * PRotMat->sk + PRotMat->co * PRotMat->ck;
    PRotMat->m23 = PRotMat->co * PRotMat->sp * PRotMat->sk + PRotMat->so * PRotMat->ck;
    PRotMat->m31 = PRotMat->sp;
    PRotMat->m32 = -PRotMat->so * PRotMat->cp;
    PRotMat->m33 = PRotMat->co * PRotMat->cp;
}




void RotationMatrixTSA(double tilt, double swing, double azimuth,
					   RotationMatrixType *PRotMat)
/* Calculate rotation matrix based on tile, swing, and azimuth */
{
    /* Compute trig functions */
    PRotMat->st = sin(tilt);
    PRotMat->ct = cos(tilt);
    PRotMat->ss = sin(swing);
    PRotMat->cs = cos(swing);
    PRotMat->sa = sin(azimuth);
    PRotMat->ca = cos(azimuth);

    /* Compute rotation matrix elements */
    PRotMat->m11 = -PRotMat->ca * PRotMat->cs - PRotMat->sa * PRotMat->ct * PRotMat->ss;
    PRotMat->m12 = PRotMat->sa * PRotMat->cs - PRotMat->ca * PRotMat->ct * PRotMat->ss;
    PRotMat->m13 = -PRotMat->st * PRotMat->ss;
    PRotMat->m21 = PRotMat->ca * PRotMat->ss - PRotMat->sa * PRotMat->ct * PRotMat->cs;
    PRotMat->m22 = -PRotMat->sa * PRotMat->ss - PRotMat->ca * PRotMat->ct * PRotMat->cs;
    PRotMat->m23 = -PRotMat->st * PRotMat->cs;
    PRotMat->m31 = -PRotMat->sa * PRotMat->st;
    PRotMat->m32 = -PRotMat->ca * PRotMat->st;
    PRotMat->m33 = PRotMat->ct;
}




double FormNormals(double *norm, double *rhs, RotationMatrixType RotMat,
				   CommonPointType *ComPoints, int numcom, double *params)
{
	int		i, j, k, pt;
	double	AI[4][8], LI[4], sumres2=0.0;

	/* Zero out normal equations */
	for (i=1; i<=7; i++) {
		rhs[i] = 0;
		for (j=i; j<=7; j++) norm[INDUT(i,j)] = 0;
	}

	for (pt=0; pt<numcom; pt++) {
		/* Form partial derivative terms for Taylor's series approximation */
		FormAI( AI, RotMat, params[1], ComPoints[pt].xarb, ComPoints[pt].yarb,
			ComPoints[pt].zarb );
		LI[1] = ComPoints[pt].xcon - params[1]*AI[1][1] - params[5];
		ComPoints[pt].xres = -LI[1];
		LI[2] = ComPoints[pt].ycon - params[1]*AI[2][1] - params[6];
		ComPoints[pt].yres = -LI[2];
		LI[3] = ComPoints[pt].zcon - params[1]*AI[3][1] - params[7];
		ComPoints[pt].zres = -LI[3];

		/* Accumulate sum of squares of residuals */
		sumres2 += LI[1]*LI[1] + LI[2]*LI[2] + LI[3]*LI[3];

		/* Add contributions from this point to normal equations */
		for (i=1; i<=7; i++) {
			for (k=1; k<=3; k++) rhs[i] += AI[k][i] * LI[k];
			for (j=i; j<=7; j++) {
				for (k=1; k<=3; k++) norm[INDUT(i,j)] += AI[k][i] * AI[k][j];
			}
		}
	}

	/* Return computed estimate of the standard error of unit weight */
	return ( sqrt(sumres2/(3*numcom-7)) );
}




void FormAI( double AI[4][8], RotationMatrixType RotMat, double scale,
			double x, double y, double z )
/* Calculate partial derivative terms for the point */
{
		AI[1][1] =	RotMat.m11*x + RotMat.m21*y + RotMat.m31*z;
		AI[2][1] =	RotMat.m12*x + RotMat.m22*y + RotMat.m32*z;
		AI[3][1] =	RotMat.m13*x + RotMat.m23*y + RotMat.m33*z;
		AI[1][2] = 0.0;
		AI[2][2] = -scale*AI[3][1];
		AI[3][2] = scale*AI[2][1];
		AI[1][3] = scale * (-RotMat.sp*RotMat.ck*x + RotMat.sp*RotMat.sk*y +
							RotMat.cp*z);
		AI[2][3] = scale * (RotMat.so*RotMat.cp*RotMat.ck*x -
			RotMat.so*RotMat.cp*RotMat.sk*y + RotMat.so*RotMat.sp*z);
		AI[3][3] = scale * (-RotMat.co*RotMat.cp*RotMat.ck*x +
			RotMat.co*RotMat.cp*RotMat.sk*y - RotMat.co*RotMat.sp*z);
		AI[1][4] = scale * (RotMat.m21*x - RotMat.m11*y);
		AI[2][4] = scale * (RotMat.m22*x - RotMat.m12*y);
		AI[3][4] = scale * (RotMat.m23*x - RotMat.m13*y);
		AI[1][5] = 1.0;
		AI[2][5] = 0.0;
		AI[3][5] = 0.0;
		AI[1][6] = 0.0;
		AI[2][6] = 1.0;
		AI[3][6] = 0.0;
		AI[1][7] = 0.0;
		AI[2][7] = 0.0;
		AI[3][7] = 1.0;
}




void AddCorrections( double *rhs, double *params )
{
	int		i;

	for (i=1; i<=7; i++) params[i] += rhs[i];

}




void PrintIter( double *rhs, double *params, CommonPointType *ComPoints, int numcom,
			   double s0, int iter, FILE *outf, char pnames[8][8] )
{
	int		i;

	fprintf(outf, "Iteration: %d\n\n", iter);
	/* Print corrections and updated approximations for the 7 parameters */
	fprintf(outf, "%5s  %10s  %10s\n", "Param", "Correction", "New Value");
	fprintf(outf, "%5s  %10.5lf  %10.5lf\n", pnames[1], rhs[1], params[1]);
	for (i=2; i<=4; i++)
		fprintf(outf, "%5s  %10.4lfd %10.4lfd\n", pnames[i],
				rhs[i]*180/M_PI, params[i]*180/M_PI);
	for (i=5; i<=7; i++)
		fprintf(outf, "%5s  %10.3lf  %10.3lf\n", pnames[i], rhs[i], params[i]);
	/* Print estimates for residuals based on current iteration */
	fprintf(outf, "\n\n\nResiduals:\n\n%8s  %7s %7s %7s\n", "Point",
		"X res", "Y res", "Z res");
	for (i=0; i<numcom; i++)
		fprintf(outf, "%8s  %7.3lf %7.3lf %7.3lf\n", ComPoints[i].name,
			ComPoints[i].xres, ComPoints[i].yres, ComPoints[i].zres );
	fprintf(outf,"\n\n");
	fprintf(outf, "Standard Error of Unit Weight: %6.5lf\n\n\n", s0);
}




void PrintResults( FILE *outf, double *norm, double s0, double *params, int numunk,
				  UnkPointType *UnkPoints, char pnames[8][8] )
{
	int					i, r, c, k;
	RotationMatrixType	RotMat;
	double				x, y, z, AI[4][8], AIQ[4][8], AIQAT[4][4];

	/* Print out transformation parameters and standard deviations */
	fprintf(outf, "Final Results:\n\n");
	fprintf(outf, "%5s  %10s  %10s\n", "Param", "Value", "Stan.Dev.");
	fprintf(outf, "%5s  %10.5lf  %10.5lf\n", pnames[1], params[1],
			s0*sqrt(norm[INDUT(1,1)]));
	for (i=2; i<=4; i++)
		fprintf(outf, "%5s  %10.4lfd %10.4lfd\n", pnames[i],
				params[i]*180/M_PI, s0*sqrt(norm[INDUT(i,i)])*180/M_PI);
	for (i=5; i<=7; i++)
		fprintf(outf, "%5s  %10.3lf  %10.3lf\n", pnames[i], params[i],
				s0*sqrt(norm[INDUT(i,i)]));

	/* Recompute rotation matrix for converged omega, phi, and kappa */
	RotationMatrixOPK( params[2], params[3], params[4], &RotMat );

	/* Compute and print transformed coordinates and standard deviations for
		additional unknown points */
	fprintf(outf, "\n\n\nTransformed Points:\n\n%8s  %9s %9s %8s %6s %6s %6s\n",
		"Point", "X    ", "Y    ", "Z    ", "SDev.X", "SDev.Y", "SDev.Z");
	for (i=0; i<numunk; i++) {
		FormAI( AI, RotMat, params[1], UnkPoints[i].x, UnkPoints[i].y, UnkPoints[i].z );
		x =	params[1]*AI[1][1] + params[5];
		y =	params[1]*AI[2][1] + params[6];
		z =	params[1]*AI[3][1] + params[7];
		for (r=1; r<=3; r++) {
			for (c=1; c<=7; c++) {
				AIQ[r][c] = 0;
				for (k=1; k<=c; k++) AIQ[r][c] += AI[r][k]*norm[INDUT(k,c)];
				for (k=c+1; k<=7; k++) AIQ[r][c] += AI[r][k]*norm[INDUT(c,k)];
			}
		}
		for (r=1; r<=3; r++) {
			for (c=r; c<=3; c++) {
				AIQAT[r][c] = 0;
				for (k=1; k<=7; k++) AIQAT[r][c] += AIQ[r][k]*AI[c][k];
			}
		}
		fprintf(outf, "%8s  %9.3lf %9.3lf %8.3lf %6.3lf %6.3lf %6.3lf\n",
			UnkPoints[i].name, x, y, z, s0*sqrt(AIQAT[1][1]), s0*sqrt(AIQAT[2][2]),
			s0*sqrt(AIQAT[3][3]) );
	}
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
This program performs a three-dimensional conformal coordinate transformation
by least squares. No weights are used for the coordinates. Data file must be
plain ASCII text and have a .dat extension. Control must be three-dimensional.
The format of the file is as follows:

  First group of lines - common points
  <point name> <x arb> <y arb> <z arb> <x control> <y control> <z control>
  Last line in first group:
  #
  Second group of lines - unknown points to be transformed
  <point name> <x arb> <y arb> <z arb>
  Last line in second group:
  #


Sample data file:


a   390.35  499.63  469.43  607.54  501.63  469.09
b   371.68  630.84   81.25  589.98  632.36   82.81
c   425.65  419.07   82.49  643.65  421.28   83.50
d   410.50  438.31   81.13  628.58  440.51   82.27
e   448.22  295.83   97.79  666.27  298.16   98.29
f   414.60  709.39  101.77  632.59  710.62  103.01
#
1   611.37  498.98  470.45
2   637.49  323.67  85.67
3   573.32  401.51  84.48
4   647.00  373.97  83.76
5   533.51  285.01  87.13
#


Output is sent to a file with the same root name and a .out extension.

*/
