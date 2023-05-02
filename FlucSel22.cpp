/* 

This program estimates the features of a multiple-allele locus in a finite population, allowing for selection and reversible mutation.
The goal is to estimate the steady-state distribution of alternative allele types over time. The equilibrium long-term result is a balance between the forces of mutation, selection, and random genetic drift. 

***** This version has two types of sites simultaneously present: major with large fitness effects, and minor with neutral effects.

***** This version allows the use of arbitrary numbers of sites and stength of selection, and has fitness function of the form exp(-s * number of deleterious alleles).

***** The selection coefficient is allowed to vary temporally, with standard deviation and mean value defined at the top.

***** In this version, a random variate for the selection coefficient is drawn each generation for BOTH alleles.



The actual population size is assumed to be equal to the effective size (N), which is constant in time.

The population experiences sequential episodes of mutation, selection, and random genetic drit.

Haploidy is assumed, and there is no recombination.

Allele designation:	the general code applies to any trait with a series of sites, each with biallelic states -/+.
	-----...., for the most deleterious allele.
	+++++...., for the most fit allele.

But because there is complete linkage, order does not matter, and a fraction of sites is allocated to the minor vs. major component.

The population is regularly censused after reproduction, at which point the population has size N.

All mutation and selection processes are treated deterministically, prior to random sampling.

Mutation:
	There are just two types of mutations at each type of site: - to +, u01 (beneficial); and + to -, u10 (deleterious); the same rates are assumed at both loci.
	
Selection:
	Fitness is determined by a function that has to be set internally. 
	
Running averages of the population features are kept track of: probabilities of being fixed in the alternative monomorphic states, and the allele-frequency distribution conditional on being polymorphic.

The run starts with an allele frequency distribution equal to the expectation under the sequential model.

After a burnin, statistics are then taken at intervals of Ne/xinc generations, with a total of ngen sampling intervals.

The rates of mutation and strength of selection are defined by input parameters entered immediately below.

NOTE THAT A SCALING FACTOR (kfac) CAN BE UTILIZED TO SCALE UP THE SELECTION COEFFICIENTS. IF THE NE IS SCALED DOWN BY THE SAME FACTOR AND THE MUTATION RATES UP BY THE SAME FACTOR,
THIS ENABLES THE OVERALL PROGRAM TO RUN FASTER, AS ALL IS A FUNCTION OF THE PRODUCTS NE*S AND NE*U. THESE ARE SET SO THAT MODIFIED NE > 1000, AND MODIFIED U AND S < 0.1.
SETTING KFAC = 1 ELIMINATES THIS SCALING. 

A SERIES OF 21 COMBINATIONS OF POPULATION SIZE (efpopn) AND MUTATION RATE (delmutrate) (USING SCALING FROM LYNCH ET AL. 2016), COVERING THE FULL RANGE OF NATURAL POPULATION VALUES, IS GIVEN INTERNALLY. 

*/



/* ********************************************************************************************************************** */

#define mutbeta	    1.00 					/* ratio of beneficial to deleterious mutation rates */

#define ellmaj		1   				    /* total number of sites for the trait with major effects */

#define ellmin		1					    /* total number of sites for the trait with minor effects */

#define scomaj		0.01				    /* mean selection coefficient for major-effect alleles */

#define sdsmaj		0.1						/* temporal standard deviation of s at the major-effect loci, + allele */

#define scomin		0.000			        /* selection coefficient for minor-effect alleles; here treated as a neutral marker with scomin = 0.0 */

#define sdsmin		0.000				    /* temporal standard deviation of s at the major-effect loci, - allele */

#define xinc		20					    /* statistics to be recorded every ne/xinc generations */

#define burnin		10000				    /* number of initial burn-in sampling increments, ignored in the statistics */

#define tintprint	20000				    /* number of sampling increments between screen printing */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_randist.h>
#include    <string.h>


/* ********************************************************************************************************** */



/* point to the output file */

FILE *stream;
char filename[100];



/* Set up the parallel runs. */

void main(int argc, char *argv[])
{

int f0, f1;
if (argc > 1)
{
    f0 = atoi(argv[1]);
    f1=f0; 
    sprintf(filename, "dataoutd21a_%d.txt", f0); 
}
else
{
    f0 = 1;
    f1 = 21; 
   
    sprintf(filename, "dataoutd21a.txt");
}


/* Set up the binomial random number generator. */

static gsl_rng* rand_new;                                                       
gsl_rng_env_setup();                                                            
if(!rand_new)                                                                   
{                                                                               
    rand_new = gsl_rng_alloc(gsl_rng_taus2);                                    
    gsl_rng_set(rand_new, time(NULL));                                          
}     


/* Set up the normal random number generator. */

static gsl_rng* rand_new2;
gsl_rng_env_setup();
if (!rand_new2)
{
	rand_new2 = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rand_new2, time(NULL));
}



   

/* ***************************************************************************************** */



/* MAIN BODY OF PROGRAM. */


int ig, igmaj, igmin;									/* counters for the classes, with 0 meaning no + alleles, and ell meaning all + alleles */

long igen;												/* generation counter */

double tint;											/* interval number for statistics */

double ngens;                                           /* time iterations in run */

double difmin, difmaj;								    /* deviation between number of + alleles and the optimum, for minor- and major-effect sites */

double meanfit;									        /* mean fitness */

double gmaxrate;									    /* maximum growth rate under additive model */

double phi;									        	/* ratio of beneficial and deleterious mutation rates per site */

double selcomaj1, selcomaj2, selcomaj, selcomin;		/* selection coefficients associated with major- and minor-effect loci */

double mut0min, mutdown0min, mutemin, mutupemin;	    /* parameters for mutational transition frequencies between classes */
double mut0maj, mutdown0maj, mutemaj, mutupemaj;

long efpopn[40];										/* effective population sizes to run -- NOTE THESE ARE SCALED DOWN ACCORDING TO KFAC TO INCREASE RUN SPEED */
double delmutrate[40];								    /* deleterious mutation rates to run -- NOTE THESE ARE SCALED UP ACCORDING TO KFAC TO INCREASE RUN SPEED */

double rlng[40];
    
int itera;												/* counter for population-size/mutation-rate iterations */

long ne;												/* effective population size, from efpopn[] */
double u10, u01;									    /* deleterious and beneficial mutation rates, from delmutrate[] */
int kfac;												/* scaling factor for speeding up runs, from scalef[] */

int newlowmin, newhighmin, newlowmaj, newhighmaj;		/* running settings for upper and lower genoytpic states for minor and major allele counts */
int lowigmin, highigmin, lowigmaj, highigmaj;
int lowmutmin, highmutmin, lowmutmaj, highmutmaj;

double sump;										    /* sum of frequencies */

int startcheckmin, startcheckmaj;						/* checkers for new low and high frequency classes */

double pp;											    /* probability associated with the binomial for drift */
long ntot;												/* integer associated with the binomial for drift */
long draw;												/* drift drawn from the binomial */
double epoi, rnum;									    /* terms for Poisson draws */

long increment;											/* increment between surveys */
long tcount;											/* counter for printing to screen to monitor simulation progress */
long counter;											/* initial number of surveys to be skipped */

double meanw, meanigmin, meanigmaj;				        /* generational means for fitness and allelic class */

double totw, totwsq;								    /* summations for grand means and variances */
double totigmin, totigsqmin, totwvarigmin;			    /* summations for grand means and variances */
double totigmaj, totigsqmaj, totwvarigmaj;

double tots, totssq;									/* summations for means and variances of selection coefficients */

double grandmeanw, varw;							    /* mean and variance of fitness */

double grandmeans, vars;								/* mean and variance of selection coefficient (*/

double grandmeanigmin, varigmin;					    /* mean and variances of numbers of minor and major alleles*/
double grandmeanigmaj, varigmaj;

double ssqigmin, wvarigmin;						        /* sum of squares and within-generation variance of trait */
double ssqigmaj, wvarigmaj;

double meanperform;								        /* mean relative performance under additive model */

double totmin0, totmaj0, meanminhet, meanmajhet;	    /* statistics for computing the long-term mean heterozygosity at minor and major sites (only applies when ellmin and ellmaj are set to 1) */
double neuthet, majhet;
double neestimate;

double fracpine;                                        /* effective population size based on pis relative to input ne */ 

double meanmin, meanmaj;                                /* grand average frequencies for minor and major alleles */

double kfacn, kfacs, kfacu;                             /* possible scaling factors based on population size, selection coefficients, and mutation rates */

double ssdismin[1002], ssdismaj[1002];					/* steady-state distributions for the + alleles at the minor (neutral and major (selected) loci */
int binmin1, binmaj1;

double foldedmin, foldedmaj;                            /* folded steady-state distribution */



double *mutransupmin = (double*)malloc((ellmin+2) * sizeof(double));
double *mutransdownmin = (double*)malloc((ellmin+2) * sizeof(double));
double *mutrans0min = (double*)malloc((ellmin+2) * sizeof(double));

double *mutransupmaj = (double*)malloc((ellmaj + 2) * sizeof(double));
double *mutransdownmaj = (double*)malloc((ellmaj + 2) * sizeof(double));
double *mutrans0maj = (double*)malloc((ellmaj + 2) * sizeof(double));

double *sumfreqmaj = (double*)malloc((ellmaj + 2) * sizeof(double)); 
double *sumfreqmin = (double*)malloc((ellmin + 2) * sizeof(double)); 
double *meanpfreqmaj = (double*)malloc((ellmaj + 2) * sizeof(double)); 
double *meanpfreqmin = (double*)malloc((ellmin + 2) * sizeof(double));


double **relw = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ relw[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **p0 = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ p0[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pmutm1 = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pmutm1[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pmutm = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pmutm[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **psel = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ psel[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }

double **pgtypexp = (double **)malloc((ellmin + 2) * (sizeof(double *)));
	for (ig = 0; ig < (ellmin + 2); ig++)
		{ pgtypexp[ig] = (double *)malloc((ellmaj + 2) * (sizeof (double))); }





/* Open the output file. */

remove("dataoutd21a.txt ");


/* Effective population sizes to run */

efpopn[21] = 100000000;
efpopn[20] = 63095734;
efpopn[19] = 39810717;
efpopn[18] = 25118864;
efpopn[17] = 15848932;
efpopn[16] = 10000000;
efpopn[15] = 6309573;
efpopn[14] = 3981072;
efpopn[13] = 2511886;
efpopn[12] = 1584893;
efpopn[11] = 1000000;
efpopn[10] = 630957;
efpopn[9] = 398107;
efpopn[8] = 251189;
efpopn[7] = 158489;
efpopn[6] = 100000;
efpopn[5] = 63096;
efpopn[4] = 39811;
efpopn[3] = 25119;
efpopn[2] = 15849;
efpopn[1] = 10000;

efpopn[1]= 1000000;




/* Associated mutation rates */

delmutrate[21] = 0.00000000025;
delmutrate[20] = 0.00000000035;
delmutrate[19] = 0.00000000048;
delmutrate[18] = 0.00000000066;
delmutrate[17] = 0.00000000091;
delmutrate[16] = 0.00000000126;
delmutrate[15] = 0.00000000174;
delmutrate[14] = 0.00000000240;
delmutrate[13] = 0.00000000331;
delmutrate[12] = 0.00000000457;
delmutrate[11] = 0.00000000631;
delmutrate[10] = 0.00000000871;
delmutrate[9] = 0.0000000120;
delmutrate[8] = 0.0000000166;
delmutrate[7] = 0.0000000229;
delmutrate[6] = 0.0000000316;
delmutrate[5] = 0.0000000436;
delmutrate[4] = 0.0000000603;
delmutrate[3] = 0.0000000832;
delmutrate[2] = 0.000000115;
delmutrate[1] = 0.000000158; 

delmutrate[1] = 0.0000000057;




/* Number of sampling increments in run; each increment is (ne/xinc) generations */

rlng[21] = 1000000.0;
rlng[20] = 1000000.0;
rlng[19] = 2000000.0;
rlng[18] = 2000000.0;
rlng[17] = 2000000.0;
rlng[16] = 4000000.0;
rlng[15] = 4000000.0;
rlng[14] = 4000000.0;
rlng[13] = 4000000.0;
rlng[12] = 4000000.0;
rlng[11] = 4000000.0;
rlng[10] = 8000000.0;
rlng[9] = 8000000.0;
rlng[8] = 8000000.0;
rlng[7] = 8000000.0;
rlng[6] = 8000000.0;
rlng[5] = 8000000.0;
rlng[4] = 8000000.0;
rlng[3] = 20000000.0;
rlng[2] = 20000000.0;
rlng[1] = 20000000.0;

rlng[1] = 5000000.0;





gmaxrate = (((double)ellmin) * scomin) + (((double)ellmaj) * scomaj);



double start, stop, time;                                                   


for (itera = f0; itera <= f1; ++itera) {							/* Start iterations over the set of population sizes and mutation rates. */

stream=fopen(filename, "a");		



/* Set the run length. */
        
    ngens = rlng[itera];

	ne = efpopn[itera];											/* effective population size */

	u10 = delmutrate[itera];									/* deleterious and beneficial mutation rates */
	u01 = u10 * mutbeta;


	/* kfac allows for the possibility of jointly scaling down Ne and scaling up s and u to speed up runs; set = 1 for eliminate this treatment. */

    kfacn = ((double) ne) / 1000.0;
    
    kfacs = 0.1 / (scomaj + (2.0*sdsmaj)); 
    
    kfacu = 0.1 / (u10 * ellmaj);

    if (kfacn < 1.0) {
        kfacn = 1.0; }
    if (kfacs < 1.0) {
        kfacs = 1.0; }
    if (kfacu < 1.0) {
        kfacu = 1.0; }

    if ((kfacn < kfacs) && (kfacn < kfacu)) {
        kfac = int(kfacn); }
    if ((kfacs < kfacn) && (kfacs < kfacu)) {
        kfac = int(kfacs); }
    if ((kfacu < kfacs) && (kfacu < kfacn)) {
        kfac = int(kfacu); }


	/* kfac = 1; */

	ne = ne / kfac;
	u10 = ((double)kfac) * u10;
	u01 = ((double)kfac) * u01;

	selcomin = ((double)kfac) * scomin;





	/* Set the mutation constants. */

	phi = u01 / u10;											/* ratio of beneficial and deleterious mutation rates per site */

	mut0min = 1.0 - (((double)ellmin) * u01);				    /* fraction remaining in minor class 0 */
	mutdown0min = u10;											/* fraction of minor class 1 degrading to minor class 0 */
	mutemin = 1.0 - (((double)ellmin) * u10);				    /* fraction remaining in minor class ellmin */
	mutupemin = u01;											/* fraction of minor class (ellmin-1) moving to minor class ellmin */

	mut0maj = 1.0 - (((double)ellmaj) * u01);				    /* fraction remaining in major class 0 */
	mutdown0maj = u10;											/* fraction of major class 1 degrading to major class 0 */
	mutemaj = 1.0 - (((double)ellmaj) * u10);				    /* fraction remaining in major class ellmaj */
	mutupemaj = u01;											/* fraction of major class (ellmaj-1) moving to major class ellmaj */

	for (igmin = 1; igmin <= (ellmin - 1); ++igmin) {
		mutransupmin[igmin] = u01 * (((double)ellmin) - ((double)(igmin - 1))); 											/* gain of a + from any of the - in next lowest class */
		mutransdownmin[igmin] = u10 * ((double)(igmin + 1));																/* loss of a + from next highest class */
		mutrans0min[igmin] = 1.0 - (u10 * ((double)igmin)) - (u01 * (((double)ellmin) - ((double)igmin))); 	}		        /* stays unchanged */

	for (igmaj = 1; igmaj <= (ellmaj - 1); ++igmaj) {
		mutransupmaj[igmaj] = u01 * (((double)ellmaj) - ((double)(igmaj - 1))); 											/* gain of a + from next lowest class */
		mutransdownmaj[igmaj] = u10 * ((double)(igmaj + 1));																/* loss of a + from next highest class */
		mutrans0maj[igmaj] = 1.0 - (u10 * ((double)igmaj)) - (u01 * (((double)ellmaj) - ((double)igmaj))); 	}		        /* stays unchanged */




	/* Set the initial genotype frequencies, and zero out the fitnesses. */

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
			p0[igmin][igmaj] = 0.0; 
			relw[igmin][igmaj] = 0.0;
		} }

	p0[0][ellmaj] = 1.0;									/* frequency of + allele before burnin = 1.0 */





	/* Initiate the allele frequencies, counters, and test statistics. */

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
			pmutm[igmin][igmaj] = 0.0;						/* zero the various allele-frequency counters */
			psel[igmin][igmaj] = 0.0;
			pgtypexp[igmin][igmaj] = 0.0;
			pmutm1[igmin][igmaj] = 0.0;  } }

	for (ig = 0; ig <= 1001; ++ig) {
		ssdismin[ig] = 0.0;
		ssdismaj[ig] = 0.0; }
	
	igen = 0;
	tcount = 0;
	tint = 0.0;
	counter = 0;
	totw = 0.0;
	totwsq = 0.0;

	meanminhet = 0.0; 
	meanmajhet = 0.0;
	
	totigmin = 0.0;
	totigsqmin = 0.0;
	totwvarigmin = 0.0;

	totigmaj = 0.0;
	totigsqmaj = 0.0;
	totwvarigmaj = 0.0;
	
	newhighmin = ellmin;
	newlowmin = 0;
	newhighmaj = ellmaj;
	newlowmaj = 0;

	increment = ne / xinc;						/* increment in generations between statistic calculations (set as a fraction of Ne). */



	/* ******************************************************************************************************************************************* */


	/* Iterate the recursion equations to obtain the equilibrium expectations. */

	while (tint < ngens)  										/* iterate until the stopping criterion has been met. */
	{
		igen = igen + 1;

    start = (double)clock()/CLOCKS_PER_SEC;                                     


		/* Set the running upper and lower boundaries to the allelic count classes. */

		lowigmin = newlowmin - 1;								
		highigmin = newhighmin + 1;

		if (lowigmin < 0) { lowigmin = 0; }
		if (highigmin > ellmin) { highigmin = ellmin; }
		if (newlowmin < 0) { newlowmin = 0; }
		if (newhighmin > ellmin) { newhighmin = ellmin; }

		lowigmaj = newlowmaj - 1;
		highigmaj = newhighmaj + 1;

		if (lowigmaj < 0) { lowigmaj = 0; }
		if (highigmaj > ellmaj) { highigmaj = ellmaj; }
		if (newlowmaj < 0) { newlowmaj = 0; }
		if (newhighmaj > ellmaj) { newhighmaj = ellmaj; }


		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {		/* zero the frequencies of the classes */
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				psel[igmin][igmaj] = 0.0; } }



		/* Draw the random selection coefficients for both alleles, and set the fitnesses for the allelic classes for this generation. */
		/* The first (0,0) and final (ellmin,ellmaj) indices are the worst and best classes. */

		selcomaj1 = ((double)kfac) * gsl_ran_gaussian(rand_new2, sdsmaj);
		selcomaj2 = ((double)kfac) * gsl_ran_gaussian(rand_new2, sdsmaj);
		selcomaj = selcomaj1 - selcomaj2 + (((double)kfac) * scomaj);
		

		for (igmin = newlowmin; igmin <= newhighmin; ++igmin) {							/* genotypic fitnesses */

			sumfreqmin[igmin] = 0.0;

			for (igmaj = newlowmaj; igmaj <= newhighmaj; ++igmaj) {

				sumfreqmaj[igmaj] = 0.0;

				difmin = ((double)ellmin) - ((double)igmin);		/* number of mismatches at loci with minor effects */
				difmaj = ((double)ellmaj) - ((double)igmaj);		/* number of mismatches at loci with major effects */

				relw[igmin][igmaj] = exp(-selcomin * difmin) * exp(-selcomaj * difmaj);		}}



		/* Impose selection. */

		meanfit = 0.0;

		for (igmin = newlowmin; igmin <= newhighmin; ++igmin) {									/* calculate mean fitness */
			for (igmaj = newlowmaj; igmaj <= newhighmaj; ++igmaj) {
				meanfit = meanfit + (p0[igmin][igmaj] * relw[igmin][igmaj]); } }

		for (igmin = newlowmin; igmin <= newhighmin; ++igmin) {									/* weight the prior genotype frequencies by relative fitness */
			for (igmaj = newlowmaj; igmaj <= newhighmaj; ++igmaj) {
				psel[igmin][igmaj] = p0[igmin][igmaj] * relw[igmin][igmaj] / meanfit; } }



		/* Impose mutation on the genotypic classes. */

		sump = 0.0;

		if (lowigmin == 0) { lowmutmin = 1; }
		else { lowmutmin = lowigmin; }
		if (highigmin == ellmin) { highmutmin = ellmin - 1; }
		else { highmutmin = highigmin; }

		if (lowigmaj == 0) { lowmutmaj = 1; }
		else { lowmutmaj = lowigmaj; }
		if (highigmaj == ellmaj) { highmutmaj = ellmaj - 1; }
		else { highmutmaj = highigmaj; }


		for (igmin = lowmutmin; igmin <= highmutmin; ++igmin) {									/* first, do the minor-allele classes */
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[igmin][igmaj] = (mutransupmin[igmin] * psel[igmin - 1][igmaj]) + (mutransdownmin[igmin] * psel[igmin + 1][igmaj]) + (mutrans0min[igmin] * psel[igmin][igmaj]); }}

		if (lowigmin == 0) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[0][igmaj] = (mut0min * psel[0][igmaj]) + (mutdown0min * psel[1][igmaj]);	}}

		if (highigmin == ellmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pmutm1[ellmin][igmaj] = (mutemin * psel[ellmin][igmaj]) + (mutupemin * psel[ellmin - 1][igmaj]); }}




		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {									/* next, do the major-allele classes */
			for (igmaj = lowmutmaj; igmaj <= highmutmaj; ++igmaj) {
				pmutm[igmin][igmaj] = (mutransupmaj[igmaj] * pmutm1[igmin][igmaj - 1]) + (mutransdownmaj[igmaj] * pmutm1[igmin][igmaj + 1]) + (mutrans0maj[igmaj] * pmutm1[igmin][igmaj]);
				sump = sump + pmutm[igmin][igmaj]; 	}}

		if (lowigmaj == 0) {
			for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
				pmutm[igmin][0] = (mut0maj * pmutm1[igmin][0]) + (mutdown0maj * pmutm1[igmin][1]);
				sump = sump + pmutm[igmin][0]; }}

		if (highigmaj == ellmaj) {
			for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
				pmutm[igmin][ellmaj] = (mutemaj * pmutm1[igmin][ellmaj]) + (mutupemaj * pmutm1[igmin][ellmaj - 1]);
				sump = sump + pmutm[igmin][ellmaj]; }}



		/* Reset the next generation's expected genotype frequencies, and ensure that they sum to 1.0. */

		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {
				pgtypexp[igmin][igmaj] = pmutm[igmin][igmaj] / sump; 
			    p0[igmin][igmaj] = 0.0; 		}}



		
    stop = (double)clock()/CLOCKS_PER_SEC;                                      


		/* Sample the population for new genotype frequencies. */

		ntot = ne;
		sump = 0.0;
		
		newlowmin = ellmin;
		newhighmin = 0;
		newlowmaj = ellmaj;
		newhighmaj = 0;
		
		for (igmin = lowigmin; igmin <= highigmin; ++igmin) {
			for (igmaj = lowigmaj; igmaj <= highigmaj; ++igmaj) {

			if ((pgtypexp[igmin][igmaj] > 0.0) && (ntot > 0))  {
				pp = pgtypexp[igmin][igmaj] / (1.0 - sump);											/* this is the remaining frequency to sample */


				if (pp >= 1.0000000000000) {														
					draw = ntot;
					p0[igmin][igmaj] = ((double)draw) / ((double)ne); 	}

				else if ( (ntot > 100) && (((double) ntot) * pp < 5.0)  ) {
					draw = gsl_ran_poisson(rand_new, (((double) ntot) * pp) );
 					p0[igmin][igmaj] = ((double)draw) / ((double)ne); 	}

				else if ( (ntot > 100) && ((((double) ntot) * (1.0-pp)) < 5.0) ) {
					draw = gsl_ran_poisson(rand_new, (((double) ntot) * (1.0-pp)) );
					draw = ntot - draw;
 					p0[igmin][igmaj] = ((double)draw) / ((double)ne); }

				else { 
					draw = gsl_ran_binomial_tpe(rand_new, pp, ntot);
					p0[igmin][igmaj] = ((double)draw) / ((double)ne); 	}


				ntot = ntot - draw;
				sump = sump + pgtypexp[igmin][igmaj];

				if (p0[igmin][igmaj] > 0.0) {
					if (igmin < newlowmin) { newlowmin = igmin; }
					if (igmin > newhighmin) { newhighmin = igmin; }
					if (igmaj < newlowmaj) { newlowmaj = igmaj; }
					if (igmaj > newhighmaj) { newhighmaj = igmaj; } 	}
								
		}}}



    start = (double)clock()/CLOCKS_PER_SEC;                                      

		

		/* Calculate the summary statistics if the sampling interval is completed. */

		if (igen == increment) {
			igen = 0;
			counter = counter + 1;

			if (counter > burnin) {
				meanw = 0.0;
				meanigmin = 0.0;
				ssqigmin = 0.0;
				meanigmaj = 0.0;
				ssqigmaj = 0.0;


				for (igmin = 0; igmin <= ellmin; ++igmin) {
					for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
						
						sumfreqmin[igmin] = sumfreqmin[igmin] + p0[igmin][igmaj];
						sumfreqmaj[igmaj] = sumfreqmaj[igmaj] + p0[igmin][igmaj];

						meanw = meanw + (p0[igmin][igmaj] * (exp(-(selcomin / ((double)kfac))  * difmin) * exp(-(selcomaj / ((double)kfac)) * difmaj)));

						meanigmin = meanigmin + (p0[igmin][igmaj] * ((double)igmin));
						ssqigmin = ssqigmin + (p0[igmin][igmaj] * pow(((double)igmin), 2.0));

						meanigmaj = meanigmaj + (p0[igmin][igmaj] * ((double)igmaj));
						ssqigmaj = ssqigmaj + (p0[igmin][igmaj] * pow(((double)igmaj), 2.0)); } }
				
				wvarigmin = (ssqigmin - (meanigmin * meanigmin)) ;
				wvarigmaj = (ssqigmaj - (meanigmaj * meanigmaj)) ;

				tots = tots + selcomaj;
				totssq = totssq + pow(selcomaj, 2.0);

				totw = totw + meanw;
				totwsq = totwsq + pow(meanw, 2.0);
				
	            grandmeans = (tots / tint) / ((double)kfac);
	            vars = ((totssq / tint) - pow(grandmeans, 2.0)) / pow(((double)kfac), 2.0);

				totigmin = totigmin + meanigmin;
				totigsqmin = totigsqmin + pow(meanigmin, 2.0);
				totwvarigmin = totwvarigmin + wvarigmin;

				totigmaj = totigmaj + meanigmaj;
				totigsqmaj = totigsqmaj + pow(meanigmaj, 2.0);
				totwvarigmaj = totwvarigmaj + wvarigmaj;

				tint = tint + 1.0;
				tcount = tcount + 1;

				grandmeanigmin = totigmin / tint;
				varigmin = totwvarigmin / tint;

				grandmeanigmaj = totigmaj / tint;
				varigmaj = totwvarigmaj / tint;

				totmin0 = 0.0;											/* Calculates the heterozygosity at the minor and major site. */
				totmaj0 = 0.0;

				for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {				/* Is only relevant if there is one such site. */
					totmin0 = totmin0 + pgtypexp[0][igmaj]; 	}

				meanminhet = meanminhet + (2.0 * totmin0 * (1.0 - totmin0));

				for (igmin = 0; igmin <= ellmin; ++igmin) {				/* Is only relevant if there is one such site. */
					totmaj0 = totmaj0 + pgtypexp[igmin][0]; 	}

				meanmajhet = meanmajhet + (2.0 * totmaj0 * (1.0 - totmaj0));

				binmin1 = int(250.0 * (1.0 - totmin0)) + 1;			/* Separate the running allele frequencies into the steady-state distribution bins. */
				binmaj1 = int(250.0 * (1.0 - totmaj0)) + 1;

				ssdismin[binmin1] = ssdismin[binmin1] + 1.0;
				ssdismaj[binmaj1] = ssdismaj[binmaj1] + 1.0;
				
				
				if (tcount > tintprint) {
					printf("%9d, %9d, %9.1f, %9.0f, %7.6f, %10.6f, %7.6f, %6d, %4d, %4d, %6d, %4d, %4d, %10.5f, %7.5f, %12.6f, %9.6f, %14.9f, %12.9f\n", (ne*kfac), kfac, rlng[itera], tint, (totw / tint), (grandmeanigmin / ((double)ellmin)), (grandmeanigmaj / ((double)ellmaj)),
						newlowmin, newhighmin, (newhighmin - newlowmin), newlowmaj, newhighmaj, (newhighmaj - newlowmaj), pow(varigmin, 0.5), pow(varigmaj, 0.5), (meanminhet/tint), (meanmajhet/tint), grandmeans, pow(vars,0.5));
					
					tcount = 0; }

			}
		}								/* ends the summary statistic analysis for this point */

	}									/* ends the loop for generations. */



	/* Calculate the final statistics. */

	for (igmaj = 0; igmaj <= ellmaj; ++igmaj) {
		meanpfreqmaj[ig] = sumfreqmaj[ig] / tint; }

	for (igmin = 0; igmin <= ellmin; ++igmin) {
		meanpfreqmin[ig] = sumfreqmin[ig] / tint; }

	grandmeanw = totw / tint;

	grandmeanigmin = totigmin / tint;
	grandmeanigmaj = totigmaj / tint;

	varw = (totwsq / tint) - pow(grandmeanw, 2.0);
	varigmin = totwvarigmin / tint;
	varigmaj = totwvarigmaj / tint;

	grandmeans = tots / tint;
	vars = (totssq / tint) - pow(grandmeans, 2.0);

	meanperform = ((((double)ellmin) * scomin * grandmeanigmin / ((double)ellmin)) + (((double)ellmaj) * scomaj * grandmeanigmaj / ((double)ellmaj))) / gmaxrate;

	meanminhet = meanminhet / tint;
	meanmajhet = meanmajhet / tint;

	neuthet = 4.0 * ne * u01 / (1.0 + mutbeta + (ne * u10 * (1.0 + (6.0*mutbeta) + pow(mutbeta,2.0))));
	
	neestimate = meanminhet * (1 + mutbeta) / u10;
	
	neestimate = neestimate / (     (4.0*mutbeta) - (meanminhet * (1.0 + (6.0*mutbeta) + pow(mutbeta,2.0)))       );
	
	fracpine = ((double) neestimate) / ((double) ne);
	
	meanmin = grandmeanigmin / ((double)ellmin);
	meanmaj = grandmeanigmaj / ((double)ellmaj);

	

	fprintf(stream, " %11d, %11d, %11d,, %12.11f, %12.11f, %12.11f,, %12.11f, %12.11f,, %12.11f, %12.11f ,,  %9d ,, %17.0f, %11d, %17.0f ,, %12.11f, %12.11f ,, %12.11f, %12.11f ,, %12.11f, %12.11f ,, %12.11f,, %12.11f, %12.11f ,, %12.11f ,, %12.0f ,, %10.8f\n  ",
		(ne*kfac), ellmin, ellmaj,
		scomin, scomaj, sdsmaj, (grandmeans / ((double)kfac)), (pow(vars,0.5) / ((double)kfac)), 
		(u10 / ((double)kfac)), mutbeta, kfac, 
		ngens, burnin, (tint*((double)increment)), grandmeanw, varw, 
		meanmin, meanmaj, pow(varigmin, 0.5), pow(varigmaj, 0.5), 
		meanperform, meanminhet, meanmajhet, neuthet, (((double)neestimate)*kfac), fracpine);



/*
	for (ig = 1; ig <= 251; ++ig) {
		fprintf(stream, " %12.11f, %12.11f, %12.11f\n", ((((double) (ig-1)) / 250.0) + 0.002), (ssdismin[ig] / tint), (ssdismaj[ig] / tint) ); }
*/


	for (ig = 1; ig <= 125; ++ig) {
	    foldedmin = (ssdismin[ig] + ssdismin[251-ig]) / tint;
	    foldedmaj = (ssdismaj[ig] + ssdismaj[251-ig]) / tint;
	    fprintf(stream, " %12.11f, %12.11f, %12.11f\n", ((((double) (ig-1)) / 250.0) + 0.002), foldedmin, foldedmaj ); }







	printf("\n");

	fclose(stream);

}									/* End the set of iterations over all population sizes and mutation rates. */


exit(0);

}





