/**************************************************************************/
/* CodonW codon usage analysis package                                    */
/* Copyright (C) 2005            John F. Peden                            */
/* This program is free software; you can redistribute                    */
/* it and/or modify it under the terms of the GNU General Public License  */
/* as published by the Free Software Foundation; version 2 of the         */
/* License,                                                               */
/*                                                                        */
/* This program is distributed in the hope that it will be useful, but    */
/* WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           */
/* GNU General Public License for more details.                           */
/* You should have received a copy of the GNU General Public License along*/
/* with this program; if not, write to the Free Software Foundation, Inc.,*/
/* 675 Mass Ave, Cambridge, MA 02139, USA.                                */
/*                                                                        */
/*                                                                        */
/* The author can be contacted by email (jfp#hanson-codonw@yahoo.com Anti-*/
/* Spam please change the # in my email to an _)                          */
/*                                                                        */
/* For the latest version and information see                             */
/* http://codonw.sourceforge.net 					  */
/**************************************************************************/
/*                                                                        */
/* -----------------------        codon_us.C     ------------------------ */
/* This file contains most of the codon usage analysis subroutines        */
/* except for the COA analysis                                            */
/* Internal subroutines and functions                                     */
/* initilize_point    assigns genetic code dependent parameters to structs*/
/* initilize_coa      decides which cod/AA to include in a COA by default */
/* codon_usage_tot    Counts codon and amino acid usage                   */
/* ident_codon        Converts codon into a numerical value in range 1-64 */
/* codon_usage_out    Write out Codon Usage to file                       */
/* codon_error        Called after all codons read, checks data was OK    */
/* rscu_usage_out     Write out RSCU                                      */
/* raau_usage_out     Write out normalised amino acid usage               */
/* aa_usage_out       Write out amino acid usage                          */
/* how_synon          Calculates how synonymous each codon is             */
/* how_synon_aa       Calculates how synonymous each AA is                */
/* clean_up           Re-zeros various internal counters and arrays       */
/* base_sil_us_out    Write out base composition at silent sites          */
/* cai_out            Write out CAI usage                                 */
/* cbi_out            Write out codon bias index                          */
/* fop_out            Write out Frequency of Optimal codons               */
/* enc_out            Write out Effective Number of codons                */
/* gc_out             Writes various analyses of base usage               */
/* dot(,X)            prints a period every X times it is called          */
/* get_aa             converts a three base codon into a 1 or 3 letter AA */
/* cutab_out          Write a nice tabulation of the RSCU+CU+AA           */
/* dinuc_count        Count the dinucleotide usage                        */
/* dinuc_out          Write out dinucleotide usage                        */
/* coa_raw_out        Write out raw codon usage for use by COA analysis   */
/* sorted_by_axis1    Sorts genes according to their axis one position    */
/* gen_cusort_fop     COA specific, write out cu of genes by axis1 posit. */
/* highlow            Used sorted cu to calculate high_low chi sq. contin */
/* hydro_out          Write out Protein hydropathicity                    */
/* aromo_out          Write out Protein aromaticity                       */
/*                                                                        */
/*                                                                        */
/* External subroutines to codon_us.c                                     */
/* my_exit            Controls exit from CodonW closes any open files     */
/* tidy               reads the input data                                */
/* output             called from tidy to decide what to do with the data */
/* toutput            handles the reformatting and translation of seqs    */
/* output_long        if sequence is very long then process what we know  */
/*                    and write sequence to disk in fragments             */
/* open_file          Open files, checks for existing files               */
/* fileclose          Closes files and returns a NULL pointer or exits    */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include "codonW.h"
/********************* Initilize Pointers**********************************/
/* Various pointers to structures are assigned here dependent on the      */
/* genetic code chosen.                                                   */
/* paa                points to a struct containing Amino Acid names      */
/* pap                points to amino acid properties                     */
/* pcai               points to Adaptation values used to calc CAI        */
/* pfop               points to a struct describing optimal codons        */
/* pcbi               points to the same structure as pfop                */
/* pcu                points to data which has the translation of codons  */
/* ds                 is a struct describing how synonymous a codon is    */
/* da                 is a struct describing the size of each AA family   */
/* pcoa               points to a struct that describes columns to be     */
/*                    included/excluded from any COA analysis             */
/**************************************************************************/
int initilize_point(char code, char fop_species, char cai_species)
{
   paa = &amino_acids;
   pap = &amino_prop;
   pcai = &cai[cai_species];
   pfop = &fop[fop_species];
   pcbi = &fop[fop_species];
   pcu = &cu[code];
   ds = how_synon();                          
   da = how_synon_aa();                                     
   pcoa = &coa;

   printf ("\n");   
   if (pm->codonW)
     printf ("Genetic code is currently set to %s %s\n\n",pcu->des,pcu->typ);

   return 1;
}
/*******************How Synonymous is this codon  *************************/
/* This function discovers at run time how synonymous a codon is by check-*/
/* ing all other codons to see if they encode the same AA                 */
/* This saves a lot of time when new genetic codes are added              */ 
/**************************************************************************/
int *how_synon(void)
{
   static int      dds[65];
   int x,i;
   
   for (x = 0; x < 65; x++)
      dds[x] = 0;

   for (x = 1; x < 65; x++)
      for (i = 1; i < 65; i++)
     if (pcu->ca[x] == pcu->ca[i])
        dds[x]++;
   return dds;                             /* return a structure          */
}
/*******************How Synonymous is this AA     *************************/
/* This function discovers at run time how synonymous an amino acid is by */
/* checking all codons to see if they encode this same AA                 */
/* This saves a lot of time when new genetic codes are added              */ 
/**************************************************************************/
int *how_synon_aa(void)
{
   static int      dda[22];
   int x;
   
   for (x = 0; x < 22; x++)
      dda[x] = 0;

   for (x = 1; x < 65; x++)
      dda[pcu->ca[x]]++;
   return dda;                             /* return a structure          */
}
/********************* Initialise COA     *********************************/
/* Decides which codons or amino acids are to be included in a COA if only*/
/* the default choice is used. For an amino acid COA, only stops are excl */
/* but for a codon usage COA stop codons and non-synonymous codons are    */
/* excluded                                                               */
/* pcoa               points to a struct that describes columns to be     */
/*                    included/excluded from any COA analysis             */
/*                    structure contains AA and Codon information         */
/**************************************************************************/
int initilize_coa(char code)
{
   static char     initilized;
   static char     oldcode;
   int i;
    
   /* if called a second time return unless the genetic code has changed  */
   if (initilized && (oldcode == code)) return 1;

   for (i = 0; i < 22; i++)         /* for each amino acid                */
      if (i == 11 || i == 0)        /* stop codons have the value 11      */
     pcoa->amino[i] = FALSE;        /* see RECODING file for more details */
      else
     pcoa->amino[i] = TRUE;

   for (i = 0; i < 65; i++)         /* for each codon                     */
      if (*(ds + i) == 1 || pcu->ca[i] == 11 || i == 0) 
     pcoa->codons[i] = FALSE;
      else
     pcoa->codons[i] = TRUE;

   initilized = TRUE;               /* we have been called  ...           */
   return 1;
}
/****************** Codon Usage Counting      *****************************/
/* Counts the frequency of usage of each codon and amino acid this data   */
/* is used throughout CodonW                                              */
/* pcu->ca contains codon to amino acid translations for the current code */
/* and is assigned in initialise point                                    */
/**************************************************************************/
int codon_usage_tot(char *seq, long int how_many)
{
   char            codon[4];
   int             icode;
   int             i;
   
   for (i = 0; i < how_many - 2; i += 3) {
      strncpy(codon, (seq + i), 3);
      icode = ident_codon(codon);
      ncod[icode]++;                          /*increment the codon count */  
      naa[pcu->ca[icode]]++;                  /*increment the AA count    */ 
      codon_tot++;                            /*increment the codon total */
   }

   if (how_many % 3) {                        /*if last codon was partial */
      icode = 0;                              /*set icode to zero and     */
      ncod[0]++;                              /*increment untranslated    */ 
   }                                          /*codons                    */
   return icode;                              /*return the last codon     */
}

/****************** Ident codon               *****************************/
/* Converts each codon into a numerical array (codon) and converts this   */
/* array into a numerical value in the range 0-64, zero is reserved for   */
/* codons that contain at least one unrecognised base                     */
/*                                                                        */
/**************************************************************************/
int ident_codon(char *codon)
{
   int             icode = 0;
   int x;
   
   for (x = 0; x < 3; x++) {
      switch (codon[x]) {
      case 'T':
      case 't':
      case 'U':
      case 'u':
     codon[x] = (char) 1;
     continue;
      case 'C':
      case 'c':
     codon[x] = (char) 2;
     continue;
      case 'A':
      case 'a':
     codon[x] = (char) 3;
     continue;
      case 'G':
      case 'g':
     codon[x] = (char) 4;
     continue;
      case '\0':
     return 0;
      default:
     codon[x] = (char) 0;
     break;
      }
   }
   if (codon[0] * codon[1] * codon[2] != 0)
      icode = (codon[0] - 1) * 16 + codon[1]
     + (codon[2] - 1) * 4;
   else
      icode = 0;

   return icode;
}

/****************** Codon error               *****************************/
/* Does some basic error checking for the input data, it can be called    */
/* using different error levels, thus generating different types of       */
/* messages. Basically checks for start, stop codons and internal stop    */
/* codons. As well as non-translatable and partial codons                 */
/**************************************************************************/
long int codon_error(int x, int y, char *ttitle, char error_level)
{
   long int             ns = 0;                  /* number of stops       */
   long int        loc_cod_tot = 0;
   static int      error_lines = 0;
   int i;

   for (i = 1, ns = 0; i < 65; i++) {
     loc_cod_tot += ncod[i];
     if (pcu->ca[i] == 11)
       ns +=  ncod[i];                           /*count  stop codons     */
   }

   switch (error_level) {
     case 1:                                    /*internal stop codons    */
       ns = ns - valid_stops;           
       /* a stop was a valid_stop if it was the last codon of a sequence  */

       if ( ! valid_start && pm->warn ) {
           dot(0,10);   
           fprintf(pm->my_err, "\nWarning: Sequence %3li \"%-20.20s\" does "
               "not begin with a recognised start codon\n"
		     ,num_sequence,ttitle);
           error_lines++;
       }

       if (ns && pm->warn ) {
	        dot(0,10);  
	        if (pm->totals && pm->warn)
	         fprintf(pm->my_err,"\nWarning: some sequences had internal stop"
		     " codons (found %li such codons)\n", ns);
	        else
	         fprintf(pm->my_err, "\nWarning: Sequence %3li \"%-20.20s\" has "
             "%li internal stop codon(s)\n", num_sequence, ttitle, ns);
	        num_seq_int_stop++;
	        error_lines++;
       }
       break;
   case 2:                                
     dot(0,10);  
     if (ncod[0] == 1 && pcu->ca[x] != 11 && pm->warn){ /*  last codon was partial */
       fprintf(pm->my_err, 
	       "\nWarning: Sequence %3li \"%-20.20s\" last codon was partial\n"
	       ,num_sequence, ttitle);
       error_lines++;
     }else {
       if (ncod[0] && pm->warn){                        /* non translatable codons */
	    if (pm->totals)
	      fprintf(pm->my_err, 
		   "\nWarning: some sequences had non translatable"
		   " codons (found %li such codons)\n",  ncod[0]);
	    else
	      fprintf(pm->my_err, 
		   "\nWarning: sequence %3li \"%-20.20s\" has %li non translatable"
		   " codon(s)\n", num_sequence, ttitle, ncod[0]);
	    error_lines++; 
       }
       if (pcu->ca[x] != 11 && pm->warn ) {
	    if (!pm->totals){
	       fprintf(pm->my_err, 
		   "\nWarning: Sequence %3li \"%-20.20s\" is not terminated by"
		   " a stop codon\n", num_sequence, ttitle);
	       error_lines++;
            }     
       }
     }  
     break;
   case 3: 
                                   /* Nc error routines see codon_us      */
     dot(0,10);                    /* dot resetting internal counter      */
     if (x==3) x=4;                /* if x=3 there are no 3 or 4 fold AA  */ 
     fprintf(pm->my_err, 
	     "\nSequence %li \"%-20.20s\" contains ",num_sequence, ttitle);
     (y) ? fprintf(pm->my_err, "only %i ", (int) y) : 
       fprintf(pm->my_err, "no ");
     fprintf(pm->my_err, "amino acids with %i synonymous codons\n", x);
     fprintf(pm->my_err, "\t--Nc was not calculated \n");
     error_lines+=2;      
     break;
   case 4:                         /* run silent                          */
     break;
   default:
      my_exit(99,"Programme error in codon_error\n");
   }
   if ((((error_lines + 2) * 2) > pm->term_length) && pm->verbose 
       && pm->my_err == stderr ) {
     error_lines = 0;              /* count lines of errors               */
     dot(0,10);                     
     pause;
   }       
   return loc_cod_tot;             /* Number of codons counted            */
}

/****************** Codon Usage Out           *****************************/
/* Writes codon usage output to file. Note this subroutine is only called */
/* when machine readable output is selected, otherwise cutab_out is used  */
/**************************************************************************/
int codon_usage_out(FILE * fblkout, long int *nncod, int last_aa, 
                    int vvalid_stops, char *ttitle)
{
  long int ccodon_tot = 0;
  int x;
  char sp=pm->seperator;
  
  ccodon_tot = codon_error(last_aa, vvalid_stops, ""   , (char) 4); /*dummy*/

  /*example of output                                                     */ 
  /*0,0,0,0,3,2,2,0,0,0,0,0,0,3,0,0,                                      */
  /*0,0,0,4,3,4,1,7,0,0,0,0,3,1,3,1,Codons=100                              */       
  /*0,0,0,0,10,6,3,0,0,0,0,0,1,1,12,0,Universal Genetic code              */
  /*0,0,0,3,7,5,7,9,0,1,1,1,8,4,5,0,MLSPCOPER.PE1                         */

  for (x = 1; x < 65; x++) {
    
    fprintf(fblkout, "%i%c",nncod[x],sp);
    
    switch (x) {
    case 16:
      fprintf(fblkout, "\n");
      break;
    case 32:
	fprintf(fblkout, "Codons=%ld\n",ccodon_tot);
      break;
    case 48:
	fprintf(fblkout, "%.30s\n", pcu->des);
      break;
    case 64:
      fprintf(fblkout, "%.20s\n",ttitle);
      break;
    default:
      break;
    }
  }
  return 1;
}
/******************  RSCU  Usage out          *****************************/
/* Writes Relative synonymous codon usage output to file. Note this subrou*/
/* tine is only called if machine readable output is selected             */
/* If human readable format was selected then what the user really wanted */
/* was cutab so this is automatically selected in codons.c                */
/* RSCU values are genetic codon dependent                                */
/**************************************************************************/
int rscu_usage_out(FILE * fblkout, long *nncod, long *nnaa)
{  
 int x;
 char sp=pm->seperator;

 /* ds points to an array[64] of synonym values i.e. how synon its AA is  */

 for (x = 1; x < 65; x++) {
   if (nnaa[pcu->ca[x]] != 0)
     fprintf(fblkout, "%5.3f%c",
	     ( (float) nncod[x] / (float) nnaa[pcu->ca[x]])
	     *  ((float) *(ds + x)), sp );
   else
     fprintf(fblkout, "0.000%c",sp);

   if (x == 64)
     fprintf(fblkout, "%-20.20s", title);

   if (!(x % 16))
     fprintf(fblkout, "\n");
   }
   return 1;
}
/******************   RAAU output             *****************************/
/* Writes Relative amino acid usage output to file. Amino Acid usage is   */
/* normalised for gene length                                             */
/**************************************************************************/
int raau_usage_out(FILE * fblkout, long *nnaa)
{
   long int        aa_tot = 0;
   static char     first_line = TRUE;
   int i,x;
   char sp;

   if (pm->seq_format=='M')                     /*  if machine readable  */
      sp = pm->seperator;
   else
      sp = '\t';   

   if (first_line) {                            /* if true write a header*/
     if ( pm->seq_format=='M')
	 fprintf(fblkout, "%s", "Gene_name");
       else
	 fprintf(fblkout, "%-20.20s", "Gene name");

      for (i = 0; i < 22; i++)
	if ( pm->seq_format=='M')
	    fprintf(fblkout, "%c%s", sp,paa->aa3[i]);/* three letter AA names*/
	  else
	    fprintf(fblkout, "%c %-6.6s", sp,paa->aa3[i]);
      fprintf(fblkout, "\n");
      first_line = FALSE;
   }
   for (i = 1; i < 22; i++)
     if (i != 11)
       aa_tot += nnaa[i];                       /* total No. of AAs      */
   
   if ( pm->seq_format=='M')
     fprintf(fblkout, "%.30s", title);
   else
     fprintf(fblkout, "%-20.20s", title);       /* don't waste spaces    */
   
   for (x = 0; x < 22; x++)
     if (x == 11)
       fprintf(fblkout, "%c0.0000",sp);         /* report 0 for stops    */
     else if (aa_tot)
       if (  pm->seq_format=='M')
	   fprintf(fblkout, "%c%.4f",sp,
		   (double) nnaa[x] / (double) aa_tot);
	 else
	   fprintf(fblkout, "%c%7.4f",sp,
		   (double) nnaa[x] / (double) aa_tot);
     else                                       /*What no AminoAcids!!!! */
       if (  pm->seq_format=='M')
	 fprintf(fblkout, "%c%c",sp,sp);            
	 else
	   fprintf(fblkout, "%c ***** ",sp);        

   fprintf(fblkout, "\n",sp);
   return 1;
}
/******************   AA usage output         *****************************/
/* Writes amino acid usage output to file.                                */
/**************************************************************************/
int aa_usage_out(FILE * fblkout, long *nnaa)
{
  static char     first_line = TRUE;
  int i;
  char sp=pm->seperator;

  if (first_line) {
    (pm->seq_format=='M')?
      fprintf(fblkout, "%s", "Gene_name"):
      fprintf(fblkout, "%-20.20s ", "Gene name");
    
    for (i = 0; i < 22; i++)
      (pm->seq_format=='M')?
	fprintf(fblkout, "%c%s", sp,paa->aa3[i]):    /* 3 letter AA code     */
      fprintf(fblkout, "%-5.5s", paa->aa3[i]);
    fprintf(fblkout, "\n");
    first_line = FALSE;
  }
  (pm->seq_format=='M')?
    fprintf(fblkout, "%.20s", title):
    fprintf(fblkout, "%-20.20s ", title);
  
  for (i = 0; i < 22; i++){
    (pm->seq_format=='M')?
      fprintf(fblkout, "%c%li", sp,nnaa[i]):
      fprintf(fblkout, "%-5li",nnaa[i]);
  }

  fprintf(fblkout, "\n");
  return 1;
}
/******************  Base Silent output     *******************************/
/* Calculates and write the base composition at silent sites              */
/* normalised as a function of the possible usage at that silent site with*/
/* changing the amino acid composition of the protein. It is inspired by  */
/* GC3s but is much more complicated to calculate as not every AA has the */
/* option to use any base at the third position                           */
/* All synonymous AA can select between a G or C though                   */
/**************************************************************************/
void base_sil_us_out(FILE * foutput, long *nncod, long *nnaa)
{
   int             id,i,x,y,z;
   long            bases_s[4];     /* synonymous GCAT bases               */
                                      
   long            cb[4];          /* codons that could have been GCAT    */
   int             done[4];
   char sp=  (char) (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

   for (x = 0; x < 4; x++) {
     cb[x] = 0;
     bases_s[x] = 0;
   }                               /* blank the arrays                    */

   for (x = 1; x < 5; x++)
     for (y = 1; y < 5; y++)
       for (z = 1; z < 5; z++) {   /* look at all 64 codons               */
	 id = (x - 1) * 16 + y + (z - 1) * 4;

	 if (*(ds + id) == 1 || pcu->ca[id] == 11)
           continue;              /* if no synon skip to next       codon */
	 bases_s[z - 1] += nncod[id]; /* count No. codon ending in base X     */
       }     
			
   for (i = 1; i < 22; i++) {
     for (x = 0; x < 4; x++)      /* don't want to count bases in 6 fold  */
         done[x] = FALSE;         /* sites twice do we so we remember     */   

     if (i == 11 || *(da + i) == 1)
       continue;                  /* if stop codon skip, or AA not synony */

      for (x = 1; x < 5; x++)    /* else add aa to could have ended count */
     for (y = 1; y < 5; y++)
        for (z = 1; z < 5; z++) {
           id = (x - 1) * 16 + y + (z - 1) * 4; 
           /* assign codon values in range 1-64                           */
           if (pcu->ca[id] == i && done[z - 1] == FALSE) {
	   /* encode AA i which we know to be synon so add could_be_x ending*/
         /* by the Number of that amino acid                              */
	     cb[z - 1] += nnaa[i];    
	     done[z - 1] = TRUE;     /* don't look for any more or we might   */
                                 /* process leu+arg+ser twice             */
           }                       
        }
   }

   /* Now the easy bit ... just output the results to file                */      
   for (i = 0; i < 4; i++) {
      if (cb[i] > 0)
     fprintf(foutput, "%6.4f%c", (double) bases_s[i]/(double)cb[i], sp);
      else
     fprintf(foutput, "0.0000%c",sp);
   }
   return;
}
/******************  Clean up               *******************************/
/* Called after each sequence has been completely read from disk          */
/* It re-zeros all the main counters, but is not called when concatenating*/
/* sequences together                                                     */
/**************************************************************************/
int clean_up(long int *nncod, long int *nnaa)
{
   int x;
   int i;
   
   for (x = 0; x < 65; x++)
      nncod[x] = 0;
   for (x = 0; x < 23; x++)
      nnaa[x] = 0;
                                    /* dinucleotide count remembers the   */                                     
   dinuc_count(" ", 1);             /* last_base from the last fragment   */
                                    /* this causes the last base to be "" */
   for (x = 0; x < 3; x++)
      for (i = 0; i < 16; i++)
         din[x][i] = 0;

   dinuc_count(" ", 1);
   master_ic = tot = 
   non_std_char = AT_TOT = GC_TOT = AA_TOT = GAP_TOT = IUBC_TOT = 0; 
   long_seq = FALSE;
   valid_stops = valid_start = codon_tot = tot = fram = 0;                   
   return 1;
}
/*****************Codon Adaptation Index output   *************************/
/* Codon Adaptation Index (CAI) (Sharp and Li 1987). CAI is a measurement */
/* of the relative adaptiveness of the codon usage of a gene towards the  */
/* codon usage of highly expressed genes. The relative adaptiveness (w) of*/
/* each codon is the ratio of the usage of each codon, to that of the most*/
/* abundant codon for the same amino acid. The relative adaptiveness of   */
/* codons for albeit a limited choice of species, can be selected from the*/
/* Menu. The user can also input a personal choice of values. The CAI     */
/* index is defined as the geometric mean of these relative adaptiveness  */
/* values. Non-synonymous codons and termination codons (genetic code     */
/* dependent) are excluded. To aid computation, the CAI is calculated as  */
/* using a natural log summation, To prevent a codon having a relative    */
/* adaptiveness value of zero, which could result in a CAI of zero;       */
/* these codons have fitness of zero (<.0001) are adjusted to 0.01        */
/**************************************************************************/
int cai_out(FILE * foutput, long int *nncod)
{
   long int        totaa = 0;
   double          sigma;
   float           ftemp;
   int x;
   char sp=  (char) (pm->seq_format=='H')? 
       (char) '\t': 
       (char) pm->seperator;
   static char       cai_ttt = FALSE;
   static char       description[61];
   static char       reference[61];
  
   static CAI_STRUCT user_cai;


   if (!cai_ttt ) {                       /* have we been called already   */     
      user_cai.des = description;         /* assign an array to a pointer  */
      user_cai.ref = reference;           /* as above                      */
      
      if ( pm->caifile==NULL && pm->verbose==TRUE 
	   && pm->menu==TRUE && (pcai == cai )){
          /* this is false                                                 */
	  /* if personal caifile is on commandline or                      */
          /* in non-interactive mode or -silent option                     */
          /* or cai values are not the default values                      */
	  

	  printf("\nDo you wish to input a personal choice of CAI"
          " values (y/n) [n] ");
      gets(pm->junk);

      /* This allows a user defined choice of CAI values to be selected    */ 
      if ('Y' == (char) toupper( (int) pm->junk[0])) {
          /* tell the user a little about what we are looking for          */
          printf("\nInput file must contain 64 CAI values\n"
                 "ranging from 0.00 to 1.00\n"
                 "values must be separated by spaces\n");
         /* open the CAI adaptiveness values file                          */
           if (!(pm->caifile = open_file("file with CAI values"
                       ,"cai.coa", "r", 0))) my_exit(6,"cai_out");
      
      }
      }                                          /* matched if pm->caifile=*/
     if (pm->caifile){  
       rewind (pm->caifile);        /* unlikely unless fopfile = caifile   */
       x = 0;
       strcpy(user_cai.des,"User supplied CAI adaptation values ");
       strcpy(user_cai.ref,"No reference");
       user_cai.cai_val[x++] = (float) 0.0;

     while ((fscanf(pm->caifile, "%f ", &ftemp)) != EOF) {
                                    /* if any bad CAI values are read EXIT*/
         if (ftemp < 0 || ftemp > 1.0) {
           printf("\nError CAI %f value out of range\nEXITING",ftemp);
           my_exit(99,"cai_out");
        }                                        
        user_cai.cai_val[x++] = ftemp;                    /* assign value */
     }                                                    /* end of while */
     if (x != 65) {                 /*             wrong number of codons */
        fprintf(pm->my_err, "\nError in CAI file, found %i values"
            " expected 64 values EXITING\n", x - 1);
        my_exit(99,"cai_out");
     }
     pcai = &user_cai;              /* assigns pointer to user CAI values */
      }                             /*        matches if( pm->caifile...  */

    
     printf ("Using %s (%s) w values to calculate "
	        		      "CAI \n",pcai->des,pcai->ref);
     cai_ttt = TRUE;                /*stops this "if" from being entered  */

    }                              /* matches if (!cai_ttt )             */
   
   for (x = 1, sigma = 0; x < 65; x++) {
      if (pcu->ca[x] == 11 || *(ds + x) == 1) continue;
      if (pcai->cai_val[x] < 0.0001)/* if value is effectively zero       */
            pcai->cai_val[x] = (float) 0.01;               /* make it .01 */
      sigma += (double) *(nncod + x) * log((double) pcai->cai_val[x]);
      totaa += *(nncod + x);
   }

   if (totaa) {                     /* catch floating point overflow error*/
      sigma = sigma / (double) totaa;
      sigma = exp(sigma);
   } else
      sigma = 0;

   fprintf(foutput, "%5.3f%c", sigma,sp);
   return 1;
}
/*****************Codon Bias Index output        **************************/
/* Codon bias index is a measure of directional codon bias, it measures   */
/* the extent to which a gene uses a subset of optimal codons.            */
/* CBI = ( Nopt-Nran)/(Nopt-Nran) Where Nopt = number of optimal codons;  */
/* Ntot = number of synonymous codons; Nran = expected number of optimal  */
/* codons if codons were assigned randomly. CBI is similar to Fop as used */
/* by Ikemura, with Nran used as a scaling factor. In a gene with extreme */
/* codon bias, CBI will equal 1.0, in a gene with random codon usage CBI  */
/* will equal 0.0. Note that it is possible for Nopt to be less than Nran.*/
/* This results in a negative value for CBI.                              */
/* ( Bennetzen and Hall 1982 )                                            */
/**************************************************************************/
int cbi_out(FILE * foutput, long int *nncod, long int *nnaa )
{
   long int        tot_cod  = 0;
   long int        opt      = 0; 
   float           exp_cod  = (float) 0.0; 
   float           fcbi;
   int             c,x;
   char            str[2];
   char sp=  (pm->seq_format=='H')? 
       (char) '\t':
       (char) pm->seperator;


   static char       description[61];
   static char       reference[61];
   static char       first_call_cbi  = TRUE;
   static char       has_opt_info[22];
   static FOP_STRUCT user_cbi;

   if (first_call_cbi) {                 /* have we been called already   */

     user_cbi.des = description;         /* assign a pointer to array     */
     user_cbi.ref = reference;    
      
      if ( pm->cbifile == NULL && pm->verbose==TRUE 
	  && pm->menu==TRUE && ( pcbi == fop )){ 
          /* this is false                                                 */
	  /* if personal fopfile is on commandline or                      */
          /* in non-interactive mode or -silent option                     */
          /* or fop values are not the default values                      */

      printf("\nDo you wish to input a personal choice of CBI"
         " values (y/n) [n] ");

      gets(pm->junk);

      if ('Y' == (char) toupper( (int) pm->junk[0])) {

     printf("\nInput file must contain 64 CBI values\n"
        " 1= rare codon\n 2= common codon\n 3= optimal codon\n");

     if (!(pm->cbifile = open_file("file with CBI values"
                       ,"cbi.coa", "r", 0)))
        my_exit(6,"cai_out");
          }                         /* matches if Y==                     */
     }                              /* matches if pm->cbifile==NULL       */


     if ( pm->cbifile ){
       rewind (pm->cbifile);        /* fopfile can be the same as cbifile */
       strcpy(user_cbi.des,"User supplied choice");
       strcpy(user_cbi.ref,"No reference");    
       x = 0;
       user_cbi.fop_cod[x++] = 0;

       while ((c = fgetc(pm->cbifile)) != EOF && x <=66) {
       sprintf (str,"%c",c);	
	 if (isdigit(c) && atoi(str) >= 0 
	     && atoi(str) <= 3) {
           user_cbi.fop_cod[x++] = (char) atoi(str);
	   
	 }                          /*                             isdigit */
       }                            /*                        end of while */

     if (x != 65) {                /*              wrong number of codons */
        sprintf(pm->messages, "\nError in CBI file %i digits found,  "
            "expected 64 EXITING\n", x - 1);
        my_exit(99,pm->messages);
     }                        
       pcbi = (&user_cbi);
    }                              /*             matches if(pm->cbifile)  */

    
     printf ("Using %s (%s) \noptimal codons to calculate "
	        		      "CBI\n",pcbi->des,pcbi->ref);


				   /* initilise has_opt_info             */			      
     for (x = 1; x < 22; x++) has_opt_info[x]=0;
     
     for (x = 1; x < 65; x++)     {
        if (pcu->ca[x] == 11 || *(ds + x) == 1) 
		continue;			      			      
        if (pcbi->fop_cod[x] == 3 ) 
		has_opt_info[pcu->ca[x]]++;        
     }  



     first_call_cbi = FALSE;       /*      this won't be called again      */
   }                               /*          matches if (first_call_cbi) */


   for (x = 1; x < 65; x++) {
      if (! has_opt_info[pcu->ca[x]])      continue;
      switch ((int) pcbi->fop_cod[x]) {
      case 3:
        opt     += nncod[x];
        tot_cod += nncod[x];
        exp_cod += (float) nnaa[pcu->ca[x]]/ (float) da[pcu->ca[x]]; 
      break;
      case 2:
      case 1:
        tot_cod += *(nncod + x);
        break;
      default:
         sprintf(pm->messages, " Serious error in CBI information found"
          " an illegal CBI value of %f for codon %i"
          " permissible values are \n 1 for non-optimal"
          " codons\n 2 for common codons\n"
          " 3 for optimal codons\n" " EXITING ",
          pcbi->fop_cod[x], x);
	 
          my_exit(99,pm->messages);
          break;
      }                             /*                   end of switch     */
   }                                /*                   for (    )        */                     

   if( tot_cod - exp_cod)
     fcbi= (opt - exp_cod) / (tot_cod - exp_cod);     
   else  
     fcbi= (float) 0.0; 
    
   fprintf(foutput, "%5.3f%c", fcbi,sp);                /* CBI     QED     */

   return 1;
}

/****************** Frequency of OPtimal codons output  ********************/
/* Frequency of Optimal codons (Fop) (Ikemura 1981). This index, is ratio  */
/* of optimal codons to synonymous codons (genetic code dependent). Optimal*/
/* codons for several species are in-built and can be selected using Menu 3*/
/* By default, the optimal codons of E. coli are assumed. The user may also*/
/* enter a personal choice of optimal codons. If rare synonymous codons    */
/* have been identified, there is a choice of calculating the original Fop */
/* index or a modified index. Fop values for the original index are always */
/* between 0 (where no optimal codons are used) and 1 (where only optimal  */
/* codons are used). When calculating the modified Fop index, any negative */
/* values are adjusted to zero.                                            */
/***************************************************************************/
int fop_out(FILE * foutput, long int *nncod)
{
   long int        nonopt = 0;
   long int        std = 0;
   long int        opt = 0;
   float           ffop;
   int             c,x;
   char            nonopt_codons = FALSE;
    
   char            str[2];


   char sp=  (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

   static char     first_call = TRUE;
   static char     description[61];
   static char     reference[61];
   static char     asked_about_fop = FALSE;
   static char     factor_in_rare = FALSE;
   static char     has_opt_info[22];
   static FOP_STRUCT user_fop;

   if (first_call) {                /* have I been called previously      */
     user_fop.des = description;
     user_fop.ref = reference;
     if ( pm->fopfile == NULL && pm->verbose==TRUE 
	  && pm->menu == TRUE && (pfop == fop )) {
          /* this is false                                                 */
	  /* if personal fopfile is on commandline or                      */
          /* in non-interactive mode or -silent option                     */
          /* or fop values are not the default values                      */

         printf("\nDo you wish to input a personal choice of Fop"
	          " values (y/n) [n] ");
         gets(pm->junk);
         if ('Y' == (char) toupper( (int) pm->junk[0])) {
          printf("\nInput file must contain 64 Fop values\n"
                 " 1= rare codon\n 2= common codon\n 3= optimal codon\n");

          if (!(pm->fopfile = open_file("file with Fop values"
                       ,"fop.coa", "r", 0))) my_exit(6,"fop_out");

         }                           /*                         if 'Y' == */
      }                              /* if (pm->fopfile == NULL........ ) */
  
 
    if ( pm->fopfile ) {
      rewind (pm->fopfile);          /*    possible for fopfile = cbifile */
      strcpy(user_fop.des,"User supplied choice");
      strcpy(user_fop.ref,"No reference");
      x = 0;
      user_fop.fop_cod[x++] = 0;
      
      while ((c = fgetc(pm->fopfile)) != EOF && x <=66) {
        sprintf (str,"%c",c);
      
        if (isdigit(c) && atoi(str) >= 0 
            && atoi(str) <= 3) {
	        user_fop.fop_cod[x++] = (char) atoi(str);	
        }                           /*                       test isdigit */
     }                              /*                       end of while */

     if (x != 65) {                 /*             wrong number of codons */
        sprintf(pm->messages, "\nError in Fop file %i values found,  "
            "expected 64 EXITING\n", x - 1);
        my_exit(99,pm->messages);
     }
     pfop = &user_fop;              /*  assigns pointer to user fop values*/
    }
     

     printf ("Using %s (%s)\noptimal codons to calculate "
	        		      "Fop\n",pfop->des,pfop->ref);
	
	
				   /* initilise has_opt_info             */			      
     for (x = 1; x < 22; x++) has_opt_info[x]=0;
        
     for (x = 1; x < 65; x++)     {
        if (pcu->ca[x] == 11 || *(ds + x) == 1) 
		continue;			      			      
        if (pfop->fop_cod[x] == 3 ) 
		has_opt_info[pcu->ca[x]]++;
	
	if (pfop->fop_cod[x] == 1 ){
	   if (!asked_about_fop && pm->verbose) {
             printf("\nIn the set of optimal codons you have selected,\n"
        	  "non-optimal codons have been identified\nThey can be "
        	  "used in the calculation of a modified Fop, "
        	  "(Fop=(opt-rare)/total)\n else the original formulae "
        	  "will be used (Fop=opt/total)\n\n\t\tDo you wish "
        	  "calculate a modified fop (y/n) [n] ");
	     gets(pm->junk);
	     if ( 'Y' == (char) toupper( (int)pm->junk[0]))
	       factor_in_rare = TRUE;
	     asked_about_fop = TRUE;
           }
	   
	   if ( factor_in_rare == TRUE )
	            has_opt_info[pcu->ca[x]]++;
        }  
    }                                 /*    matches for (x=1           */
   first_call = FALSE;
   }                                  /*    matches if ( !first_call ) */
   
   
   
   for (x = 1; x < 65; x++) {
      if (!has_opt_info[pcu->ca[x]] ) 
       continue;
      
      switch ((int) pfop->fop_cod[x]) {
      case 3:
     opt += *(nncod + x);
     break;
      case 2:
     std += *(nncod + x);
     break;
      case 1:
     nonopt_codons = TRUE;
     nonopt += *(nncod + x);
     break;
      default:                      
     sprintf(pm->messages, " Serious error in fop information found"
         " an illegal fop value of %f for codon %l"
         " permissible values are \n 1 for non-optimal"
         " codons\n 2 for common codons\n"
         " 3 for optimal codons\n" " EXITING ",
         pfop->fop_cod[x], x);
	 printf ("opt %l, std %l, nonopt %l\n",opt,std,nonopt); 
     my_exit(99,pm->messages);
     break;
      }
   }
                                    /* only ask this once  ...            */


   if (factor_in_rare && (opt + nonopt + std) )
      ffop = (float) (opt - nonopt) / (float) (opt + nonopt + std);
   else if ((opt + nonopt + std))
      ffop = (float) opt / (float) (opt + nonopt + std);
   else   
      ffop=0.0;


   fprintf(foutput, "%5.3f%c", ffop,sp);

   return 1;
}

/***************  Effective Number of Codons output   *********************/
/* The effective number of codons (NC) (Wright 1990). This index is a     */
/* simple measure of overall codon bias and is analogous to the effective */
/* number of alleles measure used in population genetics. Knowledge of the*/
/* optimal codons or a reference set of highly expressed genes is not     */
/* needed when calculating this index. Initially the homozygosity for each*/
/* amino acid is estimated from the squared codon frequencies.            */
/**************************************************************************/
float enc_out(FILE * foutput, long int *nncod, long int *nnaa) {
   int             numaa[9];
   int             fold[9];
   int             error_t = FALSE;
   int             i,z,x;
   double          totb[9];
   double          averb = 0, bb = 0, k2 = 0, s2 = 0;
   float           enc_tot = 0.0F;
   char sp=  (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

/* don't assume that 6 is the largest possible amino acid family assume 9*/
   for (i = 0; i < 9; i++) {            
      fold[i] = 0;              /* initialise arrays to zero             */
      totb[i] = 0.0;
      numaa[i] = 0;
   }

   for (i = 1; i < 22; i++) {   /* for each amino acid                  */
      if (i == 11)
     continue;                  /* but not for stop codons              */

      if (*(nnaa + i) <= 1)     /* if this aa occurs once then skip     */
     bb = 0;
      else {
     for (x = 1, s2 = 0; x < 65; x++) { 
         /* Try all codons but we are only looking for those that encode*/
         /* amino amid i, saves having to hard wire in any assumptions  */
        if (pcu->ca[x] != i) continue;           /* skip is not i       */


        if (*(nncod + x) == 0)  /* if codons not used then              */
           k2 = 0.0;            /* k2 = 0                               */
        else
           k2 = pow(((double) *(nncod + x) / (double) *(nnaa + i)),
            (double) 2);

        s2 += k2;               /* sum of all k2's for aa i             */
     }
     bb = (((double) *(nnaa + i) * s2) - 1.0) /  /* homozygosity        */
        (double) (*(nnaa + i) - 1.0);
      }

      if (bb > 0.0000001) {
     totb[*(da + i)] += bb;         /* sum of all bb's for amino acids  */
                                    /* which have z alternative codons  */
     numaa[*(da + i)]++;            /* where z = *(da+i)                */
      }
                                    /* numaa is no of aa that were z    */
      fold[*(da + i)]++;            /* fold z=4 can have 9 in univ code */
   }                                /* but some aa may be absent from   */
                                    /* gene therefore numaa[z] may be 0 */
   enc_tot = (float) fold[1];

   for (z = 2, averb = 0, error_t = FALSE; z <= 8; z++) {   
                                   /* look at all values of z if there  */
      if (fold[z]) {               /* are amino acids that are z fold   */
     if (numaa[z] && totb[z] > 0)
        averb = totb[z] / numaa[z];
     else if (z==3 && numaa[2] && numaa[4] && fold[z]==1 )   
                                   /* special case                      */
        averb = (totb[2] / numaa[2] + totb[4] / numaa[4]) * 0.5;
     else {                        /* write error to stderr             */
        codon_error( z, numaa[z], title, 3 );  
        error_t = TRUE;            /* error catch for strange genes     */
        break;
        }
     enc_tot += (float) fold[z] / (float) averb;    
                                   /* the calculation                   */
      }
   }

   if (error_t)
      fprintf(foutput, "*****%c",sp);
   else if (enc_tot <= 61)
      fprintf(foutput, "%5.2f%c", enc_tot,sp);
   else
      fprintf(foutput, "61.00%c",sp);

   return enc_tot;
}

/*******************   G+C output          *******************************/
/* This function is a real work horse, initially it counts base composit */
/* ion in all frames, length of gene, num synonymous codons, number of   */
/* non synonymous codons. Then dependent on the value for which used in  */
/* switch statement. We return various analyses of this data             */
/* if which ==1 then the output is very detailed, base by base etc.      */
/* if which ==2 then the output is for GC content only                   */
/* if which ==3 then the output is for GC3s (GC at synonymous 3rd posit) */
/* if which ==4 then the output is for L_sym                             */
/* if which ==5 then the output is for L_aa                              */
/* The output from this subroutine is in a tabular format if human read- */
/* able output is selected, and in columns if machine readable. Also the */
/* number of values reported changes as it is assumed the user has access*/
/* to a spreadsheet type programme if they are requesting tabular output */
/*************************************************************************/
void gc_out(FILE * foutput, FILE * fblkout, int which){

   long int        id;
   long int        bases[5];        /* base that are synonymous GCAT     */
   long int        base_tot[5];
   long int        base_1[5];
   long int        base_2[5];
   long int        base_3[5];
   long int        tot_s = 0;
   long int        totalaa = 0;
   static char     header = FALSE;
   int x,y,z;
   char sp=  (pm->seq_format=='H')? 
       (char) '\t': 
       (char) pm->seperator;

   typedef double lf;

   for (x = 0; x < 5; x++) {
      bases[x] = 0;                 /* initialise array values to zero    */
      base_tot[x] = 0;
      base_1[x] = 0;
      base_2[x] = 0;
      base_3[x] = 0;
   }

   for (x = 1; x < 5; x++)
      for (y = 1; y < 5; y++)
     for (z = 1; z < 5; z++) {      /* look at all 64 codons              */
        id = (x - 1) * 16 + y + (z - 1) * 4;

        if (pcu->ca[id] == 11)
           continue;                /* skip if a stop codon               */
        base_tot[x] += ncod[id];    /* we have a codon xyz therefore the  */
        base_1[x] += ncod[id];      /* frequency of each position for base*/
        base_tot[y] += ncod[id];    /* x,y,z are equal to the number of   */
        base_2[y] += ncod[id];      /* xyz codons .... easy               */
        base_tot[z] += ncod[id];    /* will be fooled a little if there   */
        base_3[z] += ncod[id];      /* non translatable codons, but these */
                                    /* are ignored when the avg is calc   */
        totalaa += ncod[id];

        if (*(ds + id) == 1)
           continue;                /* if not synon  skip codon           */

        bases[z] += ncod[id];       /* count no of codons ending in Z     */
                     
        tot_s += ncod[id];          /* count tot no of silent codons      */
                      
     }


   if (!tot_s || !totalaa) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file   \n");
      return;
   }
   switch ((int) which) {
   case 1:                          /* exhaustive output for analysis     */
      if (pm->seq_format == 'M') {  /* machine readable format            */
     if (!header) {                 /* print a first line                 */
        fprintf(fblkout,
         "Gene_description%cLen_aa%cLen_sym%cGC%cGC3s%cGCn3s%cGC1%cGC2"
         "%cGC3%cT1%cT2%cT3%cC1%cC2%cC3%cA1%cA2%cA3%cG1%cG2%cG3\n"
		,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp,sp);
        header = TRUE;
     }
                                    /* now print the information          */
     fprintf(fblkout, "%-.20s%c", title,sp); 
     fprintf(fblkout, 
	     "%ld%c%ld%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c"
         "%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f%c"
         "%5.3f%c%5.3f%c%5.3f%c%5.3f%c%5.3f\n",
	     totalaa,sp,
	     tot_s,sp,
	     (lf) (base_tot[2] + base_tot[4]) / (lf) (totalaa * 3),sp,
	     (lf) (bases[2] + bases[4]) / (lf) tot_s,sp,
	     (lf) (base_tot[2] + base_tot[4] - bases[2] - bases[4])
	     / (lf) (totalaa * 3 - tot_s),sp,
	     (lf) (base_1[2] + base_1[4]) / (lf) (totalaa),sp,
	     (lf) (base_2[2] + base_2[4]) / (lf) (totalaa),sp,
	     (lf) (base_3[2] + base_3[4]) / (lf) (totalaa),sp,
	     (lf) base_1[1] / (lf) totalaa,sp, 
	     (lf) base_2[1] / (lf) totalaa,sp, 
	     (lf) base_3[1] / (lf) totalaa,sp,
	     (lf) base_1[2] / (lf) totalaa,sp, 
	     (lf) base_2[2] / (lf) totalaa,sp, 
	     (lf) base_3[2] / (lf) totalaa,sp,
	     (lf) base_1[3] / (lf) totalaa,sp, 
	     (lf) base_2[3] / (lf) totalaa,sp, 
	     (lf) base_3[3] / (lf) totalaa,sp,
	     (lf) base_1[4] / (lf) totalaa,sp, 
	     (lf) base_2[4] / (lf) totalaa,sp, 
	     (lf) base_3[4] / (lf) totalaa);
      } else {                      /* must be human formatted output then*/
     fprintf(fblkout,               /* tabulated output                   */ 
         "Gene Name: %-69.69s\nLength   : %-ld aa"
         " \tNon_synonymous/synonymous codons (%3ld/%5ld)\n"
         " GC=%5.3f\tGC3s=%5.3f\tGC_not_GC3s=%5.3f\n"
         "base\t1\t2\t3\ttotal\t\t1\t2\t3 \ttotal\n"
         "  T\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t"
         "W\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n"
         "  C\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t"
         "S\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n"
         "  A\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t"
         "R\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n"
         "  G\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t"
         "Y\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n\n",
         title,
         totalaa,
         totalaa - tot_s,
         tot_s,
         (lf) (base_tot[2] + base_tot[4]) / (lf) (totalaa * 3),
         (lf) (bases[2] + bases[4]) / (lf) tot_s,
         (lf) (base_tot[2] + base_tot[4] - bases[2] - bases[4])
         / (lf) (totalaa * 3 - tot_s),
         (lf) base_1[1] / (lf) totalaa, (lf) base_2[1] / (lf) totalaa, 
         (lf) base_3[1] / (lf) totalaa,
         (lf) base_tot[1] / (lf) (totalaa * 3),
         (lf) (base_1[1] + base_1[3]) / (lf) totalaa,
         (lf) (base_2[1] + base_2[3]) / (lf) totalaa,
         (lf) (base_3[1] + base_3[3]) / (lf) totalaa,
         (lf) (base_tot[1] + base_tot[3]) / (lf) (totalaa * 3),
         (lf) base_1[2] / (lf) totalaa, (lf) base_2[2] / (lf) totalaa, 
         (lf) base_3[2] / (lf) totalaa,
         (lf) base_tot[2] / (lf) (totalaa * 3),
         (lf) (base_1[2] + base_1[4]) / (lf) totalaa,
         (lf) (base_2[2] + base_2[4]) / (lf) totalaa,
         (lf) (base_3[2] + base_3[4]) / (lf) totalaa,
         (lf) (base_tot[2] + base_tot[4]) / (lf) (totalaa * 3),
         (lf) base_1[3] / (lf) totalaa, (lf) base_2[3] / (lf) totalaa, 
         (lf) base_3[3] / (lf) totalaa,
         (lf) base_tot[3] / (lf) (totalaa * 3),
         (lf) (base_1[3] + base_1[4]) / (lf) totalaa,
         (lf) (base_2[3] + base_2[4]) / (lf) totalaa,
         (lf) (base_3[3] + base_3[4]) / (lf) totalaa,
         (lf) (base_tot[3] + base_tot[4]) / (lf) (totalaa * 3),
         (lf) base_1[4] / (lf) totalaa, (lf) base_2[4] / (lf) totalaa, 
         (lf) base_3[4] / (lf) totalaa,
         (lf) base_tot[4] / (lf) (totalaa * 3),
         (lf) (base_1[1] + base_1[2]) / (lf) totalaa,
         (lf) (base_2[1] + base_2[2]) / (lf) totalaa,
         (lf) (base_3[1] + base_3[2]) / (lf) totalaa,
         (lf) (base_tot[1] + base_tot[2]) / (lf) (totalaa * 3));
                                    /* What hit me, did anyone see a bus  */
      }
      break;
   case 2:                          /* a bit more simple ... GC content   */
      fprintf(foutput, "%5.3f%c", (lf) ((base_tot[2] + base_tot[4]) / (lf) 
          (totalaa * 3)),sp);
      break;
   case 3:                          /* GC3s                               */
      fprintf(foutput, "%5.3f%c", (lf) (bases[2] + bases[4]) / 
          (lf) tot_s,sp);
      break;
   case 4:                          /* Number of synonymous codons        */
      fprintf(foutput, "%3li%c", tot_s,sp);
      break;
   case 5:                          /* Total length in translatable AA    */
      fprintf(foutput, "%3li%c", totalaa,sp);
      break;

#ifdef DEBUG
   default:
      fprintf(stderr, " Programming error in GC_out which (%i) is out of "
          "valid range\n"
          ,(int) which);
      my_exit(99, "gc out");
      break;
#endif
   }
   return;
}

/********************     DOT    ******************************************/
/*   Indicates the progress of a search                                   */
/**************************************************************************/

void dot(int y, long int period)
{
   static long int xx;
   static char     dott=0;

   if (!y) dott = 0;                /* re-zero the width counter          */

   if (++xx % period == 0){         /* every period calls print a .       */
             fprintf(stderr,".");
             dott++;
            }            
   if ( dott == 50) {               /* every 50 dots wrap the line        */
             fprintf(stderr,"\n");
             dott=0;
            }
   return;
}
/**********************  get_aa    *****************************************/
/* get_aa converts a numeric codon value (range 0-64 ) into a amino acid   */
/* and returns that amino acid number                                      */
/* pcu->ca converts the codon number into amino acid number                */
/* paa->aa1 converts  amino acid code into letters                         */
/***************************************************************************/

char *get_aa(int which, char *codon)
{
   char           *amino=NULL;

   if (strlen(codon) == 3) {
      if (which == 1)
     amino = paa->aa1[pcu->ca[ident_codon(codon)]];
      else
     amino = paa->aa3[pcu->ca[ident_codon(codon)]];
   } else {
      amino = amino;
      amino = paa->aa1[0];
   }
   return amino;
}
/**********************   cutab_out     ***********************************/
/* Generates a formatted table of codon, RSCU and amino acid usage        */
/* ds points to an array[64] of synonymous values                         */
/* it reveals how many synonyms there are for each aa                     */
/**************************************************************************/
int cutab_out(FILE * fblkout, long *nncod, long *nnaa)
{
   int             last_row[4];
   int             x;
   char            sp;

   if (pm->seq_format=='M')
      sp = pm->seperator;
   else
      sp = '\t';
                         
   for (x = 0; x < 4; x++)
      last_row[x] = 0;

   codon_tot = codon_error(1, 1, "", (char) 4); /*  dummy*/

   for (x = 1; x < 65; x++) {
     if (last_row[x % 4] != pcu->ca[x]){
       (pm->seq_format=='M')?
	 fprintf(fblkout, "%s%c%s%c", paa->aa3[pcu->ca[x]], sp, paa->cod[x], sp):
	 fprintf(fblkout, "%s %s"  , paa->aa3[pcu->ca[x]], paa->cod[x]);
     }
     else{
       (pm->seq_format=='M')?
	 fprintf(fblkout, "%c%s%c", sp, paa->cod[x], sp):
	 fprintf(fblkout, "    %s",   paa->cod[x]);
     }
/* Sample of output *******************************************************/      
/*Phe UUU    0 0.00 Ser UCU    1 0.24 Tyr UAU    1 0.11 Cys UGU    1 0.67 */ 
/*    UUC   22 2.00     UCC   10 2.40     UAC   17 1.89     UGC    2 1.33 */ 
/*Leu UUA    0 0.00     UCA    1 0.24 TER UAA    0 0.00 TER UGA    1 3.00 */ 
/*    UUG    1 0.12     UCG    6 1.44     UAG    0 0.00 Trp UGG    4 1.00 */
/**************************************************************************/
   (pm->seq_format=='M')?
      fprintf(fblkout, "%i%c%.2f%c",
          (int) nncod[x],
          sp, (nncod[x]) ?
          ((float) nncod[x] / (float) nnaa[pcu->ca[x]])
          * (float) (*(ds + x)):0,sp):           /* end of fprintf        */
      fprintf(fblkout, "%5i%5.2f ",
          (int) nncod[x],
          (nncod[x]) ?
          ((float) nncod[x] / (float) nnaa[pcu->ca[x]])
          * (float) (*(ds + x)):0);              /* end of fprintf        */

      last_row[x % 4] = pcu->ca[x];

      if (!(x % 4))
     fprintf(fblkout, "\n");
      if (!(x % 16))
     fprintf(fblkout, "\n");
   }
   fprintf(fblkout, "%li codons in %16.16s (used %22.22s)\n\n", 
       (long int) codon_tot, title, pcu->des);
   return 1;
}
/********************  Dinuc_count    *************************************/
/* Count the frequency of all 16 dinucleotides in all three possible      */
/* reading frames. This one of the few functions that does not use the    */
/* codon and amino acid usage arrays ncod and naa to measure the parameter*/
/* rather they use the raw sequence data                                  */
/**************************************************************************/
int dinuc_count(char *seq, long int ttot)
{
   static char     a = 0;
   int i;
   
   for (i = 0; i < ttot; i++) {
      last_base = a;
      switch (seq[i]) {
      case 't':
      case 'T':
      case 'u':
      case 'U':
     a = 1;
     break;
      case 'c':
      case 'C':
     a = 2;
     break;
      case 'a':
      case 'A':
     a = 3;
     break;
      case 'g':
      case 'G':
     a = 4;
     break;
      default:
     a = 0;
     break;
      }
      if (!a || !last_base)
     continue;                      /* true if either of the base is not  */
                                    /* a standard UTCG, or the current bas*/
                                    /* is the start of the sequence       */
      din[fram][((last_base - 1) * 4 + a) - 1]++;
      if (++fram == 3) fram = 0;    /* resets the frame to zero           */
   }
   return 1;
}
/***************** Dinuc_out           ************************************/
/* Outputs the frequency of dinucleotides, either in fout rows per seq    */
/* if the output is meant to be in a human readable form, each row repre- */
/* senting a reading frame. The fourth row is the total of the all the    */
/* reading frames. Machine readable format writes all the data into a     */
/* single row                                                             */
/**************************************************************************/
int dinuc_out(FILE * fblkout, char *ttitle)
{
   static char     called = FALSE;
   char            bases[5] = {'T', 'C', 'A', 'G'};
   char            sp = pm->seperator;
   long            dinuc_tot[4];
   int i,x,y;


   for ( x=0 ; x<4 ; x ++)  dinuc_tot[x]=0;
 

   for ( x=0 ; x<3 ; x++ )
       for ( i=0 ; i<16 ; i++ ){
         dinuc_tot[x]+=din[x][i];   /* count dinuc usage in each frame   */
         dinuc_tot[3]+=din[x][i];   /* and total dinuc usage,            */
       }

   if (pm->seq_format=='H' ) sp = ' ';

   if (!called) {                   /* write out the first row as a header*/
      called = TRUE;

      if  (pm->seq_format=='H' ) {
	   fprintf(fblkout,"%-13.13s%cframe%c","title", sp,sp);
	   for (x = 0; x < 4; x++)
	    for (i = 0; i < 4; i++) 
	     fprintf(fblkout,"%c%c%4.4c",bases[x],bases[i],sp);        
      }else{
	   fprintf(fblkout, "%s","title");
        for (y = 0; y < 4; y ++){
	     fprintf(fblkout, "%c%s",sp,"frame");
	     for (x = 0; x < 4; x++) 
	      for (i = 0; i < 4; i++) 
          fprintf(fblkout,"%c%c%c",sp, bases[x],bases[i]);          
          }
          }
      fprintf(fblkout, "\n");
       }                            /* matches if (!called)               */ 

/*Sample output   truncated  **********************************************/
/*title         frame TT    TC    TA    TG    CT    CC    CA    CG    AT  */
/*MLSPCOPER.PE1__ 1:2 0.024 0.041 0.016 0.008 0.049 0.041 0.033 0.098 ... */
/*MLSPCOPER.PE1__ 2:3 0.000 0.195 0.000 0.098 0.000 0.138 0.008 0.073 ... */
/*MLSPCOPER.PE1__ 3:1 0.008 0.016 0.000 0.033 0.033 0.107 0.172 0.262 ... */
/*MLSPCOPER.PE1__ all 0.011 0.084 0.005 0.046 0.027 0.095 0.071 0.144 ... */
/*MLSPCOPER.PE2__ 1:2 0.026 0.026 0.009 0.009 0.053 0.035 0.053 0.061 ... */  
/**************************************************************************/
   for (x = 0; x < 4; x++) {
   if ( pm->seq_format == 'H' || x == 0 )   
     fprintf(fblkout,  (pm->seq_format=='H') ?
	     "%-15.15s%c":"%-.15s%c", ttitle, sp);

     switch (x) {
     case 0:
         fprintf(fblkout, "1:2%c", sp);
       break;
     case 1:
         fprintf(fblkout, "2:3%c", sp);
       break;
     case 2:
         fprintf(fblkout, "3:1%c", sp);
       break;
     case 3:
         fprintf(fblkout, "all%c", sp);
       break;
     }

     if ( x == 3 ){ 
       for (i = 0; i < 16; i++)
          if ( dinuc_tot[x] )
	        fprintf(fblkout,"%5.3f%c",
              (float)(din[0][i]+din[1][i]+din[2][i])/
              (float)dinuc_tot[x], sp);     
          else
            fprintf(fblkout,"%5.3f%c",0.00, sp);              
     }
     else{
       for (i = 0; i < 16; i++)
         if ( dinuc_tot[x] )	
           fprintf(fblkout, "%5.3f%c", 
             (float) din[x][i]/(float)dinuc_tot[x], sp);
           else
               fprintf(fblkout,"%5.3f%c", 0.00, sp);              
     }

     if ( pm->seq_format == 'H' || x == 3) 
       fprintf(fblkout, "\n");
   }
   return 1;
}
/************* Coa_raw_out            *************************************/
/* Write out codon usage in a format compatible with the format required  */
/* by text2bin, i.e. part of the COA analysis suite of subroutines        */
/* rather than storing this data in memory, we first write raw codon usage*/
/* to disk, and then read it in as necessary, the file handle for this    */
/* data is passed via the fcoaout pointer. By default it writes to the    */
/* files coa_raw and coa1_raw                                             */
/**************************************************************************/
char coa_raw_out(FILE * fcoaout, long *nncod, long *nnaa, char *ttitle)
{

   static int      count = 0;
   int i;
      
   for (i = 0; i < (int) strlen(ttitle); i++)  /* don't take any chances  */
      if (isspace( (int) *(ttitle + i)))    *(ttitle + i) = '_';

   strncpy(pm->junk, ttitle, 20);              /* sequence name           */
   fprintf(fcoaout, "%i_%s ", ++count, pm->junk);

   switch (pm->coa) {
   case 'c':
   case 'r':                                  /* if rscu or codon usage   */ 
      for (i = 1; i < 65; i++)
     fprintf(fcoaout, "%i\t", (int) nncod[i]);
      fprintf(fcoaout, "\n");
      break;
   case 'a':                                  /* if amino acid usage      */
      for (i = 1; i < 22; i++)
     fprintf(fcoaout, "%i\t", (int) nnaa[i]);
      fprintf(fcoaout, "\n");
      break;
#ifdef DEBUG                                  /* Debugging code           */
   default:
      fprintf(pm->my_err, " Error in coa_out_raw\n");
#endif
   }
   return 1;
}
/**********  sorted_by_axis1    *******************************************/
/* COA specific routine, after the position of the genes on the first axis*/
/* has been computed the genes are sorted according to there ordination   */
/* this allows us to identify gene positioned at either end of the first  */
/* trend. Then the codon usage of these genes is used to determine the CU */
/* of these two groups. This information is used to identify optimal codon*/
/* calculated putative CAI adaptive values and for the Chi squared con-   */
/* tingency test, used to identify the optimal and non-optimal codons     */
/* The position of each gene on axis 1 is passed via the ax1 pointer      */
/* The integer rank of each sequence is stored in sortax1                 */
/* The number of genes is passed by the integer value lig                 */
/**************************************************************************/
void sorted_by_axis1(double *ax1, int *sortax1, int lig)
{
   double          min;
   int             nmin, *tagged;
   int             i,j;
   
   /* allocated an array such that we can record which genes have been    */
   /* processed already, and are in sortax1                               */
   if ((tagged = (int *) calloc(lig + 1, sizeof(int))) == NULL)
      my_exit(3, "sorted by axis 1");

   /* blank the array, shouldn't have to do this for ANSI C compilers     */
   for (i = 1; i <= lig; i++)
      tagged[i] = FALSE;

   /* for each gene                                                       */
   for (j = 1; j <= lig; j++) {
      i = 0;
      while (tagged[++i]);          /* find the first gene not in sortax1 */
      min = ax1[i];                 /* assign it value to min             */  
      nmin = i;                     /* assign it ordination to nmin       */

      for (i = 1; i <= lig; i++) {  /* for each gene                      */
       if (tagged[i]) continue;     /* gene is already in sortax1 .. next */
       if (ax1[i] < min) {          /* find the min value among the rest  */
        min = ax1[i];               /* assign it value to min             */ 
        nmin = i;                   /* assign it ordination to nmin       */
       }
      }
      sortax1[j] = nmin;            /* gene with lowest ax1 position is   */
      tagged[nmin] = TRUE;          /* assigned to sorax1 and tagged      */
   }
   free(tagged);
}
/***********  gen_cusort_fop                 ******************************/
/* COA specific routine, takes the sorted array of axis 1 positions from  */
/* sort_by_axis1 and passed via the sortax1 pointer. The array contains   */
/* the genes in order of occurrence in the original input file, but the   */
/* ranked order of each gene is recorded as the array value               */
/* This allows us to identify genes position at either end of the main    */
/* trend. Then the codon usage of these genes is used to write out a file */
/* with the genes in a axis1 position order                               */
/* the codon usage of the two groups at either end of the principle axis  */
/* are also counted. This information is then passed to highlow()         */
/* The position of each gene on axis 1 is passed via the ax1 pointer      */
/* The integer rank of each sequence is stored in sortax1                 */
/* The number of genes is passed by the interger value lig                */
/**************************************************************************/
void gen_cusort_fop(int *sortax1, int lig, FILE * fnam, FILE *ssummary)
{
   int             stops;
   long int       *low, *high;
   int             min, max, i ;
   float           v2;
   FILE           *fcusort = NULL;
   int            j;

   
   /* first open the original raw codon usage file                        */
   if ((fcusort = open_file("", "cusort.coa", "w", FALSE)) == NULL)
      my_exit(1, "gen_cusort_fop");                       

   /* calloc enough memory for the codon usage of the low group of genes  */
   if ((low = (long int *) calloc(65, sizeof(long int))) == NULL)
      my_exit(3, "low gen_cusort_fop");
   /* calloc enought memory for the codon usage of the high group of genes*/
   if ((high = (long int *) calloc(65, sizeof(long int))) == NULL)
      my_exit(3, "high gen_cusort_fop");

   /*pcoa->fop_gene is set in the advanced correspondence menu and is used*/
   /*to set the No of genes at either end of the principle axis that are  */
   /*to be used to create the low and high codon bias subsets of genes    */
   if (pcoa->fop_gene < 0) {        /* the number represent a percentage  */
      min = (int) ((float) lig * ((float) pcoa->fop_gene * -0.01));
      max = lig - (int) ((float) lig * ((float) pcoa->fop_gene * -0.01));
   } else {                        /*  the value is an absolute number    */
      min = pcoa->fop_gene;
      max = lig - pcoa->fop_gene;
   }

   if (min <= 0) {                 /* error catch in case % is too low    */
      min = 1;                     /* or fop_gene is set too high         */
      fprintf(pm->my_err, "Problems with the number genes used for"
          " fop adjusting to 1 gene\n");
   }
   if (max <= 0) {                 /* ditto                               */
      max = 1;
      fprintf(pm->my_err, "Problems with the number genes used for"
          " fop adjusting to one gene\n");
   }
   for (j = 1; j < 65; j++) {      /* initialise the blank array          */
      low[j] = 0;
      high[j] = 0;
   }

   /* write explanation about what we are doing to summary.coa            */ 
   fprintf(ssummary, "\ncusort.coa (not shown here) contains CU of "
       "genes sorted by their\n"
       "ordination on the principle axis or factor\n"
       "Genes used to calculate fop were 1 to %i and %i to %i\n"
       "these gene numbers REFER ONLY to the file cusort.coa\n"
       ,min, max + 1, pcoa->rows);

   for (i = 1; i <= lig; i++) {     /* foreach gene                       */
      rewind(fnam);                 /* go to start of codon_raw           */
      clean_up(ncod, naa);          /* blank the codon usage array        */    
      j = 1;
      while (j++ != sortax1[i])     /* find the rank of gene i            */ 
       fgets(pm->junk, BUFSIZ,fnam);/* by scanning for lines of CU in     */ 
      fscanf(fnam, "%s", pm->junk); /* now we know the name of seq i      */

      for (j = 1; j < 64; j++) {    /* now read in the cu of each codon   */
       fscanf(fnam, "%f", &v2);     /* assign it initially to v2          */ 
       ncod[j] = (long int) v2;     /* then place this value in ncod      */
     if (min >= i)                  /* remember the codon usage of the    */ 
        low[j] += (long int) v2;    /* two groups of genes at either end  */
     if (max < i)                   /* of the axis, containing min and    */
        high[j] += (long int) v2;   /* max genes                          */
      }

      fscanf(fnam, "%f\n", &v2);    /* now read the last codon in         */
      ncod[64] = (long int) v2;
      if (min >= i)
       low[64] += (long int) v2;
      if (max < i)
        high[64] += (long int) v2;  /* as above                           */

      /* we want to use codon_us_out to write out the sorted list of CU   */
      /* to cusort.coa. But if we have any internal stops etc, it will    */
      /* generate error messages, but we have already seen this messages  */
      /* on the first pass, so we fool it by saying all the stops are     */
      /* valid stops and not to complain again                            */
      for (j = 1, stops = 0; j < 65; j++)   
                 if (pcu->ca[j] == 11)
                        stops += (int) ncod[j];
      dot( 1 , 10 );  
      codon_usage_out(fcusort, ncod, 11, stops, pm->junk);
   }
   fileclose(&fcusort);              
   highlow(low, high, ssummary);        /* now we call highlow           */
                                        /* to use the sorted cu output   */
   free(low);                           /* release the memory to the OS  */
   free(high);
}

/************ highlow          ********************************************/
/* The codon usage of the two groups on either end of the axis is assigned*/
/* to low and high ... perhaps these would be better called left and right*/
/* as when they are passed to this function it is not know which group is */
/* lowly or highly biased. This is decided within highlow, by calculating */
/* the enc (a measure of bias) for each group and assigning the group with*/
/* the lowest enc as the higher biased genes. This works if the trend     */
/* represented by axis1 is truly selection for optimal translation        */
/* IT'S THE USERS RESPONSIBILITY TO ASSERTAIN IF THIS IS VALID            */
/* This information is used to identify optimal codons, as well as        */
/* calculate  putative CAI adaptive values and for the Chi squared con-   */
/* tingency test, used to identify the optimal and non-optimal codons     */
/**************************************************************************/

void highlow(long int *low, long int *high, FILE * ssummary)
{

   int            *last_row, icode, outer,i,j,x ;

   long int       *aa_low, *aa_high, *left, *right, *left_aa, *right_aa;
   long int       *highest_x;
   long int        right_tot = 0, left_tot = 0;

   float           enc_low, enc_high;
   float           a, b, c, d, e, f, g, h, total, hr, br, *x2;
   float           w;
   char           *flag, sp;

   FILE           *fcai=NULL,*fhilo = NULL, *ffop = NULL;
   FILE           *fcbi=NULL;

   /*calloc to the pointers the required storage                          */
   if ((fhilo = open_file("", "hilo.coa", "w", FALSE)) == NULL)
      my_exit(1, "hilo.coa");
   if ((ffop = open_file("", "fop.coa", "w", FALSE)) == NULL)
      my_exit(1, "fop.coa");
   if ((aa_low = (long int *) calloc(22, sizeof(long int))) == NULL)
      my_exit(3, "aa_low");
   if ((aa_high = (long int *) calloc(22, sizeof(long int))) == NULL)
      my_exit(3, "aa_high");
   if ((highest_x = (long int *) calloc(22, sizeof(long int))) == NULL)
      my_exit(3, "last_row");      
   if ((x2 = (float *) calloc(65, sizeof(float))) == NULL)
      my_exit(3, "x2");
   if ((flag = (char *) calloc(65, sizeof(char))) == NULL)
      my_exit(3, "flag");
   if ((last_row = (int *) calloc(65, sizeof(int))) == NULL)
      my_exit(3, "last_row");
   
  
   if (pm->seq_format=='M')
      sp = pm->seperator;
   else
      sp = '\t';

   /* initialize the various arrays                                       */
   for (x = 0; x < 4; x++) last_row[x] = 0;

   for (x = 0; x < 22; x++){
      highest_x[x]=0;
      aa_low   [x]=0;
      aa_high  [x]=0;
      }
   for (x = 0; x <65 ; x++) {
      x2      [x]= (float) 0.0;
      flag    [x]=0;   
      last_row[x]=0;
      }
      
      
   /*count the amino acid usage for the two datasets, initially we only   */
   /*have the codon usage of the two groups                               */
   for (i = 1; i < 65; i++) {
      aa_low[pcu->ca[i]] += low[i];
      aa_high[pcu->ca[i]] += high[i];
      flag[i] = ' ';                /*flag is used to identify opt codons */
   }

   enc_low = enc_out(fhilo, low, aa_low);         /*calc enc for each  of */
   enc_high = enc_out(fhilo, high, aa_high);      /*datasets              */
   fprintf(fhilo, "\n");

   fprintf(ssummary, "\nenc_left %f enc_right %f\n", enc_low, enc_high);

   for (i = 1; i < 65; i++) {
      if (*(ds + i) == 1 || pcu->ca[i] == 11)     /*skip stop and nonsynon*/
     continue;

      if (enc_low < enc_high) {                  /*decide which is more   */
        left = low;                              /*biased                 */
        right = high;                            /*left and right refer   */
        left_aa = aa_low;                        /*the columns of outputed*/
        right_aa = aa_high;                      /*hilow table            */
        a = (float) low[i];
        b = (float) high[i];
        g = (float) aa_low[pcu->ca[i]];
        h = (float) aa_high[pcu->ca[i]];
      } else {
        left = high;
        right = low;
        left_aa = aa_high;
        right_aa = aa_low;
        a = (float) high[i];
        b = (float) low[i];
        g = (float) aa_high[pcu->ca[i]];
        h = (float) aa_low[pcu->ca[i]];
     }
      /* calculate the chi squared contingency value                      */
      c = g - a;
      d = h - b;
      e = a + b;
      f = c + d;
      total = a + b + c + d;
      if (e * f * h * g)
     x2[i] = ((a * d - c * b) * (a * d - c * b)) * total / (e * f * g * h);
      else
     x2[i] = (float) -99.0;                   /*if 0 assign nonsense value*/

      if (g * h) {
     hr = a / g;
     br = b / h;
     if (hr > br && x2[i] > 6.635)            /* if significant at p<.99  */
        flag[i] = '*';
     else if (hr > br && x2[i] > 3.841)       /* if significant at p<0.05 */
        flag[i] = '@';
      }
   }
   fprintf(ssummary, "Chi squared contingency test of genes from both\n"
                     "extremes of axis 1\n");
/* this created the hi-low codon usage table                              */
/* Sample output truncated (***********************************************/
/*Asp   GAU   0.10 ( 10) 1.68 ( 53)   Gly   GGU   0.21 ( 12) 0.85 ( 11)   */   
/*      GAC*  1.90 (184) 0.32 ( 10)         GGC*  3.13 (176) 2.00 ( 26)   */   
/*Glu   GAA   0.00 (  0) 1.34 ( 55)         GGA   0.05 (  3) 0.69 (  9)   */  
/*      GAG*  2.00 (255) 0.66 ( 27)         GGG   0.60 ( 34) 0.46 (  6)   */   
/*                                                                        */
/*                                                                        */
/*        Number of codons in high bias dataset 2825                      */
/*        Number of codons in low  bias dataset 1194                      */
/*Note: high bias was assigned to the dataset with the lower average Nc   */
/*NO Chi could be calculated for UGU                                      */
/*Codon UUC (Phe) chi value was 70.175                                    */
/*Codon UCC (Ser) chi value was 48.030                                    */
/*Codon UAC (Tyr) chi value was 86.069                                    */
/**************************************************************************/ 

   for (outer = 1; outer <= 3; outer += 2) {
      for (x = 1; x < 5; x++) {
      for (j = 1; j < 5; j++) {
        icode = ((x - 1) * 16) + ((j - 1) * 4) + outer; 


        for (i = icode; i <= icode + 1; i++) {   /*loop twice             */
            /* if the previous entry in this column codes for the same AA */
            if (last_row[i % 2] != pcu->ca[i]) {
	          fprintf(fhilo, "%s%c%s%c%c", paa->aa3[pcu->ca[i]],
		              sp, paa->cod[i], flag[i], sp);
	          fprintf(ssummary, "%s%c%s%c%c", paa->aa3[pcu->ca[i]],
		              sp, paa->cod[i], flag[i], sp);
	        } else {                       
	           fprintf(fhilo, "%c%s%c%c", sp, paa->cod[i], flag[i], sp);
	           fprintf(ssummary, "   %c%s%c%c",sp,paa->cod[i],flag[i],sp);
	        }
            /* write out Codon usage, RSCU and significance for both data */
	       fprintf(fhilo, "%4.2f (%3i) %4.2f (%3i)%c",
		       (left[i]) ?
		       ((float) left[i] / (float) left_aa[pcu->ca[i]])
		       * (float) (*(ds + i))
		       : 0.0,
		       (int) left[i],
		       (right[i]) ?
		       ((float) right[i] / (float) right_aa[pcu->ca[i]])
		       * (float) (*(ds + i))
		       : 0.0,
		       (int) right[i],sp);               /*       end of fprintf  */
	      fprintf(ssummary, "%4.2f (%3i) %4.2f (%3i)%c",
		       (left[i]) ?
		       ((float) left[i] / (float) left_aa[pcu->ca[i]])
		       * (float) (*(ds + i))
		       : 0.0,
		       (int) left[i],
		       (right[i]) ?
		       ((float) right[i] / (float) right_aa[pcu->ca[i]])
		       * (float) (*(ds + i))
		       : 0.0,
		       (int) right[i],sp);               /*        end of fprintf */
          last_row[i % 2] = pcu->ca[i];          /* remember the last row */
        }
        fprintf(fhilo, "\n");
        fprintf(ssummary, "\n");
       }
       fprintf(ssummary, "\n");
       fprintf(fhilo, "\n");
      }
      fprintf(ssummary, "\n");
      fprintf(fhilo, "\n");
   }

   for (i = 1; i < 65; i++) {                    /* count both datasets   */
      right_tot += right[i];
      left_tot += left[i];
   }


   fprintf(fhilo, 
       "\tNumber of codons in high bias dataset %li\n", left_tot);
   fprintf(fhilo, 
       "\tNumber of codons in low  bias dataset %li\n", right_tot);
   fprintf(fhilo, 
       "Note: high bias was assigned to the dataset with the lower"
       " average Nc\n");

   fprintf(ssummary, 
       "\tNumber of codons in high bias dataset %li\n", left_tot);
   fprintf(ssummary, 
       "\tNumber of codons in low  bias dataset %li\n", right_tot);
   fprintf(ssummary, 
       "Note high bias was assigned to the genes with the lower"
       " overall Nc\n");

   /* now printout the Chi Squared values for each significant comparison */
   for (i = 1; i < 65; i++) {
      if (flag[i] == '*' || flag[i] == '@') {
     fprintf(fhilo, "Codon %s (%s) chi value was %.3f\n", paa->cod[i],
         paa->aa3[pcu->ca[i]], x2[i]);
     fprintf(ssummary, "Codon %s (%s) chi value was %.3f\n", paa->cod[i],
         paa->aa3[pcu->ca[i]], x2[i]);
      }
      if (x2[i] == -99)       /* there were no codons in one of the groups*/
     fprintf(fhilo, "NO Chi could be calculated for %s\n", paa->cod[i]);
   }
   fprintf(fhilo, "\n");
   fprintf(ssummary, "\n");

   /* now write out the optimal codons as PUTATIVELY identified by codonW */
   fprintf(ssummary, "These are the PUTATIVE optimal codons\n"
     "This is the format required for Menu 4 option 2 (Fop) "
     "and option 3 (CBI)\n"
     "This data is also duplicated in the files \"fop.coa\" "
     "and \"cbi.coa\"\n"
     "The format of these files is that required for input "
     "as a personal choice\n"
     "of optimal codons for these indexes\n");

   for (i = 1; i < 65; i++) {
      if( left[i] > highest_x[pcu->ca[i]])    /* used for calculating CAI */
                           highest_x[pcu->ca[i]]=left[i]; 
      
      if (*(ds + i) == 1 || pcu->ca[i] == 11) {
     fprintf(ffop, "2");
     fprintf(ssummary, "2");
      } else if (flag[i] == '*') {
     fprintf(ffop, "3");
     fprintf(ssummary, "3");
      } else if (((left[i]) ?
          ((float) left[i] / (float) left_aa[pcu->ca[i]])
          * (float) (*(ds + i))
          : 0.0) < 0.1) {                        /* if RSCU <0.1 its rare */
     fprintf(ffop, "1");
     fprintf(ssummary, "1");
      } else {
     fprintf(ffop, "2");
     fprintf(ssummary, "2");
      }

      if (!(i % 16)) {                           /* handle line wrapping  */ 
     fprintf(ffop, "\n");
     fprintf(ssummary, "\n");
      } else {
     fprintf(ffop, ",");
     fprintf(ssummary, ",");
      }
   }
   fileclose(&ffop);                              /*   close the Fop file  */
  
   if ((fcbi = open_file("", "cbi.coa", "w", FALSE)) == NULL)
      my_exit(1, "cbi.coa");                     /*    open cbi.coa       */
      
 for (i = 1; i < 65; i++) {                      /* write values 2 cbi.coa*/

  if (flag[i] == '*')                       /* Only report optimal codons */
     fprintf(fcbi, "3");
  else
     fprintf(fcbi, "2");                    /* ignore non optimal codons  */

  if (!(i % 16)) 
     fprintf(fcbi, "\n");
  else
     fprintf(fcbi, ",");
    
 }
   
  fileclose(&fcbi);   
   
   fprintf(ssummary, "\n\n");
    
   /* now calculate and write out CAI adaptiveness values                 */
   fprintf(ssummary, "These are PUTATIVE CAI adaptiveness values "
     "identified by this programme\n"
     "This data is also duplicated in the file \"cai.coa\"\n"
     "The format of this file is compatible with the format\n"
     "of the file used to input a personal selection of CAI values\n"
     "That is, the format required for Menu 4 option 1\n"
     "cai.coa\tinput file to be used for CAI calculations\n"
     "\n\nCod AA    Xi\tWi\t\tCod AA    Xi\tWi\n"); 
  
   
   if ((fcai = open_file("", "cai.coa", "w", FALSE)) == NULL)
      my_exit(1, "cai.coa"); 
  
   for (i = 1, x = TRUE ; i < 65 && x ; i++) {
    
    /* if a stop or a non-synonymous codon w = 1                          */
    if (*(ds + i) == 1 || pcu->ca[i] == 11) {  
                    fprintf(fcai, "1.0000000 \n");
                    fprintf(ssummary,"%s %s %6.1f %9.7f\t", 
                      paa->cod[i], 
                      paa->aa3[pcu->ca[i]],
                      (float) left[i], 1.0000000); 
    } else  if ( highest_x[pcu->ca[i]] ) {
      
      /* if a codon is absent then adjust its frequecy to 0.5             */
      if ( left[i] ) 
       w= (float) left[i]/ (float) highest_x[pcu->ca[i]];
      else
       w= (float) 0.5   / (float) highest_x[pcu->ca[i]];
      fprintf(fcai, "%9.7f \n", w);                    /* output CAI W    */
      fprintf(ssummary,"%s %s %6.1f %9.7f\t", 
             paa->cod[i], paa->aa3[pcu->ca[i]],
             (left[i]) ? (float) left[i]:0.5 , w); 
    /* either strange amino acid composition or data sets where too small */               
    } else {                            
      fprintf(pm->my_err, 
          "WARNING An attempt to calculate CAI relative "
          "adaptivnesss FAILED\n no %s amino acids found"
          " in the high bias dataset \n",paa->aa3[pcu->ca[i]]);    
      fprintf(ssummary, 
          "\nWARNING An attempt to calculate CAI relative adaptiveness "
          "FAILED\n no %s amino acids found in the high bias dataset \n",
          paa->aa3[pcu->ca[i]]);
      x=FALSE;
   }  
   if( !(i%2)) fprintf (ssummary  , "\n");
   } /* matches for (i = 1, x = TRUE ; i < 65 && x ; i++)                 */
     
   fileclose(&fcai);                              /* close files           */
   fileclose(&fhilo);
   free(aa_low);                                 /* free memory           */
   free(aa_high);
   free(highest_x);
   free(x2);
   free(flag);
   free(last_row);
   return;
}
/*********************  hydro_out        **********************************/
/* The general average hydropathicity or (GRAVY) score, for the hypothet- */
/* ical translated gene product. It is calculated as the arithmetic mean  */
/* of the sum of the hydropathic indices of each amino acid. This index   */
/* was used to quantify the major COA trends in the amino acid usage of   */
/* E. coli genes (Lobry, 1994).                                           */
/* Calculates and outputs total protein hydropathicity based on the Kyte  */
/* and Dolittle Index of hydropathicity (1982)                            */
/* nnaa               Array with frequency of amino acids                 */
/* paa                points to a struct containing Amino Acid values     */
/* pap->hydro         Pointer to hydropathicity values for each AA        */
/**************************************************************************/
int hydro_out(FILE * foutput, long int *nnaa)
{
   long int        a2_tot = 0;
   float           hydro = (float) 0.0;
   int i;
   char sp=  (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

   for (i = 1; i < 22; i++)
      if (i != 11) a2_tot += nnaa[i];

   if (!a2_tot) {           /* whow   .. no amino acids what happened     */
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file   \n", title);
      return 1;
   }
   
   for (i = 1; i < 22; i++)
      if (i != 11)
     hydro += ((float) nnaa[i] / (float) a2_tot) * (float) pap->hydro[i];

   fprintf(foutput, "%8.6f%c", hydro,sp );

   return 1;
}
/**************** Aromo_out ***********************************************/
/* Aromaticity score of protein. This is the frequency of aromatic amino  */
/* acids (Phe, Tyr, Trp) in the hypothetical translated gene product      */
/* nnaa               Array with frequency of amino acids                 */
/* paa                points to a struct containing Amino Acid values     */
/* pap->aromo         Pointer to aromaticity values for each AA           */
/**************************************************************************/
int aromo_out(FILE * foutput, long int *nnaa)
{
   long int        a1_tot = 0;
   float           aromo = (float) 0.0;
   int i;
   char sp=  (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

   for (i = 1; i < 22; i++)
      if (i != 11)
     a1_tot += nnaa[i];


   if (!a1_tot) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file   \n", title);
      return 1;
   }
   for (i = 1; i < 22; i++)
      if (i != 11)
     aromo += ((float) nnaa[i] / (float) a1_tot) * (float) pap->aromo[i];

   fprintf(foutput, "%8.6f%c", aromo,sp);
   return 1;
}


