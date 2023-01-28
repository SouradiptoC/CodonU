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
/* -----------------------        Codons.C       ------------------------ */
/* This file contains main() function and drives CodonW.                  */
/*                                                                        */
/* External subroutines and functions                                     */
/* clearscr           screen clearing Macro defined in CodonW.h           */
/* proc_comm_line     process command line arguments                      */
/* initilize_point    assigns genetic code dependent parameters to structs*/
/* initilize_coa      selects the default codons to exclude from the      */
/*                    Correspondence Analysis                             */
/* main_menu          The interactive menu system                         */
/* clean_up           Re-zeros various internal counters and arrays       */
/* open_file          Open files, checks for existing files               */
/* fileclose          Closes files and returns a NULL pointer or exits    */
/* textbin            Converts codon usage to binary data file            */
/* dot(,X)            prints a period every X times it is called          */
/* PrepAFC            Prepare for the COA                                 */
/* DiagoRC            This routine generates the COA                      */
/* colmout            write the output from COA to file                   */
/* rowout             save as above except records the gene information   */
/* inertialig         analyse row inertia and records the results to file */
/* inertiacol         analyse column inertia and record the results       */
/* suprow             add supplementary genes into COA                    */
/* get_aa             converts a three base codon into a 1 or 3 letter AA */
/* codon_error        Called after all codons read, checks data was OK    */
/* rscu_usage_out     Write out RSCU                                      */
/* codon_usage_out    Write out Codon Usage                               */
/* raau_usage_out     Write out normalised amino acid usage               */
/* dinuc_count        Count the dinucleotide usage                        */
/* dinuc_out          Write out dinucleotide usage                        */
/* aa_usage_out       Write out amino acid usage                          */
/* gc_out             Writes various analyses of base usage               */
/* cutab_out          Write a nice tabulation of the RSCU+CU+AA           */
/* base_sil_us_out    Write out base composition at silent sites          */
/* cai_out            Write out CAI usage                                 */
/* cbi_out            Write out codon bias index                          */
/* fop_out            Write out Frequency of Optimal codons               */
/* enc_out            Write out Effective Number of codons                */
/* hydro_out          Write out Protein hydropathicity                    */
/* aromo_out          Write out Protein aromaticity                       */
/* coa_raw_out        Write out raw codon usage for use by COA analysis   */
/*                                                                        */
/*                                                                        */
/* Internal subroutines to Codon.c                                        */
/* my_exit            Controls exit from CodonW closes any open files     */
/* tidy               reads the input data                                */
/* output             called from tidy to decide what to do with the data */
/* toutput            handles the reformatting and translation of seqs    */
/* output_long        if sequence is very long then process what we know  */
/*                    and write sequence to disk in fragments             */
/* file_close         Closes open files                                   */
/* c_help             Generates help informatio                           */
/* WasHelpCalled      Checks strings to see if help was requested         */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#define  ORIG_DEFS    
      /* used to decide whether declarations are external or not          */
      /* Master Header file                                               */
#include "codonW.h"     
#undef   ORIG_DEFS


#if defined(__MWERKS__)
#include <console.h> 
#endif     

/**************************   MAIN   **************************************/
/* The main function processes commandline arguments to decide whether    */
/* CodonW is running in an interactive mode, if so then the menu is called*/
/* CodonW also had the less documented feature of imitating other useful  */
/* codon usage and sequence manipulation program.    If the program is    */
/* called by a recognised name (see proc_comm_line for a list) such as    */
/* rscu then pm->codons is false and it only performs the required tasks  */
/* bypassing the menu system.                                             */
/* Main then calls tidy() to read in the data files, and count codon usage*/
/* depending on the requested output options toutput calls various subrou */
/* tines. If COA has been requested it also calls these subroutines and   */
/* recording useful information to summary.coa.                           */
/**************************************************************************/

int main(int argc, char *argv[])
{
  FILE       *finput = NULL, *foutput = NULL, *fblkout = NULL;
  FILE       *fcoaout = NULL;
  FILE       *fsummary= NULL;
  int        num_seq  = 0;

  num_sequence     = 0;
  num_seq_int_stop = 0;
  valid_stops      = 0;
  last_aa          = 0;
 
#if defined(__MWERKS__) /* Macintosh code-warrior */
  argc=ccommand(&argv); 
#endif

  pm = &Z_menu;
  pm->totals = FALSE;
  pm->my_err = stderr;


  initilize_point(pm->code, pm->f_type, pm->c_type);
  initilize_coa(pm->code);

  proc_comm_line(&argc, &argv);


  
  /******************** main loop ****************************/

  do {
    if (pm->codonW) {  
                                  /* If the program   chosen is codons    */
      printf(" Welcome to CodonW %.*s for Help type h\n\n",
              (int) strlen(Revision) - 11, Revision +10 ); 
                                  /* Now Run the main menu interface      */      
      if (pm->menu) main_menu(0); 
    } 

          /* if users select human readable output they want nice tables  */
    if (pm->bulk == 'C' && pm->seq_format == 'H')  pm->bulk = 'O';
    if (pm->bulk == 'S' && pm->seq_format == 'H')  pm->bulk = 'O';

     pm->analysis_run = TRUE;       /* codons has started an analysis this*/
                                    /* parameter is checked by my_exit    */

    if (pm->inputfile != NULL)      /* rewind various input files in case */
      rewind(pm->inputfile);        /* this is a second analysis run      */
    if (pm->fopfile != NULL)
      rewind(pm->fopfile);
    if (pm->cbifile != NULL)
      rewind(pm->cbifile);  
    if (pm->caifile != NULL)
      rewind(pm->caifile);

    /* num_sequence                 number of sequences read              */
    /* num_seq_int_stop             number with internal stop codons      */
    /* valid_stops                  No.terminated with a stop codon       */
    /* tot                          total number of codons read           */

    num_sequence = num_seq_int_stop = valid_stops = tot = 0;

    clean_up(ncod, naa);            /*re-zero count of amino and codons   */
    finput = pm->inputfile;
    foutput = pm->outputfile;
    fblkout = pm->tidyoutfile;

    fileclose(&pm->fcoa_out);
    if (pm->coa) 
      if ((pm->fcoa_out = open_file("", "coa_raw", "w", FALSE)) == NULL)
	my_exit(1, "coa_raw");          /*controlled exit from CodonW     */
    fcoaout = pm->fcoa_out;

/*  Tidy                                                                  */
/*  reads input data, returns the number of sequences read in             */
/*  num_sequence is global so I don't really have to assign it here       */    
    num_sequence = tidy(finput, foutput, fblkout, fcoaout);

    fprintf(pm->my_err,"\n\n\t\tNumber of sequences: %i\n",
        num_sequence);

/* num_seq_int_stop  value is calculated in codon_usage_out               */ 
    if (num_seq_int_stop > 0 && pm->warn ) {   
      if (pm->totals && (num_seq_int_stop >= valid_stops ))
	   fprintf(pm->my_err, "\tWARNING\t At least one sequence in your"
                    " input file has\ninternal stop codons (found %i"
                    " internal stops) \tWARNING\n",num_seq_int_stop);
      else
	   fprintf(pm->my_err, "\tWARNING\t %i sequences had internal "
                           "stop codons \tWARNING\n",num_seq_int_stop);
    }
/* don't wait for a pause if no_menu has been set                          */
    if ( pm->codonW && pm->menu ) pause; 

    if ( pm->coa && pm->totals)          /* idiots error catch             */ 
      my_exit(99,"A COA analysis of concatenated sequences is nonsensical\n"
      "I have completed any other requests but not the COA");

/* if COA has been requested then open summary.coa and start the analysis  */
    if (pm->coa) {
     if (fsummary == NULL)
      if ((fsummary = open_file("", "summary.coa", "w", FALSE)) == NULL)
       my_exit(1, "summary.coa");
/* set the number of genes in the analysis to the number read in by tidy   */
     pcoa->rows = num_sequence;
     fileclose(&fcoaout);
/* if COA has been selected then during the reading in phase raw codon usag*/
/* will have been written to the file coa_raw                              */
/* text bin converts this to binary data for the COA analysis program      */
      textbin("coa_raw", "cbrawin");
      printf("Generating correspondence analysis\n"); 
      dot(0,10);
      
         
      fprintf(fsummary, "\t\tSummary of Correspondence Analysis \n\n"
                        "The input file was %s it contained %i genes\n"        
                        "The number of axes generated was %i\n"
                        "A COA was requested of %s%s usage\n\n\n"
                        "Most of the output presented in this file "
                        "has also been written to separate files\n"
                        "genes.coa\tThe position of the genes on the "
                        "first %i axis\n" 
                        "%s.coa\tThe position of the %i %s on the %i "
                        "principle axes\n\n\n",
                        pm->curr_infilename, 
                        pcoa->rows,
                        ((pcoa->rows<pcoa->colm)?pcoa->rows:pcoa->colm)-1,
                        (pm->coa == 'r')        ?"relative synonymous ":"", 
                        (pm->coa == 'a')        ?"amino acid" : "codon",
                        pcoa->axis,
	                    (pm->coa == 'a')        ?"amino" : "codon",
	                    pcoa->colm,
	                    (pm->coa == 'a')        ?"amino acids":"codons", 
	                    pcoa->axis);
/* allocate memory for the rows and columns, scale both, and write out the*/
/* resulting matrix to the file cbrawin                                   */

      PrepAFC("cbrawin");

/* Now do the analysis, calculate the data inertia and all the vectors    */

      DiagoRC(fsummary);

/* colmout records the position of the columns on each of the factors/axes*/

      if (pm->coa == 'a')
          colmout("cbfcco", "amino.coa", paa, fsummary);
      else
          colmout("cbfcco", "codon.coa", paa, fsummary);

/* rowout records the position of the genes on each of the axis           */

      rowout("cbfcli", "genes.coa", "coa_raw", fsummary);

/* pcoa->level == e for exhaustive analysis of inertia                    */  
      if (pcoa->level == 'e') {     

       fprintf(fsummary, "\n\n\nYou requested detailed output from the COA"
           "\n\nThe absolute and relative inertia "
           "of each gene and %s (see also inertia.coa)\n",
           (pm->coa == 'a') ? "amino acids" : "codons");
/* inertialig must preceed inertiacol, records inertia of genes to file   */
/* it opens the raw codon usage file and loads the raw data to memory     */
       inertialig("inertia.coa", "coa_raw" ,fsummary);
/* uses the preloaded raw codon usage, to calculate inertia and other data*/
/* such as contribution of each column to each factor and to the extent   */
/* each column is explained by each factor and what the residual variation*/
/* is                                                                     */
       inertiacol("inertia.coa", fsummary);
      }
      
/* if pcoa->add_row is real string, then it will be the name of the file  */
/* containing additional sequence data, that will be excluded from the COA*/
/* but factored in, using the original COA vectors and then all other     */
/* calculation can proceed as with the original data                      */
      if (strlen(pcoa->add_row)) {
          if ((finput = open_file("", pcoa->add_row, "r", FALSE)) 
              == NULL) my_exit(6, "add_row");
          if ((foutput = tmpfile()) == NULL)
              my_exit(1, "temp file foutput");
          if ((fblkout = tmpfile()) == NULL)
              my_exit(1, "temp file fblkout");

    if ((fcoaout = open_file("", "coa1_raw", "w", FALSE)) == NULL)
      my_exit(1, "coa1_raw");

    clean_up(ncod, naa);
    num_sequence =num_seq_int_stop=valid_stops=tot = 0;
/* load the additional data file and process as normal                    */
/* but don't calculate any indices or write the data to the normal output */
/* files, rather write them to tmp files which will be deleted at end of  */
/* program execution                                                      */
    num_seq = tidy(finput, foutput, fblkout, fcoaout);

/* close the files now we are finished                                    */
    fileclose(&fcoaout);
    fileclose(&foutput);
    fileclose(&fblkout);
    fileclose(&finput);

/* covert to binary, use additional raw data file, note not coa_raw this  */
    textbin("coa1_raw", "cb1raw");
/* now call the routine suprow and add these additional genes, we will    */
/* process this data for inertia and append the gene and col. coordinates */
/* to the original gene.coa and codon.coa (or amino.coa)                  */
    suprow(num_seq, "cbfcvp", "cb1raw", "genes.coa", "coa1_raw", fsummary);
    
/* close these files now that we have finished with them and the COA      */

    fileclose(&foutput);
    fileclose(&fblkout);
    fileclose(&fcoaout);
      }
    }
   printf("\n");
  } while (pm->codonW && pm->menu ); /* OK now we loop back to main_menu  */
/* though only if we are in interactive mode and running as CodonW        */
  my_exit(0,"");                     /* last call to my_exit              */
  return 0;                          /* dummy return to keep pedantic but */
                                     /* brain dead compilers happy        */
}

/**********************  END of MAIN()   **********************************/


/**********************  Subroutines     **********************************/
/* Tidy                                                                   */
/*  reads input data from a sequence file containing fasta like formatted */
/*  sequence discards numbers, but keeps other characters                 */
/*  Each sequence must begin with title line must start with > or ;       */
/*  any following descriptive lines must begin with ; or >.Sequence start */
/*  is the first alphabetic character on the line following the headers   */
/*  There is no limit to sequence length or number of sequences but       */
/*  input lines should be less than 200 char in width                     */
/**************************************************************************/

int tidy(FILE * finput, FILE * foutput, FILE * fblkout, FILE * fcoaout)
{
  char            seq[MAX_GENE + LINE_LENGTH + 1];
  char            in[LINE_LENGTH + 1];
  int             first_line = TRUE, ic = 0;
  int             ii = 0;
  int             i,x;
  long            ic_orig = 0;
/* while still able to read data from the input file keep reading         */
  while ((fgets(in, LINE_LENGTH, finput) != NULL)) {

/* idiot error check to see if the file looks like fasta or PIR format    */ 
    if (!num_sequence && in[0] != ';' && in[0] != '>') {
      fprintf(stderr, "\n Error input file not in a recognised format \n"
                      " you must convert it into FASTA/Pearson format" 
                      " EXITING\n");
      my_exit(99, "input file not in a recognised format:tidy");
    }

    if (in[0] == ';' || in[0] == '>') { /* if true them this is a header   */
      if (first_line) {                 /* if true this is the first header*/

	first_line = FALSE;                 /* will only be reset when reread  */
                                        /* the next sequence               */
	if (num_sequence) {                 /* wait till we have read the first*/
                                        /* before writing to disk          */
/* now if we are concatenating sequence data we need will handle it thus   */
        if (pm->totals) {	    

/* first if translating or reformatting the input file flush the read      */
/* data to the disk                                                        */

         if (strchr("RNT",(int)pm->bulk)!=NULL) output_long(fblkout, seq);
         if (tot) {
   /* if something we have sequence read in, then we need to process this  */
   /* check whether the last codon of the sequence was was a stop          */
            last_aa = codon_usage_tot(seq, tot);
    	    if (pcu->ca[last_aa] == 11) valid_stops++;
	     }
/* rather re-setting everything to zero, we will just blank the array seq  */
         tot = 0;  
	  } else {     
/* else matches if tot; if sequences are not being concatenated we call    */
/* output to decide what to do with all the read data                      */
/* then we blank all the data from memory and start again                  */
        output(seq, foutput, fblkout, fcoaout);
	    clean_up(ncod, naa);      
	  }                 
	}                                    /* matches if(num_sequence)       */

/* If we get here we have read a header line, this then needs to be proc'ed*/
/* first the header is tested to see does it contain spaces the string is  */
/* converted from the first non space character to the title array         */

 for (ii = 1; isspace( (int) in[ii]) && ii < (int) strlen(in); ii++)
     ;
 strncpy(title, in + ii, 99);

/* Titles are cleaned up by removing newline characters and the delimiting */
/* character p,->seperater and also null terminating the title string      */

 for (i = 0; i < (int) strlen(title); i++) {
	  if (title[i] == '\n')	
          title[i] = '\0';          /* chops new line off                  */
	  else if (title[i] == pm->seperator )
	    title[i] = '_';             /* removes the separator if present    */
	  else if (i == (int) (strlen(title) - 1))
	    title[i] = '\0';            /* if we have reached end of title     */
	}

/* if we are reformatting the data, we print a friendly dot just in-case   */
if (strchr("RNT", (int)pm->bulk) ==NULL || pm->totals)
    dot((int) num_sequence, 5);       
/* we have now finished processing our first header line and are reading   */
/* our sequence data                                                       */
num_sequence++; 
      }                            /* matches if first line                */
      continue;                    /* read another line ie. jump to while()*/
    }                              /* if (in[0] == ';' || in[0] == '>')    */
    else{                          /* this must be a line containing seq   */
	first_line = TRUE;             /* so reset the first_line variable     */
      }

/* at this point we have read in the header lines and have been or about to*/
/* process the input data, now we test how much we have read into the array*/
/* seq, tot is equivalent to the last element in the array                 */
/* if tot is greater than or equal to MAX_GENE then the array is quite full*/
/* luckily we made the array seq to be MAX_GENE plus LINE_LENGTH +1        */

    if (tot >= MAX_GENE) {         /* sequence is larger than seq          */
      master_ic += MAX_GENE;       /* now remember how many bases we are   */
      ic_orig = tot;               /* going to write to disk               */
                                   /* and what size the array was to start */

      if (strchr("RNT", (int) pm->bulk) != NULL)
	     output_long(fblkout, seq);/* flush to disk and then continue      */
      else if (pm->bulk == 'D')
       	dinuc_count(seq, tot);     /* then we had better count the dinucs  */

/* Debugging code in-case we are asking for something that we can't handle */
#ifdef DEBUG
      else if (strchr("OCASLDBX", (int) pm->bulk) != NULL) ; /* dummy      */ 
      else if (pm->bulk)
	fprintf(stderr, "ERROR-22 %c pm->bulk undefined\n", pm->bulk);

      if (pm->cai || pm->fop || pm->cbi || pm->enc || pm->gc ||
	  pm->gc3s || pm->sil_base || pm->bulk ||
	  pm->coa);
      else
	fprintf(stderr, "Programming error");
#endif


/* Now count first MAX_GENE bases, luckily MAX_GENE is always a multiple of*/
/* 3, we count the bases and amino acids in codon_usage_tot                */

      last_aa = codon_usage_tot(seq, MAX_GENE);

/* now we move all unprocessed/written/counted bases to the front of seq   */

      for (i = MAX_GENE, x = 0; i < ic_orig; i++, x++)
       seq[x] = seq[i];            /* i is pointing near the end of array  */
      tot = x;                     /* x the front of the array             */
    }                              /* Matches if (tot >= MAX_GENE)         */

    ic = 0;                        /* first base of the input file         */
    while (in[ic] != '\0') {       /* scan input line till we see a Null   */
      if (isalpha((int)in[ic])) ;  /* do nothing if a alpha                */
      else if (pm->bulk == 'R' && in[ic] == '-'); /* do nothing            */
      else if (in[ic] == '*' || in[ic] == '.') ;  /* do nothing            */

      else {
    ic++;                          /* is not one above skip to next letter */
    continue;                     
      }                            /* while( in[ic] != '\0')               */


   in[ic] = (char)toupper((int)in[ic]);/*   converts2capitals              */
   if (strrchr("CG", (int) in[ic]) != NULL)
	GC_TOT++;                          /* is it a G or C                   */
   else if (strrchr("ATU", (int) in[ic]) != NULL)
	AT_TOT++;                          /* is it an A or T                  */
   else if ( in[ic] == '-' )
	GAP_TOT++;                         /* is it a gap character            */
   else
	non_std_char++;                    /* then it isn't a standard base    */

   if (strrchr("ABCDEFGHIKLMNPQRSTVWYZX"
		  ,(int) in[ic]) != NULL)
	AA_TOT++;                          /* it might be an amino acid        */
      if (strrchr("MRWSYKVHDBXN" , (int) in[ic]) != NULL)
	IUBC_TOT++;                        /* it might be a IUBC code          */

  seq[tot] = in[ic];                       /* move base into seq array         */
  seq[tot + 1] = '\0';                     /* make sure array is null term'ed  */


 /* now we test that the first codon is a valid start codon                    */

  if ( tot == 0 && master_ic == 0 ) {
      	
      in[1] = (char)toupper((int)in[1]);  /* Uppercase the first codon         */
      in[2] = (char)toupper((int)in[2]);  

       if ( in[1] == 'T' && (in[0] == 'A' || in[2] == 'G' ))
	  valid_start=TRUE;                /* Yeup it could be a start codon   */
	else
	  valid_start=FALSE;               /* Nope it doesn't seem to be one   */
  }
      ic++;                            /* total No. of sequence bases read */
      tot++;                           /* total currently stored in memory */
    }
  }                                    /* reached end of input file        */

/* Idiot error catch, this file is empty, at least it looks empty to codonW*/

  if ( !num_sequence ) my_exit(99,"The input file was empty");

/* better make sure to write anything left in seq to disk before returning */

  output(seq, foutput, fblkout, fcoaout);
  return (int) num_sequence;
}

/************************  TOUTPUT       **********************************/
/* toutput                                                                */
/*                                                                        */
/* This subroutine is very similar to output_long, basically it reformats */
/* or translates sequences less than MAX_GENE in length as a single read  */
/* It writes in reader format "ACG ATT ATC" i.e writes the sequence in    */
/* codons. Because it works with output_long it needs to know whether     */
/* the sequence being written to disk is a fragment or a complete gene    */
/**************************************************************************/
int toutput(FILE * fblkout, char *seq) {
  long int        ic = 0;
  int             space = 3;
  char            codon[4];
  int i,x;
  
  if (long_seq == FALSE) {         /* then this must be a complete genes  */
    switch (pm->bulk) {
    case 'T':                      /* tidy or fasta formatted header      */      
      fprintf(fblkout, ">%-20.20s%6li\n",
          title, (long int) tot + master_ic);
      break;
    case 'R':                      /* reader header .. don't ask          */ 
      fprintf(fblkout, ">%6li %-70.70s\n",
          (long int) tot + master_ic, title);
      break;
    case 'N':                      /* Conceptually translated DNA header  */
      fprintf(fblkout, ">%-20.20s%6li\n", 
           title, (long int) ((tot + master_ic) / 3));
      break;
    default:                       /* whoops                              */
      printf("\nProgramming error type A2 check code \n");
      my_exit(99, "toutput");
      break;
    }
  } else {                        
      
/* then long_seq must be true, this means we are about to finish writing a*/
/* sequence that has already been written in MAX_GENE chunks to disk)     */
/* when we wrote the original header line, we didn't know the size of the */
/* sequence, but now we do so we are going to update that bit of info     */
/* luckily remembered to record where the header line is in the file      */
/* its at fl_pos_start                                                    */

    fl_pos_curr = ftell(fblkout);   /* record where we are at present     */
    fseek(fblkout, fl_pos_start, 0);/* find the header line for this seq  */
    switch (pm->bulk) {
    case 'T':                       /* Now update the info                */
      fprintf(fblkout, ">%-20.20s%6li",     
          title, (long int) tot + master_ic);   
      break;
    case 'R':                  
      fprintf(fblkout, ">%6li %-70.70s",    
          (long int) tot + master_ic, title);   
      break;
    case 'N':
      fprintf(fblkout, ">%-20.20s%6li", title, 
              (long int) ((tot + master_ic) / 3));
      break;
    default:
      printf("\nProgramming error type A3 check code \n");
      my_exit(99, "output");
    }
    fseek(fblkout, fl_pos_curr, 0);/* now we move back to where we were   */
  }


  while (ic < tot) {               /* keep writing till the array is empty*/
    switch (pm->bulk) {           
    case 'T':
      fprintf(fblkout, "%c", seq[ic++]);
      reg++;
      break;
    case 'R':
      if (space == 3) {            /* Its reader format so print a space   */                                
    fprintf(fblkout, " ");         /* every third base                     */
    space = 0;                      
      } else {                     /* not the 3rd base yet so just print   */
    fprintf(fblkout, "%c", seq[ic++]);
    space++;
    reg++;
      }
      break;
    case 'N':
      for (i = (int) ic, x = 0; i < (int) ic + 3 && i < tot; i++, x++)  
          codon[x] = *(seq + i);   /* get the next three bases if there    */
      codon[x] = '\0';             /* null terminate the codon array       */
      ic += 3;                     /* remember that we have read 3 bases   */
      /* use the function get_aa to return the amino acid for the codon    */
      /* 1 = is for the one letter code of the codon                       */
      fprintf(fblkout, "%c", *get_aa(1, codon));   
      reg++;
      break;
    }
    if (!(reg % 61)) {             /* every 60 bases print a new line char */
      reg = 1;                  
      fprintf(fblkout, "\n");
    }
  }

  if (reg != 1) {                  /* reached the end of sequence so we    */
    fprintf(fblkout, "\n");        /* print a \n char unless we just did   */
    reg = 1;                       /* reset  number of bases printed       */
  }

/* Now that we have finished writing this sequence to disk lets have a     */
/* closer look at it, and do a few diagnostics about the bases used        */

  if (AT_TOT + GC_TOT > AA_TOT*0.5) {/* Assume its DNA then                */
    fprintf(pm->my_err, "%3li>\t%6li %-40.40s\tDNA\tGC%"                   
        " =%5.3f\n"                 /* with G+C content and length of gene */
        ,num_sequence
        ,(long int) tot + master_ic, title
        ,(float) GC_TOT / (GC_TOT + AT_TOT));

    if (non_std_char - IUBC_TOT && pm->warn )  /* any non IUBC characters */
      fprintf(pm->my_err, "\t\t WARNING %d non IUBC standard characters "
          "in sequence %i\n"
          ,non_std_char - IUBC_TOT
          ,num_sequence);
  } else {                         /* if not DNA then it must be a protein */
    fprintf(pm->my_err, "\t%3i>\t%6li %-40.40s\tPROTEIN\n"
        ,num_sequence
        ,(long int) tot + master_ic
        ,title);
    if ( (tot+master_ic)-AA_TOT && pm->warn)  /* non IUBC AA chars        */ 
      fprintf(pm->my_err, "\t\t WARNING %d non "                               
          "standard AA characters "
          "in sequence %i\n"
          ,non_std_char
          ,num_sequence);
  }
  return 1;                        /* return to calling function           */
}


/************************* output_long   **********************************/
/* called to write a block of a sequence that has exceeded the MAX_GENE   */
/* limit. If this is the first time it has been called for this sequence  */
/* (ie. long_seq is false) it write a dummy header line which is updated  */
/* by toutput when the last fragment of the sequence is written to disk   */
/**************************************************************************/

int output_long(FILE * fblkout, char *seq)
{
  long int        ic = 0;
  char            space = 3;
  char            codon[4];
  int i,x;
  
  if (long_seq == FALSE) {         
/* First call to output_long for seq. So record where the header line is  */
/* and then write the dummy header line.                                  */

      fl_pos_start = ftell(fblkout);
    if (pm->bulk == 'R')
      fprintf(fblkout, ">%6s %-72.72s\n", "      ", title); 
    else                   
      fprintf(fblkout, ">%-20.20s%9s\n", title, "    ");    
    long_seq = TRUE;               
  }
/* see toutput for explanation of the switch statement                    */                       
  while (ic < MAX_GENE && ic < tot) {     
    switch (pm->bulk) {
    case 'T':
      fprintf(fblkout, "%c", seq[ic++]);
      reg++;
      break;
    case 'R':
      if (space == 3) {
    fprintf(fblkout, " ");
    space = 0;
      } else {
    fprintf(fblkout, "%c", seq[ic++]);
    space++;
    reg++;
      }
      break;
    case 'N':
      for (i = (int) ic, x = 0; i < (int) ic + 3 && i < tot; i++, x++)
	codon[x] = *(seq + i);
      codon[x] = '\0';
      fprintf(fblkout, "%c", *get_aa(1, codon));
      ic += 3;
      reg++;
      break;
    default:
      printf("\nProgramming error type A1 check code \n");
      my_exit(99, "output_long");
    }
    if (!(reg % 61)) {
      reg = 1;
      fprintf(fblkout, "\n");
    }
  }
  return 1;                        /* return to tidy                      */
}


/*************************  output       **********************************/
/* Called from after subroutine tidy has read the sequence into memory    */
/* or  more accurately counted the codon and amino acid usage. This sub-  */
/* routine, via a switch checks which parameters and indices have been    */
/* requested and write these to file, it handles all output except for COA*/
/**************************************************************************/


void output(char *seq, FILE * foutput, FILE * fblkout, FILE * fcoaout)
{
  char sp;
 
  /* set the column delimiter to something shorter than pm->seperator     */
  sp = (char) (pm->seq_format=='H')? (char) '\t': (char) pm->seperator;

  if (tot) {                       /* still data in array seq..           */
    last_aa = codon_usage_tot(seq, tot);
    if (pcu->ca[last_aa] == 11)
      valid_stops++;                /* check the last codon was a stop    */
  }
  
  /* codon_error, if 4th parameter is 1, then checks for valid start and  */
  /* internal stop codon, if 4th parmater is 2, checks that the last codon*/
  /* is a stop or was partial, and for non-translatable codons            */
  codon_error(last_aa, valid_stops, title, (char) 1); 
  codon_error(last_aa, valid_stops, title, (char) 2);

  /* if we are concatenating sequences then change the title to avger_of  */ 
  if(pm->totals)                   
    (pm->seq_format=='M')?
      strcpy(title, "Average_of_genes"):
      strcpy(title, "Average of genes");
  
                
  if (strchr("RNT", (int) pm->bulk) != NULL) {
    /* better write the remaing sequence in seq to disk                   */
    toutput(fblkout, seq);                              
  } else if (strchr("OCASDLDBX", (int) pm->bulk) != NULL) {

/* These subroutines are self explanatory (see the top of this file)      */
/* are called such that only one can be called for each sequence read     */
/* all these calls are written to the bulk output file                    */

    switch ((int) pm->bulk) {
    case 'S':
      rscu_usage_out(fblkout, ncod, naa);
      break;
    case 'C':
      codon_usage_out(fblkout, ncod, last_aa, valid_stops, title);
      break;
    case 'L':
      raau_usage_out(fblkout, naa);
      break;
    case 'D':
      dinuc_count(seq, tot);
      dinuc_out(fblkout, title);
      break;
    case 'A':
      aa_usage_out(fblkout, naa);
      break;
    case 'B':
      gc_out(foutput, fblkout, 1);
      break;
    case 'O':
      cutab_out(fblkout, ncod, naa);
      break;
    case 'X':
            /* X is no bulk output written to file */
      break;
    default:
      fprintf(stderr, "ERROR-23 %s bulk undefined\n",  pm->prog);
      my_exit(99, "output");
      break;
    }
  } else if (pm->bulk) {            /* just a programming error catch     */
    fprintf(stderr, "ERROR-24 %s -prog undefined\n", pm->prog);
    my_exit(99, "output");
  }
  
  
  /* if an index has been requested then this is true                     */                    
  if (pm->sil_base || pm->cai || pm->fop   || pm->enc  || pm->gc3s ||
            pm->gc || pm->cbi || pm->L_sym || pm->L_aa || pm->coa  || 
            pm->hyd|| pm->aro) {
      /* if this is the first sequence then write a header line           */

    if (num_sequence == 1 || pm->totals) {

      fprintf(foutput, (pm->seq_format == 'H')?
	      "%-25.25s%c":"%-.25s%c"
	      ,"title",sp);
      if (pm->sil_base)
	fprintf(foutput, "%s%c%s%c%s%c%s%c", "T3s",sp,"C3s",sp,"A3s",sp,
"G3s",sp);
      if (pm->cai)
	fprintf(foutput, "%s%c", "CAI",sp);
      if (pm->cbi)
	fprintf(foutput, "%s%c", "CBI",sp);
      if (pm->fop)
	fprintf(foutput, "%s%c", "Fop",sp);
      if (pm->enc)
	fprintf(foutput, "%s%c", "Nc",sp);
      if (pm->gc3s)
	fprintf(foutput, "%s%c", "GC3s" ,sp);
      if (pm->gc)
	fprintf(foutput, "%s%c", "GC"   ,sp);
      if (pm->L_sym)
	fprintf(foutput, "%s%c", "L_sym",sp);
      if (pm->L_aa)
	fprintf(foutput, "%s%c", "L_aa" ,sp);
      if (pm->hyd)
	fprintf(foutput, "%s%c", "Gravy",sp);
      if (pm->aro)
	fprintf(foutput, "%s%c", "Aromo",sp);

      fprintf(foutput, "\n");
    }
    
    /* if output format is human readable print the fixed width sequence  */
    /* name, else print only the name of the sequence                     */
    fprintf(foutput, (pm->seq_format == 'H')?
	      "%-25.25s%c":"%-.25s%c"
	      ,title,sp);

    /*Need to use if statements as we allow more than one index to be calc*/
    /* per sequence read in                                               */
    if (pm->sil_base)
      base_sil_us_out(foutput, ncod, naa);
    if (pm->cai)
      cai_out(foutput, ncod);
    if (pm->cbi)
      cbi_out(foutput, ncod, naa);  
    if (pm->fop)
      fop_out(foutput, ncod);
    if (pm->enc)
      enc_out(foutput, ncod, naa);
    if (pm->gc3s)
      gc_out(foutput, fblkout, 3);
    if (pm->gc)
      gc_out(foutput, fblkout, 2);
    if (pm->L_sym)
      gc_out(foutput, fblkout, 4);
    if (pm->L_aa)
      gc_out(foutput, fblkout, 5);
    if (pm->hyd)
      hydro_out(foutput, naa);
    if (pm->aro)
      aromo_out(foutput, naa);
    if (pm->coa)
      coa_raw_out(fcoaout, ncod, naa, title);

    fprintf(foutput, "\n");

  }
  return;
}

/************************* my_exit       **********************************/
/* Called to clean up open files and generate an intelligent exit message */
/* Also warns if no analysis has been run, the user did not select R from */
/* the main menu. If COA was selected then it reminds the user to look    */
/* at the file summary.coa, and deletes any stray binary files            */
/**************************************************************************/

int my_exit(int error_num, char *message)
{

  fileclose(&pm->inputfile);
  
  /* if we are masuquarading as another program we assign both outputfile */
  /* and tidyout the same filehandle (we don't want to close this twice   */
  if ( pm->outputfile == pm->tidyoutfile ){
   fileclose(&pm->outputfile);
  }else{
   fileclose(&pm->outputfile);
   fileclose(&pm->tidyoutfile);
  }
  
  fileclose(&pm->cuout);
  fileclose(&pm->fopfile);
  fileclose(&pm->cbifile);
  fileclose(&pm->caifile);
  fileclose(&pm->logfile);
  fileclose(&pm->fcoa_in);
  fileclose(&pm->fcoa_out);

  if (pm->inputfile = fopen("cbrawin", "r")) {
    fclose(pm->inputfile);
    deletefile("cbrawin");
  }
  if (pm->inputfile = fopen("cbfcco", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcco");
  }
  if (pm->inputfile = fopen("cbfcli", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcli");
  }
  if (pm->inputfile = fopen("cbfcpc", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcpc");
  }
  if (pm->inputfile = fopen("cbfcpl", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcpl");
  }
  if (pm->inputfile = fopen("cbfcta", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcta");
  }
  if (pm->inputfile = fopen("cbfcvp", "r")) {
    fclose(pm->inputfile);
    deletefile("cbfcvp");
  }
  if (pm->inputfile = fopen("cb1rawin", "r")) {
    fclose(pm->inputfile);
    deletefile("cb1rawin");
  }
  if (error_num == 2 || error_num == 0 ) {
    if (pm->analysis_run) {
      fprintf(stderr, "Files used:\n");
      if (strlen(pm->curr_infilename))
	fprintf(pm->my_err, " Input file was\t %s \n", 
            	pm->curr_infilename);

      if (strlen(pm->curr_outfilename)){
	fprintf(pm->my_err, " Output file was\t %s %s",
		pm->curr_outfilename,
		(pm->codonW) ? " (codon usage indices, e.g. gc3s)\n":"\n");
      }

      if (strlen(pm->curr_tidyoutname)){
	fprintf(pm->my_err, " Output file was\t %s %s",
		pm->curr_tidyoutname,
		(pm->codonW) ? " (bulk output e.g. raw codon usage)\n":"\n");	
      }

      if (pm->coa)
	fprintf(pm->my_err, " For more information about the COrrespondence "
		"Analysis see summary.coa\n");
    } else if ( pm->codonW )          
      fprintf(stderr, " \n\n WARNING You are exiting before codonW has generated any results\n"
	      "  Select 'r' from the main menu to run\n");
  }

  if ( pm->codonW )  printf("\n CodonW has finished\n");

  switch ((int) error_num) {

  case 0:
    /* silent exit */
    exit(0);
    break;
  case 1:
    printf("failed to open file for output <%s>\n", message);
    exit(1);
    break;
  case 2:
    printf("user requested exit <%s>\n", message);
    exit(0);
    break;
  case 3:
    printf("failed to allocate memory <%s>\n", message);
    exit(1);
    break;
  case 4:
    printf("Write to disk failed ! <%s>\n", message);
    exit(1);
    break;
  case 5:
    printf("Read from disk failed! <%s>\n", message );
    exit(1);
    break;
  case 6:
    printf("failed to open file for reading <%s>\n", message);
    exit(1);
    break;
  case 7: 
    printf("failed to close file <%s>\n", message);
    exit(1);
  case 99:
    printf(" Controlled exit <%s>\n",message);
    exit(0);
    break;
  default:
    printf("for unknown reason\n");
    exit(1);
    break;
  }
  return 0;
}

/************************** file_close   **********************************/
/* Fileclose function checks whether the filepointer is open, if so it    */
/* attempts to close the open file handle and assigns a null pointer      */
/* to that  handle                                                        */
/**************************************************************************/

int  fileclose(FILE ** file_pointer)
{
   if (*file_pointer != NULL ) {
     if (fclose(*file_pointer) == EOF ) {
       fprintf(stderr,"Failed to close file %i \n",errno);
       perror ("Unexpected condition in fileclose");
       exit(7);
     }    
     *file_pointer = NULL;           /* make sure file_pointer is null*/
   }
  return 1;
}

/************************** Chelp    **************************************/
/* Chelp scans opens the help file and returns text associated with that  */
/* help keyword. Help keywords are surrounded by hashs, starting in the   */
/* first column of the ASCII help file and are terminated by //           */
/**************************************************************************/

int chelp ( char *help_keyword )
{
 char helplib [MAX_FILENAME_LEN]="";
 char *p=NULL, inhelp=FALSE;
 char QueryString[120];            /* limit for help phrase is 120 chars  */
 char HelpMessage[121];
 int  line_counter=2;              /* assume 2 blank lines to start with  */
 FILE *hfp=NULL;
/* Inital steps is to locate help file                                    */
/* First check if CODONW_H has been set as an environment variable        */
/* If not then assume that the help file is in the current directory      */       
   
 p=getenv( "CODONW_H" );
 if ( p != NULL ) 
     strcpy ( helplib , p );   
 else {
     strcpy ( helplib , "codonW.hlp");
  }

 hfp=open_file("",helplib, "r", FALSE);

/* if we can't open the help file then explain what we where trying to do  */

 if ( hfp == NULL ) {
           fprintf ( stderr , 
                  "Could not open help file codonw.hlp\n" 
                   "Expected to find this file in %s\n" 
                   "This can be overridden by setting the" 
                   "environmental variable\n"
                   "CODONW_H to the help file location\n", 
                   helplib);
                   pause;           /* make sure they Ack. the error mesg  */
                   return 0;        /* abort                               */
     }   
/* Now that we have opened the help file, assemble the help keyword string */

strcpy (QueryString , "#");
strcat (QueryString , help_keyword );
strcat (QueryString , "#");
fprintf(stderr,"\n\n");             

/* now scan the help file looking for this keyword                         */ 

while ( fgets ( HelpMessage, 120, hfp ) )   {

 if ( strstr (HelpMessage,QueryString) != NULL ) 
     inhelp=TRUE;                                          /* we found it  */

 else if  ( inhelp && strstr ( HelpMessage  , "//") )  {   /* found the end*/
           fileclose(&hfp );                            
	       if ( line_counter )pause;
           return 1;                                     
     }
 /* if inhelp is true we have found the help keyword but not reached EOF   */
 else if  ( inhelp ) {       
     if ( strchr(HelpMessage,'\n') )
        fprintf ( stderr, "%s",HelpMessage );    
                                           /*stderr,it must be interactive */
     else
        fprintf ( stderr, "%s\n",HelpMessage );
                                           /*make sure there are line feeds*/ 


  /* count how many lines I have printed to the terminal and compare it    */
  /* with the length of the terminal screen as defined by pm->term_length  */

  if (line_counter++ >= pm->term_length-3 && line_counter ) {
    line_counter=0;
    pause;
    fprintf(stderr, "%s",HelpMessage);
   }
 }   
}   

/* Error catches for problems with help file                              */
if ( HelpMessage == NULL && inhelp == FALSE ){
   fprintf ( stderr ," Error in help file, %s not found ", QueryString);
   pause;
  }
else {
    fprintf (stderr , "Premature end of help file ...  \n");
    pause;
 }   
return 0;                           /* failed for some reason             */
}


/******************** WasHelpCalled     ***********************************/
/* Checks the string input to see if the user asked for help              */
/**************************************************************************/

char WasHelpCalled ( char * input ) {
    char ans = FALSE;

    if ( strlen ( input) == 1 && (char)toupper((int)input[0]) == 'H') 
        ans = TRUE;
    else if ( !strcmp ( input , "help") ) 
        ans = TRUE;
    else if ( !strcmp ( input , "HELP") ) 
        ans = TRUE;

    return ans;
}
