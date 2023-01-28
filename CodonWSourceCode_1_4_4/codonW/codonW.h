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


#define ARB_UNIT          100                  /* used to define the array*/
#define MAX_GENE          (ARB_UNIT*3)         /* seq, which holds readin */
#define LINE_LENGTH       (ARB_UNIT+100)       /* sequence data           */
#define GARG_EXACT     0x800                   /* used in function gargs  */
#define GARG_NEXT      0x1000                  /* used in function gargs  */
#define GARG_THERE     0x2000                  /* used in function gargs  */
#define GARG_SUBSQ     0x4000                  /* used in function gargs  */
#define MAX_ARGS       100                     /* used in function gargs  */
/*                                                debugging code          */
#define debug_         printf("Got to %i\n",debugger++); 
#define debug(x)       printf( #x " = %d", x);
/*                                                defile the macro pause  */
#define pause {fprintf(stderr,"\nPress return or enter to continue -> ");gets(pm->junk);}
#define MAX_FILENAME_LEN 90                   /* max filename             */

/* define the structures used within codonW                               */
typedef struct {
  char   *des;
  char   *typ;
  int    ca[65];
} GENETIC_CODE_STRUCT;                        /* genetic code information */  

typedef  struct { 
  char *aa1[22];                              /* 1 letter AA code         */ 
  char *aa3[22];                              /* 3 letter AA code         */  
  char *cod[65];                              /* 3 letter name of codons  */
} AMINO_STRUCT;                               

typedef struct {
  float hydro[22];                            /* hydropathicity values    */ 
  int   aromo[22];                            /* aromaticity values       */
} AMINO_PROP_STRUCT;
  
typedef struct  {
  char   *des;                                /* store a description      */
  char   *ref;                                /*       a reference        */
  char   fop_cod[65];                         /* the optimal codons       */
} FOP_STRUCT; 

typedef struct {
  char   *des;                                /* store a description      */
  char   *ref;                                /*       a reference        */
  float  cai_val[65];                         /* the CAI w values         */
} CAI_STRUCT;

typedef struct {  
char   level;                                 /* either expert or standard*/
int   axis;                                   /* how many axis to generate*/
int   rows;                                   /* how many genes in dataset*/
int   colm;                                   /* how many columns in data */
int   fop_gene;                   /* No of genes to use to ident opt codon*/
char  add_row[MAX_FILENAME_LEN];              /* file with supp sequences */
float inertia;                                /* total data inertia       */
char  codons[65];                             /* codon to be analysed     */
char  amino [22];                             /* amino acids to be COA'ed */
} COA_STRUCT;

typedef struct {
  char prog;                                  /* used to ident which prog */
  char bulk;                                  /* used to ident blk output */  
  char verbose;                          /* don't overwrite files    */
  char totals;                                /* concatenate genes ?      */
  char menu;                                  /* show a menu       ?      */
  char warn;                                  /* show sequence warning    */

  char codonW;                                /* am I codonW              */ 
  char fop;                                   /* calc index fop           */
  char cai;                                   /* calc index CAI           */
  char cbi;                                   /* calc index CBI           */
  char bases;                                 /* calc base composition    */
  char gc3s;                                  /* calc gc at sil.3rd base  */
  char gc;                                    /* calc gc                  */
  char enc;                                   /* calc enc                 */
  char sil_base;                              /* calc silent base compo   */
  char L_sym;                                 /* No of synonymous codons  */
  char L_aa;                                  /* No of amino acids        */
  char hyd;                                   /* calc hydropathicity      */
  char aro;                                   /* calc aromaticity         */
  
  char seperator;                             /* column separator         */
  char coa;                                   /* calculate a COA or not ? */
  
  char code;                                  /* which genetic code       */
  char f_type;                                /* which predefined fop val */
  char c_type;                                /* which predefined CAI val */
  
  char seq_type;                              /* DNA or Protein or CU     */
  char seq_format;                            /* Human or machine readable*/
  char curr_infilename [MAX_FILENAME_LEN];    /* input filename           */  
  char curr_outfilename[MAX_FILENAME_LEN];    /* .out filename            */   
  char curr_tidyoutname[MAX_FILENAME_LEN];    /* .blk filename            */ 
  char fop_filen[MAX_FILENAME_LEN];           /* user fop filename        */
  char cai_filen[MAX_FILENAME_LEN];           /* user CAI filename        */
  char cbi_filen[MAX_FILENAME_LEN];           /* user CBI filename        */
  char curr_logfilename[MAX_FILENAME_LEN];    /* used for logging errors  */

  char junk      [BUFSIZ+1];                  /* used to store char info  */
  char messages  [300];                       /* used to constuct messgs  */
  char analysis_run;                          /* has CodonW actually run  */

  int  term_length;                           /* how many lines are there */
                                              /* file pointers            */
  FILE *inputfile;                            /* input file               */ 
  FILE *outputfile;                           /* .out file                */
  FILE *tidyoutfile;                          /* .blk file                */
  FILE *cuout;                                /* codon usage output       */
  FILE *fopfile;                              /* fop input values         */
  FILE *caifile;                              /* cai input values         */  
  FILE *cbifile;                              /* cbi input values         */
  FILE *logfile;                              /* log file name            */
  FILE *my_err;                               /* pointer for err stream   */
  
  FILE *fcoa_in;                                 
  FILE *fcoa_out;
} MENU_STRUCT ;


#ifndef DECOSF
#define DEBUG                                 /* include debug  code      */ 
#endif

#ifndef TRUE
#define TRUE 1                                /* for dumb compilers       */
#endif

#ifndef FALSE
#define FALSE 0                               /* for dumb compilers       */
#endif


/* these handle how to delete files, and blank the screen                 */
#if defined _WINDOWS || defined _WIN32
# define deletefile(x) _unlink(x)
# define clearscr(x) {int n; for(n=0; n<x ;n++) printf("\n");}
#elif defined  _DOS
# define deletefile(x) _unlink(x)
# define clearscr(x) system("cls");
#else
# define deletefile(x)  remove(x)
#if defined DEBUG
# define clearscr(x) {int n; for(n=0; n<x ;n++) printf("\n");}
#else
# define clearscr(x) system("clear");
#endif
#endif

#ifdef ORIG_DEFS                                 /* declare only once     */ 
char Revision[] = "1.4.4";                       /* version               */
char Update[]   = "$Date: 2005/05/11 21:43:49 $";/* date                  */
char Author[]   = "$Author: johnfpeden $";       /* author                */
char  title[100];                                /* sequence description  */
char  long_seq;                                  /* length of seq title   */ 
char  last_base;
long int ncod[65];
long int naa[23];
long int din[3][16];
long int codon_tot;
long int master_ic;
long int fl_pos_start;        
long int fl_pos_curr;
long int GC_TOT;
long int AT_TOT;
long int AA_TOT;
long int IUBC_TOT;
long int GAP_TOT; 
long int num_sequence;
long int num_seq_int_stop;
long int non_std_char;
long int tot;
int last_aa = 0;
int reg = 1;
int valid_stops;
int valid_start;
int fram;              
int *da;
int *ds;

AMINO_STRUCT         *paa;                       /* pointer to structs   */
GENETIC_CODE_STRUCT  *pcu;
FOP_STRUCT           *pfop;
FOP_STRUCT           *pcbi;
CAI_STRUCT           *pcai;
MENU_STRUCT          *pm;
COA_STRUCT           *pcoa;
AMINO_PROP_STRUCT    *pap;


                                                /* declare default values */           
COA_STRUCT coa={
'n',                                            /* level                  */
4,                                              /* axis                   */
0,                                              /* rows  or genes         */
64,                                             /* colms                  */
-5,               /* fop_gene (if number is negative implies a percentage)*/ 
"",                                             /* add_row                */
(float) 0.00                                    /* inertia                */
};       

int NumGeneticCodes=8;                          /* used in menu.c         */
                                                /* No. of predefined codes*/

                                                /* define genetic codes   */    
GENETIC_CODE_STRUCT  cu[] = { 
  "Universal Genetic code",
  "TGA=* TAA=* TAG=*",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 11,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  3,  8, 15, 20,  4,  8, 15, 20,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
  "Vertebrate Mitochondrial code",
  "AGR=* ATA=M TGA=W",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 19,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  4,  8, 15, 11,  4,  8, 15, 11,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
  "Yeast Mitochondrial code",
  "CTN=* ATA=M TGA=W",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 19,  2,  6, 11, 19,
  8,  7, 12, 20,  8,  7, 12, 20,  8,  7, 13, 20,  8,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  4,  8, 15, 20,  4,  8, 15, 20,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
  "Filamentous fungi Mitochondrial code",
  "TGA=W",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 19,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  3,  8, 15, 20,  4,  8, 15, 20,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
  "Insects and Plathyhelminthes Mitochondrial code",
  "ATA=M TGA=W AGR=S",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 19,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  4,  8, 15,  6,  4,  8, 15,  6,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
 "Nuclear code of Cilitia",
 "UAA=Q=Gln  UAG=Q",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 13, 11,  2,  6, 13, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  3,  8, 15, 20,  4,  8, 15, 20,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
 "Nuclear code of Euplotes",
 "UGA=C",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 18,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  3,  8, 15, 20,  4,  8, 15, 20,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21,
  "Mitochondrial code of Echinoderms",
  "UGA=W AGR=S AAA=N",
  0,
  1,  6, 10, 18,  1,  6, 10, 18,  2,  6, 11, 19,  2,  6, 11, 19,
  2,  7, 12, 20,  2,  7, 12, 20,  2,  7, 13, 20,  2,  7, 13, 20,
  3,  8, 14,  6,  3,  8, 14,  6,  3,  8, 14,  6,  4,  8, 15,  6,
  5,  9, 16, 21,  5,  9, 16, 21,  5,  9, 17, 21,  5,  9, 17, 21
};
                                            /* define amino acid info     */
AMINO_STRUCT amino_acids ={
 "X",
 "F","L","I","M","V",
 "S","P","T","A","Y",
 "*","H","Q","N","K",
 "D","E","C","W","R","G",
 "UNK",
 "Phe","Leu","Ile","Met","Val",
 "Ser","Pro","Thr","Ala","Tyr",
 "TER","His","Gln","Asn","Lys",
 "Asp","Glu","Cys","Trp","Arg","Gly",
 "BAD",
 "UUU","UCU","UAU","UGU",
 "UUC","UCC","UAC","UGC",
 "UUA","UCA","UAA","UGA",
 "UUG","UCG","UAG","UGG",
 "CUU","CCU","CAU","CGU",
 "CUC","CCC","CAC","CGC",
 "CUA","CCA","CAA","CGA",
 "CUG","CCG","CAG","CGG",
 "AUU","ACU","AAU","AGU",
 "AUC","ACC","AAC","AGC",
 "AUA","ACA","AAA","AGA",
 "AUG","ACG","AAG","AGG",
 "GUU","GCU","GAU","GGU",
 "GUC","GCC","GAC","GGC",
 "GUA","GCA","GAA","GGA",
 "GUG","GCG","GAG","GGG"
};

int NumFopSpecies=8;                             /* again used in menu.c  */
                                                 /* predefined fop info   */   
FOP_STRUCT  fop[] = { 
  "Escherichia coli",
  "Ikemura (1985) Mol. Biol. Evol. 2:13-34 (updated by INCBI 1991)",
0,2,3,2,2,3,3,3,3,2,2,2,2,2,2,2,2,
  2,2,2,3,2,2,3,3,2,2,2,2,3,3,3,2, 
  2,3,2,2,3,3,3,3,2,2,3,2,2,2,2,2,
  3,3,2,3,2,2,3,3,2,2,3,2,2,3,2,2,
  "Bacillus subtilis ",
  "Sharp et al (1990) Genetics & Biotech of Bacilli vol3 pp89-98",
0,2,3,2,2,3,1,3,2,2,2,2,2,2,1,2,2,
  3,3,2,3,2,1,2,3,2,3,3,1,2,2,2,1, 
  2,3,2,2,3,1,3,2,1,2,3,2,2,2,2,1,
  3,3,2,3,2,1,3,2,3,2,3,2,2,2,2,1,
  "Dictyostelium discoideum ",
  "Sharp and Devine (1989) Nucl. Acids Res 17:5029-5039)",
0,2,2,2,2,3,2,3,2,2,2,2,2,2,2,2,2,
  2,2,2,3,3,2,3,2,2,3,3,2,2,2,2,2, 
  2,2,2,2,3,3,3,2,2,2,2,2,2,2,3,2,
  2,2,2,3,3,3,2,2,2,2,3,2,2,2,2,2,
  "Aspergillus nidulans ",
  "Lloyd and Sharp (1991) Mol. Gen. Genet 230: 288-294",
0,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,
  2,2,2,3,3,3,3,3,2,2,2,2,2,2,3,2, 
  2,2,2,2,3,3,3,2,2,2,2,2,2,2,3,2,
  2,3,2,3,3,3,3,2,2,2,2,2,2,2,3,2,
  "Saccharomyces cerevisiae ",
  "Sharp and Cowe (1991) Yeast 7:657-678",
0,2,3,2,3,3,3,3,2,2,2,2,2,3,2,2,2,
  2,2,2,2,2,2,3,2,2,3,3,2,2,2,2,2, 
  3,3,2,2,3,3,3,2,2,2,2,3,2,2,3,2,
  3,3,2,3,3,2,3,2,2,2,3,2,2,2,2,2,
  "Drosophila melanogaster",
  "Shields et al. (1988) Mol Biol Evol 5: 704-716",
0,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,2,
  2,2,2,3,2,3,3,3,2,2,2,2,3,2,3,2, 
  2,2,2,2,3,3,3,2,2,2,2,2,2,2,3,2,
  2,2,2,2,3,3,3,3,2,2,2,2,3,2,3,2,
 "Caenorhabditis elegans",
 "Stenico, Lloyd and Sharp Nuc. Acids Res. 22: 2437-2446(1994)",
0,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,2,
  3,2,2,3,3,2,3,3,2,3,2,2,2,2,2,2, 
  2,2,2,2,3,3,3,2,2,2,2,2,2,2,3,2,
  2,3,2,2,3,3,3,2,2,2,2,3,2,2,3,2,
 "Neurospora crassa",
 "Lloyd and Sharp (1993)",
0,2,3,2,2,3,3,3,3,2,2,2,2,2,2,2,2,
  2,2,2,3,3,3,3,3,2,2,2,2,2,2,3,2,
  2,3,2,2,3,3,3,2,2,2,2,2,2,2,3,2,
  2,2,2,3,3,3,3,3,2,2,2,2,2,2,3,2
};

int NumCaiSpecies=3;                              /* used in menu.c       */
CAI_STRUCT cai[]= {                               /* array of cai structs */
  "Escherichia coli",
  "No reference",
  0.000F,
    0.296F,1.000F,0.239F,0.500F,1.000F,0.744F,1.000F,1.000F,
    0.020F,0.077F,0.000F,0.000F,0.020F,0.017F,0.000F,1.000F,
    0.042F,0.070F,0.291F,1.000F,0.037F,0.012F,1.000F,0.356F,
    0.007F,0.135F,0.124F,0.004F,1.000F,1.000F,1.000F,0.004F,
    0.185F,0.965F,0.051F,0.085F,1.000F,1.000F,1.000F,0.410F,
    0.003F,0.076F,1.000F,0.004F,1.000F,0.099F,0.253F,0.002F,
    1.000F,1.000F,0.434F,1.000F,0.066F,0.122F,1.000F,0.724F,
    0.495F,0.586F,1.000F,0.010F,0.221F,0.424F,0.259F,0.019F,
  "Bacillus subtilis",
  "No reference",
  0.00F,
     0.571F,1.000F,0.500F,1.000F,1.000F,0.021F,1.000F,1.000F,
     1.000F,0.458F,0.000F,0.000F,0.036F,0.021F,0.000F,1.000F,
     0.857F,1.000F,1.000F,1.000F,0.143F,0.071F,0.083F,0.609F,
     0.500F,0.714F,1.000F,0.022F,0.071F,0.143F,0.214F,0.043F,
     0.500F,1.000F,0.417F,0.125F,1.000F,0.033F,1.000F,0.208F,
     0.071F,0.867F,1.000F,0.435F,1.000F,0.200F,0.097F,0.022F,
     1.000F,1.000F,0.417F,0.955F,0.188F,0.025F,1.000F,0.773F,
     0.750F,0.275F,1.000F,1.000F,0.438F,0.125F,0.412F,0.045F,
  "Saccharomyces cerevisiae",
  "Sharp and Cowe (1991) Yeast 7:657-678",
  0.00F,
    0.113F,1.000F,0.071F,1.000F,1.000F,0.693F,1.000F,0.077F,
    0.117F,0.036F,0.000F,0.000F,1.000F,0.005F,0.000F,1.000F,
    0.006F,0.047F,0.245F,0.137F,0.003F,0.009F,1.000F,0.002F,
    0.039F,1.000F,1.000F,0.002F,0.003F,0.002F,0.007F,0.002F,
    0.823F,0.921F,0.053F,0.021F,1.000F,1.000F,1.000F,0.031F,
    0.003F,0.012F,0.135F,1.000F,1.000F,0.006F,1.000F,0.003F,
    1.000F,1.000F,0.554F,1.000F,0.831F,0.316F,1.000F,0.020F,
    0.002F,0.015F,1.000F,0.002F,0.018F,0.001F,0.016F,0.004F
};


AMINO_PROP_STRUCT amino_prop={                   /* amino acid properties */
  0.00F,                
  2.80F,3.80F,4.50F,1.90F,4.20F,                 /* hydropathicity values */
  -0.8F,-1.6F,-0.7F,1.80F,-1.3F,   
  1.00F,-3.2F,-3.5F,-3.5F,-3.9F,
  -3.5F,-3.5F,2.50F,-0.9F,-4.5F,
  -0.4F,
  0,
  1,0,0,0,0,                                     /* am i aromatic ?       */ 
  0,0,0,0,1,
  0,0,0,0,0,
  0,0,0,1,0,0 
};


MENU_STRUCT Z_menu={  /* define all manner of default values              */
  FALSE,              /* prog                                             */
  'X',                /*This default is set in proc_commline to CU        */
  TRUE ,              /*verbose                                      */
  FALSE,              /*totals                                            */
  TRUE,               /*menu interface                                    */
  TRUE,               /*warnings about sequence data are to be displayed  */
  FALSE,              /*codons                                            */  
  FALSE,              /*fop                                               */  
  FALSE,              /*cai                                               */ 
  FALSE,              /*cbi                                               */ 
  FALSE,              /*bases                                             */  
  FALSE,              /*gc3s                                              */
  FALSE,              /*gc                                                */  
  FALSE,              /*enc                                               */
  FALSE,              /* silent base                                      */ 
  FALSE,              /* Length silent codons                             */
  FALSE,              /* length in codons                                 */
  FALSE,              /* hydrophobicity                                   */
  FALSE,              /* aromaticity                                      */
    
  ' ',                /* default seperator                                */
   
  FALSE,              /* coa                                              */
   
  0,                  /* genetic code                                     */
  0,                  /* type of fop_species                              */ 
  0,                  /* type of cai_species                              */
  
  FALSE,              /* sequence type                                    */
  'H',                /* Sequence format                                  */
  "",                 /* current input file name                          */
  "",                 /* current output file name                         */
  "",                 /* current tidy outfile name                        */
  "",                 /* current fop input file name                      */
  "",                 /* current cai input file name                      */
  "",                 /* current sbi input file name                      */
  "",                 /* log all stderr output to a file                  */
  "",                 /* Null the string junk                             */  
  "",                 /* Null the string messages                         */
   
  FALSE,              /* was analysis run                                 */     
  24,                 /* current number of lines (height of ) screen      */ 
  
  NULL,               /* Null pointer input file                          */   
  NULL,               /* Null pointer outputfile                          */   
  NULL,               /* Null pointer tidyout file                        */   
  NULL,               /* Null codon usage file                            */
  NULL,               /* Null pointer fopfile                             */
  NULL,               /* Null pointer caifile                             */ 
  NULL,               /* Null pointer cbifile                             */ 
  NULL,               /* Null pointer the logfile name                    */
  NULL,               /* assign NULL pointer to my_err                    */
  NULL,               /* Null pointer fcoa_in                             */
  NULL                /* Null pointer fcoa_out                            */
};


#else                 /* already been defined so declare as externals     */      

extern AMINO_STRUCT         *paa;
extern GENETIC_CODE_STRUCT  *pcu;
extern FOP_STRUCT           *pfop;
extern FOP_STRUCT           *pcbi;
extern CAI_STRUCT           *pcai;
extern MENU_STRUCT          *pm;
extern COA_STRUCT           *pcoa;
extern AMINO_PROP_STRUCT    *pap;

#if defined (_WINDOWS) || defined (_DOS) 
 extern   CAI_STRUCT          /*_near*/ cai[];       /* some MS compilers  */
 extern   GENETIC_CODE_STRUCT /*_near*/ cu[];        /* want these to be   */ 
 extern   FOP_STRUCT          /*_near*/ fop[];       /* declared as _near  */
#else                
 extern CAI_STRUCT                cai[];
 extern GENETIC_CODE_STRUCT        cu[];
 extern FOP_STRUCT                fop[];
#endif
 extern COA_STRUCT                coa;
 extern AMINO_STRUCT              amino_acids;
 extern AMINO_PROP_STRUCT         amino_prop;
 extern MENU_STRUCT               Z_menu; 

 extern char Revision[];                             /* version string    */
 extern char Update[];
 extern char Author[];
 extern char  title[100];
 extern char  long_seq;
 extern char  last_base;

 extern long int ncod[65];
 extern long int naa[23];
 extern long int din[3][16];
 extern long int codon_tot;
 extern long int master_ic;
 extern long int fl_pos_start;
 extern long int fl_pos_curr;
 extern long int GC_TOT;
 extern long int AT_TOT;
 extern long int AA_TOT;
 extern long int IUBC_TOT;
 extern long int GAP_TOT; 
 extern long int num_sequence;
 extern long int num_seq_int_stop;
 extern long int non_std_char;
 extern long int tot;
 extern int last_aa;
 extern int reg;
 extern int valid_stops;
 extern int valid_start;
 extern int fram;      
 extern int *da;
 extern int *ds;
 extern int NumGeneticCodes;
 extern int NumFopSpecies;
 extern int NumCaiSpecies;
#endif

/****************** Function type declarations *****************************/

FILE *open_file    ( char *info, char *default_name, char *mode, 
                     int  verbose );

int*  how_synon    ( void );
int*  how_synon_aa ( void );
int*  how_synon    ( void );
int*  how_synon_aa ( void );

int codon_usage_tot( char *seq, long int how_many); 
int ident_codon    ( char *codon );
int codon_usage_out( FILE *fblkout, long int *ncod,int last_aa,
                     int valid_stops, char *info);
int rscu_usage_out ( FILE *fblkout, long int *ncod,long int *naa);
int raau_usage_out ( FILE *fblkout, long int *naa );
int aa_usage_out   ( FILE *fblkout, long int *naa );
int cai_out        ( FILE *foutput, long int *ncod); 
int cbi_out        ( FILE * foutput, long int *ncod, long int *naa );
int fop_out        ( FILE *foutput, long int *ncod);
int hydro_out      ( FILE *foutput, long int *naa );
int aromo_out      ( FILE *foutput, long int *naa );
int toutput        ( FILE *fblkout, char *seq );
int output_long    ( FILE *fblkout, char *seq );
int cutab_out      ( FILE *fblkout, long *ncod, long *naa);
int dinuc_out      ( FILE *fblkout, char *title  );
int fileclose      ( FILE **file_pointer );
int clean_up       ( long int *ncod,long int *naa );
int initilize_point( char code , char fop_type, char cai_type );
int initilize_coa  ( char code );
int proc_comm_line ( int *argc, char ***arg_list);
int my_exit        ( int exit_value, char *message );          
int printinfo      ( void); 

int dinuc_count    ( char *seq , long int tot );
int tidy           ( FILE *finput , FILE *foutput , FILE *fblkout, 
                     FILE *fcoaout ) ;  
int chelp ( char *help );

long int codon_error( int last_aa, int valid_stops, char *title, 
                      char error_level);

float  enc_out      ( FILE *foutput, long int *ncod, long int *naa);
double inertot      ( void);

char* get_aa        ( int one_or_3_letter , char* the_dna_word);
char* garg          ( int argc, char *argv[], const char *targ, int mode);
char  coa_raw_out   ( FILE *fcoaout, long *ncod, long *naa, char *title);
char  WasHelpCalled ( char * input); 

void sorted_by_axis1( double *ax1, int *sortax1, int lig);
void highlow        ( long int *low , long int *high ,FILE *summ );
void menu_1         ( void);
void menu_2         ( void);
void menu_3         ( void);
void menu_4         ( void);
void menu_5         ( void);
void menu_6         ( void); 
void menu_7         ( void);
void menu_8         ( void);
void menu_coa       ( void);
void welcome        ( void);
void menu_initial   ( void);

void asummary       ( void);
void tester         ( void);
void vecalloc       ( double **vec, int n);
void vecalloc       ( double **vec, int n);
void writevec       ( double *v1, FILE *fic);
void lecmat         ( double **tab, char *nfic);
void freetab        ( double **tab);
void freevec        ( double *vec);
void taballoc       ( double ***tab, int l1, int c1);
void lecvec         ( double *v1, char *nfic);
void ecrmat         ( double **tab, char *nfic);
void ecrvec         ( double *v1, char *nfic);
void scalmat        ( double **tab, double r);
void scalvec        ( double *v1, double r);
void sqrvec         ( double *v1);
void prodmatAAtB    ( double **a, double **b);
void prodmatABC     ( double **a, double **b, double **c);
void prodmatAtAB    ( double **a, double **b);
void ecrmatred      ( double **tab, int c1, char *nfic);
void readvec        ( double *v1, FILE *fic);
void lecvalpro      ( double *v1, char *nfic);
void writescal      ( double r,  FILE *fic);
void editvalpro     ( FILE *ficlist, double *vp, int n, double s);
void DiagoRC        ( FILE *summary);
void gc_out         ( FILE *foutput, FILE *fblkout, int which);
void base_sil_us_out( FILE *foutput, long int *ncod,long int *naa);
void bintext        ( char *nfice , char *nfics);
void select_coa     ( char choice); 
void textbin        ( char *filein , char *fileout);
void colmout        ( char *nfice, char *nfics,AMINO_STRUCT *paa,
                                   FILE *summary);   
void output         ( char *seq ,  FILE *foutput , FILE* fblkout ,
                                   FILE *fcoaout);
void rowout         ( char *nfice, char *nfics, char *ncout, FILE *summary);
void PrepAFC        ( char *nfic);
void inertialig     ( char *inertia_out, char *filen, FILE *summary);
void inertiacol     ( char *inertia_out, FILE *summary);
void selectcol      ( char *nfic , double *col, int numcol);
void gen_cusort_fop ( int *sortax1, int lig , FILE *fnam ,FILE *summ ); 
void dot            ( int y    ,  long int period ); 
void DiagoComp      ( int n0, double **w, double *d, int *rang);
void suprow         ( int num_seq,char *nficvp,char *nfictasup,
                      char *nficlisup,char *option, FILE *summary);
void main_menu      ( int c );



































