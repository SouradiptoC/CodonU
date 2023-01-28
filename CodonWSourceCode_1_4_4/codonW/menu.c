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
#include <stdlib.h> 
#ifdef _WINDOWS
#include <process.h>
#endif
#include <stdio.h>    
#include <string.h>
#include <ctype.h>
#include "codonW.h" 

/************** Main menu   **********************************************/
/* Drives the menu system                                                */
/*************************************************************************/
void    main_menu ( int menu )
{
    switch ( menu ) {                                 /* go to menu X    */
    case 0:
        menu_initial();
        break;
    case 1:
        menu_1();
        break;
    case 2:
        menu_2();
        break;
    case 3:
        menu_3();
        break;
    case 4:
        menu_4();
        break;
    case 5:
        menu_5();
        break;
    case 6:
        menu_6();
        break;
    case 7:
        menu_7();
        break;
    case 8:
        menu_8();
        break;
    case 9:
        printinfo();
        welcome();
        pause;
        clearscr(pm->term_length);
        break;                     
    default:
        fprintf ( stderr,"ERROR: Unrecognised menu in main_menu\n");
        break;
    }
}


/* This is the first menu presented when running CodonW                   */ 

void menu_initial (void)
{
    int loop = TRUE;
    int c;
    
    while (loop) {                                             /* loop    */
        printf (" Initial Menu \n");
        printf (" Option\n\t (1) Load sequence file\n"); 
        
/*      printf ("\t (2) Check sequence file for redundancy\n");           */
        printf ("\t ( )\n");
        printf ("\t (3) Change defaults\n");
        printf ("\t (4) Codon usage indices\n");
        printf ("\t (5) Correspondence analysis\n");

/*      printf ("\t (6) Basic statistics\n");                             */
        printf ("\t ( ) \n");

        printf ("\t (7) Teach yourself codon usage\n");
        printf ("\t (8) Change the output written to file\n");
        printf ("\t (9) About C-codons\n");
        printf ("\t (R) Run C-codons \n");
        printf ("\t (Q) Quit \n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");

        gets(pm->junk);
        
        if (isalpha((int)pm->junk[0])) {
            c = toupper( (int) pm->junk[0]);

            switch  (c) {
            case 'Q':       
                my_exit(2,"main menu");
                break;                
            case 'R':    
                /* test that all the required files are opened             */
                if ( pm->inputfile && pm->outputfile && pm->tidyoutfile)
                    loop = FALSE;
                else {                  
                    printf("Not all required files are open\n");
		            printf("About to open input and output files\n");
		            pause;
		            main_menu(1);   
		            loop = FALSE;
                 } 
		    break;
            case 'H':                                              /* help */     
                 chelp ( "main_menu" );       
                 break;   
            default:
                fprintf( stderr, "The answer %s is not valid\n", pm->junk);
                pause;
                break;
           }                                            /* end of switch c */     
        } else if (isdigit((int) pm->junk[0])) {
            c =    atoi( pm->junk);
            if (c > 0 && c <= 9 )  
                main_menu( (int) c );
            else 
                fprintf( stderr, "The answer %s is not valid\n", pm->junk);
        }
        clearscr(pm->term_length);
    }
    return;
}

/************************* menu_1 ******************************************/
/* Opens input and output files                                            */
/* It tests if a sequence file is already in memory                        */
/* if so you have the option to reopen the same file when loaded the       */
/* pm->file_loaded is set to true and the 20 characters of the new filename*/
/* are stored                      							               */
/***************************************************************************/
void    menu_1 (void)
{
    char    root[MAX_FILENAME_LEN];
    int n;
    
    clearscr(pm->term_length);
    printf (" Loading sequence menu (type h for help)\n");

    if ( strlen(pm->curr_infilename) ) {
        printf ( "The current active file is \"%s\"\n",pm->curr_infilename);
        fileclose(&pm->inputfile);
        if (!(pm->inputfile = open_file("input sequence file",
            pm->curr_infilename, "r", FALSE)))
            my_exit(1,"menu 1");
    } else {
        printf( " No sequence file is currently loaded\n");
        if (!(pm->inputfile = open_file("input sequence file\t",
            "input.dat", "r", FALSE)))
            my_exit(1,"menu 1");
    }
    /* copies the filename into pm->curr_infilename                        */
    /* next finds the root of this filename                                */
    /* which is used to construct other filenames                          */


    strncpy(pm->curr_infilename, pm->junk, MAX_FILENAME_LEN - 1);
    strncpy(root, pm->curr_infilename    , MAX_FILENAME_LEN - 5);

    /* open the .out filename                                              */
    for (n = (int) strlen(root); n && root[n] != '.' ; --n);
    if (n)        root[n] = '\0';           /* define root of the filename */

    if ( strlen(pm->curr_outfilename)) {
        printf( "\nThe previous  output file was \"%s\"\n", 
            pm->curr_outfilename );
        fclose( pm->outputfile);
    } 
        if (!(pm->outputfile = open_file("output sequence file\t",
            strcat(root, ".out"), "w", (int)pm->verbose)))
            my_exit(1,"output menu1");

    /* open the .blk filename                                              */

    strncpy(pm->curr_outfilename, pm->junk, MAX_FILENAME_LEN - 1);
    strncpy(root, pm->curr_infilename     , MAX_FILENAME_LEN - 5);
 
    for (n = (int) strlen(root); n && root[n]!='.'  ; --n);
    if ( n  ) root[n] = '\0';                   /* find root of filename  */

    if ( strlen(pm->curr_tidyoutname)) {
        printf( "\nThe previous bulk output file was \"%s\"\n", 
            pm->curr_tidyoutname );
        fclose( pm->tidyoutfile);
   }
        if (!(pm->tidyoutfile = open_file("bulk output file\t",
            strcat(root, ".blk"), "w", (int) pm->verbose)))
            my_exit(1,"tidyout menu1");
   
    strncpy(pm->curr_tidyoutname, pm->junk, MAX_FILENAME_LEN - 1);

    clearscr(pm->term_length);
    return;
}

/************************* menu_2 ******************************************/
/* Not currently implemented                                               */
/***************************************************************************/
void    menu_2 (void)
{
    int loop = TRUE;
    int  c;
    
    clearscr(pm->term_length);
    while ( loop ) {
        printf (" Menu 2 \n");
        printf (" Purifying sequences menu\n");
        printf ("\t ( ) Sorry currently unimplemented \n");
        printf ("\t (X) Exit this menu\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");
        gets(pm->junk);
        clearscr(pm->term_length);

        if (isalpha((int)pm->junk[0]) || pm->junk[0]=='\0' ) {
            c = toupper( (int) pm->junk[0]);
            switch ( c ) {
            case 'Q':
                my_exit(2,"menu 2");
                break;
            case 'X':
            case '\0':
                return;
            case 'H':
                chelp("menu_2");
                break;
            default:
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                pause;
                break;
            }
        } 
    }
    return;
}

/************************* menu_3 ******************************************/
/* To improve flexibility, many of the default values used internally by   */
/* CodonW (defined in the header file codonW.h) can be altered at runtime  */
/* using this menu. Ten default values can be customised.                  */
/***************************************************************************/
void    menu_3 (void)
{
    int loop = TRUE;
    int i;
    int c;
    
    clearscr(pm->term_length);
    while (loop) {
        printf (" Changing defaults\n");
        printf (" Options\n");
        printf (" %-40.40s", "(1) Change the ASCII delimiter in output");
        printf ("{%s}\n", 
            (pm->seperator == ' ' ) ? "space" : 
            (pm->seperator == '\t') ? "tab" : 
            (pm->seperator == ',' ) ? "," : 
            "ERROR" );

        printf (" %-40.40s", "(2) Run silently, No Warnings");
        printf ("{%s}\n", (pm->verbose) ? "FALSE" : "TRUE");
        printf (" %-40.40s", "(3) Log warnings/information to a file");
        printf ("{%s}\n", (strlen(pm->curr_logfilename) > 1) ? "TRUE" : 
                "FALSE");
        printf (" %-40.40s", "(4) Number of lines on screen");
        printf ("{%d}\n", pm->term_length);
        printf (" %-40.40s", "(5) Change the genetic code");
        printf ("{%s}\n", cu[pm->code].des);
        printf (" %-40.40s", "(6) Change the Fop/CBI values");
        printf ("{%s}\n", fop[pm->f_type].des);
        printf (" %-40.40s", "(7) Change the CAI values");
        printf ("{%s}\n", cai[pm->c_type].des);
        printf (" %-40.40s", "(8) Output Human or Computer readable");
        printf ("{%s readable}\n", (pm->seq_format == 'M') ? "Computer" : 
                "Human"); 
        printf (" %-40.40s", "(9) Concatenate or individual genes");
        printf ("{%s genes}\n", (pm->totals == TRUE ? "concatenate":
                "individual"));      
        printf (" %s", "(10) Correspondence analysis defaults\n");
    
        printf (" (X) Return to previous menu\n");
        printf ("Choices enclosed with curly brackets are the current "
                "defaults\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");
        gets(pm->junk);
        clearscr(pm->term_length);

        if (isalpha((int) pm->junk[0])|| pm->junk[0]=='\0') {
            switch (c = toupper((int) pm->junk[0])){
              case 'Q': 
                my_exit(2,"menu 3");           /* decided to quit program  */
                break;
              case 'H':
                chelp("menu_3");
                break;
              case 'X':
              case '\0':
                return; /*     way out of loop is X or blank line          */
                break;
              default:
                fprintf(stderr,"The answer %s is not a valid\n", pm->junk);
                pause;
                continue;
                break;
                }
        }

        c=0;
        if (isdigit((int)pm->junk[0]))
            c = atoi(pm->junk);
        if ( c <= 0 && c > 10 ) {
            fprintf( stderr, "The answer %s is not valid\n", pm->junk);
            continue;
        }

        switch ((int) c) {
        case 1:
            clearscr(pm->term_length);
            printf (" The current separator is  \"%s\"\n",  
                (pm->seperator == ' ' ) ? "space" : 
                (pm->seperator == '\t') ? "tab" : 
                (pm->seperator == ',' ) ? "," : 
                "ERROR" );
            printf (" Please select a new separator \t:");
            gets(pm->junk);
            c = pm->junk[0];             /* take first character of string */

            if ( strchr ("\t, ", (int)c) == NULL || c == '\0' ) {
                                     /* remember the \0 is in every string */
                printf( "WARNING: The chosen separator %s is unsuitable\n", 
                        pm->junk);
                printf( "\tSeparator is unchanged try comma,tab "
                        "or space\n\n");
            } else
                pm->seperator = (char) c;  /* specify the column separator */

            break;
        case 2:                            /* warn about overwriting files?*/
            clearscr(pm->term_length);
            pm->verbose = (char) ((pm->verbose) ? FALSE : TRUE);
            pm->warn         = (char) ((pm->warn        ) ? FALSE : TRUE);
            break;
        case 3:                            /* redirect errors to a file    */
            if ( strlen(pm->curr_logfilename) > 1 ) {
                strcpy(pm->curr_logfilename , "" );   /* blank logfilename */
                pm->my_err = stderr;                  /* redirects errors  */
                                                      /* to stderr         */
                fclose(pm->logfile);                  /* close logfile     */ 
            } else {
                                       /* open logfile and redirect stderr */
                if (!(pm->logfile = open_file("log filename        \t",
                    "warning.log", "w", (int) pm->verbose)))
                    my_exit(1," open log file menu 3");
                pm->my_err = pm->logfile;
                strncpy(pm->curr_logfilename, pm->junk, MAX_FILENAME_LEN-1);
            }                                                 /* end of if */
            break;

        case 4:                                       /* No of line on term*/
            printf("Please give the new height of the screen [%i] ", 
                    pm->term_length);
            gets(pm->junk);
            if ( isdigit( (int) pm->junk[0]))
                pm->term_length = atoi(pm->junk) ;
            break;

        case 5:                                      /*Change genetic code */
            clearscr(pm->term_length);
            printf(" Genetic codes currently supported are\n");
           /* NumGeneticCodes is given in codonW.h                         */
           for ( i = 0 ; i < NumGeneticCodes ; i++) {
                (pm->code == i) ? printf ( " (%i) {%-45.45s %-17.17s}", i, 
                    cu[i].des, cu[i].typ) : 
                    printf ( " (%i)  %-45.45s %-17.17s ", i, cu[i].des, 
                        cu[i].typ) ;
                printf("\n");
            }
            printf("Choice enclosed with curly brackets is "
                   "the current code\n");
            printf("Please select a new code [no change]\n");
            gets(pm->junk);
            if ( isdigit( (int) pm->junk[0]) ) {
                c = (char)atoi(pm->junk);
                if ( c > 0 && c < NumGeneticCodes && pm->code!= (char) c ){ 
                    pm->code = (char) c;
                    initilize_point(pm->code,pm->f_type, pm->c_type);  
                    }
            }
            break;

        case 6:                                     /*Change optimal codons*/
            clearscr(pm->term_length);
            printf(" Fop values pre-loaded are\n");
            /* NumFopSpecies  defined with the Fop_struct in codonW.h      */
            for ( i = 0 ; i < NumFopSpecies ; i++) {
                (pm->f_type == i) ? printf (" (%i) {%-25.25s %-40.40s}", 
                    i, fop[i].des, fop[i].ref) : 
                    printf (" (%i)  %-25.25s %-40.40s ", i, fop[i].des, 
                        fop[i].ref) ;
                printf("\n");
            }
            printf ("Choice enclosed with curly brackets is the current "
                "selection\n");
            printf ("Please select a type [no change]\n");
            gets(pm->junk);
            if ( isdigit( (int) pm->junk[0]) ) {
                c = (char)atoi(pm->junk);
                if ( c > 0 && c < NumFopSpecies && pm->f_type!=(char) c) {
                        pm->f_type = (char) c;  
                        initilize_point(pm->code,pm->f_type, pm->c_type);
                }
            }
            break;

        case 7:                                      /*Change CAI w values */
            clearscr(pm->term_length);
            printf(" CAI types currently supported are\n");

            /*  NumCaiSpecies currently defined in codonW.h                */
            for ( i = 0 ; i < NumCaiSpecies ; i++) {
                (pm->c_type == i) ? printf (" (%i) {%-25.25s %-40.40s}", 
                    i, cai[i].des, cai[i].ref) : 
                    printf (" (%i)  %-25.25s %-40.40s ", i, cai[i].des, 
                        cai[i].ref) ;
                printf("\n");
            }
            printf ("Choice enclosed with curly brackets is the current "
                "selection\n");
            printf ("Please chose a new CAI [no change]\n");
            gets(pm->junk);
            if ( isdigit( (int) pm->junk[0]) ) {
                c = (char)atoi( pm->junk);

                /* if valid value and different from the current choice    */
                if (  c > 0 && c < NumCaiSpecies && pm->c_type!=(char) c){
                    pm->c_type = (char) c;
                    initilize_point(pm->code,pm->f_type, pm->c_type);
                    }
            }
            break;
       case 8:                       /* machine or human readable format  */
             clearscr(pm->term_length);
             pm->seq_format = 
                (char) (  pm->seq_format == 'M' ? 'H' : 'M'); /*toggle    */
             break;
      case 9:                        /* concatenate genes?                */
            clearscr(pm->term_length);
            pm->totals    = (char) (pm->totals == TRUE ? FALSE : TRUE); 
            break;
     case 10:                       /* change COA default then go to menu5*/
           clearscr(pm->term_length);
           if( !pm->coa ) 
                menu_5();
           else 
                menu_coa();           
           break;
     default:
            fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
            break;
     }
    }
  return;
}

/************************* menu_4 ******************************************/
/* Select which indices to calculate                                       */
/***************************************************************************/
void    menu_4 (void)
{
    char    loop = TRUE;
    char    *choices[] = {
        " ",
        "Codon Adaptation Index       (CAI)",
        "Frequency of OPtimal codons  (Fop)", 
        "Codon bias index             (CBI)",
        "Effective Number of Codons   (ENc)",
        "GC content of gene           (G+C)",
        "GC of silent 3rd codon posit.(GC3s)",
        "Silent base composition",    
        "Number of synonymous codons  (L_sym)",
        "Total number of amino acids  (L_aa )",
        "Hydrophobicity of protein    (Hydro)",
        "Aromaticity of protein       (Aromo)",
        "Select all"
    }; 
    int i,NumChoices;
    int c;
    
    
    NumChoices = (char) 12;                      /* size of choices array */

    clearscr(pm->term_length);
    while (loop) {
        printf (" Codon usage indices\n");
        printf (" Options\n");

        for (i = 1; i <= NumChoices; i++) {
            printf(" (%2i) ", i);
            switch ((int) i) {
            case 1:
                (pm->cai) ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 2:
                (pm->fop) ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 3:
                (pm->cbi) ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;         
            case 4:
                (pm->enc) ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 5:
                (pm->gc)  ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 6:
                (pm->gc3s)? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 7:
                (pm->sil_base) ?  printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 8:
                (pm->L_sym) ?  printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 9:
                (pm->L_aa)?  printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 10:
                (pm->hyd ) ? printf ("{%-45.45s}", choices[i]) : 
                printf (" %s ", choices[i]);
                break;
            case 11:    
                (pm->aro ) ? printf ("{%-45.45s}", choices[i]): 
                printf (" %s ", choices[i]);
                break;
            case 12:
                printf (" %s ", choices[i]);
                break;                      
            default:
                fprintf(stderr, "programming error \n");
                my_exit(99, "menu 4");
                break;
            }
            printf("\n");
        }
        printf (" (X)  Return to previous menu\n");
        printf ("Choices enclosed with curly brackets are the current"
                " selections\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");


        gets(pm->junk);

        if (isalpha( (int) pm->junk[0]) || pm->junk[0]=='\0') {
            switch (c = toupper( (int) pm->junk[0])){
                case 'Q': 
                 my_exit(2,"menu 4");     /* User decides to quit programme*/
                 break;
                case 'X':
                case '\0':
                    return;               /* <-back to previous menu->     */
                    break;
                case 'H':
                    chelp("menu_4");
                    continue;
                    break;
                default:
                    fprintf( stderr, "The answer %s is not a valid choice\n",
                    pm->junk);
                    continue;
                    break;
                }
        } else if (isdigit ( (int) pm->junk[0] ) ) {
            c =  atoi(pm->junk);
            switch ((int) c) {
            /* User wants to calculate CAI then we explain that it is     */
            /* dependent on the choice of CAI adaptiveness values         */
            case 1: 
                pm->cai = (char) ((pm->cai)   ? FALSE : TRUE);    
                if( pm->cai){
                clearscr(pm->term_length);
                printf("\nTo calculate CAI a reference set of highly ");
                printf("expressed genes \nmust be selected\n\n");
                printf("The reference set currently selected is that of "
                    "%s\n\n",cai[pm->c_type].des);  
                printf("See the menu 'Change defaults' to change this "
                       "selection\n\n");  
                printf("If you wish to use a personal choice of CAI "
                       "vaules\n");
                printf("\tplease continue and you will be prompted for"
                       " input\n\n");  
                pause;
                }
                break ;
            case 2: 
            /* User wants to calculate Fop then we explain that it is     */
            /* dependent on the choice of optimal codons                  */
                pm->fop = (char) ((pm->fop)   ? FALSE : TRUE); 
                if(pm->fop){   
                clearscr(pm->term_length);              
                printf("\n\nYou have chosen to calculate Fop\n\n");
                printf("To calculate Fop a set of optimal "
                       "codons must be selected\n");
                printf("The optimal codons of %s are the current selection"
                       "\n\n",fop[pm->f_type].des);  
                printf("See the menu 'Change defaults' to change Fop "
                       "selection\n\n");
                printf("If you wish to use a personal choice of Fop "
                       "vaules\n");
                printf("\tplease continue and you will be prompted for "
                       "input\n\n");            
                pause;
                }
                break ; 
           case 3: 
            /* User wants to calculate CBI then we remind then that it is */
            /* dependent on the choice of optimal codons                  */
                pm->cbi = (char) ((pm->cbi)   ? FALSE : TRUE); 
                if(pm->cbi){   
                clearscr(pm->term_length);              
                printf("\n\nYou have chosen to calculate CBI\n\n");
                printf("To calculate CBI a set of optimal "
                       "codons must be selected\n");
                printf("The optimal codons of %s are the current selection"
                       "\n\n",fop[pm->f_type].des);  
                printf("See the menu 'Change defaults' to change CBI "
                       "selection\n\n");
                printf("If you wish to use a personal choice of CBI "
                       "vaules\n");
                printf("\tplease continue and you will be prompted for "
                       "input\n\n");               
                pause;
                }
                break ;                
            case 4:                                      /* calc Nc       */
                pm->enc = (char) ( (pm->enc)   ? FALSE : TRUE);    
                break ;
            case 5:                                      /* calc GC       */
                pm->gc =  (char) ((pm->gc )   ? FALSE : TRUE);    
                break ;
            case 6:                                      /* calc GC3s     */   
                pm->gc3s =(char) ( (pm->gc3s) ? FALSE : TRUE);    
                break ;
            case 7:                                      /* calc sil base */   
                pm->sil_base = (char) ((pm->sil_base) ? FALSE : TRUE); 
                break ; 
            case 8:                                      /* No. synonyms  */
                pm->L_sym = (char) ((pm->L_sym) ? FALSE : TRUE); 
                break ; 
            case 9:                                      /* No. AminoAcids*/   
                pm->L_aa  = (char) ((pm->L_aa)  ? FALSE : TRUE); 
                break ; 
            case 10:                                     /* hydropathicity*/
                pm->hyd   =(char) ( (pm->hyd )  ? FALSE : TRUE);
                break;
            case 11:                                     /* aromatic      */
                pm->aro   = (char) ((pm->aro )  ? FALSE : TRUE);                                         
                break;
            case 12:                                     /* all the above */
                pm->cai   = (char)  TRUE;    
                pm->fop   = (char)  TRUE;
                pm->cbi   = (char)  TRUE;    
                pm->enc   = (char)  TRUE;    
                pm->gc    = (char)  TRUE;    
                pm->gc3s  = (char)  TRUE;    
                pm->sil_base 
                          = (char)  TRUE; 
                pm->L_sym = (char)  TRUE; 
                pm->L_aa  = (char)  TRUE;
                pm->hyd   = (char)  TRUE;
                pm->aro   = (char)  TRUE;          
                break ;             
            default:
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                break;
            }
        } else
            fprintf( stderr, "The answer %s is not a valid choice\n", 
            pm->junk);
    }
    return;
}

/************************* menu_5 ******************************************/
/* Select what type of COA                                                 */
/***************************************************************************/
void    menu_5 (void)
{ 
    char    *choices[] = {
        "",
        "COA on codon usage",
        "COA on RSCU",
        "COA on Amino Acid usage",
        "Do not perform a COA"
    };
    int loop = TRUE;
    int i,c,NumChoices;

    NumChoices = 4;

    clearscr(pm->term_length);

    while ( loop ) {
        printf (" Menu 5  Correspondence analysis\n");
        printf ("  Correspondence analysis (COA) \n");

        for (i = 1; i <= NumChoices; i++) {
            printf(" (%i) ", i);
            switch ((int) i) {
            case 1:
                (pm->coa=='c') ? printf ("{%-45.45s}", choices[1]):
                printf (" %s ", choices[1]);
                break;
            case 2:
                (pm->coa=='r') ? printf ("{%-45.45s}", choices[2]):
                printf (" %s ", choices[2]);
                break;
            case 3:
                (pm->coa=='a') ? printf ("{%-45.45s}", choices[3]):

                printf (" %s ", choices[3]);
                break;
            case 4:
                (pm->coa== 0 ) ? printf ("{%-45.45s}", choices[4]):
                printf (" %s ", choices[4]);
                break;                  
            default:
                fprintf(stderr, "programming error \n");
                my_exit(99,"menu 5");
                break;
            }
            printf("\n");
        }
        printf (" (X) Exit this menu\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");
        gets(pm->junk);                                              
        clearscr(pm->term_length);

        if (isalpha( (int) pm->junk[0]) ||   pm->junk[0]=='\0') {
            c =  toupper( (int) pm->junk[0]);
            switch ( c ) {
            case 'Q':
                my_exit(2,"menu 5");
                break;
            case 'X':
            case '\0':
                return;
                break;
            case 'H':
                chelp("menu_5_coa");
                continue;
                break;
            default:
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                break;
            }
        } else {
            c =  atoi(pm->junk);
            if ( c > 0 && c <= 4 ) {
            switch ((int) c){
            case 1: 
                pm->coa = 'c';                              /* COA of CU  */
                break ;
            case 2: 
                pm->coa = 'r';                              /* COA of RSCU*/
                break ;
            case 3: 
                pm->coa = 'a';                              /* COA of AA  */
                break ;
            case 4: 
                pm->coa = FALSE;        
                break;
#ifdef DEBUG
            default:
                fprintf(pm->my_err,"Error in switch in coa_raw_out\n");
#endif          
              }
            } else {
                fprintf(stderr,"The answer %s is not a valid\n", pm->junk);
                break;
            }
        }
 
     if ( pm->coa ) {  
         printf( " Do you wish to see the advanced COA menu (Y/N) [N] ");
         gets( pm->junk );

        /* Select the default codon/AAs to analyse, based on genetic code */
         initilize_coa  (pm->code);
         
         if ( (char) toupper( (int) pm->junk[0]) == 'Y' ) menu_coa(); 
         }
        
    } /* while loop */ 
    return;
}

/************************* menu_6 ******************************************/
/* Originally designed for the calculation of correlations and             */
/* other simple stats. This code is currently implemented as a perl module */
/* and is waiting to be ported to C hence the menu is unimplemented        */
/***************************************************************************/

void    menu_6 (void)
{
    int loop = TRUE;
    int c;
    
    clearscr(pm->term_length);
    while ( loop ) {
        printf (" Menu 6-Basic Stats\n");
        printf ("\n");
        printf ("\t ( ) Sorry currently unimplemented \n");
        printf ("\t (X) Exit this menu\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");
        gets(pm->junk);
        clearscr(pm->term_length);

        if (isalpha( (int) pm->junk[0])|| pm->junk[0] == '\0') {
            c =  toupper( (int) pm->junk[0]);
            switch ( c ) {
            case 'Q':
                my_exit(2,"menu 6");
                break;
            case 'X': 
            case '\0':
                return;
            case 'H':
                 chelp("menu_6");
                 break;
            default:
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                pause;
                break;
            }
        } else {
            c =  atoi(pm->junk);
            if ( c > 0 && c <= 9 )
                main_menu((int) c);
            else {
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                continue;
            }
        }
    }
    return;
}

/************************* menu_7 ******************************************/
/* This selection generates random questions about the genetic code that   */
/* has been selected. For more information see tester.c                    */
/***************************************************************************/
void    menu_7 (void)
{
    int loop = TRUE;
    int c;
    
    clearscr(pm->term_length);
    while ( loop ) {
        printf (" Menu 7 A Bit of fun \n");
        printf ("\n");
        printf (" (1) Test your knowledge of the genetic code \n");
        printf (" (X) Exit this menu\n");
        printf (" Select a menu choice, (Q)uit or (H)elp -> ");
        gets(pm->junk);
        clearscr(pm->term_length);

        if (isalpha( (int) pm->junk[0]) || pm->junk[0]=='\0') {
            c =  toupper( (int) pm->junk[0]);
            switch ( c ) {
            case 'Q':
                my_exit(2,"menu 7");
                break;
            case 'X': case '\0':
                return;

            case 'H':
                chelp("menu_7");
                continue;
                break;

            default:
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                pause;
                break;
            }
        } else {
            c = atoi(pm->junk);
            if ( c == 1 )
                tester();        /****** call tester () ********************/
            else {
                fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
                continue;
            }
        }
    }
    return;
}

/************************* menu_8 ******************************************/
/* This menu allows the selection of the output to be written to the file  */
/* .blk. Only one selection can be made at a time. However CodonW can be   */
/* rerun with the same input file but with different output options. To    */
/* make this easier each time this menu is selected the user is given the  */
/* choice of changing the output file                                      */
/***************************************************************************/


void    menu_8 (void)
{
  struct multi {                    /* struct of menu items                */
    char    *string;                /* description string                  */
    char    prog;                   /* programme name                      */
  };
  char  loop = TRUE;
  int  c;
  int  ans1,NumChoices;

  struct multi aii[] = {
    " ", ' ',              /* Initialise a single value of choices in menu */
    "Fasta format output of DNA sequence", 'T',
    "Reader format output of DNA sequence",'R',
    "Translate input file to AA sequence", 'N',
    "Codon Usage"                        , 'C',
    "Amino acid usage"                   , 'A',
    "RSCU values"                        , 'S',
    "Relative Amino Acid usage"          , 'L',
    "Dinucleotide frequencies"           , 'D',
    "Exhaustive base compostion analysis", 'B',
    "No output written to file"          , 'X' };
  
  NumChoices = 10;                            /* Number of choices in Menu */
  
            /* if there is already an output file available the user may   */
            /* select to change it                                         */

  clearscr(pm->term_length);

  /* because only one type of bulk option is permitted each time 
     codonw runs, it may be necessary to rerun with the same data
     file but changing the blk output options, if so the user
     is prompted with the choice of changing the blk filename             */

  if ( pm->analysis_run  ) {
    printf (" The current bulk output file is %s do you "
            "wish to change this (y/n) [n] ", pm->curr_tidyoutname);
    gets(pm->junk);
   
    if ( toupper( (int) pm->junk[0]) == 'Y') {
      fileclose(&pm->tidyoutfile);
    
      if (!(pm->tidyoutfile = open_file("codon usage output file",
               pm->curr_tidyoutname, "w",(int)pm->verbose)))
               my_exit(1, "menu 8");
      strncpy(pm->curr_tidyoutname, pm->junk, MAX_FILENAME_LEN - 1);
    }        /* matches  if ( !strlen (pm->junk) || toupper= ............. */
  
  } else {   /* matches  if( strlen( pm->curr_cufilename)  )               */    
    printf("Note: No output file has been selected !\n");
  }
  

  while ( loop ) {
    printf (" Menu 8\n");
    printf (" This output will be saved to %s\n\n", pm->curr_tidyoutname);
    
    for ( ans1 = 1; ans1 <= NumChoices; ans1++) {
      if (aii[ans1].prog != (char) pm->bulk)
         printf("\n\t (%2d) %s", ans1, aii[ans1].string);
      else
         printf("\n\t{(%2d) %-45.45s\t\t}", ans1, aii[ans1].string);
    }

    printf ("\n\t ( X) To return to previous menu\n");
    
    printf ("Values enclosed with curly{} brackets are the current "
            "selection\n");
    printf (" Select a menu choice, (Q)uit or (H)elp -> ");
    gets(pm->junk);
    clearscr(pm->term_length);
    
    if (isalpha( (int) pm->junk[0]) || pm->junk[0]=='\0') {
      switch (c =  toupper( (int) pm->junk[0])){
        case 'Q':  
            my_exit(2,"menu 8");         /* User decides to quit          */
            break;
        case 'X':
        case '\0':
            return;                      /* <-back to previous menu->     */      
        case 'H':
            chelp("menu_8_blk");
            continue;
            break;
        default:
          fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
          pause;
          break;
      }
    } else {
      c = atoi(pm->junk);
      if ( c > 0 && c <= NumChoices )
         pm->bulk = aii[c].prog;
      else
         fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
    }
  }                                      /* match while                  */
  return;
}


/*********************** menu_coa   ***************************************/
/* This is the advanced correspondence menu, this menu is optional, when a*/
/* a correspondence analysis is chosen, then the user is given a choice of*/
/* entering this menu                                                     */
/**************************************************************************/
void menu_coa (void)
{
  int   loop = TRUE;
  char *p;
  int c;
  int i;
  
  clearscr(pm->term_length);
  while ( loop ) {
    printf ("Advanced Correspondence Analysis\n");         
    printf (" (1) (Un)Select %s\n", (pm->coa=='a')? "amino acids": "codons");         
    printf (" (2) Change the number of axis (factors) recorded to file\n");
    printf (" (3) Add additional genes after COA\n");       
    printf (" (4) Toggle level of COA output  [%s]\n",
                  (pcoa->level=='e')? "Exhaustive":"Normal");

    if(pm->coa != 'a' )
    printf (" (5) No. genes used to identify optimal codons [%i%s]\n", 
        (pcoa->fop_gene <0)? (pcoa->fop_gene*-1): pcoa->fop_gene,
            (pcoa->fop_gene <0)? "%"                : " genes");

    printf (" (X) Exit this menu\n");
    printf (" Select a menu choice, (Q)uit or (H)elp -> ");
    gets(pm->junk);
    clearscr(pm->term_length);
    
    if (isalpha( (int) pm->junk[0]) || pm->junk[0]=='\0' ) {
      c =        toupper( (int) pm->junk[0]);
      switch ( c ) {
       case 'Q':
         my_exit(2, "menu coa");
         break;
       case 'X' :
       case '\0':
        return;
       case 'H':
        chelp("menu_coa");
        continue;
        break;
      default:
        fprintf( stderr, "The answer %s is not a valid\n", pm->junk);
        pause;
        break;
      }
    }else{
      c =       atoi(pm->junk);
      switch ( (int) c ) {
      case 1:
    select_coa( pm->coa );                /*  select what to analysis     */
    break;
      case 2:                             /* Num of axis to record        */
    printf ( "Changing the number of axis generated from %i " 
        "Please input new value [%i]", (int)pcoa->axis,(int)pcoa->axis);
    gets(pm->junk);
    if ( !strlen(pm->junk)   ) break;   
    if ( isalpha( (int) pm->junk[0])) break;
    i = (char)atoi(pm->junk);
    if ( pm->coa == 'a' && (i > 20 || i<0)  || ( i<0 || i>59 )) { 
      fprintf(pm->my_err,"Value is out of range adjusting to max value\n");
      if ( pm->coa == 'a' ) pcoa->axis = 20;
      else                  pcoa->axis = 59;
    } else {  
      pcoa->axis = i;
    } 
    break;
      case 3:                              /* Add additional genes          */
    printf("You have elected to add genes after the initial COA is complete\n"
           "these will not affect the generation of axis (factors) but can\n"
           "identify were these additional genes fall based on the trends \n"
           "identified among the original genes\n"
           "You must have a separate file containing sequence(s) that are\n"
           "to be added (these genes must be DNA in fasta format)\n"
           "Please input filename [cancel this option]: ");
    gets(pm->junk);
    if ( !strlen(pm->junk) ) break;
    strncpy(pcoa->add_row,pm->junk,MAX_FILENAME_LEN-1);
    break;
      case 4:                             /* report analysis of inertia     */
	pcoa->level = (char) ( (pcoa->level=='n')? 'e':'n'); 
	break;
      case 5:                             /* how to identify optimal codons */
        printf ("You have elected to alter the number of genes used \n"
                "to identify the optimal codons\n"
                "You can input either an absolute number of genes or a\n"
                "percentage (example 10%%)\n "  
                "\tPlease input your choice []");
    gets ( pm->junk);
    if( !strlen(pm->junk) ) continue;
    if( (p=strchr ( pm->junk,'%')) != NULL) {
          *p='\0';
      pcoa->fop_gene=atoi(pm->junk)*-1;
          if ( pcoa->fop_gene == 0 || pcoa->fop_gene < 50 ) { /* err_catch */
        printf ( " Limits are >0%% and less than 50%%\n");
        pcoa->fop_gene= (-10);                           /* assume default */
      }
    }else {
      pcoa->fop_gene=atoi(pm->junk);                      /* set No. genes */
        }
    break;
      default :
    fprintf(pm->my_err,"Answer out of range\n");
    break;  
      }              
    }    
  }
  return;
}

/*********************** select_coa ****************************************/
/* This menu is called if the user wants to change the default codons/AA   */
/* to be analysised in the COA. It is called from menu_coa                 */
/***************************************************************************/
void select_coa ( char choice ) 
{
 int   loop = TRUE;
 int   last_row[4];
 int   toggle;
 int   x;
 
 char  *startpoint, *endpoint;
 
 clearscr(pm->term_length);

 while ( loop ) { 
   if ( choice == 'a' ) {                   /* if AA analysis then         */
     for ( x = 1 ; x < 22 ; x++ ) {     
       if (!pcoa->amino[x] ) 
     printf("[(%2i)_%s_%s] ", x, paa->aa3[x],paa->aa1[x] );
       else
     printf(" (%2i)_%s_%s  ", x, paa->aa3[x],paa->aa1[x] );
       
       if ( !(x % 4) ) printf( "\n");
     }
     printf( "\n");

/*************** Sample of aa choice output    ****************************/
/* ( 1)_Phe_F   ( 2)_Leu_L   ( 3)_Ile_I   ( 4)_Met_M                      */
/* ( 5)_Val_V   ( 6)_Ser_S   ( 7)_Pro_P   ( 8)_Thr_T                      */
/* ( 9)_Ala_A   (10)_Tyr_Y  [(11)_TER_*]  (12)_His_H                      */
/* (13)_Gln_Q   (14)_Asn_N   (15)_Lys_K   (16)_Asp_D                      */
/* (17)_Glu_E   (18)_Cys_C   (19)_Trp_W   (20)_Arg_R                      */
/* (21)_Gly_G                                                             */

   }else {
     printf ( "Using %s \n", pcu->des ); 
     for ( x = 1 ; x < 65 ; x++ ) { 
       
       if ( !pcoa->codons[x] )          printf("[");
       else                         printf(" ");
       
       if (last_row[x%4] != pcu->ca[x] )
     printf( "(%2i) %s\t%s", x,paa->aa3[pcu->ca[x]],paa->cod[x]);
       else
     printf( "(%2i)    \t%s", x,paa->cod[x]);
       
       if ( !pcoa->codons[x] )         printf("]");
       else                            printf(" ");    
       
       last_row[x%4] = pcu->ca[x];
       
       if ( !(x % 4) ) 
     printf( "\n");
       if ( !(x % 16)) 
     printf( "\n");
     }
   }
      
/*************** Sample of codon choice output      ***********************/
/*   Using Universal Genetic code                                         */
/* ( 1) Phe       UUU  ( 2) Ser   UCU  ( 3) Tyr   UAU  ( 4) Cys   UGU     */
/* ( 5)           UUC  ( 6)       UCC  ( 7)       UAC  ( 8)       UGC     */
/* ( 9) Leu       UUA  (10)       UCA [(11) TER   UAA][(12) TER   UGA]    */
/* (13)           UUG  (14)       UCG [(15)       UAG][(16) Trp   UGG]    */

   printf("%s bracketed will be excluded from the COA. ", 
      (pm->coa == 'a')? "Amino Acids": "Codons" );
   printf("Select number(s) that\nidentify the %s you wish to toggle "
          "(X to exit, H for help) [X] ",
      (pm->coa == 'a')? "Amino Acids": "Codons" );
   
   gets(pm->junk);
   
   if ( !strlen(pm->junk) || toupper( (int) pm->junk[0]) == 'X' ) {
     loop=FALSE;  
     continue;
   }
    
   if ( toupper( (int) pm->junk[0]) == 'H' ) {
       chelp("select");
       continue;
       }


   endpoint   = pm->junk;
   startpoint = pm->junk;
   
   /* now toggle the codons and amino acids to be analysed                */

   while ( toggle = (int) strtol(startpoint,&endpoint,10) ) {
     if(endpoint == startpoint )    break;
     startpoint = endpoint;
     
     if (pm->coa == 'a' )  {
       if ( toggle>21 || toggle<1 ) continue;      /* check value is valid */    
       pcoa->amino [toggle]= (char)((pcoa->amino [toggle])?FALSE:TRUE);
     }else{
       if ( toggle>64 || toggle<1 ) continue;      /* check value is valid */       
       pcoa->codons[toggle]= (char)((pcoa->codons[toggle])?FALSE:TRUE);
     }
   } 
 }
 return;
}

/************************* Welcome *****************************************/
/* Prints a Banner                                                         */
/* the \'s are a problem as they must be escaped                           */
/***************************************************************************/
void    welcome ( void )
{
 printf ("\n\n");
 printf ("  //   \\   //    \\  |I    \\   //    \\  |I\\    I  /     \n");
 printf (" |I       |I      I |I     I |I      I |I\\\\   I  \\___  \n");
 printf (" |I       |I      I |I     I |I      I |I \\\\  I      \\ \n");
 printf (" |I       |I      I |I     I |I      I |I  \\\\ I       |\n");
 printf ("  \\\\___/   \\\\____/  |I____/   \\\\____/  |I   \\\\I  \\___/\n");
}

/********************** printinfo  *****************************************/
/* Prints a summary about this programme, date, version and author of code */
/* whether a debug version                                                 */
/***************************************************************************/
int  printinfo(void) {  
# if defined (__FILE__ )
  printf("\n\tSource   : %s", __FILE__);
# endif 
# if defined  (DEBUG)
  printf("(Debug version)");
# endif

  printf("\n\tAuthor   : John Peden\n");
  printf("\tVersion  : %.*s\n", strlen(Revision) , Revision ); 
  printf("\tRevised  :%.*s %s %.*s\n",(int) strlen(Update) - 7, Update + 6,
	 (*(Update + 7) ? "\n\t     by  :" : ""),
	 (int) strlen(Author) - 10, Author + 9);
  
#if defined(__DATE__ ) && defined(__TIME__)
  printf("\n\tCompiled : %s %s\n", __DATE__, __TIME__);
#endif
  
  printf("\n\t-------------------------------\n\n");
  
  printf(" All sequences must be in a single file separated by title "
      " lines whose\n first character is either ; or > \n\t any number"
      " or length of genes is acceptable\n\n");
  return 1;
}


