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
/******** Tester      *****************************************************/
/* This function is used to teach the genetic code, it generates a random */
/* series of questions about the selected genetic code.                   */
/* The questions include                                                  */
/*  1 and 3 letter amino acid names                                       */
/*  The translation of each codon                                         */
/*  The size of each amino acid family                                    */
/**************************************************************************/

#define rand_num(z) (int)((((float)rand()/((long)RAND_MAX))*(float)z)+1)

#ifdef _WINDOWS
#define beeep Beep(150,150)
#include <time.h>
#include <conio.h>
#else
#define beeep printf("\007")
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>  
#include <time.h>
#include <ctype.h>
#include "codonW.h"

/* The accuracy of the answers are recorded using these three variable     */
int num_questions = 0;                      
int num_cheats = 0;
int num_wrong = 0;


void    tester ( void ) {
    char    loop;
    char    main_loop=TRUE;
    char    tmp_AA [4];
    char    tmp_AA2[4];

    srand( (unsigned)time( NULL ) );          /* initialise random num gen */

    printf(" Welcome to TESTER \n(which just tests your "
	     "knowledge of the Genetic code)\n"
           " The genetic code used is dependant on\n what"
	     " code is selected in menu 3\n"
           " The current code is %s %s\n"
           "\n If you get stuck try typing ? for a hint\n"
           " To leave type exit or quit\n", pcu->des, pcu->typ);

    /*******************  main loop            ****************************/
    while ( main_loop )   {
      int i,x;

        i = rand_num(10);           /*  random number to between 1 and 10 */

        printf("Type Help for help:");
        /* the switch biases the questions so their freq is not equal     */
        switch (i) {
        case 1: 
        case 2:                     /*  amino acid question               */
            i = rand_num(21);
            loop = TRUE;
            while ( loop ) {
                printf("\nWhat is the three letter equivalent for the AA"
                    " %s ", paa->aa1[i]);
                gets( pm->junk ) ;
                strcpy ( tmp_AA, paa->aa3[i] );
                for ( x = 0 ; x < (int)strlen(tmp_AA); x++) 
                    tmp_AA[x] = (char) toupper( (int) tmp_AA[x]);
                for ( x = 0 ; x < (int)strlen(pm->junk  ); x++) 
                    pm->junk  [x] =  (char) toupper(  (int) pm->junk[x]);
                if ( !strcmp ( pm->junk, "QUIT" ) || 
                     !strcmp ( pm->junk, "EXIT" )) {
                    asummary();                 
                    main_loop = FALSE;
                    break;
                }

                if ( !strcmp ( pm->junk,"HELP")) {
                    chelp("fun");
                    continue;
                    }

                if ( !strcmp (pm->junk, "?" ) ) {
                    printf( "Cheat %s", paa->aa3[i]);
                    num_cheats++;             /*     The user cheated     */
                    continue; 
                }
                if ( !strcmp (pm->junk  , tmp_AA )) {
                    loop = FALSE;
                } else {
                    num_wrong++;              /*     Wrong answer       */
                    printf("Wrong answer (try ?)\n");
                }
            }
            break;
        case 3:                             /* How big is this AA family*/
            i = rand_num(21);
            loop = TRUE;
            while ( loop ) {
                printf("\nHow many codons encode the Amino Acid %s ",
                        paa->aa1[i]);
                gets( pm->junk ) ;
                for ( x = 0 ; x < (int)strlen(pm->junk); x++) 
                    pm->junk[x] = (char) toupper( (int) pm->junk[x]);
  
                if ( !strcmp ( pm->junk, "QUIT" ) || 
                     !strcmp ( pm->junk, "EXIT" )) {
                    asummary();
                    main_loop = FALSE;
                    break;
                }
  
                if ( !strcmp ( pm->junk,"HELP")) {
                    chelp("fun");
                    continue;
                    }

                if ( !strcmp (pm->junk, "?" ) ) {
                    printf( "Cheat %i\n", *(da + i) );
                    num_cheats++;
                    continue;

               }
                
               

                if ( atoi(pm->junk) == *(da + i) )
                    loop = FALSE;
                
                else {
                    num_wrong++;
                    printf("Wrong answer (try ?)\n");
                }
            }
            break;
        case 4:                                 /* 60% of the time ask    */
        case 5:                                 /* ask questions about    */     
        case 6:                                 /* codon to aa translation*/ 
        case 7: 
        case 8:
        case 9: 
        case 10:
            i = rand_num(64);
            loop = TRUE;
            while ( loop ) {
                printf("\nName the Amino Acid encoded by the codon %s ", paa->cod[i]);
                gets( pm->junk );
                for ( x = 0 ; x < (int)strlen(pm->junk ); x++) 
                    pm->junk[x] = (char) toupper( (int) pm->junk[x]);
                if ( !strcmp ( pm->junk, "QUIT" ) || 
                     !strcmp ( pm->junk, "EXIT" )) {
                    asummary();
                    main_loop = FALSE;
                    break;
                }

                
                if ( !strcmp ( pm->junk,"HELP")) {
                    chelp("fun");
                    continue;
                    }

                if ( !strcmp (pm->junk, "?" ) ) {
                    printf( "Cheat %s (%s)", paa->aa1[pcu->ca[i]]
                        , paa->aa3[pcu->ca[i]]);
                    num_cheats++;             /* tell me the answer      */
                    continue;
                }
                /* allow 1 or 3 letter amino acid code as the ans        */
                strcpy ( tmp_AA, paa->aa1[pcu->ca[i]] );
                strcpy ( tmp_AA2, paa->aa3[pcu->ca[i]] );

                /* uppercase everything, the AA names and the answer     */
                for ( x = 0 ; x < (int)strlen(tmp_AA); x++) 
                    tmp_AA[x] = (char)toupper( (int) tmp_AA[x]);
                for ( x = 0 ; x < (int)strlen(tmp_AA2); x++) 
                    tmp_AA2[x] = (char)toupper((int) tmp_AA2[x]);
                for ( x = 0 ; x < (int)strlen(pm->junk  ); x++) 
                    pm->junk  [x] = (char)toupper((int) pm->junk[x]);

                if ( !strcmp(tmp_AA, pm->junk) || 
                     !strcmp(tmp_AA2,pm->junk)  ) {         
                    loop = FALSE;
                } else {
                    printf("Wrong answer (try ?)\n");
                    num_wrong++;
                }
            }
            break;
        default:
            printf("mistake == %i \n", i);
            exit(0);                             /* error catch            */ 
            break;
        }                                        /* end of switch          */
        num_questions++;

    }                                            /* end of while           */

    return;
}                                                /* end of main            */

/*********** Asummary ******************************************************/
/* Write out a summary of the users results                                */
/***************************************************************************/
void    asummary (void) {
    printf ( " You answered\n \t %5i questions\n", num_questions);
    printf ( " \t %5i answers were wrong\n", num_wrong);
    printf ( " \t %5i times you had to ask for a hint\n", num_cheats);
    printf ( " \t  %3.0f%c accuracy \n", (float) ( (num_questions) ?                 
        (float)100 * (num_questions - num_wrong) / 
        (float)num_questions : 0 ),'%');
    pause;
    return;
}


