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

/* This is a general subroutine, so we might as well redefine TRUE & FALSE*/
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* What to do if we can't locate the file we where asked to open          */
/* On most systems we will try and be nice and show a choice of filenames */

#ifdef _DOS
#define no_file_found() system("dir/w");
#elif BSD || SYSV            
#define no_file_found() system("ls -F");
#elif defined (WIN32) || defined (_WIN) 
#define no_file_found() system("dir/w");
#else
#define no_file_found() printf("This would have presented a list of files\n\tbut I do not know howto your operating system\n");
#endif

/* Include header files                                                   */
#include <stdio.h>             
#include <limits.h> 
#include <stdlib.h>  
#include <string.h>
#include <ctype.h>         
#include "codonW.h"

/************** open_file **************************************************/
/* This subroutine is a front end to fopen open. It takes four parameters  */
/* the parameters are used to generate a user prompt for the               */
/* filename, and to give a suggested filename, to give the write perms     */
/* for the file, and whether or not to overwrite existing files.           */
/* File_needed is just a description of the file being opened. It is       */
/* assumed that if this descriptor is missing the file is to be opened     */
/* without further user input. If default_filename is blank then there is  */
/* no default_filename                                                     */
/* write_perm sets up the type of file being opened                        */
/* verbose tells this function whether to check if there is a         */
/* previous version of any file being opened for writing                   */
/***************************************************************************/

FILE *open_file(char *file_needed, char *default_filename, 
char *write_perm, int  verbose )
{
    char   infile_name[MAX_FILENAME_LEN]="";
    FILE  *input=NULL;
    char   temp[4];
    char  *answer = pm->junk;

    /**********************************************************************/
    /* If a string has been given for file_needed it is assumed           */
    /* that the user will have a choice of file_names to choose           */
    /* therefore (s)he will be prompted for a name                        */
    /* if a default filename was supplied by the calling function this    */
    /* will be suggested as well, otherwise there is no default           */
    /**********************************************************************/

    if ( strlen(file_needed)) {
        while (!strlen(infile_name) )  {
            printf("\nName of %s (h for help) [%s] ", 
                      file_needed,default_filename);
            gets(infile_name);                          /* get filename   */
            
            if ( WasHelpCalled ( infile_name ) ) {
                     chelp("open_file_query");          /* Help ....      */
                     infile_name[0]='\0';
                     continue;
                }

            if ( !strlen(infile_name) && default_filename )
                strcpy(infile_name, default_filename);
        }                                         /* end of get filename  */
    } else if ( strlen(default_filename) )        /* use default filename */
        strcpy(infile_name, default_filename);
    else {                                        /* not enough info      */
        fprintf(stderr, "Programming error: no filename supplied\n");
        my_exit (0,"open file");
    }


    /**********************************************************************/
    /* At this point infile_name contains a possible filename             */
    /* Depending on the mode (write_perm) string this is tested 3 ways    */
    /*                                                                    */
    /* (r or r+) Test if the file exists if not, all the files in the     */
    /* current directory are listed and the the user is prompted for      */
    /* an alternative name or they may quit the programme                 */
    /*                                                                    */
    /* (a, a+) Not tested, just open the file                             */
    /*                                                                    */
    /* (w, w+) If the variable verbose = FALSE then no test          */
    /* If verbose == TRUE then the file is checked to see if         */
    /* it already exsists, if it does then the user is prompted for       */
    /* either for permission to overwrite this file or to                 */
    /* suggest an alternative file_name which is then tested as well      */
    /* the user can type q to quit at any stage of this prompting process */
    /**********************************************************************/

    if ( !strcmp(write_perm, "r") || !strcmp(write_perm, "r+") 
       ||!strcmp(write_perm, "rb") ){
        while ( !(input = fopen (infile_name , write_perm )))  {
            fprintf(stderr,"\nThese are the files in the current directory "
                "I cannot find %.*s \n\n",strlen(infile_name),infile_name);
            no_file_found();
            fprintf(stderr, "\n\nPlease enter another filename, "
                " (Q)uit, (H)elp [%s] ",infile_name);
            gets(answer);
            
            if (strlen (answer)==1 && 
                   ((char)toupper((int)answer[0])=='Q'))
                   my_exit(2,"open_file");
            else if (WasHelpCalled ( infile_name )){
                   chelp ("File_not_found");
                }
            else if (strlen (answer))
                strcpy (infile_name, answer);  
		}                                     /* end of while loop */       
        strcpy ( answer,infile_name);           /* allow transfer    */
        return input;
    }                                               

    /************************* Append  ***********************************/
    else if ( !strcmp(write_perm, "a") || !strcmp(write_perm, "a+")
           || !strcmp(write_perm, "ab") ) {
        input = fopen (infile_name, write_perm);
        strcpy ( answer,infile_name);      
        return input;
    }                                              
    /************************* Write    **********************************/
    else if ( !strcmp( write_perm, "w") || !strcmp(write_perm, "w+") 
            ||!strcmp( write_perm, "wb") ) {

         while ( verbose == TRUE ) {
            if ( (input = fopen (infile_name , "r")) ) {
                fclose(input);                  /* close the filehandle  */
                fprintf(stderr, "\nWarning :File %.*s "
                    "exists already \n\tDo you wish to"
                    " overwrite ? (y/n/h/q)\t [y] ",
                    strlen(infile_name), infile_name);
                fgets(temp, 3, stdin);

                switch (toupper( (int) temp[0])) {
                case 'Y': 
                case '\0': 
                case '\n':
                    verbose = FALSE;
                    continue;
                case 'Q':
                    my_exit(2,"open_file2");
                    break;
                case 'H':
                    chelp("file_exists");
                    continue;
                    break;
                default:
                    fprintf(stderr, 
                        "\nYou decided not to overwrite, please enter\n"
                        " another filename, (q)uit, (a)ppend, (h)elp \n"
                        " (a/q/h/filename)\t[a] ");
                    gets(answer);
                }

                /* if the answer is 'a' then the default file is opened  */
                /* as appendable else if 'q' then the programme exits    */
                /* anything else is taken as a file name                 */

                if ( strlen(answer) <= 1 ) {
                    switch (toupper( (int) answer[0])) {
                    case 'Q':
                        return (NULL);
                    case 'A':
                    case '\0':
                    case'\n':
                        verbose = FALSE;  /* leave the while loop   */
                        strcpy(write_perm, "a+");
                        break;
                    case 'H':
                     chelp("file_append");        
                     continue;
                     break;
                    default:
                        continue;
                    };                                /* end of switch   */
                }               
            } else                              /* filename is unique    */
              verbose = FALSE;             /* exit the while loop   */
              }                                 /* match while preserve  */
        input = fopen (infile_name,write_perm); /* opens filehandle      */
        strcpy ( answer,infile_name);         
        return input;
    }                                           /* matchs if w or w+     */
    return (NULL);
}

/************** Main just for testing purposes ***************************/
/* uncomment to test function as a standalone subroutine                 */
/* will also need to replace my_exit with exit calls                     */
/*************************************************************************/
/*  main ()
 {
 FILE *test=NULL;
 if( test = open_file( "test file","","r",NULL))
    printf( "Success\n");
 else 
    printf( "Failed\n");
 } */
/*************************************************************************/



