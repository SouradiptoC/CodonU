

CodonW is a package for codon usage analysis. It was designed to simplify
Multivariate Analysis (MVA) of codon usage. The MVA method employed in
CodonW is correspondence analysis (COA) (the most popular MVA method for
codon usage analysis). CodonW can gen erate a COA for codon usage,
relative synonymous codon usage or amino acid usage. Additional analyses
of codon usage include investigation of optimal codons, codon and
dinucleotide bias, and/or base composition. 

CodonW also has the capacity to analysis sequences encoded by genetic
codes other than the universal code.

Why call it codonW? 

Well first you must realise that "clustal" (a very popular multiple
alignment program by Des Higgins) was originally written in Paul's lab in
Trinity College Dublin. Clustal has since been rewritten from FORTRAN into
C and undergone several name changes c lustal-> clustalv-> clustalw ->
clustalx. There was also a program called "codons" written in FORTRAN by
Andrew Lloyd (a post-doc in Paul's lab), this was the original inspiration
for codonW. An early version of codonW, written in C, was called codonv.
Wh en the code was enhanced to include multivariate analysis, what better
name than codonW. 


CodonW version 1.3 June 1997 
================= 

The source code for CodonW can be obtained from
ftp://molbiol.ox.ac.uk/cu/codonW.tar.Z. Binaries for a number of platforms
are also available at this site see ftp://molbiol.ox.ac.uk/cu.


To Install and Build on UNIX Platforms 
================= 

Get the source code from ftp://molbiol.ox.ac.uk/cu/codonW.tar.Z Change
directory to the directory where you intend to install CodonW. 

uncompress codonW.tar.Z 
tar -xvf codonW.tar 
cd codonw 
./codonWinstall all     (this writes a makefiel and then builds codonw) 

This will ask a few questions regarding 'make' and 'cc' and then configure
the installation and compile the programs. If you don't understand the
questions, just accept the default by pressing the return key and the
installation should be OK using the defaults. The install script also
creates a number of links to the compiled executable codonW.  These links
allow codonW to emulate other useful codon usage analysis and sequence
manipulation software by passing the menu interface (for more informa tion
see README.links). Alternatively you can just elect to only build the main
program, and not install the linked programs. 

./codonWinstall codonw (compile only the executable codonw) 

Once you have successfully built codonw, try these commands to get you
started.  

./codonw -help (for commandline summary)  
./codonw        (menu interface) 

There is also a short tutorial. 


For the most recent documentation on codonW see
http://www.molbiol.ox.ac.uk/cu/


To Set the Codonw Help Environment:  
================= 

CodonW has an in-built help system, the help file is called codonW.hlp and
should be located in the same directory as the executable codonw.
Alternatively the help file can be pointed to by the environment variable
CODONW_H, if you are using a C shell you
 can add something similar to this to your .login script. 

setenv CODONW_H file_path 

Where file_path is the fully defined path name for codonW.hlp. 

Additional Files:
=================

README.indices - explanation about the various codon usage indices that
codonW calculates.  

README.coa- explanation about the output files from the correspondence
analysis. 

README.links- explanation about the auxiliary programmes created during
the making of codonw. 

Tutorial- A quick tutorial on the analysis of codon usage of the open
reading frames from Saccharomyces cerevisiae chromosome III. 

input.dat- An input file containing 167 open reading frames from
Saccharomyces cerevisiae chromosome III. (see Tutorial). 

Recoding - A quick explanation about how amino acids and codons have are
represented internally within codonW. 


Bugs 

This is a beta version of codonW, therefore there may be bugs within the
code. If you do find or notice anything strange please e-mail bug
reports/complaints/suggestions to johnp@molbiol.ox.ac.uk. Remember to
include an example of the input file (and outp ut files) and the options
selected that generated the error, don't forget to tell me the make of
computer and operating system it was running under.


