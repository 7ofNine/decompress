/********************************************************************/
/*  Voyager Image Decompression Program - C Version for PC, VAX,    */
/*  UNIX and Macintosh systems.                                     */
/*                                                                  */
/*  Decompresses images using Kris Becker's subroutine DECOMP.C     */
/*  which is included in this program in a shortened version.       */
/*                                                                  */
/*  Reads a variable length compressed VOYAGER image and outputs a  */
/*  fixed length uncompressed image file in PDS format with         */
/*  labels, image histogram, engineering table and 800 lines of     */
/*  836 bytes (800 samples, 36 engineering bytes); or an 800 by     */
/*  800 array with FITS, VICAR or no labels.  If used on a non-     */
/*  byte-swapped machine the image histogram is un-swapped.         */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/*  Use the following commands to compile and link to produce an    */
/*  executable file:                                                */
/*                                                                  */
/*  On an IBM PC (using Microsoft C Version 5.x)                    */
/*                                                                  */
/*    cl /c cdcomp.c                                                */
/*    link  cdcomp/stack:10000;                                     */
/*                                                                  */
/*  On a VAX:                                                       */
/*                                                                  */
/*    cc   cdcomp                                                   */
/*    $define lnk$library sys$library:vaxcrtl.olb                   */
/*    link cdcomp                                                   */
/*                                                                  */
/*  On a Unix host (Sun, Masscomp)                                  */
/*                                                                  */
/*    cc -o cdcomp cdcomp.c                                         */
/*                                                                  */
/*  On a Macintosh (using Lightspeed C)                             */
/*                                                                  */
/*    link with the following libraries:                            */
/*    stdio, storage, strings, unix and MacTraps, with MacTraps     */
/*    in a separate segment.                                        */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/*  Use the following command to run the program:                   */
/*                                                                  */
/*    CDCOMP [infile] [outfile] [output format]                     */
/*                                                                  */
/*       infile        - name of compressed image file.             */
/*       outfile       - name of uncompressed output file.          */
/*       output format - selected from the following list:          */
/*                                                                  */
/*          1  SFDU/PDS format [DEFAULT].                           */
/*          2  FITS format.                                         */
/*          3  VICAR format.                                        */
/*          4  Unlabelled binary array.                             */
/*                                                                  */
/*  On the VAX computer you will need to 'install' the program to   */
/*  be able to use command line arguments using the following       */
/*  command:                                                        */
/*                                                                  */
/*  $ vaxdcomp :== $DISKNAME:[DIRECTORY]vaxcomp.exe                 */
/*                                                                  */
/*  where DISKNAME is the disk drive and DIRECTORY is the           */
/*  directory where VAXDCOMP.EXE is stored.                         */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/* LIMS                                                             */
/*  This program has been tested on a VAX 780 (VMS 4.6), SUN        */
/*  Workstation (UNIX 4.2, release 3.4), an IBM PC                  */
/*  (MICROSOFT 5.1 compiler) and Macintosh IIx using Lightspeed C   */
/*  version 3.0.  When converting to other systems, check for       */
/*  portability conflicts.                                          */
/*                                                                  */
/* HIST                                                             */
/*  AUG89 Added code to get command line arguments for filenames    */
/*  and output format; routines to free memory used by the Huffman  */
/*  tree); fixed the SFDU label output length; and modified the     */
/*  I/O routines so that the open for Host type 2 uses binary I/O.  */
/*  JUN89 Fixed READVAR, to get length on 16-bit unswapped hosts.   */
/*  JUL88 C driver to decompress standard Voyager Compressed images */
/*  by Mike Martin 1989/06/10                                       */
/*                                                                  */
/*  Inputs   - Input file to be decompressed.                       */
/*                                                                  */
/*  Outputs  - Output file containing decompressed image.           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>

#define TRUE                  1
#define FALSE                 0
#define RECORD_BYTES        836

                                    /* pc i/o defines               */
#define O_RDONLY         0x0000     /* open for reading only        */
#define O_BINARY         0x8000     /* file mode is binary          */

                                    /* vax i/o defines              */
#define RECORD_TYPE      "rfm=fix"  /* VAX fixed length output      */
#define CTX              "ctx=bin"  /* no translation of \n         */
#define FOP          "fop=cif,sup"  /* file processing ops          */

typedef struct leaf
  {
   struct leaf *right;
   short int dn;
   struct leaf *left;
  } NODE;

/*************************************************************************
 Declare the tree pointer. This pointer will hold the root of the tree
 once the tree is created by the accompanying routine huff_tree.
**************************************************************************/

  NODE *tree;

/* subroutine definitions                                           */

void               pds_labels();
void               fits_labels();
void               vicar_labels();
void               no_labels();
int                check_host();
int                get_files();
long               swap_long();
void               decompress();
void               decmpinit();
void               free_tree();

/* global variables                                                 */

int                infile;
FILE               *outfile;
char               inname[80],outname[80];
int                output_format;

main(argc,argv)

int  argc;
char **argv;

{
char          ibuf[2048],obuf[2048];
unsigned char blank=32;
short         host,length,total_bytes,line,i;
long          long_length,hist[511];
int           out_bytes = RECORD_BYTES;
int           count;

/*********************************************************************/
/*                                                                   */
/* get host information and input and output files                   */
/*                                                                   */
/*********************************************************************/


   strcpy(inname,"   ");  
   strcpy(outname,"   ");
   output_format = 0;

   if (argc == 1);                     /* prompt user for parameters */
   else if (argc == 2 && (strncmp(argv[1],"help",4) == 0 || 
                          strncmp(argv[1],"HELP",4) == 0 ||
                          strncmp(argv[1],"?",1)    == 0))
     {
      printf("Voyager Image Decompression Program.  Command line format:\n\n");
      printf("CDCOMP [infile] [outfile] [format code]\n");
      printf("   infile        - name of compressed image file.     \n");
      printf("   outfile       - name of uncompressed output file.  \n");
      printf("   output format - selected from the following list:  \n");   
      printf("                                                      \n");   
      printf("     1  SFDU/PDS format [DEFAULT].                    \n");   
      printf("     2  FITS format.                                  \n");   
      printf("     3  VICAR format.                                 \n");   
      printf("     4  Unlabelled binary array.                    \n\n");   
      exit(1);
     }  
   else 
     {
      strcpy(inname,argv[1]);  
      if (argc >= 3) strcpy(outname,argv[2]); 
      if (argc == 3) output_format = 1;
      if (argc == 4) sscanf(argv[3],"%d",&output_format); 
     }
   host = check_host();
   host = get_files(host); /* may change host if VAX */

/*********************************************************************/
/*                                                                   */
/* read and edit compressed file labels                              */
/*                                                                   */
/*********************************************************************/

   switch (output_format)
     {
       case 1: pds_labels(host);
               break;
       case 2: fits_labels(host);
               break;
       case 3: vicar_labels(host);
               break;
       case 4: no_labels(host);
     }
/*********************************************************************/
/*                                                                   */
/* process the image histogram                                       */
/*                                                                   */
/*********************************************************************/

/* need to know record_bytes,hist_count,hist_item_type,item_count.*/
   total_bytes = 0;
   length = read_var((char *)hist,host);
   total_bytes = total_bytes + length;
   length = read_var((char *)hist+836,host);
   total_bytes = total_bytes + length;

   if (host == 2 || host == 5)             /* If non-byte swapped    */
     for (i=0;i<256;i++)                   /* host, swap bytes in    */
       hist[i] = swap_long(hist[i]);       /* the output histogram   */

   if (output_format == 1)
     {
      fwrite((char *)hist,    836,1,outfile);
      fwrite((char *)hist+836,length,1,outfile);

      /*  pad out the histogram to a multiple of RECORD_BYTES */
      for (i=total_bytes;i<RECORD_BYTES*2;i++) fputc(blank,outfile);
     }
/*********************************************************************/
/*                                                                   */
/* process the encoding histogram                                    */
/* don't have to byte-swap because DECOMP.C does it for us           */
/*                                                                   */
/*********************************************************************/

   length = read_var((char *)hist,host);
   length = read_var((char *)hist+836,host);
   length = read_var((char *)hist+1672,host);

/*********************************************************************/
/*                                                                   */
/* process the engineering summary                                   */
/*                                                                   */
/*********************************************************************/

   total_bytes = 0;
   length = read_var(ibuf,host);

   if (output_format == 1)
     {
      fwrite(ibuf,length,1,outfile);
      total_bytes = total_bytes + length;

      /*  pad out engineering to multiple of 836 */
      for (i=total_bytes;i<RECORD_BYTES;i++) fputc(blank,outfile);
     }
/*********************************************************************/
/*                                                                   */
/* initialize the decompression                                      */
/*                                                                   */
/*********************************************************************/

    printf("\nInitializing decompression routine...");
	decmpinit(hist);

/*********************************************************************/
/*                                                                   */
/* decompress the image                                              */
/*                                                                   */
/*********************************************************************/

	printf("\nDecompressing data...");
    line=0;
    do
      {
       length = read_var(ibuf,host);
       if (length <= 0) break;
       long_length = (long)length;
       line += 1;
       decompress(ibuf, obuf,&long_length, &out_bytes);
       if (output_format == 1)
         {
          count = fwrite(obuf,RECORD_BYTES,1,outfile);
          if (count != 1)
            {
             printf("\nError writing output file.  Aborting program.");
             printf("\nCheck disk space or for duplicate file name on VAX.");
             exit(1);
            }
         }
       else fwrite(obuf,800,1,outfile);

       if (line % 100 == 0) printf("\n line %d",line);
      } while (length > 0 && line < 800);

 /*  pad out FITS file to a multiple of 2880 */
 if (output_format == 2)
   for (i=0;i<2240;i++) fputc(blank,outfile);

 printf("\n");
 free_tree(&long_length);
 close(infile);
 fclose(outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine get_files - get input filenames and open               */
/*                                                                   */
/*********************************************************************/

int get_files(host)
int host;

{
short   shortint;

  if (inname[0] == ' ')
    {
     printf("\nEnter name of file to be decompressed: ");
     gets (inname);
    }
  if (host == 1 | host == 2)
    {
     if ((infile = open(inname,O_RDONLY | O_BINARY)) <= 0)
       {
        printf("\ncan't open input file: %s\n",inname);
        exit(1);
       }
    }
  else if (host == 3 | host == 5)
    {
     if ((infile = open(inname,O_RDONLY)) <= 0)
       {
        printf("\ncan't open input file: %s\n",inname);
        exit(1);
       }

  /****************************************************************/
  /* If we are on a vax see if the file is in var length format.  */
  /* This logic is in here in case the vax file has been stored   */
  /* in fixed or undefined format.  This might be necessary since */
  /* vax variable length files can't be moved to other computer   */
  /* systems with standard comm programs (kermit, for example).   */
  /****************************************************************/

     if (host == 3)
       {
        read(infile,&shortint,2);
        if (shortint > 0 && shortint < 80)
          {
           host = 4;              /* change host to 4                */
           printf("This is not a VAX variable length file.");
          }
        else printf("This is a VAX variable length file.");
        lseek(infile,0,0);        /* reposition to beginning of file */
       }
    }

  if (output_format == 0)
  do
    {
     printf("\nEnter a number for the output format desired:\n");
     printf("\n  1.  SFDU/PDS format.");
     printf("\n  2.  FITS format.");
     printf("\n  3.  VICAR format.");
     printf("\n  4.  Unlabelled binary array.\n");
     printf("\n  Enter format number:");
     gets(inname);
     output_format = atoi(inname);
    } while (output_format < 1 || output_format > 4);

  if (outname[0] == ' ')
    {
     printf("\nEnter name of uncompressed output file: ");
     gets (outname);
    }

  if (host == 1 || host == 2 || host == 5)
    {
     if ((outfile = fopen(outname,"wb"))==NULL)
       {
        printf("\ncan't open output file: %s\n",outname);
        exit(1);
       }
    }
  else if (host == 3 || host == 4)
    {
     if (output_format == 1)      /* write PDS format blocks */
       {
        if ((outfile=fopen(outname,"w",
                           "mrs=836",FOP,CTX,RECORD_TYPE))==NULL)
         {
          printf("\ncan't open output file: %s\n",outname);
          exit(1);
         }
       }
     else if (output_format == 2) /* write FITS format blocks */
       {
        if ((outfile=fopen(outname,"w",
                           "mrs=2880",FOP,CTX,RECORD_TYPE))==NULL)
         {
          printf("\ncan't open output file: %s\n",outname);
          exit(1);
         }
       }
     else                         /* write 800 byte records */
       {
        if ((outfile=fopen(outname,"w",
                           "mrs=800",FOP,CTX,RECORD_TYPE))==NULL)
         {
          printf("\ncan't open output file: %s\n",outname);
          exit(1);
         }
       }
    }

  return(host);  /* In case its been updated */
}

/*********************************************************************/
/*                                                                   */
/* subroutine pds_labels - edit PDS labels and write to output file  */
/*                                                                   */
/*********************************************************************/

void pds_labels(host)
int host;

{
char          outstring[80],ibuf[2048];
unsigned char cr=13,lf=10,blank=32;
short         length,nlen,total_bytes,line,i;


total_bytes = 0;
do
  {
   length = read_var(ibuf,host);
   ibuf[length]=NULL;

  /******************************************************************/
  /*  edit labels which need to be changed                          */
  /******************************************************************/

   if      ((i = strncmp(ibuf,"NJPL1I00PDS1",12)) == 0)
   /*****************************************************************/
   /* add the output file length to the sfdu label                  */
   /*****************************************************************/
     {
      strcpy(outstring,ibuf);
      strcpy(outstring+12,"00673796");
      strcpy(outstring+20,ibuf+20);
      fwrite(outstring,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"RECORD_TYPE",11)) == 0)
   /*****************************************************************/
   /* change the record_type value from variable to fixed           */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"FIXED_LENGTH");
      length = length - 3;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"FILE_RECORDS",12)) == 0)
   /*****************************************************************/
   /* change the file_records count to 806                          */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"806");
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"LABEL_RECORDS",13)) == 0)
   /*****************************************************************/
  /* change the label_records count from 56 to 3                    */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"3");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"^IMAGE_HISTOGRAM",16)) == 0)
   /*****************************************************************/
   /* change the location pointer of image_histogram to record 4    */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"4");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"^ENCODING_HISTOGRAM)",19)) == 0);
   /*****************************************************************/
   /* delete the encoding_histogram location pointer                */
   /*****************************************************************/
   else if ((i = strncmp(ibuf,"^ENGINEERING_TABLE",18)) == 0)
   /*****************************************************************/
   /* change the location pointer of engineering_summary to record 6*/
   /*****************************************************************/
     {
      strcpy(ibuf+35,"6");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"^IMAGE",6)) == 0)
   /*****************************************************************/
   /* change the location pointer of image to record 7              */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"7");
      length = length -1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,
             "OBJECT                           = ENCODING",43)) == 0)
   /*****************************************************************/
   /* delete the 4 encoding histogram labels                        */
   /*****************************************************************/
     {
      for (i=0;i<4;i++)   /* ignore these labels */
        {
         length = read_var(ibuf,host);
        }
     }
   else if ((i = strncmp(ibuf," ENCODING",9)) == 0);
   /*****************************************************************/
   /* delete the encoding type label in the image object            */
   /*****************************************************************/
   else if ((host == 2 || host == 5) && (i = strncmp(ibuf,
             " ITEM_TYPE                       = VAX_INTEGER",46)) == 0)
   /*****************************************************************/
   /* change the record_type value from variable to fixed           */
   /*****************************************************************/
     {
      strcpy(ibuf+35,"INTEGER");
      length = length - 4;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }



   /*****************************************************************/
   /* if none of above write out the label to the output file       */
   /*****************************************************************/
   else
     {
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"%c%c",cr,lf);
      total_bytes = total_bytes + length + 2;
     }
   /*****************************************************************/
   /* test for the end of the PDS labels                            */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (i=total_bytes;i<RECORD_BYTES*3;i++) fputc(blank,outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine fits_labels - write FITS header to output file */
/*                                                                   */
/*********************************************************************/

void fits_labels(host)
int host;

{
char          ibuf[2048],outstring[80];
unsigned char cr=13,lf=10,blank=32;
short         length,nlen,total_bytes,line,i;

do
  {
   length = read_var(ibuf,host);
   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

total_bytes = 0;

strcpy(outstring,
"SIMPLE  =                    T                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

strcpy(outstring,
"BITPIX  =                    8                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

strcpy(outstring,
"NAXIS   =                    2                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

strcpy(outstring,
"NAXIS1  =                  800                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

strcpy(outstring,
"NAXIS2  =                  800                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

strcpy(outstring,
"END                                                                           ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 80;

/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (i=total_bytes;i<2880;i++) fputc(blank,outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine vicar_labels - write vicar labels to output file       */
/*                                                                   */
/*********************************************************************/

void vicar_labels(host)
int host;

{
char          ibuf[2048],outstring[80];
unsigned char cr=13,lf=10,blank=32;
short         length,nlen,total_bytes,line,i;

do
  {
   length = read_var(ibuf,host);
   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

total_bytes = 0;

strcpy(outstring,
"LBLSIZE=800             FORMAT='BYTE'  TYPE='IMAGE'  BUFSIZ=800  DIM=2  ");
fwrite(outstring,72,1,outfile);
total_bytes = total_bytes + 72;
strcpy(outstring,
"EOL=0  RECSIZE=800  ORG='BSQ'  NL=800  NS=800  NB=1  N1=0  N2=0  N3=0  ");
total_bytes = total_bytes + 71;
fwrite(outstring,71,1,outfile);
strcpy(outstring,
"N4=0  NBB=0  NLB=0");
fwrite(outstring,18,1,outfile);
fprintf(outfile,"%c%c",cr,lf);
total_bytes = total_bytes + 20;

/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (i=total_bytes;i<800;i++) fputc(blank,outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine no_labels - read past the pds labels                   */
/*                                                                   */
/*********************************************************************/

void no_labels(host)
int host;

{
char          ibuf[2048],outstring[80];
unsigned char cr=13,lf=10,blank=32;
short         length,nlen,total_bytes,line,i;

do
  {
   length = read_var(ibuf,host);
   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

}


/*********************************************************************/
/*                                                                   */
/* subroutine read_var - read variable length records from input file*/
/*                                                                   */
/*********************************************************************/

read_var(ibuf,host)
char  *ibuf;
int   host;
{
int   length,result,nlen;
char  temp;
union /* this union is used to swap 16 and 32 bit integers          */
  {
   char  ichar[4];
   short slen;
   long  llen;
  } onion;

  switch (host)
    {
     case 1: /*******************************************************/
             /* IBM PC host                                         */
             /*******************************************************/
       length = 0;
       result = read(infile,&length,2);
       nlen =   read(infile,ibuf,length+(length%2));
       return (length);
       break;

     case 2: /*******************************************************/
             /* Macintosh host                                      */
             /*******************************************************/
       length = 0;
       result = read(infile,onion.ichar,2);
       /*     byte swap the length field                            */
       temp   = onion.ichar[0];
       onion.ichar[0]=onion.ichar[1];
       onion.ichar[1]=temp;
       length = onion.slen;       /* left out of earlier versions   */
       /* printf("length=%04x,result=%d\n",length,result);          */
       nlen =   read(infile,ibuf,length+(length%2));
       return (length);
       break;

     case 3: /*******************************************************/
             /* VAX host with variable length support               */
             /*******************************************************/
       length = read(infile,ibuf,RECORD_BYTES);
       return (length);

     case 4: /*******************************************************/
             /* VAX host, but not a variable length file            */
             /*******************************************************/
       length = 0;
       result = read(infile,&length,2);
       nlen =   read(infile,ibuf,length+(length%2));

       /* check to see if we crossed a vax record boundary          */
       while (nlen < length)
         nlen += read(infile,ibuf+nlen,length+(length%2)-nlen);
       return (length);
       break;

     case 5: /*******************************************************/
             /* Unix workstation host (non-byte-swapped 32 bit host)*/
             /*******************************************************/
       length = 0;
       result = read(infile,onion.ichar,2);
       /*     byte swap the length field                            */
       temp   = onion.ichar[0];
       onion.ichar[0]=onion.ichar[1];
       onion.ichar[1]=temp;
       length = onion.slen;
       /* printf("length=%04x,result=%d\n",length,result); */
       nlen =   read(infile,ibuf,length+(length%2));
       return (length);
       break;
    }
}

/*********************************************************************/
/*                                                                   */
/* subroutine check_host - find out what kind of machine we are on   */
/*                                                                   */
/*********************************************************************/

int check_host()
{
/*  This subroutine checks the attributes of the host computer and
    returns a host code number.
*/
char hostname[80];

int swap,host,bits,var;
union
  {
   char  ichar[2];
   short ilen;
  } onion;

if (sizeof(var) == 4) bits = 32;
                 else bits = 16;

onion.ichar[0] = 1;
onion.ichar[1] = 0;

if (onion.ilen == 1) swap = TRUE;
else                 swap = FALSE;

if (bits == 16 && swap == TRUE)
  {
   host = 1; /* IBM PC host  */
   strcpy(hostname,
          "Host 1 - 16 bit integers with swapping, no var len support.");
  }

if (bits == 16 && swap == FALSE)
  {
   host = 2; /* Non byte swapped 16 bit host  */
   strcpy(hostname,
          "Host 2 - 16 bit integers without swapping, no var len support.");
  }

if (bits == 32 && swap == TRUE)
 { host = 3; /* VAX host with var length support */
   strcpy(hostname,
          "Host 3,4 - 32 bit integers with swapping.");
 }

if (bits == 32 && swap == FALSE)
  {
   host = 5; /* OTHER 32-bit host  */
   strcpy(hostname,
          "Host 5 - 32 bit integers without swapping, no var len support.");
  }

printf("%s\n",hostname);
return(host);
}

long swap_long(inval)  /* swap 4 byte integer                       */
long inval;
{
union /* this union is used to swap 16 and 32 bit integers          */
  {
   char  ichar[4];
   short slen;
   long  llen;
  } onion;
  char   temp;

  /* byte swap the input field                                      */
  onion.llen   = inval;
  temp   = onion.ichar[0];
  onion.ichar[0]=onion.ichar[3];
  onion.ichar[3]=temp;
  temp   = onion.ichar[1];
  onion.ichar[1]=onion.ichar[2];
  onion.ichar[2]=temp;
  return (onion.llen);
}

 void decompress(ibuf,obuf,nin,nout)
/****************************************************************************
*_TITLE decompress - decompresses image lines stored in compressed format   *
*_ARGS  TYPE       NAME      I/O        DESCRIPTION                         */
        char       *ibuf;  /* I         Compressed data buffer              */
        char       *obuf;  /* O         Decompressed image line             */
        long int   *nin;   /* I         Number of bytes on input buffer     */
        long int   *nout;  /* I         Number of bytes in output buffer    */

  {
 /* The external root pointer to tree */
    extern NODE *tree;

 /* Declare functions called from this routine */
    void dcmprs();

/*************************************************************************
  This routine is fairly simple as it's only function is to call the
  routine dcmprs.
**************************************************************************/

    dcmprs(ibuf,obuf,nin,nout,tree);

    return;
  }

void decmpinit(hist)
/***************************************************************************
*_TITLE decmpinit - initializes the Huffman tree                           *
*_ARGS  TYPE       NAME      I/O        DESCRIPTION                        */
        long int   *hist;  /* I         First-difference histogram.        */

{
  extern NODE *tree;          /* Huffman tree root pointer */

  /* Specify the calling function to initialize the tree */
  NODE *huff_tree();

/****************************************************************************
  Simply call the huff_tree routine and return.
*****************************************************************************/

  tree = huff_tree(hist);

  return;
 }

NODE *huff_tree(hist)
/****************************************************************************
*_TITLE huff_tree - constructs the Huffman tree; returns pointer to root    *
*_ARGS  TYPE          NAME        I/O   DESCRIPTION                         */
        long int     *hist;     /* I    First difference histogram          */

  {
  /*  Local variables used */
    long int freq_list[512];      /* Histogram frequency list */
    NODE **node_list;             /* DN pointer array list */

    register long int *fp;        /* Frequency list pointer */
    register NODE **np;           /* Node list pointer */

    register long int num_freq;   /* Number non-zero frequencies in histogram */
    long int sum;                 /* Sum of all frequencies */

    register short int num_nodes; /* Counter for DN initialization */
    register short int cnt;       /* Miscellaneous counter */

    short int znull = -1;         /* Null node value */

    register NODE *temp;          /* Temporary node pointer */

  /* Functions called */
    void sort_freq();
    NODE *new_node();
    char *malloc();

/***************************************************************************
  Allocate the array of nodes from memory and initialize these with numbers
  corresponding with the frequency list.  There are only 511 possible
  permutations of first difference histograms.  There are 512 allocated
  here to adhere to the FORTRAN version.
****************************************************************************/

   fp = freq_list;
   node_list = (NODE **) malloc(sizeof(temp)*512);
   if (node_list == NULL)
    {
      printf("\nOut of memory in huff_tree!\n");
      exit(1);
    }
   np = node_list;

   for (num_nodes=1, cnt=512 ; cnt-- ; num_nodes++)
     {
/**************************************************************************
    The following code has been added to standardize the VAX byte order
    for the "long int" type.  This code is intended to make the routine
    as machine independant as possible.
***************************************************************************/
        unsigned char *cp = (unsigned char *) hist++;
        unsigned long int j;
        short int i;
        for (i=4 ; --i >= 0 ; j = (j << 8) | *(cp+i));

/* Now make the assignment */
        *fp++ = j;
        temp = new_node(num_nodes);
        *np++ = temp;
     }

     (*--fp) = 0;         /* Ensure the last element is zeroed out.  */

/***************************************************************************
  Now, sort the frequency list and eliminate all frequencies of zero.
****************************************************************************/

  num_freq = 512;
  sort_freq(freq_list,node_list,num_freq);

  fp = freq_list;
  np = node_list;

  for (num_freq=512 ; (*fp) == 0 && (num_freq) ; fp++, np++, num_freq--);


/***************************************************************************
  Now create the tree.  Note that if there is only one difference value,
  it is returned as the root.  On each interation, a new node is created
  and the least frequently occurring difference is assigned to the right
  pointer and the next least frequency to the left pointer.  The node
  assigned to the left pointer now becomes the combination of the two
  nodes and it's frequency is the sum of the two combining nodes.
****************************************************************************/

  for (temp=(*np) ; (num_freq--) > 1 ; )
    {
        temp = new_node(znull);
        temp->right = (*np++);
        temp->left = (*np);
        *np = temp;
        *(fp+1) = *(fp+1) + *fp;
        *fp++ = 0;
        sort_freq(fp,np,num_freq);
    }

  return temp;
 }

NODE *new_node(value)
/****************************************************************************
*_TITLE new_node - allocates a NODE structure and returns a pointer to it   *
*_ARGS  TYPE        NAME        I/O     DESCRIPTION                         */
        short int   value;    /* I      Value to assign to DN field         */

  {
    NODE *temp;         /* Pointer to the memory block */

  char *malloc();       /* Memory allocation function */

/***************************************************************************
  Allocate the memory and intialize the fields.
****************************************************************************/

  temp = (NODE *) malloc(sizeof(NODE));

  if (temp != NULL)
    {
      temp->right = NULL;
      temp->dn = value;
      temp->left = NULL;
    }
  else
    {
       printf("\nOut of memory in new_node!\n");
       exit(1);
    }

   return temp;
  }

 void sort_freq(freq_list,node_list,num_freq)
/****************************************************************************
*_TITLE sort_freq - sorts frequency and node lists in increasing freq. order*
*_ARGS  TYPE       NAME            I/O  DESCRIPTION                         */
        long int   *freq_list;   /* I   Pointer to frequency list           */
        NODE       **node_list;  /* I   Pointer to array of node pointers   */
        long int   num_freq;     /* I   Number of values in freq list       */

  {
    /* Local Variables */
    register long int *i;       /* primary pointer into freq_list */
    register long int *j;       /* secondary pointer into freq_list */

    register NODE **k;          /* primary pointer to node_list */
    register NODE **l;          /* secondary pointer into node_list */

    long int temp1;             /* temporary storage for freq_list */
    NODE *temp2;                /* temporary storage for node_list */

    register long int cnt;      /* count of list elements */


/************************************************************************
  Save the current element - starting with the second - in temporary
  storage.  Compare with all elements in first part of list moving
  each up one element until the element is larger.  Insert current
  element at this point in list.
*************************************************************************/

   if (num_freq <= 0) return;      /* If no elements or invalid, return */

   for (i=freq_list, k=node_list, cnt=num_freq ; --cnt ; *j=temp1, *l=temp2)
     {
        temp1 = *(++i);
        temp2 = *(++k);

        for (j = i, l = k ;  *(j-1) > temp1 ; )
          {
            *j = *(j-1);
            *l = *(l-1);
            j--;
            l--;
            if ( j <= freq_list) break;
          }

     }
  return;
  }

 void dcmprs(ibuf,obuf,nin,nout,root)
/****************************************************************************
*_TITLE dcmprs - decompresses Huffman coded compressed image lines          *
*_ARGS  TYPE       NAME       I/O       DESCRIPTION                         */
        char       *ibuf;   /* I        Compressed data buffer              */
        char       *obuf;   /* O        Decompressed image line             */
        long int   *nin;    /* I        Number of bytes on input buffer     */
        long int   *nout;   /* I        Number of bytes in output buffer    */
        NODE       *root;   /* I        Huffman coded tree                  */

  {
    /* Local Variables */
    register NODE *ptr = root;        /* pointer to position in tree */
    register unsigned char test;      /* test byte for bit set */
    register unsigned char idn;       /* input compressed byte */

    register char odn;                /* last dn value decompressed */

    char *ilim = ibuf + *nin;         /* end of compressed bytes */
    char *olim = obuf + *nout;        /* end of output buffer */

/**************************************************************************
  Check for valid input values for nin, nout and make initial assignments.
***************************************************************************/

    if (ilim > ibuf && olim > obuf)
       odn = *obuf++ = *ibuf++;
    else
       {
           printf("\nInvalid byte count in dcmprs!\n");
           exit(1);
       }

/**************************************************************************
  Decompress the input buffer.  Assign the first byte to the working
  variable, idn.  An arithmatic and (&) is performed using the variable
  'test' that is bit shifted to the right.  If the result is 0, then
  go to right else go to left.
***************************************************************************/

    for (idn=(*ibuf) ; ibuf < ilim  ; idn =(*++ibuf))
     {
        for (test=0x80 ; test ; test >>= 1)
           {
            ptr = (test & idn) ? ptr->left : ptr->right;

            if (ptr->dn != -1)
              {
                if (obuf >= olim) return;
                odn -= ptr->dn + 256;
                *obuf++ = odn;
                ptr = root;
              }
          }
     }
   return;
  }


void free_tree(nfreed)
/****************************************************************************
*_TITLE free_tree - free memory of all allocated nodes                      *
*_ARGS  TYPE       NAME       I/O        DESCRIPTION                        */
        long int   *nfreed;  /* O        Return of total count of nodes     *
*                                        freed.                             */

/*
*_DESCR This routine is supplied to the programmer to free up all the       *
*       allocated memory required to build the huffman tree.  The count     *
*       of the nodes freed is returned in the parameter 'nfreed'.  The      *
*       purpose of the routine is so if the user wishes to decompress more  *
*       than one file per run, the program will not keep allocating new     *
*       memory without first deallocating all previous nodes associated     *
*       with the previous file decompression.                               *

*_HIST  16-AUG-89 Kris Becker   USGS, Flagstaff Original Version            *
*_END                                                                       *
****************************************************************************/

{
	long int total_free = 0;

	extern NODE *tree;      /* Huffman tree root pointer */

/* Specify the function to free the tree */
	long int free_node();
/****************************************************************************
  Simply call the free_node routine and return the result.
*****************************************************************************/

	*nfreed = free_node(tree,total_free);

	return;
}

long int free_node(pnode,total_free)
/***************************************************************************
*_TITLE free_node - deallocates an allocated NODE pointer
*_ARGS  TYPE     NAME          I/O   DESCRIPTION                           */
        NODE     *pnode;       /* I  Pointer to node to free               */
        long int total_free;   /* I  Total number of freed nodes           */

/*
*_DESCR  free_node will check both right and left pointers of a node       *
*        and then free the current node using the free() C utility.        *
*        Note that all nodes attached to the node via right or left        *
*        pointers area also freed, so be sure that this is the desired     *
*        result when calling this routine.                                 *

*        This routine is supplied to allow successive calls to the         *
*        decmpinit routine.  It will free up the memory allocated          *
*        by previous calls to the decmpinit routine.  The call to free     *
*        a previous huffman tree is:  total = free_node(tree,(long) 0);    *
*        This call must be done by the programmer application routine      *
*        and is not done by any of these routines.                         *
*_HIST   16-AUG-89  Kris Becker U.S.G.S  Flagstaff Original Version        */
{
	if (pnode == (NODE *) NULL) return(total_free);
	
	if (pnode->right != (NODE *) NULL)
		total_free = free_node(pnode->right,total_free);
	if (pnode->left != (NODE *) NULL)
		total_free = free_node(pnode->left,total_free);

	free((char *) pnode);
	return(total_free + 1);
}
