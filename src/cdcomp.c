/********************************************************************/
/*  Voyager Image Decompression Program - C Version for PC          */
/*                                                                  */
/*  Decompresses images using Kris Becker's subroutine DECOMP.C     */
/*  which is included in this program in a shortened version.       */
/*                                                                  */
/*  Reads a variable length compressed VOYAGER image and outputs a  */
/*  fixed length uncompressed image file in PDS format with         */
/*  labels, image histogram, engineering table and 800 lines of     */
/*  836 bytes (800 samples, 36 engineering bytes); or an 800 by     */
/*  800 array with FITS, VICAR or no labels.                        */
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <io.h>
#include <fcntl.h>
#include <share.h>
#include <sys/stat.h>

#define RECORD_BYTES        836

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

static NODE *tree = NULL;

/* subroutine definitions                                           */

void               pds_labels();
void               fits_labels();
void               vicar_labels();
void               no_labels();
void               get_files();
void               dcmprs(char*, char*, long int*, long int*, NODE*);
void               decmpinit();
void               free_tree(long int *);
int                read_var(char*);


/* global variables     */

#define NAME_SIZE 80
#define BUFF_SIZE 2048

#define FORMAT_UNDEFINED 0
#define FORMAT_PDS 1
#define FORMAT_FITS 2
#define FORMAT_VICAR 3
#define FORMAT_RAW 4

#define ENCODING_SIZE 511 // Size of encoding histogram



int                infile;
FILE               *outfile;
char               inname[NAME_SIZE];
char               outname[NAME_SIZE];
int                output_format;



void main(int argc, char **argv)


{
char          ibuf[BUFF_SIZE];
char          obuf[BUFF_SIZE];
short         length, total_bytes, line, i;
long          long_length;
long          hist[ENCODING_SIZE] = { 0 };
int           out_bytes = RECORD_BYTES;
size_t        count;   


/*********************************************************************/
/*                                                                   */
/* get input and output files                                        */
/*                                                                   */
/*********************************************************************/


   strcpy_s(inname, NAME_SIZE, "   ");
   strcpy_s(outname, NAME_SIZE, "   ");
   output_format = FORMAT_UNDEFINED;

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
      strcpy_s(inname, NAME_SIZE,argv[1]);
      if (argc >= 3) strcpy_s(outname, NAME_SIZE, argv[2]);
      if (argc == 3) output_format = FORMAT_PDS;
      int dummy = -1;
      if (argc == 4) dummy = sscanf_s(argv[3], "%d", &output_format); 
     }

   get_files(); 

/*********************************************************************/
/*                                                                   */
/* read and edit compressed file labels                              */
/*                                                                   */
/*********************************************************************/

   switch (output_format)
     {
       case FORMAT_PDS : pds_labels();
               break;
       case FORMAT_FITS : fits_labels();
               break;
       case FORMAT_VICAR: vicar_labels();
               break;
       case FORMAT_RAW: no_labels();
     }
/*********************************************************************/
/*                                                                   */
/* process the image histogram                                       */
/*                                                                   */
/*********************************************************************/

/* need to know record_bytes,hist_count,hist_item_type,item_count.*/
   total_bytes = 0;
   length = read_var((char *)hist);
   total_bytes = total_bytes + length;
   length = read_var((char *)hist + RECORD_BYTES);
   total_bytes = total_bytes + length;

  
   if (output_format == FORMAT_PDS)
     {
      fwrite((char *)hist, RECORD_BYTES,1,outfile);
      fwrite((char *)hist + RECORD_BYTES,length,1,outfile);

      /*  pad out the histogram to a multiple of RECORD_BYTES */
      for (int i = total_bytes; i < RECORD_BYTES*2; i++) fputc(' ', outfile);
     }
/*********************************************************************/
/*                                                                   */
/* process the encoding histogram                                    */
/* don't have to byte-swap because DECOMP.C does it for us           */
/*                                                                   */
/*********************************************************************/

   length = read_var((char *)hist);
   length = read_var((char *)hist + RECORD_BYTES);
   length = read_var((char *)hist + 2* RECORD_BYTES);

/*********************************************************************/
/*                                                                   */
/* process the engineering summary                                   */
/*                                                                   */
/*********************************************************************/

   total_bytes = 0;
   length = read_var(ibuf);

   if (output_format == FORMAT_PDS )
     {
      fwrite(ibuf,length,1,outfile);
      total_bytes = total_bytes + length;

      /*  pad out engineering to multiple of RECORD_BYTES (836) */
      for (i=total_bytes; i<RECORD_BYTES; i++) fputc(' ', outfile);
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
       length = read_var(ibuf);
       if (length <= 0) break;
       long_length = (long)length;
       line += 1;
       dcmprs(ibuf, obuf,&long_length, &out_bytes, tree);
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
 if (output_format == FORMAT_FITS)
   for (int i = 0; i < 2240; i++) fputc(' ', outfile);

 printf("\n");
 free_tree(&long_length);
 _close(infile);
 fclose(outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine get_files - get input filenames and open               */
/*                                                                   */
/*********************************************************************/

void get_files()

{
//short   shortint;

  if (inname[0] == ' '){

     printf("\nEnter name of file to be decompressed: ");
     gets_s (inname, 79);
  }
  _set_errno(0);
    if (( _sopen_s(&infile,inname,_O_RDONLY | _O_BINARY, _SH_DENYWR, _S_IREAD)) != 0){

        printf("\ncan't open input file: %s\n",inname);
        exit(1);
    }
  
  if (output_format == FORMAT_UNDEFINED)
  do
  {
     printf("\nEnter a number for the output format desired:\n");
     printf("\n  1.  SFDU/PDS format.");
     printf("\n  2.  FITS format.");
     printf("\n  3.  VICAR format.");
     printf("\n  4.  Unlabelled binary array.\n");
     printf("\n  Enter format number:");
     gets_s(inname, 79);
     output_format = atoi(inname);
  } while (output_format < FORMAT_PDS || output_format > FORMAT_RAW);

  if (outname[0] == ' ')
  {
     printf("\nEnter name of uncompressed output file: ");
     gets_s (outname,79);
  }

  _set_errno(0);
  if (fopen_s(&outfile, outname,"wb") != 0){

       printf("\ncan't open output file: %s\n",outname);
       exit(1);
   }
  
   return;  
}



/*********************************************************************/
/*                                                                   */
/* subroutine pds_labels - edit PDS labels and write to output file  */
/*                                                                   */
/*********************************************************************/

void pds_labels()

{
char          outstring[NAME_SIZE];
char          ibuf[BUFF_SIZE];
short         length;
short         total_bytes;
short         i;


total_bytes = 0;
do
  {
   length = read_var(ibuf);
   ibuf[length] = '\0';// NULL;

  /******************************************************************/
  /*  edit labels which need to be changed                          */
  /******************************************************************/

   if      ((i = strncmp(ibuf,"NJPL1I00PDS1",12)) == 0)
   /*****************************************************************/
   /* add the output file length to the sfdu label                  */
   /*****************************************************************/
     {
      strcpy_s(outstring, NAME_SIZE,ibuf);
      strcpy_s(outstring+12, NAME_SIZE, "00673796");
      strcpy_s(outstring+20, NAME_SIZE, ibuf+20);
      fwrite(outstring,length,1,outfile);
      fprintf(outfile,"\r\n");
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"RECORD_TYPE",11)) == 0)
   /*****************************************************************/
   /* change the record_type value from variable to fixed           */
   /*****************************************************************/
     {
      strcpy_s(ibuf + 35, BUFF_SIZE, "FIXED_LENGTH");
      length = length - 3;
      fwrite(ibuf, length, 1, outfile);
      fprintf(outfile, "\r\n");
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"FILE_RECORDS",12)) == 0)
   /*****************************************************************/
   /* change the file_records count to 806                          */
   /*****************************************************************/
     {
      strcpy_s(ibuf+35, BUFF_SIZE, "806");
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"LABEL_RECORDS",13)) == 0)
   /*****************************************************************/
  /* change the label_records count from 56 to 3                    */
   /*****************************************************************/
     {
      strcpy_s(ibuf+35, BUFF_SIZE, "3");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"^IMAGE_HISTOGRAM",16)) == 0)
   /*****************************************************************/
   /* change the location pointer of image_histogram to record 4    */
   /*****************************************************************/
     {
      strcpy_s(ibuf+35, BUFF_SIZE, "4");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
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
      strcpy_s(ibuf+35, BUFF_SIZE, "6");
      length -= 1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
      total_bytes = total_bytes + length + 2;
     }
   else if ((i = strncmp(ibuf,"^IMAGE",6)) == 0)
   /*****************************************************************/
   /* change the location pointer of image to record 7              */
   /*****************************************************************/
     {
      strcpy_s(ibuf+35, BUFF_SIZE, "7");
      length = length -1;
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
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
         length = read_var(ibuf);
        }
     }
   else if ((i = strncmp(ibuf," ENCODING",9)) == 0);
  
   /*****************************************************************/
   /* if none of above write out the label to the output file       */
   /*****************************************************************/
   else
     {
      fwrite(ibuf,length,1,outfile);
      fprintf(outfile,"\r\n");
      total_bytes = total_bytes + length + 2;
     }
   /*****************************************************************/
   /* test for the end of the PDS labels                            */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (i=total_bytes;i<RECORD_BYTES*3;i++) fputc(' ', outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine fits_labels - write FITS header to output file */
/*                                                                   */
/*********************************************************************/

void fits_labels()

{
char        ibuf[BUFF_SIZE];
char        outstring[NAME_SIZE];
short       length,total_bytes,i;

do
  {
   length = read_var(ibuf);
   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

total_bytes = 0;

strcpy_s(outstring, NAME_SIZE, 
"SIMPLE  =                    T                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

strcpy_s(outstring, NAME_SIZE, 
"BITPIX  =                    8                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

strcpy_s(outstring, NAME_SIZE, 
"NAXIS   =                    2                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

strcpy_s(outstring, NAME_SIZE, 
"NAXIS1  =                  800                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

strcpy_s(outstring, NAME_SIZE, 
"NAXIS2  =                  800                                                ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

strcpy_s(outstring, NAME_SIZE, 
"END                                                                           ");
fwrite(outstring,78,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 80;

/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (int i = total_bytes; i < 2880; i++) fputc(' ', outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine vicar_labels - write vicar labels to output file       */
/*                                                                   */
/*********************************************************************/

void vicar_labels()
{
char          ibuf[BUFF_SIZE];
char          outstring[NAME_SIZE];
short         length,total_bytes,i;

do
  {
   length = read_var(ibuf);
   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp(ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

total_bytes = 0;

strcpy_s(outstring, NAME_SIZE, 
"LBLSIZE=800             FORMAT='BYTE'  TYPE='IMAGE'  BUFSIZ=800  DIM=2  ");
fwrite(outstring,72,1,outfile);
total_bytes = total_bytes + 72;

strcpy_s(outstring, NAME_SIZE, 
"EOL=0  RECSIZE=800  ORG='BSQ'  NL=800  NS=800  NB=1  N1=0  N2=0  N3=0  ");
total_bytes = total_bytes + 71;
fwrite(outstring,71,1,outfile);

strcpy_s(outstring, NAME_SIZE, 
"N4=0  NBB=0  NLB=0");
fwrite(outstring,18,1,outfile);
fprintf(outfile,"\r\n");
total_bytes = total_bytes + 20;


/* pad out the labels with blanks to multiple of RECORD_BYTES */
   for (int i = total_bytes; i < 800; i++) fputc(' ', outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine no_labels - read past the pds labels                   */
/*                                                                   */
/*********************************************************************/

void no_labels()

{
char          ibuf[BUFF_SIZE];
short         length,i;

do
  {
   length = read_var(ibuf);
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

read_var(char * ibuf)
{
int   length,result,nlen;

    length = 0;
    result = _read(infile,&length,2);
    nlen =   _read(infile,ibuf,length + (length%2));

    return length;
}

static void decmpinit(long int *hist)
/***************************************************************************
*_TITLE decmpinit - initializes the Huffman tree                           *
*_ARGS  TYPE       NAME      I/O        DESCRIPTION                        *
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

static NODE *huff_tree(long int *hist)
/****************************************************************************
*_TITLE huff_tree - constructs the Huffman tree; returns pointer to root    *
*_ARGS  TYPE          NAME        I/O   DESCRIPTION                         *
        long int     *hist;     /* I    First difference histogram          */

  {
  /*  Local variables used */
    long int freq_list[512];      /* Histogram frequency list */
    NODE **node_list;             /* DN pointer array list */

    long int *fp;        /* Frequency list pointer */
    NODE **np;           /* Node list pointer */

    long int num_freq;   /* Number non-zero frequencies in histogram */
    //long int sum;                 /* Sum of all frequencies */

    short int num_nodes; /* Counter for DN initialization */
    short int cnt;       /* Miscellaneous counter */

    short int znull = -1;         /* Null node value */

    NODE *temp;          /* Temporary node pointer */

  /* Functions called */
    void sort_freq();
    NODE *new_node();
    //char *malloc();

/***************************************************************************
  Allocate the array of nodes from memory and initialize these with numbers
  corresponding with the frequency list.  There are only ENCODING_SIZE (511) possible
  permutations of first difference histograms.  There are 512 allocated
  here to adhere to the FORTRAN version.
****************************************************************************/

   fp = freq_list;
   node_list = (NODE **) malloc(sizeof(temp)*512);  // list of node pointers
   if (node_list == NULL)
    {
      printf("\nOut of memory in huff_tree!\n");
      exit(1);
    }
   np = node_list;

   for (num_nodes = 1, cnt = 512 ; cnt-- ; num_nodes++)
     {
/**************************************************************************
    The following code has been added to standardize the VAX byte order
    for the "long int" type.  This code is intended to make the routine
    as machine independant as possible.
***************************************************************************/
        unsigned char *cp = (unsigned char *) hist++;
        unsigned long int j;
        short int i;
        for (i = 4 ; --i >= 0 ; j = (j << 8) | *(cp + i));

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
  sort_freq(freq_list, node_list, num_freq);

  fp = freq_list;
  np = node_list;

  for (num_freq = 512 ; (*fp) == 0 && (num_freq) ; fp++, np++, num_freq--);


/***************************************************************************
  Now create the tree.  Note that if there is only one difference value,
  it is returned as the root.  On each interation, a new node is created
  and the least frequently occurring difference is assigned to the right
  pointer and the next least frequency to the left pointer.  The node
  assigned to the left pointer now becomes the combination of the two
  nodes and it's frequency is the sum of the two combining nodes.
****************************************************************************/

  for (temp = (*np) ; (num_freq--) > 1 ; )
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

  //char *malloc();       /* Memory allocation function */

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

 static void sort_freq(long int *freq_list, NODE **node_list, long int num_freq)
/****************************************************************************
*_TITLE sort_freq - sorts frequency and node lists in increasing freq. order*
*_ARGS  TYPE       NAME            I/O  DESCRIPTION                         *
        long int   *freq_list;   /* I   Pointer to frequency list           *
        NODE       **node_list;  /* I   Pointer to array of node pointers   *
        long int   num_freq;     /* I   Number of values in freq list       */

  {
    /* Local Variables */
    long int *i;       /* primary pointer into freq_list */
    long int *j;       /* secondary pointer into freq_list */

    NODE **k;          /* primary pointer to node_list */
    NODE **l;          /* secondary pointer into node_list */

    long int temp1;             /* temporary storage for freq_list */
    NODE *temp2;                /* temporary storage for node_list */

    long int cnt;      /* count of list elements */


/************************************************************************
  Save the current element - starting with the second - in temporary
  storage.  Compare with all elements in first part of list moving
  each up one element until the element is larger.  Insert current
  element at this point in list.
*************************************************************************/

   if (num_freq <= 0) return;      /* If no elements or invalid, return */

   for (i = freq_list, k = node_list, cnt = num_freq ; --cnt ; *j=temp1, *l=temp2)
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

 static void dcmprs(char *ibuf, char *obuf, long int *nin, long int *nout, NODE *root)
/****************************************************************************
*_TITLE dcmprs - decompresses Huffman coded compressed image lines          *
*_ARGS  TYPE       NAME       I/O       DESCRIPTION                         *
        char       *ibuf;   /* I        Compressed data buffer              *
        char       *obuf;   /* O        Decompressed image line             *
        long int   *nin;    /* I        Number of bytes on input buffer     *
        long int   *nout;   /* I        Number of bytes in output buffer    *
        NODE       *root;   /* I        Huffman coded tree                  */

  {
    /* Local Variables */
    NODE *ptr = root;        /* pointer to position in tree */
    unsigned char test;      /* test byte for bit set */
    unsigned char idn;       /* input compressed byte */

    char odn;                /* last dn value decompressed */

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

    for (idn = (*ibuf) ; ibuf < ilim  ; idn =(*++ibuf))
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


static void free_tree(long int *nfreed)
/****************************************************************************
*_TITLE free_tree - free memory of all allocated nodes                      *
*_ARGS  TYPE       NAME       I/O        DESCRIPTION                        *
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

	*nfreed = free_node(tree, total_free);

	return;
}


static long int free_node(NODE *pnode, long int total_free)
/***************************************************************************
*_TITLE free_node - deallocates an allocated NODE pointer
*_ARGS  TYPE     NAME          I/O   DESCRIPTION                           *
        NODE     *pnode;       /* I  Pointer to node to free               *
        long int total_free;   /* I  Total number of freed nodes           *

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
	if (pnode ==  NULL) return(total_free);
	
	if (pnode->right !=  NULL)
		total_free = free_node(pnode->right, total_free);
	if (pnode->left !=  NULL)
		total_free = free_node(pnode->left, total_free);

	free(pnode);
	return(total_free + 1);
}
