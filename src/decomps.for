      SUBROUTINE DECOMPRESS (IBUF, OBUF, NIN, NOUT)
C**********************************************************************
C_TITLE DECOMPRESS - Decompresses image lines stored in a compressed 
C                    format. (alternate fortran souce code which will
C                    run on a SUN workstation)
C
C_ARGS
      INTEGER*4       NIN             
C                                       Input - number of bytes in
C                                       the compressed buffer.
      INTEGER*4       NOUT            
C                                       Input - known number of bytes
C                                       in decompressed output line
      CHARACTER       IBUF(1)         
C                                       Input - line of data to 
C                                       decompress. The first byte is 
C                                       the actual pixel value, the 
C                                       rest of the array is the bit 
C                                       stream of coded "first
C                                       difference" values.
      CHARACTER       OBUF(1)         
C                                       Output - decompressed line of 
C                                       data. The routine assumes 
C                                       original image composed of 8 
C                                       bit pixels.  This is the 
C                                       restored image line - 
C                                       compression and first 
C                                       differences undone - of the 
C                                       line passed in ibuf.
C
C_DESC  This routine decompresses Huffman encoded compressed data.  
C       (Huffman compression encoding can be found in any standard 
C       data compression reference.  One such book is: "Data 
C       Compression, Techniques and Applications, Hardware and 
C       Software Considerations" by Gilbert Held, John Wiley & Sons, 
C       publishers.)  Huffman coding uses variable number of bits to 
C       encode different values of the original data.  The compression 
C       results because the most frequently occurring values have
C       the smaller number of bits for their codes.
C
C       The routine, used in conjuction with DECMPINIT,  provides 
C       a common data base from which to call the processing
C       subroutines. DECMPINIT builds the Huffman tree from the first-
C       difference histogram and is called only once per image.  
C       DECOMPRESS processes one compressed input line per call and 
C       returns it completely restored.
C
C_CALLS This routine calls two subroutines DCMPRS and HUFFTREE.
C
C
C_LIMS  The fortran code found in the DECOMPS.FOR file was adapted
C       to run on a SUN work station using fortran version 3.4
C
C_HIST  DEC87, DMcMacken, ISD, USGS, Flagstaff, Original version
C
C_END
C***********************************************************************
      INTEGER*2       TREE(5,1024)    
C                                       Huffman tree
      INTEGER*2       VALUE
      INTEGER*2       RIGHT
      INTEGER*2       LEFT
      INTEGER*2       PARENT
      INTEGER*2       BRANCH
      PARAMETER       (VALUE=1, RIGHT=2, LEFT=3, PARENT=4, BRANCH=5)
      COMMON /DECOM/ TREE
C***********************************************************************
C Decompress an input line
C***********************************************************************
      CALL DCMPRS(IBUF, OBUF, NIN, NOUT, TREE)
      RETURN
      END
      SUBROUTINE DECMPINIT(HIST)
C***********************************************************************
C_TITLE DECMPINIT - initialize the decompression TREE from the first
C                   difference histogram
C_ARGS
      INTEGER*4       HIST(512)       
C                                       Input, difference histogram 
C                                       table. The histogram of the 
C                                       whole image after first 
C                                       difference has been run on 
C                                       each line.  In a "first 
C                                       difference" line the first pixel
C                                       retains its actual value; all 
C                                       other pixels are obtained by 
C                                       subtracting the actual value of 
C                                       the pixel from the actual value 
C                                       of the preceding pixel and 
C                                       adding 256 to provide a 
C                                       positive number.  The i-th 
C                                       element of the array HIST 
C                                       should be the number of pixels 
C                                       in the image with value i.
C
C Some computer hardware systems may need to swap the byte order
C of the 32-bit words which make up the first-difference historgram.  
C The histogram  is configured on the CDROM in "least significant byte 
C first" order. This is the order for integer values used by VAX and 
C IBM PC computer systems. Users of other computer architectures
C (IBM Mainframes, MacIntosh, SUN, and Apollo) will need to swap
C byte pairs 1 and 4, and 2 and 3 before passing the histogram to
C the DECMPINIT subroutine.
C                                     
C
C_DESC This routine initializes the Huffman tree using the first
C      difference histogram passed by the calling program. This 
C      subroutine must be called before the DECOMPRESS subroutine
C      can be utilized.
C
C_HIST DEC87, DMcMacken, ISD, USGS, Flagstaff, Original Version
C
C_END
C***********************************************************************
      INTEGER*2       TREE(5,1024)    
      COMMON /DECOM/ TREE

      DO 10 I = 1,1024
      DO 20 J = 1,5
      TREE(J,I) = 0
   20 CONTINUE
   10 CONTINUE
      CALL HUFFTREE (HIST, TREE)
      RETURN
      END
      SUBROUTINE HUFFTREE(HIST, TREE)
C**********************************************************************
C_TITLE HUFFTREE creates a Huffman code tree from input histogram
C
C_ARGS
      INTEGER*4       HIST(512)               
C                                       In - image histogram
      INTEGER*2       TREE(5,1024)            
C                                       Out - Huffman code tree
C
C_DESC  This routine creates a Huffman code tree from the input first
C       difference histogram of an image.  It starts by making all valid
C       (non-zero) densities for the image leaf nodes.  These are sorted
C       in increasing order by histogram frequency.  The two nodes with
C       smallest frequencies are connected by branches to a new node 
C       which is given a frequency equal to the sum of the two 
C       combining nodes. The new node takes the place of the two nodes 
C       and the frequency list is resorted.  The process is repeated 
C       until the frequency set is reduced to one node.  The tree 
C       consists of a double dimensioned array whose first dimension 
C       gives five parameters for each node. The first parameter gives 
C       a value to the node.  This is a density value for leaf nodes 
C       and a -1 for all others.  The second parameter is a pointer to 
C       the next node on the right branch from this node. The third 
C       points to the node on the left branch.  The fourth points
C       to the parent node from which this node branches.  The fifth
C       parameter notes whether this node is on the right or left 
C       branch of the parent node.
C
C_CALLS This routine calls the specially provided sort routine
C               SORTFREQ.
C
C_HIST  Fall86, DMcMacken, ISD, USGS, Flagstaff, Original version
C
C_END
C**********************************************************************
      INTEGER*4       FREQLIST(512)          
C                                       Histogramm frequency list
      INTEGER*2       NODELIST(512)          
C                                       Tree node list
      INTEGER*2       VALUE
      INTEGER*2       RIGHT
      INTEGER*2       LEFT
      INTEGER*2       PARENT
      INTEGER*2       BRANCH
      PARAMETER       (VALUE=1,RIGHT=2,LEFT=3,PARENT=4,BRANCH=5)
      INTEGER*4       NUMNODES               
C                                       Node counter
      INTEGER*4       I                       
C                                       Do loop index
      INTEGER*4       PTR                     
C                                       Current node pointer
      INTEGER*4       NUMFREQ                
C                                       # non-zero frequencies
      INTEGER*4       INDEX                   
C                                       Frequency list pointer
C**********************************************************************
C Create leaf nodes, and pointer tables
C**********************************************************************
      NUMNODES = 0
      DO 10 I = 1, 512
         FREQLIST(I) = HIST(I)
         NUMNODES = NUMNODES + 1
         PTR = NUMNODES
         NODELIST(I) = PTR
         TREE(VALUE, PTR) = I
         TREE(PARENT, PTR) = 0
   10 CONTINUE
C**********************************************************************
C Make sure the last element is always 0. For the first difference
C histogram, there can only be 511 values.
C**********************************************************************
      FREQLIST(512) = 0
C**********************************************************************
C Sort freqlist in increasing order; skip over zero frequencies
C**********************************************************************
      NUMFREQ = 512
      CALL SORTFREQ (FREQLIST, NODELIST, NUMFREQ)
      INDEX = 1
   20 IF (FREQLIST(INDEX) .EQ. 0) THEN
         INDEX = INDEX + 1
         GOTO 20
      ENDIF

      NUMFREQ = NUMFREQ - INDEX + 1
   30 IF (NUMFREQ .GT. 1) THEN
         NUMNODES = NUMNODES + 1
         PTR = NUMNODES
         TREE(RIGHT, PTR) = NODELIST(INDEX)
         TREE(LEFT, PTR) = NODELIST(INDEX+1)
         TREE(PARENT, TREE(RIGHT, PTR)) = PTR
         TREE(BRANCH, TREE(RIGHT, PTR)) = RIGHT
         TREE(PARENT, TREE(LEFT, PTR)) = PTR
         TREE(BRANCH, TREE(LEFT, PTR)) = LEFT
         TREE(VALUE, PTR) = -1

         FREQLIST(INDEX + 1) = FREQLIST(INDEX)
     1   + FREQLIST(INDEX + 1)
         NODELIST(INDEX + 1) = PTR
         FREQLIST(INDEX) = 0
         INDEX = INDEX + 1
         NUMFREQ = NUMFREQ - 1
         CALL SORTFREQ (FREQLIST(INDEX), NODELIST(INDEX),
     1   NUMFREQ)
         GOTO 30
      ENDIF
      TREE(VALUE, 1024) = NUMNODES
      RETURN
      END
      SUBROUTINE SORTFREQ(FREQLIST, NODELIST, NUMFREQ)
C**********************************************************************
C_TITLE SORTFREQ sorts frequency and node lists in increasing freq. 
C       order.
C
C_ARGS
      INTEGER*4       FREQLIST(512)          
C                                       input - frequency list
      INTEGER*2       NODELIST(512)          
C                                       input - node list
      INTEGER*4       NUMFREQ                
C                                       input - # values in freq list
C
C_DESC  This routine uses an insertion sort to reorder a frequency list
C       in order of increasing frequency.  The corresponding elements
C       of the node list are reordered to maintain correspondence.
C
C_HIST  Fall86, DMcMacken, ISD, USGS, Flagstaff, Original version
C
C_END
C***********************************************************************
      INTEGER*4       I                       
C                                       Do loop index
      INTEGER*4       J                       
C                                       List position pointer
      INTEGER*4       TEMP1                   
C                                       Temporary storage - freq.
      INTEGER*4       TEMP2                   
C                                       Temporary storage - node
C***********************************************************************
C Save current element - starting with second - in temporary
C storage.  Compare with all elements in first part of list -
C moving each up one element - until  element is larger.
C Insert current element at this point in list
C***********************************************************************
      DO 20 I = 2, NUMFREQ
         TEMP1 = FREQLIST(I)
         TEMP2 = NODELIST(I)
         J = I
   10    IF (FREQLIST(J-1) .GT. TEMP1) THEN
            FREQLIST(J) = FREQLIST(J-1)
            NODELIST(J) = NODELIST(J-1)
            J = J - 1
            IF (J .GT. 1) GOTO 10
         ENDIF
         FREQLIST(J) = TEMP1
         NODELIST(J) = TEMP2
   20 CONTINUE
      RETURN
      END
      SUBROUTINE DCMPRS(IBUF, OBUF, NIN, NOUT, TREE)
C***********************************************************************
C_TITLE DCMPRS decompresses Huffman coded compressed image lines
C       (Remove this subroutine from the source code if you are
C       planning to use the VAX/VMS DCMPRS.MAR macro routine in
C       place of this fortran code. The DCMPRS macro version is more 
C       than twice as fast as the fortran version)
C
C_ARGS
      INTEGER*4       NIN             
C                                       Input, # bytes in input buffer
      INTEGER*4       NOUT            
C                                       Input, # bytes in output line
      INTEGER*2       TREE(5,1024)    
C                                       Input, Huffman code tree
      CHARACTER       IBUF(1)         
C                                       Input, compressed data buffer
      CHARACTER       OBUF(1)         
C                                       Output, decompressed line
C
C_DESC  This subroutine follows a path from the root of the Huffman tree
C       to one of its leaves.  The choice at each branch is decided by
C       the successive bits of the compressed input stream, left for 1,
C       right for 0.  Only leaf nodes have a value other than -1.  The
C       routine traces a path through the tree until it finds a node 
C       with a value not equal to -1 (a leaf node).  The value at the 
C       leaf node is subtracted from the preceding pixel value plus 256 
C       to restore the the umcompressed pixel.  This algorithm is then 
C       repeated until the entire line has been processed.
C
C_CALLS This routine uses the VAX/VMS implicit function
C       BTEST (VAR, IBIT) to test whether the bit number IBIT is set 
C       in the variable VAR. A fortran version of this routine
C       is available for BTEST.
C
C_HIST  DEC86, DMcMacken, ISD, USGS, Flagstaff, Original version
C
C_END
C**********************************************************************
      INTEGER*4       PTR             
C                                       Pointer to position in tree
      INTEGER*4       K               
C                                       Do loop index
      INTEGER*4       BIT             
C                                       Pointer to current bit in word
      INTEGER*4       NBYTE           
C                                       Counter, # bytes decompressed
      INTEGER*4       TEMP            
C                                       Working variable
      INTEGER*2       VALUE           
C                                       Value index in tree
      INTEGER*2       RIGHT           
C                                       Right branch index in tree
      INTEGER*2       LEFT            
C                                       Left branch index in tree
      INTEGER*2       PARENT          
C                                       Parent pointer index in tree
      INTEGER*2       BRANCH          
C                                       Parent branch index in tree
      PARAMETER       (VALUE=1, RIGHT=2, LEFT=3, PARENT=4, BRANCH=5)
C
      INTEGER*2       DN              
C                                       Current density value
      INTEGER*2       DNL             
C                                       Last density value
      LOGICAL       END             
C                                       End of line flag
      LOGICAL       BTEST
C
C**********************************************************************
C Start at root of tree
C**********************************************************************
      PTR = TREE(1,1024)
C**********************************************************************
C Just copy first byte it is uncompressed.
C**********************************************************************
      OBUF(1) = IBUF(1)
C**********************************************************************
C Count starts here
C**********************************************************************
      NBYTE = 1
C**********************************************************************
C Initialize current and last density values
C**********************************************************************
      DN = ICHAR(IBUF(1))
      IF (DN.LT.0) DN = 256 + DN
      DNL = DN
C**********************************************************************
C Reset end of line flag
C**********************************************************************
      END = .FALSE.
C**********************************************************************
C K points to current input byte loop is done when k exceeds # 
C input bytes.
C**********************************************************************
      K = 2
      IF (K .GT. NIN) END = .TRUE.
C**********************************************************************
C This is the processing loop
C**********************************************************************
   10 IF (.NOT. END) THEN
         TEMP = ICHAR(IBUF(K))
         BIT = 7                         
   20    IF (BIT .GE. 0 .AND. .NOT. END) THEN 
            IF (BTEST(TEMP, BIT)) THEN
               PTR = TREE(LEFT, PTR)
            ELSE
               PTR = TREE(RIGHT, PTR)
            ENDIF
C**********************************************************************
C We are at end leaf when value is not -1
C**********************************************************************
            IF (TREE(VALUE, PTR) .NE. -1) THEN
C**********************************************************************
C Compute value of pixel, reset or increment pointers and counters
C**********************************************************************
               DN = TREE(VALUE, PTR)
               DN = DNL - DN + 256
               DNL = DN
               NBYTE = NBYTE + 1
               OBUF(NBYTE) = CHAR(DN)
               PTR = TREE(1, 1024)
               IF (NBYTE .EQ. NOUT) END = .TRUE.
            ENDIF
C**********************************************************************
C Process next bit in byte
C**********************************************************************
            BIT = BIT - 1
            GOTO 20
         ENDIF
C**********************************************************************
C Set index to next input byte
C**********************************************************************
         K = K + 1
         IF (K .GT. NIN) END = .TRUE.    
         GOTO 10
      ENDIF
C**********************************************************************
C Go back to caller with decompressed line
C**********************************************************************
      RETURN
      END
