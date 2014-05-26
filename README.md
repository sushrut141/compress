compress
========
JPEG image compression algorithm has been partly implemented using opencv liraries.
Huffman encoding has not been implemented.
The output of run level encoding of each macroblock has been commented out.
The program takes an input image (preferable large), compresses it by eliminating high spatial frequencies and 
shows the output image.
There is also a provision to reorder macroblocks after calculating Discrete Cosine Transform and carry out Run level encoding to generate tuples for each block.
These tuples can be written to a file using the print _out function.
The functions for reordering, run level encoding, printing tuples have been commented out.
