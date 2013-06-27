To build html documentation run

      doxygen
	
from the command line in this directory, or otherwise run 
doxygen using the Doxygen config file here.  The root of
the documentation is then

      html/index.html
	
To turn on the latex (and thus pdf) output, in Doxyfile set

      GENERATE_LATEX = YES
	
and run pdflatex from the latex directory.
	