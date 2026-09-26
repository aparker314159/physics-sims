# Numerical Methods Self-Study
I've been going through Paolo Giannozzi's book/lecture notes for numerical methods in quantum mechanics.
Unfortunately, the links to the code associated with the book seem to be broken, so I have had to write my own implementations 
from the book's descriptions. I figured I'd make my code public, in case it's useful for anyone else.

I've wrote my code in Rust, and it uses the plotters library with a cairo backend to display plots. 

I might add some other numerical methods physics simulators later, hence the very general repository name.

## Chapter 1
The solver seems to work for several symmetric potentials, including the harmonic one. The blue and red are the real and imaginary parts respectively, while the green is the probability amplitude.
The code here roughly corresponds to harmonic1.c (I'd assume - the book's code is unavailable). I added some of the routines described in the book's Appendix B, to allow for the wavefunction to evolve 
over time. The code is still very messy, since it's mostly for self-study rather than a real coding project.

