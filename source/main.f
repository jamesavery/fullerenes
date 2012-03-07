      PROGRAM Fullerene
      IMPLICIT REAL*8 (A-H,O-Z)
C---------------------------------------------------------------------- 
C
C  PROGRAM FULLERENE FOR TOPOLOGICAL ANALYSIS OF FULLERENES
C    A program to create cartesian coordinates for fullerenes
C    isomers and to perform a topological analysis.
C    The results can be used for plotting fullerene graphs
C    (Schlegel diagrams)i and structures, and as a starting point
C    for further quantum theoretical treatment. 
C    Version 4 now incorporates C++ routines linked to the
C    original Fortran program using improved algorithms.
C
C---------------------------------------------------------------------- 
C  R U N N I N G    T H E    P R O G R A M:
C---------------------------------------------------------------------- 
C* This program works with gfortran and gc compilers
C    You need to use the make file included in fullerene.zip 
C    The executable "fullerene" runs on a mac intel as:
C      ./fullerene <inp >out
C    If you type      make tests      it runs all the input jobs and
C      puts them into *.out
C    If you type      make clean      all the object files are deleted
C    It currently compiles in a 64 bit version, but you can change to
C      32 bits in the Makefile if necessary
C    All fortran files are in the directory    source
C    All C-files are in the directory          libgraph
C
C----------------------------------------------------------------------
C  G E N E R A L   D E S C R I P T I O N:
C---------------------------------------------------------------------- 
C    The program is written in standard Fortran and C++ (~10,000 lines)
C    Some standard routines from Mathematical Recepies were modified and
C    are used here for matrix diagonalization and geometry optimization.
C
C* Function: To perform a topological analysis of a regular fullerene 
C    (i.e. consisting of pentagons and hexagons) fulfilling Euler's theorem. 
C    The program calculates the volume and surface area of a fullerene 
C    (irregular or not). It further constructs the structure obtaining 
C    cartesian coordinates from the canonical ring spiral pentagon indices
C    through either Tutte embedding or Hueckel eigenvector method. Note that 
C    there is no unique definition for the volume of a fullerene for
C    nonplanar 5- or 6-rings on the fullerene surface, but there is no 
C    reason why any other definition than the tesselation algorithm or
C    convex hull chosen should be preferred. The Wu force-field and 
C    geometry optimization using a Fletcher-Reeves-Polak-Ribiere
C    minimization with analytical gradients is also implemented, providing 
C    good a initial guess for cartesian coordinates.
C
C    Note: This program works for any (distorted or not) regular fullerene
C     (i.e. a fullerene of genus 0 consisting of pentagons and hexagons only).
C     The spiral algorithm of Fowler and Manolopoulus is not
C     restricted to starting from a pentagon.
C     For a general list of fullerenes see "The House of graphs" at 
C     https://hog.grinvin.org/Fullerenes. 
C
C    Lit.: 1) P. Schwerdtfeger, J. Avery, "Topological Analysis of Fullerenes - 
C             A Fortran and C++ Program (Version 4.0)", Massey University Albany, 
C             Auckland, New Zealand (2012).
C          2) P. W. Fowler and D. E. Manopoulus, "An Atlas of Fullerenes" 
C             (Dover Publ., New York, 2006).
C          3) D. Babic, "Nomenclature and Coding of Fullerenes",
C             J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).
C          4) Z. C. Wu, D. A. Jelski, T. F. George, "Vibrational Motions of
C             Buckminsterfullerene", Chem. Phys. Lett. 137, 291-295 (1987).
C
C      Further reading:
C          5) D. E. Manopoulus and P. W. Fowler, "Molecular graphs, point groups, 
C             and fullerenes", J. Chem. Phys. 96, 7603-7614 (1992).
C          6) G. B. Adams, M. O'Keefe, and R. S. Ruoff, "Van der Waals Surface Areas
C             and Volumes of Fullerenes", J. Phys. Chem. 98, 9465-9469 (1994).
C          7) W. O. J. Boo, "An Introduction to Fullerene Structures",
C             J. Chem. Ed. 69, 605-609 (1992).
C          8) D. Babic, D. J. Klein and C. H. Sah, "Symmetry of fullerenes",
C             Chem. Phys. Lett. 211 (1993) 235-241.
C          9) R. Tonner, G. Frenking, M. Lein, and P. Schwerdtfeger, 
C             "Packed to the Rafters – Filling up C60 with Rare Gas Atoms", 
C             Chem. Phys. Chem. 12, 2081-2084 (2011).
C         10) G. Brinkmann, J. Goedgebeur, B. D. McKay, "The smallest fullerene 
C             without a spiral", Chem. Phys. Lett. 522, 54–55 (2012).
C
C      There is a paper in preparation which explains most of the features
C        in this program. You may ask for a preliminary copy.
C
C      If you use the CYLview program, the input file is written into cylview.dat.
C        CYLview plots molecules, is written by C. Y. Legault and is available
C        under http://www.cylview.org/Home.html
C
C      If you use QMGA for plotting fullerene graphs the input file is written
C        in qmga.dat. The program is available at http://qmga.sourceforge.net/
C        In this case you need to cite A. T. Gabriel, T. Meyer, G. Germano, 
C        "Molecular graphics of convex body fluids", J. Chem. Theory Comput. 
C        4, 468-476 (2008). There will be soon another algorith availble to
C        plot fullerene graphs.
C
C   -----------------------------------------------------------------------------
C   |  NOTE: You may do whatever you like with this program, but if you use it  |
C   |  and publish data please cite at least references 1-4 given above.        |
C   |  The book by Fowler and Manopoulus is highly recommended.                 |
C   |  It helps understanding how this program works (the book is fun to read!).|
C   |  Many of the concepts used in this program can be found in this book.     |
C   -----------------------------------------------------------------------------
C
C   -> Many definitions depend on the use of Angstroems, so please use this unit.
C
C   Important steps in the program are:
C--------------------------------------
C - Create the structure of fullerene:
C    Read in Cartesian coordinates for a fullerene (see files c20.inp to c540.inp), 
C    or construct them for the Ih isomer (files ico.inp, icoExp.inp, icoideal.inp)
C    of C60, or get cartesian coordinates from ring spiral pentagon indices
C    by the Fowler-Manopoulus algorithm (as defined in refs.2 and 5) or Tutte
C    embedding algorithm (used by files pentagon1.inp to pentagon25.inp). 
C    Barycenter of the fullerene is set to the origin. See ref.2 for the 
C    construction of cartesian coordinates from ring spiral pentagon indices and 
C    the use of Hueckel P-type eigenvectors. Note it is critical to get the right
C    vectors for the construction of the cartesian coordinates. The 3 P-type
C    vectors may need to be read in (see file pentagon8.inp for such an example).
C    It is important that the end-product is viewed by molecular visualization.
C    program. We recommend CYLview by Claude Legault (see http://www.cylview.org).
C    For this purpose a file is written out to cylview.xyz to be used as an input 
C    file. If this algorithm fails you can use our Tutte embedding algorithm 
C    (which to our opinion is less trouble and may be used as the standard method).
C    Using these coordinates the program calculates the smallest and largest cage 
C    diameters which gives already a measure for distortion from spherical symmetry
C    (Subroutine DIAMETER). It produces the distance matrix if print level is set
C    to high (Subroutine DISTMATRIX). The diameters already indicate if the
C    fullerene is heavily distorted from spherical symmetry.
C
C - Use program SPIRAL of Fowler and Manopoulus for creating a ring spiral.
C    It also produces canonical ring spiral pentagon indices 
C    (Subroutine SPIRALSEARCH) if cartesian coordinate input is chosen.
C
C-  Print all isomers and perform analysis introduced mostly in the book 
C    written by Fowler and Manopoulus (ref.2), e.g. pentagon indices and 
C    pentagon number, hexagon indices and strain parameter, NMR information
C    and number of distinct Hamiltonian cycles if required.
C
C - Perform Hueckel analysis (Subroutine HUECKEL). 
C    This gives you a good hint if the fullerene is open or closed shell.
C    Note: It is known that for fullerenes the Hueckel analysis is not
C      very reliable as hybridization with the C(2s) orbitals occur due
C      to nonplanarity. Hence the sigma-pi separation breaks down.
C      Nevertheless, we adopt  alpha=-0.21 au  and  beta=-0.111 au obtained
C      from the exp. ionization potential (7.58 eV) and excitation energy
C      (3.02 eV) of C60 (note the electron goes into the second LUMO of 
C      t1g symmetry). This gives also orbital energies for C60 in reasonable 
C      agreement with DFT Kohn-Sham orbital energies. The fullerene might
C      also Jahn-Teller distort adopting a singlet ground state instead
C      of one of a higher multiplicity. Such electronic effects are not
C      captured in this program, and one needs to perform a proper quantum
C      theoretical calculation.
C
C - Use program HAMILTON of Babic (ref.3) for Hamiltonian cycles and IUPAC 
C    Nomenclature. The number of Hamiltonian cycles given has been checked 
C    against a second algorithm, so it should work. Note, the number gives 
C    all distinct Hamiltonian cycles and left-right cycles counted as the 
C    same. Although finding all Hamiltonian cycles by the back-track algorithm
C    of Babich used here is a NP-complete problem, it works fine up to 
C    about C100. After that it becomes computationally very demanding. 
C    Note also that the existence of Hamiltonian cycles for fullerenes is 
C    only conjectured, and only for layered fullerenes (e.g. fullerene 
C    nanotubes) it has been proven to exist, our calculations show that
C    they exist for all fullerene isomers up to C100.
C
C - Establish connectivities between atoms either from given cartesian
C    coordinates or from adjacency matrix if ring spiral as input is used:
C    1) Bond between atoms, 2) Vertices (Subroutine CONNECT).
C    Identify all closed 5- and 6-ring systems (Subroutine RING).
C    This routine also determines if Euler's theorem is fulfilled.
C
C - Determine the center for each 5- and 6-ring (Subroutine RINGC)
C    This is required for the trigonal pyramidal tessellation to obtain
C    the volume and surface.  This routine also analyzes all 2- and 3-ring fusions
C    It further gives the Rhagavachari-Fowler-Manoupoulos neighboring pentagon 
C    and hexagon indices as described in the Fowler and Manolopoulos book
C    (ref.2). From the hexagon indices one derives if the fullerne fulfills the
C    IPR or not.
C
C - Fletcher-Reeves-Polak-Ribiere geometry optimization using analytical 
C    gradients for the Wu force field (Subroutine OPTFF).
C    It is very fast, even for C840. Note that the force field optimization
C    might distort the fullerene from the ideal point group symmetry.
C    On the other hand, the construction of the fullerene by using
C    pentagon indices leads to a more spherical arrangement in both
C    algorithms (Fowler-Manoupoulos or Tutte), e.g. barrels instead of
C    nanotubes.
C
C - Calculate the volume of the fullerene by summing over all
C    tetrahedrons spanned by the three vectors (Subroutine VOLUME).
C     CM-CR  (center of cage to the center of ring)
C     CM-CA1 (center of cage to atom 1 in ring)
C     CM-CA2 (center of cage to atom 2 in ring)
C    There are 5 such tetrahedrons in a 5-ring and 6 in a 6-ring
C    Note that CM is already in the origin
C    Let CR=(X1,Y1,Z1) , CA1=(X2,Y2,Z2) , and CA2=(X2,Y2,Z2)
C    Then the volume V for a irregular tetrahedron is given by 
C    the determinant
C
C                                 | X1 Y1 Z1 |
C     V = abs(Vdet)  ,   V =  1/6 | X2 Y2 Z2 |
C                                 | X3 Y3 Z3 |
C
C   Calculate the surface area A and the area/volume ratio (Subroutine VOLUME)
C
C                                        2              2              2
C                             | Y1 Z1 1 |    | Z1 X1 1 |    | X1 Y1 1 |   
C     A = 1/2 d**0.5 ,   d =  | Y2 Z2 1 | +  | Z2 X2 1 | +  | X2 Y2 1 |
C                             | Y3 Z3 1 |    | Z3 X3 1 |    | X3 Y3 1 |
C
C   Note that the ideal C60 coordinates can be constructed from scratch 
C    (Subroutine COORDC60). This routine was constructed to test the program 
C    for the case of an ideal capped icosahedron, where the analytical formula 
C    is well known, i.e.   V=[(125+43*sqrt(5))*R**3]/4
C    and R is the distance between the atoms (all the same)
C    Setting R=Rmin (Rmin is the smallest distance in C60) this gives a
C    lower bound for the volume, while the volume of the covering central sphere 
C    gives the upper bound.
C    For two different bond distances, R5 for the 5-ring and R6 for the 6-ring
C    joining another 6-ring, the volume can be determined as well after some tedious
C    algebraic manipulations:
C    i.e.   V=5[(3+sqrt(5))*(2R5+R6)**3]/12-[(5+sqrt(5))*R5**3]/2
C    For C20 (ideal dodecahedron) we have V=[(15+7sqrt(5))*R5**3]/4
C    For C20 and C60 the result of these formulae are also printed.
C   This method gives sensible results for convex fullerenes.
C
C - Calculate the minimum covering sphere (MCS) of the cage molecule 
C    (Subroutine MINCOVSPHERE): The MCS in m-dimensional space exists, is unique 
C     and can be expressed as a convex combination of at most (m+1) points, hence
C     our algorithm stops when 4 points are left over in the iteration process.
C
C    The problem can be reduced to
C
C    min(c) max(i) || p(I) - Cmcs ||
C
C    where ||..|| is the Eucledian norm in m dimensions and Cmcs is the center 
C    of the MCS.
C
C    Note: The spherical central cover SCC is not the minimum covering sphere MCS
C    (except if all distances from the center of points CM are the same as in the
C    ideal capped icosahedron). The spherical central cover is taken from the
C    CM point with radius Rmax (longest distance to one vertex).
C    The minimum covering sphere is calculated at the end of the
C    program using the algorithm of E. A. Yildirim, SIAM Journal on Optimization
C    Vol. 19(3),1368-1391 (2008) and the test by T. H. Hopp and C. P. Reeve,
C    NIST, US Department of Commerce (1996)'). Note that the much simpler algorithm 
C    by F. Lu and W. He, Global Congress on Intelligent Systems (2009),
C    DOI 10:1109/GCIS:2009:381, simply does not work for more than m+1 points on a 
C    surface or close by as their linear equation becomes linearly dependent.
C    Note also that the function Psi in Yildirim's algorithm is really the function
C    Phi defined earlier in his paper in section 2 and corrected in the SIAM paper.
C    His easier to program algorithm 1 was also tested, but is much slower.
C    If the value given in the iteration as "convergence" is close to
C    zero (equal zero), the iteration stops (if it falls below epsilon).
C    You can change the epsilon parameter in subroutine Sphere.
C    You can also try algorithm 1 of Yildirim through subroutine Sphere1
C    which is included in an file called algorithm1.f (although this file has not
C    been updated and further developed). Note that we changed the first
C    condition in this algorithm by choosing the furthest point from CM.
C    In the final statistics there should be 0 points outside the sphere
C    and at least 1 point on the sphere.
C    At the end the Van der Waals radius of carbon (1.415 Angstroems) is added to the
C    radius of the minimum covering sphere (note input coordinates for this need to be
C    in Angstroems otherwise change the program), and the volume of
C    an ideal fcc solid is calculated. The Van der Waals radius is chosen such that
C    for C60 the solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
C    are reproduced. The definition of the distortion parameter D from the MCS or for the
C    the isoperimetric quotient IPQ is
C
C    IPQ=36Pi(V^2/A^3)
C
C    D=[100/(N*Rmin)]* sum(i=1,N) {Rmcs - ||pi-Cmcs|| }    (N=MAtom)
C
C - Calculate the minimum distance sphere (MDS) of the fullerene.
C    The MCS definition for the distortion is biased for the case that few atoms stick 
C    out on a sphere and the MDS measure may be more appropriate for a measure
C    from spherical distortion. The MDS is defined as
C    The problem can be reduced to
C
C    min(Cmds) 1/N sum(i=1,N) | Rmds - || p(I) - Cmds || |
C
C    where ||..|| is the Eucledian norm in m dimensions. Cmds has to lie within the
C    convex hull. The MDS may not be uniquely defined, as there can be many 
C    (even degenerate) local minima, but for most spherical fullerenes it should 
C    just be fine. Analogous to the MCS there will be a measure for distortion 
C    from spherical symmetry.
C
C    D=[100/(N*Rmin)]* sum(i=1,N) | Rmds - || p(I) - Cmds || |
C
C - Calculate the maximum inner sphere (MCS) of the cage molecule
C
C    max(Cmds) min(i) || p(I) - Cmds ||
C
C    The maximum inner sphere is important for evaluating how much space
C    there is in a fullerene for encapsulating atoms and molecules. For
C    this the radius and volume is printed out with the Van der Waals
C    radius of carbon taken off Rmds. 
C
C - Produce the (X,Y) coordinates of a fullerene graph (Subroutine SCHLEGEL).
C    Here the points are rotated (if in input I1,I2, and I3 are given) so to
C    put the selected vertex, edge or ring center on top of the z-axis as the
C    point of projection (otherwise the point (0,0,zmax) is chosen with zmax
C    being the point with maximum z-value from the original input coordinates).
C    Points are then sorted in descending order according to their z-values.
C    The circumfences for atoms and rings down the z-axis are determined.
C    The Schlegel projection is created giving as output the projected (X,Y)
C    coordinates. The connections between the points are already written out
C    earlier in the output such that the fullerene graph can be drawn.
C    There are two choices for the projection, depending if you choose the
C    outer ring or the center part of the fullerene graph as a starting point:
C    1) The cone projection, i.e. points are projected out to an enveloping cone
C       and then down to a plane below the fullerene. The input I1,I2,I3 
C       defines the center of the fullerene graph. The last ring center should 
C       be at the bottom of the fullerene and if detected, will not be projected 
C       out, or if not will have a large scale factor (this center may be ignored 
C       in the drawing). Also, the last points on the outer ring in the fullerene
C       graph are scaled in distance by 1.2 in order to make the outer rings 
C       more visible. This also moves the outer centers within the ring.
C    2) The perspective projection, i.e. points are projected down a plane from
C       a set projection point. In this case the input I1,I2,I3 defines the outer
C       ring of the Schlegel diagram. 
C    From the output you can easily construct the name of the fullerene. 
C    At the end a rough printout of the fullerene graph is
C    produced. Note that this is o.k for fullerenes up to about C100, beyond it
C    it becomes too crowded and a proper plotting program should be used.
C    Nevertheless, it serves for a first rough picture. 
C    Furthermore, for large fullerenes it becomes critical to correctly set the
C    projection point or point of the cone. If for example the projection
C    point is too far away from the fullerene, edges may cross. 
C    
C-----------------------------------------------------------------------------------
C  G E N E R A L    I N F O R M A T I O N
C-----------------------------------------------------------------------------------
C Input and output files are in the folders  input  and   output  respectively.
C This program has been tested for the ideal capped icosahedron (input file ico.inp)
C   and for many other fullerenes which are found in the following input files:
C       C20 (c20.inp), C24 (c24.inp), C26 (c26.inp), C28 (c28.inp), C30 (c30.inp),
C       C36 (c36.inp), C50 (c50.inp), C60 (c60.inp), C70 (c70.inp), C72 (c72.inp),
C       C74 (c74.inp), C78 (c78.inp), C80 (c80.inp), C92 (c92.inp), C100 (c100.inp),
C       C180 (c180.inp), C320 (c320.inp), and C540 (c540.inp), pentagon1.inp, ...,
C       pentagon25.inp
C   The coordinates are mostly B3LYP aug-cc-pVDZ optimized up to C60, and
C    cc-pVDZ up to C180, and 6-31G for the rest and all for singlet states
C    (except of course for the ones where the pentagon indices input is
C    chosen). Note that for some of the fullerene coordinates the singlet
C    state chosen may not be the electronic ground state.
C Number of atoms is set in NATOMS currently at 900, so change this parameter
C   if you do not have enough RAM, alternatively the distances need to be
C   calculated directly and the DistMat(natomL) Matrix removed 
C   (which I can do if somebody insists). Also, maxit=2000000, which sets the
C   number of isomers or the number of Hamiltonian cycles to this value in order
C   for the program to not run forever.
C NB: The algorithm for locating all 5-and 6-rings might not be the smartest
C      one, but as this program requires only a second or so to run it was
C      not important to find a better algorithm.
C You should also be aware of program fullgen for generating nonisomorphic fullerenes.
C      It is written by Gunnar Brinkmann (Gunnar.Brinkmann@Ugent.be) and can be
C      downloaded from Brendan McKay's website (Australian National University) at
C      http://cs.anu.edu.au/~bdm/plantri/
C
C  This program is under permanent construction. The program has been tested for bugs.
C      Nevertheless, if you have any problem or find a bug please report to:
C                   peter.schwerdtfeger@gmail.com
C  Note: This program has neither been optimized for performance nor style.
C      It is however quite fast even for fullerenes such as C840.
C
C  Not implemented yet and on the to do list (in progress) is:
C      1) A subroutine to fit the minimum outer ellipsoidal cover 
C          useful for rugby ball like fullerenes.
C      2) Volume of the convex hull to compare to the tesselation method.
C      3) Produce symmetric fullerene graphs by optimizing the coordinates.
C         to avoid Schlegel projection which can result in edge crossings.
C      4) Geometry optimization using the extended Wu-Fowler force field.
C      5) Frequency calculations from the force field optimized geometry.
C      6) Construction of leap-frog fullerenes from adjacency matrix.
C      7) Construction of non-ring-spiral isomers using the genus algorithm.
C      8) Symmetry labels for Hueckel orbital energies.
C      9) Extend to non-regular fullerenes of genus 0 (heptagons and squares).
C     10) Extend to non-regular fullerenes of genus 1.
C     11) Convex Hull.
C     12) Symmetrize coordinates to point group symmetry.
C
C For any questions concerning this program please contact P. Schwerdtfeger
C   Centre for Theoretical Chemistry and Physics (CTCP)
C   The New Zealand Institute for Advanced Study (NZIAS) 
C   Massey University Auckland, Bldg.44, Private Bag 102904
C   North Shore City, 0745 Auckland, New Zealand
C   email: peter.schwerdtfeger@gmail.com
C   http://ctcp.massey.ac.nz/   and  http://www.nzias.ac.nz/
C   --> Always open to improvements and suggestions
C
C-----------------------------------------------------------------------------------
C  A C K N O W L E D G M E N T
C-----------------------------------------------------------------------------------
C I am indebted to the Alexander von Humboldt Foundation (Bonn) for financial support 
C in terms of a Humboldt Research Award, and to both Prof. Gernot Frenking and 
C Dr. Ralf Tonner (Marburg) for support during my extended stay in Marburg where I
C started to write this program. I acknowledge also the help of Darko Babich, Patrick 
C W. Fowler and David E. Manopoulus to allow the free distribution of their Fortran 
C subroutines.
C
C-----------------------------------------------------------------------------------
C  I N P U T:
C-----------------------------------------------------------------------------------
C Input (Either Namelist or in Free Format):   (use Angstroem for distances)
C
C 1) Text card of 80 characters (A80 format)
C
C 2) Input to create cartesian coordinates and main flags for the program
C    &Coord options /       (e.g. &Coord IC=20, IOPT=1, R6=1.42 /)
C    list of options: NA,IC,IP,IV1,IV2,IV3,icyl,TolR,R5,R6
C    NA= Number of Atoms (Default: 60)
C    IC= Flag for construction of cartesian coordinates (Default: 0)
C    IP= Print option (Default: 0)
C    IV1= Number for Hueckel P-type eigenvector for Fowler-Manopoulus algorithm (Default: 2)
C    IV2= Number for Hueckel P-type eigenvector for Fowler-Manopoulus algorithm (Default: 3)
C    IV3= Number for Hueckel P-type eigenvector for Fowler-Manopoulus algorithm (Default: 4)
C    icyl= Flag for producing input file for CYLview (Default: 0)
C    TolR= Tolerance in % (Default: 33)
C    R5= pentagon bond distance (Default: 1.455)
C    R6= hexagon  bond distance (Default: 1.391)
C    In detail:
C      If IC = 0 No coordinate input required, cordinates are constructed 
C              for the IPR isomer of C60
C           In this case only one card is read in:
C              R5,R6        (arbitrary units, e.g. Angstroms)
C              R5: Bond lengths in the pentagons 
C              R6: Bond length of the bonds connecting hexagons
C              If R5=R6 chosen then the ideal capped icosahedron is obtained
C      If IC = 1 Cartesian Coordinates expected as input
C         In this case N lines with    Z, X, Y, Z  in free format are expected.
C         (Z= Nuclear Charge, X,Y,C Cartesian Coordinates for Atom).
C         NB: Z is not really needed, but you can copy Gaussian output
C          directly into the input file
C      If IC = 2 or 3 Cartesian Coordinates are created from pentagon 
C             ring spiral list. Extra input required (free format):
C
C           IP(I),I=1,12 
C
C           R6 is taken as the smallest bond distance in the fullerene and IP(I)
C           identify the locations of the pentagons as described in detail
C           in ref.2, i.e. IP is the pentagon ring spiral numbering scheme. 
C           Note this only works for ring spiral fullerenes as described in 
C           detail by P. W. Fowler and D. E. Manopoulus. Use the canonical
C           pentagon ring indices if possible (transformation to the canonical
C           from should work as well).
C           IC=3: Fowler-Manopoulus algorithm using P-type eigenvectors:
C            If problem with eigenvectors are found to construct the
C            cartesian coordinates, i.e. the identification of P-type
C            eigenvectors, three integer values IV1, IV2, IV3 can be specified
C            identifying the eigenvectors to be chosen. pentagon8.inp is such an example.
C            In this case a severe warning occurs which means you should carefully
C            check the eigenvectors used and cartesian coordinates produced.
C            Otherwise coordinates are obtained which are useless. This is more
C            often the case as you might expect. 
C           IC=4: Tutte embedding algorithm. This should alway work, although
C            the initial fullerene might be too spherical.
C           Examples are given in the input files starting with 'pentagon'.
C           Please use Angstroems.
C      If iprintf>0 larger output produced, i.e. the full distance matrix, all
C          Hamiltonian cycles and all 3-ring connections.
C      Connectivities are found for atoms with distances between
C         R6   and   R6*(1+TolR/100)   if cartesian coordinate input is chosen.
C      If TolR=0. default value of 33% is used. 
C         NB: If this parameter is set at a value too large, unwanted connectivities
C          are produced resulting in smaller polygons. This parameter
C          should reflect the maximum deviation in distance from the
C          smallest distance found.
C
C 3) Option for force-field optimization:
C      &Opt options /        (e.g. &Opt Iopt=1 /)
C      list of options: Iopt
C      Iopt= Flag for force-field optimization (Default: 0)
C      In detail:
C       If Iopt=1  then fullerene is optimized using the force field method
C         of Wu et al within a Fletcher-Reeves-Polak-Ribiere algorithm:
C         Z.C.Wu, D.A.Jelski, T.F.George, Chem. Phys. Lett. 137, 291-295 (1987).
C         Note that there are more sophisticated force fields available,
C         but for these general parameters for fullerenes are not
C         yet available, and the Wu force field does the job to create good initial
C         cartesian coordinates for further refinement using more sophisticated
C         QM methods. Note that a converged energy much greater than zero implies
C         that the set distances and angles in the Wu force field cannot be
C         reached for all atoms and rings. For further reading see: 
C         A.Ceulemans, B.C.Titeca, L.F.Chibotaru, I.Vos, P.W.Fowler, 
C         J. Phys. Chem. A 105, 8284-8295 (2001).
C
C 4) Option for calculating Hamiltonian cycles and IUPAC numbers 
C      &Hamilton options /        (e.g. &Hamilton IHam=1 IUPAC=1 /)
C      list of options: IHam,IUPAC
C      In detail:
C      If IHam>0 Then Subroutine HAMILTON is called.     (Default: 0)
C         IHam=1 Routine will stop after 1 million Hamiltonian cycles are found
C         IHam>1 it will run up to 10**IHam Hamiltonian cycles.
C      If IUPAC=1 IUPAC numbers are produced.     (Default: 0)
C
C 5) Option for producing list of isomers and properties.
C      &Isomers options /        (e.g.&Isomers IPR=1, IPH=1 /)
C      list of options: IPR,IPH,istop     (Default 0 for all options)
C      In detail:
C      If IPR>0 then the ring spiral subroutine of Fowler and Manopoulus is used.
C         This sub-program catalogues fullerenes with a given number of
C         vertices using the spiral algorithm and a uniqueness test
C         based on equivalent spirals. The required input is IPR.
C         IPR=1 for isolated-pentagon isomers from C60 onwards.
C         IPR=2 for general isomers (note that this generates a lot of output and
C            takes some computer time for fullerenes larger than C60).
C      If IPH=1 then number of distinct Hamiltonian cycles is calculated for
C          every isomer (note this is computer time extensive).
C      If istop=1 program stops after calling this routine.
C      The resulting output is a catalogue of the isomers found containing
C         their idealized point groups, canonical spirals, and NMR patterns
C         (see ref.2).
C
C 6) Option for producing coordinates for fullerene graphs (Schlegel diagrams).
C      &Graph options /      (e.g. &Graph IG=1, ISO1=1, ISO2=3, ISO3=7 /)
C      list of options: IG,ISO1,ISO2,ISO3,PS
C      In detail:
C      If IG>0 Use Program Schlegel for generating the Schlegel projection.
C         IG=1 Use the perspective projection method.
C         IR=2 Use the cone projection method.
C      If ISO1=0 Use the input coordinates for the construction of 
C               the Schlegel diagram.
C      If ISO1.ne.0 Specifying the ring center, edge or
C              vertex through which the z-axis goes at the top of
C              the projection (under construction). This rotates
C              the fullerene around the origin of points into the
C              right position for the Schlegel diagram. Since 3 atoms
C              uniquely define the ring, two the edge, and one the vertex,
C              the input is three integers with the corresponding
C              atom numbers  IO1, IO2, and IO3, i.e. for these values
C                    1 0 0 is a vertex with atom 1 on the top of the z-axis;
C                    7 8 0 is an edge between atom 7 and 8 and z-axis goes
C                          the middle of the bond;
C                    12 13 50 is a ring (either pentagon or hexagon determined
C                          by the program) and z-axis goes through the center
C                          of the ring;
C              NB: You can run the program first with 0 0 0, check the positions
C                  and run it again with the appropriate values.
C              NB2: For IG=1 the input (if chosen) requires to be a ring, i.e.
C                  IO1, IO2 and IO3 are required, otherwise they are chosen by the
C                  program using the largest z-coordinate.
C      If PS.eq.0. Projection angle of 45 degrees is chosen (default value)
C              otherwise to be used as an input parameter in the cone projection
C              method. The angle is the cone angle from the point of projection 
C              to the projection plane, which touches the point with the smallest 
C              z-value (opposite the projection point). Note that angle is reset 
C              to 45 deg if chosen larger than 89.
C              In the case of the perspective projection PS is the distance
C              between the focal point and the ring center underneath, Note this
C              is chosen by the program, but if nonzero this parameter has to be
C              carefully chosen. The picture produced gives a good idea if the
C              ParamS is correctly chosen or not. Note that the ring centers get
C              an extra boost of 10% in the scaling factors such that they appear
C              more in the center of the rings produced by the Schlegel projection.
C              
C---------------------------------------------------------------------- 

C    Set the dimensions for the distance matrix
      PARAMETER (natom=900)      !  Change natom if RAM is not sufficient
      PARAMETER (nat11=natom*11)  
      PARAMETER (msrs=56+1)      !  Size of Schlegel output matrix
      PARAMETER (natom2=natom*natom)
      PARAMETER (natomL=(natom*(natom-1))/2)
      PARAMETER (Nfaces=natom/2+2)
      PARAMETER (NSpScale=12)
      PARAMETER (NSpirals=Nfaces*NSpScale)
      PARAMETER (Nedges=3*natom/2)
      PARAMETER (maxit=2000000) 
      DIMENSION CRing5(3,Nfaces),CRing6(3,Nfaces),cmcs(3),CR(3,Nfaces)
      DIMENSION DistMat(natomL),Dist(3,natom),DistCM(3)
      DIMENSION A(NAtom,NAtom),evec(Natom),df(Natom)
      DIMENSION N5MEM(Nfaces,5),N6MEM(Nfaces,6),Iring(Nfaces)
      DIMENSION Icon2(natom2),distP(natom),IDA(Natom,Natom)
      DIMENSION IATOM(natom),IC3(natom,3),Nring(Nfaces)
      DIMENSION NringA(Nedges),NringB(Nedges)
      DIMENSION NringC(Nedges),NringD(Nedges)
      DIMENSION NringE(Nedges),NringF(Nedges)
      DIMENSION IDual(Nfaces,Nfaces)
      DIMENSION Symbol(Nfaces)
      Real*4 TimeX
CG77  CHARACTER CDAT*9,CTIM*8
      CHARACTER CDAT*8,CTIM*10,Zone*5
      CHARACTER*1  Symbol
      CHARACTER*2 El(99)
      CHARACTER*13 routine
      Integer Values(8)
      DATA El/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE','NA',
     1 'MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V','CR',
     1 'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',    
     1 'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',   
     1 'IN','SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',  
     1 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF', 
     1 'TA',' W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO',
     1 'AT','RN','FR','RA','AC','TH','PA',' U','NP','PU','AM','CM',   
     1 'BK','CF','ES'/                                               
      DATA Tol,anglew,Rdist/0.33d0,45.d0,1.391d0/
C     Van der Waals radius of carbon, adjusted approximately to the
C     solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
      DATA RVdWC/1.415d0/
      IN=5
      IOUT=6

C  You might like to comment these 2 lines out 
C  (and same at the end of this routine) or substitute them with your
C  compiler specific option. Next two are g77 options
CG77    CALL Date(CDAT)
CG77    CALL Time(CTIM)
        call date_and_time(CDAT,CTIM,zone,values)
        TIMEX=0.d0
        CALL Timer(TIMEX)
C       WRITE(IOUT,1000) CDAT,CTIM,natom
        WRITE(IOUT,1000) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7),natom
C  INPUT and setting parameters for running the subroutines
      routine='DATAIN       '
      Write(Iout,1008) routine
        CALL Datain(IN,IOUT,NAtom,MAtom,Icart,Iopt,iprintf,IHam,IPR,
     1   ISchlegel,IS1,IS2,IS3,IER,istop,IOPD,iupac,Ipent,iprintham,
     1   IV1,IV2,IV3,icyl,ParamS,TolX,R5,R6,Rdist)

C  Stop if error in input
      If(IER.ne.0) go to 99
C  Only do isomer statistics
      if(istop.ne.0) go to 98

C Options for Input coordinates
      If(Icart-1) 10,20,30 
C  Cartesian coordinates produced for Ih C60
   10 routine='COORDC60     '
      Write(Iout,1008) routine
      CALL CoordC60(Natom,IN,Iout,IAtom,R5,R6,Dist)
      Matom=60
      Do I=1,60
      IAtom(I)=6
      enddo
      Go to 40
C Input Cartesian coordinates for fullerenes
   20 Do J=1,MAtom
      Read(IN,*) IAtom(J),(Dist(I,J),I=1,3)
      enddo
      Go to 40
C Cartesian coordinates produced ring from spiral pentagon list
C currently using the Fowler-Manopoulus algorithm to
C identify P-type eigenvectors and construct the 3D fullerene
   30 Ipent=1
      routine='COORDPENT    '
      Write(Iout,1008) routine
      Do I=1,Matom
      IAtom(I)=6
      enddo
      CALL CoordPent(Natom,NFaces,Nedges,MAtom,IN,Iout,IDA,IDual,
     1 Icart,IV1,IV2,IV3,A,evec,df,Dist,distp,Rdist)
   40 WRITE(IOUT,1001) MAtom,TolX*100.d0

C Some general infos on isomers and spiral routine
C of Fowler and Manopoulus. Set parameter IPR for independent
C pentagon rule as full list beyond C60 is computer time 
C intensive
  98  routine='ISOMERS      '
      Write(Iout,1008) routine
      CALL Isomers(NAtom,NFaces,Nedges,MAtom,IPR,IOUT,
     1 maxit,iprintham,IDA,A)
      if(istop.ne.0) go to 99

C Move carbon cage to Atomic Center
      routine='MOVECM       '
      Write(Iout,1008) routine
      Iprint=1
      Call MoveCM(Natom,Matom,Iout,Iprint,IAtom,Dist,DistCM,El)

C Calculate largest and smallest atom-to-atom diameters
C Also get moment of inertia (to be implemented)
      routine='DIAMETER     '
      Write(Iout,1008) routine
      CALL Diameter(NAtom,MAtom,IOUT,Dist,distp)

C Calculate the distance Matrix and print out distance Matrix
      routine='DISTMATRIX   '
      Write(Iout,1008) routine
      CALL Distmatrix(NAtom,natomL,MAtom,IOUT,iprintf,Iopt,
     1 Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)

C Establish Connectivities
      routine='CONNECT      '
      Write(Iout,1008) routine
      CALL Connect(NAtom,natomL,Natom2,MCon2,MAtom,Ipent,IOUT,
     1 Icon2,IC3,IDA,TolX,DistMat,Rmin)

C Hueckel matrix and eigenvalues
      if(ipent.eq.0) then
      routine='HUECKEL      '
      Write(Iout,1008) routine
      CALL Hueckel(NAtom,MAtom,IOUT,IC3,IDA,A,evec,df)
      endif

C Generate IUPAC name and locate Hamiltonian cycles. 
C Routine written by D. Babic. Note routine
C is called only if IPR>0 as computer time is extensive beyond
C C100 (PN-hard problem). Last routine uses the adjaceny matrix
C to calculate the number of all distinct paths between 
C adjacent vertices
      routine='HAMILTON     '
      Write(Iout,1008) routine
       maxiter=maxit
      if(IHam.gt.1.and.IHam.le.10) then
       maxiter=10**IHam
      endif
      if(IHam.ne.0) then
       if(iupac.ne.0) then
         CALL Hamilton(NAtom,MAtom,Iout,iprintf,maxiter,IC3)
       else
         CALL HamiltonCyc(NAtom,MAtom,maxiter,IDA,Nhamilton)
         WRITE(Iout,1010) Nhamilton
       endif
      endif
      CALL Paths(NAtom,Nedges,MAtom,IOUT,IDA,A,evec,df)

C Establish all closed ring systems
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(NAtom,Nfaces,natomL,Natom2,MCon2,MAtom,IOUT,N5Ring,
     1 N6Ring,IC3,Icon2,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,DistMat)

C Produce the dual of the fullerene
C     if(IOPD.eq.1) then
C     routine='DualFullerene'
C     Write(Iout,1008) routine
C     WRITE(Iout,1033) 
C1033 Format(/1X,'Produce the adjacency matrix for the dual fullerene')
C     CALL DUALFull(NMAX,MAtom,NFaces,IDA,IDual,
C    1 N5Ring,N6Ring,N5MEM,N6MEM)
C     IOPD=0
C     go to 98
C     endif

C Optimize Geometry through force field method
      If(Iopt.ne.0) then
      routine='OPTFF        '
      Write(Iout,1008) routine
      CALL OptFF(Natom,NFaces,MAtom,Iout,IDA,N5Ring,N6Ring,
     1 N5MEM,N6MEM,Dist,Rdist)
      Iprint=0
      Call MoveCM(Natom,Matom,Iout,Iprint,IAtom,Dist,DistCM,El)
      routine='DISTMATRIX   '
      Write(Iout,1008) routine
      CALL Distmatrix(NAtom,natomL,MAtom,IOUT,Iprintf,0,
     1 Dist,DistMat,Rmin,Rmax,VolSphere,ASphere)
      routine='DIAMETER     '
      Write(Iout,1008) routine
      CALL Diameter(NAtom,MAtom,IOUT,Dist,distp)

      endif
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(NAtom,Nfaces,natomL,Natom2,MCon2,MAtom,IOUT,N5Ring,
     1 N6Ring,IC3,Icon2,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,DistMat)

C Calculate the center for each ring system
      routine='RINGC        '
      Write(Iout,1008) routine
      CALL RingC(NAtom,Nfaces,Nedges,NAtom2,Matom,nat11,Iout,iprintf,
     1 N5MEM,N6MEM,N5Ring,N6Ring,NRing,Iring5,Iring6,Iring56,NringA,
     1 NringB,NringC,NringD,NringE,NringF,DIST,CRing5,CRing6)

C Now produce clockwise spiral ring pentagon count a la Fowler and Manopoulus
      if(ipent.eq.0) then
      routine='SPIRALSEARCH '
      Write(Iout,1008) routine
      CALL SpiralSearch(NAtom,Nfaces,Nedges,Nspirals,MAtom,Iout,
     1 Iring5,Iring6,Iring56,NringA,NringB,NringC,NringD,NringE,NringF)
      endif

C Calculate the volume
      routine='VOLUME       '
      Write(Iout,1008) routine
      CALL Volume(NAtom,Nfaces,NAtom2,Matom,Iout,N5MEM,N6MEM,
     1 N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,ASphere,
     2 Atol,VTol,Rmin5,Rmin6,Rmax5,Rmax6)

C Print out Coordinates used as input for CYLview
      if(icyl.ne.0) then
      Open(unit=3,file='cylview.xyz',form='formatted')
      routine='PRINTCOORD   '
      Write(Iout,1008) routine
      WRITE(Iout,1002) 
      WRITE(3,1003) MATOM 
      Do J=1,MAtom
      IM=IAtom(J)      
      Write(3,1007) El(IM),(Dist(I,J),I=1,3)
      enddo
      endif

C Calculate the minimum distance sphere
C     routine='CONVEXHULL'
C     Write(Iout,1008) routine
C     CALL ConvexHull(NAtom,MAtom,Dist,VolumeCH,AreaCH)

C Calculate the minimum covering sphere and volumes
C     MinCovSphere1 contains algorithm 1 and MinCovSphere2 algorithm2
C     MinCovSphere2 is more efficient and contains more useful information
C     CALL MinCovSphere1(NAtom,MAtom,IOUT,Dist,
C    1 Rmin,Rmax,VolSphere,ASphere,Atol,VTol,cmcs,rmcs,RVdWC)
      routine='MINCOVSPHERE2'
      Write(Iout,1008) routine
      CALL MinCovSphere2(NAtom,MAtom,IOUT,Dist,Rmin,Rmax,
     1 VolSphere,ASphere,Atol,VTol,distP,cmcs,rmcs,RVdWC)

C Calculate the minimum distance sphere
      routine='MINDISTSPHERE'
      Write(Iout,1008) routine
      CALL MinDistSphere(NAtom,MAtom,IOUT,Dist,
     1 Rmin,Rmax,distP,cmcs,rmcs)

C Calculate the maximum inner sphere
      routine='MAXINSPHERE'
      Write(Iout,1008) routine
      CALL MaxInSphere(NAtom,MAtom,IOUT,Dist,cmcs,RVdWC)

C Calculate Schlegel diagram
      if(ISchlegel.ne.0) then
      routine='SCHLEGEL     '
      Write(Iout,1008) routine
      if(ISchlegel.eq.2) then
       if(ParamS.le.1.d0.or.ParamS.gt.8.9d1) then
       ParamS=anglew
       WRITE(IOUT,1006) ParamS
       endif
      else
       ParamS=dabs(ParamS)
      endif
      CALL Schlegel(NAtom,Nfaces,Nedges,MAtom,msrs,IOUT,IS1,IS2,IS3,
     1 N5MEM,N6MEM,N5Ring,N6Ring,NRing,Iring,Ischlegel,IC3,Dist,ParamS,
     1 Rmin,TolX,CR,CRing5,CRing6,Symbol)
      endif

C  E N D   O F   P R O G R A M
CG77 99  CALL TIME(CTIM)
  99  call date_and_time(CDAT,CTIM,zone,values)
        WRITE(IOUT,1004) Values(3),Values(2),Values(1),Values(5),
     1    Values(6),Values(7)
      CALL Timer(TIMEX)
      Hours=TIMEX/3.6d3
      WRITE(IOUT,1009) TIMEX,Hours
      Close(unit=3)
C Formats 
 1000 FORMAT(
     1  1X,'*********************************************************',
     1 /1X,'**      P R O G R A M   F U L L E R E N E              **',
     1 /1X,'**  Fortran Program for the topological analysis       **',
     1 /1X,'**    of Fullerenes C(20+2H) (H: hexagonal faces)      **',
     1 /1X,'**  Written by P.Schwerdtfeger and J.Avery             **',
     1 /1X,'**    with routines from Fowler, Manopoulus and Babic  **',
     1 /1X,'**  Massey University,  Auckland,  New Zealand         **',
     1 /1X,'**  First version                     08/06/10         **',
     1 /1X,'**  Version 4.0                       05/03/12         **',
     1 /1X,'*********************************************************',
CG77 1 /1X,'DATE: ',A9,10X,'TIME: ',A8,/1X,'Limited to ',I6,' Atoms',
     1 /1X,'Date: ',I2,'/',I2,'/',I4,10X,'Time: ',I2,'h',I2,'m',I2,'s',
     1 /1X,'Limited to ',I6,' Atoms',
     1 //1X,'Literature (please cite refs.1-4):',/1X,
     1 '1) P. Schwerdtfeger, J. Avery, Topological Aspects of ',
     1 'Fullerenes - A Fortran Program',/9X,'(Version 4.0, Massey ',
     1 'University Albany, Auckland, New Zealand, 2012).',/1X,
     1 '2) P. W. Fowler, D. E. Manopoulus, An Atlas of Fullerenes',
     1 ' (Dover Publ., New York, 2006).',/1X,
     1 '3) D. Babic, Nomenclature and Coding of Fullerenes,',
     1 ' J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).',/1X,
     1 '4) Z. C. Wu, D. A. Jelski, T. F. George, Vibrational Motions of'
     1 ' Buckminsterfullerene, Chem. Phys. Lett. 137, 291-295 (1987).',
     1 /1X,' Further reading:',/1X,
     1 '5) D. E. Manopoulus, P. W. Fowler, Molecular graphs, point '
     1 'groups, and fullerenes,',
     1 ' J. Chem. Phys. 96, 7603-7614 (1992).',/1X,
     1 '6) D. Babic, D. J. Klein, C. H. Sah, Symmetry of fullerenes,',
     1 ' Chem. Phys. Lett. 211, 235-241 (1993).',/1X,
     1 '7) G. B. Adams, M. O Keefe, R. S. Ruoff, Van der Waals Surface'
     1 ' Areas and Volumes of Fullerenes, J. Phys. Chem. 98, 9465-9469'
     1 '(1994).',/1X,
     1 '8) W. O. J. Boo, An Introduction to Fullerene Structures,',
     1 ' J. Chem. Ed. 69, 605-609 (1992).',/1X,
     1 '9) R. Tonner, G. Frenking, M. Lein, P.',
     1 ' Schwerdtfeger, Packed to the Rafters – Filling up C60 with',
     1 ' Rare Gas Atoms,',/9X,'Chem. Phys. Chem. 12, 2081-2084 (2011).',
     1 /1X)
 1001 FORMAT(/1X,'Number of Atoms: ',I4,', and distance tolerance: ',
     1 F12.2,'%')
 1002 FORMAT(/1X,'Input coordinates to be used for program CYLview ',
     1 'by C.Y.Legault:',/1X,'Output written into cylview.xyz')
 1003 FORMAT(I5,/,'Cartesian coordinate input for CYLview')
CG77 1004 FORMAT(1X,124(1H-),/1X,6HTIME: ,A8)
 1004 FORMAT(132(1H-),/1X,'DATE: ',I2,'/',I2,'/',I4,10X,
     1 'TIME: ',I2,'h',I2,'m',I2,'s')
 1006 FORMAT(/1X,'Angle for Schlegel diagram reset to ',
     1 F10.4,' degrees')
 1007 FORMAT(A2,6X,3(E18.12,2X))
 1008 FORMAT(132('-'),/1x,'--> Enter Subroutine ',A13)
 1009 FORMAT(1x,'CPU Seconds: ',F14.2,', CPU Hours: ',F12.5)
 1010 FORMAT(1X,'Number of Hamiltonian cycles: ',I10)
      STOP 
      END

      SUBROUTINE TIMER(TIMEX)
      Real TA(2)
      Call DTIME(TA,time)
      TIMEX=TIME
      RETURN
      END
