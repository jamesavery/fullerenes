!
!           P R O G R A M     F U L L E R E N E
!
!----------------------------------------------------------------------
!
!          A PROGRAM FOR TOPOLOGICAL ANALYSIS OF FULLERENES
!    The program creates cartesian coordinates for fullerenes isomers
!        and performs a topological/graph theoretical analysis.
!        The results can be used for plotting fullerene graphs
!     (Schlegel diagrams) and structures, and as a starting point
!            for further quantum theoretical treatment. 
!        Version 4 now incorporates C++ routines linked to the
!         original Fortran program using improved algorithms.
!
!----------------------------------------------------------------------------
!| Important Copyright Message: Hey there is none! If you like to know why  |
!| read Darrel C. Ince1, Leslie Hatton, John Graham-Cumming, Nature 482,    |
!| p.485 (2012). So: You may do whatever you like with this program, but if |
!| you use it and publish data please cite at least the references given    |
!| below. Before you use the program for commercial purposes you nedd to    |
!| contact the authors.                                                     |
!| The book of Fowler and Manopoulus is highly recommended. It helps        |
!| understanding how this program functions (the book is fun to read as     |
!| well. Many of the concepts used in this program can be found in this     | 
!| book. A standard book on graph theory helps as well.                     |
!----------------------------------------------------------------------------
!
!---------------------------------------------------------------------- 
!  R U N N I N G    T H E    P R O G R A M:
!---------------------------------------------------------------------- 
!* The program is LINUX/UNIX based and works fine with gfortran and 
!      gc compilers
!    You need to use the Makefile included in fullerene.zip 
!    The executable "fullerene" runs on a mac intel as:
!      ./fullerene <inp >out
!    If you type      make tests      it runs all the input jobs and
!      puts them into *.out
!    If you type      make clean      all the object files are deleted
!    It currently compiles in a 64 bit version, but you can change to
!      32 bits in the Makefile if necessary
!    All fortran files are in the directory    source
!    All C-files are in the directory          libgraph
!    If you use the database, the file needs to be in the directory
!      where the source and libgraph directories are.
!
!----------------------------------------------------------------------
!  G E N E R A L   D E S C R I P T I O N:
!---------------------------------------------------------------------- 
!    The program is written in standard Fortran and C++ (~12,000 lines)
!    and is LINUX/UNIX based with links to plotting programs.
!    Reason: I am good in old-fashioned Fortran and James is good in C++.
!    Some standard routines from Mathematical Recipies were modified and
!    are used here for matrix diagonalization and geometry optimization.
!
!* Function: To perform a topological analysis of a regular fullerene 
!    (i.e. consisting of pentagons and hexagons) fulfilling Euler's theorem. 
!    The program calculates the volume and surface area of a fullerene 
!    (irregular or not). It further constructs the structure obtaining 
!    cartesian coordinates from the canonical ring spiral pentagon indices
!    through either Tutte embedding or matrix eigenvector methods. 
!    It can also construct the n-th leapfrog fullerne. Note that 
!    there is no unique definition for the volume of a fullerene for
!    nonplanar 5- or 6-rings on the fullerene surface except for the convex
!    hull, but there is no reason why any other definition than the fast
!    tesselation algorithm should be preferred. The Wu force-field and 
!    geometry optimization using a Fletcher-Reeves-Polak-Ribiere
!    minimization with analytical gradients is also implemented, providing 
!    good a initial guess for cartesian coordinates. Also Schlegel diagrams
!    can be produced using various algorithms.
!
!    Note: This program works for any (distorted or not) regular fullerene
!     (i.e. a fullerene of genus 0 consisting of pentagons and hexagons only).
!     The spiral algorithm of Fowler and Manolopoulus is not
!     restricted to starting from a pentagon or to canonical indices.
!     For a general list of fullerenes see "The House of graphs" at 
!     https://hog.grinvin.org/Fullerenes. 
!
!    Lit.: 1) P. Schwerdtfeger, J. Avery, "Topological Analysis of Fullerenes - 
!             A Fortran and C++ Program (Version 4.0)", Massey University Albany, 
!             Auckland, New Zealand (2012).
!          2) P. W. Fowler and D. E. Manopoulus, "An Atlas of Fullerenes" 
!             (Dover Publ., New York, 2006).
!          3) D. Babic, "Nomenclature and Coding of Fullerenes",
!             J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).
!          4) Z. C. Wu, D. A. Jelski, T. F. George, "Vibrational Motions of
!             Buckminsterfullerene", Chem. Phys. Lett. 137, 291-295 (1987).
!
!      Further reading:
!          5) D. E. Manopoulus and P. W. Fowler, "Molecular graphs, point groups, 
!             and fullerenes", J. Chem. Phys. 96, 7603-7614 (1992).
!          6) G. B. Adams, M. O'Keefe, and R. S. Ruoff, "Van der Waals Surface Areas
!             and Volumes of Fullerenes", J. Phys. Chem. 98, 9465-9469 (1994).
!          7) W. O. J. Boo, "An Introduction to Fullerene Structures",
!             J. Chem. Ed. 69, 605-609 (1992).
!          8) D. Babic, D. J. Klein and C. H. Sah, "Symmetry of fullerenes",
!             Chem. Phys. Lett. 211 (1993) 235-241.
!          9) R. Tonner, G. Frenking, M. Lein, and P. Schwerdtfeger, 
!             "Packed to the Rafters - Filling up C60 with Rare Gas Atoms", 
!             Chem. Phys. Chem. 12, 2081-2084 (2011).
!         10) T. Pisanski, B. Plestenjak, A. Graovac, "NiceGraph Program and its 
!             applications in chemistry", Croatica Chemica Acta 68, 283-292 (1995).
!         11) B. Plestenjak, "An algorithm for drawing Schlegel diagrams",
!             unpublished.
!
!      There is a paper in preparation which explains most of the features
!        in this program. You may ask for a preliminary copy.
!
!      If you use the CYLview program, the input file is written into cylview.xyz
!        if specified (xyz format). This file also works with Avogadro, Jmol, Pymol
!        or any other plot program accepting standard xyz format files.
!        CYLview plots molecules, is written by C. Y. Legault and is freely
!        available from http://www.cylview.org/Home.html. Note for using these
!        programs it is important to force-field optimize them, otherwise
!        bonds cannot be idendified if structures are used directly from
!        AME, LME or 3D-TE algorithms (see below). We recommend CYLview, it
!        is more robust and works for the largest fullerenes up to 1000 atoms,
!        where Avogadro has already some difficulties.
!
!      If you use QMGA for plotting fullerene graphs the input file is written
!        in qmga.dat. The program is available at http://qmga.sourceforge.net/
!        In this case you need to cite A. T. Gabriel, T. Meyer, G. Germano, 
!        "Molecular graphics of convex body fluids", J. Chem. Theory Comput. 
!        4, 468-476 (2008). There will be soon another program availble to
!        plot fullerene graphs (Schlegel diagrams).
!
!   -> Many definitions depend on the use of Angstroems, so please use this unit.
!
!   Important steps in the program are:
!--------------------------------------
! - Create the structure of a fullerene:
!    Read in Cartesian coordinates for a fullerene (see files c20.inp to c540.inp), 
!    or construct them for the Ih isomer (files ico.inp, icoExp.inp, icoideal.inp)
!    of C60, or get cartesian coordinates from ring spiral pentagon indices
!    by the Fowler-Manopoulus algorithm also known as the AME algorithm
!    (adjacency matrix eigenvector algorithm) as defined in refs.2 and 5, or Tutte
!    embedding algorithm, 3D-TEA (used by files pentagon1.inp to pentagon25.inp). 
!    Barycenter of the fullerene is set to the origin. See ref.2 for the 
!    construction of cartesian coordinates from ring spiral pentagon indices and 
!    the use of Hueckel P-type eigenvectors. Note it is critical to get the right
!    vectors for the construction of the cartesian coordinates. The 3 P-type
!    vectors may need to be read in (see file pentagon8.inp for such an example).
!    It is important that the end-product is viewed by molecular visualization.
!    program. We recommend CYLview by Claude Legault (see http://www.cylview.org),
!    Avogadro (see http://avogadro.openmolecules.net) or JMol (see 
!    http://jmol.sourceforge.net/). This are all open source codes.
!    For this purpose a file is written out to cylview.xyz to be used as an input 
!    file. If the AME algorithm fails you can use our Tutte embedding algorithm 
!    (which to our opinion is less trouble and may be used as the standard method).
!    Using these coordinates the program calculates the smallest and largest cage 
!    diameters which gives already a measure for distortion from spherical symmetry
!    (Subroutine DIAMETER). It produces the distance matrix if print level is set
!    to high (Subroutine DISTMATRIX). The diameters already indicate if the
!    fullerene is heavily distorted from spherical symmetry.
!
! - Use program SPIRAL of Fowler and Manopoulus for creating a ring spiral.
!    It also produces canonical ring spiral pentagon indices 
!    (Subroutine SPIRALSEARCH) if cartesian coordinate input is chosen.
!
!-  Print all isomers and perform analysis introduced mostly in the book 
!    written by Fowler and Manopoulus (ref.2), e.g. pentagon indices and 
!    pentagon number, hexagon indices and strain parameter, NMR information
!    and number of distinct Hamiltonian cycles if required.
!
! - Perform Hueckel analysis (Subroutine HUECKEL). 
!    This gives you a good hint if the fullerene is open or closed shell.
!    Note: It is known that for fullerenes the Hueckel analysis is not
!      very reliable as hybridization with the C(2s) orbitals occur due
!      to nonplanarity. Hence the sigma-pi separation breaks down.
!      Nevertheless, we adopt  alpha=-0.21 au  and  beta=-0.111 au obtained
!      from the exp. ionization potential (7.58 eV) and excitation energy
!      (3.02 eV) of C60 (note the electron goes into the second LUMO of 
!      t1g symmetry). This gives also orbital energies for C60 in reasonable 
!      agreement with DFT Kohn-Sham orbital energies. The fullerene might
!      also Jahn-Teller distort adopting a singlet ground state instead
!      of one of a higher multiplicity. Such electronic effects are not
!      captured in this program, and one needs to perform a proper quantum
!      theoretical calculation.
!
! - Use program HAMILTON of Babic (ref.3) for Hamiltonian cycles and IUPAC 
!    Nomenclature. The number of Hamiltonian cycles given has been checked 
!    against a second algorithm, so it should work. Note, the number gives 
!    all distinct Hamiltonian cycles and left-right cycles counted as the 
!    same. Although finding all Hamiltonian cycles by the back-track algorithm
!    of Babich used here is a NP-complete problem, it works fine up to 
!    about C100. After that it becomes computationally very demanding. 
!    Note also that the existence of Hamiltonian cycles for fullerenes is 
!    only conjectured, and only for layered fullerenes (e.g. fullerene 
!    nanotubes) it has been proven to exist, our calculations show that
!    they exist for all fullerene isomers up to C100.
!
! - Establish connectivities between atoms either from given cartesian
!    coordinates or from adjacency matrix if ring spiral as input is used:
!    1) Bond between atoms, 2) Vertices (Subroutine CONNECT).
!    Identify all closed 5- and 6-ring systems (Subroutine RING).
!    This routine also determines if Euler's theorem is fulfilled.
!
! - Determine the center for each 5- and 6-ring (Subroutine RINGC)
!    This is required for the trigonal pyramidal tessellation to obtain
!    the volume and surface.  This routine also analyzes all 2- and 3-ring fusions
!    It further gives the Rhagavachari-Fowler-Manoupoulos neighboring pentagon 
!    and hexagon indices as described in the Fowler and Manolopoulos book
!    (ref.2). From the hexagon indices one derives if the fullerene fulfills the
!    IPR or not.
!
! - Fletcher-Reeves-Polak-Ribiere geometry optimization using analytical 
!    gradients for the Wu force field (Subroutine OPTFF).
!    It is very fast, even for C840. Note that the force field optimization
!    might distort the fullerene from the ideal point group symmetry.
!    On the other hand, the construction of the fullerene by using
!    pentagon indices leads to a more spherical arrangement in both
!    algorithms (Fowler-Manoupoulos or Tutte), e.g. barrels instead of
!    nanotubes.
!
! - Calculate the volume of the fullerene by summing over all
!    tetrahedrons spanned by the three vectors (Subroutine VOLUME).
!     CM-CR  (center of cage to the center of ring)
!     CM-CA1 (center of cage to atom 1 in ring)
!     CM-CA2 (center of cage to atom 2 in ring)
!    There are 5 such tetrahedrons in a 5-ring and 6 in a 6-ring
!    Note that CM is already in the origin
!    Let CR=(X1,Y1,Z1) , CA1=(X2,Y2,Z2) , and CA2=(X2,Y2,Z2)
!    Then the volume V for a irregular tetrahedron is given by 
!    the determinant
!
!                                 | X1 Y1 Z1 |
!     V = abs(Vdet)  ,   V =  1/6 | X2 Y2 Z2 |
!                                 | X3 Y3 Z3 |
!
!   Calculate the surface area A and the area/volume ratio (Subroutine VOLUME)
!
!                                        2              2              2
!                             | Y1 Z1 1 |    | Z1 X1 1 |    | X1 Y1 1 |   
!     A = 1/2 d**0.5 ,   d =  | Y2 Z2 1 | +  | Z2 X2 1 | +  | X2 Y2 1 |
!                             | Y3 Z3 1 |    | Z3 X3 1 |    | X3 Y3 1 |
!
!   Note that the ideal C60 coordinates can be constructed from scratch 
!    (Subroutine COORDC60). This routine was constructed to test the program 
!    for the case of an ideal capped icosahedron, where the analytical formula 
!    is well known, i.e.   V=[(125+43*sqrt(5))*R**3]/4
!    and R is the distance between the atoms (all the same)
!    Setting R=Rmin (Rmin is the smallest distance in C60) this gives a
!    lower bound for the volume, while the volume of the covering central sphere 
!    gives the upper bound.
!    For two different bond distances, R5 for the 5-ring and R6 for the 6-ring
!    joining another 6-ring, the volume can be determined as well after some tedious
!    algebraic manipulations:
!    i.e.   V=5[(3+sqrt(5))*(2R5+R6)**3]/12-[(5+sqrt(5))*R5**3]/2
!    For C20 (ideal dodecahedron) we have V=[(15+7sqrt(5))*R5**3]/4
!    For C20 and C60 the result of these formulae are also printed.
!   This method gives sensible results for convex fullerenes.
!
! - Calculate the minimum covering sphere (MCS) of the cage molecule 
!    (Subroutine MINCOVSPHERE): The MCS in m-dimensional space exists, is unique 
!     and can be expressed as a convex combination of at most (m+1) points, hence
!     our algorithm stops when 4 points are left over in the iteration process.
!
!    The problem can be reduced to
!
!    min(c) max(i) || p(I) - Cmcs ||
!
!    where ||..|| is the Euclidian norm in m dimensions and Cmcs is the center 
!    of the MCS.
!
!    Note: The spherical central cover SCC is not the minimum covering sphere MCS
!    (except if all distances from the center of points CM are the same as in the
!    ideal capped icosahedron). The spherical central cover is taken from the
!    CM point with radius Rmax (longest distance to one vertex).
!    The minimum covering sphere is calculated at the end of the
!    program using the algorithm of E. A. Yildirim, SIAM Journal on Optimization
!    Vol. 19(3),1368-1391 (2008) and the test by T. H. Hopp and C. P. Reeve,
!    NIST, US Department of Commerce (1996)'). Note that the much simpler algorithm 
!    by F. Lu and W. He, Global Congress on Intelligent Systems (2009),
!    DOI 10:1109/GCIS:2009:381, simply does not work for more than m+1 points on a 
!    surface or close by as their linear equation becomes linearly dependent.
!    Note also that the function Psi in Yildirim's algorithm is really the function
!    Phi defined earlier in his paper in section 2 and corrected in the SIAM paper.
!    His easier to program algorithm 1 was also tested, but is much slower.
!    If the value given in the iteration as "convergence" is close to
!    zero (equal zero), the iteration stops (if it falls below epsilon).
!    You can change the epsilon parameter in subroutine Sphere.
!    You can also try algorithm 1 of Yildirim through subroutine Sphere1
!    which is included in an file called algorithm1.f (although this file has not
!    been updated and further developed). Note that we changed the first
!    condition in this algorithm by choosing the furthest point from CM.
!    In the final statistics there should be 0 points outside the sphere
!    and at least 1 point on the sphere.
!    At the end the Van der Waals radius of carbon (1.415 Angstroems) is added to the
!    radius of the minimum covering sphere (note input coordinates for this need to be
!    in Angstroems otherwise change the program), and the volume of
!    an ideal fcc solid is calculated. The Van der Waals radius is chosen such that
!    for C60 the solid-state results of P.A.Heiney et al., Phys. Rev. Lett. 66, 2911 (1991)
!    are reproduced. The definition of the distortion parameter D from the MCS or for the
!    the isoperimetric quotient IPQ is
!
!    IPQ=36Pi(V^2/A^3)
!
!    D=[100/(N*Rmin)]* sum(i=1,N) {Rmcs - ||pi-Cmcs|| }    (N=MAtom)
!
! - Calculate the minimum distance sphere (MDS) of the fullerene.
!    The MCS definition for the distortion is biased for the case that few atoms stick 
!    out on a sphere and the MDS measure may be more appropriate for a measure
!    from spherical distortion. The MDS is defined as
!    The problem can be reduced to
!
!    min(Cmds) 1/N sum(i=1,N) | Rmds - || p(I) - Cmds || |
!
!    where ||..|| is the Eucledian norm in m dimensions. Cmds has to lie within the
!    convex hull. The MDS may not be uniquely defined, as there can be many 
!    (even degenerate) local minima, but for most spherical fullerenes it should 
!    just be fine. Analogous to the MCS there will be a measure for distortion 
!    from spherical symmetry.
!
!    D=[100/(N*Rmin)]* sum(i=1,N) | Rmds - || p(I) - Cmds || |
!
! - Calculate the maximum inner sphere (MCS) of the cage molecule
!
!    max(Cmds) min(i) || p(I) - Cmds ||
!
!    The maximum inner sphere is important for evaluating how much space
!    there is in a fullerene for encapsulating atoms and molecules. For
!    this the radius and volume is printed out with the Van der Waals
!    radius of carbon taken off Rmds. 
!
! - Produce the (X,Y) coordinates of a fullerene graph (Subroutine SCHLEGEL).
!    Schlegel projection (SP):
!    Here the points are rotated (if in input I1,I2, and I3 are given) so to
!    put the selected vertex, edge or ring center on top of the z-axis as the
!    point of projection (otherwise the point (0,0,zmax) is chosen with zmax
!    being the point with maximum z-value from the original input coordinates).
!    Points are then sorted in descending order according to their z-values.
!    The circumference for atoms and rings down the z-axis are determined.
!    The Schlegel projection is created giving as output the projected (X,Y)
!    coordinates. The connections between the points are already written out
!    earlier in the output such that the fullerene graph can be drawn.
!    There are two choices for the projection, depending if you choose the
!    outer ring or the center part of the fullerene graph as a starting point:
!    1) The cone projection (CSP), i.e. points are projected out to an enveloping 
!       cone and then down to a plane below the fullerene. The input I1,I2,I3 
!       defines the center of the fullerene graph. The last ring center should 
!       be at the bottom of the fullerene and if detected, will not be projected 
!       out, or if not will have a large scale factor (this center may be ignored 
!       in the drawing). Also, the last points on the outer ring in the fullerene
!       graph are scaled in distance by 1.2 in order to make the outer rings 
!       more visible. This also moves the outer centers within the ring.
!    2) The perspective projection (PSP), i.e. points are projected down a plane 
!       from a set projection point. In this case the input I1,I2,I3 defines the
!       peripheral ring of the Schlegel diagram. 
!    From the output you can easily construct the name of the fullerene. 
!    At the end a rough printout of the fullerene graph is
!    produced. Note that this is o.k for fullerenes up to about C100, beyond it
!    it becomes too crowded and a proper plotting program should be used.
!    Nevertheless, it serves for a first rough picture. 
!    Furthermore, for large fullerenes it becomes critical to correctly set the
!    projection point or point of the cone. If for example the projection
!    point is too far away from the fullerene, edges may cross.
!    Other algorithms for producing Schlegel diagrams:
!     - Tutte graph with linear scaling (2D-TGE-LS)
!     - Spring embedding with barycentric Coulomb repulsion (SE+C)
!     - Pisanski-Plestenjak-Graovac embedding (PPGA)
!     - Kamada-Kawai embedding (2D-KKE, this gives a 2D picture of a 3D structure
!        and has edge crossings).
!    
!-----------------------------------------------------------------------------------
!  G E N E R A L    I N F O R M A T I O N
!-----------------------------------------------------------------------------------
! Input and output files are in the folders  input  and   output  respectively.
! This program has been tested for the ideal capped icosahedron (input file ico.inp)
!   and for many other fullerenes which are found in the following input files:
!       C20 (c20.inp), C24 (c24.inp), C26 (c26.inp), C28 (c28.inp), C30 (c30.inp),
!       C36 (c36.inp), C50 (c50.inp), C60 (c60.inp), C70 (c70.inp), C72 (c72.inp),
!       C74 (c74.inp), C78 (c78.inp), C80 (c80.inp), C92 (c92.inp), C100 (c100.inp),
!       C180 (c180.inp), C320 (c320.inp), and C540 (c540.inp), pentagon1.inp, ...,
!       pentagon25.inp
!   The coordinates are mostly B3LYP aug-cc-pVDZ optimized up to C60, and
!    cc-pVDZ up to C180, and 6-31G for the rest and all for singlet states
!    (except of course for the ones where the pentagon indices input is
!    chosen). Note that for some of the fullerene coordinates the singlet
!    state chosen may not be the electronic ground state.
! Number of atoms is set in NATOMS currently at 900, so change this parameter
!   if you do not have enough RAM, alternatively the distances need to be
!   calculated directly and the DistMat(natomL) Matrix removed 
!   (which I can do if somebody insists). Also, maxit=2000000, which sets the
!   number of isomers or the number of Hamiltonian cycles to this value in order
!   for the program to not run forever.
! NB: The algorithm for locating all 5-and 6-rings might not be the smartest
!      one, but as this program requires only a second or so to run it was
!      not important to find a better algorithm.
! You should also be aware of program fullgen for generating nonisomorphic fullerenes.
!      It is written by Gunnar Brinkmann (Gunnar.Brinkmann@Ugent.be) and can be
!      downloaded from Brendan McKay's website (Australian National University) at
!      http://cs.anu.edu.au/~bdm/plantri/
!
!  This program is under permanent construction. The program has been tested for bugs.
!      Nevertheless, if you have any problem or find a bug please report to:
!                   peter.schwerdtfeger@gmail.com
!  Note: This program has neither been optimized for performance nor style.
!      It is however quite fast even for fullerenes such as C840.
!
!  Not implemented yet and on the to do list (in progress) is:
!      1) A subroutine to fit the minimum outer ellipsoidal cover useful for  
!          rugby ball like fullerenes and close packing of ellipsoids.
!      2) Volume of the convex hull (coming soon).
!      3) Use the Coxeter construction for fullerenes.
!      4) Geometry optimization using the extended Wu-Fowler force field.
!      5) Frequency calculations from the force field optimized geometry.
!      6) Construction of non-ring-spiral isomers using the genus algorithm.
!      7) Symmetry labels for Hueckel orbital energies.
!      8) Extend to non-regular fullerenes of genus 0 (heptagons and squares).
!      9) Extend to non-regular fullerenes of genus 1.
!     10) Symmetrize coordinates to point group symmetry.
!     11) Implement Cioslowski's scheme for enthalpy of formation.
!     12) Restart option for subroutine for Hamiltonian cycle.
!     13) Use of databank for all non-IPR fullerenes uo to C120 and IPR up to
!          C150 (coming soon).
!     14) Produce pictures of Schlegel diagrams and corresponding duals (coming soon).
!     15) Produce fullerene name from Schlegel diagram.
!     16) Use of a databases (coming soon).
!
! For any questions concerning this program please contact P. Schwerdtfeger
!   Centre for Theoretical Chemistry and Physics (CTCP)
!   The New Zealand Institute for Advanced Study (NZIAS) 
!   Massey University Auckland, Bldg.44, Private Bag 102904
!   North Shore City, 0745 Auckland, New Zealand
!   email: peter.schwerdtfeger@gmail.com
!   http://ctcp.massey.ac.nz/   and  http://www.nzias.ac.nz/
!   --> Always open to improvements and suggestions
!
!-----------------------------------------------------------------------------------
!  A C K N O W L E D G M E N T
!-----------------------------------------------------------------------------------
! PS is indebted to the Alexander von Humboldt Foundation (Bonn) for financial support 
! in terms of a Humboldt Research Award, and to both Prof. Gernot Frenking and 
! Dr. Ralf Tonner (Marburg) for support during my extended stay in Marburg where I
! started to write this program. I acknowledge also the help of Darko Babich, Patrick 
! W. Fowler and David E. Manopoulus to allow the free distribution of their Fortran 
! subroutines.
!
!-----------------------------------------------------------------------------------
!  I N P U T:
!-----------------------------------------------------------------------------------
! Input (Either Namelist or in Free Format):   (use Angstroem for distances)
!
! 1) Text card of 80 characters (A80 format)
!    (you can put as many text cards in as you want, the main title (first card)
!     is printed out extra
!
! 2) Input to create cartesian coordinates and main flags for the program
!    &Coord options /       (e.g. &Coord IC=20, IOPT=1, R6=1.42 /)
!    list of options: NA,IC,IP,IV1,IV2,IV3,ixyz,ichk,leap,isonum,IPRC,TolR,R5,R6,xyzname
!    NA= Number of Atoms (Default: 60)
!    IC= Flag for construction of cartesian coordinates (Default: 0)
!    IP= Print option (Default: 0)
!    IV1= Number for Hueckel P-type eigenvector for AME algorithm (Default: 2)
!    IV2= Number for Hueckel P-type eigenvector for AME algorithm (Default: 3)
!    IV3= Number for Hueckel P-type eigenvector for AME algorithm (Default: 4)
!    ixyz= Flag for producing input file for CYLview, Avogadro or others in standard
!      xyz format (Default: 0)
!    isonum= Isomer number according to the scheme of Fowler and Manopoulus (Default: 0)
!      If IC=2, 3 or 4 and isonum not zero, than pentagon indices are taken from the
!      isomer list contained in a database (see below). There are two databases, one
!      for the general isomers (IPRC=2) and one for the IPR isomers (IPRC=1), the
!      definition is similar to the IPR parameter below (Default: 2).
!    xyzname (max 20 characters) file name if ixyz.ne.0 (default: cylview.xyz)
!    TolR= Tolerance in % (Default: 33)
!    R5= pentagon bond distance (Default: 1.455)
!    R6= hexagon  bond distance (Default: 1.391)
!    In detail:
!      If IC = 0 No coordinate input required, cordinates are constructed 
!              for the IPR isomer of C60
!           In this case only one card is read in:
!              R5,R6        (arbitrary units, e.g. Angstroms)
!              R5: Bond lengths in the pentagons 
!              R6: Bond length of the bonds connecting hexagons
!              If R5=R6 chosen then the ideal capped icosahedron is obtained
!      If IC = 1 Cartesian Coordinates expected as input
!         In this case N lines with    Z, X, Y, Z  in free format are expected.
!         (Z= Nuclear Charge, X,Y,C Cartesian Coordinates for Atom).
!         NB: Z is not really needed, but you can copy Gaussian output
!          directly into the input file
!      If IC = 2 or 3 Cartesian Coordinates are created from pentagon 
!             ring spiral list. Extra input required (free format):
!
!           IRSP(I),I=1,12 
!
!           R6 is taken as the smallest bond distance in the fullerene and IP(I)
!           identify the locations of the pentagons as described in detail
!           in ref.2, i.e. IP is the pentagon ring spiral numbering scheme. 
!           Note this only works for ring spiral fullerenes as described in 
!           detail by P. W. Fowler and D. E. Manopoulus. Use the canonical
!           pentagon ring indices if possible (transformation to the canonical
!           from should work as well).
!           IC=2: AME algorithm using P-type eigenvectors produced from the 
!            adjacency matrix:
!            If problem with eigenvectors are found to construct the
!            cartesian coordinates, i.e. the identification of P-type
!            eigenvectors, three integer values IV1, IV2, IV3 can be specified
!            identifying the eigenvectors to be chosen. pentagon8.inp is such an example.
!            In this case a severe warning occurs which means you should carefully
!            check the eigenvectors used and cartesian coordinates produced.
!            Otherwise coordinates are obtained which are useless. This is more
!            often the case as you might expect. 
!          IC=3: Same as IC=2 but the Laplacian matrix is used instead of the
!            adjacency matrix (LME algorithm). 
!          IC=4: Tutte embedding (3D-TE) algorithm. This should always work,
!            although the initial fullerene might be too spherical. But this
!            algorithm is easier, and (in theory) should never fail.
!           Examples are given in the input files starting with 'pentagon'.
!           Please use Angstroems.
!          IC=5: Goldberg-Coxeter construction of fullerene
!      If IP>0 larger output produced, i.e. the full distance matrix, all
!          Hamiltonian cycles and all 3-ring connections.
!      if leap=n than the n-th leapfrog fullerene is generated.
!      Connectivities are found for atoms with distances between
!         R6   and   R6*(1+TolR/100)   if cartesian coordinate input is chosen.
!      If TolR=0. default value of 33% is used. 
!         NB: If this parameter is set at a value too large, unwanted connectivities
!          are produced resulting in smaller polygons. This parameter
!          should reflect the maximum deviation in distance from the
!          smallest distance found.
!
! 3) Option for force-field optimization:
!      &Opt options /        (e.g. &Opt Iopt=1 /)
!      list of options: Iopt,ftol,WuR5,WuR6,WuA5,WuA6,WufR,WufA,fCoulomb
!      Iopt= Flag for force-field optimization (Default: 0)
!      In detail:
!       If Iopt=1  then fullerene is optimized using the force field method
!         of Wu et al within a Fletcher-Reeves-Polak-Ribiere algorithm:
!         Z.C.Wu, D.A.Jelski, T.F.George, Chem. Phys. Lett. 137, 291-295 (1987).
!         Note that there are more sophisticated force fields available,
!         but for these general parameters for fullerenes are not
!         yet available, and the Wu force field does the job to create good initial
!         cartesian coordinates for further refinement using more sophisticated
!         QM methods. Note that a converged energy much greater than zero implies
!         that the set distances and angles in the Wu force field cannot be
!         reached for all atoms and rings. For further reading see: 
!         A.Ceulemans, B.C.Titeca, L.F.Chibotaru, I.Vos, P.W.Fowler, 
!         J. Phys. Chem. A 105, 8284-8295 (2001).
!       If Iopt=2  Preoptimize with input force field, then optimize with Wu
!         force field. This is especially usefull for fcoulomb input (see below).  
!         NB: Avogadro has a more sophisticated force-field which you can try out.
!       ftol: The convergence tolerance on the function value is input as ftol
!         (Default: 5.0E-8)
!       WuR5,WuR6,WuA5,WuA6,WufR,WufA: Force field parameters for Wu force field for
!         distances R5, R6, angles A5, A6, force constants for distance WufR and
!         angles WufA (see paper by Wu et al. for details).
!         Defaults: WuR5=1.455, WuR6=1.391, WuA5=1.08d2, WuA6=1.2d2, WufR=1.d6,
!                   WufA=1.d5
!       If fCoulomb>0. then add an additional repulsive Coulomb forct from the 
!         barycenter to the atoms      (Default:  fCoulomb=0.d0)
!         This is extremely useful for an initial geometry optimization to keep
!         the cage convex if for example the Tutte construction leads to not
!         a good guess of the initial structure. In such cases fCoulomb=1.d2
!         is a good choice.
!
! 4) Option for calculating Hamiltonian cycles and IUPAC numbers 
!      &Hamilton options /        (e.g. &Hamilton IHam=1 IUPAC=1 /)
!      list of options: IHam,IUPAC
!      In detail:
!      If IHam>0 Then Subroutine HAMILTON is called.     (Default: 0)
!         IHam=1 Routine will stop after 1 million Hamiltonian cycles are found
!         IHam=1000000000 Program runs forever and prints if IHam is reached
!      If IUPAC=1 IUPAC numbers are produced.   (Default: 0) 
!         Note only with this option together with IP=1 in &Coord input 
!         all Hamiltonian cycles are printed out. IP=0 only gives the
!         best cycle (see Babic, ref.3).
!         IUPAC=0 goes into a fast subroutine and prints only the total number of
!         Hamiltonian cycles.
!
! 5) Option for producing list of isomers and properties.
!      &Isomers options /        (e.g.&Isomers IPR=1, IPH=1 /)
!      list of options: IPR,IPH,IStop,IChk,chkname (Default 0 for all options
!                                                   and 'checkpoint' for chkname)
!      In detail:
!      If IPR>0 then the ring spiral subroutine of Fowler and Manopoulus is used.
!         This sub-program catalogues fullerenes with a given number of
!         vertices using the spiral algorithm and a uniqueness test
!         based on equivalent spirals. The required input is IPR.
!         IPR=1 for isolated-pentagon isomers from C60 onwards.
!         IPR=2 for general isomers (note that this generates a lot of output and
!            takes some computer time for fullerenes larger than C60).
!      If IPH=1 then number of distinct Hamiltonian cycles is calculated for
!          every isomer (note this is computer time extensive).
!      If istop=1 program stops after calling this routine.
!      If IChk=1  Restart: Isomer list is continued from previous output file 
!          called 'checkpoint' as default if not otherwise give in chkname.
!          This is a restart option from a previous run which terminated. Note
!          that the new output file does not contain the previous one.
!      The resulting output is a catalogue of the isomers found containing
!         their idealized point groups, canonical spirals, and NMR patterns
!         (see ref.2).
!
! 6) Option for producing coordinates for fullerene graphs (Schlegel diagrams).
!      &Graph options /      (e.g. &Graph IG=1, ISO1=1, ISO2=3, ISO3=7 /)
!      list of options: IG,ISO1,ISO2,ISO3,PS,SCALE,SCALEPPG
!      We recommend IG= 1, 2, 4 or 7
!      In detail:
!      If IG>0 Use Program Schlegel for generating fullerene graphs
!         IG=1 Use the perspective projection method (PSP).
!         IG=2 Use the cone projection method (CSP).
!         IG=3 Produce the Tutte graph (2D-TE)
!         IG=4 Produce the Tutte graph and perform linear scaling (2D-TE-LS)
!              (scale factor can be read in).
!         IG=5 Produce the Tutte graph and perform spring embedding optimization
!              (rij-r0)**2 with r0 set to 2.0 (2D-SE).
!         IG=6 Starting from the Tutte graph perform spring + repulsive Coulomb
!              force embedding optimization (2D-SE+C). The repulsive force is taken
!              from the barycenter to the vertex. r0 is set to 2.0. The force
!              constants are set such that the graph looks nice.
!         IG=7 Starting from the Tutte graph perform a Pisanski-Plestenjak-Graovac 
!              embedding PPGE) optimization (called Schlegel by the authors)
!              Note: This is much faster than their simulated annealing algorithm.
!         IG=8 Starting from the Tutte graph perform a Kamada-Kawai embedding 
!              optimization using the distance matrix (2D-KKE).
!              Note: This gives a 2D picture of a 3D structure, thus
!               it is not related to a Schlegel diagram and has edge crossings.
!      If ISO1=0 Use the input coordinates for the construction of 
!               the Schlegel diagram.
!      If ISO1.ne.0 Specifying the ring center, edge or
!              vertex through which the z-axis goes at the top of
!              the projection (under construction). This rotates
!              the fullerene around the origin of points into the
!              right position for the Schlegel diagram. Since 3 atoms
!              uniquely define the ring, two the edge, and one the vertex,
!              the input is three integers with the corresponding
!              atom numbers  IO1, IO2, and IO3, i.e. for these values
!                    1 0 0 is a vertex with atom 1 on the top of the z-axis;
!                    7 8 0 is an edge between atom 7 and 8 and z-axis goes
!                          the middle of the bond;
!                    12 13 50 is a ring (either pentagon or hexagon determined
!                          by the program) and z-axis goes through the center
!                          of the ring;
!              NB: You can run the program first with 0 0 0, check the positions
!                  and run it again with the appropriate values.
!              NB2: For IG=1 the input (if chosen) requires to be a ring, i.e.
!                  IO1, IO2 and IO3 are required, otherwise they are chosen by the
!                  program using the largest z-coordinate.
!      If Scale.ne.0. Scaling factor for linear scaling of Tutte graph (default  2.5)
!              2D coordinated are scaled by 1.+.5*Scale*(rmin-r)/rmin
!              r is the distance of the vertex from the barycenter
!              rmin is the smallest distance of the vertex from the barycenter 
!               belonging to the outer circumferencing 5- or 6-ring
!      If ScalePPG.ne.0. Scaling factor for exponential in Pisanski-Plestenjak-Graovac
!              algorithm     (default 1.0)
!      If PS.ne.0. 
!              - Projection angle in degrees is chosen (default 45 degrees)
!              to be used as an input parameter in the cone projection
!              method. The angle is the cone angle from the point of projection 
!              to the projection plane, which touches the point with the smallest 
!              z-value (opposite the projection point). Note that angle is reset 
!              to 45 deg if chosen larger than 89.
!              - In the case of the perspective projection PS is the distance
!              between the focal point and the ring center underneath, Note this
!              is chosen by the program, but if nonzero this parameter has to be
!              carefully chosen. The picture produced gives a good idea if the
!              ParamS is correctly chosen or not. Note that the ring centers get
!              an extra boost of 10% in the scaling factors such that they appear
!              more in the center of the rings produced by the Schlegel projection.
!-----------------------------------------------------------------------------------
!  F U L L E R E N E     I S O M E R     D A T A B A S E
!-----------------------------------------------------------------------------------
! A database is provided for general isomers up tp C100 and for IPR isomers up to
! C120 including the number of Hamiltonian cycles. The database can be copied into
! the main program folder and can be used to read the ring spiral pentagon indices.
! The numbering scheme is exactly that chosen in the book by Fowler and Manopoulus
! (ref.2), that is each isomer in the books appendix can be constructed easily
! from the database. An example is given in the input file   pentagon13.inp.
! The datafiles are formatted and can easily be read. It is our intension to
! extend the isomer list beyond C100/C12 (without Hamiltonian cycles). New lists
! will be available on our website. Note the determination of the number of
! distinct Hamiltonian cycles is NP-complete and beyond 100 (120 for IPR) 
! computationally too demanding. The longest file for our database ran for 3
! months on a single processor.
! Note: the directory database needs to be in the same directories as source or 
! libgraph.
!-----------------------------------------------------------------------------------

      PROGRAM Fullerene
      IMPLICIT REAL*8 (A-H,O-Z)
C    Set the dimensions for the distance matrix
      PARAMETER (natom=5000)      !  Change natom if RAM is not sufficient
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
      DIMENSION DistMat(natomL),Dist(3,natom),DistCM(3),Dist2D(2,natom)
      DIMENSION A(NAtom,NAtom),evec(Natom),df(Natom)
      DIMENSION forceWu(9),forceWuP(9)
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
      CHARACTER*20 xyzname
      CHARACTER*20 chkname
      Character TEXTINPUT*80
      CHARACTER*3 GROUP
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
      leapspiral=0
      Write(Iout,1008) routine
        CALL Datain(IN,IOUT,NAtom,MAtom,Icart,Iopt,iprintf,IHam,IPR,
     1  IPRC,ISchlegel,IS1,IS2,IS3,IER,istop,leap,iupac,Ipent,iprintham,
     1  IGC1,IGC2,IV1,IV2,IV3,icyl,ichk,isonum,ParamS,TolX,R5,R6,Rdist,
     1  scales,scalePPG,ftolP,forceWu,forceWuP,xyzname,chkname,
     1  TEXTINPUT)

C  Stop if error in input
      If(IER.ne.0) go to 99
C  Only do isomer statistics
      if(istop.ne.0) go to 98

C Options for Input coordinates
      go to (10,20,30,30,30,31) Icart+1
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
     1 Icart,IV1,IV2,IV3,isonum,IPR,IPRC,A,evec,df,Dist,Dist2D,distp,
     1 Rdist,GROUP)
      CALL Chiral(Iout,GROUP)
      Go to 40

C Goldberg-Coxeter construction of fullerens
   31 CALL GoldbergCoxeter(Natom,NFaces,Nedges,MAtom,IGC1,IGC2,
     1 IN,Iout,IDA,A,evec,df,Dist,Dist2D,distp,Rdist)

   40 WRITE(Iout,1001) MAtom,TolX*100.d0

C Some general infos on isomers and spiral routine
C of Fowler and Manopoulus. Set parameter IPR for independent
C pentagon rule as full list beyond C60 is computer time 
C intensive
  98  routine='ISOMERS      '
      Write(Iout,1008) routine
      CALL Isomers(NAtom,NFaces,Nedges,MAtom,IPR,IOUT,
     1 maxit,iprintham,ichk,IDA,A,chkname)
      if(istop.ne.0) go to 99

C Move carbon cage to Atomic Center
  999 routine='MOVECM       '
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

C Produce the nth leapfrog of the fullerene
      if(leap.gt.0) then
      routine='Leapfrog'
      Write(Iout,1008) routine
      CALL Leapfrog(NAtom,MAtom,Iout,leap,LeapErr,IDA,
     1 A,evec,df,Dist,Dist2D,distp,Rdist)
      Do I=1,Matom
       IAtom(I)=6
      enddo
      leap=0
      ipent=1
      leapspiral=1
      if(MAtom.gt.100) IHam=0
      if(LeapErr.eq.0) go to 999
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
      if(IHam.gt.1.and.IHam.le.9) then
       maxiter=10**IHam
      endif
      if(IHam.ne.0) then
       if(iupac.ne.0) then
         CALL Hamilton(NAtom,MAtom,Iout,iprintf,maxiter,IC3)
       else
         CALL HamiltonCyc(NAtom,MAtom,maxiter,Iout,nbatch,IDA,Nhamilton)
         WRITE(Iout,1010) Nhamilton
         if(nbatch.ne.0) WRITE(Iout,1014)
       endif
      endif
      CALL Paths(NAtom,Nedges,MAtom,IOUT,IDA,A,evec,df)

C Establish all closed ring systems
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(NAtom,Nedges,Nfaces,natomL,Natom2,MCon2,MAtom,IOUT,
     1 N5Ring,N6Ring,IC3,Icon2,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,
     1 DistMat)

C Optimize Geometry through force field method
      If(Iopt.ne.0) then
      routine='OPTFF        '
      Write(Iout,1008) routine
      if(Iopt.eq.2) then
       ftol=1.d-4
      else
       ftol=ftolP
      endif
 333  CALL OptFF(Natom,NFaces,MAtom,Iout,IDA,N5Ring,N6Ring,
     1 N5MEM,N6MEM,Dist,Rdist,ftol,forceWu)
      if(Iopt.eq.2) then
       ftol=ftolP
      Write(Iout,1003)
       Do I=1,9
        forceWu(I)=forceWuP(I)
       enddo
       Iopt=1
       go to 333
      endif
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

C Rings
      routine='RING         '
      Write(Iout,1008) routine
      CALL Ring(NAtom,Nedges,Nfaces,natomL,Natom2,MCon2,MAtom,IOUT,
     1 N5Ring,N6Ring,IC3,Icon2,N5MEM,N6MEM,Rmin5,Rmin6,Rmax5,Rmax6,
     1 DistMat)

C Calculate the center for each ring system
      routine='RINGC        '
      Write(Iout,1008) routine
      CALL RingC(NAtom,Nfaces,Nedges,NAtom2,Matom,nat11,Iout,iprintf,
     1 N5MEM,N6MEM,N5Ring,N6Ring,NRing,Iring5,Iring6,Iring56,NringA,
     1 NringB,NringC,NringD,NringE,NringF,DIST,CRing5,CRing6)

C Now produce clockwise spiral ring pentagon count a la Fowler and Manopoulus
      if(ipent.eq.0.or.leapspiral.ne.0) then
      routine='SPIRALSEARCH '
      Write(Iout,1008) routine
      CALL SpiralSearch(NAtom,Nfaces,Nedges,Nspirals,MAtom,Iout,Iring5,
     1 Iring6,Iring56,NringA,NringB,NringC,NringD,NringE,NringF,GROUP)
      CALL Chiral(Iout,GROUP)
      endif

C Calculate the volume
      routine='VOLUME       '
      Write(Iout,1008) routine
      CALL Volume(NAtom,Nfaces,NAtom2,Matom,Iout,N5MEM,N6MEM,
     1 N5Ring,N6Ring,DIST,CRing5,CRing6,VolSphere,ASphere,
     2 Atol,VTol,Rmin5,Rmin6,Rmax5,Rmax6)

C Print out Coordinates used as input for CYLview
      if(icyl.ne.0) then
      Open(unit=3,file=xyzname,form='formatted')
      routine='PRINTCOORD   '
      Write(Iout,1008) routine
      WRITE(Iout,1002) xyzname 
      if(MAtom.lt.100) WRITE(3,1011) MAtom,MAtom,TEXTINPUT 
      if(MAtom.ge.100.and.MAtom.lt.1000) 
     1 WRITE(3,1012) MAtom,MAtom,TEXTINPUT 
      if(MAtom.ge.1000.and.MAtom.lt.10000) 
     1 WRITE(3,1013) MAtom,MAtom,TEXTINPUT
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
      CALL MinDistSphere(NAtom,MAtom,IOUT,Dist,distP,cmcs,rmcs)

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
     1 N5MEM,N6MEM,N5Ring,N6Ring,NRing,Iring,Ischlegel,IC3,IDA,Dist,
     1 ParamS,Rmin,TolX,scales,scalePPG,CR,CRing5,CRing6,Symbol)
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
     1  1X,' ________________________________________________________ ',
     1 /1X,'|                                                        |',
     1 /1X,'|          P R O G R A M   F U L L E R E N E             |',
     1 /1X,'|    Fortran/C++ Program for the topological analysis    |',
     1 /1X,'|      of regular fullerenes (pentagons and hexagons)    |',
     1 /1X,'|    Written by Peter Schwerdtfeger and James Avery      |',
     1 /1X,'|      with routines from Fowler, Manopoulus and Babic   |',
     1 /1X,'|    Massey University,  Auckland,  New Zealand          |',
     1 /1X,'|    First version: 1.0:               from 08/06/10     |',
     1 /1X,'|    This  version: 4.0, last revision from 08/04/12     |',
     1 /1X,'|________________________________________________________|',
CG77 1 /1X,'DATE: ',A9,10X,'TIME: ',A8,/1X,'Limited to ',I6,' Atoms',
     1 //1X,'Date: ',I2,'/',I2,'/',I4,10X,'Time: ',I2,'h',I2,'m',I2,'s',
     1 /1X,'Limited to ',I6,' Atoms',
     1 /1X,'For citation when running this program use:',/1X,
     1 '1) P. Schwerdtfeger, J. Avery, Topological Aspects of ',
     1 'Fullerenes - A Fortran Program',/9X,'(Version 4.0, Massey ',
     1 'University Albany, Auckland, New Zealand, 2012).',/1X,
     1 '2) P. W. Fowler, D. E. Manopoulus, An Atlas of Fullerenes',
     1 ' (Dover Publ., New York, 2006).',/1X,
     1 '3) D. Babic, Nomenclature and Coding of Fullerenes,',
     1 ' J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).',/1X,
     1 'See README file for further literature and input instructions ',
     1 'concerning this program')
 1001 FORMAT(/1X,'Number of Atoms: ',I4,', and distance tolerance: ',
     1 F12.2,'%')
 1002 FORMAT(/1X,'Input coordinates to be used for program CYLview ',
     1 'by C.Y.Legault:',/1X,'Output written into ',A20)
CG77 1004 FORMAT(1X,124(1H-),/1X,6HTIME: ,A8)
 1003 FORMAT(/1X,'Reoptimization using the original Wu force field')
 1004 FORMAT(132(1H-),/1X,'DATE: ',I2,'/',I2,'/',I4,10X,
     1 'TIME: ',I2,'h',I2,'m',I2,'s')
 1006 FORMAT(/1X,'Angle for Schlegel diagram reset to ',
     1 F10.4,' degrees')
 1007 FORMAT(A2,6X,3(F15.6,2X))
 1008 FORMAT(132('-'),/1x,'--> Enter Subroutine ',A13)
 1009 FORMAT(1x,'CPU Seconds: ',F15.2,', CPU Hours: ',F13.5)
 1010 FORMAT(1X,'Number of Hamiltonian cycles: ',I10)
 1011 FORMAT(I5,/,'C',I2,'/  ',A80)
 1012 FORMAT(I5,/,'C',I3,'/  ',A80)
 1013 FORMAT(I5,/,'C',I4,'/  ',A80)
 1014 FORMAT(3X,'(Add to this batches from previous cycles!)')
      STOP 
      END

      SUBROUTINE TIMER(TIMEX)
      Real TA(2)
      Call DTIME(TA,time)
      TIMEX=TIME
      RETURN
      END
