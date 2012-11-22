<!DOCTYPE xhtml>
<xhtml>
<head>

<title> Program Fullerene webpage </title>
<!
--------------------------------------------------------------------------
This is the main webpage for the Fullerene program
Two php files are called:

fulleruf.php
fullerdl.php

fulleruf.php processes user information (iname, address and email),
  appends it to a file called fulleruser.txt and sends an email out
  with a password for downloading the fullerene program. The email is
  checked for email injections by validating the input

fullerdl.php receives the password, checks it and if correct
  the user can download the fullerene program
--------------------------------------------------------------------------
>

<!-- Main content. -->
</head>
<td valign="top" >
<body background="data/background.jpg" text=#000000>


<table width="95%" align="center">
<tr>
<td valign="middle" align="center">
<a href="data/C4860.png"><img src="data/C4860.png" width="270" alt="C4860"></a><br><br>
<a href="data/C4050.png"><img src="data/C4050.png" width="270" alt="C4050"></a><br>
<font size="-1"><i>click picture for larger image</i></font>
</td>
<td valign="middle"  align="center">
<h1> <big> Program FULLERENE </big> <br><br>
A Fortran/C++ Program for Constructing Fullerene Structures <br> 
and for Topological Analysis <br><br> </h1>
<h2>by Peter Schwerdtfeger, Lukas Wirz and James Avery<br>
<I>Centre for Theoretical Chemistry and Physics<br>
The Quantum Chemistry and Physics Group<br>
The New Zealand Institute for Advanced Study<br>
Massey University Auckland<br>
0745 Auckland, New Zealand</I></h2>
</td>
<td valign="middle" align="center">
<a href="data/C60.png"><img src="data/C60.png" width="250" alt="C60"></a><br>
<a href="data/C540.png"><img src="data/C540.png" width="270" alt="C540"></a>
</td>
</tr></table><br>


The program creates cartesian coordinates for fullerene isomers and performs a topological/graph theoretical analysis.
The results can be used for plotting 2D/3D fullerene graphs (e.g. Schlegel diagrams) and structures, and serves as a starting 
point for further quantum theoretical treatment. Version 4 now incorporates C++ routines linked to the
original Fortran program using much improved algorithms.
<br><br>

This is an open source code and free for distribution (except for commercial purposes). The following references should be used for citation:<br><br>
1) P. W. Fowler and D. E. Manolopoulos, "An Atlas of Fullerenes"
(Dover Publ., New York, 2006).<br>
2) P. Schwerdtfeger, L. Wirz and J. Avery, "Topological Analysis of Fullerenes - A Fortran/C++
Program", Version 4.0, Massey University, Auckland, New Zealand (2012).<br>
3) D. Babic, "Nomenclature and Coding of Fullerenes", J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).<br>
For further literature see the README file.<br><br>

The program has been tested for gfortran and gc compilers. You need to use the 
Makefile included in fullerene.zip. The executable runs just fine on a PC or Mac running LINUX/UNIX.
The program is written in standard Fortran and C++ and has been checked
for a large number of fullerenes up to 10000 vertices. A number of input 
and output files are included in the fullerene.zip file and a more detailed description 
is found in the User's Manual.<br><br>

<b>Download for Manual: </b><a href="http://ctcp.massey.ac.nz/data/UserManualv4.1.pdf">UserManualv4.1.pdf (23/08/2012)</a> 
<br><br>

<b>Download for Database: </b><a href="http://ctcp.massey.ac.nz/data/database.zip">database.zip (26/07/2012)</a> (compressed 30.8 MB, uncompressed 251.9 MB)
<br><br>

<b>Download for fullerene program:</b><br>
If you are a new user or if your passwd expired, please supply Title (Prof., Dr., Mr., Mrs.), 
First and Last name, Institution etc, and Email address. 
You will receive automatically an email message with the password to download the program. 
Your entry will be stored and kept confidentially in our internal user-file database.<br><br>
<I>Maximum 30 characters allowed per box (except for comment and email)</I><br>

<table>
<form action="fulleruf.php" method="get">
<tr> <td> Title:</td> <td> <input type="text" name="title" size="30" /></td> </tr>
<tr> <td> First name:</td> <td> <input type="text" name="firstname" size="30" /></td> </tr>
<tr> <td> Last name: </td> <td> <input type="text" name="lastname" size="30" /></td> </tr>
<tr> <td> Institution: </td> <td> <input type="text" name="institution" size="30" /></td> </tr>
<tr> <td> Department: </td> <td> <input type="text" name="department" size="30" /></td> </tr>
<tr> <td> City: </td> <td> <input type="text" name="city" size="30" /></td> </tr>
<tr> <td> Country: </td> <td> <input type="text" name="country" size="30" /></td> </tr>
<tr> <td> Email: </td> <td> <input type="text" name="Emailuser" size="50" /> </td> </tr>
<tr> <td> Comment: </td> <td> <input type="text" name="comment" size="50" /> </td> </tr>
<tr> <td> <input type="submit" value="Submit" /> </td> </tr>
</form>
</table> <br><br>

Once you have a password you can download the program here: (Version 4.1, last update 23/08/2012)<br>
<form action="fullerdl.php" method="get">
Password: <input type="password" name="pwd" size="10" /><br> <input type="submit" value="Download" /> 
</form>

<a href="data/C840.png"><img src="data/C840.png" width="750" alt="C840"></a><br><br>
<b>Important steps in the program are:</b><br>
<li>Read in cartesian coordinates for a fullerene.</li>
<li>Construct cartesian coordinates for the Ih isomer of C<SUB>60</SUB> or C<SUB>20</SUB>.
<li>Construct cartesian coordinates using the face-spiral pentagon indices of Fowler and Manolopoulos.</li>
<li>Construct cartesian coordinates from a Goldberg-Coxeter transformation of C<SUB>20</SUB>.</li>
<li>Use matrix-eigenvector method or Tutte embedding to obtain cartesian coordinates.</li>
<li>Goldberg-Coxeter, leapfrog and halma transformations of fullerenes.</li>
<li>Endo-Kroto 2-vertex and Yoshida-Fowler 4- and 6-vertex insertions.</li>
<li>Use program SPIRAL of Fowler and Manolopoulos for face-spiral fullerenes.</li>
<li>Use program HAMILTON of Babic for Hamiltonian cycles and IUPAC Nomenclature.</li>
<li>Establish connectivities between vertices to construct the adjacency matrix.</li>
<li>Perform a Hueckel analysis. This gives you a good hint if the fullerene is open or closed shell.</li>
<li>Identify all closed 5- and 6-ring systems. Perform a topological analysis for ring connectivities.</li>
<li>For IPR fullerenes calculate the enthalpy of formation using Cioslowski's scheme of motifs.</li>
<li>Use Fletcher-Reeves-Polak-Ribiere geometry optimization with analytical gradients for the 
Wu force field and perform geometry optimization.</li>
<li>Calculate the volume of the fullerene through the convex hull or trigonal pyramidal tesselation.</li>
<li>Calculate the surface area.<br>
<li>Calculate the minimum covering sphere (MCS) of the fullerene using the Yildirim algorithm.</li>
<li>Calculate the minimum distance sphere (MDS) and maximum inner sphere (MIS) of the cage molecule.</li>
<li>Measure of distortion from spherical symmetry.</li>
<li>Produce (x,y) coordinates of a fullerene graph (Schlegel diagram) using a variety of different algorithms.</li>
<li>The output file ending with .xyz can be used as input for CYLVIEW, AVOGADRO or PYMOL.</li>
<li>The output file ending with .dat can be used as input for drawing Schlegel diagrams.</li><br>

A database is provided for general isomers up to C<SUB>100</SUB> and for IPR isomers up to
C<SUB>120</SUB> including the number of Hamiltonian cycles. The database can be copied into
the main program folder and can be used to read the ring spiral pentagon indices.
The numbering scheme is identical to that chosen in the Atlas of Fullerenes by Fowler and Manolopoulos,
that is each isomer in the book's appendix can be constructed easily from the database.  
The datafiles are formatted and can easily be read. It is our intension to
extend the isomer list beyond C<SUB>100</SUB>/C<SUB>120</SUB> (without Hamiltonian cycles)
and to link it to the database available from the House of Graphs by Brinkmann and co-workers.<br><br>

<b>Note: </b>The program is under construction and new features are being added.
You can check the date on the right hand side of the download button to see
if you already have the newest version. 
Please report any bug to <a href="mailto:p.a.schwerdtfeger&#64massey.ac.nz">p.a.schwerdtfeger&#64massey.ac.nz</a>.<br><br>
<b>Acknowledgment:</b><br>
PS is indebted to the Alexander von Humboldt Foundation (Bonn) for financial support
in terms of a Humboldt Research Award. The first version of the program was written during
his stay at the Philipps University Marburg. We acknowledge the help of Darko Babic, 
Patrick W. Fowler and David E. Manolopoulos to kindly allow their Fortran routines to be implemented in our program.
<br>

<!-- Main content ends here. -->
</html>
