<!DOCTYPE xhtml>
<xhtml>
<head>

<title> Program Fullerene webpage </title>
<!
--------------------------------------------------------------------------
This is the main webpage for the Fullerene program 
Last update: October 23, 2012
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
<h1> <big> Program FULLERENE <br>
Version 4.2.2</big> <br><br>
A Fortran/C++ Program for Fullerene Structures <br> 
and Topological Analysis <br><br> </h1>
<h2>by Peter Schwerdtfeger, Lukas Wirz and James Avery<br>
<I>Centre for Theoretical Chemistry and Physics<br>
The Quantum Chemistry and Physics Group<br>
The New Zealand Institute for Advanced Study<br>
Massey University Auckland<br>
0745 Auckland, New Zealand (Middle Earth)</I></h2>
</td>
<td valign="middle" align="center">
<a href="data/C60.png"><img src="data/C60.png" width="250" alt="C60"></a><br>
<a href="data/C540.png"><img src="data/C540.png" width="270" alt="C540"></a>
</td>
</tr></table><br>


The program creates cartesian coordinates for fullerene isomers and performs a topological/graph theoretical analysis.
The results can be used for plotting 2D/3D fullerene graphs (e.g. Schlegel diagrams and 3D structures), and serves as a starting 
point for further quantum theoretical treatment. Version 4 now incorporates C++ routines linked to the
original Fortran program using much improved algorithms. Version 4.2 is now more robust than 4.1 or 4.0, and contains
a more extended version of topological indicators, an extended Wu force field using dihedral angles, access to both
the Yoshida and House of Graphs databases, and an improved Tutte embedder. Our own database has been extended up to 
C<SUB>122</SUB>, and up to C<SUB>170</SUB> for IPR isomers. Version 4.2.2 is a slightly improved version and contains some bug fixes.
<br><br>

This is an open-source code and free for distribution (except for commercial purposes). The following references should be used for citation:<br><br>
1) P. W. Fowler and D. E. Manolopoulos, "An Atlas of Fullerenes"
(Dover Publ., New York, 2006).<br>
2) P. Schwerdtfeger, L. Wirz, and J. Avery, "Topological Analysis of Fullerenes - 
A Fortran and C++ Program", Version 4.2.2, Massey University Albany, Auckland, New Zealand, 2012.<br>
3) D. Babic, "Nomenclature and Coding of Fullerenes", J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).<br>
For further literature see the User's Manual.<br><br>

The program has been tested for the gnu compiler collection (gfortran and g++). 
You need to use the Makefile included in fullerene.zip. The executable runs 
just fine on a PC or Mac running LINUX/UNIX.
The program is written in standard Fortran and C++ and has been checked
for a large number of fullerenes up to 10,000 vertices. A number of input 
and output files are included in the fullerene.zip file and a more detailed description 
is found in the User's Manual.<br><br>

<b>Download for Manual: </b><a href="http://ctcp.massey.ac.nz/data/UserManualv4.2.2.pdf">UserManualv4.2.2.pdf (22/11/2012)</a> (2.3 MB)
<br><br>

<b>Download for Database1: General fullerenes C<SUB>20</SUB>-C<SUB>122</SUB> </b><a href="http://ctcp.massey.ac.nz/data/ALL.zip"> ALL.zip (01/11/2012)</a> (compressed 198.9 MB, uncompressed 1.92 GB)
<br>
<b>Download for Database2: IPR fullerenes C<SUB>60</SUB>-C<SUB>170</SUB> </b><a href="http://ctcp.massey.ac.nz/data/IPR.zip"> IPR.zip (01/11/2012)</a> (compressed 169.2 MB, uncompressed 1.86 GB)
<br>
<b>Download for Database3: Yoshida selected fullerenes C<SUB>20</SUB>-C<SUB>720</SUB></b><a href="http://ctcp.massey.ac.nz/data/Yoshida.zip"> Yoshida.zip</a> (compressed 4.4 MB, uncompressed 15.4 MB)
<br>
<br>

<b>Download for fullerene program:</b><br>
If you are a new user or if your password expired, please supply Title (Prof., Dr., Mr., Mrs.), 
First and Last name, Institution etc, and Email address. 
You will receive automatically an email message with the password to download the program. 
Your entry will be stored and kept confidentially in our internal user-file database.
This has the advantage that you will receive an email message if a new major upgrade becomes available.<br><br>
<I>Maximum 30 characters allowed per box (except for comment and email)</I><br>

<table>
<form action="fulleruf.php" method="get">
<tr> <td> Title:</td> <td> <select name="Title" />
<option value="">Title...</option>
<option value="Prof. Dr.">Prof. Dr.</option>
<option value="Prof.">Prof.</option>
<option value="Dr.">Dr.</option>
<option value="Ms.">Mr.</option>
<option value="Mr.">Mrs.</option>
</select>
</td> </tr>
<tr> <td> First name:</td> <td> <input type="text" name="firstname" size="30" /></td> </tr>
<tr> <td> Last name: </td> <td> <input type="text" name="lastname" size="30" /></td> </tr>
<tr> <td> Institution: </td> <td> <input type="text" name="institution" size="30" /></td> </tr>
<tr> <td> Department: </td> <td> <input type="text" name="department" size="30" /></td> </tr>
<tr> <td> City: </td> <td> <input type="text" name="city" size="30" /></td> </tr>
<tr> <td> Country: </td> <td> <select name="country" />
<option value="">Country...</option>
<option value="Afganistan">Afghanistan</option>
<option value="Albania">Albania</option>
<option value="Algeria">Algeria</option>
<option value="American Samoa">American Samoa</option>
<option value="Andorra">Andorra</option>
<option value="Angola">Angola</option>
<option value="Anguilla">Anguilla</option>
<option value="Antigua &amp; Barbuda">Antigua &amp; Barbuda</option>
<option value="Argentina">Argentina</option>
<option value="Armenia">Armenia</option>
<option value="Aruba">Aruba</option>
<option value="Australia">Australia</option>
<option value="Austria">Austria</option>
<option value="Azerbaijan">Azerbaijan</option>
<option value="Bahamas">Bahamas</option>
<option value="Bahrain">Bahrain</option>
<option value="Bangladesh">Bangladesh</option>
<option value="Barbados">Barbados</option>
<option value="Belarus">Belarus</option>
<option value="Belgium">Belgium</option>
<option value="Belize">Belize</option>
<option value="Benin">Benin</option>
<option value="Bermuda">Bermuda</option>
<option value="Bhutan">Bhutan</option>
<option value="Bolivia">Bolivia</option>
<option value="Bonaire">Bonaire</option>
<option value="Bosnia &amp; Herzegovina">Bosnia &amp; Herzegovina</option>
<option value="Botswana">Botswana</option>
<option value="Brazil">Brazil</option>
<option value="British Indian Ocean Ter">British Indian Ocean Ter</option>
<option value="Brunei">Brunei</option>
<option value="Bulgaria">Bulgaria</option>
<option value="Burkina Faso">Burkina Faso</option>
<option value="Burundi">Burundi</option>
<option value="Cambodia">Cambodia</option>
<option value="Cameroon">Cameroon</option>
<option value="Canada">Canada</option>
<option value="Canary Islands">Canary Islands</option>
<option value="Cape Verde">Cape Verde</option>
<option value="Cayman Islands">Cayman Islands</option>
<option value="Central African Republic">Central African Republic</option>
<option value="Chad">Chad</option>
<option value="Channel Islands">Channel Islands</option>
<option value="Chile">Chile</option>
<option value="China">China</option>
<option value="Christmas Island">Christmas Island</option>
<option value="Cocos Island">Cocos Island</option>
<option value="Colombia">Colombia</option>
<option value="Comoros">Comoros</option>
<option value="Congo">Congo</option>
<option value="Cook Islands">Cook Islands</option>
<option value="Costa Rica">Costa Rica</option>
<option value="Cote DIvoire">Cote D'Ivoire</option>
<option value="Croatia">Croatia</option>
<option value="Cuba">Cuba</option>
<option value="Curaco">Curacao</option>
<option value="Cyprus">Cyprus</option>
<option value="Czech Republic">Czech Republic</option>
<option value="Denmark">Denmark</option>
<option value="Djibouti">Djibouti</option>
<option value="Dominica">Dominica</option>
<option value="Dominican Republic">Dominican Republic</option>
<option value="East Timor">East Timor</option>
<option value="Ecuador">Ecuador</option>
<option value="Egypt">Egypt</option>
<option value="El Salvador">El Salvador</option>
<option value="Equatorial Guinea">Equatorial Guinea</option>
<option value="Eritrea">Eritrea</option>
<option value="Estonia">Estonia</option>
<option value="Ethiopia">Ethiopia</option>
<option value="Falkland Islands">Falkland Islands</option>
<option value="Faroe Islands">Faroe Islands</option>
<option value="Fiji">Fiji</option>
<option value="Finland">Finland</option>
<option value="France">France</option>
<option value="French Guiana">French Guiana</option>
<option value="French Polynesia">French Polynesia</option>
<option value="French Southern Ter">French Southern Ter</option>
<option value="Gabon">Gabon</option>
<option value="Gambia">Gambia</option>
<option value="Georgia">Georgia</option>
<option value="Germany">Germany</option>
<option value="Ghana">Ghana</option>
<option value="Gibraltar">Gibraltar</option>
<option value="Great Britain">Great Britain</option>
<option value="Greece">Greece</option>
<option value="Greenland">Greenland</option>
<option value="Grenada">Grenada</option>
<option value="Guadeloupe">Guadeloupe</option>
<option value="Guam">Guam</option>
<option value="Guatemala">Guatemala</option>
<option value="Guinea">Guinea</option>
<option value="Guyana">Guyana</option>
<option value="Haiti">Haiti</option>
<option value="Hawaii">Hawaii</option>
<option value="Honduras">Honduras</option>
<option value="Hong Kong">Hong Kong</option>
<option value="Hungary">Hungary</option>
<option value="Iceland">Iceland</option>
<option value="India">India</option>
<option value="Indonesia">Indonesia</option>
<option value="Iran">Iran</option>
<option value="Iraq">Iraq</option>
<option value="Ireland">Ireland</option>
<option value="Isle of Man">Isle of Man</option>
<option value="Israel">Israel</option>
<option value="Italy">Italy</option>
<option value="Jamaica">Jamaica</option>
<option value="Japan">Japan</option>
<option value="Jordan">Jordan</option>
<option value="Kazakhstan">Kazakhstan</option>
<option value="Kenya">Kenya</option>
<option value="Kiribati">Kiribati</option>
<option value="Korea North">Korea North</option>
<option value="Korea Sout">Korea South</option>
<option value="Kuwait">Kuwait</option>
<option value="Kyrgyzstan">Kyrgyzstan</option>
<option value="Laos">Laos</option>
<option value="Latvia">Latvia</option>
<option value="Lebanon">Lebanon</option>
<option value="Lesotho">Lesotho</option>
<option value="Liberia">Liberia</option>
<option value="Libya">Libya</option>
<option value="Liechtenstein">Liechtenstein</option>
<option value="Lithuania">Lithuania</option>
<option value="Luxembourg">Luxembourg</option>
<option value="Macau">Macau</option>
<option value="Macedonia">Macedonia</option>
<option value="Madagascar">Madagascar</option>
<option value="Malaysia">Malaysia</option>
<option value="Malawi">Malawi</option>
<option value="Maldives">Maldives</option>
<option value="Mali">Mali</option>
<option value="Malta">Malta</option>
<option value="Marshall Islands">Marshall Islands</option>
<option value="Martinique">Martinique</option>
<option value="Mauritania">Mauritania</option>
<option value="Mauritius">Mauritius</option>
<option value="Mayotte">Mayotte</option>
<option value="Mexico">Mexico</option>
<option value="Midway Islands">Midway Islands</option>
<option value="Moldova">Moldova</option>
<option value="Monaco">Monaco</option>
<option value="Mongolia">Mongolia</option>
<option value="Montserrat">Montserrat</option>
<option value="Morocco">Morocco</option>
<option value="Mozambique">Mozambique</option>
<option value="Myanmar">Myanmar</option>
<option value="Nambia">Nambia</option>
<option value="Nauru">Nauru</option>
<option value="Nepal">Nepal</option>
<option value="Netherland Antilles">Netherland Antilles</option>
<option value="Netherlands">Netherlands (Holland, Europe)</option>
<option value="Nevis">Nevis</option>
<option value="New Caledonia">New Caledonia</option>
<option value="New Zealand">New Zealand</option>
<option value="Nicaragua">Nicaragua</option>
<option value="Niger">Niger</option>
<option value="Nigeria">Nigeria</option>
<option value="Niue">Niue</option>
<option value="Norfolk Island">Norfolk Island</option>
<option value="Norway">Norway</option>
<option value="Oman">Oman</option>
<option value="Pakistan">Pakistan</option>
<option value="Palau Island">Palau Island</option>
<option value="Palestine">Palestine</option>
<option value="Panama">Panama</option>
<option value="Papua New Guinea">Papua New Guinea</option>
<option value="Paraguay">Paraguay</option>
<option value="Peru">Peru</option>
<option value="Phillipines">Philippines</option>
<option value="Pitcairn Island">Pitcairn Island</option>
<option value="Poland">Poland</option>
<option value="Portugal">Portugal</option>
<option value="Puerto Rico">Puerto Rico</option>
<option value="Qatar">Qatar</option>
<option value="Republic of Montenegro">Republic of Montenegro</option>
<option value="Republic of Serbia">Republic of Serbia</option>
<option value="Reunion">Reunion</option>
<option value="Romania">Romania</option>
<option value="Russia">Russia</option>
<option value="Rwanda">Rwanda</option>
<option value="St Barthelemy">St Barthelemy</option>
<option value="St Eustatius">St Eustatius</option>
<option value="St Helena">St Helena</option>
<option value="St Kitts-Nevis">St Kitts-Nevis</option>
<option value="St Lucia">St Lucia</option>
<option value="St Maarten">St Maarten</option>
<option value="St Pierre &amp; Miquelon">St Pierre &amp; Miquelon</option>
<option value="St Vincent &amp; Grenadines">St Vincent &amp; Grenadines</option>
<option value="Saipan">Saipan</option>
<option value="Samoa">Samoa</option>
<option value="Samoa American">Samoa American</option>
<option value="San Marino">San Marino</option>
<option value="Sao Tome & Principe">Sao Tome &amp; Principe</option>
<option value="Saudi Arabia">Saudi Arabia</option>
<option value="Senegal">Senegal</option>
<option value="Seychelles">Seychelles</option>
<option value="Sierra Leone">Sierra Leone</option>
<option value="Singapore">Singapore</option>
<option value="Slovakia">Slovakia</option>
<option value="Slovenia">Slovenia</option>
<option value="Solomon Islands">Solomon Islands</option>
<option value="Somalia">Somalia</option>
<option value="South Africa">South Africa</option>
<option value="Spain">Spain</option>
<option value="Sri Lanka">Sri Lanka</option>
<option value="Sudan">Sudan</option>
<option value="Suriname">Suriname</option>
<option value="Swaziland">Swaziland</option>
<option value="Sweden">Sweden</option>
<option value="Switzerland">Switzerland</option>
<option value="Syria">Syria</option>
<option value="Tahiti">Tahiti</option>
<option value="Taiwan">Taiwan</option>
<option value="Tajikistan">Tajikistan</option>
<option value="Tanzania">Tanzania</option>
<option value="Thailand">Thailand</option>
<option value="Togo">Togo</option>
<option value="Tokelau">Tokelau</option>
<option value="Tonga">Tonga</option>
<option value="Trinidad &amp; Tobago">Trinidad &amp; Tobago</option>
<option value="Tunisia">Tunisia</option>
<option value="Turkey">Turkey</option>
<option value="Turkmenistan">Turkmenistan</option>
<option value="Turks &amp; Caicos Is">Turks &amp; Caicos Is</option>
<option value="Tuvalu">Tuvalu</option>
<option value="Uganda">Uganda</option>
<option value="Ukraine">Ukraine</option>
<option value="United Arab Erimates">United Arab Emirates</option>
<option value="United Kingdom">United Kingdom</option>
<option value="United States of America">United States of America</option>
<option value="Uraguay">Uruguay</option>
<option value="Uzbekistan">Uzbekistan</option>
<option value="Vanuatu">Vanuatu</option>
<option value="Vatican City State">Vatican City State</option>
<option value="Venezuela">Venezuela</option>
<option value="Vietnam">Vietnam</option>
<option value="Virgin Islands (Brit)">Virgin Islands (Brit)</option>
<option value="Virgin Islands (USA)">Virgin Islands (USA)</option>
<option value="Wake Island">Wake Island</option>
<option value="Wallis &amp; Futana Is">Wallis &amp; Futana Is</option>
<option value="Yemen">Yemen</option>
<option value="Zaire">Zaire</option>
<option value="Zambia">Zambia</option>
<option value="Zimbabwe">Zimbabwe</option>
</select>
</td> </tr>
<tr> <td> Email: </td> <td> <input type="text" name="Emailuser" size="50" /> </td> </tr>
<tr> <td> Comment: </td> <td> <input type="text" name="comment" size="50" /> </td> </tr>
<tr> <td> <input type="submit" value="Submit" /> </td> </tr>
</form>
</table> <br><br>

Once you have a password you can download the program here: (Version 4.2.2, last update 22/11/2012)<br>
<form action="fullerdl.php" method="get">
Password: <input type="password" name="pwd" size="10" /><br> <input type="submit" value="Download" /> 
</form>

<a href="data/C840.png"><img src="data/C840.png" width="750" alt="C840"></a><br><br>
<b>Important steps in the program are:</b><br>
<li>Read in cartesian coordinates for a fullerene.</li>
<li>Construct cartesian coordinates for the Ih isomer of C<SUB>60</SUB> or C<SUB>20</SUB>.
<li>Construct cartesian coordinates using the face-spiral pentagon indices of Fowler and Manolopoulos.</li>
<li>Construct cartesian coordinates from a Goldberg-Coxeter transformation of C<SUB>20</SUB>.</li>
<li>Use matrix-eigenvector or Tutte embedding algorithms to obtain cartesian coordinates.</li>
<li>Goldberg-Coxeter, leapfrog and halma transformations of fullerenes.</li>
<li>Endo-Kroto 2-vertex, Yoshida-Fowler 4- and 6-vertex and Brinkmann-Fowler 6-vertex insertions.</li>
<li>Use of program SPIRAL of Fowler and Manolopoulos for face-spiral fullerenes.</li>
<li>Use of program HAMILTON of Babic for Hamiltonian cycles and IUPAC Nomenclature.</li>
<li>Establish connectivities between vertices to construct the adjacency matrix.</li>
<li>Perform a Hueckel analysis. This gives you a good hint if the fullerene is open or closed shell.</li>
<li>Topological Indices: Wiener, Estrada, Balaban ...</li>
<li>Identify all pentagons and hexagons and perform a topological analysis for ring connectivities.</li>
<li>Calculate the enthalpy of formation using Martin's or Cioslowski's scheme of motifs.</li>
<li>Use Fletcher-Reeves-Polak-Ribiere geometry optimization with analytical gradients for the 
Wu- and extended Wu force fields and perform a geometry optimization.</li>
<li>Calculate the volume of the fullerene through the convex hull or trigonal pyramidal tesselation.</li>
<li>Calculate the surface area.<br>
<li>Calculate the minimum covering sphere (MCS) of the fullerene using the Yildirim algorithm.</li>
<li>Calculate the minimum distance sphere (MDS) and maximum inner sphere (MIS).</li>
<li>Measure of distortion from spherical symmetry and convexity.</li>
<li>Produce (x,y) coordinates of a fullerene graph (Schlegel diagram) using a variety of different algorithms.</li>
<li>The output file ending with .xyz or .cc1 can be used as input for CYLVIEW, AVOGADRO, VMD or PYMOL.</li>
<li>The output file ending with .dat can be used as input for drawing Schlegel diagrams.</li><br>

A database is provided for general isomers up to C<SUB>122</SUB> and for IPR isomers up to
C<SUB>170</SUB>, including the number of Hamiltonian cycles for general isomers up to C<SUB>100</SUB> 
and for IPR isomers up to C<SUB>120</SUB>. The database can be copied into
the main program folder and can be used to read the ring spiral pentagon indices.
The numbering scheme is identical to that chosen in the "Atlas of Fullerenes" by Fowler and Manolopoulos,
that is each isomer in the book's appendix can be constructed easily from the database.  
The data files are formatted and can easily be read on its own. It is our intention to
extend the isomer list beyond C<SUB>122</SUB>/C<SUB>170</SUB> (without Hamiltonian cycles).
A link to the  "House of Graphs" database by Brinkmann and co-workers is already implemented, 
and the database can be downloaded from <a href="http://hog.grinvin.org/Fullerenes">House of Graphs Database</a>.<br><br>
A link to the Yoshida database is also implemented and the database can be downloaded from
<a href="http://www.jcrystal.com/steffenweber/gallery/Fullerenes/Fullerenes.html">Yoshida Database</a>.

<b>Note: </b>The program is under construction and new features are being added.
You can check the date on the right hand side of the download button to see
if you already have the newest version. 
Please report any bug to <a href="mailto:p.a.schwerdtfeger&#64massey.ac.nz">p.a.schwerdtfeger&#64massey.ac.nz</a>.<br>
The program is in use in the following countries: 
Brazil, China, Colombia, France, Germany, Japan, India, Itali, Iran, New Zealand, Poland, Russia, Taiwan, Ukraine, United Kingdom, USA.
<br><br>
<b>Acknowledgment:</b><br>
PS is indebted to the Alexander von Humboldt Foundation (Bonn) for financial support
in terms of a Humboldt Research Award. The first version of the program was written during
his stay at the Philipps University Marburg. We acknowledge the help of Darko Babic (Zagreb), 
Patrick W. Fowler (Sheffield) and David E. Manolopoulos (Oxford) to kindly allow their Fortran routines to be 
modified and implemented in our program system. We also thank Prof. Ottorino Ori (Rome) for fruitful discussions.

<br>

<!-- Main content ends here. -->
</html>