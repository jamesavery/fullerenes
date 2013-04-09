<!DOCTYPE xhtml>
<xhtml>
<head>

<title> Program Fullerene webpage </title>
<!
--------------------------------------------------------------------------
This is the main webpage for the Fullerene program 
Last update: March 19, 2013
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
<a href="data/C250.png"><img src="data/C250.png" width="270" alt="C250"></a>
</td>
<td valign="middle"  align="center">
<h1> <big> Program FULLERENE <br>
Version 4.4</big> <br><br>
A Fortran/C++ Program for Creating Fullerene Structures <br> 
and for Performing Topological Analyses <br><br> </h1>
<h2>by Peter Schwerdtfeger, Lukas Wirz and James Avery<br>
<I>Centre for Theoretical Chemistry and Physics<br>
The Quantum Chemistry and Physics Group<br>
The New Zealand Institute for Advanced Study<br>
Massey University Auckland<br>
0745 Auckland, New Zealand (Middle Earth)</I></h2><br>
<a href="data/C4050.png"><img src="data/C4050.png" width="350" alt="C4050"></a><br>
<font size="-1"><i>click picture for larger image</i></font></td>
<td valign="middle" align="center">
<a href="data/C60.png"><img src="data/C60.png" width="270" alt="C60"></a><br>
<a href="data/C540.png"><img src="data/C540.png" width="300" alt="C540"></a>
</td>
</tr></table><br>


The program generates cartesian/internal coordinates for fullerene isomers and performs a topological/graph theoretical analysis.
The results can be used for plotting 2D/3D fullerene graphs (e.g. Schlegel diagrams and 3D structures), and serves as 
a starting point for further quantum theoretical treatment. Version 4 incorporated C++ routines for the first time
linked to the original Fortran program using much improved algorithms. Version 4.4 is more robust than previous versions, 
and contains a more extended evaluation of topological indicators, a harmonic oscillator force field including dihedral angles, 
vibrational frequencies (although the force constants have not been adjusted yet), general Goldberg-Coxeter transformations
of C<SUB>20</SUB> (now all icosahedral fullerenes can be constructed), access to both the Yoshida and House of Graphs 
databases, an improved Tutte embedder, a correct ring spiral count, a generalized ring spiral algorithm with jumps, 
tight lower and upper bounds for Hamiltonian cycles, and a neighbourhood spiral ring search algorithm. The fullerene database 
has been brought into a more compact form and extended up to C<SUB>150</SUB>, and up to C<SUB>200</SUB> for IPR isomers, 
and is now more compact.  
<br><br>

This is an open-source code and free for distribution (except for commercial purposes). Please cite the following paper if you use this program for publishing data:<br><br>
1) P. Schwerdtfeger, L. Wirz, and J. Avery, "Program Fullerene - A Software Package for Constructing and Analyzing 
Structures of Regular Fullerenes", Version 4.4, J. Comput. Chem., in press. DOI10.1002/jcc.23278.<br>
We also recommend to cite the following book (many of the concepts used can be found here):<br>
2) P. W. Fowler and D. E. Manolopoulos, "An Atlas of Fullerenes"
(Dover Publ., New York, 2006).<br>
The Hamiltonian cycle part of the program uses the code by Babic:<br>
3) D. Babic, "Nomenclature and Coding of Fullerenes", J. Chem. Inf. Comput. Sci. 35, 515-526 (1995).<br>
For further literature see the User's Manual.<br><br>

The program has been tested for the gnu compiler collection (gfortran and g++). 
You need to use the Makefile included in fullerene.zip. The executable runs 
just fine on a PC or Mac running LINUX/UNIX.
The program is written in standard Fortran and C++ and has been checked
for a large number of fullerenes up to 35,000 vertices. A number of input 
and output files are included in the fullerenev4.4.zip file and a more detailed description 
is found in the User's Manual.<br><br>

<b>Download for Manual: </b><a href="http://ctcp.massey.ac.nz/data/UserManualv4.4.pdf">UserManualv4.4.pdf (20/01/2013)</a> (2.4 MB)
<br><br>

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
<option value="Afghanistan">Afghanistan</option>
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
<option value="Curacao">Curacao</option>
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
<option value="Korea South">Korea South</option>
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
<option value="Namibia">Namibia</option>
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
<option value="Philipines">Philippines</option>
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
<option value="United Arab Emirates">United Arab Emirates</option>
<option value="United Kingdom">United Kingdom</option>
<option value="United States of America">United States of America</option>
<option value="Uruguay">Uruguay</option>
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

Once you have a password you can download the program here: (Version 4.4, last update 19/01/2013)<br>
<form action="fullerdl.php" method="get">
Password: <input type="password" name="pwd" size="10" /><br> <input type="submit" value="Download" /> 
</form><br><br>

<b>Database files:</b> The database is provided for general isomers up to C<SUB>150</SUB> and for IPR isomers up to
C<SUB>200</SUB>, including the number of Hamiltonian cycles for general isomers up to C<SUB>110</SUB> 
and for IPR isomers up to C<SUB>122</SUB>. The database can be copied into
the main program folder and can be used to read the ring spiral pentagon indices.
The numbering scheme is identical to that chosen in the "Atlas of Fullerenes" by Fowler and Manolopoulos,
that is each isomer in the book's appendix can be constructed easily from the database.  
The data files are gzipped. before usage of a specific fullerene database you need to unzip the specific file.
Once printed by the program the data can easily be read. It is our intention to
extend the isomer list beyond C<SUB>150</SUB>/C<SUB>200</SUB> in due time (without Hamiltonian cycles).
A link to the  "House of Graphs" database by Brinkmann and co-workers is already implemented, 
and the database can be downloaded from <a href="http://hog.grinvin.org/Fullerenes">House of Graphs Database</a>.<br><br>
A link to the Yoshida database is also implemented and the database can be downloaded from
<a href="http://www.jcrystal.com/steffenweber/gallery/Fullerenes/Fullerenes.html">Yoshida Database</a>.<br><br>

<b>Download for DatabaseY (Yoshida selected fullerenes):</b><br>
C<SUB>20</SUB>-C<SUB>720</SUB></b><a href="http://ctcp.massey.ac.nz/data/Yoshida.zip"> Yoshida.zip (31/10/2012)</a> (compressed 4.4 MB, uncompressed 15.4 MB)
<br>
The following database files are more compact than the ones used in the previous version 4.3. They are now obsolete. Once the new database has been downloaded,
the gzipped files require about 2.3 GB of space. In order to save disk space, only gunzip the files which are to be used in the calculation. 
Start downloading here for the various fullerene database files (in tar format):<br>
<b>Download for DatabaseAllC20-C200 (General fullerenes):</b><br> 
C<SUB>020</SUB>-C<SUB>110</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C020-110.tar"> databaseAll.C020-110.tar (20/03/2013)</a> (compressed 60.4 MB)<br>
C<SUB>112</SUB>-C<SUB>120</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C112-120.tar"> databaseAll.C112-120.tar (20/03/2013)</a> (compressed 66.3 MB)<br>
C<SUB>122</SUB>-C<SUB>130</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C122-130.tar"> databaseAll.C122-130.tar (20/03/2013)</a> (compressed 144.8 MB)<br>
C<SUB>132</SUB>-C<SUB>140</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C132-140.tar"> databaseAll.C132-140.tar (20/03/2013)</a> (compressed 298.7 MB)<br>
C<SUB>142</SUB>-C<SUB>146</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C142-146.tar"> databaseAll.C142-146.tar (20/03/2013)</a> (compressed 303.0 MB)<br>
C<SUB>148</SUB>-C<SUB>150</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseAll.C148-150.tar"> databaseAll.C148-150.tar (03/04/2013)</a> (compressed 278.3 MB)<br>
<b>Download for DatabaseIPRC60-C200 (IPR fullerenes):</b><br>
C<SUB>060</SUB>-C<SUB>160</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C060-160.tar"> databaseIPR.C060-160.tar (20/03/2013)</a> (compressed 45.8 MB)<br>
C<SUB>162</SUB>-C<SUB>170</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C162-170.tar"> databaseIPR.C162-170.tar (20/03/2013)</a> (compressed 67.1 MB)<br>
C<SUB>172</SUB>-C<SUB>180</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C172-180.tar"> databaseIPR.C172-180.tar (20/03/2013)</a> (compressed 145.9 MB)<br>
C<SUB>182</SUB>-C<SUB>186</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C182-186.tar"> databaseIPR.C182-186.tar (20/03/2013)</a> (compressed 153.5 MB)<br>
C<SUB>188</SUB>-C<SUB>190</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C188-190.tar"> databaseIPR.C188-190.tar (20/03/2013)</a> (compressed 144.3 MB)<br>
C<SUB>192</SUB>-C<SUB>194</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C192-194.tar"> databaseIPR.C192-194.tar (20/03/2013)</a> (compressed 188.0 MB)<br>
C<SUB>196</SUB>-C<SUB>198</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C196-198.tar"> databaseIPR.C196-198.tar (20/03/2013)</a> (compressed 243.6 MB)<br>
C<SUB>200</SUB> <a href="http://ctcp.massey.ac.nz/data/databaseIPR.C200.tar"> databaseIPR.C200.tar (20/03/2013)</a> (compressed 146.8 MB)<br>


<table width="95%" align="center">
<tr>
<td valign="left" align="left">
<a href="data/C840.png"><img src="data/C840.png" width="750" alt="C840"></a><br><br>
<b>Important steps in the program are:</b><br>
<li>Read in cartesian coordinates for a fullerene.</li>
<li>Establish connectivities between vertices to construct the adjacency matrix.</li>
<li>Construct cartesian coordinates for the Ih isomer of C<SUB>60</SUB> or C<SUB>20</SUB>.
<li>Construct cartesian coordinates using the ring-spiral pentagon indices of Fowler and Manolopoulos.</li>
<li>Construct cartesian coordinates for non-spiral fullerenes using generalized ring-spiral pentagon indices with jumps.</li>
<li>Construct cartesian coordinates from a Goldberg-Coxeter transformation of C<SUB>20</SUB>.</li>
<li>Use matrix-eigenvector or Tutte embedding algorithms to obtain cartesian coordinates.</li>
<li>Goldberg-Coxeter leapfrog and halma transformations of a fullerene.</li>
<li>Endo-Kroto 2-vertex, Yoshida-Fowler 4- and 6-vertex and Brinkmann-Fowler 6-vertex insertions.</li>
<li>Use of program SPIRAL of Fowler and Manolopoulos for face-spiral fullerenes.</li>
<li>Find fullerenes in the neighbourhood of given ring-spiral pentagon indices.</li>
<li>Use of program HAMILTON of Babic for Hamiltonian cycles and IUPAC Nomenclature.</li> 
<li>Tight upper and lower bounds for Hamiltonian cycles</li>
<li>Perform a Hueckel analysis. This gives you a good hint if the fullerene is open or closed shell.</li>
<li>Topological Indices: Wiener, Estrada, Balaban, Szeged, and many more.</li>
<li>Identify all pentagons and hexagons and perform a topological analysis for ring connectivities.</li>
<li>Calculate the enthalpy of formation using Martin's or Cioslowski's scheme of motifs.</li>
<li>Use Fletcher-Reeves-Polak-Ribiere geometry optimization with analytical gradients for the 
harmonic oscillator force field including dihedrals.</li>
<li>Get Hessian for the force field and calculate the vibrational frequencies.</li>
<li>Calculate the volume of the fullerene through the convex hull or trigonal pyramidal tesselation.</li>
<li>Calculate the surface area.</li>
<li>Calculate the minimum covering sphere (MCS) of the fullerene using the Yildirim algorithm.</li>
<li>Calculate the minimum distance sphere (MDS) and maximum inner sphere (MIS).</li>
<li>Measure of distortion from spherical symmetry and convexity.</li>
<li>Produce (x,y) coordinates of a fullerene graph (Schlegel diagram) using a variety of different algorithms.</li>
<li>Produce Schlegel diagram, with Hamiltonian cycle or the dual graph included, in .tex format.</li>
<li>The output file ending with .xyz or .cc1 can be used as input for CYLVIEW, AVOGADRO, VMD or PYMOL.</li>
<li>The output file ending with .dat can be used as input for drawing Schlegel diagrams.</li><br>
</td>
<td valign="right" align="right">
<a href="data/Obelix.png"><img src="data/Obelix.png" width="350" alt="Obelix"></a><br><br>
</td>
</tr></table><br>


<b>Note: </b>The program is under construction and new features are being added.
You can check the date on the right hand side of the download button to see
if you already have the newest version. The next intended version 5 will undergo major restructuring.
Please report any bug to <a href="mailto:p.a.schwerdtfeger&#64massey.ac.nz">p.a.schwerdtfeger&#64massey.ac.nz</a>.<br>
The program is used in over 30 countries worldwide.<br><br>
<b>Acknowledgment:</b><br>
PS is indebted to the Alexander von Humboldt Foundation (Bonn) for financial support
in terms of a Humboldt Research Award. The first version of the program was written during
his stay at the Philipps University Marburg. We acknowledge the help of Darko Babic (Zagreb), 
Patrick W. Fowler (Sheffield) and David E. Manolopoulos (Oxford) to kindly allow their Fortran routines to be 
modified and implemented in our program system. We also thank Prof. Ottorino Ori (Rome) for fruitful discussions.
<br>

<!-- Main content ends here. -->
</html>
