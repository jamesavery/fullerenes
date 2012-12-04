<?php
/*
The following password needs to be identical to
the one in defined in file    userdata.php
*/

$file='data/fullerenev4.3.zip';
$passwd = $_GET["pwd"];

if ($passwd == "R4TY19KVBN")
  {
        /*
        File download code here
        */
        //download.php
        //content type
        header('Content-type: application/zip');
        //open/save dialog box
        header('Content-Disposition: attachment; filename="fullerenev4.3.zip"');
        //read from server and write to buffer
        readfile($file);
  }
else
  {
   echo "<html>";
   echo "<body>";
   echo "Password incorrect. Please contact us directly if you experience problems.";
   echo "</body>";
   echo "</html>";
  }
?>

