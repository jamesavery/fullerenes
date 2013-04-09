<html>
<body>

<!
--------------------------------------------------------------------------
fulleruf.php processes user information (name, address and email),
  appends it to a file called fulleruser.txt and sends an email out
  with a password for downloading the fullerene program. The email is
  checked for email injections by validating the input
--------------------------------------------------------------------------
>

<!--- get user data-->
<h2>Your user details have been stored in our datafile:</h2>
Title: <?php echo $_GET["Title"]; ?><br>
First Name:  <?php echo $_GET["firstname"]; ?><br> 
Last Name:   <?php echo $_GET["lastname"]; ?><br> 
Institution: <?php echo $_GET["institution"]; ?><br> 
Department:  <?php echo $_GET["department"]; ?><br> 
City:        <?php echo $_GET["city"]; ?><br> 
Country:     <?php echo $_GET["country"]; ?><br> 
Comment:     <?php echo $_GET["comment"]; ?><br>
Email:       <?php echo $_GET["Emailuser"]; ?><br>
<br><br>
<!--- Check email-->
<?php
/*
Setting email variables and password
The following password needs to be identical to
the one in defined in file    fullerdl.php
*/
$passwd="R4TY19KVBN";
$subject = "Program Fullerene";
$from = "p.a.schwerdtfeger@massey.ac.nz";
$headers = "From:" . $from;
$email = $_GET["Emailuser"];
$commentS = $_GET["comment"];
$userfile="fulleruser.txt";
$separator=" & ";
$string="@";
$string1="http";

$check = "True";
/* 
--- Check for if email contains only one @ and no more to avoid compound emails --
*/
if ($_GET["firstname"]=="") $check="False";
if ($_GET["lastname"]=="") $check="False";
if ($_GET["city"]=="") $check="False";
if ($_GET["country"]=="") $check="False";

if ($email=="") $check="False";
$pos = strpos($email,$string);
if ($pos =="") $check="False";
if ($pos =="0") $check="False";
$pos1 = strrpos($email,$string);
if ($pos !="$pos1") $check="False";

/* 
--- Check for comment contains http to avoid spam --
*/
$pos2 = strpos($commentS,$string1);
if ($pos2 !="") $check="False";

/*
--- If data is valid --
*/
if ($check=="True")
  {
mail($email,$subject,$passwd,$headers);
echo "Thank you for using our mail form. Your password will be sent to you shortly."; 
echo " If you do not receive the password please email me directly.";
/* 
--- Store data in fulleruser.txt --
*/
$fh=fopen($userfile,"a+");
fwrite($fh,$_GET["Title"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["firstname"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["lastname"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["institution"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["department"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["city"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["country"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["Emailuser"]);
fwrite($fh,$separator);
fwrite($fh,$_GET["comment"]);
fwrite($fh,"\n");
fclose($fh);
  }
else
/*
--- If data is invalid --
*/
  {
echo "Your entry is not acceptable or incorrect. Please provide all necessary information."; 
echo " If you experience problems please email us directly.";
  }
?>

</body>
</html> 
