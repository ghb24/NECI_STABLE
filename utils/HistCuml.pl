#!/usr/bin/env perl
open(F1,"<$ARGV[0]")|| die("Cannot open file '$ARGV[0]'");
open(F2,"<$ARGV[1]")|| die("Cannot open file '$ARGV[1]'");
while($v1=<F1>)
{
   $v2=<F2>;
   ($num1)=($v1=~/\s*(\S*)/);
   ($num2)=($v2=~/\s*(\S*)/);
   if($num1!=0.0) {   print "",($num1-$num2)/$num1,"\n";}
   else {print "0\n";}
   
}
close F1;
close F2;


