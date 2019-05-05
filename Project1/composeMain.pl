#! /usr/bin/perl

use strict;
#use Getopt::Long;

my $source = 'Source2.cpp';
my $out = 'main.cpp';
my %used;
my $res="";
my ($what, $problemName, $problemNumber) = @ARGV; 
# $what = build, test, submit, submit-f, gen
sub gpp{
    #system("g++ -std=c++14 main.cpp -fsanitize=address");
    GetOptions('')
    if(system("g++ -std=c++14 main.cpp")) {die("compile error\n");}
}
sub expand{
    $_ = `g++ -std=c++14 -MM $_[0]`;
    my @files = split(" ", $_);
    my @includes = @files[2..$#files];
    print  "analysing: $_";
    for my $file (@includes){
        if($used{$file}==0){
            print "expanding: $file\n";
            expand($file);
        }
    }   
    print "writing file[1] $files[1]\n";
    $used{$files[1]}=1;
    $res .= `cat $files[1]`;
    $res .= "\n";
}
expand($source);
# removing ,  #include "file name ";
$res =~ s/.*#include\s+".*"//g ;
# removing pragma once
$res =~ s/.*#pragma once.*//g;
open(OUT, "> $out") or die;
print OUT $res; 

my $testCaseDir = "../TestCases";
my $workspace = "$testCaseDir/$problemName/$problemNumber";
# check online judge tool
if(@ARGV == 0){exit(0);}
elsif($what eq 'build'){gpp; print "\nbuild success!\n"; exit(0);}
elsif($problemName =~ /http.*/){
    system("rm -f ./test/*");
    print "oj dl $problemName";
    if(system("oj dl $problemName")){die("downloading test case failed\n");}
    gpp;
    if(system("oj test")) {die("test failed\n");}
    if($what eq 'submit'){
        if(system("oj login $problemName") || 
            system("oj submit $problemName main.cpp")){
            print "submittion failed\n";
        }
    }
    exit(0);
}
else{
    #create testcase if not exist
    unless((-e "$testCaseDir/$problemName") ){
        print "generating testcase\n";
        system("atcoder-tools gen $problemName --workspace $testCaseDir");
        if($what eq 'gen'){exit(0);}
    }
    `cp -f $out $workspace`;
    chdir $workspace;
    gpp;
    print "starting test\n\n";
    $what=~s/^submit$/submit -u/;
    $what=~s/^submit-f$/submit -u -f/;
    system("atcoder-tools $what");
    print "\nWork space is here: $workspace\n";
}



