#! /usr/bin/perl

use strict;
use Term::ANSIColor;
#use Getopt::Long;

my $source = 'Source2.cpp';
my $out = 'main.cpp';
my %used;
my $res="";
my ($what, $problemName, $problemNumber) = @ARGV; 
# $what = build, test, submit, submit-f, gen
sub gpp{
    #system("g++ -std=c++17 main.cpp -fsanitize=address");
    if(system("g++ -std=c++17 -O2 main.cpp")) {die("compile error\n");}
}
sub expand{
    $_ = `g++ -std=c++14 -MM $_[0]`;
    $_ =~ s/\\//;
    my @files = split(" ", $_);
    my @includes = @files[2..$#files];
    print  "analysing: $_";
    my $num =0;
    for my $file (@includes){
        if($used{$file}==0){
            print "expanding: $file\n";
            expand($file);
        }
        $num++;
    }   
    print "writing file[1] $files[1]\n";
    $used{$files[1]}=1;
    $res .= `cat $files[1]`;
    $res .= "\n";
}
expand($source);
# removing ,  #include "file name " except for boost
$res =~ s/.*#include\s+"(?!boost).*"//g ;
# removing pragma once
$res =~ s/.*#pragma once.*//g;
open(OUT, "> $out") or die;
print OUT $res; 

my $testCaseDir = "../TestCases";
my $workspace = "$testCaseDir/$problemName/$problemNumber";
# check online judge tool
if(@ARGV == 0){exit(0);}
elsif($what eq 'build'){gpp; print "\nbuild success!\n"; exit(0);}
elsif($what eq 'open'){gpp; print "\nbuild success!\n"; `notepad.exe ./$out`; exit(0);}
elsif($what eq 'code'){`code $out`; exit(0);}
elsif($what eq 'test-man'){
    gpp;
    print colored("build success!\n",'cyan');
    print colored("Input EOF\n",'cyan');
    my $EOF = <STDIN>; chomp($EOF);
    while(1){
        print colored("Input test in data: \n",'cyan');
        my @test_in;
        while(<STDIN>){
            chomp($_);
            if($_ ne $EOF){
                push(@test_in, $_);
            }
            else{ last; }
        }
        # test begin
        my $test_in = join('\n', @test_in);
        print colored("result: \n", "yellow");
        system("echo $test_in | ./a.out");
    }
}
elsif($problemName =~ /http.*/){
    system("rm -f ./test/*");
    print "oj dl $problemName";
    if(system("oj dl $problemName")){die("downloading test case failed\n");}
    gpp;
    if(system("oj test")) {die("test failed\n");}
    if($what eq 'submit'){
	$ENV{"BROWSER"} = "/mnt/c/Program Files (x86)/Google/Chrome/Application/chrome.exe";
        if(system("oj login $problemName") || 
            system(" oj submit -y $problemName main.cpp")){
            print "submittion failed\n";
        }
    }
    elsif($what eq 'test-sys'){
        if(system("oj dl $problemName --system")){die("downloading test case failed\n");}
        if(system("oj test")) {die("test failed\n");}
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
    if($what eq "submit-f"){
    	print colored("submit without check. ok? (y/n):", "red");
	my $ok = <STDIN>; chomp($ok);
	if($ok ne "y"){
		exit 1;
	}
    }
    $what=~s/^submit$/submit -u/;
    $what=~s/^submit-f$/submit -u -f/;
    system("atcoder-tools $what");
    print "\nWork space is here: $workspace\n";
}



