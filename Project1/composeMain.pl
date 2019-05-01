use strict;

my $source = 'Source2.cpp';
my $out = 'main.cpp';
my %used;
my $res="";
sub expand{
    $_ = `g++ -MM $_[0]`;
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
$res =~ s/#include\s+".*"//g ;
# removing pragma once
$res =~ s/.*#pragma once.*//g;
open(OUT, "> $out") or die;
print OUT $res; 

# select problem 
my ($what, $problemName, $problemNumber) = $ARGV; #test or submit or gen
my $testCaseDir = "../TestCases";
my $workspace = "$testCaseDir/$problemName/$problemNumber";
#create testcase if not exist
unless(-e "$testCaseDir/$problemName" ){
    print "generating testcase\n";
    system("atcoder-tools gen $problemName --workspace $testCaseDir");
    if($what eq 'gen'){exit(0);}
}
`cp -f $out $workspace`;
chdir $workspace;
if($what eq 'test'){
    system("g++ -std=c++14 main.cpp -fsanitize=address");
}
print "starting test\n\n";
system("atcoder-tools $what");
print "\nWork space is here: $workspace\n";


