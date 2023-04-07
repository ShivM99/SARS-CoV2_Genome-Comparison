open (F, "All_variants_blastn_3") or die "File not found";
$id= ();
%identity= ();
%len= ();
while (chomp($line=<F>))
{
if ($line =~ /^Query= ([A-Z0-9_\.]+) /)
{
$c++;
$id{$c}= $1;
$identity{$c}= 0;
}
elsif ($line =~ /^Length=([0-9]+)/)
{ $len{$id{$c}}= $1; }
}
print "Number of sequences: $c\n";
foreach $s(keys%id)
{ print "Sequence ID: $id{$s}\tLength: $len{$id{$s}}\n"; }
seek (F,0,0);
while (chomp($line=<F>))
{
if ($line =~ /^Query= ([A-Z0-9_\.]+) /)
{
$x= '-'x50;
print "$x\n";
print "Query sequence ID: $1\t";
}
elsif ($line =~ /^Length=([0-9]+)/)
{ print "Length of the sequence: $1\n"; }
elsif ($line =~ /^Sequences producing significant alignments:/)
{
$line= <F>;
while (chomp($line=<F>))
{
if ($line =~ /^$/)
{ last; }
elsif ($line =~ /^([A-Z0-9_\.]+) .*  ([0-9]+) +([0-9\.]+)/)
{ print "Subject sequence ID: $1\tBit score: $2\tE-value: $3\n"; }
}
}
elsif ($line =~ /^Query_[0-9]+ +[0-9]+ +[ATGCN\-]+/)
{
while (chomp($line=<F>))
{
if ($line =~ /^Effective search space used:/)
{ last; }
elsif ($line =~ /^([0-9]+) +[0-9]+ +([ATGCN\.\-]+) +/)
{ $identity{$1+1} += ($2 =~ tr/\./\./); }
}
foreach $k(keys%identity)
{
$per= ($identity{$k}*100)/$len{$id{$k}};
print "Subject sequence: $id{$k}\t%identity: $per\n";
$identity{$k}= 0;
}
}
}
close (F);