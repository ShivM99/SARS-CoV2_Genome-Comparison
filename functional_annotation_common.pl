open (F, "RefSeq_blastn_3") or die "File not found";
%db= (0=>"MW422256.1", 1=>"MW580576.1", 2=>"MW520923.1", 3=>"MW931310.1", 4=>"OM059472.1");
%seq_space= ();
%seq= ();
foreach $k(keys%db)
{
$seq_space{$k}= 0;
$seq{$k}= '';
$count++;
}
while (chomp($line=<F>))
{
if ($line =~ /^Query_[0-9]+ +([0-9]+) +[ATtGCRYSWKMBDHVN\-]+/)
{ $ref_start= $1-1; }
if ($line =~ /^(Query_[0-9]+ +[0-9]+ +)[ATtGCRYSWKMBDHVN\-]+/)
{ $start= length($1); }
elsif ($line =~ /^([0-9]+) +[0-9]+ +[ATtGCRYSWKMBDHVN\-\.]+ +/)
{
$line= reverse($line);
for ($i=1; $i<=$start; $i++)
{ chop($line); }
$line= reverse($line);
$space= ($line =~ tr/ / /);
$seq_space{$1} = $space-2;
$c++;
if ($c == $count)
{ last; }
}
}
seek (F,0,0);
while (chomp($line=<F>))
{
if ($line =~ /^Query_[0-9]+ +[0-9]+ +([ATtGCRYSWKMBDHVN\-]+)/)
{ $ref .= $1; }
elsif ($line =~ /^([0-9]+) +[0-9]+ +([ATtGCRYSWKMBDHVN\-\.]+) +/)
{ $seq{$1} .= $2; }
elsif ($line =~ /^Effective search space used:/)
{ last; }
}
@ref= split//,$ref;
foreach $k(keys%db)
{
@seq= split//,$seq{$k};
print "Total number of mutations in subject sequence $db{$k}: ";
foreach $nuc(@seq)
{
if ($nuc ne '.')
{ $mut++; }
}
print "$mut\n";
print "Ref\tRefPos\tSeq\tSeqPos\n";
$refp= $seq_space{$k};
foreach $nuc(@seq)
{
$seqp++;
$refp++;
if ($nuc eq '-')
{ $seq_gap++; }
if ($ref[$refp-1] eq '-')
{ $ref_gap++; }
if ($nuc ne '.')
{
$seq_pos= $seqp - $seq_gap;
$ref_pos= $ref_start + $refp - $ref_gap;
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}
}
$mut= 0;
$seqp= 0;
$refp= 0;
$seq_gap= 0;
$ref_gap= 0;
$x= '-'x50;
print "$x\n";
}
close (F);