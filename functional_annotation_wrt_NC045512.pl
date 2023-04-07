open (F, "RefSeq_blastn_3") or die "File not found";
open (E, ">Functional_annotation_output.txt");
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
print E "Total number of mutations in subject sequence $db{$k}: ";
foreach $nuc(@seq)
{
if ($nuc ne '.')
{ $mut++; }
}
print E "$mut\n";
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

if ($ref_pos>=1 and $ref_pos<=265)
{
if ($head == 0)
{
print E "Mutations in 5'UTR region:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 1;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=266 and $ref_pos<=805)
{
if ($head <= 1)
{
print E "Mutations in leader protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 2;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=806 and $ref_pos<=2719)
{
if ($head <= 2)
{
print E "Mutations in nsp2:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 3;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=2720 and $ref_pos<=8554)
{
if ($head <= 3)
{
print E "Mutations in nsp3:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 4;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=8555 and $ref_pos<=10054)
{
if ($head <= 4)
{
print E "Mutations in nsp4:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 5;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=10055 and $ref_pos<=10972)
{
if ($head <= 5)
{
print E "Mutations in 3C-like proteinase:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 6;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=10973 and $ref_pos<=11842)
{
if ($head <= 6)
{
print E "Mutations in nsp6:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 7;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=11843 and $ref_pos<=12091)
{
if ($head <= 7)
{
print E "Mutations in nsp7:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 8;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=12092 and $ref_pos<=12685)
{
if ($head <= 8)
{
print E "Mutations in nsp8:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 9;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=12686 and $ref_pos<=13024)
{
if ($head <= 9)
{
print E "Mutations in nsp9:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 10;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=13025 and $ref_pos<=13441)
{
if ($head <= 10)
{
print E "Mutations in nsp10:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 11;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=13442 and $ref_pos<=16236)
{
if ($head <= 11)
{
print E "Mutations in RNA-dependent RNA-polymerase:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 12;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=16237 and $ref_pos<=18039)
{
if ($head <= 12)
{
print E "Mutations in helicase:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 13;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=18040 and $ref_pos<=19620)
{
if ($head <= 13)
{
print E "Mutations in 3’ to 5’ exonuclease:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 14;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=19621 and $ref_pos<=20658)
{
if ($head <= 14)
{
print E "Mutations in EndoRNAse:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 15;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=20659 and $ref_pos<=21552)
{
if ($head <= 15)
{
print E "Mutations in 2’-O-ribose methyltransferase:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 16;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=21563 and $ref_pos<=25384)
{
if ($head <= 16)
{
print E "Mutations in surface glycoprotein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 17;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=25393 and $ref_pos<=26220)
{
if ($head <= 17)
{
print E "Mutations in ORF3a protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 18;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=26245 and $ref_pos<=26472)
{
if ($head <= 18)
{
print E "Mutations in envelope protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 19;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=26523 and $ref_pos<=27191)
{
if ($head <= 19)
{
print E "Mutations in membrane glycoprotein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 20;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=27202 and $ref_pos<=27387)
{
if ($head <= 20)
{
print E "Mutations in ORF6 protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 21;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=27394 and $ref_pos<=27759)
{
if ($head <= 21)
{
print E "Mutations in ORF7a protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 22;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=27756 and $ref_pos<=27887)
{
if ($head <= 22)
{
print E "Mutations in ORF7b protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 23;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=27894 and $ref_pos<=28259)
{
if ($head <= 23)
{
print E "Mutations in ORF8 protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 24;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=28274 and $ref_pos<=29533)
{
if ($head <= 24)
{
print E "Mutations in nucleocapsid phosphoprotein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 25;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=29558 and $ref_pos<=29674)
{
if ($head <= 25)
{
print E "Mutations in ORF10 protein:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 26;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

elsif ($ref_pos>=29675 and $ref_pos<=29903)
{
if ($head <= 26)
{
print E "Mutations in 3'UTR region:\n";
print E "Ref\tRefPos\tSeq\tSeqPos\n";
$head= 27;
}
if ($ref[$refp-1] ne '-' and $nuc ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t$seq_pos\n"; }
elsif ($ref[$refp-1] ne '-')
{ print E "$ref[$refp-1]\t$ref_pos\t$nuc\t\n"; }
elsif ($nuc ne '-')
{ print E "$ref[$refp-1]\t\t$nuc\t$seq_pos\n"; }
}

}
}
$mut= 0;
$seqp= 0;
$refp= 0;
$seq_gap= 0;
$ref_gap= 0;
$head= 0;
$x= '-'x50;
print E "$x\n";
}
close (F);