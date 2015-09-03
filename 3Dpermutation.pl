use strict;
use Getopt::Long;
use lib("/System/Library/Perl/5.10.0/darwin-thread-multi-2level/");

my %CodonChart=('TTT','F', 'TCT','S', 'TAT','Y', 'TGT','C',
                'TTC','F', 'TCC','S', 'TAC','Y', 'TGC','C',
                'TTA','L', 'TCA','S', 'TAA','X', 'TGA','X',
                'TTG','L', 'TCG','S', 'TAG','X', 'TGG','W',
                'CTT','L', 'CCT','P', 'CAT','H', 'CGT','R',
                'CTC','L', 'CCC','P', 'CAC','H', 'CGC','R',
                'CTA','L', 'CCA','P', 'CAA','Q', 'CGA','R',
                'CTG','L', 'CCG','P', 'CAG','Q', 'CGG','R',
                'ATT','I', 'ACT','T', 'AAT','N', 'AGT','S',
                'ATC','I', 'ACC','T', 'AAC','N', 'AGC','S',
                'ATA','I', 'ACA','T', 'AAA','K', 'AGA','R',
                'ATG','M', 'ACG','T', 'AAG','K', 'AGG','R',
                'GTT','V', 'GCT','A', 'GAT','D', 'GGT','G',
                'GTC','V', 'GCC','A', 'GAC','D', 'GGC','G',
                'GTA','V', 'GCA','A', 'GAA','E', 'GGA','G',
                'GTG','V', 'GCG','A', 'GAG','E', 'GGG','G');

my %AA = (ALA => "A", CYS => "C", ASP => "D", GLU => "E", PHE => "F", GLY => "G", HIS => "H",
ILE => "I", LYS => "K", LEU => "L", MET => "M", ASN => "N", PRO => "P", GLN => "Q",
ARG => "R", SER => "S", THR => "T", VAL => "V", TRP => "W", TYR => "Y");

my $mRNA_id = undef;
my $PDB_file = undef;
my $mutation_pos = undef;
my $num_mut_cutoff = 2;
my $gap_ratio_cutoff = 0.5;
my $length_ratio_cutoff = 0.5;
my $chan;
my $num_permutation;
my $reference_seqence_file;
my $gene_annotation_file;

GetOptions(
	"G=s" => \$mRNA_id,
	"I=s" => \$PDB_file,
	"C=s" => \$chan,
	"M=s" => \$mutation_pos,
	"N=i" => \$num_mut_cutoff,
	"A=i" => \$gap_ratio_cutoff,
	"L=i" => \$length_ratio_cutoff,
	"P=i" => \$num_permutation,
	"R=s" => \$reference_seqence_file,
	"A=s" => \$gene_annotation_file,
);

if(! $mRNA_id){die"-G <mRNA id (NM_001128226)> !!\n";}
if(! -f $PDB_file){die"-I <file name> !!\n";}
if(! $mutation_pos){die"-M <mutation position (100:200)> !!\n";}
if(! $chan){die"-C <chan of the target (A)> !!\n";}
if($gap_ratio_cutoff > 1 || $gap_ratio_cutoff < 0){die "-U <gap ratio 0-1> !!";}
if($length_ratio_cutoff > 1 || $length_ratio_cutoff < 0){die "-L <ratio of sequence length 0-1> !!";}
if(! $num_permutation){$num_permutation = 10000;}
if(! -f $reference_seqence_file){die "-R <Reference sequence file> !!";}
if(! -f $gene_annotation_file){die "-A <Gene annotation file> !!";}

#my $gene_data = "/home/icgc/WORK/fujimoto/src/3D/seq_gene_like.refGene_110529.txt.d";
#my $seq_gene_md = "/Users/fujimoto/Desktop/solexa/3D_structure/Github/seq_gene_like.gencode.v19.txt";
#my $samtools_seq = "/Users/fujimoto/Desktop/solexa/3D_structure/Github/All.fa";

my %gene_seq = get_gene_seq2($mRNA_id, $gene_annotation_file, $reference_seqence_file, \%CodonChart);

my $gene_symbol = $gene_seq{gene_symbol};

my %moelcule_ave_pos = &get_structure($chan, $PDB_file);

if(keys(%moelcule_ave_pos) == 0){warn "chain $chan data not found in $PDB_file\n"; exit;}

my @aligned_sequence = &align_by_mafft(\%gene_seq, \%moelcule_ave_pos);

my ($length_ratio, $gap_ratio) = &check_sequence_identity(\@aligned_sequence);

if($length_ratio < $length_ratio_cutoff || $gap_ratio < $gap_ratio_cutoff){
	my @out_low_name = ("name", "gene_symbol", "chan", "mRNA_id", "PDB_file", "length_ratio_cutoff", "gap_ratio_cutoff", "length_ratio", "gap_ratio", "num_converted_mutation", "mutation_pos", "mutation_pos_on_PDB", "mutation_num", "ave_mutation_dis", "permutaiton_P_value", "median_mutation_dis", "median_permutaiton_P_value", "sd_mutation_dis", "sd_permutaiton_P_value", "Chimera_command"),"\n";
	my @out_low = ("result", $gene_symbol, $chan, $mRNA_id, $PDB_file, $length_ratio_cutoff, $gap_ratio_cutoff, $length_ratio, $gap_ratio, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA");
	
	print join("\t", @out_low_name), "\n";
	print join("\t", @out_low), "\n";
	exit;
}

my %mutation_pos_on_PDB = &calculate_mutation_pos3(\@aligned_sequence, $mutation_pos);
my %mutation_pos_on_PDB2 = &calculate_mutation_pos(\@aligned_sequence, $mutation_pos);

my @mutation_pos_on_PDB;
my @mutation_pos2;
foreach(sort {$a <=> $b} keys %mutation_pos_on_PDB){
	push(@mutation_pos2, $_);
	push(@mutation_pos_on_PDB, $mutation_pos_on_PDB{$_});
}

my $mutation_pos_on_PDB = join(":", @mutation_pos_on_PDB);

my $mutation_pos_on_PDB_for_permutation;
my @mutation_pos_on_PDB_for_permutation;
my $num_converted_mutation = 0;

foreach my $pos_on_UCSC ((split(":", $mutation_pos))){
	if(exists($mutation_pos_on_PDB{$pos_on_UCSC}) && $mutation_pos_on_PDB2{$pos_on_UCSC} ne "NA"){
		push(@mutation_pos_on_PDB_for_permutation, $mutation_pos_on_PDB2{$pos_on_UCSC});
		$num_converted_mutation++;
	}
}

$mutation_pos_on_PDB_for_permutation = join(":", @mutation_pos_on_PDB_for_permutation);

if($num_converted_mutation < $num_mut_cutoff){
	my @out_row_name = ("gene_symbol", "chan", "mRNA_id", "PDB_file", "Length_ratio", "Gap_ratio", "Number_of_converted_mutation", "Mutation_pos", "Mutation_pos_on_PDB", "Average_mutation_distance", "Permutaiton_P_value"),"\n";
	my @out_row = ($gene_symbol, $chan, $mRNA_id, $PDB_file, $length_ratio, $gap_ratio, $num_converted_mutation, "NA", "NA", "NA", "NA");

	warn join("\t", @out_row_name), "\n";
	warn join("\t", @out_row), "\n";
	exit;
}

my ($mutation_num, $ave_mutation_dis, $median_mutation_dis, $sd_mutation_dis, $permutaiton_P_value_ave, $permutaiton_P_value_median, $permutaiton_P_value_sd) = permutation_3D_pos(\%moelcule_ave_pos, $mutation_pos_on_PDB_for_permutation, $num_permutation);
my @out_row_name = ("gene_symbol", "chain", "mRNA_id", "PDB_file", "Length_ratio", "Gap_ratio", "Number_of_mapped_mutation", "Mutation_position", "Mutation_position_on_protein", "Average_distance_between_mutations", "Permutaiton_P_value"),"\n";
my @out_row = ($gene_symbol, $chan, $mRNA_id, $PDB_file, $length_ratio, $gap_ratio, $num_converted_mutation, join(",", @mutation_pos2), $mutation_pos_on_PDB, $ave_mutation_dis, $permutaiton_P_value_ave);

print join("\t", @out_row_name), "\n";
print join("\t", @out_row), "\n";

sub check_sequence_identity{
	my $aligned_sequence = shift;
	my @aligned_sequence = @{$aligned_sequence};

	my $total_PDB = 0;
	my $total_UCSC = 0;
	my $gap_PDB = 0;
	my $gap_UCSC = 0;
	foreach(@aligned_sequence){
		my ($consensus, $pos_seq_UCSC, $seq_UCSC, $pos_seq_PDB, $pos_seq_PDB2, $seq_PDB) = @{$_};
		if($seq_UCSC =~ /[A-Z]/){$total_UCSC++;}
		if($seq_UCSC =~ /[A-Z]/ && $seq_PDB eq "-"){$gap_UCSC++;}
		if($seq_PDB =~ /[A-Z]/){$total_PDB++;}
		if($seq_PDB =~ /[A-Z]/ && $seq_UCSC eq "-"){$gap_PDB++;}
	}

	my $length_ratio;
	my $gap_ratio;
	if($total_UCSC > $total_PDB){
		$length_ratio = $total_UCSC/$total_PDB;
		$gap_ratio = ($total_PDB - $gap_PDB)/$total_PDB;
	}
	else{
		if($total_UCSC == 0){
			$length_ratio = 0;
			$gap_ratio = 1;
		}
		else{
			$length_ratio = $total_PDB/$total_UCSC;
			$gap_ratio = ($total_UCSC - $gap_UCSC)/$total_UCSC;
		}
	}
	
	return(&round($length_ratio, 3), &round($gap_ratio, 3));
}

sub get_window_seq_in_UCSC{
	my ($window_pos_UCSC, $gene_seq) = @_;
	my @window_pos_UCSC = split(":", $window_pos_UCSC);
	my %gene_seq = %{$gene_seq};
	
	my %window_seq_length_result;
	foreach (@{$gene_seq{seq_array}}){
		my @l = split("\t", $_);
		my $found = 0;
		foreach(@window_pos_UCSC){
			if($l[3] == $_){$found++;}
		}
		if($found > 0){
			$window_seq_length_result{length}++;
			if($l[6] ne $l[7]){$window_seq_length_result{nonsyn}++;}
		}
	}
	
	$window_seq_length_result{length} = $window_seq_length_result{length}/3;
	$window_seq_length_result{nonsyn} = $window_seq_length_result{nonsyn}/3;
	
	return %window_seq_length_result;
}

sub calculate_mutation_pos2{
	my ($aligned_sequence, $mutation_pos) = @_;
	my @aligned_sequence = @{$aligned_sequence};

	my %mutation_pos;
	foreach((split(":", $mutation_pos))){$mutation_pos{$_} = 1;}

	my %mutaiton_pos_on_UCSC;
	foreach(@aligned_sequence){
		my ($consensus, $pos_seq_UCSC, $seq_UCSC, $pos_seq_PDB, $pos_seq_PDB2, $seq_PDB) = @{$_};

		if($mutation_pos{$pos_seq_PDB2} == 1 && $seq_UCSC ne "-"){$mutaiton_pos_on_UCSC{$pos_seq_PDB2} = $pos_seq_UCSC;}
		
		if(keys(%mutaiton_pos_on_UCSC) == keys(%mutation_pos)){last;}
	}
	
	foreach my $pos (keys %mutation_pos){
		if(defined($mutaiton_pos_on_UCSC{$pos}) != 1 || $mutaiton_pos_on_UCSC{$pos} == 0){$mutaiton_pos_on_UCSC{$pos} = "NA";}
	}
	
	return %mutaiton_pos_on_UCSC;
}

sub calculate_mutation_pos{
	my ($aligned_sequence, $mutation_pos) = @_;
	my @aligned_sequence = @{$aligned_sequence};

	my %mutation_pos;
	foreach((split(":", $mutation_pos))){$mutation_pos{$_} = 1;}

	my %mutaiton_pos_on_PDB;
	foreach(@aligned_sequence){
		my ($consensus, $pos_seq_UCSC, $seq_UCSC, $pos_seq_PDB, $pos_seq_PDB2, $seq_PDB) = @{$_};

		if($mutation_pos{$pos_seq_UCSC} == 1 && $seq_PDB ne "-"){
			$mutaiton_pos_on_PDB{$pos_seq_UCSC} = $pos_seq_PDB2;
		}
		
		if(keys(%mutaiton_pos_on_PDB) == keys(%mutation_pos)){last;}
	}
	
	foreach my $pos (keys %mutation_pos){
		if(defined($mutaiton_pos_on_PDB{$pos}) != 1 || $mutaiton_pos_on_PDB{$pos} == 0){$mutaiton_pos_on_PDB{$pos} = "NA";}
	}
	
	return %mutaiton_pos_on_PDB;
}

sub calculate_mutation_pos3{
	my ($aligned_sequence, $mutation_pos) = @_;
	my @aligned_sequence = @{$aligned_sequence};

	my %mutation_pos;
	foreach((split(":", $mutation_pos))){$mutation_pos{$_} = 1;}

	my %mutaiton_pos_on_PDB;
	foreach(@aligned_sequence){
		my ($consensus, $pos_seq_UCSC, $seq_UCSC, $pos_seq_PDB, $pos_seq_PDB2, $seq_PDB) = @{$_};

		if($mutation_pos{$pos_seq_UCSC} == 1 && $seq_PDB ne "-"){
			$mutaiton_pos_on_PDB{$pos_seq_UCSC} = $pos_seq_PDB; #140321
		}
		
		if(keys(%mutaiton_pos_on_PDB) == keys(%mutation_pos)){last;}
	}
	
	foreach my $pos (keys %mutation_pos){
		if(defined($mutaiton_pos_on_PDB{$pos}) != 1 || $mutaiton_pos_on_PDB{$pos} == 0){$mutaiton_pos_on_PDB{$pos} = "NA";}
	}
	
	return %mutaiton_pos_on_PDB;
}

sub align_by_mafft{
	my ($gene_seq, $moelcule_ave_pos) = @_;
	my %gene_seq = %{$gene_seq};
	my %moelcule_ave_pos = %{$moelcule_ave_pos};
	
	my $seq_in_PDB;
	foreach my $pos (sort {$a <=> $b} keys %moelcule_ave_pos){
		if(defined($AA{${$moelcule_ave_pos{$pos}}[0]}) == 1){$seq_in_PDB = "$seq_in_PDB"."$AA{${$moelcule_ave_pos{$pos}}[0]}";}
		else{$seq_in_PDB = "$seq_in_PDB"."X";}
	}
	
	my %seq_UCSC;
	foreach(@{$gene_seq{seq_array}}){
		my @l = split("\t");
		$seq_UCSC{$l[0]} = $l[1];
	}
	
	my @seq_in_UCSC;
	foreach(sort {$a <=> $b} keys %seq_UCSC){push(@seq_in_UCSC, $seq_UCSC{$_});}
	my $seq_in_UCSC = join("", @seq_in_UCSC);
	
	my $input_file_for_mafft = "$PDB_file".".fas";
	open MAFFT_INPUT, ">$input_file_for_mafft" or die "Can not open $input_file_for_mafft !!";
	print MAFFT_INPUT ">$gene_seq{gene_symbol}UCSC\n$seq_in_UCSC\n>$gene_seq{gene_symbol}PDB\n$seq_in_PDB\n";
	
	my $output_file_for_mafft = "$PDB_file".".fas.align";
	my $mafft_cmd = "(mafft $input_file_for_mafft > $output_file_for_mafft) >& /dev/null";
	system("$mafft_cmd");
	
	open MAFFT_OUTPUT, "$output_file_for_mafft" or die "Can not open $output_file_for_mafft !!";
	my %seq;
	my $seq_name;
	while(<MAFFT_OUTPUT>){
		chomp;
		if($_ =~ />/ && $_ =~ /UCSC$/){
			$seq_name = "UCSC";
		}
		if($_ =~ />/ && $_ =~ /PDB$/){
			$seq_name = "PDB";
		}
		
		$seq{$seq_name} .= $_;
	}
	close(MAFFT_OUTPUT);
	
	my @seq_UCSC = split("", $seq{UCSC});
	my @seq_PDB = split("", $seq{PDB});
	
	my @pos = sort {$a <=> $b} keys %moelcule_ave_pos;
	
	my @aligned_sequence = ();
	my $consensus = 0;
	my $pos_seq_UCSC = 0;
	my $pos_seq_PDB = 0;
	my $pos_seq_PDB2 = 0;
	for(my $i = 0; $i < @seq_UCSC; $i++){
		$consensus++;
		if($seq_UCSC[$i] ne "-"){$pos_seq_UCSC++;}
		if($seq_PDB[$i] ne "-"){
			$pos_seq_PDB2 = $pos[$pos_seq_PDB];
			$pos_seq_PDB++;
		}
		
		push(@aligned_sequence, [$consensus, $pos_seq_UCSC, $seq_UCSC[$i], $pos_seq_PDB, $pos_seq_PDB2, $seq_PDB[$i]]);
	}

	return @aligned_sequence;
}

sub get_DNAseq{
        my $chr = shift;
        my $start = shift;
        my $last = shift;
        my $seqfile = shift;

        my $sequence = `samtools faidx $seqfile \"chr$chr:$start-$last\"|grep \"^>\" -v|tr -d '\n'`;

        return $sequence;
}

sub get_gene_seq2{
	my ($mRNA_id, $seq_gene_md, $samtools_seq, $CodonChart) = @_;
	
	my %CodonChart = %{$CodonChart};

	open GENE_MD, "$seq_gene_md" or die "$seq_gene_md !!";

	my $strand;
	my $seq;
	while(<GENE_MD>){
		chomp;
		my @l = split("\t", $_);
		if($l[-1] ne $mRNA_id){next}
		else{
			my ($chr, $start, $last) = ($l[2], $l[5], $l[6]);
			my $seq_tmp = &get_DNAseq($chr, $start, $last, $samtools_seq);
			$seq .= uc($seq_tmp);
			$strand = $l[3];
			$gene_symbol = $l[0];
#			last;
		}
	}

	if($strand eq "-"){
		$seq =~ tr/[ATGC]/[TACG]/;
		$seq = reverse($seq);
	}

	my %gene_seq = ("mRNA_id" => $mRNA_id, "gene_symbol" => $gene_symbol);
	my $aa_seq;
	for(my $i = 1; $i < length($seq)/3; $i++){
		my $codon_start = 3*($i - 1);
		my $aa_tmp = substr($seq, $codon_start, 3);
		$aa_seq .= $aa_tmp;
		my $tmp = "$i"."\t"."$CodonChart{$aa_tmp}";
		push(@{$gene_seq{seq_array}}, $tmp);
	}

	return %gene_seq;
}

sub get_window_summary{
	my $window = shift;
	my %window = %{$window};

	my ($num_window, $ave_mutation, $ave_size) = (0, 0, 0);

	foreach my $pos (sort {$a <=> $b} keys %window){
		$num_window++;
		$ave_size += @{${$window{$pos}}{pos}};
	}
	
	return  ($num_window, &round($ave_size/$num_window, 3));
}

sub get_mutation_window{
	my ($moelcule_ave_pos, $mutation_pos, $diameter) = @_;
	my %moelcule_ave_pos = %{$moelcule_ave_pos};
	my %mutation_pos;
	foreach((split(":", $mutation_pos))){$mutation_pos{$_} = 1;}

	my %window;
	foreach my $pos1 (sort {$a <=> $b} keys %moelcule_ave_pos){
		my ($x1, $y1, $z1) = (${$moelcule_ave_pos{$pos1}}[1], ${$moelcule_ave_pos{$pos1}}[2], ${$moelcule_ave_pos{$pos1}}[3]);
		
		foreach my $pos2 (sort {$a <=> $b} keys %moelcule_ave_pos){
			my ($x2, $y2, $z2) = (${$moelcule_ave_pos{$pos2}}[1], ${$moelcule_ave_pos{$pos2}}[2], ${$moelcule_ave_pos{$pos2}}[3]);
			my $distance = (($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2)**(1/2);
			if($distance <= $diameter){push(@{${$window{$pos1}}{pos}}, $pos2);}
		}
	}

	foreach my $pos (sort {$a <=> $b} keys %window){
		my $num = 0;
		
		foreach(@{${$window{$pos}}{pos}}){
			$num += $mutation_pos{$_}; 
			if($mutation_pos{$_} == 1){push(@{${$window{$pos}}{mut}}, $_);}
		}
		${$window{$pos}}{num} = $num;
	}
	return %window;
}

sub permutation_3D_pos{
	my ($moelcule_ave_pos, $mutation_pos, $num_permutation) = @_;
	my %moelcule_ave_pos = %{$moelcule_ave_pos};
	my @mutation_pos = split(":", $mutation_pos);
	my $num_combination = 0;
	my $total_mutation_dis = 0;
	my @combination;
	my @distance;

	for(my $i = 0; $i < @mutation_pos; $i++){
		my $pos1 = $mutation_pos[$i];
		if($pos1 =~ /[^0-9]/){next;}
		my ($x1, $y1, $z1) = (${$moelcule_ave_pos{$pos1}}[1], ${$moelcule_ave_pos{$pos1}}[2], ${$moelcule_ave_pos{$pos1}}[3]);
		for(my $j = 0; $j < @mutation_pos; $j++){
			if($i <= $j){next;}
			my $pos2 = $mutation_pos[$j];
			if($pos2 =~ /[^0-9]/){next;}
			my ($x2, $y2, $z2) = (${$moelcule_ave_pos{$pos2}}[1], ${$moelcule_ave_pos{$pos2}}[2], ${$moelcule_ave_pos{$pos2}}[3]);
			my $distance = (($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2)**(1/2);
			push(@distance, $distance);
		}
	}

	my $ave_mutation_dis = round(average(@distance), 2);
	my $median_mutation_dis = round(median(@distance), 2);
	my $sd_mutation_dis = round(sd(@distance), 2);

	my $mutation_num = 0;
	foreach(@mutation_pos){if($_ =~ /[^0-9]/){next;} $mutation_num++;}

	my @total_sim_dis_ave;
	my @total_sim_dis_median;
	my @total_sim_dis_sd;
	my @moelcule_ave_pos = keys(%moelcule_ave_pos);
	for(my $i = 1; $i <= $num_permutation; $i++){
		my $total_sim_dis = 0;
		my $num_combination2 = 0;
		my @random_pos = shuffle_array(\@moelcule_ave_pos, $mutation_num);

		my @combination;
		my @random_distance;
		for(my $j = 0; $j < $mutation_num; $j++){
			my $pos1 = $random_pos[$j];
			my ($x1, $y1, $z1) = (${$moelcule_ave_pos{$pos1}}[1], ${$moelcule_ave_pos{$pos1}}[2], ${$moelcule_ave_pos{$pos1}}[3]);
			for(my $k = 0; $k < $j; $k++){
				my $pos2 = $random_pos[$k];
				my ($x2, $y2, $z2) = (${$moelcule_ave_pos{$pos2}}[1], ${$moelcule_ave_pos{$pos2}}[2], ${$moelcule_ave_pos{$pos2}}[3]);
				my $distance = (($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2)**(1/2);
				push(@random_distance, $distance);
			}
		}
		my $ave_sim_dis = round(average(@random_distance), 2);
		my $median_sim_dis = round(median(@random_distance), 2);
		my $sd_sim_dis = round(sd(@random_distance),2);

		push(@total_sim_dis_ave, $ave_sim_dis);
		push(@total_sim_dis_median, $median_sim_dis);
		push(@total_sim_dis_sd, $sd_sim_dis);

		if(($i)%int($num_permutation/10) == 0){warn"$i", "th permutation\n";}
	}
	my $num_permutation = @total_sim_dis_ave;
	my $num_lower_than_obs_ave = 0;
	foreach(@total_sim_dis_ave){if($_ < $ave_mutation_dis){$num_lower_than_obs_ave++;}}
	
	my $num_lower_than_obs_median = 0;
	foreach(@total_sim_dis_median){if($_ < $median_mutation_dis){$num_lower_than_obs_median++;}}
	
	my $num_lower_than_obs_sd = 0;
	foreach(@total_sim_dis_sd){if($_ < $sd_mutation_dis){$num_lower_than_obs_sd++;}}

	my $permutaiton_P_value_ave = $num_lower_than_obs_ave/$num_permutation;
	my $permutaiton_P_value_median = $num_lower_than_obs_median/$num_permutation;
	my $permutaiton_P_value_sd = $num_lower_than_obs_sd/$num_permutation;

	return ($mutation_num, $ave_mutation_dis, $median_mutation_dis, $sd_mutation_dis, $permutaiton_P_value_ave, $permutaiton_P_value_median, $permutaiton_P_value_sd);
}

sub sd{
	my @data = @_;

	my $mean = 0;
	my $n = @data;

	foreach my $i (@data){
		$mean += $i/$n;
	}

	my $square_sum = 0;
	foreach my $i (@data){
		$square_sum += ( $i - $mean )**2;
	}
	
	my $sd=($square_sum/$n)**(1/2);
	
	return($sd);
}

sub average{
	my @data = @_;

	my $sum = 0;
	foreach(@data){
		$sum+=$_;
	}
	return ($sum/@data);
}

sub median{
	my @data = @_;

	@data = sort { $a <=> $b} @data;

	my $n_data = @data;
	my $median; 
	if(($n_data%2) == 1){
		$median = $data[(($n_data + 1)/2) - 1];
	}else{
		$median = (($data[$n_data/2 - 1] + $data[($n_data/2 + 1) - 1]) / 2);
	}
	return($median);
}

sub shuffle_array{
	my ($array, $number) = @_;
	my @array = @{$array};

	my @array2;
	for(my $i = 0; $i < $number; $i++){
		my $j = int rand (@array + 1);
		push(@array2, $array[$j]);
	}

	return @array2;
}

sub select_window{
	my ($window, $mutation_pos, $num_mut_in_window) = @_;
	my %window = %{$window};
	my %mutation_pos;
	foreach((split(":", $mutation_pos))){$mutation_pos{$_} = 1;}

	my %num_mutation;
	foreach my $pos (sort {$a <=> $b} keys %window){
		my $num = ${$window{$pos}}{num};
		$num_mutation{$num}++;
	}

	my @mutation;
	foreach(sort {$b <=> $a} keys %num_mutation){
		if($_ >= $num_mut_in_window){push(@mutation, $_);}
	}

	my %selected_window;
	my %mutation_in_selected_windows;
	foreach(@mutation){
		my $num_mut = $_;
		if($num_mut < $num_mut_in_window){next;}

		if(keys(%mutation_in_selected_windows) > 0){
			foreach my $pos (sort {$a <=> $b} keys %window){
				my @tmp_mut;
				foreach(@{${$window{$pos}}{mut}}){if(defined($mutation_in_selected_windows{$_}) != 1){push(@tmp_mut, $_);}}
				@{${$window{$pos}}{mut}} = @tmp_mut;
				${$window{$pos}}{num} = @tmp_mut;
			}
		}

		my %window_tmp;
		my $num_of_window = 0;
		foreach my $pos (sort {$a <=> $b} keys %window){
			if(${$window{$pos}}{num} != $num_mut){next;}
			my $aa_num = @{${$window{$pos}}{pos}};
			$window_tmp{$pos} = $window{$pos};
			$num_of_window++;
		}

		while($num_of_window > 0){
			my $min_len;
			my $selected_pos;
			my $j = 0;
			foreach my $pos (keys %window_tmp){
				$j++;
				if($j == 1){$min_len = @{${$window_tmp{$pos}}{pos}}; $selected_pos = $pos;}

				if($min_len > @{${$window_tmp{$pos}}{pos}}){
					$min_len = @{${$window_tmp{$pos}}{pos}};
					$selected_pos = $pos;
				}
			}
			${$selected_window{$selected_pos}}{num} = ${$window_tmp{$selected_pos}}{num};
			@{${$selected_window{$selected_pos}}{mut}} = @{${$window_tmp{$selected_pos}}{mut}};
			@{${$selected_window{$selected_pos}}{pos}} = @{${$window_tmp{$selected_pos}}{pos}};

			foreach(@{${$selected_window{$selected_pos}}{mut}}){$mutation_in_selected_windows{$_} = 1;}

			foreach my $pos (keys %window_tmp){
				my @tmp_mut;
				foreach(@{${$window_tmp{$pos}}{mut}}){
					my $mutation_pos = $_;
					if($mutation_in_selected_windows{$mutation_pos} != 1){push(@tmp_mut, $mutation_pos);}
				}
				${$window_tmp{$pos}}{num} = @tmp_mut;
				@{${$window_tmp{$pos}}{mut}} = @tmp_mut;
				if(@tmp_mut != $num_mut){delete($window_tmp{$pos});}
			}

			my $aa_num = @{${$selected_window{$selected_pos}}{pos}};

			$num_of_window = keys(%window_tmp);
		}
	}

	return %selected_window;
}

sub get_structure{
	my $chan = shift;
	my $PDB_file = shift;
	
	if($PDB_file =~ /gz$/){open IN, "zcat $PDB_file | " or die "Can not open $PDB_file !!";}
	else{open IN, $PDB_file or die "Can not open $PDB_file !!";}
	my %moelcule;
	while(<IN>){
		chomp;
		
		my @l = split(" ");
		
		if($l[0] ne "ATOM"){next;}
		if($l[4] ne $chan){next;}

		my @pos;
		my $AA_pos;
		if($l[4] =~ /[0-9]/){
			$AA_pos = $l[4];
			$AA_pos =~ s/[A-Z]//g;
			@pos = ($l[3], $l[5], $l[6], $l[7]);
		}
		else{
			$AA_pos = $l[5];
			@pos = ($l[3], $l[6], $l[7], $l[8]);
		}
		
		push(@{$moelcule{$AA_pos}}, \@pos);
	}

	my %moelcule_ave_pos;
	foreach my $pos (sort {$a <=> $b} keys %moelcule){
		my ($sum_x, $sum_y, $sum_z, $num, $AA) = (0, 0, 0, 0, "");
		
		foreach(@{$moelcule{$pos}}){
			my @pos = @{$_};
			$AA = $pos[0];
			$sum_x += $pos[1];
			$sum_y += $pos[2];
			$sum_z += $pos[3];
			$num++;
		}
		
		@{$moelcule_ave_pos{$pos}} = ($AA, &round($sum_x/$num, 3), &round($sum_y/$num, 3), &round($sum_z/$num, 3));
	}
	
	return %moelcule_ave_pos;
}

sub get_chan{
	my $gene_symbol = shift;
	my $PDB_file = shift;

	open IN, $PDB_file or die "Can not open $PDB_file !!";
	my %COMPND;
	my $mol_id = undef;
	my $chain = undef;
	my $synonym = undef;
	while(<IN>){
		chomp;
		if($_ =~ /MOL_ID:/){
			$mol_id = $_; $mol_id =~ s/.+MOL_ID//g; $mol_id =~ s/[^0-9]//g; ${$COMPND{$mol_id}}[0] = $mol_id;
			if(defined($COMPND{$mol_id}) != 1){@{$COMPND{$mol_id}} = ();}
		}
		if($_ =~ /CHAIN:/){$chain = $_; $chain =~ s/.+CHAIN://; $chain =~ s/\;//; ${$COMPND{$mol_id}}[1] = $chain;}
		if($_ =~ /SYNONYM:/){$synonym = $_; $synonym =~ s/.+SYNONYM: //; $synonym =~ s/\;//; ${$COMPND{$mol_id}}[2] = $synonym;}
		if($_ =~ /GENE:/){my $name = $_; $name =~ s/.+GENE: //; $name =~ s/;.+//; ${$COMPND{$mol_id}}[2] = "${$COMPND{$mol_id}}[2]"."$name";}
	}
	
	my $chan = undef;
	foreach my $id (keys %COMPND){
		my @name = split("\, ", ${$COMPND{$id}}[2]);
		foreach(@name){
			if($_ =~ /$gene_symbol/){$chan = ${$COMPND{$id}}[1];}
		}
	}
	$chan = (split(/\,/, $chan))[0];
	$chan =~ s/[^a-zA-Z]//g;
	
	if(! $chan){die"Gene symbol not found !!\n";}

	return $chan;
}

sub round{
	my $num = shift;
	my $digit = shift;
	$num = int($num*(10**$digit) + 0.5)/(10**$digit);
	return $num;
}
