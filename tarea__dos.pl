#!/usr/bin/perl -w

# prog1.1 
# Bruno Contreras-Moreira
# Nearest Neighbor dG calculator

use strict;

# global variables
my $T           = 37; # temperature(C)
my $windowL     = 15;  # window length, http://www.biomedcentral.com/1471-2105/6/1
my %NNparams    = ( 
	# SantaLucia J (1998) PNAS 95(4): 1460-1465.
	# [NaCl] 1M, 37C & pH=7 
	# H(enthalpy): kcal/mol	, S(entropy): cal/k�mol
	# stacking dinucleotides
	'AA/TT' , {'H',-7.9, 'S',-22.2},
	'AT/TA' , {'H',-7.2, 'S',-20.4},
	'TA/AT' , {'H',-7.2, 'S',-21.3},
	'CA/GT' , {'H',-8.5, 'S',-22.7},
	'GT/CA' , {'H',-8.4, 'S',-22.4},
	'CT/GA' , {'H',-7.8, 'S',-21.0},
	'GA/CT' , {'H',-8.2, 'S',-22.2},
	'CG/GC' , {'H',-10.6,'S',-27.2},
	'GC/CG' , {'H',-9.8, 'S',-24.4},
	'GG/CC' , {'H',-8.0, 'S',-19.9},
	# initiation costs
	'G'     , {'H', 0.1, 'S',-2.8 },
	'A'     , {'H', 2.3, 'S',4.1  },
	# symmetry correction
	'sym'   , {'H', 0, 'S',-1.4 } );

my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";

print "# parameters: Temperature=$T Window=$windowL\n\n";

my @DGs = ();	#lista 2D en donde guardo los deltaG de cada ventana por secuencia
my $index = 0;
		
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>)
{
	if(/^(b\d{4}) \\ ([ATGC]+)/)
	{
		my ($name,$seq) = ($1,$2); 
		#printf("sequence %s (%d nts)\n",$name,length($seq));
		my @secs = split("", $seq); #guardo mis secuencias en array	
		foreach my $ren (14..$#secs) {
			my $ini = $ren - 14;
			my @cuts = @secs;  
			my @nucs = splice(@cuts, $ini, 15); #corto ventanas de 15nuc
			my $reads = join("",@nucs); #guardo las ventanas en la variable reads
			#my $DGs = duplex_deltaG($reads,37);
			#print "$DGs\n";	
			$DGs[$index][$ini] = duplex_deltaG($reads,37);#guardolo anterios en array2D		
						}
	}
	$index++;		
};
close(SEQ);

#abro lista2D DGs	
foreach my $elem(@DGs){
		my @values = @{$elem};
		#defino los limites para carcular E1 Y E2
		foreach my $val(199..$#values) {
			my @dgsE2 = @values[$val-100..$val];#corto fracciones de 100 valore de deltag
			my @dgsE1 = @values[$val-199..$val-150];#corto fracciones de 50valores de deltag
		my $E2 = 0;				
		foreach my $reb(@dgsE2){
		$E2 += $reb;
		#sumatoria de cada fraccion
		}
		$E2 /= 100;
		#promedio de cada fracción
		my $E1 = 0;
		foreach my $reb(@dgsE1){
		$E1 += $reb;
		#sumatoria de cada fraccion
		}
		$E1 /= 50;
		#promedio de cada fracción
		#calculo D y defino los limites o cutoff 
		my $D = $E1 - $E2;
		if(($D>3.3) && ($E1> -17.1)){
		print "$D\t$E1\t$val\t$elem\n"; #imprimo n's positivas = 499 promotores, 238n's por secuencia
		}
		print "$E1\t$E2\t$D\t$val\t$elem\n";#imprimo tabla E1 y D por secuencia		
	}
}

sub duplex_deltaG {	
	my ($reads, $T) = @_; 			#parámetros
	my ($DNAstep,$nt,$dG,$total_dG) = ('','',0,0);	#variables locales
	my @sequence = split(//,uc($reads)); 	#guardo los reads en array separados por espacio
	my $tK = 273.15 + $T; 			#convierto temperatura celcius en kelvin	
	sub complement{ $_[0] =~ tr/ATGC/TACG/; return $_[0] } #saca el complemento de la secuencia
	# add dG for overlapping dinucleotides
	for(my $n=0;$n<$#sequence;$n++) {
			$DNAstep = $sequence[$n].$sequence[$n+1].'/'.
				complement($sequence[$n].$sequence[$n+1]);
			if(!defined($NNparams{$DNAstep}))
			{$DNAstep = reverse($DNAstep)}
			$dG = (((1000*$NNparams{$DNAstep}{'H'})-
					($tK*$NNparams{$DNAstep}{'S'}))/ 1000);
			$total_dG += $dG; #deltaG por dinucleotido
			}
	# add correction for helix initiation
	$nt = $sequence[0]; # first pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) } 
	$total_dG += (((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))/1000);
	#añadimos corrección por ultimo nuc								 
	$nt = $sequence[$#sequence]; # last pair
	if(!defined($NNparams{$nt})){ $nt = complement($nt) }
	$total_dG += (((1000*$NNparams{$nt}{'H'})-
					($tK*$NNparams{$nt}{'S'}))/ 1000);
	#añadimos correccion si la secuancia es simetrica				
	if(complement($reads) eq reverse(complement($reads))) {
	$total_dG += (((1000*$NNparams{'sym'}{'H'})-
					($tK*$NNparams{'sym'}{'S'}))/ 1000);
								}		
return $total_dG;	#regresa el deltaG total por ventana incluyendo las correcciones		
}	
	

