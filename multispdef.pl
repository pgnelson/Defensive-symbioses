#!/usr/bin/perl
use strict;
#all parameter names are the same as in Nelson and May 2018
my $reproduce_figure = 2;
#this will reproduce either figure 1, 2 or 3, first it will find equilibria when shared costs affect b, and then equilibria when shared costs affect m
#each row of the prinout contains all parameters values and the equilibrium v, I, f, found under those parameters

#figure 4 involved running the code for a two dimensional grid of  p and v, and manually finding values of p at which virulence trasitioned from negative to positive
#as such we are unable to provide a simple configuration for reproducing figure 4

my @ab = (4, 4); #the number of values in the parameters arrays determines the number of species
my @am = (1, 1);#all parameter arrays need to have the same number of values
my @gb = (1, 1);
my @gm = (1, 1);
my @b0 = (0.4, 0);
my @m0 = (0.1, 0.1);
my @hm = (0, 0);
my @hb = (0.1, 0.1);
my $accuracy = 0.0001;#controls length of run, the smaller this number is the longer the run
my $pb =0;
my $pm =0;
my $bigM = 1; #exponent on m
my $bigB = 0.5;#exponent on b
my $type = 1;#how defence is incorporated, "1" corresponds to the functional form used in Nelson and May 2018
my $withinspecies;#1 if shared costs are within species, 0 if only between species, 2 if both
my $numsp = scalar @ab;
my $numprintouts;#number of times it prints values during a run
my $generations = 200 / $accuracy;
my $parameter_min = 0;
my $parameter_max = 5;
open (OUTPUT, '>defdata.txt');

for(my $i = 0; $i < $numsp; $i++){
	print OUTPUT "shared\tab\tam\tb0\tm0\tpb\tpm\tgb\tgm\thb\thm\tv\tf\tI\tb\tm\tdI\t";
	print "shared,ab,am,b0,m0,pb,pm,gb,gm,hb,hm,v,f,I,b,m,dI,";
}
print OUTPUT "\n";print "\n";

if($reproduce_figure == 1){$withinspecies = 1;$numprintouts=1; $parameter_min=0; $parameter_max =5}
if($reproduce_figure == 2){$withinspecies = 0;$numprintouts=1;$parameter_min=0; $parameter_max =5}
if($reproduce_figure == 3){$withinspecies = 0;$numprintouts=200;$parameter_min=5; $parameter_max =5}

for (my $shared_costs_to_birth = 0; $shared_costs_to_birth <= 1; $shared_costs_to_birth++){
	my @initial_Inf = (0.1, 0.1);
	my @initial_v = (0.1, 0.1);
	my @initial_f = (0, 0);
	for (my $j = $parameter_min; $j <= $parameter_max; $j += 0.05){ #incrementing the parameter of interest
		#$b0[0] = $j;
		if($shared_costs_to_birth==0){$pb=$j;$pm=0;}else{$pb=0;$pm=$j;}
		my @Inf = @initial_Inf;
		my @v = @initial_v;
		my @f = @initial_f;
		for(my $i = 0; $i < $numsp; $i++){	#loop to change any other parameters
			#$am[$i] = 1 - $pm;
		}
		print $shared_costs_to_birth. "\t" .$j . "\n";
		for(my $counter = 0; $counter <= $generations; $counter++){
			my $PubB = 0;my $PubM = 0;
			if ($withinspecies == 2){
				for(my $i = 0; $i < $numsp; $i++){
					$PubM += $v[abs($i)] * $Inf[abs($i)] * $pm / $numsp;
					$PubB += $v[abs($i)] * $Inf[abs($i)] * $pb / $numsp;
				}
			}
			for(my $sp = 0; $sp < $numsp; $sp++){
				my $dI;
				if ($withinspecies == 0){
					$PubM = $v[abs($sp-1)] * $Inf[abs($sp-1)] * $pm;
					$PubB = $v[abs($sp-1)] * $Inf[abs($sp-1)] * $pb;
				}elsif ($withinspecies == 1){
					$PubM = $v[abs($sp)] * $Inf[abs($sp)] * $pm;
					$PubB = $v[abs($sp)] * $Inf[abs($sp)] * $pb;
				}
				my $c = (1 - $Inf[$sp] );
				if ($c < $accuracy){$c = 1000*$accuracy;}#probability of new infections never reaches zero
				my $dI;	my $dIvplus;my $dIvminus;my $dIfplus;my $dIfminus;my $newPubB;my $newPubM;
				for (my $j = -1; $j <2; $j++){#calculates dI for vir above below and at current values
					my $testv = $v[$sp] + $j * $accuracy;
					if ($type == 0){#adjusting public costs for defence
						$newPubM = $PubM * (1 - $gm[$sp] * $f[$sp]) ;
						$newPubB = $PubB * (1 - $gb[$sp] * $f[$sp]);
					}elsif ($type == 1){
						$newPubM = $PubM /(1 + $gm[$sp] * $f[$sp]);
						$newPubB = $PubB /(1 + $gb[$sp] * $f[$sp]);
					}elsif ($type == 2){
						$newPubM = $PubM * exp(-1*$gm[$sp] * $f[$sp]);
						$newPubB = $PubB * exp(-1*$gb[$sp] * $f[$sp]);
					}
					my $m = $testv * $am[$sp] + $newPubM + $m0[$sp] + $hm[$sp]*$f[$sp];
					if ($m > 0){ $m = $m ** $bigM;}else{$m=0;}
					my $b = $testv * $ab[$sp] - $newPubB + $b0[$sp] - $hb[$sp]*$f[$sp];
					if ($b > 0){ $b = $b ** $bigB;}else{$b=0;}
					if ($j == -1) {$dIvminus = $c*$b-$m;}
					if ($j == 0) {$dI = $Inf[$sp]*($c*$b-$m);}
					if ($j == 1) {$dIvplus = $c*$b-$m;}
				}

				for (my $k = -1; $k <2; $k++){#calculates dI for def above below and at current values
					my $testf = $f[$sp] + $k * $accuracy;
					if ($type == 0){#adjusting public costs for defence
						$newPubM = $PubM * (1 - $gm[$sp] * $testf) ;
						$newPubB = $PubB * (1 - $gb[$sp] * $testf);
					}elsif ($type == 1){
						$newPubM = $PubM /(1 + $gm[$sp] * $testf);
						$newPubB = $PubB /(1 + $gb[$sp] * $testf);
					}elsif ($type == 2){
						$newPubM = $PubM * exp(-1*$gm[$sp] * $testf);
						$newPubB = $PubB * exp(-1*$gb[$sp] * $testf);
					}
					my $m = $v[$sp] * $am[$sp] + $newPubM + $m0[$sp] + $hm[$sp]*$testf;
					if ($m > 0){ $m = $m ** $bigM;}else{$m=0;}
					my $b = $v[$sp] * $ab[$sp] - $newPubB + $b0[$sp] - $hb[$sp]*$testf;
					if ($b > 0){ $b = $b ** $bigB;}else{$b=0;}
					if ($k == -1) {$dIfminus = $c*$b-$m;}
					if ($k == 1) {$dIfplus = $c*$b-$m;}
				}

				$v[$sp] += ($dIvplus - $dIvminus)*$accuracy / 0.0005;#incrementing virulence
				$f[$sp] += ($dIfplus - $dIfminus)*$accuracy / 0.0005;#incrementing defence
				$Inf[$sp] += $accuracy * $dI;	#incrementing infection frequency
				#checking for boundary conditions
				if (($gm[$sp] * $f[$sp] > 1) and $type == 0){$f[$sp] = 1 / $gm[$sp];}
				if (($gb[$sp] * $f[$sp] > 1) and $type == 0){$f[$sp] = 1 / $gb[$sp];}
				if ($Inf[$sp] < 0){$Inf[$sp]=0;}	if ($Inf[$sp] > 1){$Inf[$sp]=1;}
				if ($f[$sp] < 0){$f[$sp]=0;}
			#This section prints out the results
				if($counter > 0 and $counter%(($generations)/$numprintouts) == 0 or $counter == $generations){
					my $PubB = 0;my $PubM = 0;my $newPubB = 0;my $newPubM = 0;
					if ($withinspecies == 2){
						for(my $i = 0; $i < $numsp; $i++){
							$PubM += $v[abs($i)] * $Inf[abs($i)] * $pm / $numsp;
							$PubB += $v[abs($i)] * $Inf[abs($i)] * $pb / $numsp;
						}
					}
					for(my $i = 0; $i < $numsp; $i++){
						if ($withinspecies == 0){
							$PubM = $v[abs($i-1)] * $Inf[abs($i-1)] * $pm;
							$PubB = $v[abs($i-1)] * $Inf[abs($i-1)] * $pb;
						}elsif ($withinspecies == 1){
							$PubM = $v[abs($i)] * $Inf[abs($i)] * $pm;
							$PubB = $v[abs($i)] * $Inf[abs($i)] * $pb;
						}
						if ($type == 0){
							$newPubM = $PubM * (1 - $gm[$i] * $f[$i]);
							$newPubB = $PubB * (1 - $gb[$i] * $f[$i]);
						}elsif ($type == 1){
							$newPubM = $PubM /(1 + $gm[$i] * $f[$i]);
							$newPubB = $PubB /(1 + $gb[$i] * $f[$i]);
						}elsif ($type == 2){
							$newPubM = $PubM * exp(-1*$gm[$i] * $f[$i]);
							$newPubB = $PubB * exp(-1*$gb[$i] * $f[$i]);
						}
						my $m = $v[$i] * $am[$i] + $newPubM + $m0[$i] + $hm[$i]*$f[$i];
						if ($m > 0){ $m = $m ** $bigM;}else{$m=0;}
						my $b = $v[$i] * $ab[$i] - $newPubB + $b0[$i] - $hb[$i]*$f[$i];
						if ($b > 0){$b = $b ** $bigB;}else{$b=0;}
						my $c = (1-$Inf[$i]);
						my $dI = $c*$b-$m;
						#if (abs($dI)>0.01){$v[$i]=0;$f[$i]=0;$Inf[$i]=0;}
						printf OUTPUT ("%.2f",$shared_costs_to_birth);	print OUTPUT "\t";
						printf OUTPUT ("%.2f",$ab[$i]);	print OUTPUT "\t";printf OUTPUT ("%.2f",$am[$i]); 	print OUTPUT "\t";
						printf OUTPUT ("%.2f",$b0[$i]);	print OUTPUT "\t";printf OUTPUT ("%.2f",$m0[$i]); 	print OUTPUT "\t";
						printf OUTPUT ("%.2f",$pb);		print OUTPUT "\t";printf OUTPUT ("%.2f",$pm); 		print OUTPUT "\t";
						printf OUTPUT ("%.2f",$gb[$i]);	print OUTPUT "\t";printf OUTPUT ("%.2f",$gm[$i]); 	print OUTPUT "\t";
						printf OUTPUT ("%.2f",$hb[$i]);	print OUTPUT "\t";printf OUTPUT ("%.2f",$hm[$i]); 	print OUTPUT "\t";
						printf OUTPUT ("%.3f",$v[$i]);  print OUTPUT "\t";printf OUTPUT ("%.3f",$f[$i]);  	print OUTPUT "\t";	printf OUTPUT ("%.3f",$Inf[$i]);print OUTPUT "\t";
						printf OUTPUT ("%.3f",$b); 		print OUTPUT "\t";printf OUTPUT ("%.3f",$m); 		print OUTPUT "\t";	printf OUTPUT ("%.3f",$dI);		print OUTPUT "\t";

						printf  ("%.2f",$shared_costs_to_birth);	print ",";
						printf  ("%.2f",$ab[$i]);	print  ",";printf  ("%.2f",$am[$i]);	print  ",";
						printf  ("%.2f",$b0[$i]);	print  ",";printf  ("%.2f",$m0[$i]);	print  ",";
						printf  ("%.2f",$pb);		print  ",";printf  ("%.2f",$pm); 		print  ",";
						printf  ("%.2f",$gb[$i]);	print  ",";printf  ("%.2f",$gm[$i]); 	print  ",";
						printf  ("%.2f",$hb[$i]);	print  ",";printf  ("%.2f",$hm[$i]); 	print  ",";
						printf  ("%.3f",$v[$i]); 	print  ",";printf  ("%.3f",$f[$i]);  	print  ",";	printf  ("%.3f",$Inf[$i]); 	print  ",";
						printf  ("%.3f",$b);  		print  ",";printf  ("%.3f",$m);  		print  ",";	printf  ("%.3f",$dI);		print  ",";
						#$initial_Inf[$i] = $v[$i];#initial infection frequency for next run
						#$initial_v[$i] = $f[$i];#initial virulence for next run
						#$initial_f[$i] = $Inf[$i];#initial defence for next run
					}
					print OUTPUT "\n";
					print "\n";
				}
			}
		}
	}
}
