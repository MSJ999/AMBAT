
# PROTEIN INSERSION BIT OF THE CODE

proc execution_pi {tx ty} {
	
	puts "	************************************************"
	puts "	THIS CODE HELPS TO ORIENT THE MEMBRANE PROTEIN IN"
	puts "	RIGHT ORIENTATION"
	puts "	************************************************"

	puts ""
	puts ""
	puts "	DEVELOPED BY TARUN KHANNA AND DR. IAN GOULD"
	puts "	IMPERIAL COLLEGE LONDON, U.K."
	puts "	FUNDED BY: EUROPEAN UNION FP7 PROGRAMME"
	
	puts ""
	puts "			#### ENTER THE NAME OF THE PDB ####"
	puts "			    (MAKE SURE TO REMOVE THE REMARK SECTION OF THE PDB)"
	set inppdb [gets stdin]
	exec ls $inppdb

	puts ""
	puts "				#### DO YOU WANT AMBAT PI TO TRY TO FIND THE BEST POSITION FOR THE PROTEIN? (Y/N) ####"
	puts ""

	set inp [gets stdin]
	set max -1000

	global bt
	set bt 30.0

	if { $inp == "y" || $inp == "Y" } {
		
		puts "			***** NOTE : BASED ON THE INPUTS YOU GIVE BELOW THE INPUT PROTEIN PDB WILL BE TRANSLATED AND ROTATED ACCORDING TO THESE VALUES *****"
		puts ""

		set inp1 $tx
		set inp2 $ty

		puts "			#### THE CODE WILL TRY A SERIES OF X AND Y ROTATION AND Z TRANSLATIONS TO INSERT THE PROTEIN ####"
		puts "			#### AT THE END CHECK 'pro_ori.pdb' TO CHECK THE CORRECTNESS OF THE ORIENTATION ####"
		puts ""

		set zstep [expr { $bt / 5 }]
		set txstep 36.0
		set tystep 36.0
		temp_grid 
		set i 0
		for {set z1 [expr { -1 * $bt }]} {$z1 < $bt} {set z1 [expr { $z1 + $zstep }]} {
			for {set tx 0.0} {$tx <= 360.0} {set tx [expr { $tx + $txstep }]} {
				for {set ty 0.0} {$ty <= 360.0} {set ty [expr { $ty + $tystep }] } {
					puts "$z1 $tx $ty"
					set inp3 $z1
					set inp4 $tx
					set inp5 $ty
					pi 0 $inppdb $inp1 $inp2 $inp3 $inp4 $inp5 0.0
					set tcheck [protein_profile_pi ref1_no.pdb 2 y $max]
					puts "		DIFF = $tcheck"
					set max $tcheck
					set arrayz1($i) $z1
					set arraytx($i) $tx
					set arraytz($i) $ty 
					set arraydiff($i) $max
					incr i
				}
			}
		}
		set npoints $i
		set max -1000
		for {set i 0} {$i < $npoints} {incr i} {
			if { $max < $arraydiff($i) } {
				set inp3 $arrayz1($i)
				set inp4 $arraytx($i)
				set inp5 $arraytz($i)
				set max $arraydiff($i)
			}
		}
		puts "			#### ENTER THE VALUE FOR Z ROTATION (NOTE: PROTEIN IS INVARIANT UNDER Z ROTATION,THIS VALUE ONLY CONTROLS THE SIMULATION SETUP) ####"
		puts ""
		set inp6 [gets stdin]
		puts "PARAMETERS :: $inp1 $inp2 $inp3 $inp4 $inp5 $inp6"
		pi 0 $inppdb $inp1 $inp2 $inp3 $inp4 $inp5 $inp6
		set tcheck [protein_profile_pi ref1_no.pdb 2 n $max]
	} else {
		puts "			***** NOTE : BASED ON THE INPUTS YOU GIVE BELOW THE INPUT PROTEIN PDB WILL BE TRANSLATED AND ROTATED ACCORDING TO THESE VALUES *****"
		puts ""

		set inp1 $tx
		set inp2 $ty

		set f [open "$inppdb" "r"]
		set data [read $f]
		close $f

		set g [open "mod_inp.pdb" "w"]

		set k 0
		set count 0
		while {$k < [llength $data] } {
			if { [lindex $data $k] == "MASTER" } {
				incr count
			}
			incr k
		}
		if { $count != 0 } {
			set strng [string first "MASTER" $data]
			set strng [expr { $strng - 1 }]
		} else {
			set strng [string first "DUM" $data]
			set strng [expr { $strng - 20 }]
		}
		set data2 [string replace $data $strng end ""]

		puts $g "$data2"
		close $g	

		set k 0

		while { [lindex $data $k] != "DUM" } {
			incr k
		}

		set z1 [lindex $data [expr { $k + 4 }]]

		set inp3 [expr { abs($z1) * -1 }]

		set inp4 0.0
		set inp5 0.0

		puts "			#### ENTER THE VALUE FOR Z ROTATION (NOTE: PROTEIN IS INVARIANT UNDER Z ROTATION, THIS VALUE ONLY CONTROLS THE SIMULATION SETUP) ####"
		puts ""
		set inp6 [gets stdin]

		temp_grid
		pi 1 "mod_inp.pdb" $inp1 $inp2 $inp3 $inp4 $inp5 $inp6
		set tcheck [protein_profile_pi ref1_no.pdb 2 n $max]
	}
}

proc pi { con inppdb inp1 inp2 inp3 inp4 inp5 inp6} {

	package require math::linearalgebra

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	set p(5) ""
	
	set p1(1) "    "
	set p1(2) "   "
	set p1(3) "  "
	set p1(4) " "
	set p1(5) ""

	set c(4) "    "
	set c(5) "   "
	set c(6) "  "
	set c(7) " "
	set c(8) ""

	set sat(1) "  "
	set sat(2) " "
	set sat(3) ""
	set sat(4) ""

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""
	
	set satn(1) "   "
	set satn(2) "  "
	set satn(3) " "
	set satn(4) " "

	# TRANSFORMATIONS

	# TRANSFORMATION 1 X AXIS

	set theta1 $inp4
	set theta1 [expr { (3.14*$theta1) / 180.0 }]

	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta1) }] [expr { -1 * (sin($theta1))}]] 
	set tmr3 [list 0.0 [expr { sin($theta1) }] [expr { cos($theta1) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta2 $inp6
	set theta2 [expr { (3.14*$theta2) / 180.0 }]

	set tm2r1 [list [expr { cos($theta2) }] [expr { -1 * sin($theta2) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta2))}] [expr { cos($theta2) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Y AXIS

	set theta3 $inp5
	set theta3 [expr { (3.14*$theta3) / 180.0 }]

	set tm3r1 [list [expr { cos($theta3) }] 0.0 [expr { sin($theta3) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta3))}] 0.0 [expr { cos($theta3) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	set f [open "$inppdb" "r"]
	set data1 [read $f]
	close $f

	set h [open "ref1_no.pdb" "w"]

	# SETTING THE ORIGIN TO THE FIRST ATOM OF FIRST RESIDUE

	set xori [lindex $data1 6]  
	set xtrans $inp1
	set yori [lindex $data1 7]
	set ytrans $inp2
	set zori [lindex $data1 8] 
	set ztrans $inp3

	set tcoord [list $xori $yori $zori]

	# TRANSFORMATION 1

	set tcoord [::math::linearalgebra::matmul $tm $tcoord]

	# TRANSFORMATION 2

	set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

	# TRANSFORMATION 3

	set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

	set dum [open "dummy" "w"]
	
	puts $dum "$tcoord"

	close $dum

	set dum [open "dummy" "r"]
	set dvar [read $dum]
	close $dum

	if { $con == 0 } {
		set xori [format "%.3f" [lindex $dvar 0]]
		set yori [format "%.3f" [lindex $dvar 1]]
		set zori [format "%.3f" [lindex $dvar 2]]
	} else {
		set xori 0.0
		set yori 0.0
		set zori 0.0
	}
				
	set k 0

	while { $k < [llength $data1] } {
		set term [lindex $data1 $k]
		set t1 [string range $term 0 5]
		if { $term == "TER" } {
			puts $h "TER"
		}
		if { [lindex $data1 $k] == "ATOM" || $t1 == "HETATM" } {
			if { $t1 == "HETATM" } {
				set sterm [string length $term]
				if { $sterm > 6 } {
					set shift 1
					set rn1 $term
					set srn1 5
					set ft ""
				} else {
						set shift 0
						set rn1 [lindex $data1 [expr { $k + 1 }]]
						set srn1 [string length $rn1]
						set ft "HETATM"
				}
			} else { 
					set shift 0
					set rn1 [lindex $data1 [expr { $k + 1 }]]
					set srn1 [string length $rn1]
					set ft "ATOM  "
			}
		
			set atype [lindex $data1 [expr { $k + 2 - $shift}]]
			set satype [string length $atype]

			if { $satype > 5} {
				set shift2 1
				set at1 [lindex $data1 [expr { $k + 2 - $shift }]]		
				set sat1 4
				set resn ""
				set sresn 4
			} else {
				set shift2 0
				set at1 [lindex $data1 [expr { $k + 2 - $shift }]]		
				set sat1 [string length $at1]
				set resn [lindex $data1 [expr { $k + 3 - $shift}]]
				set sresn [string length $resn]
			}			

			set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
			set schain_id [string length $chain_id]
			if { $schain_id > 1 } {
				set shift1 1
				set an1 [lindex $data1 [expr { $k + 4 - $shift - $shift2 }]]
				set san1 4
				set cn ""
			} else { 
				set cn [lindex $data1 [expr { $k + 4 - $shift - $shift2  }]]
				set an1 [lindex $data1 [expr { $k + 5 -$shift - $shift2 }]]
				set san1 [string length $an1]
				set shift1 0
			}

			set x1 [lindex $data1 [expr { $k + 6 - $shift - $shift1 -$shift2}]]
			set sx1 [string length $x1]
			if { $sx1 > 8 } {
				set shift3 1
				set t 0
				while { [string range $x1 $t $t] != "." } {
						incr t
				}
				set corx [string range $x1 0 [expr { $t + 3 }]]
				set cory [string range $x1 [expr { $t + 4 }] end]
				set sx1 [string length [string range $x1 0 [expr { $t + 3 }]]]
				set y1 ""
				set sy1 8
				set scory [string length $cory]
				if { $scory > 8 } {
					set y2 $cory
					set shift4 1
					set t 0
					while { [string range $y2 $t $t] != "." } {
						incr t
					}
					set cory [string range $y2 0 [expr { $t + 3 }]]
					set corz [string range $y2 [expr { $t + 4 }] end]
					set sy1 [string length [string range $y2 0 [expr { $t + 3 }]]]
					set z1 ""
					set sz1 8
				} else {
					set shift4 0
					set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
					set z1 [format "%.3f" [expr { $z1 - 0.0 }]]
					set corz $z1
					set sz1 [string length $z1]
				}
			} else { 
				set shift3 0
				set x1 [format "%.3f" [expr { $x1 - 0.0 }]]
				set corx $x1
				set sx1 [string length $x1]
				set y1 [lindex $data1 [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]]
				set cory $y1 
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set shift4 1
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set cory [string range $y1 0 [expr { $t + 3 }]]
					set corz [string range $y1 [expr { $t + 4 }] end]
					set sy1 [string length [string range $y1 0 [expr { $t + 3 }]]]
					set z1 ""
					set sz1 8
				} else {
					set shift4 0
					set y1 [format "%.3f" [expr { $y1 - 0.0 }]]
					set cory $y1
					set sy1 [string length $y1]
					set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
					set corz $z1
					set z1 [format "%.3f" [expr { $z1 - 0.0 }]]		
					set sz1 [string length $z1]
				}
			}

			if { [string length [lindex $data1 [expr { $k + 9 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]]] > 5 } {
				set shift5 1
			} else {
				set shift5 0
			}

			if { $shift4 == 0 && $shift3 == 0 } {
				set tcoord [list $x1 $y1 $z1]

				# TRANSFORMATION 1

				set tcoord [::math::linearalgebra::matmul $tm $tcoord]

				# TRANSFORMATION 2

				set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

				# TRANSFORMATION 3

				set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

				set dum [open "dummy" "w"]
	
				puts $dum "$tcoord"

				close $dum

				set dum [open "dummy" "r"]
				set dvar [read $dum]
				close $dum

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori + $xtrans }]] 
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori + $ytrans }]]
				set sy1 [string length $y1]

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori + $ztrans }]]
				set sz1 [string length $z1]
			}

			if { $shift3 != 0 } {
				set tcoord [list $corx $cory $z] 
		
				# TRANSFORMATION 1

				set tcoord [::math::linearalgebra::matmul $tm $tcoord]

				# TRANSFORMATION 2

				set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

				# TRANSFORMATION 3

				set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

				set dum [open "dummy" "w"]
	
				puts $dum "$tcoord"

				close $dum

				set dum [open "dummy" "r"]
				set dvar [read $dum]
				close $dum

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori + $xtrans}]]
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori + $ytrans }]]
				set sy1 8
	
				set x1 $x1$y1
				set y1 ""

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori + $ztrans }]]
				set sz1 [string length $z1]
			}

			if { $shift4 != 0 } {
				set tcoord [list $x $cory $corz] 
		
				# TRANSFORMATION 1

				set tcoord [::math::linearalgebra::matmul $tm $tcoord]

				# TRANSFORMATION 2

				set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

				# TRANSFORMATION 3

				set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

				set dum [open "dummy" "w"]
	
				puts $dum "$tcoord"

				close $dum

				set dum [open "dummy" "r"]
				set dvar [read $dum]
				close $dum

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori + $xtrans }]]
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori + $ytrans }]]
				set sy1 [string length $y1]

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori + $ztrans }]]
				set sz1 9

				set y1 $y1$z1
				set z1 ""
			}
		
			puts $h "$ft$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$p($sresn)$resn $cn$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data1 [expr { $k + 11 - $shift - $shift1 - $shift2 - $shift3 - $shift4 - $shift5 }]]"
		}
		incr k
	}
	#puts $h "END"
	close $h

	set h [open "ref1_no.pdb" "r"]
	set data2 [read $h]
	close $h

	set m [open "temp_grid.pdb" "r"]
	set data3 [read $m]
	close $m

	set data4 $data2$data3

	set n [open "pro_ori.pdb" "w"]
	puts $n "$data4"
	close $n
}	

proc protein_spread_pi {inp} {	

	set f [open "$inp" "r"]
	set data1 [read $f]
	close $f

	set k 0
	set minx 1000.0
	set maxx -1000.0
	set miny 1000.0
	set maxy -1000.0
	set minz 1000.0
	set maxz -1000.0

	while { $k < [llength $data1] } {
		set term [lindex $data1 $k]
		set t1 [string range $term 0 5]
		if { [lindex $data1 $k] == "ATOM" || $t1 == "HETATM" } {
			if { $t1 == "HETATM" } {
				set sterm [string length $term]
				if { $sterm > 6 } {
					set shift 1
				} else {
						set shift 0
				}
			} else { 
					set shift 0
			}
		
			set atype [lindex $data1 [expr { $k + 2 - $shift}]]
			set satype [string length $atype]

			if { $satype > 5} {
				set shift2 1
			} else {
				set shift2 0
			}			

			set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
			set schain_id [string length $chain_id]
			if { $schain_id > 1 } {
				set shift1 1
			} else { 
				set shift1 0
			}


			set x1 [lindex $data1 [expr { $k + 6 - $shift - $shift1 -$shift2}]]
			set sx1 [string length $x1]
			if { $sx1 > 8 } {
				set shift3 1
				set t 0
				while { [string range $x1 $t $t] != "." } {
					incr t
				}
				set corx [string range $x1 0 [expr { $t + 3 }]]
				set cory [string range $x1 [expr { $t + 4 }] end]
				set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3}]] 
			} else { 
				set shift3 0
				set corx $x1
				set y1 [lindex $data1 [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]] 
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set shift4 1
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set cory [string range $y1 0 [expr { $t + 3 }]]
					set corz [string range $y1 [expr { $t + 4 }] end]
				} else {
					set shift4 0
					set cory $y1
					set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
					set corz $z1
				}
			}
			if { $corx < $minx } {
				set minx $corx
			} 
			if { $cory < $miny } {
				set miny $cory
			}
			if { $corz < $minz } {	
				set minz $corz
			}
			if { $corx > $maxx } {
				set maxx $corx
			} 
			if { $cory > $maxy } {
				set maxy $cory
			}
			if { $corz > $maxz } {	
				set maxz $corz
			}
		}
		incr k
	}
	#puts "			*** THE SPREAD OF THE PROTEIN ALONG X AXIS IS [expr { $maxx - $minx }] ***"
	#puts "			*** THE SPREAD OF THE PROTEIN ALONG Y AXIS IS [expr { $maxy - $miny }] ***"
	#puts "			*** THE SPREAD OF THE PROTEIN ALONG Z AXIS IS [expr { $maxz - $minz }] ***"
	set g [open "dummy" "w"]
	puts $g " { $minx $maxx } "
	puts $g " { $miny $maxy } "
	puts $g " { $minz $maxz } "
	puts $g " { [expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }] }"
	close $g
}

proc protein_profile_pi { inp dir tcheck maxdiff } {

	protein_spread_pi $inp
	
	#puts ""
	#puts "				**** DETERIMING THE PROTEIN PROFILE ****"
	#puts ""

	set g [open "dummy" "r"]
	set dum [read $g]
	close $g

	set xext [lindex $dum 3 0]
	set yext [lindex $dum 3 1]
	set zext [lindex $dum 3 2]

	if { $dir == 0} {
		set len $xext
		set min [lindex $dum 0 0]
		set max [lindex $dum 0 1]
		set dir 0
	}
	if { $dir == 1} {
		set len $yext
		set min [lindex $dum 1 0]
		set max [lindex $dum 1 1]
		set dir 1
	}
	if { $dir == 2} {
		set len $zext
		set min [lindex $dum 2 0]
		set max [lindex $dum 2 1]
		set dir 2
	}
	#puts "						(x,y,z) = (0,1,2) "
	set f1 [open "dir1" "w"]
	puts $f1 "# dir $dir"
	
	set f [open "$inp" "r"]
	set data1 [read $f]
	close $f
	
	set step 1.0

	set f2 [open "amino_acid_dis" "w"]
	
	set f3 [open "ac_nature" "w"]

	set f4 [open "res_num" "w"]

	for {set ij $min} {$ij < $max} {set ij [expr { $ij + $step }]} {

		set oldresnum -1

		# AMINO ACID VARIABLES

		for {set i 0} {$i < 20} {incr i} {
			set ac($i) 0
		}
		for {set i 0} {$i < 5} {incr i} {
			set gp($i) 0
		}

		#puts "				**** $ij ALONG DIRECTION $dir ****"
		set k 0
		set minx 1000.0
		set maxx -1000.0
		set miny 1000.0
		set maxy -1000.0
		set minz 1000.0
		set maxz -1000.0
		while { $k < [llength $data1] } {
			set term [lindex $data1 $k]
			set t1 [string range $term 0 5]
			if { [lindex $data1 $k] == "ATOM" || $t1 == "HETATM" } {
				if { $t1 == "HETATM" } {
					set sterm [string length $term]
					if { $sterm > 6 } {
						set shift 1
					} else {
							set shift 0
					}
				} else { 
						set shift 0
				}
		
				set atype [lindex $data1 [expr { $k + 2 - $shift}]]
				set satype [string length $atype]

				if { $satype > 5} {
					set shift2 1
					set at1 [lindex $data1 [expr { $k + 2 - $shift }]]		
					set resn [string range $at1 1 3]
				} else {
					set shift2 0
					set at1 [lindex $data1 [expr { $k + 2 - $shift }]]		
					set resn [lindex $data1 [expr { $k + 3 - $shift}]]
				}			
				
				set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
				set schain_id [string length $chain_id]
				if { $schain_id > 1 } {
					set shift1 1
				} else { 
					set shift1 0
				}


				set x1 [lindex $data1 [expr { $k + 6 - $shift - $shift1 -$shift2}]]
				set sx1 [string length $x1]
				if { $sx1 > 8 } {
					set shift3 1
					set t 0
					while { [string range $x1 $t $t] != "." } {
						incr t
					}
					set corx [string range $x1 0 [expr { $t + 3 }]]
					set cory [string range $x1 [expr { $t + 4 }] end]
					set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3}]] 
				} else { 
					set shift3 0
					set corx $x1
					set y1 [lindex $data1 [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]] 
					set sy1 [string length $y1]
					if { $sy1 > 8 } {
						set shift4 1
						set t 0
						while { [string range $y1 $t $t] != "." } {
							incr t
						}
						set cory [string range $y1 0 [expr { $t + 3 }]]
						set corz [string range $y1 [expr { $t + 4 }] end]
					} else {
						set shift4 0
						set cory $y1
						set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
						set corz $z1
					}
				}
				if { $dir == 0 } {
					if { $corx >= [expr { $ij - $step }] && $corx <= $ij } {
						if { $cory < $miny } {
							set miny $cory
						}
						if { $corz < $minz } {	
							set minz $corz
						}
						if { $cory > $maxy } {
							set maxy $cory
						}
						if { $corz > $maxz } {	
							set maxz $corz
						}
					}
				}

				if { $dir == 1 } {
					if { $cory >= [expr { $ij - $step }] && $cory <= $ij } {
						if { $corx < $minx } {
							set minx $corx
						}
						if { $corz < $minz } {	
							set minz $corz
						}
						if { $corx > $maxx } {
							set maxx $corx
						}
						if { $corz > $maxz } {	
							set maxz $corz
						}
					}
				}

				if { $dir == 2 } {
					if { $corz >= [expr { $ij - $step }] && $corz <= $ij } {
						if { $corx < $minx } {
							set minx $corx
						}
						if { $cory < $miny } {	
							set miny $cory
						}
						if { $cory > $maxy } {
							set maxy $cory
						}
						if { $corx > $maxx } {	
							set maxx $corx
						}

						# DETERMING THE RESIDUE DISTRIBUTION

						if { $oldresnum != [lindex $data1 [expr { $k + 5 - $shift - $shift1 -$shift2}]] } { 

							# HYDROPHOBIC RESIDUES :: GROUP 1 :: MEMBRANE CORE

							if { $resn == "VAL" } {
								incr ac(0)
								incr gp(0)
							}
							if { $resn == "LEU" } {
								incr ac(1)
								incr gp(0)
							}
							if { $resn == "ILE" } {
								incr ac(2)
								incr gp(0)
							}
							if { $resn == "MET" } {
								incr ac(3)
								incr gp(0)
							}
							if { $resn == "PRO" } {
								incr ac(4)
								incr gp(0)
							}
							if { $resn == "PHE" } {
								incr ac(5)
								incr gp(0)
							}

							# CHARGED AND LARGE POLAR RESIDUES :: GROUP 2 :: MEMBRANE WATER INTERFACE

							if { $resn == "ARG" } {
								incr ac(6)
								incr gp(1)
							}
							if { $resn == "LYS" } {
								incr ac(7)
								incr gp(1)
							}
							if { $resn == "ASP" } {
								incr ac(8)
								incr gp(1)
							}
							if { $resn == "ASN" } {
								incr ac(9)
								incr gp(1)
							}
							if { $resn == "GLU" } {
								incr ac(10)
								incr gp(1)
							}
							if { $resn == "GLN" } {
								incr ac(11)
								incr gp(1)
							}

							# SMALL POLAR RESIDUES :: GROUP 3 :: INTERFACE AS WELL AS WITHIN BILAYER

							if { $resn == "SER" } {
								incr ac(12)
								incr gp(2)
							}
							if { $resn == "THR" } {
								incr ac(13)
								incr gp(2)
							}
						
							# OTHER CHARGED RESIDUES :: GROUP 4  :: BOUNDARY BETWEEN HYDROPHOBIC CORE AND HEADGROUP REGION

							if { $resn == "HIS" || $resn == "HIE" || $resn == "HID" } {
								incr ac(14)
								incr gp(3)
							}
							if { $resn == "TYR" } {
								incr ac(15)	
								incr gp(3)
							}
							if { $resn == "TRP" } {
								incr ac(16)
								incr gp(3)
							}

							# OTHER SMALL RESIDUES :: GROUP 5 :: FOUND THROUGHOUT THE BILAYER

							if { $resn == "GLY" } {
								incr ac(17)
								incr gp(4)
							}
							if { $resn == "ALA" } {
								incr ac(18)
								incr gp(4)
							}
							if { $resn == "CYS" } {
								incr ac(19)
								incr gp(4)
							}					
						}	
					}
				}
				set oldresnum [lindex $data1 [expr { $k + 5 - $shift - $shift1 -$shift2}]] 
			}
			incr k
		}
	puts $f1 " $ij [expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }]"

	puts $f2 " [format %.3f $ij] $ac(0) $ac(1) $ac(2) $ac(3) $ac(4) $ac(5) $ac(6) $ac(7) $ac(8) $ac(9) $ac(10) $ac(11) $ac(12) $ac(13) $ac(14) $ac(15) $ac(16) $ac(17) $ac(18) $ac(19)"

	puts $f3 " [format %.3f $ij] $gp(0) $gp(1) $gp(2) $gp(3) $gp(4)"
	}
	close $f1
	close $f2
	close $f3
	close $f4
	if { $tcheck == "y" } {
		set parm [trans]
		return $parm
	} else {
		return 0
	}
}

proc temp_grid {} {
	# TO FORM A TEMPORARY GRID

	set f [open "temp_grid.pdb" "w"]

	set tg {HETATM   76  O21 PC      1       0.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      2       0.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      3       0.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      4       0.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      5       0.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      6       0.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      7       0.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      8       0.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC      9       0.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     10       0.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     11       0.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     12       0.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     13       0.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     14       0.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     15       0.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     16       0.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     17       0.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     18       0.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     19       8.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     20       8.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     21       8.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     22       8.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     23       8.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     24       8.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     25       8.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     26       8.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     27       8.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     28       8.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     29       8.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     30       8.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     31       8.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     32       8.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     33       8.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     34       8.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     35       8.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     36       8.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     37      16.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     38      16.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     39      16.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     40      16.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     41      16.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     42      16.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     43      16.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     44      16.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     45      16.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     46      16.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     47      16.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     48      16.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     49      16.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     50      16.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     51      16.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     52      16.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     53      16.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     54      16.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     55      25.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     56      25.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     57      25.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     58      25.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     59      25.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     60      25.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     61      25.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     62      25.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     63      25.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     64      25.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     65      25.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     66      25.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     67      25.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     68      25.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     69      25.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     70      25.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     71      25.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     72      25.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     73      33.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     74      33.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     75      33.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     76      33.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     77      33.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     78      33.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     79      33.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     80      33.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     81      33.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     82      33.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     83      33.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     84      33.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     85      33.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     86      33.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     87      33.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     88      33.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     89      33.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     90      33.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     91      41.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     92      41.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     93      41.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     94      41.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     95      41.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     96      41.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     97      41.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     98      41.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC     99      41.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    100      41.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    101      41.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    102      41.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    103      41.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    104      41.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    105      41.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    106      41.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    107      41.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    108      41.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    109      50.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    110      50.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    111      50.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    112      50.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    113      50.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    114      50.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    115      50.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    116      50.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    117      50.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    118      50.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    119      50.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    120      50.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    121      50.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    122      50.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    123      50.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    124      50.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    125      50.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    126      50.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    127      58.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    128      58.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    129      58.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    130      58.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    131      58.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    132      58.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    133      58.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    134      58.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    135      58.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    136      58.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    137      58.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    138      58.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    139      58.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    140      58.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    141      58.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    142      58.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    143      58.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    144      58.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    145      66.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    146      66.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    147      66.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    148      66.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    149      66.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    150      66.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    151      66.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    152      66.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    153      66.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    154      66.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    155      66.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    156      66.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    157      66.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    158      66.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    159      66.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    160      66.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    161      66.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    162      66.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    163      75.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    164      75.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    165      75.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    166      75.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    167      75.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    168      75.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    169      75.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    170      75.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    171      75.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    172      75.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    173      75.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    174      75.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    175      75.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    176      75.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    177      75.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    178      75.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    179      75.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    180      75.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    181      83.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    182      83.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    183      83.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    184      83.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    185      83.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    186      83.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    187      83.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    188      83.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    189      83.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    190      83.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    191      83.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    192      83.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    193      83.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    194      83.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    195      83.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    196      83.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    197      83.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    198      83.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    199      91.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    200      91.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    201      91.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    202      91.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    203      91.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    204      91.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    205      91.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    206      91.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    207      91.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    208      91.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    209      91.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    210      91.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    211      91.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    212      91.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    213      91.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    214      91.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    215      91.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    216      91.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    217     100.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    218     100.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    219     100.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    220     100.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    221     100.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    222     100.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    223     100.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    224     100.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    225     100.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    226     100.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    227     100.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    228     100.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    229     100.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    230     100.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    231     100.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    232     100.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    233     100.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    234     100.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    235     108.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    236     108.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    237     108.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    238     108.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    239     108.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    240     108.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    241     108.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    242     108.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    243     108.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    244     108.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    245     108.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    246     108.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    247     108.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    248     108.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    249     108.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    250     108.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    251     108.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    252     108.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    253     116.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    254     116.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    255     116.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    256     116.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    257     116.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    258     116.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    259     116.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    260     116.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    261     116.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    262     116.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    263     116.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    264     116.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    265     116.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    266     116.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    267     116.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    268     116.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    269     116.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    270     116.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    271     125.000   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    272     125.000   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    273     125.000  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    274     125.000  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    275     125.000  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    276     125.000  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    277     125.000  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    278     125.000  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    279     125.000  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    280     125.000  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    281     125.000  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    282     125.000  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    283     125.000 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    284     125.000 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    285     125.000 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    286     125.000 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    287     125.000 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    288     125.000 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    289     133.333   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    290     133.333   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    291     133.333  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    292     133.333  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    293     133.333  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    294     133.333  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    295     133.333  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    296     133.333  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    297     133.333  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    298     133.333  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    299     133.333  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    300     133.333  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    301     133.333 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    302     133.333 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    303     133.333 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    304     133.333 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    305     133.333 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    306     133.333 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    307     141.667   0.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    308     141.667   8.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    309     141.667  16.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    310     141.667  25.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    311     141.667  33.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    312     141.667  41.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    313     141.667  50.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    314     141.667  58.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    315     141.667  66.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    316     141.667  75.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    317     141.667  83.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    318     141.667  91.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    319     141.667 100.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    320     141.667 108.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    321     141.667 116.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    322     141.667 125.000   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    323     141.667 133.333   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    324     141.667 141.667   0.000  1.00  0.00           O
ter
HETATM   76  O21 PC    325       0.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    326       0.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    327       0.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    328       0.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    329       0.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    330       0.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    331       0.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    332       0.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    333       0.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    334       0.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    335       0.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    336       0.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    337       0.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    338       0.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    339       0.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    340       0.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    341       0.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    342       0.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    343       8.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    344       8.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    345       8.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    346       8.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    347       8.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    348       8.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    349       8.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    350       8.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    351       8.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    352       8.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    353       8.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    354       8.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    355       8.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    356       8.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    357       8.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    358       8.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    359       8.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    360       8.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    361      16.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    362      16.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    363      16.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    364      16.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    365      16.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    366      16.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    367      16.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    368      16.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    369      16.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    370      16.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    371      16.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    372      16.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    373      16.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    374      16.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    375      16.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    376      16.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    377      16.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    378      16.667 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    379      25.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    380      25.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    381      25.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    382      25.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    383      25.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    384      25.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    385      25.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    386      25.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    387      25.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    388      25.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    389      25.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    390      25.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    391      25.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    392      25.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    393      25.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    394      25.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    395      25.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    396      25.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    397      33.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    398      33.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    399      33.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    400      33.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    401      33.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    402      33.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    403      33.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    404      33.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    405      33.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    406      33.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    407      33.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    408      33.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    409      33.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    410      33.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    411      33.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    412      33.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    413      33.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    414      33.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    415      41.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    416      41.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    417      41.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    418      41.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    419      41.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    420      41.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    421      41.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    422      41.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    423      41.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    424      41.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    425      41.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    426      41.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    427      41.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    428      41.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    429      41.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    430      41.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    431      41.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    432      41.667 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    433      50.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    434      50.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    435      50.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    436      50.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    437      50.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    438      50.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    439      50.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    440      50.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    441      50.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    442      50.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    443      50.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    444      50.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    445      50.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    446      50.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    447      50.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    448      50.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    449      50.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    450      50.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    451      58.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    452      58.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    453      58.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    454      58.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    455      58.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    456      58.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    457      58.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    458      58.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    459      58.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    460      58.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    461      58.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    462      58.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    463      58.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    464      58.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    465      58.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    466      58.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    467      58.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    468      58.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    469      66.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    470      66.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    471      66.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    472      66.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    473      66.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    474      66.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    475      66.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    476      66.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    477      66.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    478      66.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    479      66.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    480      66.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    481      66.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    482      66.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    483      66.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    484      66.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    485      66.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    486      66.667 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    487      75.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    488      75.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    489      75.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    490      75.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    491      75.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    492      75.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    493      75.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    494      75.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    495      75.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    496      75.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    497      75.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    498      75.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    499      75.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    500      75.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    501      75.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    502      75.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    503      75.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    504      75.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    505      83.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    506      83.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    507      83.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    508      83.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    509      83.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    510      83.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    511      83.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    512      83.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    513      83.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    514      83.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    515      83.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    516      83.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    517      83.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    518      83.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    519      83.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    520      83.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    521      83.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    522      83.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    523      91.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    524      91.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    525      91.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    526      91.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    527      91.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    528      91.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    529      91.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    530      91.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    531      91.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    532      91.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    533      91.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    534      91.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    535      91.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    536      91.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    537      91.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    538      91.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    539      91.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    540      91.667 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    541     100.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    542     100.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    543     100.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    544     100.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    545     100.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    546     100.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    547     100.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    548     100.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    549     100.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    550     100.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    551     100.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    552     100.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    553     100.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    554     100.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    555     100.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    556     100.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    557     100.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    558     100.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    559     108.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    560     108.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    561     108.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    562     108.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    563     108.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    564     108.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    565     108.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    566     108.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    567     108.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    568     108.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    569     108.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    570     108.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    571     108.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    572     108.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    573     108.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    574     108.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    575     108.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    576     108.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    577     116.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    578     116.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    579     116.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    580     116.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    581     116.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    582     116.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    583     116.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    584     116.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    585     116.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    586     116.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    587     116.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    588     116.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    589     116.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    590     116.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    591     116.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    592     116.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    593     116.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    594     116.667 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    595     125.000   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    596     125.000   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    597     125.000  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    598     125.000  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    599     125.000  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    600     125.000  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    601     125.000  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    602     125.000  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    603     125.000  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    604     125.000  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    605     125.000  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    606     125.000  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    607     125.000 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    608     125.000 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    609     125.000 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    610     125.000 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    611     125.000 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    612     125.000 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    613     133.333   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    614     133.333   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    615     133.333  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    616     133.333  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    617     133.333  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    618     133.333  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    619     133.333  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    620     133.333  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    621     133.333  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    622     133.333  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    623     133.333  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    624     133.333  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    625     133.333 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    626     133.333 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    627     133.333 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    628     133.333 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    629     133.333 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    630     133.333 141.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    631     141.667   0.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    632     141.667   8.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    633     141.667  16.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    634     141.667  25.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    635     141.667  33.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    636     141.667  41.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    637     141.667  50.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    638     141.667  58.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    639     141.667  66.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    640     141.667  75.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    641     141.667  83.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    642     141.667  91.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    643     141.667 100.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    644     141.667 108.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    645     141.667 116.667 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    646     141.667 125.000 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    647     141.667 133.333 -36.000  1.00  0.00           O
TER
HETATM   76  O21 PC    648     141.667 141.667 -36.000  1.00  0.00           O
TER
END}

	puts $f "$tg"
	close $f
}

proc trans {} {

	global bt

	set f [open "ac_nature" "r"]
	set data [read $f]
	close $f

	set k 0
	
	set i 0
	while { $k < [llength $data] } {
		incr i
		incr k 6
	}
	set num_dp $i

	# COUNTING THE HYDROPHOBIC AND HYDROPHYLIC RESIDUES

	set hpho 0
	set hphy 0
	set hphoo 0
	set hphyo 0
	set gp4 0
	for {set i 0} {$i < $num_dp} {incr i} {
		set znew [lindex $data [expr { 6 * $i }]]
		if { $znew < 0.0 && $znew > [expr { $bt * -1.0 }] } {
			set hpho [expr { $hpho + [lindex $data [expr { (6*$i) + 1 }]] }]
			set hphy [expr { $hphy + [lindex $data [expr { (6*$i) + 2 }]] }]
		} else {
			set hphoo [expr { $hphoo + [lindex $data [expr { (6*$i) + 1 }]] }]
			set hphyo [expr { $hphyo + [lindex $data [expr { (6*$i) + 2 }]] }]
		}
		if { $znew < 5.0 && $znew > 0.0 } {
			set gp4 [expr { $gp4 + [lindex $data [expr { (6*$i) + 4 }]] }]
		}
		if { $znew < [expr { $bt * -1.0 }] && $znew > [expr {($bt + 5.0) * -1.0 }] } {
			set gp4 [expr { $gp4 + [lindex $data [expr { (6*$i) + 4 }]] }]
		}
	}
	set diff [expr { $hpho - $hphy }]
	set diffo [expr { $hphyo - $hphoo }]
	#if {$diff > $max } {
		#set max $diff
	#}
	puts "GP4 == $gp4"
	return [expr { $diff + $diffo + (1*$gp4)}]
}

execution_pi 70.0 70.0
