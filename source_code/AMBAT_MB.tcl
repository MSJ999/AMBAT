proc protein_spread {} {	

	# THIS PROCEDURE WILL CALCULATE THE PROTEIN PROFILE ALONG THE AXIS OF LARGEST SPREAD

	puts "				**** CALCULATING THE PROTEIN SPREAD, LOOK FOR DIR1 FILE IN THE FOLDER TO PLOT THE SPREAD ALONG THE LONGEST AXIS ****"
	puts ""
	
	set in [open "input" "r"]
	set inp [read $in]
	close $in

	set f [open "[lindex $inp 7 1]" "r"]
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
	puts "			*** THE SPREAD OF THE PROTEIN ALONG X AXIS IS [expr { $maxx - $minx }] ***"
	puts "			*** THE SPREAD OF THE PROTEIN ALONG Y AXIS IS [expr { $maxy - $miny }] ***"
	puts "			*** THE SPREAD OF THE PROTEIN ALONG Z AXIS IS [expr { $maxz - $minz }] ***"
	set g [open "dummy" "w"]
	puts $g " { $minx $maxx } "
	puts $g " { $miny $maxy } "
	puts $g " { $minz $maxz } "
	puts $g " { [expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }] }"
	close $g
}

proc protein_profile {} {

	protein_spread
	
	puts ""
	puts "				**** DETERIMING THE PROTEIN PROFILE ****"
	puts ""

	set g [open "dummy" "r"]
	set dum [read $g]
	close $g

	set xext [lindex $dum 3 0]
	set yext [lindex $dum 3 1]
	set zext [lindex $dum 3 2]

	if { $xext > $yext && $xext > $zext} {
		set len $xext
		set min [lindex $dum 0 0]
		set max [lindex $dum 0 1]
		set dir 0
	}
	if { $yext > $xext && $yext > $zext} {
		set len $yext
		set min [lindex $dum 1 0]
		set max [lindex $dum 1 1]
		set dir 1
	}
	if { $zext > $xext && $zext > $yext} {
		set len $zext
		set min [lindex $dum 2 0]
		set max [lindex $dum 2 1]
		set dir 2
	}

	set len $zext
	set dir 2
	set min [lindex $dum 2 0]
	set max [lindex $dum 2 1]

	puts "						(x,y,z) = (0,1,2) "
	set f1 [open "dir1" "w"]
	puts $f1 "# dir $dir"
	
	set in [open "input" "r"]
	set inp [read $in]
	close $in

	set f [open "[lindex $inp 7 1]" "r"]
	set data1 [read $f]
	close $f
	
	set step 1.0

	for {set ij $min} {$ij < $max} {set ij [expr { $ij + $step }]} {
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
				if { $dir == 0 } {
					if { $corx > [expr { $ij - $step }] && $corx < [expr { $ij + $step }] } {
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
					if { $cory > [expr { $ij - $step }] && $cory < [expr { $ij + $step }] } {
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
					if { $corz > [expr { $ij - $step }] && $corz < [expr { $ij + $step }] } {
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
					}
				}

			}
			incr k
		}
	puts "				**** $ij POSITION ALONG THE BILAYER THICKNESS EXTENTION ALONG OTHER DIRECTIONS IS[expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }] (IGNORE -2000)****"
	puts $f1 " $ij [expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }]"
	}
	puts ""
	puts "				#### ALL PUT IN FILE DIR1 DON'T WORRY ####"
	puts ""
	close $f1
}

proc lip_num {} {

	# THIS PROCEDURE WILL CALCULTE THE APPROXIMATE NUMBER OF LIPIDS WHICH WILL BE REMOVED BECAUSE OF PROTEIN INSERTION

	puts "				**** CALCULATING THE NUMBER OF LIPIDS TO BE REMOVED ****"
	puts ""

	set f [open "dir1" "r"]
	set data [read $f]
	close $f

	set g [open "input" "r"]
	set data1 [read $g]
	close $g

	set h [open "output" "w"]

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]

	set k1 0
	for {set k1 0} {$k1 < 100} {incr k1} {
		set l_ul($k1) 0
		set l_ll($k1) 0
	}

	set k1 1

	# X Y PLANE

	set num_lipids_ul 0.0
	set num_lipids_ll 0.0

	while { $k1 < [llength [lindex $data1 3]] } {
		set num_lipids_ul [expr { $num_lipids_ul + [lindex $data1 3 $k1] }]
		set num_lipids_ll [expr { $num_lipids_ll + [lindex $data1 4 $k1] }]
		set l_ul([expr { $k1 - 1 }]) [lindex $data1 3 $k1] 
		set l_ll([expr { $k1 - 1 }]) [lindex $data1 4 $k1]

		incr k1
	}

	if { $num_lipids_ul > $num_lipids_ll } {
		set num_lipids $num_lipids_ul
	} else {
		set num_lipids $num_lipids_ll
	}

	set stepx [expr { $X / (sqrt($num_lipids))  }]
	set stepy [expr { $Y / (sqrt($num_lipids))  }]

	set dir [lindex $data 2]

	set th [lindex $data1 5 3]
	set th [expr { $th * 2 }]

	set k 3
	set start $k

	while { $k < [llength $data] } {	
		set xold [lindex $data $k]
		set k1 $k
		if { $xold > [expr { (-1*$th) - 4.0 }] } {
			while { $k1 < [llength $data] } {
				set xnew [lindex $data $k1]
				set delx [expr { abs($xnew-$xold) }]
			
				if { $delx < [expr { $th + 4.0 }] && $xnew <= 4.000 } {
					set k2 $k
					set maxyul 0.0
					set maxzul 0.0
					set maxyll 0.0
					set maxzll 0.0
					while { $k2 < $k1 } {
						set xch [lindex $data $k2]
						if { $dir == 0 } { 
							set y [lindex $data [expr { $k2 + 2 }]]
							set z [lindex $data [expr { $k2 + 3 }]]
						} elseif { $dir == 1 } {
							set y [lindex $data [expr { $k2 + 1 }]]
							set z [lindex $data [expr { $k2 + 3 }]]
						} elseif { $dir == 2 } {
							set y [lindex $data [expr { $k2 + 1 }]]
							set z [lindex $data [expr { $k2 + 2 }]]
						}
						# UPPER LEAFLET
						if { [expr { abs($xch-$xold) }] < $Z } {
							if { $y > $maxyul } {
								set maxyul $y
							}
							if { $z > $maxzul } {
								set maxzul $z
							}
						# LOWER LEAFLET
						} else { 
							if { $y > $maxyll } {
								set maxyll $y
							}
							if { $z > $maxzll } {
								set maxzll $z
							}
						}
						incr k2 4
					}
					puts $h "$xold	$xnew	$maxyul	$maxzul	$maxyll	$maxzll"
				}
				incr k1 4
			}
			set k [llength $data]
		}
		incr k 4
	}
	close $h

	set h [open "output" "r"]
	set data2 [read $h]	
	close $h

	# CALCULATING THE APPROXIMATE NUMBER OF LIPIDS WHICH WILL BE REMOVED ON INSERTING THE PROTEIN

	set k 0

	set maxul 0
	set maxll 0
	while { $k < [llength $data2] } {
		set yul [lindex $data2 [expr { $k + 2 }]]
		set zul [lindex $data2 [expr { $k + 3 }]]
		set yll [lindex $data2 [expr { $k + 4 }]]
		set zll [lindex $data2 [expr { $k + 5 }]]
		
		set ngyul [expr { $yul / $stepx }]
		set ngzul [expr { $zul / $stepy }]
		set ngul [expr { $ngyul * $ngzul }]

		set ngyll [expr { $yll / $stepx }]
		set ngzll [expr { $zll / $stepy }]
		set ngll [expr { $ngyll * $ngzll }]
		
		set numlul [format "%.0f" $ngul]
		set numlll [format "%.0f" $ngll]
		
		if { $numlul > $maxul } {
			set maxul $numlul
		}
		if { $numlll > $maxll } {
			set maxll $numlll
		}
		incr k 6
	}
	
	puts "			#### APPROXIMATELY $maxul NUMBER OF LIPIDS WILL BE REMOVED FROM LOWER LEAFLET####"
	puts "			#### APPROXIMATELY $maxll NUMBER OF LIPIDS WILL BE REMOVED FROM UPPER LEAFLET####"
	file delete output
}

proc spread {} {

	puts "				**** CALCULATING THE SPREAD OF EACH LIPID BASED ON THE INPUT ****"	

	set f [open "input" "r"]
	set data [read $f]
	close $f

	set k 1

	while { $k < [llength [lindex $data 0]] } {
		set filename [lindex $data 0 $k]
		set pivot [lindex $data 2 $k]

		set maxx 0.0
		set maxy 0.0
		set maxz 0.0

		set g [open "$filename" "r"]
		set data1 [read $g]
		close $g

		set k1 0

		while { [lindex $data1 $k1] != "$pivot" } {
			incr k1
		}

		set xori [lindex $data1 [expr { $k1 + 4 }]]
		set yori [lindex $data1 [expr { $k1 + 5 }]]
		set zori [lindex $data1 [expr { $k1 + 6 }]]

		set k1 0

		while { $k1 < [llength $data1] } {
			if { [lindex $data1 $k1] == "HETATM" && [lindex $data1 [expr { $k1 + 2 }]] != "$pivot" } {
				set x [lindex $data1 [expr { $k1 + 6 }]]
				set y [lindex $data1 [expr { $k1 + 7 }]]
				set z [lindex $data1 [expr { $k1 + 8 }]]

				set delx [expr { abs($xori - $x) }]
				set dely [expr { abs($yori - $y) }]
				set delz [expr { abs($zori - $z) }]

				set delx2 [expr { $delx * $delx }]
				set dely2 [expr { $dely * $dely }]
				set delz2 [expr { $delz * $delz }]
	
				if { $delx > $maxx } {
					set maxx $delx
				} 
	
				if { $dely > $maxy} {
					set maxy $dely
				}

				if { $delz > $maxz } {
					set maxz $delz
				}

			}
			incr k1
		}

		puts "		** EXTENSION OF [lindex $data 1 $k] ALONG X AXIS IS $maxx **"
		puts "		** EXTENSION OF [lindex $data 1 $k] ALONG Y AXIS IS $maxy **"
		puts "		** EXTENSION OF [lindex $data 1 $k] ALONG Z AXIS IS $maxz **"
		puts "		____________________________________________________________"
		puts ""
	
		incr k
	}
}

proc lipid_bilayer {} {

	puts ""
	puts "				**** FORMING A LIPID BILAYER ****"
	puts "				----------------------------------"
	puts ""
	
	# BUILDING LIPID BYLAYER

	package require math::linearalgebra

	set h [open "input" "r"]
	set data1 [read $h]
	close $h

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]
	set gd [lindex $data1 9 2]
	set method [lindex $data1 9 1]
	set version [lindex $data1 13 1]
	set lip_name [lindex $data1 13 2]

	# TRANSFORMATIONS

	# TRANSFORMATION 1 X AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
	set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Z TRANSLATION

	set z_trans [expr { $Z * -1 }]

	set ztrans [list 0 0 $z_trans]

	# TRANSFORMATION 4 Y AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm3r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	set g [open "lipids.pdb" "w"]

	set k1 0
	for {set k1 0} {$k1 < 100} {incr k1} {
		set l_ul($k1) 0
		set l_ll($k1) 0
	}

	set k1 1


	# X Y PLANE

	set num_lipids_ul 0.0
	set num_lipids_ll 0.0

	while { $k1 < [llength [lindex $data1 3]] } {
		set num_lipids_ul [expr { $num_lipids_ul + [lindex $data1 3 $k1] }]
		set num_lipids_ll [expr { $num_lipids_ll + [lindex $data1 4 $k1] }]
		set l_ul([expr { $k1 - 1 }]) [lindex $data1 3 $k1] 
		set l_ll([expr { $k1 - 1 }]) [lindex $data1 4 $k1]

		incr k1
	}

	if { $num_lipids_ul > $num_lipids_ll } {
		set num_lipids $num_lipids_ul
	} else {
		set num_lipids $num_lipids_ll
	}

	set div [expr { 1.0 / [expr { $k1 - 1 }] }]

	if { $method == 0 } {

		set pla [open "xy_grid" "w"]

		# RANDOM POINTS

		for {set i 0} {$i < $num_lipids} {incr i} {
			puts "				**** FINDING POINT [expr { $i + 1 }] of $num_lipids **** "
			set count 0
			set x_xy($i) [expr { rand() * $X }]
			set y_xy($i) [expr { rand() * $Y }]
			set x_xy_l($i) $x_xy($i)
			set y_xy_l($i) $y_xy($i)
			for {set j 0} {$j < $i} {incr j} {
				set diffx [expr { $x_xy($i) - $x_xy($j) }]
				set diffx2 [expr { $diffx * $diffx }]
				set diffy [expr { $y_xy($i) - $y_xy($j) }]
				set diffy2 [expr { $diffy * $diffy }] 
				set dis [expr { ($diffx2 + $diffy2) }]
				if { $dis < $gd } {
					incr count
				}
			}
			if { $count > 0 } {
				set i [expr { $i - 1 }]
			} else {
				puts $pla "$x_xy($i) $y_xy($i)"
			}
		}
		close $pla
		set n 0
	}


	if { $method == 1 } {
	
		# UNIFORM GRID + RANDOM LIPIDS UPPER LAYER

		set stepx [expr { $X / (sqrt($num_lipids))  }]
		set stepy [expr { $Y / (sqrt($num_lipids))  }]

		puts "			**** BUILDING THE UPPER LAYER ACCORDING TO METHOD1 ****"
		puts "			-------------------------------------------------------"
		puts ""

		set ij1 0
		set ij2 $l_ul(0)
		set ij3 [expr { $l_ul(0) + $l_ul(1) }]
		set ij4 [expr { $ij3 + $l_ul(2) }]
		set ij5 [expr { $ij4 + $l_ul(3) }]
		set ij6 [expr { $ij5 + $l_ul(4) }]
		set ij7 [expr { $ij6 + $l_ul(5) }]
		set ij8 [expr { $ij7 + $l_ul(6) }]
		for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
			#puts "			**** [format %.1f [expr { ($i * 100) / $X }]] % COMPLETE ****"	
			for {set j 0.0} {$j< $Y} { set j [expr { $j + $stepy }] } {		
				set pointer [expr { rand() }]
				if { $ij1 == $l_ul(0) && $ij2 == [expr { [expr { $l_ul(0) + $l_ul(1) }] }] && $ij3 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) }] && $ij4 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) }] && $ij5 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) }] && $ij6 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) }] && $ij7 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) + $l_ul(6)}]  && $ij8 == [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) + $l_ul(6) + $l_ul(7) }]} {
					set j [expr { $j + $stepy }] 
				} else { 
					set div1 [expr { $l_ul(0) / $num_lipids_ul }]
					if { $pointer <= $div1 } {
						if { $ij1 < $l_ul(0) } {
							set x_xy($ij1) $i 
							set y_xy($ij1) $j
							incr ij1
						} else { 
							set j [expr { $j - $stepy }]
						}
					}
				
					if { [expr { 2 * $div }] <= 1.0 } {
						set div2 [expr { $l_ul(1) / $num_lipids_ul }]
						if { $pointer > $div1 && $pointer <= [expr { $div1 + $div2 }]} {
							if { $ij2 < [expr { $l_ul(0) + $l_ul(1) }] } {
								set x_xy($ij2) $i
								set y_xy($ij2) $j
								incr ij2
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 3 * $div }] <= 1.0 } {
						set div3 [expr { $l_ul(2) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 }] && $pointer <= [expr { $div1 + $div2 + $div3 }] } {
							if { $ij3 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) }] } {
								set x_xy($ij3) $i
								set y_xy($ij3) $j
								incr ij3
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 4 * $div }] <= 1.0 } {
						set div4 [expr { $l_ul(3) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 + $div3 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 }] } {
							if { $ij4 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) }] } {
								set x_xy($ij4) $i
								set y_xy($ij4) $j
								incr ij4
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 5 * $div }] <= 1.0 } {
						set div5 [expr { $l_ul(4) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 }] } {
							if { $ij5 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) }] } {
								set x_xy($ij5) $i
								set y_xy($ij5) $j
								incr ij5
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 6 * $div }] <= 1.0 } {
						set div6 [expr { $l_ul(5) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5}] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 }] } {
							if { $ij6 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) }] } {
								set x_xy($ij6) $i
								set y_xy($ij6) $j
								incr ij6
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 7 * $div }] <= 1.0 } {
						set div7 [expr { $l_ul(6) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6}] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 }] } {
							if { $ij7 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) + $l_ul(6) }] } {
								set x_xy($ij7) $i
								set y_xy($ij7) $j
								incr ij7
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 8 * $div }] <= 1.0 } {
						set div8 [expr { $l_ul(7) / $num_lipids_ul }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 + $div8 }] } {
							if { $ij8 < [expr { $l_ul(0) + $l_ul(1) + $l_ul(2) + $l_ul(3) + $l_ul(4) + $l_ul(5) + $l_ul(6) + $l_ul(7)}] } {
								set x_xy($ij8) $i
								set y_xy($ij8) $j
								incr ij8
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}
				}
			}
		}

		# UNIFORM GRID + RANDOM LIPIDS LOWER LAYER

		puts "			**** BUILDING THE LOWER LAYER ACCORDING TO METHOD1 ****"
		puts "			-------------------------------------------------------"
		puts ""

		set ij1 0
		set ij2 $l_ll(0)
		set ij3 [expr { $l_ll(0) + $l_ll(1) }]
		set ij4 [expr { $ij3 + $l_ll(2) }]
		set ij5 [expr { $ij4 + $l_ll(3) }]
		set ij6 [expr { $ij5 + $l_ll(4) }]
		set ij7 [expr { $ij6 + $l_ll(5) }]
		set ij8 [expr { $ij7 + $l_ll(6) }]
		for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
			#puts "			**** [format %.1f [expr { ($i * 100) / $X }]] % COMPLETE ****"				
			for {set j 0.0} {$j< $Y} {set j [expr { $j + $stepy }] } {
				set pointer [expr { rand() }]
				if { $ij1 == $l_ll(0) && $ij2 == [expr { [expr { $l_ll(0) + $l_ll(1) }] }] && $ij3 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) }] && $ij4 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) }] && $ij5 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) }] && $ij6 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) }] && $ij7 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) + $l_ll(6)}]  && $ij8 == [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) + $l_ll(6) + $l_ll(7) }]} {
					set j [expr { $j + $stepy }] 
				} else {
					set div1 [expr { $l_ll(0) / $num_lipids_ll }]
					if { $pointer <= $div1 } {
						if { $ij1 < $l_ll(0) } {
							set x_xy_l($ij1) $i 
							set y_xy_l($ij1) $j
							incr ij1
						} else { 
							set j [expr { $j - $stepy }]
						}
					}

					if { [expr { 2 * $div }] <= 1.0 } {
						set div2 [expr { $l_ll(1) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 }] && $pointer <= [expr { $div1 + $div2 }] } {
							if { $ij2 < [expr { $l_ll(1) + $l_ll(0) }] } {
								set x_xy_l($ij2) $i
								set y_xy_l($ij2) $j
								incr ij2
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 3 * $div }] <= 1.0 } {
						set div3 [expr { $l_ll(2) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 }] && $pointer <= [expr { $div1 + $div2 + $div3 }] } {
							if { $ij3 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) }] } {
								set x_xy_l($ij3) $i
								set y_xy_l($ij3) $j
								incr ij3
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 4 * $div }] <= 1.0 } {
						set div4 [expr { $l_ll(3) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 + $div3 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 }] } {
							if { $ij4 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) }] } {
								set x_xy_l($ij4) $i
								set y_xy_l($ij4) $j
								incr ij4
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 5 * $div }] <= 1.0 } {
						set div5 [expr { $l_ll(4) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 }] } {
							if { $ij5 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) }] } {
								set x_xy_l($ij5) $i
								set y_xy_l($ij5) $j
								incr ij5
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 6 * $div }] <= 1.0 } {
						set div6 [expr { $l_ll(5) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 }] } {
							if { $ij6 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) }] } {
								set x_xy_l($ij6) $i
								set y_xy_l($ij6) $j
								incr ij6
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 7 * $div }] <= 1.0 } {
						set div7 [expr { $l_ll(6) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 }] } {
							if { $ij7 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) + $l_ll(6) }] } {
								set x_xy_l($ij7) $i
								set y_xy_l($ij7) $j
								incr ij7
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}

					if { [expr { 8 * $div }] <= 1.0 } {
						set div8 [expr { $l_ll(7) / $num_lipids_ll }]
						if { $pointer > [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 }] && $pointer <= [expr { $div1 + $div2 + $div3 + $div4 + $div5 + $div6 + $div7 + $div8 }] } {
							if { $ij8 < [expr { $l_ll(0) + $l_ll(1) + $l_ll(2) + $l_ll(3) + $l_ll(4) + $l_ll(5) + $l_ll(6) + $l_ll(7)}] } {
								set x_xy_l($ij8) $i
								set y_xy_l($ij8) $j
								incr ij8
							} else { 
								set j [expr { $j - $stepy }]
							}
						}
					}
				}
			}
		}
	}

	set num_residues [format "%0.0f" $num_lipids]
	set m 0
	if { $method == 1 } {
		set n 0
	} 

	set k1 1

	puts "				****PUTTING THE LIPIDS INSIDE ****"
	puts "				-----------------------------------"
	puts ""

	while { $k1 < [llength [lindex $data1 0]] } {
		set pdb [lindex $data1 0 $k1]
		set resname [lindex $data1 1 $k1]
		set pivot [lindex $data1 2 $k1]
		set num_lu [lindex $data1 3 $k1]
		set num_ll [lindex $data1 4 $k1]

		puts "		**** PDB :: $pdb ****"

		set f [open "$pdb" "r"]
		set data [read $f]
		close $f

		set k 0
		 
		while { $k < [llength $data] } {
			if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
				set xori [lindex $data [expr { $k + 6 }]]
				set yori [lindex $data [expr { $k + 7 }]]
				set zori [lindex $data [expr { $k + 8 }]]
			}
			incr k
		}

		# BUILDING THE NEW PDB WITH ORIGIN AT PIVOT ELEMENT

		# space variables

		set p(1) "   "
		set p(2) "  "
		set p(3) " "
		set p(4) ""
	
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

		set sat(1) "   "
		set sat(2) "  "
		set sat(3) " "
		set sat(4) " "

		set ic(1) " "
		set ic(2) " "
		set ic(3) " "
		set ic(4) ""

		for {set i 0} {$i < $num_lu} {incr i} {

			set x_shift $x_xy($m)
			set y_shift $y_xy($m)
			incr m
			#set z_shift $zori
			set z_shift 0.0

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					#set an1 [lindex $data [expr { $k + 5 }]]
					set an1 $m
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 6 }]]
					set x1 [format "%.3f" [expr { $x1 - $xori + $x_shift }]]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 7 }]] 
					set y1 [format "%.3f" [expr { $y1 - $yori + $y_shift}]]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 8 }]] 
					set z1 [format "%.3f" [expr { $z1 - $zori + $z_shift}]]
					set sz1 [string length $z1]

					#set sresname [string length [lindex $data [expr { $k + 3 }]]]
					set sresname [llength $resname]

					if { $zori < 0.0 } {

						set coord [list [expr { $x1 - $x_shift }] [expr { $y1 - $y_shift }] [expr { $z1 - $z_shift }]]

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						set dum [open "dummy" "w"]
	
						puts $dum "$tcoord"

						close $dum

						set dum [open "dummy" "r"]
						set dvar [read $dum]
						close $dum

						set x1 [format "%.3f" [lindex $dvar 0]]
						set x1 [format "%.3f" [expr { $x1 + $x_shift}]]
						set sx1 [string length $x1]

						set y1 [format "%.3f" [lindex $dvar 1]]
						set y1 [format "%.3f" [expr { $y1 + $y_shift}]]
						set sy1 [string length $y1]

						set z1 [format "%.3f" [lindex $dvar 2]]
						set z1 [format "%.3f" [expr { $z1 + $z_shift}]]
						set sz1 [string length $z1]
					}
					if { $method == 1 } {
						if { $at1 == $pivot } {
							puts $g "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 11 }]]"
						}
					} else { 
						puts $g "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 11 }]]"
					}
				}
			incr k
			}
		puts $g "TER"
		}

		# BUILDING ITS IMAGE ON THE OTHER LAYER

		for {set i 0} {$i < $num_ll} {incr i} {

			set x_shift $x_xy_l($n)
			set y_shift $y_xy_l($n)

			incr n
			
			#set z_shift $zori
			set z_shift 0.0

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					#set an1 [expr { 1 + [lindex $data [expr { $k + 5 }]] }]
		
					set an1 [expr { $num_residues + $n }]
				
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 6 }]]
					set x1 [format "%.3f" [expr { $x1 - $xori }]]
	
					set y1 [lindex $data [expr { $k + 7 }]] 
					set y1 [format "%.3f" [expr { $y1 - $yori }]]

					set z1 [lindex $data [expr { $k + 8 }]] 
					set z1 [format "%.3f" [expr { $z1 - $zori }]]

					#set sresname [string length [lindex $data [expr { $k + 3 }]]]
					set sresname [llength $resname]

					set coord [list $x1 $y1 $z1]

					if { $zori > 0.0 } {

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						# TRANSFORMATION 3

						set tcoord [::math::linearalgebra::add $tcoord $ztrans]
					} else { 
						set tcoord [::math::linearalgebra::add $coord $ztrans]
					}

					set dum [open "dummy" "w"]
	
					puts $dum "$tcoord"

					close $dum

					set dum [open "dummy" "r"]
					set dvar [read $dum]
					close $dum

					set x1 [format "%.3f" [lindex $dvar 0]]
					set x1 [format "%.3f" [expr { $x1 + $x_shift}]]
					set sx1 [string length $x1]

					set y1 [format "%.3f" [lindex $dvar 1]]
					set y1 [format "%.3f" [expr { $y1 + $y_shift}]]
					set sy1 [string length $y1]

					set z1 [format "%.3f" [lindex $dvar 2]]
					set z1 [format "%.3f" [expr { $z1 + $z_shift}]]
					set sz1 [string length $z1]
					if { $method == 1 } {
						if { $at1 == $pivot } {
							puts $g "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 11 }]]"
						}
					} else { 
						puts $g "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 11 }]]"
					}
				}
			incr k 
			}
			puts $g "TER"
		}
		incr k1
	}
	puts $g "END"
	close $g

	set h [open "input" "r"]
	set data1 [read $h]
	close $h

	set version [lindex $data1 13 1]
	set num_dl [lindex $data1 13 [expr { [llength [lindex $data1 13]] - 1 } ]]

	set k 2

	if { $method == 1 } {
		if { $version == 3.0 } {
			puts "				#### USING VERSION 3.0 ####"
			if { $num_dl > 1 } {
				grid_mapping_assy
			} else {
				grid_mapping "[lindex $data1 13 $k]"
			}		
		} else {
			puts "				#### USING VERSION 2.0 ####" 
			lipid_growth
		}
	} else { 
		lipid_overlap
	}
}

proc lipid_growth {} {	

	puts "				**** REMOVING THE OVERLAPS IN BOTH LAYERS ACCORDING TO LIPID GROWTH ALGORITHM ****"

	package require math::linearalgebra

	set f [open "lipids.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "input" "r"]
	set data1 [read $g]
	close $g

	set h [open "lipid_order.pdb" "w"]

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]

	# TRANSFORMATIONS

	# TRANSFORMATION 1 X AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
	set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Z TRANSLATION

	set z_trans [expr { $Z * -1 }]

	set ztrans [list 0 0 $z_trans]

	# TRANSFORMATION 4 Y AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm3r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set k1 1
	set num_lipids_ul 0.0
	set num_lipids_ll 0.0

	while { $k1 < [llength [lindex $data1 3]] } {
		set num_lipids_ul [expr { $num_lipids_ul + [lindex $data1 3 $k1] }]
		set num_lipids_ll [expr { $num_lipids_ll + [lindex $data1 4 $k1] }]
		if { [lindex $data1 3 $k1] > [lindex $data1 4 $k1] } {
			set l([expr { $k1 - 1 }]) [lindex $data1 3 $k1]
		} else { 
			set l([expr { $k1 - 1 }]) [lindex $data1 4 $k1]
		}
		incr k1
	}

	if { $num_lipids_ul > $num_lipids_ll } {
		set num_lipids $num_lipids_ul
	} else {
		set num_lipids $num_lipids_ll
	}

	set stepx [expr { $X / (sqrt($num_lipids))  }]
	set stepy [expr { $Y / (sqrt($num_lipids))  }]

	# UPPER LAYER
	set z 0.000

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j< $Y} {set j [expr { $j + $stepy }] } {
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
	
					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "ter"
					}
				}
				incr k
			}
		}
	}

	# LOWER LAYER

	set z [expr { -1 * $Z }]
	set z [format "%.3f" $z]

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j < $Y} { set j [expr { $j + $stepy }] } {	
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]

					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "TER"
					}
				}
				incr k
			}
		}
	}
	puts $h "END"
	close $h


	# BUILDING THE WHOLE LIPID BILAYER WITH NO OVERLAP

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n
	global rn1n
	global at1n
	global resnamen
	global atomnamen
	global xshift
	global yshift
	global zshift

	set h [open "lipid_order.pdb" "r"]
	set data2 [read $h]
	close $h

	set g1 [open "lipids_no.pdb" "w"]

	# REMOVING THE OVERLAPS IN UPPER LAYER

	puts "				**** REMOVING THE OVERLAPS IN UPPER LAYER **** "

	# DETERMING THE NEIGBOURS OF EACH RESIDUE

	set nei 1

	set h1 [open "neighbours" "w"]

	set k 0

	while { $k < [llength $data2] } {
		if { [lindex $data2 $k] == "HETATM" } {
			set x [lindex $data2 [expr { $k + 5 }]]
			set y [lindex $data2 [expr { $k + 6 }]]
			set z [lindex $data2 [expr { $k + 7 }]]
			set rid [lindex $data2 [expr { $k + 4 }]]
			puts $h1 "\{ $rid"
			for {set i [expr { -1*$nei }]} {$i <= $nei} {incr i} {
				set icord [format "%.3f" [expr { $x - ($i*$stepx) }]]
				for {set j [expr { -1*$nei } ]} {$j <= $nei } {incr j} {
					set jcord [format "%.3f" [expr { $y - ($j*$stepy) }]]
					set k1 0
					while {$k1 < [llength $data2]} {
						if { [lindex $data2 $k1] == "HETATM" } {
							set x1 [lindex $data2 [expr { $k1 + 5 }]]
							set delx1 [expr { abs($x-$x1) }]
			
							set y1 [lindex $data2 [expr { $k1 + 6 }]]
							set dely1 [expr { abs($y-$y1) }]
							
							set z1 [lindex $data2 [expr { $k1 + 7 }]]
							set rid [lindex $data2 [expr { $k1 + 4 }]]
							if { $x1 != $x || $y1 != $y } {
								if { [format "%.0f" $x1]  == [format "%.0f" $icord] && [format "%.0f" $y1] == [format "%.0f" $jcord] && $z1 == $z } {
									puts $h1 "$rid"
								}
							}
						}
						incr k1
					}
				}
			}
			puts $h1 "\}"
		}
		incr k
	}

	close $h1

	set z 0.000
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)

					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					if { $zori < 0.0 } {

						set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						set dum [open "dummy" "w"]

						puts $dum "$tcoord"

						close $dum

						set dum [open "dummy" "r"]
						set dvar [read $dum]
						close $dum

						set x1n($n) [format "%.3f" [lindex $dvar 0]]
						set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
						set sx1 [string length $x1n($n)]

						set y1n($n) [format "%.3f" [lindex $dvar 1]]
						set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
						set sy1 [string length $y1n($n)]

						set z1n($n) [format "%.3f" [lindex $dvar 2]]
						set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
						set sz1 [string length $z1n($n)]
					}
				incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1

	set num_atom $n
	set cum_n $n
	set an1n(-1) 0
	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				set otx 0.0
				set oty 0.0
				set otz -3.14
				for {set tx 0.0} {$tx < [expr { $stepx / 2.0 }] } {set tx [expr { $tx + 0.5 }] } {
					for {set ty 0.0} {$ty < [expr { $stepy / 2.0 }] } {set ty [expr { $ty + 0.5 }] } {
						for {set tz -3.14} {$tz <= 3.14} { set tz [expr { $tz + [format "%.2f" [expr { 3.14 / 30.0 }]] }] } {
							#puts "				Tx=$tx	Ty=$ty	Tz=$tz"
							transformations 3 [expr { $tx - $otx }] $an1n($n)
							transformations 4 [expr { $ty - $oty }] $an1n($n)
							transformations 2 [expr { $tz - $otz }] $an1n($n)
							set otz $tz
							set oty $ty
							set otx $tx
							set check [overlap $an1n($n)] 
							if { $check == 0 } {
								puts "				FINAL VALUE :: Tx=$tx	Ty=$ty	Tz=$tz"
								set tx [expr {  $stepx / 2.0  }]
								set ty [expr {  $stepy / 2.0  }]
								#set tz [expr { (10.0 * 3.14) /180.0 }]
								set tz 4.00
							}
						}
						set otz -3.14
					}
				}
			}
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}
		}
	}

	# REMOVING THE OVERLAPS IN LOWER LAYER
	puts ""
	puts "				**** REMOVING THE OVERLAPS IN LOWER LAYER **** "

	set z [expr { -1 * $Z }]
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)
					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

					if { $zori > 0.0 } {

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						# TRANSFORMATION 3

						set tcoord [::math::linearalgebra::add $tcoord $ztrans]
					} else { 
						set tcoord [::math::linearalgebra::add $coord $ztrans]
					}

					set dum [open "dummy" "w"]

					puts $dum "$tcoord"

					close $dum

					set dum [open "dummy" "r"]
					set dvar [read $dum]
					close $dum

					set x1n($n) [format "%.3f" [lindex $dvar 0]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
					set sx1 [string length $x1n($n)]

					set y1n($n) [format "%.3f" [lindex $dvar 1]]
					set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
					set sy1 [string length $y1n($n)]

					set z1n($n) [format "%.3f" [lindex $dvar 2]]
					set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
					set sz1 [string length $z1n($n)]

					incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1
	set num_atom $n
	set cum_n $n
	set an1n(-1) 0
	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				set otx 0.0
				set oty 0.0
				set otz -3.14
				for {set tx 0.0} {$tx < [expr { $stepx / 2.0 }] } {set tx [expr { $tx + 0.5 }] } {
					for {set ty 0.0} {$ty < [expr { $stepy / 2.0 }] } {set ty [expr { $ty + 0.5 }] } {
						for {set tz -3.14} {$tz <= 3.14} { set tz [expr { $tz + [format "%.2f" [expr { 3.14 / 30.0 }]] }] } {
							#puts "				Tx=$tx	Ty=$ty	Tz=$tz"
							transformations 3 [expr { $tx - $otx }] $an1n($n)
							transformations 4 [expr { $ty - $oty }] $an1n($n)
							transformations 2 [expr { $tz - $otz }] $an1n($n)
							set otz $tz
							set oty $ty
							set otx $tx
							set check [overlap $an1n($n)] 
							if { $check == 0 } {
								puts "				FINAL VALUE :: Tx=$tx	Ty=$ty	Tz=$tz"
								set tx [expr {  $stepx / 2.0  }]
								set ty [expr {  $stepy / 2.0  }]
								#set tz [expr { (10.0 * 3.14) /180.0 }]
								set tz 4.0
							}
						}
						set otz -3.14
					}
				}
			}
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}
		}
	}

	#puts $g1 "END"	
	close $g1
}

proc grid_mapping { inpres } {	

	puts "				**** REMOVING THE OVERLAPS IN BOTH LAYERS ACCORDING TO LIPID GROWTH ALGORITHM ****"

	package require math::linearalgebra

	set f [open "lipids.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "input" "r"]
	set data1 [read $g]
	close $g

	set h [open "lipid_order.pdb" "w"]

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]

	# TRANSFORMATIONS

	# TRANSFORMATION 1 X AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
	set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Z TRANSLATION

	set z_trans [expr { $Z * -1 }]

	set ztrans [list 0 0 $z_trans]

	# TRANSFORMATION 4 Y AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm3r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set k1 1
	set num_lipids_ul 0.0
	set num_lipids_ll 0.0

	while { $k1 < [llength [lindex $data1 3]] } {
		set num_lipids_ul [expr { $num_lipids_ul + [lindex $data1 3 $k1] }]
		set num_lipids_ll [expr { $num_lipids_ll + [lindex $data1 4 $k1] }]
		if { [lindex $data1 3 $k1] > [lindex $data1 4 $k1] } {
			set l([expr { $k1 - 1 }]) [lindex $data1 3 $k1]
		} else { 
			set l([expr { $k1 - 1 }]) [lindex $data1 4 $k1]
		}
		incr k1
	}

	if { $num_lipids_ul > $num_lipids_ll } {
		set num_lipids $num_lipids_ul
	} else {
		set num_lipids $num_lipids_ll
	}

	set stepx [expr { $X / (sqrt($num_lipids))  }]
	set stepy [expr { $Y / (sqrt($num_lipids))  }]

	# UPPER LAYER
	set z 0.000

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j< $Y} {set j [expr { $j + $stepy }] } {
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
	
					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "ter"
					}
				}
				incr k
			}
		}
	}

	# LOWER LAYER

	set z [expr { -1 * $Z }]
	set z [format "%.3f" $z]

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j < $Y} { set j [expr { $j + $stepy }] } {	
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]

					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "TER"
					}
				}
				incr k
			}
		}
	}
	puts $h "END"
	close $h


	# BUILDING THE WHOLE LIPID BILAYER WITH NO OVERLAP

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n
	global rn1n
	global at1n
	global resnamen
	global atomnamen
	global xshift
	global yshift
	global zshift

	set h [open "lipid_order.pdb" "r"]
	set data2 [read $h]
	close $h

	set g1 [open "lipids_no.pdb" "w"]

	# REMOVING THE OVERLAPS IN UPPER LAYER

	puts "				**** REMOVING THE OVERLAPS IN UPPER LAYER **** "

	# DETERMING THE NEIGBOURS OF EACH RESIDUE

	set nei 1

	set h1 [open "neighbours" "w"]

	set k 0

	while { $k < [llength $data2] } {
		if { [lindex $data2 $k] == "HETATM" } {
			set x [lindex $data2 [expr { $k + 5 }]]
			set y [lindex $data2 [expr { $k + 6 }]]
			set z [lindex $data2 [expr { $k + 7 }]]
			set rid [lindex $data2 [expr { $k + 4 }]]
			puts $h1 "\{ $rid"
			for {set i [expr { -1*$nei }]} {$i <= $nei} {incr i} {
				set icord [format "%.3f" [expr { $x - ($i*$stepx) }]]
				for {set j [expr { -1*$nei } ]} {$j <= $nei } {incr j} {
					set jcord [format "%.3f" [expr { $y - ($j*$stepy) }]]
					set k1 0
					while {$k1 < [llength $data2]} {
						if { [lindex $data2 $k1] == "HETATM" } {
							set x1 [lindex $data2 [expr { $k1 + 5 }]]
							set delx1 [expr { abs($x-$x1) }]
			
							set y1 [lindex $data2 [expr { $k1 + 6 }]]
							set dely1 [expr { abs($y-$y1) }]
							
							set z1 [lindex $data2 [expr { $k1 + 7 }]]
							set rid [lindex $data2 [expr { $k1 + 4 }]]
							if { $x1 != $x || $y1 != $y } {
								if { [format "%.0f" $x1]  == [format "%.0f" $icord] && [format "%.0f" $y1] == [format "%.0f" $jcord] && $z1 == $z } {
									puts $h1 "$rid"
								}
							}
						}
						incr k1
					}
				}
			}
			puts $h1 "\}"
		}
		incr k
	}

	close $h1

	set z 0.000
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)

					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					if { $zori < 0.0 } {

						set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						set dum [open "dummy" "w"]

						puts $dum "$tcoord"

						close $dum

						set dum [open "dummy" "r"]
						set dvar [read $dum]
						close $dum

						set x1n($n) [format "%.3f" [lindex $dvar 0]]
						set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
						set sx1 [string length $x1n($n)]

						set y1n($n) [format "%.3f" [lindex $dvar 1]]
						set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
						set sy1 [string length $y1n($n)]

						set z1n($n) [format "%.3f" [lindex $dvar 2]]
						set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
						set sz1 [string length $z1n($n)]
					}
				incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1

	set num_atom $n
	set cum_n $n
	set an1n(-1) 0

	for {set n 0} {$n < 10} {incr n} {
		set tra($n) 0.0
	}
	set g10 [open "grid" "r"]
	set gd [read $g10]
	close $g10

	set k10 0

	while { [lindex $gd $k10] != $inpres } {
		incr k10
	}  
	set gp $k10
	incr k10

	set k11 1

	set frs [expr { $X / $stepx }]
	set frs [format "%.0f" $frs]
	set fcs [expr { $Y / $stepy }]
	set fcs [format "%.0f" $fcs]

	set resad 0
	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"

			set ulgrid [lindex $gd $k10]

			set gtype [lindex $ulgrid 0]

			set k11 [expr { ($resad / $frs) + 1 }]

			if { $gtype == "A" } {
				if { $k11 == 1} {
					set numft [llength [lindex $ulgrid $k11]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid $k11 $k12]
						incr k12
					}
				}
				if { $k11 >= 2 && $k11 < $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
			} elseif { $gtype == "S" } {
				if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs } {
					set numft [llength [lindex $ulgrid 1]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 1 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
			}	elseif { $gtype == "S1" } {
				if { [expr { $k11 % 2 }] != 0 } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 3 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 1 $k12]
							incr k12
						}
					}
				}
				if { [expr { $k11 % 2 }] == 0 } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 4 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 2]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 $k12]
							incr k12
						}
					}
				}
			} elseif { $gtype == "S2" } {
				if { $k11 == 1 } {
					set numft [llength [lindex $ulgrid 1]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 1 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs && $k11 != 1} {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 4]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 4 $k12]
						incr k12
					}
				}
			} elseif { $gtype == "A1" } {
				if { $k11 == 1} {
					if { $resad == 0 } {
						set numft [llength [lindex $ulgrid $k11 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid $k11 0 $k12]
							incr k12
						}
					} else { 
						set numft [llength [lindex $ulgrid $k11 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid $k11 1 $k12]
							incr k12
						}
					}
				}
				if { $k11 >= 2 && $k11 < $fcs} {
					if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
						set numft [llength [lindex $ulgrid 2 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 0 $k12]
							incr k12
						} 
						} else {
						set numft [llength [lindex $ulgrid 2 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 1 $k12]
							incr k12
						}
					}
				}
				if { $k11 == $fcs } {
					if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
						set numft [llength [lindex $ulgrid 3 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 3 0 $k12]
							incr k12
						} 
						} else {
						set numft [llength [lindex $ulgrid 3 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 3 1 $k12]
							incr k12
						}
					}
				}
			} elseif { $gtype == "S3" } {
				if { [expr { $k11 % 2 }] != 0 && $k11 != $frs } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 3 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 1 $k12]
							incr k12
						}
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $frs} {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 4 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 2]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 $k12]
							incr k12
						}
					}
				}
				if { $k11 == $frs } {
					set numft [llength [lindex $ulgrid 5]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 5 $k12]
						incr k12
					}
				}	
			}
				
			if { $k12 >= $numft } {
				set k12 0
			} else {
				incr k12
			} 

			set tx [lindex $tra($k12) 0]
			set ty [lindex $tra($k12) 1]
			set tz [lindex $tra($k12) 2]

			if { $ty != 0 } {
				for {set tyn 0.0} {$tyn < $ty} {set tyn [expr { $tyn + 0.5}]} {
					transformations 3 0.0 $an1n($n)
					transformations 4 $tyn $an1n($n)
					transformations 2 6.2 $an1n($n)
				}
			}
				
			transformations 3 $tx $an1n($n)
			transformations 4 $ty $an1n($n)
			transformations 2 $tz $an1n($n)
						
			#set check [overlap $an1n($n)]
			set check 0
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}
			incr resad
		}
	}

	# REMOVING THE OVERLAPS IN LOWER LAYER
	puts ""
	puts "				**** REMOVING THE OVERLAPS IN LOWER LAYER **** "

	set z [expr { -1 * $Z }]
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)
					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

					if { $zori > 0.0 } {

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						# TRANSFORMATION 3

						set tcoord [::math::linearalgebra::add $tcoord $ztrans]
					} else { 
						set tcoord [::math::linearalgebra::add $coord $ztrans]
					}

					set dum [open "dummy" "w"]

					puts $dum "$tcoord"

					close $dum

					set dum [open "dummy" "r"]
					set dvar [read $dum]
					close $dum

					set x1n($n) [format "%.3f" [lindex $dvar 0]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
					set sx1 [string length $x1n($n)]

					set y1n($n) [format "%.3f" [lindex $dvar 1]]
					set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
					set sy1 [string length $y1n($n)]

					set z1n($n) [format "%.3f" [lindex $dvar 2]]
					set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
					set sz1 [string length $z1n($n)]

					incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1
	set num_atom $n
	set cum_n $n
	set an1n(-1) 0

	set resad 0
	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"

			set ulgrid [lindex $gd [expr { $k10 + 1 }]]

			set gtype [lindex $ulgrid 0]

			set k11 [expr { ($resad / $frs) + 1 }]

			if { $gtype == "A" } {
				if { $k11 == 1} {
					set numft [llength [lindex $ulgrid $k11]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid $k11 $k12]
						incr k12
					}
				}
				if { $k11 >= 2 && $k11 < $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
			} elseif { $gtype == "S" } {
				if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs } {
					set numft [llength [lindex $ulgrid 1]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 1 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
			}	elseif { $gtype == "S1" } {
				if { [expr { $k11 % 2 }] != 0 } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 3 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 1 $k12]
							incr k12
						}
					}
				}
				if { [expr { $k11 % 2 }] == 0 } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 4 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 2]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 $k12]
							incr k12
						}
					}
				}
			} elseif { $gtype == "S2" } {
				if { $k11 == 1 } {
					set numft [llength [lindex $ulgrid 1]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 1 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs && $k11 != 1} {
					set numft [llength [lindex $ulgrid 3]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 3 $k12]
						incr k12
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
					set numft [llength [lindex $ulgrid 2]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 2 $k12]
						incr k12
					}
				}
				if { $k11 == $fcs } {
					set numft [llength [lindex $ulgrid 4]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 4 $k12]
						incr k12
					}
				}
			}  elseif { $gtype == "A1" } {
				if { $k11 == 1} {
					if { $resad == 0 } {
						set numft [llength [lindex $ulgrid $k11 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid $k11 0 $k12]
							incr k12
						}
					} else { 
						set numft [llength [lindex $ulgrid $k11 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid $k11 1 $k12]
							incr k12
						}
					}
				}
				if { $k11 >= 2 && $k11 < $fcs} {
					if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
						set numft [llength [lindex $ulgrid 2 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 0 $k12]
							incr k12
						} 
						} else {
						set numft [llength [lindex $ulgrid 2 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 1 $k12]
							incr k12
						}
					}
				}
				if { $k11 == $fcs } {
					if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
						set numft [llength [lindex $ulgrid 3 0]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 3 0 $k12]
							incr k12
						} 
						} else {
						set numft [llength [lindex $ulgrid 3 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 3 1 $k12]
							incr k12
						}
					}
				}
			} elseif { $gtype == "S3" } {
				if { [expr { $k11 % 2 }] != 0 && $k11 != $frs } {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 3 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 1]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 1 $k12]
							incr k12
						}
					}
				}
				if { [expr { $k11 % 2 }] == 0 && $k11 != $frs} {
					if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
						set numft 1
						set k12 0
						set tra($k12) [lindex $ulgrid 4 0]
						incr k12
					} else {
						set numft [llength [lindex $ulgrid 2]]
						set k12 0
						while { $k12 < $numft } {
							set tra($k12) [lindex $ulgrid 2 $k12]
							incr k12
						}
					}
				}
				if { $k11 == $frs } {
					set numft [llength [lindex $ulgrid 5]]
					set k12 0
					while { $k12 < $numft } {
						set tra($k12) [lindex $ulgrid 5 $k12]
						incr k12
					}
				}	
			}
			
			if { $k12 >= $numft } {
				set k12 0
			} else {
				incr k12
			} 

			set tx [lindex $tra($k12) 0]
			set ty [lindex $tra($k12) 1]
			set tz [lindex $tra($k12) 2]

			if { $ty != 0 } {
				for {set tyn 0.0} {$tyn < $ty} {set tyn [expr { $tyn + 0.5}]} {
					transformations 3 0.0 $an1n($n)
					transformations 4 $tyn $an1n($n)
					transformations 2 6.2 $an1n($n)
				}
			}
		
			transformations 3 $tx $an1n($n)
			transformations 4 $ty $an1n($n)
			transformations 2 $tz $an1n($n)
			
			#set check [overlap $an1n($n)]
			set check 0

			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}	
		incr resad
		}
	}

	#puts $g1 "END"	
	close $g1
}

proc grid_mapping_assy { } {	

	puts "				**** REMOVING THE OVERLAPS IN BOTH LAYERS ACCORDING TO LIPID GROWTH ALGORITHM ****"

	package require math::linearalgebra

	set f [open "lipids.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "input" "r"]
	set data1 [read $g]
	close $g

	set h [open "lipid_order.pdb" "w"]

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]

	set ndlip [lindex $data1 13 [expr { [llength [lindex $data1 13]] - 1 }]]

	# TRANSFORMATIONS

	# TRANSFORMATION 1 X AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
	set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Z TRANSLATION

	set z_trans [expr { $Z * -1 }]

	set ztrans [list 0 0 $z_trans]

	# TRANSFORMATION 4 Y AXIS

	set theta [expr { (1.0 * 3.14159) / 1.0 }]

	set tm3r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set k1 1
	set num_lipids_ul 0.0
	set num_lipids_ll 0.0

	while { $k1 < [llength [lindex $data1 3]] } {
		set num_lipids_ul [expr { $num_lipids_ul + [lindex $data1 3 $k1] }]
		set num_lipids_ll [expr { $num_lipids_ll + [lindex $data1 4 $k1] }]
		if { [lindex $data1 3 $k1] > [lindex $data1 4 $k1] } {
			set l([expr { $k1 - 1 }]) [lindex $data1 3 $k1]
		} else { 
			set l([expr { $k1 - 1 }]) [lindex $data1 4 $k1]
		}
		incr k1
	}

	if { $num_lipids_ul > $num_lipids_ll } {
		set num_lipids $num_lipids_ul
	} else {
		set num_lipids $num_lipids_ll
	}

	set stepx [expr { $X / (sqrt($num_lipids))  }]
	set stepy [expr { $Y / (sqrt($num_lipids))  }]

	# UPPER LAYER
	set z 0.000

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j< $Y} {set j [expr { $j + $stepy }] } {
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
	
					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "ter"
					}
				}
				incr k
			}
		}
	}

	# LOWER LAYER

	set z [expr { -1 * $Z }]
	set z [format "%.3f" $z]

	for {set i 0.0} {$i < $X} {set i [expr { $i + $stepx  }] } {
		for {set j 0.0} {$j < $Y} { set j [expr { $j + $stepy }] } {	
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]

					set at1 [lindex $data [expr { $k + 2 }]]
					set sat1 [string length $at1]

					set an1 [lindex $data [expr { $k + 4 }]]
					set san1 [string length $an1]

					set x1 [lindex $data [expr { $k + 5 }]]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

					set y1 [lindex $data [expr { $k + 6 }]] 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]

					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]

					set sresname [string length [lindex $data [expr { $k + 3 }]]]

					if { $x1 == [format "%.3f" $i] && $y1 == [format "%.3f" $j] && $z1 == $z} {
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
						puts $h "TER"
					}
				}
				incr k
			}
		}
	}
	puts $h "END"
	close $h


	# BUILDING THE WHOLE LIPID BILAYER WITH NO OVERLAP

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n
	global rn1n
	global at1n
	global resnamen
	global atomnamen
	global xshift
	global yshift
	global zshift

	set h [open "lipid_order.pdb" "r"]
	set data2 [read $h]
	close $h

	set g1 [open "lipids_no.pdb" "w"]

	# REMOVING THE OVERLAPS IN UPPER LAYER

	puts "				**** REMOVING THE OVERLAPS IN UPPER LAYER **** "

	# DETERMING THE NEIGBOURS OF EACH RESIDUE

	set nei 1

	set h1 [open "neighbours" "w"]

	set k 0

	while { $k < [llength $data2] } {
		if { [lindex $data2 $k] == "HETATM" } {
			set x [lindex $data2 [expr { $k + 5 }]]
			set y [lindex $data2 [expr { $k + 6 }]]
			set z [lindex $data2 [expr { $k + 7 }]]
			set rid [lindex $data2 [expr { $k + 4 }]]
			puts $h1 "\{ $rid"
			for {set i [expr { -1*$nei }]} {$i <= $nei} {incr i} {
				set icord [format "%.3f" [expr { $x - ($i*$stepx) }]]
				for {set j [expr { -1*$nei } ]} {$j <= $nei } {incr j} {
					set jcord [format "%.3f" [expr { $y - ($j*$stepy) }]]
					set k1 0
					while {$k1 < [llength $data2]} {
						if { [lindex $data2 $k1] == "HETATM" } {
							set x1 [lindex $data2 [expr { $k1 + 5 }]]
							set delx1 [expr { abs($x-$x1) }]
			
							set y1 [lindex $data2 [expr { $k1 + 6 }]]
							set dely1 [expr { abs($y-$y1) }]
							
							set z1 [lindex $data2 [expr { $k1 + 7 }]]
							set rid [lindex $data2 [expr { $k1 + 4 }]]
							if { $x1 != $x || $y1 != $y } {
								if { [format "%.0f" $x1]  == [format "%.0f" $icord] && [format "%.0f" $y1] == [format "%.0f" $jcord] && $z1 == $z } {
									puts $h1 "$rid"
								}
							}
						}
						incr k1
					}
				}
			}
			puts $h1 "\}"
		}
		incr k
	}

	close $h1

	set z 0.000
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)

					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					if { $zori < 0.0 } {

						set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						set dum [open "dummy" "w"]

						puts $dum "$tcoord"

						close $dum

						set dum [open "dummy" "r"]
						set dvar [read $dum]
						close $dum

						set x1n($n) [format "%.3f" [lindex $dvar 0]]
						set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
						set sx1 [string length $x1n($n)]

						set y1n($n) [format "%.3f" [lindex $dvar 1]]
						set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
						set sy1 [string length $y1n($n)]

						set z1n($n) [format "%.3f" [lindex $dvar 2]]
						set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
						set sz1 [string length $z1n($n)]
					}
				incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1

	set num_atom $n
	set cum_n $n
	set an1n(-1) 0

	for {set n 0} {$n < 10} {incr n} {
		set tra($n) 0.0
	}
	set g10 [open "grid" "r"]
	set gd [read $g10]
	close $g10

	
	for {set ndl 0} {$ndl < $ndlip} {incr ndl} {
		
		puts ""
		puts "			#### IMPOSING THE SET OF TRASFORMATIONS FROM THE PURE LIPID BILAYER OF LIPID [expr { $ndl + 1 }] ####" 

		set inpres [lindex $data1 13 [expr { 2 + $ndl }]]
		puts "$inpres"
		set k10 0

		while { [lindex $gd $k10] != $inpres } {
			incr k10
		}  
		set gp $k10
		incr k10

		set k11 1

		set frs [expr { $X / $stepx }]
		set frs [format "%.0f" $frs]
		set fcs [expr { $Y / $stepy }]
		set fcs [format "%.0f" $fcs]


		set resad 0
		for {set n 0} {$n < $num_atom} {incr n} {
			if { $an1n([expr { $n -1 }]) != $an1n($n) } {
				puts ""
				puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
				puts " 				_______________________________________"

				set ulgrid [lindex $gd $k10]

				set gtype [lindex $ulgrid 0]

				set k11 [expr { ($resad / $frs) + 1 }]

				set check [overlap $an1n($n)]
				if { $check != 0 } {
					if { $gtype == "A" } {
						if { $k11 == 1} {
							set numft [llength [lindex $ulgrid $k11]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid $k11 $k12]
								incr k12
							}
						}
						if { $k11 >= 2 && $k11 < $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
					} elseif { $gtype == "S" } {
						if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs } {
							set numft [llength [lindex $ulgrid 1]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 1 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
					}	elseif { $gtype == "S1" } {
						if { [expr { $k11 % 2 }] != 0 } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 3 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 1 $k12]
									incr k12
								}
							}
						}
						if { [expr { $k11 % 2 }] == 0 } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 4 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 2]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 $k12]
									incr k12
								}
							}
						}
					} elseif { $gtype == "S2" } {
						if { $k11 == 1 } {
							set numft [llength [lindex $ulgrid 1]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 1 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs && $k11 != 1} {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 4]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 4 $k12]
								incr k12
							}
						}
					} elseif { $gtype == "A1" } {
						if { $k11 == 1} {
							if { $resad == 0 } {
								set numft [llength [lindex $ulgrid $k11 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid $k11 0 $k12]
									incr k12
								}
							} else { 
								set numft [llength [lindex $ulgrid $k11 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid $k11 1 $k12]
									incr k12
								}
							}
						}
						if { $k11 >= 2 && $k11 < $fcs} {
							if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
								set numft [llength [lindex $ulgrid 2 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 0 $k12]
									incr k12
								} 
								} else {
								set numft [llength [lindex $ulgrid 2 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 1 $k12]
									incr k12
								}
							}
						}
						if { $k11 == $fcs } {
							if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
								set numft [llength [lindex $ulgrid 3 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 3 0 $k12]
									incr k12
								} 
								} else {
								set numft [llength [lindex $ulgrid 3 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 3 1 $k12]
									incr k12
								}
							}
						}
					} elseif { $gtype == "S3" } {
						if { [expr { $k11 % 2 }] != 0 && $k11 != $frs } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 3 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 1 $k12]
									incr k12
								}
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $frs} {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 4 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 2]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 $k12]
									incr k12
								}
							}
						}
						if { $k11 == $frs } {
							set numft [llength [lindex $ulgrid 5]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 5 $k12]
								incr k12
							}
						}	
					}
				
					if { $k12 >= $numft } {
						set k12 0
					} else {
						incr k12
					} 

					set tx [lindex $tra($k12) 0]
					set ty [lindex $tra($k12) 1]
					set tz [lindex $tra($k12) 2]

					if { $ty != 0 } {
						for {set tyn 0.0} {$tyn < $ty} {set tyn [expr { $tyn + 0.5}]} {
							transformations 3 0.0 $an1n($n)
							transformations 4 $tyn $an1n($n)
							transformations 2 6.2 $an1n($n)
						}
					}
					transformations 3 $tx $an1n($n)
					transformations 4 $ty $an1n($n)
					transformations 2 $tz $an1n($n)
				}
				incr resad		
			}
		}
	}
	
	puts ""
	puts "				#### TRYING NEW SET OF TRANSFORMATIONS ####"
	puts ""

	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				set otx 0.0
				set oty 0.0
				set otz -3.14
				for {set tx 0.0} {$tx < [expr { $stepx / 2.0 }] } {set tx [expr { $tx + 0.5 }] } {
					for {set ty 0.0} {$ty < [expr { $stepy / 2.0 }] } {set ty [expr { $ty + 0.5 }] } {
						for {set tz -3.14} {$tz <= 3.14} { set tz [expr { $tz + [format "%.2f" [expr { 3.14 / 30.0 }]] }] } {
							#puts "				Tx=$tx	Ty=$ty	Tz=$tz"
							transformations 3 [expr { $tx - $otx }] $an1n($n)
							transformations 4 [expr { $ty - $oty }] $an1n($n)
							transformations 2 [expr { $tz - $otz }] $an1n($n)
							set otz $tz
							set oty $ty
							set otx $tx
							set check [overlap $an1n($n)] 
							if { $check == 0 } {
								puts "				FINAL VALUE :: Tx=$tx	Ty=$ty	Tz=$tz"
								set tx [expr {  $stepx / 2.0  }]
								set ty [expr {  $stepy / 2.0  }]
								#set tz [expr { (10.0 * 3.14) /180.0 }]
								set tz 4.00
							}
						}
						set otz -3.14
					}
				}
			}
			#set check [overlap $an1n($n)]
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}
		}
	}

	# REMOVING THE OVERLAPS IN LOWER LAYER
	puts ""
	puts "				**** REMOVING THE OVERLAPS IN LOWER LAYER **** "

	set z [expr { -1 * $Z }]
	set nr 1
	set n 0
	set cum_n $n
	set k1 0
	while { $k1 < [llength $data2] } {
		set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
		if { [lindex $data2 $k1] == "HETATM" && $zshift($n) == $z } {
			set res [lindex $data2 [expr {$k1 + 3 }]]
			set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
			set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
			set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]
			set pivot [lindex $data2 [expr { $k1 + 2 }]] 

			set k 0 
			while { [lindex $data1 1 $k] != $res } {
				incr k
			}
			set pdb [lindex $data1 0 $k]
			puts "				**** READING PDB $pdb ****"
			set m [open "$pdb" "r"]
			set data [read $m]
			close $m 

			set k 0

			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 6 }]]
					set yori [lindex $data [expr { $k + 7 }]]
					set zori [lindex $data [expr { $k + 8 }]]
				}
				incr k
			}

			set k 0
			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "HETATM" } {

					set xshift($n) [lindex $data2 [expr { $k1 + 5 }]]
					set yshift($n) [lindex $data2 [expr { $k1 + 6 }]]
					set zshift($n) [lindex $data2 [expr { $k1 + 7 }]]

					set rn1n($n) [lindex $data [expr { $k + 1 }]]

					set at1n($n) [lindex $data [expr { $k + 2 }]]
					
					set an1n($n) [lindex $data2 [expr { $k1 + 4 }]]
					set nr $an1n($n)
					#set an1($n) $nr
						

					set x1n($n) [lindex $data [expr { $k + 6 }]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) - $xori + $xshift($n) }]]
						

					set y1n($n) [lindex $data [expr { $k + 7 }]] 
					set y1n($n) [format "%.3f" [expr { $y1n($n) - $yori + $yshift($n)}]]

					set z1n($n) [lindex $data [expr { $k + 8 }]] 
					set z1n($n) [format "%.3f" [expr { $z1n($n) - $zori + $zshift($n)}]]
						

					set resnamen($n) [lindex $data [expr { $k + 3 }]]

					set atomnamen($n) [lindex $data [expr { $k + 11 }]]
					set satomnamen [string length $atomnamen($n)]

					set coord [list [expr { $x1n($n) - $xshift($n) }] [expr { $y1n($n) - $yshift($n) }] [expr { $z1n($n) - $zshift($n) }]]

					if { $zori > 0.0 } {

						# TRANSFORMATION 1

						set tcoord [::math::linearalgebra::matmul $tm $coord]

						# TRANSFORMATION 2

						set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

						# TRANSFORMATION 3

						set tcoord [::math::linearalgebra::add $tcoord $ztrans]
					} else { 
						set tcoord [::math::linearalgebra::add $coord $ztrans]
					}

					set dum [open "dummy" "w"]

					puts $dum "$tcoord"

					close $dum

					set dum [open "dummy" "r"]
					set dvar [read $dum]
					close $dum

					set x1n($n) [format "%.3f" [lindex $dvar 0]]
					set x1n($n) [format "%.3f" [expr { $x1n($n) + $xshift($n)}]]
					set sx1 [string length $x1n($n)]

					set y1n($n) [format "%.3f" [lindex $dvar 1]]
					set y1n($n) [format "%.3f" [expr { $y1n($n) + $yshift($n)}]]
					set sy1 [string length $y1n($n)]

					set z1n($n) [format "%.3f" [lindex $dvar 2]]
					set z1n($n) [format "%.3f" [expr { $z1n($n) + $zshift($n)}]]
					set sz1 [string length $z1n($n)]

					incr n
				}
			incr k
			}
		}
		incr k1
	}
	set an1n($n) -1
	set num_atom $n
	set cum_n $n
	set an1n(-1) 0

	for {set ndl 0} {$ndl < $ndlip} {incr ndl} {
		
		puts ""
		puts "			#### IMPOSING THE SET OF TRASFORMATIONS FROM THE PURE LIPID BILAYER OF LIPID [expr { $ndl + 1 }] ####" 

		set inpres [lindex $data1 13 [expr { 2 + $ndl }]]

		set k10 0

		while { [lindex $gd $k10] != $inpres } {
			incr k10
		}  
		set gp $k10
		incr k10

		set k11 1

		set frs [expr { $X / $stepx }]
		set frs [format "%.0f" $frs]
		set fcs [expr { $Y / $stepy }]
		set fcs [format "%.0f" $fcs]

		set resad 0
		for {set n 0} {$n < $num_atom} {incr n} {
			if { $an1n([expr { $n -1 }]) != $an1n($n) } {
				puts ""
				puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
				puts " 				_______________________________________"

				set ulgrid [lindex $gd [expr { $k10 + 1 }]]

				set gtype [lindex $ulgrid 0]

				set k11 [expr { ($resad / $frs) + 1 }]

				set check [overlap $an1n($n)]
		
				if { $check != 0 } {
					if { $gtype == "A" } {
						if { $k11 == 1} {
							set numft [llength [lindex $ulgrid $k11]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid $k11 $k12]
								incr k12
							}
						}
						if { $k11 >= 2 && $k11 < $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
					} elseif { $gtype == "S" } {
						if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs } {
							set numft [llength [lindex $ulgrid 1]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 1 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
					}	elseif { $gtype == "S1" } {
						if { [expr { $k11 % 2 }] != 0 } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 3 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 1 $k12]
									incr k12
								}
							}
						}
						if { [expr { $k11 % 2 }] == 0 } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 4 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 2]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 $k12]
									incr k12
								}
							}
						}
					} elseif { $gtype == "S2" } {
						if { $k11 == 1 } {
							set numft [llength [lindex $ulgrid 1]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 1 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] != 0 && $k11 != $fcs && $k11 != 1} {
							set numft [llength [lindex $ulgrid 3]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 3 $k12]
								incr k12
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $fcs} {
							set numft [llength [lindex $ulgrid 2]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 2 $k12]
								incr k12
							}
						}
						if { $k11 == $fcs } {
							set numft [llength [lindex $ulgrid 4]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 4 $k12]
								incr k12
							}
						}
					} elseif { $gtype == "A1" } {
						if { $k11 == 1} {
							if { $resad == 0 } {
								set numft [llength [lindex $ulgrid $k11 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid $k11 0 $k12]
									incr k12
								}
							} else { 
								set numft [llength [lindex $ulgrid $k11 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid $k11 1 $k12]
									incr k12
								}
							}
						}
						if { $k11 >= 2 && $k11 < $fcs} {
							if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
								set numft [llength [lindex $ulgrid 2 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 0 $k12]
									incr k12
								} 
								} else {
								set numft [llength [lindex $ulgrid 2 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 1 $k12]
									incr k12
								}
							}
						}
						if { $k11 == $fcs } {
							if { $resad != 0 && [expr { $resad % $frs }] == 0 } {
								set numft [llength [lindex $ulgrid 3 0]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 3 0 $k12]
									incr k12
								} 
								} else {
								set numft [llength [lindex $ulgrid 3 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 3 1 $k12]
									incr k12
								}
							}
						}
					} elseif { $gtype == "S3" } {
						if { [expr { $k11 % 2 }] != 0 && $k11 != $frs } {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 3 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 1]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 1 $k12]
									incr k12
								}
							}
						}
						if { [expr { $k11 % 2 }] == 0 && $k11 != $frs} {
							if { [expr { $frs * $k11 }] == [expr { $resad + 1 }] } {
								set numft 1
								set k12 0
								set tra($k12) [lindex $ulgrid 4 0]
								incr k12
							} else {
								set numft [llength [lindex $ulgrid 2]]
								set k12 0
								while { $k12 < $numft } {
									set tra($k12) [lindex $ulgrid 2 $k12]
									incr k12
								}
							}
						}
						if { $k11 == $frs } {
							set numft [llength [lindex $ulgrid 5]]
							set k12 0
							while { $k12 < $numft } {
								set tra($k12) [lindex $ulgrid 5 $k12]
								incr k12
							}
						}	
					}
				
					if { $k12 >= $numft } {
						set k12 0
					} else {
						incr k12
					} 

					set tx [lindex $tra($k12) 0]
					set ty [lindex $tra($k12) 1]
					set tz [lindex $tra($k12) 2]

					if { $ty != 0 } {
						for {set tyn 0.0} {$tyn < $ty} {set tyn [expr { $tyn + 0.5}]} {
							transformations 3 0.0 $an1n($n)
							transformations 4 $tyn $an1n($n)
							transformations 2 6.2 $an1n($n)
						}
					}
				
					transformations 3 $tx $an1n($n)
					transformations 4 $ty $an1n($n)
					transformations 2 $tz $an1n($n)
				}
			incr resad		
			}
		}
	}

	puts ""
	puts "				#### TRYING NEW SET OF TRANSFORMATIONS ####"
	puts ""

	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"
			set check [overlap $an1n($n)]
			if { $check != 0 } {
				set otx 0.0
				set oty 0.0
				set otz -3.14
				for {set tx 0.0} {$tx < [expr { $stepx / 2.0 }] } {set tx [expr { $tx + 0.5 }] } {
					for {set ty 0.0} {$ty < [expr { $stepy / 2.0 }] } {set ty [expr { $ty + 0.5 }] } {
						for {set tz -3.14} {$tz <= 3.14} { set tz [expr { $tz + [format "%.2f" [expr { 3.14 / 30.0 }]] }] } {
							#puts "				Tx=$tx	Ty=$ty	Tz=$tz"
							transformations 3 [expr { $tx - $otx }] $an1n($n)
							transformations 4 [expr { $ty - $oty }] $an1n($n)
							transformations 2 [expr { $tz - $otz }] $an1n($n)
							set otz $tz
							set oty $ty
							set otx $tx
							set check [overlap $an1n($n)] 
							if { $check == 0 } {
								puts "				FINAL VALUE :: Tx=$tx	Ty=$ty	Tz=$tz"
								set tx [expr {  $stepx / 2.0  }]
								set ty [expr {  $stepy / 2.0  }]
								#set tz [expr { (10.0 * 3.14) /180.0 }]
								set tz 4.00
							}
						}
						set otz -3.14
					}
				}
			}
			#set check [overlap $an1n($n)]
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g1 "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g1 "TER"
			}
		}
	}
}

proc lipid_overlap {} {

	set g [open "input" "r"]
	set data1 [read $g]
	close $g
	

	puts ""
	puts "				****REMOVING THE OVERLAPS ****"
	puts "				------------------------------"
	puts ""
 

	# space variables

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set f [open "lipids.pdb" "r"]
	set data [read $f]
	close $f

	set h [open "lipids_no.pdb" "w"]

	set h1 [open "overlap" "w"]

	set k1 0
	set i 0
	set l 0
	set count 0
	set oldanc 1
	while { $k1 < [llength $data] } {
		if { [lindex $data $k1] == "HETATM" } {
			set xc1($l) [lindex $data [expr { $k1 + 5 }]]
			set yc1($l) [lindex $data [expr { $k1 + 6 }]]
			set zc1($l) [lindex $data [expr { $k1 + 7 }]] 
			set anc1($l) [lindex $data [expr { $k1 + 4 }]]
			set atom [lindex $data [expr { $k1 + 1 }]]
	

			puts "				**** RESIDUE $anc1($l) ATOM $atom ****"
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" } {
					set anc2 [lindex $data [expr { $k + 4 }]]
					if { $anc2 != $anc1($l) } {
						set xc2 [lindex $data [expr { $k + 5 }]]
						set yc2 [lindex $data [expr { $k + 6 }]]
						set zc2 [lindex $data [expr { $k + 7 }]]
						set delx [expr { $xc1($l) - $xc2 }]
						set delx2 [expr { $delx * $delx }]
						set dely [expr { $yc1($l) - $yc2 }]
						set dely2 [expr { $dely * $dely }]
						set delz [expr { $zc1($l) - $zc2 }]
						set delz2 [expr { $delz * $delz }]
						set del [expr { $delx2 + $dely2 + $delz2 }]
						if { $del < 1.0 } {
							puts $h1 "$anc1($l) $anc2"
						}	
					}
				}
				incr k
			}
			incr l
		}
		incr k1 
	}

	close $h1

	set h1 [open "overlap" "r"]
	set ov [read $h1]
	close $h1

	set h2 [open "overlap_sorted" "w"]

	set k 0
	set i 0
	set t1(-1) -1
	set count1 0
	after 1000
	while { $k < [llength $ov] } {
		set t1($i) [lindex $ov $k]
		set t2($i) [lindex $ov [expr { $k + 1 }]]
		if { $t1([expr { $i-1 }]) != $t1($i) } {
			if { $i > 0 } {
				puts $h2 "OVERLAP_PAIR $count1" 
			}
			set count1 0
		}
		set count 0
		for {set j 0} {$j < $i} {incr j} {
			if { $t1($i) == $t1($j) && $t2($i) == $t2($j) } {
				incr count
			}
		}
		if { $count == 0} {
			puts $h2 "$t1($i) $t2($i)"
			incr count1

		}
		incr k 2
		incr i
	}
	puts $h2 "OVERLAP_PAIR $count1"
	close $h2

	set h2 [open "overlap_sorted" "r"]
	set ov [read $h2]
	close $h2

	set k 0
	set i 0
	set j 0
	set rem(0) -1
	
	while { $k < [llength $ov] } {
		set count 0
		set t1c [lindex $ov $k]
		set t2c [lindex $ov [expr { $k + 1 }]]
		for {set ij 0} {$ij < $j} {incr ij} {
			if { $t1c == $rem($ij) || $t2c == $rem($ij) } {
				incr count
			}
		}
		if { $count == 0 && $t1c != "OVERLAP_PAIR" } {
			set t2c [lindex $ov [expr { $k + 1 }]]

			set k3 $k
			while { [lindex $ov $k3] != "OVERLAP_PAIR" } {	
				incr k3
			}	
			incr k3
			set pc1 [lindex $ov $k3]
		
			set k1 0
			while { $k1 < [llength $ov] } {
				set t1c1 [lindex $ov $k1]
				if { $t1c1 == $t2c } {
					set k4 $k1
					while { [lindex $ov $k4] != "OVERLAP_PAIR" } {	
						incr k4
					}
					incr k4
					set pc2 [lindex $ov $k4] 
					set k1 [llength $ov]
				}
				incr k1 2
			}
			if { $pc2 < $pc1 } {
				set rem($j) $t1c	
				incr j	
				incr i	
			} else {
				set rem($j) $t2c
				incr i
				incr j
			}
		}
		incr k 2
	}

	set num_removal $j
	puts  "$num_removal"
	for {set i 0} {$i < $j} {incr i} {
		puts "REMOVING	$rem($i)"
	}

	set k 0
	while { [lindex $data $k] != "END" } {
		if { [lindex $data $k] == "TER" && $count == 0} {
			puts $h "TER"
		}
		set count 0
		if { [lindex $data $k] == "HETATM" } {
			set rn1 [lindex $data [expr { $k + 1 }]]
			set srn1 [string length $rn1]

			set at1 [lindex $data [expr { $k + 2 }]]
			set sat1 [string length $at1]

			set an1 [lindex $data [expr { $k + 4 }]]
			set san1 [string length $an1]

			set x1 [lindex $data [expr { $k + 5 }]]
			set x1 [format "%.3f" $x1]
			set sx1 [string length $x1]

			set y1 [lindex $data [expr { $k + 6 }]] 
			set y1 [format "%.3f" $y1]
			set sy1 [string length $y1]

			set z1 [lindex $data [expr { $k + 7 }]] 
			set z1 [format "%.3f" $z1]
			set sz1 [string length $z1]

			set sresname [string length [lindex $data [expr { $k + 3 }]]]

			for {set i 0} {$i < $num_removal} {incr i} { 
				if { $an1 == $rem($i) } {
					incr count
				}
			}
			if { $count == 0 } {
				puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
			}
			incr k 11
		} else { 
			incr k
		}	
	}
	#puts $h "END"
	close $h
}

proc lipid_hole {} {

	set g [open "input" "r"]
	set data1 [read $g]
	close $g
	
	set hole_radiix [lindex $data1 6 1]
	set hole_radiiy [lindex $data1 6 2]

	set X [lindex $data1 5 1]
	set Y [lindex $data1 5 2]
	set Z [lindex $data1 5 3]

	set xmin [expr { ($X / 2.0) - $hole_radiix }]
	set ymin [expr { ($Y / 2.0) - $hole_radiiy }]
	set xmax [expr { ($X / 2.0) + $hole_radiix }]
	set ymax [expr { ($Y / 2.0) + $hole_radiiy }]

	puts ""
	puts "				**** CARVING OUT A HOLE WITH $xmin < X < $xmax && $ymin < Y < $ymax ****"
	puts "				-------------------------------------------------------------------------"
	puts ""
 

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set f [open "lipids_no.pdb" "r"]
	set data [read $f]
	close $f

	set h [open "lip_hole.pdb" "w"]

	set k1 1
	set i 0

	while { $k1 < [llength [lindex $data1 0]] } {
		set pivot [lindex $data1 2 $k1] 
			set k 0
			while { $k < [llength $data] } {
				if { [lindex $data $k] == "HETATM" && [lindex $data [expr { $k + 2 }]] == "$pivot" } {
					set xori [lindex $data [expr { $k + 5 }]]
					set yori [lindex $data [expr { $k + 6 }]]
					set tagres [lindex $data [expr { $k + 4 }]]
					if { $xori > $xmin && $xori < $xmax && $yori > $ymin && $yori < $ymax } {
						set tag($i) $tagres
						incr i
					}
				}
				incr k
			}
		incr k1
	}
	set num_tag $i

	set k 0
	while { $k < [llength $data] } {
		if { [lindex $data $k] == "TER" && $count == 0} {
			puts $h "TER"
		}
		set count 0
		if { [lindex $data $k] == "HETATM" } {
			set rn1 [lindex $data [expr { $k + 1 }]]
			set srn1 [string length $rn1]

			set at1 [lindex $data [expr { $k + 2 }]]
			set sat1 [string length $at1]

			set an1 [lindex $data [expr { $k + 4 }]]
			set san1 [string length $an1]

			set x1 [lindex $data [expr { $k + 5 }]]
			set x1 [format "%.3f" $x1]
			set sx1 [string length $x1]

			set y1 [lindex $data [expr { $k + 6 }]] 
			set y1 [format "%.3f" $y1]
			set sy1 [string length $y1]

			set z1 [lindex $data [expr { $k + 7 }]] 
			set z1 [format "%.3f" $z1]
			set sz1 [string length $z1]

			set sresname [string length [lindex $data [expr { $k + 3 }]]]

			for {set i 0} {$i < $num_tag} {incr i} { 
				if { $an1 == $tag($i) } {
					incr count
				}
			}
			if { $count == 0 } {
				puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
			}
			incr k 11
		} else { 
			incr k
		}	
	}
	puts $h "END"
	close $h
}

proc protein_insersion {} {	

	package require math::linearalgebra

	set inp [open "input" "r"]
	set mn [read $inp]
	close $inp

	set pdb [lindex $mn 7 1]
	
	puts ""
	puts "				**** INSERTING A PROTEIN WITH PDB NAME $pdb INSIDE THE LIPID BILAYER ****"
	puts "				---------------------------------------------------------------------------"
	puts ""

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

	set theta [lindex $mn 8 4]


	set tmr1 [list 1.0 0.0 0.0] 
	set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
	set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
	set tm [list $tmr1 $tmr2 $tmr3]

	# TRANSFORMATION 2 Z AXIS

	set theta [lindex $mn 8 6]

	set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
	set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
	set tm2r3 [list 0.0 0.0 1.0]
	set tm2 [list $tm2r1 $tm2r2 $tm2r3]

	# TRANSFORMATION 3 Y AXIS

	set theta [lindex $mn 8 5]

	set tm3r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
	set tm3r2 [list 0.0 1.0 0.0]
	set tm3r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
	set tm3 [list $tm3r1 $tm3r2 $tm3r3]

	if { [lindex $mn 10 1] == 1 } {
		set f [open "lipids_no.pdb" "r"]
		set data [read $f]
		close $f
	} else { 
		set f [open "lip_hole.pdb" "r"]
		set data [read $f]
		close $f
	}

	set g [open "$pdb" "r"]
	set data1 [read $g]
	close $g

	set h [open "lip_pro.pdb" "w"]

	set hh [open "only_mod_protein.pdb" "w"]

	# FIXING THE COORDINATE TO THE FIRST ATOM OF LIPID FILE

	if { [lindex $mn 12 1] == 1 || [lindex $mn 12 1] == 2 || [lindex $mn 12 1] == 3 || [lindex $mn 12 1] == 4 } {
		set xoris [lindex $data1 6]
		set yoris [lindex $data1 7]
		set zoris [lindex $data1 8]
	} else {
		set xoris 0.0
		set yoris 0.0
		set zoris 0.0
	}

	# CHANGING THE COORDINATES OF LIPIDS AND WRITING INTO A NEW FILE

	set k 0

	while { $k < [llength $data] } {
		if { [lindex $data $k] == "TER" } {
			puts $h "TER"
		}
		if { [lindex $data $k] == "HETATM" } {
			set rn1 [lindex $data [expr { $k + 1 }]]
			set srn1 [string length $rn1]

			set at1 [lindex $data [expr { $k + 2 }]]
			set sat1 [string length $at1]

			set an1 [lindex $data [expr { $k + 4 }]]
			set san1 [string length $an1]

			set x1 [lindex $data [expr { $k + 5 }]]
			set x1 [format "%.3f" [expr { $x1 - 0.0 }]]
			set sx1 [string length $x1]

			set y1 [lindex $data [expr { $k + 6 }]] 
			set y1 [format "%.3f" [expr { $y1 - 0.0 }]]
			set sy1 [string length $y1]

			set z1 [lindex $data [expr { $k + 7 }]] 
			set z1 [format "%.3f" [expr { $z1 - 0.0 }]]
			set sz1 [string length $z1]

			set sresname [string length [lindex $data [expr { $k + 3 }]]]
			set sresname [expr { $sresname + 1 }]

			puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$satn($sat1)[lindex $data [expr { $k + 3 }]] $p($sresname)a$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 10 }]]"
			incr k 11
		} else { 
			incr k
		}	
	}

	set xori [expr { $xoris + [lindex $mn 8 1] }]
	set yori [expr { $yoris + [lindex $mn 8 2] }]
	set zori [expr { $zoris + [lindex $mn 8 3] }]

	set tcoord [list $xori $yori $zori]

	# TRANSFORMATION 1

	#set tcoord [::math::linearalgebra::matmul $tm $tcoord]

	# TRANSFORMATION 2

	#set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

	# TRANSFORMATION 3

	#set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

	set dum [open "dummy" "w"]
	
	puts $dum "$tcoord"

	close $dum

	set dum [open "dummy" "r"]
	set dvar [read $dum]
	close $dum

	set xori [format "%.3f" [lindex $dvar 0]]
	set yori [format "%.3f" [lindex $dvar 1]]
	set zori [format "%.3f" [lindex $dvar 2]]
				
	set k 0

	while { $k < [llength $data1] } {
		set term [lindex $data1 $k]
		set t1 [string range $term 0 5]
		if { $term == "TER" } {
			puts $h "TER"
			puts $hh "TER"
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

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori }]]
				set x1 
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori }]]
				set sy1 [string length $y1]

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori }]]
				set sz1 [string length $z1]
			}

			if { $shift3 != 0 } {
				set tcoord [list $corx $cory $z1] 
		
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

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori }]]
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori }]]
				set sy1 8
	
				set x1 $x1$y1
				set y1 ""

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori }]]
				set sz1 [string length $z1]
			}

			if { $shift4 != 0 } {
				set tcoord [list $x1 $cory $corz] 
		
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

				set x1 [format "%.3f" [expr { [lindex $dvar 0] - $xori }]]
				set sx1 [string length $x1]

				set y1 [format "%.3f" [expr { [lindex $dvar 1] - $yori }]]
				set sy1 [string length $y1]

				set z1 [format "%.3f" [expr { [lindex $dvar 2] - $zori }]]
				set sz1 8

				set y1 $y1$z1
				set z1 ""
			}
		
			puts $h "$ft$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$p($sresn)$resn $cn$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data1 [expr { $k + 11 - $shift - $shift1 - $shift2 - $shift3 - $shift4 - $shift5 }]]"
			puts $hh "$ft$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$p($sresn)$resn $cn$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data1 [expr { $k + 11 - $shift - $shift1 - $shift2 - $shift3 - $shift4 - $shift5 }]]"
		}
		incr k
	}
	puts $h "TER"
	puts $h "END"
	close $h
	puts $hh "TER"
	puts $hh "END"
	close $hh
}	

proc replace {} {

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set check [open "overlap_atom" "w"]

	set g [open "input" "r"]
	set data1 [read $g]
	close $g

	set pro [open "[lindex $data1 7 1]" "r"]
	set pr [read $pro]
	close $pro

	set k 0

	while { [lindex $pr $k] != "ATOM" && [lindex $pr $k] != "HETATM"} {
		incr k
	}

	set resc [lindex $pr [expr { $k + 3 }]] 

	# IN THIS PROCEDURE WE INSERT THE PROTEIN FIRST AND THEN REMOVES ANY LIPID WHICH IS 1A FROM THE PROTEIN

	set f [open "lip_pro.pdb" "r"]
	set data [read $f]
	close $f

	set lp1 [llength $data]
	set lp2 [llength $pr]
	set lp [expr { $lp1 - $lp2 }]

	set h [open "lip_pro_final.pdb" "w"]

	# READING THE END OF THE PDB WHICH CONTAIN INFORMATION ABOUT THE LIPIDS


	set k 0
	set old_rn -1	

	while { [lindex $data [expr { $k + 3 }]] != "$resc" } {
		if { [lindex $data $k] == "HETATM" } {
			set rn1 [lindex $data [expr { $k + 5 }]]
			set slrn1 [string length $rn1]
			if { $slrn1 >= 4 } {
				set shift1 1
				set rn1 [string range [lindex $data [expr { $k + 4 }]] 1 4]
			} else {
				set shift1 0
			}
			set resn1 [lindex $data [expr { $k + 3 }]]
			if { $rn1 != $old_rn } {
				puts "				**** CHECKING RESIDUE $rn1 ****"
				set old_rn $rn1
				set k1 $k
				set count 0
				while { [lindex $data $k1] != "TER" } {
					if { [lindex $data $k1] == "HETATM" } {
						set xl [lindex $data [expr { $k1 + 6 - $shift1}]]
						set yl [lindex $data [expr { $k1 + 7 - $shift1}]]
						set zl [lindex $data [expr { $k1 + 8 - $shift1}]]
						set k2 $lp
						while { $k2 < [llength $data] } {

							if { [lindex $data $k2] == "HETATM" || [lindex $data $k2] == "ATOM" } {

								set rn2 [lindex $data [expr { $k2 + 5 }]]
								set slrn2 [string length $rn2]
								if { $slrn2 >= 4 } {
									set shift2 1
									set rn2 [string range [lindex $data [expr { $k2 + 4 }]] 1 4]
								} else {
									set shift2 0
								}
								set resn2 [lindex $data [expr { $k2 + 3 }]]

								if { $rn1 != $rn2 } {
									set tm 1
									set ct 0
									while { $tm < [llength [lindex $data1 11]] } {
										if { $resn2 == [lindex $data1 11 $tm] } {
											incr ct
										}
										incr tm
									}
									if { $ct == 0 } {

										set xc [lindex $data [expr { $k2 + 6 - $shift2}]]
										set sxc [string length $xc]

										if { $sxc > 8 } {
											set t 0
											while { [string range $xc $t $t] != "." } {
												incr t
											}
											set xcc [string range $xc 0 [expr { $t + 3 }]]
											set ycc [string range $xc [expr { $t + 4 }] end]
											set zcc [lindex $data [expr { $k + 4 }]]
										} else { 
											set xcc $xc
											set yc [lindex $data [expr { $k2 + 7 - $shift2}]]
											set syc [string length $yc]
											if { $syc > 8 } {
												set t 0
												while { [string range $yc $t $t] != "." } {
													incr t
												}
												set ycc [string range $yc 0 [expr { $t + 3 }]]
										 	 	set zcc [string range $yc [expr { $t + 4 }] end]
											}	else {
												set ycc [lindex $data [expr { $k2 + 7 - $shift2}]]
												set zcc [lindex $data [expr { $k2 + 8 - $shift2}]]
											}
										}
					
										set xc $xcc
										set yc $ycc
										set zc $zcc

										set delx [expr { $xc - $xl }]
										set delx [expr { $delx * $delx }]
										set dely [expr { $yc - $yl }]
										set dely [expr { $dely * $dely }]
										set delz [expr { $zc - $zl }]
										set delz [expr { $delz * $delz }]
		
										# CALCULATING THE DISTANCE BETWEEN THE LIPID ATOM AND THE ATOM FROM THE PDB
			
										set del [expr { $delx + $dely + $delz }]

										if { $del < 1.44 } {
											puts ""
											puts "			#### RESIDUE $resn1 AND $resn2 ARE OVERLAPING ####"
											puts "" 
											puts $check "$rn1 $rn2"
											puts ""
											puts "				**** REMOVING RESIDUE $rn1 ****"
											incr count
											set k2 [llength $data]
											while { [lindex $data $k1] != "TER" } {
												incr k1
											}
											set k $k1
											set k1 [expr { $k1 - 1 }]
										}
									}
								}
							}
							incr k2
						}
					}
					incr k1
				}
			}
			if { $count == 0 } {
				while { [lindex $data $k] != "TER" } {
					if { [lindex $data $k] == "HETATM" } {
						set rn1 [lindex $data [expr { $k + 1 }]]
						set srn1 [string length $rn1]

						set at1 [lindex $data [expr { $k + 2 }]]
						set sat1 [string length $at1]

						set cn [lindex $data [expr { $k + 4 }]]

						set an1 [lindex $data [expr { $k + 5 }]]
						set san1 [string length $an1]
						if { $san1 >= 4 } {
							set shift 1
							set an1 [string range [lindex $data [expr { $k + 4 }]] 1 4]
							set san1 4
							set cn a
							} else {
								set shift 0
							}

						set x1 [lindex $data [expr { $k + 6 - $shift}]]
						set x1 [format "%.3f" $x1]
						set sx1 [string length $x1]

						set y1 [lindex $data [expr { $k + 7 - $shift}]] 
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]

						set z1 [lindex $data [expr { $k + 8 -$shift }]] 
						set z1 [format "%.3f" $z1]
						set sz1 [string length $z1]

						set resn [lindex $data [expr { $k + 3 }]]
						set sresn [string length [lindex $data [expr { $k + 3 }]]]
						set sresn [expr { $sresn + 1 }]
				
						puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$satn($sat1)$resn $p($sresn)a$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data [expr { $k + 11 - $shift }]]"
					}
				incr k
				}
				puts $h "TER"
			}
		}
		incr k
	}
	#puts $h "TER"
	#puts $h "END"
	close $h

	set h1 [open "lip_pro_final.pdb" "r"]
	set dat1 [read $h1]
	close $h1

	set h2 [open "only_mod_protein.pdb" "r"]
	set dat2 [read $h2]
	close $h2

	set dat3 $dat1$dat2

	set h3 [open "1_final_struct.pdb" "w"]

	puts $h3 "$dat3"
	close $h3
}


proc transformations { flag theta nr} {

	package require math::linearalgebra

	#puts "transformations"

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n
	global xshift
	global yshift
	global zshift
	global resnamen

	set k 0
	while { $an1n($k) != $nr } { 
		incr k
	}
	set k1 $k

		
	# TRANSFORMATIONS

	if { $flag == 0 } {
		
		# TRANSFORMATION 1 X AXIS

		set tmr1 [list 1.0 0.0 0.0] 
		set tmr2 [list 0.0 [expr { cos($theta) }] [expr { -1 * (sin($theta))}]] 
		set tmr3 [list 0.0 [expr { sin($theta) }] [expr { cos($theta) }]]
		set tm [list $tmr1 $tmr2 $tmr3]

		while { $an1n($k) == $nr } {

			set tcoord [list [expr { $x1n($k) - $xshift($k) }] [expr { $y1n($k) - $yshift($k) }] [expr { $z1n($k) - $zshift($k) }]]
			set tcoord [::math::linearalgebra::matmul $tm $tcoord]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set x1n($k) [format "%.3f" [expr { $x1n($k) + $xshift($k)}]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set y1n($k) [format "%.3f" [expr { $y1n($k) + $yshift($k)}]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set z1n($k) [format "%.3f" [expr { $z1n($k) + $zshift($k)}]]
			set sz1 [string length $z1n($k)]

			incr k
		}
	}

	set k $k1

	if { $flag == 1 } {

		# TRANSFORMATION 2 	Y AXIS

		set tm2r1 [list [expr { cos($theta) }] 0.0 [expr { sin($theta) }]] 
		set tm2r2 [list 0.0 1.0 0.0]
		set tm2r3 [list [expr { -1 * (sin($theta))}] 0.0 [expr { cos($theta) }]]
		set tm2 [list $tm2r1 $tm2r2 $tm2r3]

		while { $an1n($k) == $nr } {

			set tcoord [list [expr { $x1n($k) - $xshift($k) }] [expr { $y1n($k) - $yshift($k) }] [expr { $z1n($k) - $zshift($k) }]]
			set tcoord [::math::linearalgebra::matmul $tm2 $tcoord]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set x1n($k) [format "%.3f" [expr { $x1n($k) + $xshift($k)}]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set y1n($k) [format "%.3f" [expr { $y1n($k) + $yshift($k)}]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set z1n($k) [format "%.3f" [expr { $z1n($k) + $zshift($k)}]]
			set sz1 [string length $z1n($k)]

			incr k
		}
	}

	set k $k1 
			
	if { $flag == 2 } {

		# TRANSFORMATION 3 Z AXIS

		set tm3r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
		set tm3r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
		set tm3r3 [list 0.0 0.0 1.0]
		set tm3 [list $tm3r1 $tm3r2 $tm3r3]

		while { $an1n($k) == $nr } {

			set tcoord [list [expr { $x1n($k) - $xshift($k) }] [expr { $y1n($k) - $yshift($k) }] [expr { $z1n($k) - $zshift($k) }]]
			set tcoord [::math::linearalgebra::matmul $tm3 $tcoord]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set x1n($k) [format "%.3f" [expr { $x1n($k) + $xshift($k)}]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set y1n($k) [format "%.3f" [expr { $y1n($k) + $yshift($k)}]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set z1n($k) [format "%.3f" [expr { $z1n($k) + $zshift($k)}]]
			set sz1 [string length $z1n($k)]

			incr k
		}
	}
	set k $k1

	if { $flag == 3 } {
		set xtrans [list $theta 0 0]

		while { $an1n($k) == $nr } {

			set tcoord [list [expr { $x1n($k) - $xshift($k) }] [expr { $y1n($k) - $yshift($k) }] [expr { $z1n($k) - $zshift($k) }]]

			set tcoord [::math::linearalgebra::add $tcoord $xtrans]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set x1n($k) [format "%.3f" [expr { $x1n($k) + $xshift($k)}]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set y1n($k) [format "%.3f" [expr { $y1n($k) + $yshift($k)}]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set z1n($k) [format "%.3f" [expr { $z1n($k) + $zshift($k)}]]
			set sz1 [string length $z1n($k)]

			incr k

		}
	}
	set k $k1

	if { $flag == 4 } {

		set ytrans [list 0 $theta 0]

		while { $an1n($k) == $nr } {

			set tcoord [list [expr { $x1n($k) - $xshift($k) }] [expr { $y1n($k) - $yshift($k) }] [expr { $z1n($k) - $zshift($k) }]]

			set tcoord [::math::linearalgebra::add $tcoord $ytrans]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set x1n($k) [format "%.3f" [expr { $x1n($k) + $xshift($k)}]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set y1n($k) [format "%.3f" [expr { $y1n($k) + $yshift($k)}]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set z1n($k) [format "%.3f" [expr { $z1n($k) + $zshift($k)}]]
			set sz1 [string length $z1n($k)]

			incr k
		}
	}
}

proc overlap {nr} {

	#puts "overlap"

	set f [open "neighbours" "r"]
	set data [read $f]
	close $f

	set k 0
	while { [lindex $data $k 0] != $nr } {
		incr k
	}
	
	set pos $k
	set num_nei [llength [lindex $data $k]]

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n

	set k 0
	while { $an1n($k) != $nr } { 
		incr k
	}

	while { $an1n($k) == $nr } {
		set x1c $x1n($k)
		set y1c $y1n($k)
		set z1c $z1n($k)
		for {set i 0} {$i < $num_nei} {incr i} {
			set k1 0
			while { $an1n($k1) != [lindex $data $pos $i] } {
				#puts "[lindex $data $pos $i] $an1n($k1)"
				incr k1
			}
			set nr1 $an1n($k1)
			while { $an1n($k1) == $nr1 } { 
				if { $an1n($k1) != $nr } {
					set x2c $x1n($k1)
					set y2c $y1n($k1)
					set z2c $z1n($k1)

					set delx [expr { $x2c - $x1c }]
					set delx2 [expr { $delx * $delx }]

					set dely [expr { $y2c - $y1c }]
					set dely2 [expr { $dely * $dely }]

					set delz [expr { $z2c - $z1c }]
					set delz2 [expr { $delz * $delz }]

					set del [expr { $delx2 + $dely2 + $delz2 }]

					if { $del < 1.0 } {
						return -1
					}
				}
			incr k1
			}
		}
	incr k
	}
	return 0
}

proc lip_comp {fil} {
	
	# THIS PROCEDURE CALCULATES THE COMPOSITION OF THE LIPID BILAYER AFTER THE INSERSION OF THE PROTEIN
	puts ""
	puts "			**** CHECKING THE FINAL COMPOSITION OF LIPIDS ****"

	set f [open "$fil" "r"]
	set data [read $f]
	close $f

	set in [open "input" "r"]
	set inp [read $in]
	close $in

	set pro [open "[lindex $inp 7 1]" "r"]
	set pr [read $pro]
	close $pro

	set k 0

	while { [lindex $pr $k] != "ATOM" && [lindex $pr $k] != "HETATM"} {
		incr k
	}

	set resc [lindex $pr [expr { $k + 3 }]] 

	set k 1

	while { $k < [llength [lindex $inp 0]] } {
		incr k
	}
	set k [expr { $k - 1 }]

	set num_lips $k

	# LIPID NAMES

	set j 0
	for {set i 0} {$i < $k} {incr i} {
		set lp($i,$j) [lindex $inp 1 [expr { $i + 1}]]
	}

	# PIVOT ATOM

	for {set i 0} {$i < $k} {incr i} {
		set p($i) [lindex $inp 2 [expr { $i + 1}]]
	}

	# INITIAL LIPID COMPOSITION UPPER LEAFLET

	set j 1
	for {set i 0} {$i < $k} {incr i} {
		set lp($i,$j) [lindex $inp 3 [expr { $i + 1}]]
	}

	# INITIAL LIPID COMPOSITION LOWER LEAFLET

	set j 2
	for {set i 0} {$i < $k} {incr i} {
		set lp($i,$j) [lindex $inp 4 [expr { $i + 1}]]
	}

	# FINAL LIPID COMPOSITION IN THE UPPER LEAFLET

	set j 3
	for {set i 0} {$i < $k} {incr i} {
		set lp($i,$j) 0
	}

	# FINAL LIPID COMPOSITION IN THE UPPER LEAFLET

	set j 4
	for {set i 0} {$i < $k} {incr i} {
		set lp($i,$j) 0
	}

	# LIPID COMPOSITION AFTER INSERSION OF PROTEIN

	set k 0

	while { [lindex $data [expr { $k + 3 }]] != $resc } {
		if { [lindex $data $k] == "HETATM" } {
			set resnum [lindex $data [expr { $k + 5 }]]
			set slrn [string length $resnum]
			if { $slrn > 4 } {
				set shift 1
			} else { 
				set shift 0
			}
			set res [lindex $data [expr { $k + 3 }]]
			set count 0
			for {set i 0} {$i < $num_lips} {incr i} {
				if { $res == $lp($i,0) } {
					incr count
				}
			}
			if { $count != 0 } {
				set k1 0
				while { $res != $lp($k1,0) } {
					incr k1
				}

				set k2 $k
			
				while { [lindex $data $k2] != $p($k1) } {
					incr k2
				}
				set z [lindex $data [expr { $k2 + 6 - $shift }]]

				if { $z < 0.0 } {
					set lp($k1,4) [expr { $lp($k1,4) + 1 }]
				} else { 
					set lp($k1,3) [expr { $lp($k1,3) + 1 }]
				}
			
				set k2 $k
				while { [lindex $data $k2] != "TER" } {
					incr k2
				}
				set k $k2
			}
		}
		incr k
	}

	puts ""
	puts ""
	puts "			**** CHECK CAREFULLY ****"
	puts "			#### LIPID COMPOSITION BEFORE AND AFTER PROTEIN INSERTION ####"
	puts "			    LIPID		INITIAL(UL)		FINAL(UL)		INITIAL(LL)		FINAL(LL)"
	for {set i 0} {$i < $num_lips} {incr i} {
		puts "		  	    $lp($i,0)			$lp($i,1)			$lp($i,3)			$lp($i,2)			$lp($i,4)"
	}

}

proc mem_builder {} {

	puts "	************************************************"
	puts "	THIS IS AMBER BASED LIPID BILAYER BUILDER FOR"
	puts "	LIPID BILAYER AND TRANSMEMBRANE PROTEIN "
	puts "	SIMULATIONS SETUP"
	puts "	************************************************"

	puts ""
	puts ""
	puts "	DEVELOPED BY TARUN KHANNA AND DR. IAN GOULD"
	puts "	IMPERIAL COLLEGE LONDON, U.K."
	puts ""
	puts "	**********************************************************"
	puts "	The research leading to AMBAT has received funding from" 
	puts "	the People Programme (Marie Curie Actions) of the European" 
	puts "	Union’s	Seventh Framework Programme (FP7/2007-2013) " 
	puts "	under REA grant agreement no 607466."
	puts "	**********************************************************"


	# PARAMETERS 

	# VARIOUS LIPIDS AND THEIR HEAD AND TAIL
	set p1 { DLPC DMPC DOPC DPPC POPC DLPE DMPE DOPE DPPE POPE DOPG DOPS POPS POPG CHL} 
	# HEAD
	set p2 { PC PC PC PC PC PE PE PE PE PE PG PS PS PG CHL}
	# TAIL
	set p3 { LA MY OL PA PA LA MY OL PA PA OL OL PA PA CHL}
	# TAIL
	set p4 { LA MY OL PA OL LA MY OL PA OL OL OL OL OL CHL}

	# LINKER ATOMS
	set p5 { O21 O21 O11 O21 O21 O21 O21 O21 O21 O21 O11 O11 O21 O21 O1}

	# AREA PER LIPID
	set p6 { 67.24 67.1 72.7 68.0 71.3 65.0 63.9 67.4 63.0 62.8 70.4 74.6 64.4 66.0 45.0}

	# Z DIMENSIONS
	set p7 { 15.0 18.0 19.0 18.0 19.0 15.0 18.0 19.0 18.0 19.0 19.0 19.0 19.0 19.0 16.0}
	 
	# USER INPUTS

	# VALID INPUTS
	set vopt1 { 0 1 2 3 4 5 6 7 8}
	set vopt2 { 0 1 2 3 4 5 6 7 8 }

	puts ""
	puts ""

	puts "		### ENTER THE TASK YOU WANT TO PERFORM"
	puts ""
	puts "	CHOOSE ANY OF THE BELOW OPTIONS"

	puts "	0 = BUILD A PURE LIPID BILAYER"
	puts "	1 = INSERT A MOLECULE INSIDE AND BUILD A BILAYER"
	puts "	2 = INSERT A MOLECULE OUTSIDE AND BUILD A BILAYER"
	puts "	3 = INSERT A MOLECULE INSIDE A PREBUILD BILAYER"
	puts "	4 = INSERT A MOLECULE OUTSIDE A PREBUILD BILAYER"
	puts "	5 = INSERT A PROTEIN AND BUILD A BILAYER AROUND IT"
	puts "	6 = INSERT A PROTEIN IN A PREBUILD BILAYER"
	puts "	7 = TO INITIATE THE VESICLE BUILDER"
	puts "	8 = TO INITIATE THE MEMBRANE BUILDER FROM A INPUT FILE AND A LIPID PDB FILE"
	puts "	9 = JUST SOLVATE THE MEMEBRANE SYSTEM "
	puts "	IMP NOTE: PRE BUILD BILAYER SHOULD BE IN THE AMBER FORMAT AND IN THE SAME FOLDER WHERE THE CODE IS EXECUTED"

	set opt1 [gets stdin]

	if { $opt1 == 7 } {
		execution
	} else {
		if { $opt1 != 9 } {
			if { $opt1 != 8 } {
				puts ""

				puts "		### ENTER THE NUMBER OF DIFFERENT LIPIDS IN UPPER LEAFLET"
				set opt2 [gets stdin]
				puts "	(FOLLOW LIPID 14 NOMENCLATURE)"

				for {set i 0} {$i < $opt2} {incr i} {
					puts "	ENTER THE HEAD GROUP OF LIPID [expr { $i + 1 }]"
					set opt2H($i) [gets stdin]

					puts ""
					puts " ENTER FIRST TAIL GROUP OF LIPID [expr { $i + 1 }]"
					set opt2T1($i) [gets stdin]

					puts ""
					puts " ENTER SECOND TAIL GROUP OF LIPID [expr { $i + 1 }]"
					set opt2T2($i) [gets stdin]
					puts ""

					set k 0
					set count 0
					while { $k < [llength $p2] } {
						if { [lindex $p2 $k] == $opt2H($i) } {
							if { [lindex $p3 $k] == $opt2T1($i) && [lindex $p4 $k] == $opt2T2($i) } {
								set lipul($i) [lindex $p1 $k]
								set linkul($i) [lindex $p5 $k]
								incr count
							} elseif { [lindex $p3 $k] == $opt2T2($i) && [lindex $p4 $k] == $opt2T1($i) } {
								set lipul($i) [lindex $p1 $k]
								set linkul($i) [lindex $p5 $k]
								incr count
							}
						}
						incr k
					}
					if { $count == 0 } {
						puts "			### ERROR : LIPID NOT FOUND IN THE DATABASE"
						exit
					}
					puts "		### ENTER THE NUMBER OF $lipul($i) LIPIDS IN UPPER LEAFLET"
					set opt2n($i) [gets stdin] 
				}
				set lipul($i) -1

				puts "		### ENTER THE NUMBER OF DIFFERENT LIPIDS IN LOWER LEAFLET"
				set opt3 [gets stdin]
				puts "  (FOLLOW LIPID 14 NOMENCLATURE)"

				for {set i 0} {$i < $opt3} {incr i} {
					puts "	ENTER THE HEAD GROUP OF LIPID [expr { $i + 1 }]"
					set opt3H($i) [gets stdin]

					puts ""
					puts " ENTER FIRST TAIL GROUP OF LIPID [expr { $i + 1 }]"
					set opt3T1($i) [gets stdin]

					puts ""
					puts " ENTER SECOND TAIL GROUP OF LIPID [expr { $i + 1 }]"
					set opt3T2($i) [gets stdin]
					puts ""

					set k 0
					set count 0
					while { $k < [llength $p2] } {
						if { [lindex $p2 $k] == $opt3H($i) } {
							if { [lindex $p3 $k] == $opt3T1($i) && [lindex $p4 $k] == $opt3T2($i) } {
								set lipll($i) [lindex $p1 $k]
								set linkll($i) [lindex $p5 $k]
								incr count
							} elseif { [lindex $p3 $k] == $opt3T2($i) && [lindex $p4 $k] == $opt3T1($i) } {
								set lipll($i) [lindex $p1 $k]
								set linkll($i) [lindex $p5 $k]
								incr count
							}
						}
						incr k
					}
					if { $count == 0 } {
						puts "			### ERROR : LIPID NOT FOUND IN THE DATABASE"
						exit
					}
					puts "		### ENTER THE NUMBER OF $lipll($i) LIPIDS IN LOWER LEAFLET"
					set opt3n($i) [gets stdin] 
				} 

				set lipll($i) -1

				set count1 0
				for {set i 0} {$i < [llength $vopt1]} {incr i} {
					if { $opt1 == [lindex $vopt1 $i] } {
						incr count1
					}
				}

				set count2 0
				set count3 0

				for {set i 0} {$i < [llength $vopt2]} {incr i} {
					if { $opt2 == [lindex $vopt2 $i] } {
						incr count2
					}
					if { $opt3 == [lindex $vopt2 $i] } {
						incr count3
					}
				}

				if { $count1 == 0 || $count2 == 0 || $count3 == 0 } {
					puts "	### ERROR: INVALID INPUT ###"
					puts "		EXITING THE BUILDER"
					exit
				}

				# USING A PRE BUILT BILAYER

				if { $opt1 == 3 || $opt1 == 4 || $opt1 == 6 } {
					puts ""
					puts "### ENTER THE NAME OF THE PRE BUILT BILAYER"
					set opt5 [gets stdin]
					exec ls $opt5

					# SHIFTING THE ORIGIN TO THE PIVOT ATOM OF THE MOST NEGATIVE RESIDUE

					set f [open "$opt5" "r"]
					set data [read $f]
					close $f
	
					set k 0
					set xmin 1000.0
					set ymin 1000.0
					set zmax 0.000
					while { $k < [llength $data] } {
						if { [lindex $data $k] == "O21" } {
							set x1 [lindex $data [expr { $k + 3 }]]
							set sx1 [string length $x1]

							if { $sx1 > 8 } {
								set t 0
								while { [string range $x1 $t $t] != "." } {
									incr t
								}
								set xori [string range $x1 0 [expr { $t + 3 }]]
								set yori [string range $x1 [expr { $t + 4 }] end]
								set zori [lindex $data [expr { $k + 4 }]]
							} else { 
								set xori $x1
								set y1 [lindex $data [expr { $k + 4 }]]
								set sy1 [string length $y1]
								if { $sy1 > 8 } {
									set t 0
									while { [string range $y1 $t $t] != "." } {
										incr t
									}
									set yori [string range $y1 0 [expr { $t + 3 }]]
							 	 	set zori [string range $y1 [expr { $t + 4 }] end]
								}	else {
									set yori [lindex $data [expr { $k + 4 }]]
									set zori [lindex $data [expr { $k + 5 }]]
								}
							}
							if { $xori < $xmin } {
								set xmin $xori
							} 
							if { $yori < $ymin } {
								set ymin $yori
							}
							if { $zori > $zmax } {
								set zmax $zori
							}
						}
						incr k
					}

					mod1 $opt5 $xmin $ymin $zmax "lipids_no.pdb"
				}

		
				# FORMING THE INPUT FILE FOR THE BUILDER
				set count 0
				if { $opt3 > $opt2 } {
					set num $opt3
					for {set i 0} {$i < $num} {incr i} {
						set tmp_count 0
						set in1($i) $lipll($i)
						set k 0
						if { $lipul(0) != -1 } {
							while { $lipul($k) != -1 } {
								if { $lipul($k) == $lipll($i) } {
									set in4($i) $opt2n($k)
									incr count
									incr tmp_count
								} elseif { $tmp_count == 0 } { 
									set in4($i) 0
								}
								incr k
							}
						} else { 
							set in4($i) 0
						}
						set in2($i) $opt3H($i)
						set in3($i) $linkll($i)
						set in5($i) $opt3n($i)
					}
					set numl $i
					if { $count != $opt2 } {
						set k 0
						while { $k < $opt2 } {
							set dum $lipul($k)
							set count 0
							for {set i 0} {$i < $num } {incr i} {
								if { $lipll($i) == $dum } {
									incr count
								}
							}
							if { $count == 0 } {
								set in1($numl) $dum
								set in2($numl) $opt2H($k)
								set in3($numl) $linkul($k)
								set in4($numl) $opt2n($k)
								set in5($numl) 0
								incr numl
							}
							incr k
						}
					}		
				} else { 
					set num $opt2
					for {set i 0} {$i < $num} {incr i} {
						set tmp_count 0
						set in1($i) $lipul($i)
						set k 0
						if { $lipll(0) != -1 } {
							while { $lipll($k) != -1 } {
								if { $lipll($k) == $lipul($i) } {
									set in5($i) $opt3n($k)
									incr count
									incr tmp_count
								} elseif { $tmp_count == 0 } { 
									set in5($i) 0
								}
								incr k
							}
						} else { 
							set in5($i) 0
						}
							set in2($i) $opt2H($i)
							set in3($i) $linkul($i)
							set in4($i) $opt2n($i)
					}
					set numl $i
					if { $count != $opt3 } {
						set k 0
						while { $k < $opt3 } {
							set dum $lipll($k)
							set count 0
							for {set i 0} {$i < $num } {incr i} {
								if { $lipul($i) == $dum } {
									incr count
								}
							}
							if { $count == 0 } {
								set in1($numl) $dum
								set in2($numl) $opt2H($k)
								set in3($numl) $linkul($k)
								set in4($numl) 0
								set in5($numl) $opt2n($k)
								incr numl
							}
							incr k
						}
					}		
				}

				set num $numl
				set f [open "input" "w"]

				set max_apl 0.0
				set tnlul 0.0
				set tnlll 0.0
				set in6z 0.0
				for {set i 0} {$i < $num} {incr i} {
					set k 0
					while { $in1($i) != [lindex $p1 $k] } {
						incr k
					}
					if { [lindex $p7 $k] > $in6z } {
						set in6z [lindex $p7 $k]
					}
					set apl [lindex $p6 $k]
					if { $apl > $max_apl } {
						set max_apl $apl
						set ml_id $k
					}
					set tnlul [expr { $in4($i) + $tnlul }]
					set tnlll [expr { $in5($i) + $tnlll }]
				}
				set area_ul [expr { $max_apl * $tnlul }]
				set area_ll [expr { $max_apl * $tnlll }]

				if { $area_ul > $area_ll } {
					set apl $area_ul
				} else { 
					set apl $area_ll
				}

				if { $num > 1 } {
					# 20 PECENT INCREASE IN AREA PER LIPID FOR ASYMMETRIC LIPID BILAYERS
					set apl [expr { $apl * 1.2 }]
				}

				set in6xy [expr { sqrt($apl) }]	

				# INPUT1
				set inp_lip ""
				puts $f "\{ LIPIDS_PDB:" 
				for {set i 0} {$i < $num} {incr i} {
					set inp_lip [linsert $inp_lip $i $in1($i)]
					puts $f "$in1($i)_A.pdb"
				}
				puts $f "\}"

				set inp_lip [linsert $inp_lip $i "NA"]
		
				input_PDB $inp_lip

				# INPUT2
				set k 0
				puts $f "\{ LIPIDS_PDB_NAME:" 
				for {set i 0} {$i < $num} {incr i} {
					set count 0
					for {set j 0} {$j < $i} {incr j} {
						if { $in2($i) == $in2($j) } { 
							incr count
						}
					}
					if { $count == 0 } { 
						puts $f "$in2($i)"
					} else { 
						puts $f "$in2($i)$k"
						incr k
					}
				}
				puts $f "\}"

				# INPUT3
				puts $f "\{ HEAD_TAIL_LINKAGE:" 
				for {set i 0} {$i < $num} {incr i} {
					puts $f "$in3($i)"
				}
				puts $f "\}"

				# INPUT4
				puts $f "\{ LIPID_COMPOSITION_UL:" 
				for {set i 0} {$i < $num} {incr i} {
					puts $f "$in4($i)"
				}
				puts $f "\}"

				# INPUT5
				puts $f "\{ LIPID_COMPOSITION_LL:" 
				for {set i 0} {$i < $num} {incr i} {
					puts $f "$in5($i)"
				}
				puts $f "\}"

				# INPUT6
				puts $f "\{ BOX_DIMENSIONS: $in6xy $in6xy $in6z \}" 

				# INPUT7
				puts $f "\{ HOLE_DIMENSIONS: \}"

				# INPUT8
				if { $opt1 == 5 || $opt1 == 6 } {
					execution_pi [expr { $in6xy / 2.0 }] [expr { $in6xy / 2.0 }]
					file_delete
					files
					puts ""
					puts "***** CHECK CAREFULLY THE FILES GENERATED BEFORE CONTINUING *****"
					puts ""
					set opt4 "ref1_no.pdb"
					exec ls $opt4
				} elseif { $opt1 != 0 } {
					puts "# ENTER THE NAME OF THE PDB CONTATING THE MOLECULE"
					set opt4 [gets stdin]
					exec ls $opt4 
				} elseif { $opt1 == 0 } {
					set opt4 ""
				}
				puts $f "\{ PROTEIN_PDB: $opt4 \}"

				# INPUT9
				if { $opt1 == 1 || $opt1 == 3 } {
					puts ""
					puts "	### PUTING THE MOLECULE AT THE CENTRE OF THE BILAYER ### "
					puts ""
					set xtrans [expr { (-1 * $in6xy) / 2.0 }]
					set ytrans [expr { (-1 * $in6xy) / 2.0 }]
					set ztrans [expr { $in6z }]
					puts $f "\{ PROTEIN_SHIFT: $xtrans $ytrans $ztrans 0.0 0.0 0.0 \}"
				} 

				if { $opt1 == 2 || $opt1 == 4 } {
					puts ""
					puts "	### HOW OUT YOU WANT THE MOLECULE FROM THE TOP LEAFLET? ###"
					set out [gets stdin]
					puts ""
					set xtrans [expr { (-1 * $in6xy) / 2.0 }]
					set ytrans [expr { (-1 * $in6xy) / 2.0 }]
					set ztrans [expr { -5.0 - $out }]
					puts $f "\{ PROTEIN_SHIFT: $xtrans $ytrans $ztrans 0.0 0.0 0.0 \}"
				}
				 
				if { $opt1 == 0 || $opt1 == 5 || $opt1 == 6 } {
					puts $f "\{ PROTEIN_SHIFT: 0.0 0.0 0.0 0.0 0.0 0.0 \}"
				} 

				# INPUT10
				puts $f "\{ LIPID_BILAYER_METHOD: 1 36.0 \}"

				# INPUT11
				puts $f "\{ HOLE_METHOD: 1 \}"

				# INPUT12
				puts $f "\{ LIPID14_DIVISION: \}"

				# INPUT13
				puts $f "\{ CODE_EXECUTION: $opt1 \}"

				# INPUT14
				puts ""
				puts	"				WHICH VERSION DO YOU WANT TO EXECUTE? (VERSION 3.0 IS THE LATEST ONE) "
				set ve [gets stdin]
				puts ""

				puts $f "\{ VERSION: $ve" 
				for {set i 0} {$i < $num} {incr i} {
					puts $f "$in1($i)"
				}
				puts $f "$num \}"

				close $f
			}
		
			inp_grid

			puts "				#### FILE 'input' CONTAINS THE INPUT PARAMETERS FOR THE MEMBRANE BUILDER ####"

			set f [open "input" "r"]
			set data [read $f]
			close $f

			# CONTROL TO EXECUTE VARIOUS PART OF THE CODE

			spread
			if { [lindex $data 12 1] == 0 } {
				puts "				##### ONLY LIPID BILAYER WILL BE FORMED #####"
				puts ""	
				puts ""
				puts "				#### DO YOU WANT TO FORM A NON-RANDOM GRID BASED ON SIMPLE == AND != RULES? (Y/N)####"
				set gridinp [gets stdin]
				puts ""
				if { $gridinp == "Y" || $gridinp == "y" } {
					non_random_grid
				} else {
					lipid_bilayer
				}				
			} elseif { [lindex $data 12 1] == 3 || [lindex $data 12 1] == 4 || [lindex $data 12 1] == 6 } {
				puts "				##### INSERSION PART OF THE CODE WILL BE EXECUTED #####"
				if { [lindex $data 10 1] == 0 } {
					protein_profile
					lip_num
					lipid_hole
					protein_insersion
					lip_comp lip_pro.pdb
				} else {
					protein_profile
					lip_num
					protein_insersion
					replace
					lip_comp 1_final_struct.pdb
				}
			} elseif { [lindex $data 12 1] == 1 || [lindex $data 12 1] == 2 || [lindex $data 12 1] == 5} {
				puts "				##### BOTH LIPID BILAYER AND INSERSION PART OF THE CODE WILL BE EXECUTED #####"
				puts "				#### DO YOU WANT TO FORM A NON-RANDOM GRID BASED ON SIMPLE == AND != RULES? (Y/N) ####"
				set gridinp [gets stdin]
				puts ""
				if { $gridinp == "Y" || $gridinp == "y" } {
					non_random_grid
				} else {
					lipid_bilayer
				}				
				if { [lindex $data 10 1] == 0 } {
					protein_profile
					lip_num
					lipid_hole
					protein_insersion
					lip_comp lip_pro.pdb
				} else {
					protein_profile
					lip_num
					protein_insersion
					replace
					lip_comp 1_final_struct.pdb
				} 
			} elseif { [lindex $data 12 1] == 7 } {
				puts "				#### DETERMING THE PROTEIN PROFILE AND THE NUMBER OF LIPIDS TO BE INCREASED ####"
				protein_profile
				lip_num
			}
		} else { 
			set f [open "input" "r"]
			set data [read $f]
			close $f
		}

		puts ""
		puts "				##### DO YOU WANT TO SOLVATE THE SYSTEM? (REQUIRES AMBERTOOLS) (Y/N) #####"
		set sa [gets stdin]
		if { $sa == "Y" || $sa == "y" } {	
			puts ""
			puts "				##### DOES THE SYSTEM CONTAINS MORE THAN 100,000 RESIDUES (EXCLUDING WATER)? (Y/N) #####"
			set sol_ans [gets stdin]

			if { $sol_ans == "n" || $sol_ans == "N" } {
				if { [lindex $data 12 1] == 0 } {
					puts "				#### DO YOU WANT TO STRICTLY CONTROL THE NUMBER OF WATER MOLECULES? (NOTE: ONLY TO BE USED FOR PURE BILAYER SYSTEMS) (Y/N) ####"
					set sans [gets stdin]
					if { $sans == "Y" || $sans == "y" } {
						solvation_nc
					} else {
						solvation [lindex $data 12 1]
					}
				} else {
					solvation [lindex $data 12 1]
				}
			} else {
				puts "				#### FOR MORE THAN 100,000 RESIDUES WE PREFER THE APPROACH OF NORMALLY SOLVATING THE SYSTEM USING LEAP AND THEN USING CPPTRAJ TO SELECTIVELY REMOVE THE WATER FROM INSIDE"
				puts "				     AND ALONG THE BILAYER THICKNESS ####"
				puts ""
				puts "				#### THE CODE WILL GIVE YOU THE UNWANTED RESIDUE IDS OF THE WATER, USE MANUALLY CPPTRAJ TO STRIP THEM FROM THE PRMTOP AND INPCRD ####"
				puts ""
				solvation_100000 [lindex $data 12 1]
				puts "				#### DONE :: 1wat_resid_sorted CONTAINS THE RESID'S OF THE UNWATED WATER MOLECULES"
				puts "				     STRIP THEM FROM mol1.inpcrd AND mol1.prmtop USING CPPTRAJ ####"
			}
		}
	}
}

#################################################################### VESCICLE BUILDER ##############################################################################################################################

proc execution {} {

	puts "	************************************************"
	puts "	THIS IS AMBER BASED VESICLE BUILDER"
	puts "	************************************************"

	puts ""
	puts ""
	puts "	DEVELOPED BY TARUN KHANNA AND DR. IAN GOULD"
	puts "	IMPERIAL COLLEGE LONDON, U.K."

	# PARAMETERS 

	set gridsize 5.5
	set max_layer 20.0
	set gridsize_w 1.0

	# EMPTY WATER PORE FILE 

	set em [open "water_pore.pdb" "w"]
	close $em


	# VARIOUS LIPIDS AND THEIR HEAD AND TAIL
	set p1 { DLPC DMPC DOPC DPPC POPC DLPE DMPE DOPE DPPE POPE DOPG DOPS } 
	
	# LINKER ATOMS
	set p2 { O21 O21 O11 O21 O21 O21 O21 O21 O21 O21 O11 O11 }

	# AREA PER LIPID
	set p3 { 67.24 67.1 72.7 67.0 71.3 64.8 63.9 67.4 63.0 62.8 70.4 74.6 }

	# Z DIMENSIONS
	set p4 { 15.0 18.0 19.0 18.0 19.0 15.0 18.0 19.0 18.0 19.0 19.0 19.0 }

	# HEAD
	set p5 { PC PC PC PC PC PE PE PE PE PE PG PS}

	# LIPID INPUT
	
	puts ""
	puts "				#### NOTE: THIS VERSION ONLY SUPPORTS THE BUILDING OF THE SYMMETRIC LIPID VESICLES ####"
	puts ""

	puts "				#### ENTER THE NAME OF THE LIPID YOU WANT TO FORM VESICLE FROM ####"
	set lip [gets stdin]
	input_PDB [list $lip "NA"]
	puts ""

	set k 0

	set count 0
	while { $k < [llength $p1] } {
		if { [lindex $p1 $k] == $lip } {
			set ln $k
			incr count
		}
		incr k
	}

	if { $count == 0 } {	
		puts "				#### ERROR :: THE LIPID NOT FOUND IN THE DATABASE ####"
	}

	set lt [lindex $p4 $ln]

	set liptype [list [lindex $p5 $ln]]

	set patom1 [lindex $p2 $ln]
	set patom2 "C112"
	set patom [list $patom1 $patom2]

	set apl [lindex $p3 $ln]

	puts "				#### ENTER THE RADIUS OF THE INNER LEAFLET IN ANGSTOM ####"
	set R_LL [gets stdin]
	puts ""

	puts "				#### DO YOU WANT TO FORM VESICLE? #### (y/n) ####"
	set inp [gets stdin]
	
	set f [open "input" "w"]

	puts $f "\{ LIPIDS_PDB: $lip\_A.pdb \}"
	puts $f "\{ LIPID_PDB_NAME: [lindex $p5 $ln] \}"
	puts $f "\{ LIPID_ELEM_OLD: $patom1 \}"
	puts $f "\{ LIPID_ELEM_NEW: $patom1 \}"
	puts $f "\{ PROTEIN_PDB: lipids.pdb \}"
	puts $f "\{ LIPID_TAILS_OLD: $patom2 \}"
	puts $f "\{ LIPID_TAILS_NEW: $patom2 \}"

	close $f

	if { $inp == "Y" || $inp == "y" } { 
		puts ""
		puts "				#### DO YOU WANT TO PUT WATER PORES FOR FAST EQUILIBRATION BETWEEN INNER AND OUTER LEAFLET (NOT TESTED RIGOROUSLY) ? (y/n) ####"
		set inpwp [gets stdin]
		set R_LL_OLD $R_LL
		if { $inpwp == "Y" || $inpwp == "y" } {
			set R_UL [expr { $R_LL + (2*$lt) + (2.0) }]
			set rpore [expr { $R_UL * 0.1 }]
			set rpore [format "%.1f" $rpore]
			set R_LL [expr { ((4*$R_LL_OLD*$R_LL_OLD) + (6*$rpore * $rpore)) / 4 }]
			set R_LL [expr { sqrt($R_LL) }]
			set R_LL [format "%.1f" $R_LL]
			spherical_grid $gridsize $lt $apl $R_LL $liptype $patom $rpore
			set twat [water_pore $lt $R_LL $rpore [expr { $gridsize_w + 0.4 }]]
		} else {
			spherical_grid $gridsize $lt $apl $R_LL $liptype $patom NA
		}
		puts ""
		puts "				#### DO YOU WANT TO PUT THE SOLVATING WATER? (y/n) ####"
		set inp2 [gets stdin]
		if { $inp2 == "Y" || $inp2 == "y" } {
			if { $inpwp == "Y" || $inpwp == "y" } {
				water_grid $max_layer $gridsize_w $lt $R_LL $R_LL_OLD $rpore $twat
			} else {
				water_grid $max_layer $gridsize_w $lt $R_LL $R_LL_OLD NA 0
			}
			puts ""
			puts "				***** sol_lipid.pdb IS THE FINAL SOLVATED LIPID VESICLE *****"
			puts ""
			box_dimensions "sol_lipid.pdb"
			
			if { $inpwp == "Y" || $inpwp == "y" } {
				puts ""
				puts "		*******************************************************************************"
				puts "		#### IMP NOTE ::: $twat NUMBER OF WATERS AFTER THE LIPIDS NEED TO BE FIXED ####"
				puts "		*******************************************************************************"
				puts ""
			}
		} else {
			puts "				***** lipids_wo.pdb IS THE FINAL DRY VESICLE *****"
			box_dimensions "lipids_wo.pdb"
		}
	} else {
		puts ""
		puts "				#### DO YOU WANT TO FORM DIBS? #### (y/n) ####"
		set inp [gets stdin]
		if { $inp == "Y" || $inp == "y" } {
			spherical_grid_dibs $gridsize [expr { $lt - 4.0 }] $apl $R_LL $liptype $patom
			puts ""
			puts "				#### DO YOU WANT TO PUT THE SOLVATING WATER? (y/n) ####"
			set inp2 [gets stdin]
			if { $inp2 == "Y" || $inp2 == "y" } {
				water_grid_dibs $max_layer $gridsize_w [expr { $lt - 4.0 }] $R_LL
				puts ""
				puts "				***** sol_lipid.pdb IS THE FINAL SOLVATED DIBs *****"
				box_dimensions "sol_lipid.pdb"
			} else {
				puts "				***** lipids_wo.pdb IS THE FINAL DRY DIBs *****"
				box_dimensions "lipids_wo.pdb"
			} 
		} 
	}

	set bd [open "box_dimensions" "r"]
	set datbd [read $bd]
	close $bd

	puts ""
	puts "			#### IMPORTANT NOTE : BOX DIMENSIONS OF THE SYSTEM IS [lindex $datbd 0] [lindex $datbd 1] [lindex $datbd 2] ####"
	puts ""

}

proc spherical_grid_dibs { gridsize lt apl R_LL liptype patom } {

	# FORMING A CONCENTRIC SPHERICAL GRID

	global nlipid_OL
	global nlipid_IL

	set h [open "spherical_coord" "w"]
	set h1 [open "spherical_coord_ll" "w"]

	set R_UL [expr { $R_LL + (2*$lt) + (2.0) }]

	set ndl 1

	# OUR CRITERION

	set nlipid [expr { (4 * 3.14 * $R_UL * $R_UL) / $apl }]
	set npoints [expr { sqrt($nlipid) }]
	set npoints [format "%.0f" $npoints]
	set nlipid_OL [expr { $npoints * $npoints }]

	set nlipid [expr { (4 * 3.14 * $R_LL * $R_LL) / $apl }]
	set npoints [expr { sqrt($nlipid) }]
	set npoints [format "%.0f" $npoints]
	set nlipid_IL [expr { $npoints * $npoints }]

	# CHARMM CRITERION

	#set a 0.205
	#set b 30.38
	#set c 1569.8

	#set nlipid [expr { ($a* $R_UL * $R_UL) + ($b*$R_UL) - ($c) }]
	#set npoints [expr { sqrt($nlipid) }]
	#set npoints [format "%.0f" $npoints]
	#set nlipid_OL [expr { $npoints * $npoints }]

	#set nlipid [expr { ($a* $R_LL * $R_LL) + ($b*$R_LL) - ($c) }]
	#set npoints [expr { sqrt($nlipid) }]
	#set npoints [format "%.0f" $npoints]
	#set nlipid_IL [expr { $npoints * $npoints }]

	puts "			#### PUTTING IN $nlipid_IL IN THE FORM OF DIBS OF $R_LL Ang RADIUS####"
	puts ""

	#set Ustep [expr { 1.0 / $npoints }]
	#set Ustep [format "%.2f" $Ustep]
	#set Vstep [expr { 1.0 / $npoints }]
	#set Vstep [format "%.2f" $Vstep]

	# LEAFLET 1

	puts "				  LEAFLET 1"
	puts "				*************"
	puts ""

	set k 0
	while { $k < $nlipid_IL } {
		set u [expr { rand() }]
		set v [expr { rand() }]
		set u [format "%.2f" $u]
		set v [format "%.2f" $v]
		set phi [expr { acos((2*$v)-1) }]
		set theta [expr { 2 * 3.14 * $u }]
		set x($k) [expr { sin($phi) * cos($theta) }]
		set y($k) [expr { sin($phi) * sin($theta) }]
		set z($k) [expr { cos($phi) }]
		set count 0
		for {set i 0} {$i < $k} {incr i} {
			set xc [expr { $x($k) - $x($i) }]
			set yc [expr { $y($k) - $y($i) }]
			set zc [expr { $z($k) - $z($i) }]			
			set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
			set dis [expr { $dis * $R_LL * $R_LL }]
			if { $dis < [expr { $gridsize * $gridsize }] } {
				incr count
			}
		}
		if { $count == 0 } {
			#puts "				#### PUTTING IN LIPID [expr { $k + 1 }] ####"	
			puts $h1 "[expr { $R_LL * $x($k) }] [expr { $R_LL * $y($k) }] [expr { $R_LL * $z($k) }]"
			puts $h1 "[expr { ($R_LL + $lt) * $x($k) }] [expr { ($R_LL + $lt) * $y($k) }] [expr { ($R_LL + $lt) * $z($k) }]"
			incr k
		}
	}

	# LEAFLET 2

	puts "				  LEAFLET 2"
	puts "				*************"
	puts ""

	set k 0
	while { $k < $nlipid_IL } {
		set u [expr { rand() }]
		set v [expr { rand() }]
		set u [format "%.2f" $u]
		set v [format "%.2f" $v]
		set phi [expr { acos((2*$v)-1) }]
		set theta [expr { 2 * 3.14 * $u }]
		set x($k) [expr { sin($phi) * cos($theta) }]
		set y($k) [expr { sin($phi) * sin($theta) }]
		set z($k) [expr { cos($phi) }]
		set count 0
		for {set i 0} {$i < $k} {incr i} {
			set xc [expr { $x($k) - $x($i) }]
			set yc [expr { $y($k) - $y($i) }]
			set zc [expr { $z($k) - $z($i) }]			
			set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
			set dis [expr { $dis * $R_LL * $R_LL }]
			if { $dis < [expr { $gridsize * $gridsize }] } {
				incr count
			}
		}
		if { $count == 0 } {
			#puts "				#### PUTTING IN LIPID [expr { $k + 1 }] ####"	
			# SHIFTING X BY 2*RIL + lt
			set x1 [expr { $R_LL * $x($k) }]
			set x1 [expr { $x1 + (2*$R_LL) + (2*$lt) }]
			set x2 [expr { (($R_LL  + $lt) * $x($k)) }]
			set x2 [expr { $x2 + (2*$R_LL) + (2*$lt) }]
			puts $h "$x1 [expr { $R_LL * $y($k) }] [expr { $R_LL * $z($k) }]"
			puts $h "$x2	[expr { ($R_LL + $lt) * $y($k) }] [expr { ($R_LL + $lt) * $z($k) }]"
			incr k
		}
	}
	close $h
	close $h1
	
	vesicles $nlipid_IL $nlipid_IL $liptype $patom
}		

proc spherical_grid { gridsize lt apl R_LL liptype patom rpore } {

	# FORMING A CONCENTRIC SPHERICAL GRID

	global nlipid_OL
	global nlipid_IL

	set h [open "spherical_coord" "w"]
	set h1 [open "spherical_coord_ll" "w"]

	set R_UL [expr { $R_LL + (2*$lt) + (2.0) }]

	# CHECKING THE WATER PORE REQUIREMENT

	if { $rpore != "NA" } {
		set pparm 1
	} else { 
		set pparm 0
	}

	set ndl 1

	# OUR CRITERION

	set nlipid [expr { (4 * 3.14 * $R_UL * $R_UL) / $apl }]
	set npoints [expr { sqrt($nlipid) }]
	set npoints [format "%.0f" $npoints]
	set nlipid_OL [expr { $npoints * $npoints }]

	set nlipid [expr { (4 * 3.14 * $R_LL * $R_LL) / $apl }]
	set npoints [expr { sqrt($nlipid) }]
	set npoints [format "%.0f" $npoints]
	set nlipid_IL [expr { $npoints * $npoints }]

	# CHARMM CRITERION

	#set a 0.205
	#set b 30.38
	#set c 1569.8

	#set nlipid [expr { ($a* $R_UL * $R_UL) + ($b*$R_UL) - ($c) }]
	#set npoints [expr { sqrt($nlipid) }]
	#set npoints [format "%.0f" $npoints]
	#set nlipid_OL [expr { $npoints * $npoints }]

	#set nlipid [expr { ($a* $R_LL * $R_LL) + ($b*$R_LL) - ($c) }]
	#set npoints [expr { sqrt($nlipid) }]
	#set npoints [format "%.0f" $npoints]
	#set nlipid_IL [expr { $npoints * $npoints }]

	puts "			#### PUTTING IN $nlipid_OL ON THE OUTER LEAFLET AND $nlipid_IL ON THE INNER LEAFLET IN A LIPID VESICLE STRUCTURE WITH A RADIUS OF $R_UL Ang ####"
	puts ""

	#set Ustep [expr { 1.0 / $npoints }]
	#set Ustep [format "%.2f" $Ustep]
	#set Vstep [expr { 1.0 / $npoints }]
	#set Vstep [format "%.2f" $Vstep]

	# INNER LAYER

	puts "				INNER LEAFLET"
	puts "				*************"
	puts ""

	set k 0
	while { $k < $nlipid_IL } {
		set u [expr { rand() }]
		set v [expr { rand() }]
		set u [format "%.2f" $u]
		set v [format "%.2f" $v]
		set phi [expr { acos((2*$v)-1) }]
		set theta [expr { 2 * 3.14 * $u }]
		set x($k) [expr { sin($phi) * cos($theta) }]
		set y($k) [expr { sin($phi) * sin($theta) }]
		set z($k) [expr { cos($phi) }]

		set count 0
		# CHECK PORE THE WATER PORE 
		if { $pparm == 1 } {
			set xc1 [expr { ($x($k) - 1.0) * ($x($k) - 1.0) }]
			set xc2 [expr { ($x($k) + 1.0) * ($x($k) + 1.0) }]
			set xc3 [expr { $x($k) * $x($k) }]
			set yc1 [expr { ($y($k) - 1.0) * ($y($k) - 1.0) }]
			set yc2 [expr { ($y($k) + 1.0) * ($y($k) + 1.0) }]
			set yc3 [expr { $y($k) * $y($k) }]
			set zc1 [expr { ($z($k) - 1.0) * ($z($k) - 1.0) }]
			set zc2 [expr { ($z($k) + 1.0) * ($z($k) + 1.0) }]
			set zc3 [expr { $z($k) * $z($k) }]

			set check1 [expr { $xc1 + $yc3 + $zc3 }]
			set check2 [expr { $xc2 + $yc3 + $zc3 }]
			set check3 [expr { $yc1 + $zc3 + $xc3 }]
			set check4 [expr { $yc2 + $zc3 + $xc3 }]
			set check5 [expr { $zc1 + $yc3 + $xc3 }]
			set check6 [expr { $zc2 + $yc3 + $xc3 }]
			set rsd [expr { $rpore / $R_LL }]
			set rsd [expr { $rsd * $rsd }]

			if { $check1 <= $rsd || $check2 <= $rsd || $check3 <= $rsd || $check4 <= $rsd || $check5 <= $rsd || $check6 <= $rsd } {
				incr count
			}
		} 

		for {set i 0} {$i < $k} {incr i} {
			set xc [expr { $x($k) - $x($i) }]
			set yc [expr { $y($k) - $y($i) }]
			set zc [expr { $z($k) - $z($i) }]			
			set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
			set dis [expr { $dis * $R_LL * $R_LL }]
			if { $dis < [expr { $gridsize * $gridsize }] } {
				incr count
				set i $k
			}
		}
		if { $count == 0 } {
			#puts "				#### PUTTING IN LIPID [expr { $k + 1 }] ####"	
			puts $h1 "[expr { $R_LL * $x($k) }] [expr { $R_LL * $y($k) }] [expr { $R_LL * $z($k) }]"
			puts $h1 "[expr { ($R_LL + $lt) * $x($k) }] [expr { ($R_LL + $lt) * $y($k) }] [expr { ($R_LL + $lt) * $z($k) }]"
			incr k
		}
	}

	# OUTER LAYER

	puts "				*************"
	puts "				OUTER LEAFLET"
	puts "				*************"
	puts ""

	set k 0
	while { $k < $nlipid_OL } {
		set u [expr { rand() }]
		set v [expr { rand() }]
		set u [format "%.2f" $u]
		set v [format "%.2f" $v]
		set phi [expr { acos((2*$v)-1) }]
		set theta [expr { 2 * 3.14 * $u }]
		set x($k) [expr { sin($phi) * cos($theta) }]
		set y($k) [expr { sin($phi) * sin($theta) }]
		set z($k) [expr { cos($phi) }]
		set count 0
		# CHECK PORE THE WATER PORE 
		
		if { $pparm == 1 } {
			set xc1 [expr { ($x($k) - 1.0) * ($x($k) - 1.0) }]
			set xc2 [expr { ($x($k) + 1.0) * ($x($k) + 1.0) }]
			set xc3 [expr { $x($k) * $x($k) }]
			set yc1 [expr { ($y($k) - 1.0) * ($y($k) - 1.0) }]
			set yc2 [expr { ($y($k) + 1.0) * ($y($k) + 1.0) }]
			set yc3 [expr { $y($k) * $y($k) }]
			set zc1 [expr { ($z($k) - 1.0) * ($z($k) - 1.0) }]
			set zc2 [expr { ($z($k) + 1.0) * ($z($k) + 1.0) }]
			set zc3 [expr { $z($k) * $z($k) }]

			set check1 [expr { $xc1 + $yc3 + $zc3 }]
			set check2 [expr { $xc2 + $yc3 + $zc3 }]
			set check3 [expr { $yc1 + $zc3 + $xc3 }]
			set check4 [expr { $yc2 + $zc3 + $xc3 }]
			set check5 [expr { $zc1 + $yc3 + $xc3 }]
			set check6 [expr { $zc2 + $yc3 + $xc3 }]
			set rsd [expr { $rpore / $R_UL }]
			set rsd [expr { $rsd * $rsd }]

			if { $check1 <= $rsd || $check2 <= $rsd || $check3 <= $rsd || $check4 <= $rsd || $check5 <= $rsd || $check6 <= $rsd } {
				incr count
			}
		} 

		for {set i 0} {$i < $k} {incr i} {
			set xc [expr { $x($k) - $x($i) }]
			set yc [expr { $y($k) - $y($i) }]
			set zc [expr { $z($k) - $z($i) }]			
			set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
			set dis [expr { $dis * $R_UL * $R_UL }]
			if { $dis < [expr { $gridsize * $gridsize }] } {
				incr count
				set i $k
			}
		}
		if { $count == 0 } {	
			#puts "				#### PUTTING IN LIPID [expr { $k + 1 }] ####"	
			puts $h "[expr { $R_UL * $x($k) }] [expr { $R_UL * $y($k) }] [expr { $R_UL * $z($k) }]"
			puts $h "[expr { ($R_UL - $lt) * $x($k) }] [expr { ($R_UL - $lt) * $y($k) }] [expr { ($R_UL - $lt) * $z($k) }]"
			incr k
		}
	}
	close $h
	close $h1
	
	vesicles $nlipid_IL $nlipid_OL $liptype $patom
}		

proc vesicles { nlipid_IL nlipid_OL liptype patom } {

# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "	
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set h [open "spherical_coord" "r"]
	set data [read $h]
	close $h

	set h1 [open "spherical_coord_ll" "r"]
	set data1 [read $h1]
	close $h1

	set f [open "lipids.pdb" "w"]

	set ndl [llength $liptype]

	# PUTTING THE LIPIIDS IN THE OUTER LEAFLET
	
	set k 0
	set j 1
	for {set i 0} {$i < $nlipid_OL} {incr i} {
		set rn1 $j
		set srn1 [string length $rn1]
		incr j

		set rn11 $j
		set srn11 [string length $rn11]
		incr j

		set at1 [lindex $patom 0]
		set sat1 [string length $at1]

		set at11 [lindex $patom 1]
		set sat11 [string length $at11]

		set an1 [expr { $i + 1 }] 
		set san1 [string length $an1]

		set resname [lindex $liptype 0]
		set sresname [string length $resname]

		set x1 [lindex $data $k]
		set x1 [format "%.3f" $x1]
		set sx1 [string length $x1]

		set y1 [lindex $data [expr { $k + 1 }]]
		set y1 [format "%.3f" $y1]
		set sy1 [string length $y1]

		set z1 [lindex $data [expr { $k + 2 }]]
		set z1 [format "%.3f" $z1]
		set sz1 [string length $z1]
		
		set x2 [lindex $data [expr { $k + 3 }]]
		set x2 [format "%.3f" $x2]
		set sx2 [string length $x2]

		set y2 [lindex $data [expr { $k + 4 }]]
		set y2 [format "%.3f" $y2]
		set sy2 [string length $y2]

		set z2 [lindex $data [expr { $k + 5 }]]
		set z2 [format "%.3f" $z2]
		set sz2 [string length $z2]

		incr k 6
		
		puts $f "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname  a$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
		puts $f "HETATM$p1($srn11)$rn11 $ic($sat11)$at11$sat($sat1)$resname  a$p($san1)$an1    $c($sx2)$x2$c($sy2)$y2$c($sz2)$z2  1.00  0.00           C"
		puts $f "TER"
	}

	# PUTTING THE LIPIIDS IN THE INNER LEAFLET
	
	set k 0
	for {set i $nlipid_OL} {$i < [expr { $nlipid_IL + $nlipid_OL }]} {incr i} {
		set rn1 $j
		set srn1 [string length $rn1]
		incr j

		set rn11 $j
		set srn11 [string length $rn11]
		incr j

		set at1 [lindex $patom 0]
		set sat1 [string length $at1]

		set at11 [lindex $patom 1]
		set sat11 [string length $at11]

		set an1 [expr { $i + 1 }] 
		set san1 [string length $an1]

		set resname [lindex $liptype 0]
		set sresname [string length $resname]

		set x1 [lindex $data1 $k]
		set x1 [format "%.3f" $x1]
		set sx1 [string length $x1]

		set y1 [lindex $data1 [expr { $k + 1 }]]
		set y1 [format "%.3f" $y1]
		set sy1 [string length $y1]

		set z1 [lindex $data1 [expr { $k + 2 }]]
		set z1 [format "%.3f" $z1]
		set sz1 [string length $z1]
		
		set x2 [lindex $data1 [expr { $k + 3 }]]
		set x2 [format "%.3f" $x2]
		set sx2 [string length $x2]

		set y2 [lindex $data1 [expr { $k + 4 }]]
		set y2 [format "%.3f" $y2]
		set sy2 [string length $y2]

		set z2 [lindex $data1 [expr { $k + 5 }]]
		set z2 [format "%.3f" $z2]
		set sz2 [string length $z2]

		incr k 6
		
		puts $f "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname  a$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
		puts $f "HETATM$p1($srn11)$rn11 $ic($sat11)$at11$sat($sat1)$resname  a$p($san1)$an1    $c($sx2)$x2$c($sy2)$y2$c($sz2)$z2  1.00  0.00           C"
		puts $f "TER"
	}
	close $f

	lipid_insert
}

proc lipid_insert {} {

	set patom [list "O21" "C112"]

	puts ""
	puts "				#### INSERTING THE LIPIDS IN THE RIGHT ORIENTATION ####"

	package require math::linearalgebra

	# REMOVE AND THEN REPLACE

	package require math::linearalgebra

	set inp [open "input" "r"]
	set mn [read $inp]
	close $inp

	set pdb [lindex $mn 4 1]

	set f [open "$pdb" "r"]
	set data1 [read $f]
	close $f

	set h [open "lipids_wo.pdb" "w"]

	# space variables

	set old_res -99
	set new_res 99

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set k 0

	set next 0
	
	set start 1
	while { $k < [llength $data1] } {
		set term [lindex $data1 $k]
		set t1 [string range $term 0 5]
		if { [lindex $data1 $k] == "ATOM" || $t1 == "HETATM" } {
			if { $t1 == "HETATM" } {
				set sterm [string length $term]
				if { $sterm > 6 } {
					set shift 1
					set rn1 $term 
					set srn1 [string length $rn1]
					set rnum [string range $rn1 6 $srn1]
					set srn1 5
					set ft ""
				} else {
						set shift 0
						set rn1 [lindex $data1 [expr { $k + 1 }]]
						set rnum $rn1
						set srn1 [string length $rn1]
						set ft "HETATM"
				}
			} else { 
					set shift 0
					set rn1 [lindex $data1 [expr { $k + 1 }]]
					set srn1 [string length $rn1]
					set rnum $rn1
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

				set new_res $an1

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

				set t 0
				set count 0
				set countt 0

				while { $t < [llength [lindex $mn 1]] } {
					if { $resn == [lindex $mn 1 $t] && $at1 == [lindex $mn 2 $t]} {
						incr count
						set place $t
					}
					incr t
				}

				set t 0
				while { $t < [llength [lindex $mn 1]] } {
					if { $resn == [lindex $mn 1 $t] && $at1 == [lindex $mn 5 $t]} {
						incr countt
					}
					incr t
				}
		
				if { $count != 0 } {
					set headx $corx
					set heady $cory
					set headz $corz
					incr next
				}

				if { $countt != 0 } {	
					set tailx $corx
					set taily $cory
					set tailz $corz
					incr next
				}

				if { $next == 2 } {
					set next 0

					#set start $rnum
	
					#puts "			**** REPLACING THE LIPID ****"

					set npdb "[lindex $mn 0 $place]"
					set m [open "$npdb" "r"]
					set rpdb [read $m]
					close $m

					set k2 0

					while { $k2 < [llength $rpdb] } {
						if { [lindex $rpdb $k2] == [lindex $mn 3 $place] } {
							set xr [lindex $rpdb [expr { $k2 + 4 }]]
							set yr [lindex $rpdb [expr { $k2 + 5 }]]
							set zr [lindex $rpdb [expr { $k2 + 6 }]]
						}
						incr k2 
					}

					set xorigin [expr { $xr - $headx }]
					set yorigin [expr { $yr - $heady }]
					set zorigin [expr { $zr - $headz }]
			
					# GETTING THE RIGHT ORIENTATION OF THE LIPID :: SIMILAR TO THE LIPID BEING REPLACED

					set old_vec [list [expr { $headx - $tailx }] [expr { $heady - $taily }] [expr { $headz - $tailz }]]

					set k2 0 
	
					while { $k2 < [llength $rpdb] } {
						if { [lindex $rpdb $k2] == "HETATM" } {

							set at1 [lindex $rpdb [expr { $k2 + 2 }]]

							if { $at1 == [lindex $mn 3 $place] } {
								set nheadx [lindex $rpdb [expr { $k2 + 6 }]]
								set nheadx [format "%.3f" [expr { $nheadx - $xorigin }]]
							
								set nheady [lindex $rpdb [expr { $k2 + 7 }]] 
								set nheady [format "%.3f" [expr { $nheady - $yorigin }]]
								

								set nheadz [lindex $rpdb [expr { $k2 + 8 }]] 
								set nheadz [format "%.3f" [expr { $nheadz - $zorigin }]]
							}
								
							if { $at1 == [lindex $mn 6 $place] } {
								set ntailx [lindex $rpdb [expr { $k2 + 6 }]]
								set ntailx [format "%.3f" [expr { $ntailx - $xorigin }]]
							
								set ntaily [lindex $rpdb [expr { $k2 + 7 }]] 
								set ntaily [format "%.3f" [expr { $ntaily - $yorigin }]]
								

								set ntailz [lindex $rpdb [expr { $k2 + 8 }]] 
								set ntailz [format "%.3f" [expr { $ntailz - $zorigin }]]
							}
						}
						incr k2
					}
					set new_vec [list [expr { $nheadx - $ntailx }] [expr { $nheady - $ntaily }] [expr { $nheadz - $ntailz }]]
				
					set old_vec [::math::linearalgebra::unitLengthVector $old_vec]
					set new_vec [::math::linearalgebra::unitLengthVector $new_vec]
				
					set angle [::math::linearalgebra::dotproduct $old_vec $new_vec]
					set normal [::math::linearalgebra::crossproduct $old_vec $new_vec]
					set normal [::math::linearalgebra::unitLengthVector $normal]

					set dum1 [open "dummy1" "w"]
					puts $dum1 "$normal"
					close $dum1

					set dum1 [open "dummy1" "r"]
					set du [read $dum1]
					close $dum1

					# ROTATION OF THE TAIL ABOUT THE ORIGIN AND SO IS THE OTHER ATOMS IN THE RESIDUE

					set u [lindex $du 0]
					set v [lindex $du 1]
					set w [lindex $du 2]
					set u2 [expr { $u * $u }]
					set v2 [expr { $v * $v }]
					set w2 [expr { $w * $w }]
					set l [expr { $u2 + $v2 + $w2 }]

					set rl [expr { sqrt($l) }]
	
					#set theta [expr { ($angle * 3.14 * -1.0) / 180.0 }]
					set theta [expr { acos($angle) }]

					set a11 [expr { ($u2 + (($v2 + $w2) * cos($theta))) / $l } ]
					set a12 [expr { ((($u*$v)*(1-cos($theta))) - ($w * $rl * sin($theta))) / $l }]
					set a13 [expr { ((($u*$w)*(1-cos($theta))) + ($v * $rl * sin($theta))) / $l }]
					set a21 [expr { ((($u*$v)*(1-cos($theta))) + ($w * $rl * sin($theta))) / $l }]
					set a22 [expr { ($v2 + (($u2 + $w2) * cos($theta))) / $l }]
					set a23 [expr { ((($v*$w)*(1-cos($theta))) - ($u * $rl * sin($theta))) / $l }]
					set a31 [expr { ((($u*$w)*(1-cos($theta))) - ($v * $rl * sin($theta))) / $l }]
					set a32 [expr { ((($v*$w)*(1-cos($theta))) + ($u * $rl * sin($theta))) / $l }]
					set a33 [expr { ($w2 + (($u2 + $v2) * cos($theta))) / $l }]

					set row1 [list $a11 $a12 $a13]
					set row2 [list $a21 $a22 $a23]
					set row3 [list $a31 $a32 $a33]
					set tr [list $row1 $row2 $row3]

					#puts "$normal"
					#puts "$angle" 
						
					set k2 0

					while { $k2 < [llength $rpdb] } {
						if { [lindex $rpdb $k2] == "HETATM" } {
							if { $start == 99999 } {
								set start 1
							}
							set rn1 $start
							set srn1 [string length $rn1]

							set at1 [lindex $rpdb [expr { $k2 + 2 }]]
							set sat1 [string length $at1]

							#set an1 [lindex $rpdb [expr { $k2 + 5 }]]
							#set san1 [string length $an1]

							set x1 [lindex $rpdb [expr { $k2 + 6 }]]
							set x1 [expr { $x1 - $xorigin }]
							set x1 [expr { $nheadx - $x1 }]

							set y1 [lindex $rpdb [expr { $k2 + 7 }]] 
							set y1 [expr { $y1 - $yorigin }]
							set y1 [expr { $nheady - $y1 }]

							set z1 [lindex $rpdb [expr { $k2 + 8 }]] 
							set z1 [expr { $z1 - $zorigin }]
							set z1 [expr { $nheadz - $z1 }]
		
							#set R [expr { sqrt(($x1*$x1) + ($y1*$y1) + ($z1*$z1)) }]

							set coord [list $x1 $y1 $z1]
							#set coord [::math::linearalgebra::unitLengthVector $coord]
							set tcoord [::math::linearalgebra::matmul $coord $tr]

							set dum [open "dummy" "w"]
	
							puts $dum "$tcoord"

							close $dum

							set dum [open "dummy" "r"]
							set dvar [read $dum]
							close $dum

							set x1 [format "%.3f" [expr { $nheadx - [lindex $dvar 0] }]]
							#set x1 [format "%.3f" [lindex $dvar 0]]
							set sx1 [string length $x1]

							set y1 [format "%.3f" [expr { $nheady - [lindex $dvar 1] }]]
							#set y1 [format "%.3f" [lindex $dvar 1]]
							set sy1 [string length $y1]

							set z1 [format "%.3f" [expr { $nheadz - [lindex $dvar 2] }]]
							#set z1 [format "%.3f" [lindex $dvar 2]]
							set sz1 [string length $z1]

							set resn [lindex $rpdb [expr { $k2 + 3 }]]
							set sresn [string length [lindex $rpdb [expr { $k2 + 3 }]]]
							incr start
				
							puts $h "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$p($sresn)$resn $cn$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $rpdb [expr { $k2 + 11 }]]"
						}
				incr k2
				}
			puts $h "TER"
			set old_res $new_res
			}
		}
		incr k
	}
	#puts $h "END"
	close $h	

	#lipid_overlap_VES $patom
}

proc lipid_overlap_VES { patom } {

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	puts ""
	puts "				#### REMOVING THE OVERLAPS ####"
	puts ""

	# REMOVING THE OVERLAPS BY ROTATING THE LIPIDS ABOUT THE VECTOR CONNECTING THE PIVOT ATOM AND THE END OF THE TAIL

	global x1n
	global y1n
	global z1n
	global an1n
	global rn1n
	global at1n
	global nlipid_OL
	global nlipid_IL
	
	# NEIGHBOUR DEFINATION IS 10A FROM THE PIVOT OXYGEN ATOM
	
	set nei 12.0

	set pivotO [lindex $patom 0]

	set f [open "lipids_wo.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "lipids_no.pdb" "w"]

	# DETERMING THE NEIGHBOURS

	set h [open "neighbours" "w"]

	set k 0
	while { $k < [llength $data] } {
		if { [lindex $data $k] == "$pivotO" } {
			set x1 [lindex $data [expr { $k + 4 }]]
			set sx1 [string length $x1]

			if { $sx1 > 8 } {
				set t 0
				while { [string range $x1 $t $t] != "." } {
					incr t
				}
				set xori [string range $x1 0 [expr { $t + 3 }]]
				set yori [string range $x1 [expr { $t + 4 }] end]
				set syori [string length $yori]
				if { $syori > 8 } {
					set y1 $yori
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set yori [string range $y1 0 [expr { $t + 3 }]]
			 	 	set zori [string range $y1 [expr { $t + 4 }] end]
				} else {
					set zori [lindex $data [expr { $k + 5 }]]
				}
			} else { 
				set xori $x1
				set y1 [lindex $data [expr { $k + 5 }]]
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set yori [string range $y1 0 [expr { $t + 3 }]]
			 	 	set zori [string range $y1 [expr { $t + 4 }] end]
				}	else {
					set yori [lindex $data [expr { $k + 5 }]]
					set zori [lindex $data [expr { $k + 6 }]]
				}
			}
			set xr $xori
			set yr $yori
			set zr $zori

			set r2r [expr { ($xr*$xr) + ($yr*$yr) + ($zr*$zr) }]

			set rid [lindex $data [expr { $k + 3 }]]

			puts $h "\{ $rid"

			set k2 0

			while { $k2 < [llength $data] } {
				if { [lindex $data $k2] == "$pivotO" } {
					set x1 [lindex $data [expr { $k2 + 4 }]]
					set sx1 [string length $x1]

					if { $sx1 > 8 } {
						set t 0
						while { [string range $x1 $t $t] != "." } {
							incr t
						}
						set xori [string range $x1 0 [expr { $t + 3 }]]
						set yori [string range $x1 [expr { $t + 4 }] end]
						set syori [string length $yori]
						if { $syori > 8 } {
							set y1 $yori
							set t 0
							while { [string range $y1 $t $t] != "." } {
								incr t
							}
							set yori [string range $y1 0 [expr { $t + 3 }]]
					 	 	set zori [string range $y1 [expr { $t + 4 }] end]
						} else {
							set zori [lindex $data [expr { $k2 + 5 }]]
						}
					} else { 
						set xori $x1
						set y1 [lindex $data [expr { $k2 + 5 }]]
						set sy1 [string length $y1]
						if { $sy1 > 8 } {
							set t 0
							while { [string range $y1 $t $t] != "." } {
								incr t
							}
							set yori [string range $y1 0 [expr { $t + 3 }]]
					 	 	set zori [string range $y1 [expr { $t + 4 }] end]
						}	else {
							set yori [lindex $data [expr { $k2 + 5 }]]
							set zori [lindex $data [expr { $k2 + 6 }]]
						}
					}
					set xc $xori
					set yc $yori
					set zc $zori

					set ridc [lindex $data [expr { $k2 + 3 }]]

					set r2c [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]

					set x2 [expr { ($xc - $xr) * ($xc - $xr) }]
					set y2 [expr { ($yc - $yr) * ($yc - $yr) }]
					set z2 [expr { ($zc - $zr) * ($zc - $zr) }]
				
					set dis [expr { $x2 + $y2 + $z2 }]

					if { $dis < [expr { $nei * $nei }] && [expr { abs($r2c-$r2r) }] < 2.0 && $rid != $ridc } {
						puts $h "$ridc"
					}
				}
				incr k2
			}
			puts $h "\}"
		}
		incr k
	}
	close $h			

	# GETTING THE PDB INFORMATION IN THE FORM OF GLOBAL ARRAY FOR OUTER LEAFLET

	set k 0
	set n 0
	set cum_n $n
	set i 0
	set count 0
	set old_count -1

	while { $k < [llength $data] } {
		set term [lindex $data $k]
		set t1 [string range $term 0 5]
		if { [lindex $data $k] == "ATOM" || $t1 == "HETATM" } {
			if { $t1 == "HETATM" } {
				set sterm [string length $term]
				if { $sterm > 6 } {
					set shift 1
					set rn1 $term 
					set srn1 [string length $rn1]
					set rnum [string range $rn1 6 $srn1]
					set srn1 5
					set ft ""
				} else {
						set shift 0
						set rn1 [lindex $data [expr { $k + 1 }]]
						set rnum $rn1
						set srn1 [string length $rn1]
						set ft "HETATM"
				}
			} else { 
					set shift 0
					set rn1 [lindex $data [expr { $k + 1 }]]
					set srn1 [string length $rn1]
					set rnum $rn1
					set ft "ATOM  "
			}
		
			set atype [lindex $data [expr { $k + 2 - $shift}]]
			set satype [string length $atype]

			if { $satype > 5} {
				set shift2 1
				set at1 [lindex $data [expr { $k + 2 - $shift }]]		
				set sat1 4
				set resn ""
				set sresn 4
			} else {
				set shift2 0
				set at1 [lindex $data [expr { $k + 2 - $shift }]]		
				set sat1 [string length $at1]
				set resn [lindex $data [expr { $k + 3 - $shift}]]
			}
			set chain_id [lindex $data [expr { $k + 4 - $shift - $shift2}]]
			set schain_id [string length $chain_id]
				if { $schain_id > 1 } {
					set shift1 1
					set an1 [lindex $data [expr { $k + 4 - $shift - $shift2 }]]
					set san1 4
					set cn ""
				} else { 
					set cn [lindex $data [expr { $k + 4 - $shift - $shift2  }]]
					set an1 [lindex $data [expr { $k + 5 -$shift - $shift2 }]]
					set san1 [string length $an1]
					set shift1 0
				}

			set new_res $an1

			set x1 [lindex $data [expr { $k + 6 - $shift - $shift1 -$shift2}]]
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
					set z1 [lindex $data [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
					set z1 [format "%.3f" [expr { $z1 - 0.0 }]]
					set sz1 [string length $z1]
				}
			} else { 
				set shift3 0
				set x1 [format "%.3f" [expr { $x1 - 0.0 }]]
				set corx $x1
				set sx1 [string length $x1]
				set y1 [lindex $data [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]]
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
					set z1 [lindex $data [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
					set corz $z1
					set z1 [format "%.3f" [expr { $z1 - 0.0 }]]
					
					set sz1 [string length $z1]
				}
			}
			set rn1n($n) $rn1
			set at1n($n) $atype
			set an1n($n) $an1
			set resnamen($n) $resn
			set atomnamen($n) [lindex $data [expr { $k + 11 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]]
			set x1n($n) $corx
			set y1n($n) $cory
			set z1n($n) $corz

			# LINE OF ROTATION

			if { $atype == [lindex $patom 0] } {
				set v1x $corx
				set v1y $cory
				set v1z $corz
				incr count
			}

			if { $atype == [lindex $patom 1] && $resn == "OL" } {
				set v2x $corx
				set v2y $cory
				set v2z $corz
				incr count
			}

			if { $count != 0 && [expr { $count % 2 }] == 0 && $old_count != $count} {
				set vecx($i) [expr { $v1x - $v2x }]
				set vecy($i) [expr { $v1y - $v2y }]
				set vecz($i) [expr { $v1z - $v2z }]
				set vecpx($i) $v2x
				set vecpy($i) $v2y
				set vecpz($i) $v2z
				set old_count $count
				incr i
			}

			#if { $an1 > $nlipid_OL } {
				#set k [llength $data]
			#}	else {
				incr n
			#}
		}
		incr k
	}
	set an1n($n) -1

	set num_atom $n
	set cum_n $n
	set an1n(-1) 0
	set i 0
	set step [format "%.2f" [expr { 6.28 / 60.0 }]]
	for {set n 0} {$n < $num_atom} {incr n} {
		if { $an1n([expr { $n -1 }]) != $an1n($n) } {
			puts ""
			puts "				**** PUTTING IN RESIDUE $an1n($n) ****"
			puts " 				_______________________________________"
			set check [overlap_VES $an1n($n)]
			if { $check != 0 } {
				set otz 0.0
				for {set tz 0.0} {$tz <= 6.28} { set tz [expr { $tz + $step }]} {
					set tz [format "%.2f" $tz]
					#puts " Tz=$tz ABOUT $vecx($i) $vecy($i) $vecz($i)"
					transformations 0 [expr { $tz - $otz }] $an1n($n) $vecx($i) $vecy($i) $vecz($i) $vecpx($i) $vecpy($i) $vecpz($i)
					set otz $tz
					set check [overlap_VES $an1n($n)] 
					if { $check == 0 } {
						puts "				FINAL VALUE :: Tz=$tz ABOUT $vecx($i) $vecy($i) $vecz($i)"
						set tz 7.00
					}
				}
				set otz -3.14
			}
			set check [overlap_VES $an1n($n)]
			if { $check != 0 } {
				puts "RESIDUE $an1n($n) REMOVED"
			} else {
				set k3 0
				set k 0
				while { $an1n($k3) != $an1n($n) } { 
					incr k3
				}
				while { $an1n($k3) == $an1n($n) } {
					set srn1 [string length $rn1n($k3)]
					set sat1 [string length $at1n($k3)]
					set san1 [string length $an1n($k3)]
					set sx1 [string length $x1n($k3)]
					set sy1 [string length $y1n($k3)]
					set sz1 [string length $z1n($k3)]
					set sresname [string length $resnamen($k3)]
					puts $g "HETATM$p1($srn1)$rn1n($k3) $ic($sat1)$at1n($k3)$sat($sat1)$resnamen($k3)$p($sresname)$p1($san1)$an1n($k3)    $c($sx1)$x1n($k3)$c($sy1)$y1n($k3)$c($sz1)$z1n($k3)  1.00  0.00           $atomnamen($k3)"
					incr k3
				}
				puts $g "TER"
			}
		incr i
		}
	}
}

proc overlap_VES {nr} {

	#puts "overlap"

	set f [open "neighbours" "r"]
	set data [read $f]
	close $f

	set k 0
	while { [lindex $data $k 0] != $nr } {
		incr k
	}
	
	set pos $k
	set num_nei [llength [lindex $data $k]]

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n

	set k 0
	while { $an1n($k) != $nr } { 
		incr k
	}

	while { $an1n($k) == $nr } {
		set x1c $x1n($k)
		set y1c $y1n($k)
		set z1c $z1n($k)
		for {set i 0} {$i < $num_nei} {incr i} {
			set k1 0
			while { $an1n($k1) != [lindex $data $pos $i] } {
				#puts "[lindex $data $pos $i] $an1n($k1)"
				incr k1
			}
			set nr1 $an1n($k1)
			while { $an1n($k1) == $nr1 } { 
				if { $an1n($k1) != $nr } {
					set x2c $x1n($k1)
					set y2c $y1n($k1)
					set z2c $z1n($k1)

					set delx [expr { $x2c - $x1c }]
					set delx2 [expr { $delx * $delx }]

					set dely [expr { $y2c - $y1c }]
					set dely2 [expr { $dely * $dely }]

					set delz [expr { $z2c - $z1c }]
					set delz2 [expr { $delz * $delz }]

					set del [expr { $delx2 + $dely2 + $delz2 }]

					if { $del < 1.0 } {
						return -1
					}
				}
			incr k1
			}
		}
	incr k
	}
	return 0
}

proc transformations_VES { flag theta nr vecx vecy vecz vecpx vecpy vecpz} {

	package require math::linearalgebra

	#puts "transformations"

	global cum_n
	global x1n
	global y1n
	global z1n
	global an1n
	global xshift
	global yshift
	global zshift
	global resnamen

	set k 0
	while { $an1n($k) != $nr } { 
		incr k
	}
	set k1 $k

	# TRANSFORMATIONS (ROTATION AROUND A LINE WITH THE POSITION VECTOR (u,v,w))

	if { $flag == 0 } {

		set u $vecx
		set v $vecy
		set w $vecz
		set a $vecpx
		set b $vecpy
		set cc $vecpz
	
		set u2 [expr { $u * $u }]
		set v2 [expr { $v * $v }]
		set w2 [expr { $w * $w }]
		set l [expr { $u2 + $v2 + $w2 }]

		set rl [expr { sqrt($l) }]
		
		set a11 [expr { ($u2 + (($v2 + $w2) * cos($theta))) / $l } ]
		set a12 [expr { ((($u*$v)*(1-cos($theta))) - ($w * $rl * sin($theta))) / $l }]
		set a13 [expr { ((($u*$w)*(1-cos($theta))) + ($v * $rl * sin($theta))) / $l }]
		set a21 [expr { ((($u*$v)*(1-cos($theta))) + ($w * $rl * sin($theta))) / $l }]
		set a22 [expr { ($v2 + (($u2 + $w2) * cos($theta))) / $l }]
		set a23 [expr { ((($v*$w)*(1-cos($theta))) - ($u * $rl * sin($theta))) / $l }]
		set a31 [expr { ((($u*$w)*(1-cos($theta))) - ($v * $rl * sin($theta))) / $l }]
		set a32 [expr { ((($v*$w)*(1-cos($theta))) + ($u * $rl * sin($theta))) / $l }]
		set a33 [expr { ($w2 + (($u2 + $v2) * cos($theta))) / $l }]
		set a41 [expr { (((($a * ($v2 + $w2)) - ($u * (($b*$v) + ($cc*$w)))) * (1 - cos($theta))) + (((($b*$w) - ($cc*$v))) * ($rl * sin($theta)))) / $l }]
		set a42 [expr { (((($b * ($u2 + $w2)) - ($v * (($a*$u) + ($cc*$w)))) * (1 - cos($theta))) + (((($cc*$u) - ($a*$w))) * ($rl * sin($theta)))) / $l }]
		set a43 [expr { (((($cc * ($u2 + $v2)) - ($w * (($a*$u) + ($b*$v)))) * (1 - cos($theta))) + (((($a*$v) - ($b*$u))) * ($rl * sin($theta)))) / $l }]

		set row1 [list $a11 $a12 $a13]
		set row2 [list $a21 $a22 $a23]
		set row3 [list $a31 $a32 $a33]

		set tm [list $row1 $row2 $row3]

		while { $an1n($k) == $nr } {

			set tcoord [list $x1n($k) $y1n($k) $z1n($k)]
			set tcoord [::math::linearalgebra::matmul $tm $tcoord]

			set dum [open "dummy" "w"]
	
			puts $dum "$tcoord"

			close $dum

			set dum [open "dummy" "r"]
			set dvar [read $dum]
			close $dum

			set x1n($k) [format "%.3f" [lindex $dvar 0]]
			set sx1 [string length $x1n($k)]

			set y1n($k) [format "%.3f" [lindex $dvar 1]]
			set sy1 [string length $y1n($k)]

			set z1n($k) [format "%.3f" [lindex $dvar 2]]
			set sz1 [string length $z1n($k)]

			incr k
		}
	}
}

proc water_grid { max_layer gridsize_w lt R_LL R_LL_OLD rpore twat } {

	# FORMING A CONCENTRIC SPHERICAL GRID OF WATER MOLECULES INSIDE AND OUTSIDE THE VESICLE

	# MAINTAINING THE NUMBER DENSITY OF 0.033 / A^3

	# CHECKING THE WATER PORE REQUIREMENT

	if { $rpore != "NA" } {
		set pparm 1
	} else { 
		set pparm 0
	}

	set gridsize $gridsize_w

	set h [open "spherical_coord" "w"]
	set h1 [open "spherical_coord_ll" "w"]

	set R_UL [expr { $R_LL + (2*$lt) + (2.0) }]
	set height [expr { $R_UL * 2 }]

	set nwat_IL [expr { 4 * 3.14 * $R_LL_OLD * $R_LL_OLD * $R_LL_OLD * 0.011 }] 
	set nwat_IL [format "%.0f" $nwat_IL]
	set nwat_IL_AD [expr { $nwat_IL - ($twat / 2) }]

	set R_UL_OLD [expr { $R_LL_OLD + (2*$lt) + (2.0) }] 
	set R_M [expr { $R_UL_OLD + $max_layer }]

	set nwat_OL [expr { ($R_M * $R_M * $R_M * $nwat_IL) / ($R_LL_OLD * $R_LL_OLD * $R_LL_OLD) }]
	set nwat_OL [format "%.0f" $nwat_OL]

	puts "			#### PUTTING IN $nwat_OL WATER OUTSIDE AND $nwat_IL_AD INSIDE THE VESICLE ####"
	puts ""
	after 1000


	# INNER LAYER

	puts "				INNER LEAFLET"
	puts "				*************"
	puts ""

	#set k 0
	#while { $k < $nwat_IL } {
		#set u [expr { rand() }]
		#set v [expr { rand() }]
		#set u [format "%.2f" $u]
		#set v [format "%.2f" $v]
		#set r [expr { rand() * $R_LL }]
		#set phi [expr { acos((2*$v)-1) }]
		#set theta [expr { 2 * 3.14 * $u }]
		#set x($k) [expr { sin($phi) * cos($theta) }]
		#set y($k) [expr { sin($phi) * sin($theta) }]
		#set z($k) [expr { cos($phi) }]
		#set count 0
		#for {set i 0} {$i < $k} {incr i} {
			#set xc [expr { $x($k) - $x($i) }]
			#set yc [expr { $y($k) - $y($i) }]
			#set zc [expr { $z($k) - $z($i) }]			
			#set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
			#set dis [expr { $dis * $r * $r }]
			#if { $dis < [expr { $gridsize * $gridsize }] } {
				#incr count
				#set i $k
			#}
		#}
		#if { $count == 0 } {
			#puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat_IL ####"	
			#puts $h1 "[expr { $r * $x($k) }] [expr { $r * $y($k) }] [expr { $r * $z($k) }]"
			#incr k
		#}
	#}

	# CUBIC WATER LAYER

	set k 0
	while { $k < $nwat_IL_AD } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($R_LL*$m) }]
					set y($k) [expr { rand() * ($R_LL*$j) }]
					set z($k) [expr { rand() * ($R_LL*$l) }]
			
					set r [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r < [expr { $R_LL * $R_LL }] } {

						# CHECK PORE THE WATER PORE 

						if { $pparm == 1 } {
							set r1 [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) }]
							set r2 [expr { ($x($k)*$x($k)) + ($z($k)*$z($k)) }]
							set r3 [expr { ($z($k)*$z($k)) + ($y($k)*$y($k)) }]

							if { $r1 <= [expr { $rpore + 1.4 }] && $z($k) >= [expr { -1.0 * $height / 2.0 }] && $z($k) < 0 } {
								incr count
							}

							if { $r1 <= [expr { $rpore + 1.4 }] && $z($k) <= [expr { $height / 2.0 }] && $z($k) > 1.4 } {
								incr count
							}
			
							if { $r2 <= [expr { $rpore + 1.4 }] && $y($k) >= [expr { -1.0 * $height / 2.0 }] && $y($k) < 0 } {
								incr count
							}

							if { $r2 <= [expr { $rpore + 1.4 }] && $y($k) <= [expr { $height / 2.0 }] && $y($k) > 1.4 } {
								incr count
							}

							if { $r3 <= [expr { $rpore + 1.4 }] && $x($k) >= [expr { -1.0 * $height / 2.0 }] && $x($k) < 0 } {
								incr count
							}
							if { $r3 <= [expr { $rpore + 1.4 }] && $x($k) <= [expr { $height / 2.0 }] && $x($k) > 1.4 } {
								incr count
							}
						}

						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat_IL_AD ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	# OUTER LAYER

	puts "				*************"
	puts "				OUTER LEAFLET"
	puts "				*************"
	puts ""

	#set k 0
	#while { $k < $nwat_OL } {
		#set u [expr { rand() }]
		#set v [expr { rand() }]
		#set u [format "%.2f" $u]
		#set v [format "%.2f" $v]
		#set r [expr { rand() * ($R_UL + $max_layer) }]
		#if { $r > [expr { $R_UL + $gridsize }] } {
			#set phi [expr { acos((2*$v)-1) }]
			#set theta [expr { 2 * 3.14 * $u }]
			#set x($k) [expr { sin($phi) * cos($theta) }]
			#set y($k) [expr { sin($phi) * sin($theta) }]
			#set z($k) [expr { cos($phi) }]
			#set count 0
			#for {set i 0} {$i < $k} {incr i} {
				#set xc [expr { $x($k) - $x($i) }]
				#set yc [expr { $y($k) - $y($i) }]
				#set zc [expr { $z($k) - $z($i) }]			
				#set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
				#set dis [expr { $dis * $R_UL * $R_UL }]
				#if { $dis < [expr { $gridsize * $gridsize }] } {
					#incr count
					#set i $k
				#}
			#}
			#if { $count == 0 } {	
				#puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat_OL ####"	
				#puts $h "[expr { $r * $x($k) }] [expr { $r * $y($k) }] [expr { $r * $z($k) }]"
				#incr k
			#}
		#}
	#}

	# CUBIC WATER LAYER

	set k 0
	while { $k < $nwat_OL } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * (($R_UL + $max_layer)*$m) }]
					set y($k) [expr { rand() * (($R_UL + $max_layer)*$j) }]
					set z($k) [expr { rand() * (($R_UL + $max_layer)*$l) }]
			
					set r [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r > [expr { $R_UL * $R_UL }] } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat_OL ####"	
							puts $h "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}
	close $h
	close $h1

	water $nwat_IL_AD $nwat_OL
	
}		

proc water_grid_dibs { max_layer gridsize_w lt R_LL } {

	# FORMING A CONCENTRIC SPHERICAL GRID OF WATER MOLECULES INSIDE AND OUTSIDE THE VESICLE

	# MAINTAINING THE NUMBER DENSITY OF 0.033 / A^3

	set gridsize $gridsize_w

	set h [open "spherical_coord" "w"]
	set h1 [open "spherical_coord_ll" "w"]

	set nwat_IL [expr { 4 * 3.14 * $R_LL * $R_LL * $R_LL * 0.011 }] 
	set nwat_IL [format "%.0f" $nwat_IL]

	puts "			#### PUTTING IN $nwat_IL INSIDE THE DIBs ####"
	puts ""
	after 1000


	# INNER LAYER OF BOTH

	puts "				INNER LEAFLET"
	puts "				*************"
	puts ""

	# CUBIC WATER LAYER

	set k 0
	while { $k < $nwat_IL } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($R_LL*$m) }]
					set y($k) [expr { rand() * ($R_LL*$j) }]
					set z($k) [expr { rand() * ($R_LL*$l) }]
			
					set r [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r < [expr { $R_LL * $R_LL }] } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat_IL ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							puts $h "[expr { $x($k) + (2*($R_LL + $lt)) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	close $h
	close $h1

	water $nwat_IL $nwat_IL
	
}		

proc water { nwat_IL nwat_OL } {

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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set h [open "spherical_coord" "r"]
	set data [read $h]
	close $h

	set h1 [open "spherical_coord_ll" "r"]
	set data1 [read $h1]
	close $h1

	set f [open "sol_water.pdb" "w"]


	# PUTTING THE WATER IN THE OUTER LEAFLET
	
	set k 0
	set j 1
	set ij 1
	for {set i 0} {$i < $nwat_OL} {incr i} {
		if { $j > 99999 } {
			set j 1
		}
		set rn1 $j
		set srn1 [string length $rn1]
		incr j

		set at1 "O"
		set sat1 [string length $at1]

		set an1 [expr { $i + 1 }] 
		set san1 [string length $an1]

		if { $an1 > 99999 } {	
			while { $an1 > 99999 } {
				set an1 [expr { $an1 - 99999 }]
			}
		}
		set san1 [string length $an1]

		set resname "WAT"
		set sresname [string length $resname]

		set x1 [lindex $data $k]
		set x1 [format "%.3f" $x1]
		set sx1 [string length $x1]

		set y1 [lindex $data [expr { $k + 1 }]]
		set y1 [format "%.3f" $y1]
		set sy1 [string length $y1]

		set z1 [lindex $data [expr { $k + 2 }]]
		set z1 [format "%.3f" $z1]
		set sz1 [string length $z1]

		incr k 3
		
		puts $f "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname $p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
		puts $f "TER"
	}

	# PUTTING WATER IN THE INNER LEAFLET
	
	set k 0
	for {set i $nwat_OL} {$i < [expr { $nwat_IL + $nwat_OL }]} {incr i} {

		if { $j > 99999 } {
			set j 1
		}

		set rn1 $j
		set srn1 [string length $rn1]
		incr j

		set at1 "O"
		set sat1 [string length $at1]

		set an1 [expr { $i + 1 }] 
		set san1 [string length $an1]

		if { $an1 > 99999 } {	
			while { $an1 > 99999 } {
				set an1 [expr { $an1 - 99999 }]
			}
		}
		set san1 [string length $an1]

		set resname "WAT"
		set sresname [string length $resname]

		set x1 [lindex $data1 $k]
		set x1 [format "%.3f" $x1]
		set sx1 [string length $x1]

		set y1 [lindex $data1 [expr { $k + 1 }]]
		set y1 [format "%.3f" $y1]
		set sy1 [string length $y1]

		set z1 [lindex $data1 [expr { $k + 2 }]]
		set z1 [format "%.3f" $z1]
		set sz1 [string length $z1]
	
		incr k 3
		
		puts $f "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname $p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
		puts $f "TER"
	}
	close $f
	
	comb
}

proc water_pore { lt R_LL rpore gridsize_w} {

	# FORMING A CONCENTRIC CYLINDRICAL GRID OF WATER PORE

	# MAINTAINING THE NUMBER DENSITY OF 0.015 / A^3

	set gridsize $gridsize_w

	set h1 [open "spherical_coord_wp" "w"]

	set R_UL [expr { $R_LL + (2*$lt) + (2.0) }]

	set height [expr { $R_UL * 2 }]

	set nwat [expr { 3.14 * $rpore * $rpore * $height * 0.015 }] 
	set nwat [format "%.0f" $nwat]

	puts "			#### FORMING A CYLINDRICAL WATER PORE OF RADIUS $rpore Ang AND HEIGHT $height ####"
	puts ""
	after 1000


	# WATER PORE 1,2,3,4,5,6,

	puts ""
	puts "				WATER PORE 1 ALONG Z"
	puts "				********************"
	puts ""

	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($rpore*$m) }]
					set y($k) [expr { rand() * ($rpore*$j) }]
					set z($k) [expr { rand() * ($height*$l) }]
			
					set r [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) }]
					set count 0
					if { $r <= $rpore && $z($k) >= [expr { -1.0 * $height / 2.0 }] && $z($k) < 0 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	puts ""
	puts "				WATER PORE 2 ALONG Z"
	puts "				********************"
	puts ""


	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($rpore*$m) }]
					set y($k) [expr { rand() * ($rpore*$j) }]
					set z($k) [expr { rand() * ($height*$l) }]
			
					set r [expr { ($x($k)*$x($k)) + ($y($k)*$y($k)) }]
					set count 0
					if { $r <= $rpore && $z($k) <= [expr { $height / 2.0 }] && $z($k) > 1.4 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	puts ""
	puts "				WATER PORE 3 ALONG Y"
	puts "				********************"
	puts ""

	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($rpore*$m) }]
					set y($k) [expr { rand() * ($height*$l) }]
					set z($k) [expr { rand() * ($rpore*$j) }]
			
					set r [expr { ($x($k)*$x($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r <= $rpore && $y($k) >= [expr { -1.0 * $height / 2.0 }] && $y($k) < 0 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	puts ""
	puts "				WATER PORE 4 ALONG Y"
	puts "				********************"
	puts ""

	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($rpore*$m) }]
					set y($k) [expr { rand() * ($height*$l) }]
					set z($k) [expr { rand() * ($rpore*$j) }]
			
					set r [expr { ($x($k)*$x($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r <= $rpore && $y($k) <= [expr { $height / 2.0 }] && $y($k) > 1.4 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	puts ""
	puts "				WATER PORE 5 ALONG X"
	puts "				********************"
	puts ""

	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($height*$m) }]
					set y($k) [expr { rand() * ($rpore*$j) }]
					set z($k) [expr { rand() * ($rpore*$l) }]
			
					set r [expr { ($y($k)*$y($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r <= $rpore && $x($k) >= [expr { -1.0 * $height / 2.0 }] && $x($k) <  0 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}

	puts ""
	puts "				WATER PORE 6 ALONG X"
	puts "				********************"
	puts ""

	set k 0
	while { $k < $nwat } {
		for {set m -1} {$m <= 1} {incr m 2} {
			for {set j -1} {$j<= 1} {incr j 2} {
				for {set l -1} {$l <= 1 } {incr l 2} {
					set x($k) [expr { rand() * ($height*$m) }]
					set y($k) [expr { rand() * ($rpore*$j) }]
					set z($k) [expr { rand() * ($rpore*$l) }]
			
					set r [expr { ($y($k)*$y($k)) + ($z($k)*$z($k)) }]
					set count 0
					if { $r <= $rpore && $x($k) <= [expr { $height / 2.0 }] && $x($k) > 1.4 } {
						for {set i 0} {$i < $k} {incr i} {
							set xc [expr { $x($k) - $x($i) }]
							set yc [expr { $y($k) - $y($i) }]
							set zc [expr { $z($k) - $z($i) }]			
							set dis [expr { ($xc*$xc) + ($yc*$yc) + ($zc*$zc) }]
							if { $dis < [expr { $gridsize * $gridsize }] } {
								incr count
								set i $k
							}
						}
						if { $count == 0 } {	
							puts "				#### PUTTING IN WATER [expr { $k + 1 }] OF $nwat ####"	
							puts $h1 "[expr { $x($k) }] [expr { $y($k) }] [expr { $z($k) }]"
							incr k
						}
					}
				}
			}
		}	
	}
	close $h1
	
	set twat [expr { 6 * $nwat }]
	water_wp $twat
	return $twat
}		

proc water_wp { twat } {

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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set h1 [open "spherical_coord_wp" "r"]
	set data1 [read $h1]
	close $h1

	set f [open "water_pore.pdb" "w"]


	# PUTTING THE WATER IN THE OUTER LEAFLET
	
	set k 0
	set j 1
	set ij 1
	for {set i 0} {$i < $twat} {incr i} {
		if { $j > 99999 } {
			set j 1
		}
		set rn1 $j
		set srn1 [string length $rn1]
		incr j

		set at1 "O"
		set sat1 [string length $at1]

		set an1 [expr { $i + 1 }] 
		set san1 [string length $an1]

		if { $an1 > 99999 } {	
			while { $an1 > 99999 } {
				set an1 [expr { $an1 - 99999 }]
			}
		}
		set san1 [string length $an1]

		set resname "WAT"
		set sresname [string length $resname]

		set x1 [lindex $data1 $k]
		set x1 [format "%.3f" $x1]
		set sx1 [string length $x1]

		set y1 [lindex $data1 [expr { $k + 1 }]]
		set y1 [format "%.3f" $y1]
		set sy1 [string length $y1]

		set z1 [lindex $data1 [expr { $k + 2 }]]
		set z1 [format "%.3f" $z1]
		set sz1 [string length $z1]

		incr k 3
		
		puts $f "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname $p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
		puts $f "TER"
	}
	close $f
}

proc comb {} {
	set f [open "lipids_wo.pdb" "r"]
	set data [read $f]
	close $f

	set g1 [open "water_pore.pdb" "r"]
	set data3 [read $g1]
	close $g1

	set g [open "sol_water.pdb" "r"]
	set data1 [read $g]
	close $g

	set h [open "sol_lipid.pdb" "w"]

	set data2 $data$data3$data1

	puts $h "$data2"

	close $h
}

#################################################################### VESCICLE BUILDER END ##############################################################################################################################


proc input_PDB {inp} {

	set p1 { DLPC DMPC DOPC DPPC POPC DLPE DMPE DOPE DPPE POPE DOPG DOPS POPS POPG CHL} 

	# DLPE

set DLPE {HETATM    1  C12 LA  a   2      -8.388   3.147   8.309  1.00  0.00           C
HETATM    2  H2R LA  a   2      -9.342   3.693   8.149  1.00  0.00           H
HETATM    3  H2S LA  a   2      -8.695   2.080   8.271  1.00  0.00           H
HETATM    4  C13 LA  a   2      -7.435   3.433   7.116  1.00  0.00           C
HETATM    5  H3R LA  a   2      -7.762   2.754   6.300  1.00  0.00           H
HETATM    6  H3S LA  a   2      -6.355   3.385   7.372  1.00  0.00           H
HETATM    7  C14 LA  a   2      -7.607   4.877   6.461  1.00  0.00           C
HETATM    8  H4R LA  a   2      -7.115   5.604   7.143  1.00  0.00           H
HETATM    9  H4S LA  a   2      -8.699   5.061   6.373  1.00  0.00           H
HETATM   10  C15 LA  a   2      -7.026   5.095   5.025  1.00  0.00           C
HETATM   11  H5R LA  a   2      -5.923   5.224   5.044  1.00  0.00           H
HETATM   12  H5S LA  a   2      -7.519   6.018   4.649  1.00  0.00           H
HETATM   13  C16 LA  a   2      -7.354   3.902   4.110  1.00  0.00           C
HETATM   14  H6R LA  a   2      -8.386   3.572   4.355  1.00  0.00           H
HETATM   15  H6S LA  a   2      -6.630   3.080   4.296  1.00  0.00           H
HETATM   16  C17 LA  a   2      -7.349   4.241   2.618  1.00  0.00           C
HETATM   17  H7R LA  a   2      -8.064   5.055   2.373  1.00  0.00           H
HETATM   18  H7S LA  a   2      -7.569   3.454   1.866  1.00  0.00           H
HETATM   19  C18 LA  a   2      -5.913   4.732   2.135  1.00  0.00           C
HETATM   20  H8R LA  a   2      -5.286   3.826   1.995  1.00  0.00           H
HETATM   21  H8S LA  a   2      -5.441   5.385   2.900  1.00  0.00           H
HETATM   22  C19 LA  a   2      -5.951   5.612   0.819  1.00  0.00           C
HETATM   23  H9R LA  a   2      -4.965   6.110   0.702  1.00  0.00           H
HETATM   24  H9S LA  a   2      -6.762   6.343   1.024  1.00  0.00           H
HETATM   25 C110 LA  a   2      -6.417   4.787  -0.442  1.00  0.00           C
HETATM   26 H10R LA  a   2      -7.504   4.608  -0.295  1.00  0.00           H
HETATM   27 H10S LA  a   2      -5.830   3.848  -0.360  1.00  0.00           H
HETATM   28 C111 LA  a   2      -6.190   5.448  -1.735  1.00  0.00           C
HETATM   29 H11R LA  a   2      -5.133   5.761  -1.592  1.00  0.00           H
HETATM   30 H11S LA  a   2      -7.037   6.159  -1.840  1.00  0.00           H
HETATM   31 C112 LA  a   2      -6.339   4.450  -2.877  1.00  0.00           C
HETATM   32 H12R LA  a   2      -7.404   4.154  -2.769  1.00  0.00           H
HETATM   33 H12S LA  a   2      -5.615   3.608  -2.837  1.00  0.00           H
HETATM   34 H12T LA  a   2      -6.117   5.003  -3.814  1.00  0.00           H
HETATM   35  N31 PE  a   2      -6.541   5.647  17.156  1.00  0.00           N
HETATM   36 HN1A PE  a   2      -5.769   5.902  17.803  1.00  0.00           H
HETATM   37 HN1B PE  a   2      -6.541   4.630  16.943  1.00  0.00           H
HETATM   38 HN1C PE  a   2      -7.392   6.029  17.616  1.00  0.00           H
HETATM   39  C32 PE  a   2      -6.110   6.347  15.832  1.00  0.00           C
HETATM   40  H2A PE  a   2      -6.098   7.429  16.085  1.00  0.00           H
HETATM   41  H2B PE  a   2      -6.902   6.275  15.056  1.00  0.00           H
HETATM   42  C31 PE  a   2      -4.784   5.666  15.312  1.00  0.00           C
HETATM   43  H1A PE  a   2      -3.974   6.060  15.963  1.00  0.00           H
HETATM   44  H1B PE  a   2      -4.636   5.702  14.212  1.00  0.00           H
HETATM   45  P31 PE  a   2      -3.930   3.254  14.802  1.00  0.00           P
HETATM   46  O33 PE  a   2      -3.963   1.943  15.547  1.00  0.00           O
HETATM   47  O34 PE  a   2      -2.601   3.866  14.428  1.00  0.00           O
HETATM   48  O31 PE  a   2      -4.606   3.054  13.477  1.00  0.00           O
HETATM   49  O32 PE  a   2      -4.848   4.270  15.608  1.00  0.00           O
HETATM   50  C3  PE  a   2      -4.141   2.092  12.546  1.00  0.00           C
HETATM   51  HA  PE  a   2      -4.439   1.078  12.890  1.00  0.00           H
HETATM   52  HB  PE  a   2      -3.037   2.102  12.420  1.00  0.00           H
HETATM   53  C2  PE  a   2      -4.771   2.358  11.120  1.00  0.00           C
HETATM   54  HX  PE  a   2      -4.382   1.572  10.438  1.00  0.00           H
HETATM   55  O21 PE  a   2      -4.511   3.711  10.696  1.00  0.00           O
HETATM   56  C21 PE  a   2      -4.034   3.897   9.463  1.00  0.00           C
HETATM   57  O22 PE  a   2      -3.780   3.029   8.675  1.00  0.00           O
HETATM   58  C1  PE  a   2      -6.282   2.156  11.150  1.00  0.00           C
HETATM   59  HR  PE  a   2      -6.601   2.982  11.821  1.00  0.00           H
HETATM   60  HS  PE  a   2      -6.466   1.161  11.609  1.00  0.00           H
HETATM   61  O11 PE  a   2      -6.942   2.230   9.880  1.00  0.00           O
HETATM   62  C11 PE  a   2      -7.644   3.328   9.664  1.00  0.00           C
HETATM   63  O12 PE  a   2      -7.758   4.243  10.426  1.00  0.00           O
HETATM   64  C12 LA  a   2      -3.959   5.368   9.202  1.00  0.00           C
HETATM   65  H2R LA  a   2      -3.351   5.837  10.005  1.00  0.00           H
HETATM   66  H2S LA  a   2      -4.998   5.762   9.212  1.00  0.00           H
HETATM   67  C13 LA  a   2      -3.318   5.621   7.791  1.00  0.00           C
HETATM   68  H3R LA  a   2      -3.748   5.005   6.973  1.00  0.00           H
HETATM   69  H3S LA  a   2      -2.306   5.179   7.908  1.00  0.00           H
HETATM   70  C14 LA  a   2      -3.491   7.162   7.554  1.00  0.00           C
HETATM   71  H4R LA  a   2      -3.185   7.676   8.490  1.00  0.00           H
HETATM   72  H4S LA  a   2      -4.593   7.290   7.490  1.00  0.00           H
HETATM   73  C15 LA  a   2      -2.687   7.698   6.347  1.00  0.00           C
HETATM   74  H5R LA  a   2      -1.873   6.948   6.245  1.00  0.00           H
HETATM   75  H5S LA  a   2      -2.177   8.684   6.388  1.00  0.00           H
HETATM   76  C16 LA  a   2      -3.527   7.668   4.996  1.00  0.00           C
HETATM   77  H6R LA  a   2      -4.260   8.502   5.028  1.00  0.00           H
HETATM   78  H6S LA  a   2      -4.174   6.768   5.066  1.00  0.00           H
HETATM   79  C17 LA  a   2      -2.799   7.833   3.669  1.00  0.00           C
HETATM   80  H7R LA  a   2      -3.416   8.228   2.834  1.00  0.00           H
HETATM   81  H7S LA  a   2      -1.961   8.562   3.639  1.00  0.00           H
HETATM   82  C18 LA  a   2      -2.131   6.538   3.329  1.00  0.00           C
HETATM   83  H8R LA  a   2      -1.139   6.414   3.812  1.00  0.00           H
HETATM   84  H8S LA  a   2      -2.760   5.671   3.624  1.00  0.00           H
HETATM   85  C19 LA  a   2      -1.885   6.348   1.758  1.00  0.00           C
HETATM   86  H9R LA  a   2      -1.195   7.159   1.439  1.00  0.00           H
HETATM   87  H9S LA  a   2      -2.901   6.512   1.340  1.00  0.00           H
HETATM   88 C110 LA  a   2      -1.439   4.932   1.399  1.00  0.00           C
HETATM   89 H10R LA  a   2      -2.345   4.361   1.104  1.00  0.00           H
HETATM   90 H10S LA  a   2      -0.936   4.447   2.263  1.00  0.00           H
HETATM   91 C111 LA  a   2      -0.462   5.000   0.226  1.00  0.00           C
HETATM   92 H11R LA  a   2      -0.931   5.553  -0.616  1.00  0.00           H
HETATM   93 H11S LA  a   2      -0.186   3.957  -0.036  1.00  0.00           H
HETATM   94 C112 LA  a   2       0.866   5.616   0.698  1.00  0.00           C
HETATM   95 H12R LA  a   2       0.980   5.628   1.803  1.00  0.00           H
HETATM   96 H12S LA  a   2       1.030   6.651   0.330  1.00  0.00           H
HETATM   97 H12T LA  a   2       1.732   5.118   0.212  1.00  0.00           H
END}

	set f [open "DLPE_A.pdb" "w"]
	puts $f "$DLPE"
	close $f

	# DLPC

	set DLPC {HETATM    1  C12 LA  a   4      -8.388   3.147   8.309  1.00  0.00           C
HETATM    2  H2R LA  a   4      -9.342   3.693   8.149  1.00  0.00           H
HETATM    3  H2S LA  a   4      -8.695   2.080   8.271  1.00  0.00           H
HETATM    4  C13 LA  a   4      -7.435   3.433   7.116  1.00  0.00           C
HETATM    5  H3R LA  a   4      -7.762   2.754   6.300  1.00  0.00           H
HETATM    6  H3S LA  a   4      -6.355   3.385   7.372  1.00  0.00           H
HETATM    7  C14 LA  a   4      -7.607   4.877   6.461  1.00  0.00           C
HETATM    8  H4R LA  a   4      -7.115   5.604   7.143  1.00  0.00           H
HETATM    9  H4S LA  a   4      -8.699   5.061   6.373  1.00  0.00           H
HETATM   10  C15 LA  a   4      -7.026   5.095   5.025  1.00  0.00           C
HETATM   11  H5R LA  a   4      -5.923   5.224   5.044  1.00  0.00           H
HETATM   12  H5S LA  a   4      -7.519   6.018   4.649  1.00  0.00           H
HETATM   13  C16 LA  a   4      -7.354   3.902   4.110  1.00  0.00           C
HETATM   14  H6R LA  a   4      -8.386   3.572   4.355  1.00  0.00           H
HETATM   15  H6S LA  a   4      -6.630   3.080   4.296  1.00  0.00           H
HETATM   16  C17 LA  a   4      -7.349   4.241   2.618  1.00  0.00           C
HETATM   17  H7R LA  a   4      -8.064   5.055   2.373  1.00  0.00           H
HETATM   18  H7S LA  a   4      -7.569   3.454   1.866  1.00  0.00           H
HETATM   19  C18 LA  a   4      -5.913   4.732   2.135  1.00  0.00           C
HETATM   20  H8R LA  a   4      -5.286   3.826   1.995  1.00  0.00           H
HETATM   21  H8S LA  a   4      -5.441   5.385   2.900  1.00  0.00           H
HETATM   22  C19 LA  a   4      -5.951   5.612   0.819  1.00  0.00           C
HETATM   23  H9R LA  a   4      -4.965   6.110   0.702  1.00  0.00           H
HETATM   24  H9S LA  a   4      -6.762   6.343   1.024  1.00  0.00           H
HETATM   25 C110 LA  a   4      -6.417   4.787  -0.442  1.00  0.00           C
HETATM   26 H10R LA  a   4      -7.504   4.608  -0.295  1.00  0.00           H
HETATM   27 H10S LA  a   4      -5.830   3.848  -0.360  1.00  0.00           H
HETATM   28 C111 LA  a   4      -6.190   5.448  -1.735  1.00  0.00           C
HETATM   29 H11R LA  a   4      -5.133   5.761  -1.592  1.00  0.00           H
HETATM   30 H11S LA  a   4      -7.037   6.159  -1.840  1.00  0.00           H
HETATM   31 C112 LA  a   4      -6.339   4.450  -2.877  1.00  0.00           C
HETATM   32 H12R LA  a   4      -7.404   4.154  -2.769  1.00  0.00           H
HETATM   33 H12S LA  a   4      -5.615   3.608  -2.837  1.00  0.00           H
HETATM   34 H12T LA  a   4      -6.117   5.003  -3.814  1.00  0.00           H
HETATM   35  N31 PC  a   4      -7.390   5.865  15.217  1.00  0.00           N
HETATM   36  C33 PC  a   4      -7.342   6.274  13.751  1.00  0.00           C
HETATM   37  H3A PC  a   4      -8.194   5.836  13.252  1.00  0.00           H
HETATM   38  H3B PC  a   4      -6.413   6.046  13.250  1.00  0.00           H
HETATM   39  H3C PC  a   4      -7.549   7.329  13.646  1.00  0.00           H
HETATM   40  C34 PC  a   4      -8.584   6.623  15.827  1.00  0.00           C
HETATM   41  H4A PC  a   4      -9.470   6.121  15.469  1.00  0.00           H
HETATM   42  H4B PC  a   4      -8.636   7.630  15.440  1.00  0.00           H
HETATM   43  H4C PC  a   4      -8.493   6.689  16.901  1.00  0.00           H
HETATM   44  C35 PC  a   4      -7.686   4.338  15.221  1.00  0.00           C
HETATM   45  H5A PC  a   4      -6.925   3.772  14.705  1.00  0.00           H
HETATM   46  H5B PC  a   4      -8.693   4.228  14.845  1.00  0.00           H
HETATM   47  H5C PC  a   4      -7.751   4.080  16.268  1.00  0.00           H
HETATM   48  C32 PC  a   4      -6.054   6.318  15.810  1.00  0.00           C
HETATM   49  H2A PC  a   4      -6.056   6.150  16.877  1.00  0.00           H
HETATM   50  H2B PC  a   4      -5.969   7.351  15.505  1.00  0.00           H
HETATM   51  C31 PC  a   4      -4.784   5.666  15.312  1.00  0.00           C
HETATM   52  H1A PC  a   4      -3.974   6.060  15.963  1.00  0.00           H
HETATM   53  H1B PC  a   4      -4.636   5.702  14.212  1.00  0.00           H
HETATM   54  P31 PC  a   4      -3.930   3.254  14.802  1.00  0.00           P
HETATM   55  O33 PC  a   4      -3.963   1.943  15.547  1.00  0.00           O
HETATM   56  O34 PC  a   4      -2.601   3.866  14.428  1.00  0.00           O
HETATM   57  O31 PC  a   4      -4.606   3.054  13.477  1.00  0.00           O
HETATM   58  O32 PC  a   4      -4.848   4.270  15.608  1.00  0.00           O
HETATM   59  C3  PC  a   4      -4.141   2.092  12.546  1.00  0.00           C
HETATM   60  HA  PC  a   4      -4.439   1.078  12.890  1.00  0.00           H
HETATM   61  HB  PC  a   4      -3.037   2.102  12.420  1.00  0.00           H
HETATM   62  C2  PC  a   4      -4.771   2.358  11.120  1.00  0.00           C
HETATM   63  HX  PC  a   4      -4.382   1.572  10.438  1.00  0.00           H
HETATM   64  O21 PC  a   4      -4.511   3.711  10.696  1.00  0.00           O
HETATM   65  C21 PC  a   4      -4.034   3.897   9.463  1.00  0.00           C
HETATM   66  O22 PC  a   4      -3.780   3.029   8.675  1.00  0.00           O
HETATM   67  C1  PC  a   4      -6.282   2.156  11.150  1.00  0.00           C
HETATM   68  HR  PC  a   4      -6.601   2.982  11.821  1.00  0.00           H
HETATM   69  HS  PC  a   4      -6.466   1.161  11.609  1.00  0.00           H
HETATM   70  O11 PC  a   4      -6.942   2.230   9.880  1.00  0.00           O
HETATM   71  C11 PC  a   4      -7.644   3.328   9.664  1.00  0.00           C
HETATM   72  O12 PC  a   4      -7.758   4.243  10.426  1.00  0.00           O
HETATM   73  C12 LA  a   4      -3.959   5.368   9.202  1.00  0.00           C
HETATM   74  H2R LA  a   4      -3.351   5.837  10.005  1.00  0.00           H
HETATM   75  H2S LA  a   4      -4.998   5.762   9.212  1.00  0.00           H
HETATM   76  C13 LA  a   4      -3.318   5.621   7.791  1.00  0.00           C
HETATM   77  H3R LA  a   4      -3.748   5.005   6.973  1.00  0.00           H
HETATM   78  H3S LA  a   4      -2.306   5.179   7.908  1.00  0.00           H
HETATM   79  C14 LA  a   4      -3.491   7.162   7.554  1.00  0.00           C
HETATM   80  H4R LA  a   4      -3.185   7.676   8.490  1.00  0.00           H
HETATM   81  H4S LA  a   4      -4.593   7.290   7.490  1.00  0.00           H
HETATM   82  C15 LA  a   4      -2.687   7.698   6.347  1.00  0.00           C
HETATM   83  H5R LA  a   4      -1.873   6.948   6.245  1.00  0.00           H
HETATM   84  H5S LA  a   4      -2.177   8.684   6.388  1.00  0.00           H
HETATM   85  C16 LA  a   4      -3.527   7.668   4.996  1.00  0.00           C
HETATM   86  H6R LA  a   4      -4.260   8.502   5.028  1.00  0.00           H
HETATM   87  H6S LA  a   4      -4.174   6.768   5.066  1.00  0.00           H
HETATM   88  C17 LA  a   4      -2.799   7.833   3.669  1.00  0.00           C
HETATM   89  H7R LA  a   4      -3.416   8.228   2.834  1.00  0.00           H
HETATM   90  H7S LA  a   4      -1.961   8.562   3.639  1.00  0.00           H
HETATM   91  C18 LA  a   4      -2.131   6.538   3.329  1.00  0.00           C
HETATM   92  H8R LA  a   4      -1.139   6.414   3.812  1.00  0.00           H
HETATM   93  H8S LA  a   4      -2.760   5.671   3.624  1.00  0.00           H
HETATM   94  C19 LA  a   4      -1.885   6.348   1.758  1.00  0.00           C
HETATM   95  H9R LA  a   4      -1.195   7.159   1.439  1.00  0.00           H
HETATM   96  H9S LA  a   4      -2.901   6.512   1.340  1.00  0.00           H
HETATM   97 C110 LA  a   4      -1.439   4.932   1.399  1.00  0.00           C
HETATM   98 H10R LA  a   4      -2.345   4.361   1.104  1.00  0.00           H
HETATM   99 H10S LA  a   4      -0.936   4.447   2.263  1.00  0.00           H
HETATM  100 C111 LA  a   4      -0.462   5.000   0.226  1.00  0.00           C
HETATM  101 H11R LA  a   4      -0.931   5.553  -0.616  1.00  0.00           H
HETATM  102 H11S LA  a   4      -0.186   3.957  -0.036  1.00  0.00           H
HETATM  103 C112 LA  a   4       0.866   5.616   0.698  1.00  0.00           C
HETATM  104 H12R LA  a   4       0.980   5.628   1.803  1.00  0.00           H
HETATM  105 H12S LA  a   4       1.030   6.651   0.330  1.00  0.00           H
HETATM  106 H12T LA  a   4       1.732   5.118   0.212  1.00  0.00           H
END}

	set f [open "DLPC_A.pdb" "w"]
	puts $f "$DLPC"
	close $f

	# DMPC

	set DMPC {HETATM    1  C12 MY  a  39     -14.281 -22.942  10.676  1.00  0.00           C
HETATM    2  H2R MY  a  39     -14.466 -23.933  11.145  1.00  0.00           H
HETATM    3  H2S MY  a  39     -13.394 -22.432  11.108  1.00  0.00           H
HETATM    4  C13 MY  a  39     -14.256 -23.097   9.150  1.00  0.00           C
HETATM    5  H3R MY  a  39     -13.345 -23.695   8.933  1.00  0.00           H
HETATM    6  H3S MY  a  39     -13.998 -22.137   8.654  1.00  0.00           H
HETATM    7  C14 MY  a  39     -15.390 -23.921   8.564  1.00  0.00           C
HETATM    8  H4R MY  a  39     -16.342 -23.363   8.689  1.00  0.00           H
HETATM    9  H4S MY  a  39     -15.404 -24.907   9.076  1.00  0.00           H
HETATM   10  C15 MY  a  39     -15.207 -24.291   7.102  1.00  0.00           C
HETATM   11  H5R MY  a  39     -16.203 -24.196   6.620  1.00  0.00           H
HETATM   12  H5S MY  a  39     -14.855 -25.343   7.042  1.00  0.00           H
HETATM   13  C16 MY  a  39     -14.337 -23.298   6.348  1.00  0.00           C
HETATM   14  H6R MY  a  39     -13.293 -23.360   6.722  1.00  0.00           H
HETATM   15  H6S MY  a  39     -14.622 -22.227   6.274  1.00  0.00           H
HETATM   16  C17 MY  a  39     -14.170 -23.763   4.853  1.00  0.00           C
HETATM   17  H7R MY  a  39     -13.721 -24.769   4.988  1.00  0.00           H
HETATM   18  H7S MY  a  39     -13.316 -23.171   4.461  1.00  0.00           H
HETATM   19  C18 MY  a  39     -15.420 -23.756   3.967  1.00  0.00           C
HETATM   20  H8R MY  a  39     -15.615 -22.683   3.755  1.00  0.00           H
HETATM   21  H8S MY  a  39     -16.288 -24.113   4.561  1.00  0.00           H
HETATM   22  C19 MY  a  39     -15.126 -24.518   2.639  1.00  0.00           C
HETATM   23  H9R MY  a  39     -15.102 -25.629   2.667  1.00  0.00           H
HETATM   24  H9S MY  a  39     -14.183 -24.088   2.238  1.00  0.00           H
HETATM   25 C110 MY  a  39     -16.022 -24.102   1.457  1.00  0.00           C
HETATM   26 H10R MY  a  39     -15.945 -22.999   1.572  1.00  0.00           H
HETATM   27 H10S MY  a  39     -17.098 -24.339   1.603  1.00  0.00           H
HETATM   28 C111 MY  a  39     -15.509 -24.533   0.136  1.00  0.00           C
HETATM   29 H11R MY  a  39     -15.297 -25.620   0.212  1.00  0.00           H
HETATM   30 H11S MY  a  39     -14.519 -24.061  -0.042  1.00  0.00           H
HETATM   31 C112 MY  a  39     -16.588 -24.412  -0.936  1.00  0.00           C
HETATM   32 H12R MY  a  39     -16.916 -23.351  -0.954  1.00  0.00           H
HETATM   33 H12S MY  a  39     -17.405 -25.109  -0.650  1.00  0.00           H
HETATM   34 C113 MY  a  39     -16.261 -24.842  -2.377  1.00  0.00           C
HETATM   35 H13R MY  a  39     -17.193 -24.921  -2.976  1.00  0.00           H
HETATM   36 H13S MY  a  39     -15.837 -25.829  -2.094  1.00  0.00           H
HETATM   37 C114 MY  a  39     -15.183 -23.958  -3.049  1.00  0.00           C
HETATM   38 H14R MY  a  39     -15.647 -23.119  -3.611  1.00  0.00           H
HETATM   39 H14S MY  a  39     -14.526 -24.562  -3.710  1.00  0.00           H
HETATM   40 H14T MY  a  39     -14.444 -23.456  -2.389  1.00  0.00           H
HETATM   41  N31 PC  a  39     -12.452 -25.364  17.439  1.00  0.00           N
HETATM   42  C33 PC  a  39     -11.812 -26.657  17.938  1.00  0.00           C
HETATM   43  H3A PC  a  39     -11.530 -27.317  17.131  1.00  0.00           H
HETATM   44  H3B PC  a  39     -12.474 -27.146  18.638  1.00  0.00           H
HETATM   45  H3C PC  a  39     -10.995 -26.370  18.584  1.00  0.00           H
HETATM   46  C34 PC  a  39     -12.762 -24.511  18.609  1.00  0.00           C
HETATM   47  H4A PC  a  39     -12.245 -24.860  19.491  1.00  0.00           H
HETATM   48  H4B PC  a  39     -13.842 -24.538  18.611  1.00  0.00           H
HETATM   49  H4C PC  a  39     -12.450 -23.532  18.277  1.00  0.00           H
HETATM   50  C35 PC  a  39     -11.402 -24.626  16.722  1.00  0.00           C
HETATM   51  H5A PC  a  39     -10.931 -25.262  15.986  1.00  0.00           H
HETATM   52  H5B PC  a  39     -10.653 -24.232  17.392  1.00  0.00           H
HETATM   53  H5C PC  a  39     -11.934 -23.783  16.307  1.00  0.00           H
HETATM   54  C32 PC  a  39     -13.662 -25.757  16.559  1.00  0.00           C
HETATM   55  H2A PC  a  39     -14.338 -26.177  17.290  1.00  0.00           H
HETATM   56  H2B PC  a  39     -13.416 -26.585  15.910  1.00  0.00           H
HETATM   57  C31 PC  a  39     -14.236 -24.742  15.529  1.00  0.00           C
HETATM   58  H1A PC  a  39     -15.266 -25.070  15.272  1.00  0.00           H
HETATM   59  H1B PC  a  39     -13.663 -24.652  14.582  1.00  0.00           H
HETATM   60  P31 PC  a  39     -15.688 -22.889  16.951  1.00  0.00           P
HETATM   61  O33 PC  a  39     -15.489 -21.415  17.244  1.00  0.00           O
HETATM   62  O34 PC  a  39     -16.084 -23.902  17.991  1.00  0.00           O
HETATM   63  O31 PC  a  39     -16.751 -22.812  15.827  1.00  0.00           O
HETATM   64  O32 PC  a  39     -14.369 -23.451  16.189  1.00  0.00           O
HETATM   65  C3  PC  a  39     -16.623 -22.167  14.596  1.00  0.00           C
HETATM   66  HA  PC  a  39     -15.651 -22.494  14.168  1.00  0.00           H
HETATM   67  HB  PC  a  39     -16.617 -21.056  14.628  1.00  0.00           H
HETATM   68  C2  PC  a  39     -17.772 -22.729  13.661  1.00  0.00           C
HETATM   69  HX  PC  a  39     -18.684 -22.318  14.146  1.00  0.00           H
HETATM   70  O21 PC  a  39     -17.855 -24.149  13.789  1.00  0.00           O
HETATM   71  C21 PC  a  39     -18.809 -24.736  13.067  1.00  0.00           C
HETATM   72  O22 PC  a  39     -19.792 -24.142  12.635  1.00  0.00           O
HETATM   73  C1  PC  a  39     -17.523 -22.332  12.147  1.00  0.00           C
HETATM   74  HR  PC  a  39     -17.391 -21.233  12.238  1.00  0.00           H
HETATM   75  HS  PC  a  39     -18.399 -22.460  11.475  1.00  0.00           H
HETATM   76  O11 PC  a  39     -16.332 -22.930  11.651  1.00  0.00           O
HETATM   77  C11 PC  a  39     -15.498 -22.195  10.959  1.00  0.00           C
HETATM   78  O12 PC  a  39     -15.550 -21.004  10.653  1.00  0.00           O
HETATM   79  C12 MY  a  39     -18.526 -26.207  13.020  1.00  0.00           C
HETATM   80  H2R MY  a  39     -19.160 -26.708  13.782  1.00  0.00           H
HETATM   81  H2S MY  a  39     -17.466 -26.532  13.088  1.00  0.00           H
HETATM   82  C13 MY  a  39     -18.924 -26.731  11.648  1.00  0.00           C
HETATM   83  H3R MY  a  39     -20.022 -26.641  11.506  1.00  0.00           H
HETATM   84  H3S MY  a  39     -18.625 -27.800  11.702  1.00  0.00           H
HETATM   85  C14 MY  a  39     -18.194 -26.160  10.387  1.00  0.00           C
HETATM   86  H4R MY  a  39     -17.094 -26.133  10.535  1.00  0.00           H
HETATM   87  H4S MY  a  39     -18.531 -25.107  10.275  1.00  0.00           H
HETATM   88  C15 MY  a  39     -18.343 -27.001   9.110  1.00  0.00           C
HETATM   89  H5R MY  a  39     -19.410 -27.309   9.108  1.00  0.00           H
HETATM   90  H5S MY  a  39     -17.793 -27.954   9.265  1.00  0.00           H
HETATM   91  C16 MY  a  39     -18.039 -26.201   7.771  1.00  0.00           C
HETATM   92  H6R MY  a  39     -16.928 -26.215   7.752  1.00  0.00           H
HETATM   93  H6S MY  a  39     -18.244 -25.115   7.885  1.00  0.00           H
HETATM   94  C17 MY  a  39     -18.600 -26.921   6.541  1.00  0.00           C
HETATM   95  H7R MY  a  39     -19.563 -27.351   6.890  1.00  0.00           H
HETATM   96  H7S MY  a  39     -17.934 -27.744   6.205  1.00  0.00           H
HETATM   97  C18 MY  a  39     -18.860 -25.969   5.359  1.00  0.00           C
HETATM   98  H8R MY  a  39     -17.928 -25.411   5.122  1.00  0.00           H
HETATM   99  H8S MY  a  39     -19.557 -25.168   5.683  1.00  0.00           H
HETATM  100  C19 MY  a  39     -19.375 -26.629   4.110  1.00  0.00           C
HETATM  101  H9R MY  a  39     -20.273 -27.180   4.461  1.00  0.00           H
HETATM  102  H9S MY  a  39     -18.737 -27.372   3.586  1.00  0.00           H
HETATM  103 C110 MY  a  39     -19.860 -25.633   3.085  1.00  0.00           C
HETATM  104 H10R MY  a  39     -19.107 -24.938   2.654  1.00  0.00           H
HETATM  105 H10S MY  a  39     -20.670 -25.076   3.602  1.00  0.00           H
HETATM  106 C111 MY  a  39     -20.438 -26.209   1.749  1.00  0.00           C
HETATM  107 H11R MY  a  39     -21.264 -26.934   1.914  1.00  0.00           H
HETATM  108 H11S MY  a  39     -19.679 -26.664   1.079  1.00  0.00           H
HETATM  109 C112 MY  a  39     -21.195 -25.133   0.961  1.00  0.00           C
HETATM  110 H12R MY  a  39     -21.862 -24.582   1.658  1.00  0.00           H
HETATM  111 H12S MY  a  39     -21.874 -25.552   0.187  1.00  0.00           H
HETATM  112 C113 MY  a  39     -20.320 -24.118   0.222  1.00  0.00           C
HETATM  113 H13R MY  a  39     -19.624 -24.604  -0.495  1.00  0.00           H
HETATM  114 H13S MY  a  39     -19.660 -23.533   0.897  1.00  0.00           H
HETATM  115 C114 MY  a  39     -21.029 -22.967  -0.447  1.00  0.00           C
HETATM  116 H14R MY  a  39     -20.265 -22.312  -0.918  1.00  0.00           H
HETATM  117 H14S MY  a  39     -21.703 -22.364   0.199  1.00  0.00           H
HETATM  118 H14T MY  a  39     -21.613 -23.447  -1.261  1.00  0.00           H
END}

	set f [open "DMPC_A.pdb" "w"]
	puts $f "$DMPC"
	close $f

	# DOPC

	set DOPC {HETATM    1  C12 OL  a   1      -9.326 -10.224  15.460  1.00  0.00           C
HETATM    2  H2R OL  a   1      -9.322  -9.321  16.108  1.00  0.00           H
HETATM    3  H2S OL  a   1      -8.285 -10.237  15.071  1.00  0.00           H
HETATM    4  C13 OL  a   1     -10.340  -9.864  14.369  1.00  0.00           C
HETATM    5  H3R OL  a   1     -11.344  -9.949  14.838  1.00  0.00           H
HETATM    1  H3S OL  a   1     -10.516 -10.639  13.593  1.00  0.00           H
HETATM    7  C14 OL  a   1     -10.195  -8.458  13.685  1.00  0.00           C
HETATM    8  H4R OL  a   1     -10.024  -7.636  14.412  1.00  0.00           H
HETATM    9  H4S OL  a   1     -11.091  -8.197  13.081  1.00  0.00           H
HETATM   10  C15 OL  a   1      -8.890  -8.341  12.877  1.00  0.00           C
HETATM   11  H5R OL  a   1      -8.055  -8.658  13.538  1.00  0.00           H
HETATM   12  H5S OL  a   1      -8.795  -7.250  12.689  1.00  0.00           H
HETATM   13  C16 OL  a   1      -8.734  -9.123  11.503  1.00  0.00           C
HETATM   14  H6R OL  a   1      -7.715  -8.923  11.108  1.00  0.00           H
HETATM   15  H6S OL  a   1      -8.616 -10.209  11.707  1.00  0.00           H
HETATM   16  C17 OL  a   1      -9.682  -8.746  10.316  1.00  0.00           C
HETATM   17  H7R OL  a   1     -10.054  -7.703  10.406  1.00  0.00           H
HETATM   18  H7S OL  a   1      -9.134  -8.784   9.350  1.00  0.00           H
HETATM   19  C18 OL  a   1     -10.960  -9.623  10.257  1.00  0.00           C
HETATM   20  H8R OL  a   1     -11.768  -9.345   9.547  1.00  0.00           H
HETATM   21  H8S OL  a   1     -11.435  -9.566  11.259  1.00  0.00           H
HETATM   22  C19 OL  a   1     -10.487 -10.987   9.903  1.00  0.00           C
HETATM   23  H9R OL  a   1     -10.016 -11.533  10.734  1.00  0.00           H
HETATM   24 C110 OL  a   1     -10.577 -11.608   8.703  1.00  0.00           C
HETATM   25 H10R OL  a   1     -10.135 -12.613   8.769  1.00  0.00           H
HETATM   26 C111 OL  a   1     -11.176 -11.145   7.405  1.00  0.00           C
HETATM   27 H11R OL  a   1     -11.920 -10.328   7.526  1.00  0.00           H
HETATM   28 H11S OL  a   1     -11.790 -11.950   6.948  1.00  0.00           H
HETATM   29 C112 OL  a   1     -10.115 -10.796   6.325  1.00  0.00           C
HETATM   30 H12R OL  a   1      -9.496 -11.705   6.167  1.00  0.00           H
HETATM   31 H12S OL  a   1      -9.557  -9.928   6.736  1.00  0.00           H
HETATM   32 C113 OL  a   1     -10.849 -10.446   5.004  1.00  0.00           C
HETATM   33 H13R OL  a   1     -11.564 -11.205   4.619  1.00  0.00           H
HETATM   34 H13S OL  a   1     -11.290  -9.452   5.228  1.00  0.00           H
HETATM   35 C114 OL  a   1      -9.679 -10.204   3.946  1.00  0.00           C
HETATM   36 H14R OL  a   1      -8.897  -9.605   4.461  1.00  0.00           H
HETATM   37 H14S OL  a   1      -9.253 -11.202   3.708  1.00  0.00           H
HETATM   38 C115 OL  a   1     -10.269  -9.522   2.700  1.00  0.00           C
HETATM   39 H15R OL  a   1     -10.491  -8.475   2.997  1.00  0.00           H
HETATM   40 H15S OL  a   1     -11.206 -10.014   2.362  1.00  0.00           H
HETATM   41 C116 OL  a   1      -9.279  -9.503   1.587  1.00  0.00           C
HETATM   42 H16R OL  a   1      -8.764 -10.473   1.423  1.00  0.00           H
HETATM   43 H16S OL  a   1      -9.736  -9.207   0.619  1.00  0.00           H
HETATM   44 C117 OL  a   1      -8.154  -8.500   1.824  1.00  0.00           C
HETATM   45 H17R OL  a   1      -8.409  -7.617   2.449  1.00  0.00           H
HETATM   46 H17S OL  a   1      -7.413  -8.909   2.544  1.00  0.00           H
HETATM   47 C118 OL  a   1      -7.458  -8.036   0.512  1.00  0.00           C
HETATM   48 H18R OL  a   1      -8.191  -8.026  -0.322  1.00  0.00           H
HETATM   49 H18S OL  a   1      -6.913  -7.080   0.665  1.00  0.00           H
HETATM   50 H18T OL  a   1      -6.656  -8.800   0.424  1.00  0.00           H
HETATM   51  N31 PC  a   1      -6.184  -7.750  17.865  1.00  0.00           N
HETATM   52  C32 PC  a   1      -4.867  -7.862  18.684  1.00  0.00           C
HETATM   53  H2A PC  a   1      -4.484  -6.870  18.873  1.00  0.00           H
HETATM   54  H2B PC  a   1      -4.042  -8.317  18.156  1.00  0.00           H
HETATM   55  C33 PC  a   1      -5.973  -6.592  16.907  1.00  0.00           C
HETATM   56  H3A PC  a   1      -5.632  -5.745  17.484  1.00  0.00           H
HETATM   57  H3B PC  a   1      -5.231  -6.961  16.215  1.00  0.00           H
HETATM   58  H3C PC  a   1      -6.829  -6.455  16.263  1.00  0.00           H
HETATM   59  C34 PC  a   1      -7.331  -7.388  18.747  1.00  0.00           C
HETATM   60  H4A PC  a   1      -8.146  -7.130  18.087  1.00  0.00           H
HETATM   61  H4B PC  a   1      -7.641  -8.234  19.344  1.00  0.00           H
HETATM   62  H4C PC  a   1      -6.972  -6.676  19.476  1.00  0.00           H
HETATM   63  C35 PC  a   1      -6.382  -8.974  17.004  1.00  0.00           C
HETATM   64  H5A PC  a   1      -6.287  -9.834  17.651  1.00  0.00           H
HETATM   65  H5B PC  a   1      -5.616  -8.850  16.254  1.00  0.00           H
HETATM   66  H5C PC  a   1      -7.339  -9.002  16.503  1.00  0.00           H
HETATM   67  C31 PC  a   1      -4.901  -8.707  19.962  1.00  0.00           C
HETATM   68  H1A PC  a   1      -3.835  -8.574  20.248  1.00  0.00           H
HETATM   69  H1B PC  a   1      -5.538  -8.389  20.815  1.00  0.00           H
HETATM   70  P31 PC  a   1      -4.445 -11.321  19.967  1.00  0.00           P
HETATM   71  O33 PC  a   1      -4.160 -11.217  21.421  1.00  0.00           O
HETATM   72  O34 PC  a   1      -3.277 -11.380  19.079  1.00  0.00           O
HETATM   73  O32 PC  a   1      -5.335 -10.025  19.624  1.00  0.00           O
HETATM   74  O31 PC  a   1      -5.401 -12.594  19.739  1.00  0.00           O
HETATM   75  C3  PC  a   1      -6.705 -12.389  19.347  1.00  0.00           C
HETATM   76  HA  PC  a   1      -6.970 -11.315  19.251  1.00  0.00           H
HETATM   77  HB  PC  a   1      -7.313 -12.838  20.162  1.00  0.00           H
HETATM   78  C2  PC  a   1      -6.905 -13.282  18.056  1.00  0.00           C
HETATM   79  HX  PC  a   1      -6.746 -14.361  18.267  1.00  0.00           H
HETATM   80  O21 PC  a   1      -5.900 -12.759  17.113  1.00  0.00           O
HETATM   81  C21 PC  a   1      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   82  O22 PC  a   1      -6.234 -14.582  15.821  1.00  0.00           O
HETATM   83  C1  PC  a   1      -8.359 -13.112  17.525  1.00  0.00           C
HETATM   84  HR  PC  a   1      -9.050 -13.396  18.347  1.00  0.00           H
HETATM   85  HS  PC  a   1      -8.565 -13.815  16.689  1.00  0.00           H
HETATM   86  O11 PC  a   1      -8.573 -11.766  17.067  1.00  0.00           O
HETATM   87  C11 PC  a   1      -9.675 -11.555  16.259  1.00  0.00           C
HETATM   88  O12 PC  a   1     -10.638 -12.267  16.242  1.00  0.00           O
HETATM   89  C12 OL  a   1      -5.038 -12.692  14.886  1.00  0.00           C
HETATM   90  H2R OL  a   1      -3.962 -12.929  15.034  1.00  0.00           H
HETATM   91  H2S OL  a   1      -5.223 -11.604  15.010  1.00  0.00           H
HETATM   92  C13 OL  a   1      -5.565 -13.291  13.565  1.00  0.00           C
HETATM   93  H3R OL  a   1      -6.674 -13.323  13.501  1.00  0.00           H
HETATM   94  H3S OL  a   1      -5.106 -14.290  13.402  1.00  0.00           H
HETATM   95  C14 OL  a   1      -5.066 -12.562  12.388  1.00  0.00           C
HETATM   96  H4R OL  a   1      -5.609 -12.965  11.507  1.00  0.00           H
HETATM   97  H4S OL  a   1      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   98  C15 OL  a   1      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   99  H5R OL  a   1      -4.722 -10.479  12.952  1.00  0.00           H
HETATM  100  H5S OL  a   1      -6.378 -10.838  12.716  1.00  0.00           H
HETATM  101  C16 OL  a   1      -5.328 -10.465  10.922  1.00  0.00           C
HETATM  102  H6R OL  a   1      -5.283  -9.355  10.939  1.00  0.00           H
HETATM  103  H6S OL  a   1      -6.083 -10.781  10.172  1.00  0.00           H
HETATM  104  C17 OL  a   1      -4.072 -10.915  10.088  1.00  0.00           C
HETATM  105  H7R OL  a   1      -3.803 -11.985  10.224  1.00  0.00           H
HETATM  106  H7S OL  a   1      -3.228 -10.288  10.447  1.00  0.00           H
HETATM  107  C18 OL  a   1      -4.096 -10.707   8.561  1.00  0.00           C
HETATM  108  H8R OL  a   1      -3.023 -10.549   8.317  1.00  0.00           H
HETATM  109  H8S OL  a   1      -4.539  -9.771   8.158  1.00  0.00           H
HETATM  110  C19 OL  a   1      -4.543 -11.872   7.693  1.00  0.00           C
HETATM  111  H9R OL  a   1      -4.056 -12.833   7.917  1.00  0.00           H
HETATM  112 C110 OL  a   1      -5.517 -11.984   6.700  1.00  0.00           C
HETATM  113 H10R OL  a   1      -5.665 -12.996   6.296  1.00  0.00           H
HETATM  114 C111 OL  a   1      -6.391 -10.916   6.211  1.00  0.00           C
HETATM  115 H11R OL  a   1      -6.147  -9.972   6.743  1.00  0.00           H
HETATM  116 H11S OL  a   1      -7.457 -11.176   6.387  1.00  0.00           H
HETATM  117 C112 OL  a   1      -6.043 -10.766   4.699  1.00  0.00           C
HETATM  118 H12R OL  a   1      -4.989 -10.419   4.634  1.00  0.00           H
HETATM  119 H12S OL  a   1      -6.774 -10.059   4.253  1.00  0.00           H
HETATM  120 C113 OL  a   1      -6.272 -11.989   3.860  1.00  0.00           C
HETATM  121 H13R OL  a   1      -7.325 -12.174   4.161  1.00  0.00           H
HETATM  122 H13S OL  a   1      -5.670 -12.857   4.208  1.00  0.00           H
HETATM  123 C114 OL  a   1      -5.995 -11.819   2.357  1.00  0.00           C
HETATM  124 H14R OL  a   1      -4.903 -11.898   2.168  1.00  0.00           H
HETATM  125 H14S OL  a   1      -6.282 -10.759   2.192  1.00  0.00           H
HETATM  126 C115 OL  a   1      -6.855 -12.761   1.560  1.00  0.00           C
HETATM  127 H15R OL  a   1      -7.917 -12.453   1.672  1.00  0.00           H
HETATM  128 H15S OL  a   1      -6.820 -13.811   1.919  1.00  0.00           H
HETATM  129 C116 OL  a   1      -6.528 -12.913   0.035  1.00  0.00           C
HETATM  130 H16R OL  a   1      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM  131 H16S OL  a   1      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM  132 C117 OL  a   1      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM  133 H17R OL  a   1      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM  134 H17S OL  a   1      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM  135 C118 OL  a   1      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM  136 H18R OL  a   1      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM  137 H18S OL  a   1      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM  138 H18T OL  a   1      -5.769 -10.553  -2.572  1.00  0.00           H
END}
	
	set f [open "DOPC_A.pdb" "w"]
	puts $f "$DOPC"
	close $f

	# DOPG

	set DOPG {HETATM    1  C12 OL  a   1      -9.326 -10.224  15.460  1.00  0.00           C
HETATM    2  H2R OL  a   1      -9.322  -9.321  16.108  1.00  0.00           H
HETATM    3  H2S OL  a   1      -8.285 -10.237  15.071  1.00  0.00           H
HETATM    4  C13 OL  a   1     -10.340  -9.864  14.369  1.00  0.00           C
HETATM    5  H3R OL  a   1     -11.344  -9.949  14.838  1.00  0.00           H
HETATM    6  H3S OL  a   1     -10.516 -10.639  13.593  1.00  0.00           H
HETATM    7  C14 OL  a   1     -10.195  -8.458  13.685  1.00  0.00           C
HETATM    8  H4R OL  a   1     -10.024  -7.636  14.412  1.00  0.00           H
HETATM    9  H4S OL  a   1     -11.091  -8.197  13.081  1.00  0.00           H
HETATM   10  C15 OL  a   1      -8.890  -8.341  12.877  1.00  0.00           C
HETATM   11  H5R OL  a   1      -8.055  -8.658  13.538  1.00  0.00           H
HETATM   12  H5S OL  a   1      -8.795  -7.250  12.689  1.00  0.00           H
HETATM   13  C16 OL  a   1      -8.734  -9.123  11.503  1.00  0.00           C
HETATM   14  H6R OL  a   1      -7.715  -8.923  11.108  1.00  0.00           H
HETATM   15  H6S OL  a   1      -8.616 -10.209  11.707  1.00  0.00           H
HETATM   16  C17 OL  a   1      -9.682  -8.746  10.316  1.00  0.00           C
HETATM   17  H7R OL  a   1     -10.054  -7.703  10.406  1.00  0.00           H
HETATM   18  H7S OL  a   1      -9.134  -8.784   9.350  1.00  0.00           H
HETATM   19  C18 OL  a   1     -10.960  -9.623  10.257  1.00  0.00           C
HETATM   20  H8R OL  a   1     -11.768  -9.345   9.547  1.00  0.00           H
HETATM   21  H8S OL  a   1     -11.435  -9.566  11.259  1.00  0.00           H
HETATM   22  C19 OL  a   1     -10.487 -10.987   9.903  1.00  0.00           C
HETATM   23  H9R OL  a   1     -10.016 -11.533  10.734  1.00  0.00           H
HETATM   24 C110 OL  a   1     -10.577 -11.608   8.703  1.00  0.00           C
HETATM   25 H10R OL  a   1     -10.135 -12.613   8.769  1.00  0.00           H
HETATM   26 C111 OL  a   1     -11.176 -11.145   7.405  1.00  0.00           C
HETATM   27 H11R OL  a   1     -11.920 -10.328   7.526  1.00  0.00           H
HETATM   28 H11S OL  a   1     -11.790 -11.950   6.948  1.00  0.00           H
HETATM   29 C112 OL  a   1     -10.115 -10.796   6.325  1.00  0.00           C
HETATM   30 H12R OL  a   1      -9.496 -11.705   6.167  1.00  0.00           H
HETATM   31 H12S OL  a   1      -9.557  -9.928   6.736  1.00  0.00           H
HETATM   32 C113 OL  a   1     -10.849 -10.446   5.004  1.00  0.00           C
HETATM   33 H13R OL  a   1     -11.564 -11.205   4.619  1.00  0.00           H
HETATM   34 H13S OL  a   1     -11.290  -9.452   5.228  1.00  0.00           H
HETATM   35 C114 OL  a   1      -9.679 -10.204   3.946  1.00  0.00           C
HETATM   36 H14R OL  a   1      -8.897  -9.605   4.461  1.00  0.00           H
HETATM   37 H14S OL  a   1      -9.253 -11.202   3.708  1.00  0.00           H
HETATM   38 C115 OL  a   1     -10.269  -9.522   2.700  1.00  0.00           C
HETATM   39 H15R OL  a   1     -10.491  -8.475   2.997  1.00  0.00           H
HETATM   40 H15S OL  a   1     -11.206 -10.014   2.362  1.00  0.00           H
HETATM   41 C116 OL  a   1      -9.279  -9.503   1.587  1.00  0.00           C
HETATM   42 H16R OL  a   1      -8.764 -10.473   1.423  1.00  0.00           H
HETATM   43 H16S OL  a   1      -9.736  -9.207   0.619  1.00  0.00           H
HETATM   44 C117 OL  a   1      -8.154  -8.500   1.824  1.00  0.00           C
HETATM   45 H17R OL  a   1      -8.409  -7.617   2.449  1.00  0.00           H
HETATM   46 H17S OL  a   1      -7.413  -8.909   2.544  1.00  0.00           H
HETATM   47 C118 OL  a   1      -7.458  -8.036   0.512  1.00  0.00           C
HETATM   48 H18R OL  a   1      -8.191  -8.026  -0.322  1.00  0.00           H
HETATM   49 H18S OL  a   1      -6.913  -7.080   0.665  1.00  0.00           H
HETATM   50 H18T OL  a   1      -6.656  -8.800   0.424  1.00  0.00           H
HETATM   51  C33 PG  a   1      -7.424  -9.612  21.287  1.00  0.00           C
HETATM   52  H3A PG  a   1      -7.571 -10.051  20.277  1.00  0.00           H
HETATM   53  H3B PG  a   1      -6.837 -10.366  21.854  1.00  0.00           H
HETATM   54  O36 PG  a   1      -8.615  -9.495  22.014  1.00  0.00           O
HETATM   55 HO6A PG  a   1      -9.225  -9.882  21.381  1.00  0.00           H
HETATM   56  C32 PG  a   1      -6.682  -8.201  21.140  1.00  0.00           C
HETATM   57  H2A PG  a   1      -6.846  -7.649  22.091  1.00  0.00           H
HETATM   58  O35 PG  a   1      -5.261  -8.330  21.137  1.00  0.00           O
HETATM   59 HO5A PG  a   1      -4.979  -7.482  21.488  1.00  0.00           H
HETATM   60  C31 PG  a   1      -4.901  -8.707  19.962  1.00  0.00           C
HETATM   61  H1A PG  a   1      -3.835  -8.574  20.248  1.00  0.00           H
HETATM   62  H1B PG  a   1      -5.538  -8.389  20.815  1.00  0.00           H
HETATM   63  P31 PG  a   1      -4.445 -11.321  19.967  1.00  0.00           P
HETATM   64  O33 PG  a   1      -4.160 -11.217  21.421  1.00  0.00           O
HETATM   65  O34 PG  a   1      -3.277 -11.380  19.079  1.00  0.00           O
HETATM   66  O32 PG  a   1      -5.335 -10.025  19.624  1.00  0.00           O
HETATM   67  O31 PG  a   1      -5.401 -12.594  19.739  1.00  0.00           O
HETATM   68  C3  PG  a   1      -6.705 -12.389  19.347  1.00  0.00           C
HETATM   69  HA  PG  a   1      -6.970 -11.315  19.251  1.00  0.00           H
HETATM   70  HB  PG  a   1      -7.313 -12.838  20.162  1.00  0.00           H
HETATM   71  C2  PG  a   1      -6.905 -13.282  18.056  1.00  0.00           C
HETATM   72  HX  PG  a   1      -6.746 -14.361  18.267  1.00  0.00           H
HETATM   73  O21 PG  a   1      -5.900 -12.759  17.113  1.00  0.00           O
HETATM   74  C21 PG  a   1      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   75  O22 PG  a   1      -6.234 -14.582  15.821  1.00  0.00           O
HETATM   76  C1  PG  a   1      -8.359 -13.112  17.525  1.00  0.00           C
HETATM   77  HR  PG  a   1      -9.050 -13.396  18.347  1.00  0.00           H
HETATM   78  HS  PG  a   1      -8.565 -13.815  16.689  1.00  0.00           H
HETATM   79  O11 PG  a   1      -8.573 -11.766  17.067  1.00  0.00           O
HETATM   80  C11 PG  a   1      -9.675 -11.555  16.259  1.00  0.00           C
HETATM   81  O12 PG  a   1     -10.638 -12.267  16.242  1.00  0.00           O
HETATM   82  C12 OL  a   1      -5.038 -12.692  14.886  1.00  0.00           C
HETATM   83  H2R OL  a   1      -3.962 -12.929  15.034  1.00  0.00           H
HETATM   84  H2S OL  a   1      -5.223 -11.604  15.010  1.00  0.00           H
HETATM   85  C13 OL  a   1      -5.565 -13.291  13.565  1.00  0.00           C
HETATM   86  H3R OL  a   1      -6.674 -13.323  13.501  1.00  0.00           H
HETATM   87  H3S OL  a   1      -5.106 -14.290  13.402  1.00  0.00           H
HETATM   88  C14 OL  a   1      -5.066 -12.562  12.388  1.00  0.00           C
HETATM   89  H4R OL  a   1      -5.609 -12.965  11.507  1.00  0.00           H
HETATM   90  H4S OL  a   1      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   91  C15 OL  a   1      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   92  H5R OL  a   1      -4.722 -10.479  12.952  1.00  0.00           H
HETATM   93  H5S OL  a   1      -6.378 -10.838  12.716  1.00  0.00           H
HETATM   94  C16 OL  a   1      -5.328 -10.465  10.922  1.00  0.00           C
HETATM   95  H6R OL  a   1      -5.283  -9.355  10.939  1.00  0.00           H
HETATM   96  H6S OL  a   1      -6.083 -10.781  10.172  1.00  0.00           H
HETATM   97  C17 OL  a   1      -4.072 -10.915  10.088  1.00  0.00           C
HETATM   98  H7R OL  a   1      -3.803 -11.985  10.224  1.00  0.00           H
HETATM   99  H7S OL  a   1      -3.228 -10.288  10.447  1.00  0.00           H
HETATM  100  C18 OL  a   1      -4.096 -10.707   8.561  1.00  0.00           C
HETATM  101  H8R OL  a   1      -3.023 -10.549   8.317  1.00  0.00           H
HETATM  102  H8S OL  a   1      -4.539  -9.771   8.158  1.00  0.00           H
HETATM  103  C19 OL  a   1      -4.543 -11.872   7.693  1.00  0.00           C
HETATM  104  H9R OL  a   1      -4.056 -12.833   7.917  1.00  0.00           H
HETATM  105 C110 OL  a   1      -5.517 -11.984   6.700  1.00  0.00           C
HETATM  106 H10R OL  a   1      -5.665 -12.996   6.296  1.00  0.00           H
HETATM  107 C111 OL  a   1      -6.391 -10.916   6.211  1.00  0.00           C
HETATM  108 H11R OL  a   1      -6.147  -9.972   6.743  1.00  0.00           H
HETATM  109 H11S OL  a   1      -7.457 -11.176   6.387  1.00  0.00           H
HETATM  110 C112 OL  a   1      -6.043 -10.766   4.699  1.00  0.00           C
HETATM  111 H12R OL  a   1      -4.989 -10.419   4.634  1.00  0.00           H
HETATM  112 H12S OL  a   1      -6.774 -10.059   4.253  1.00  0.00           H
HETATM  113 C113 OL  a   1      -6.272 -11.989   3.860  1.00  0.00           C
HETATM  114 H13R OL  a   1      -7.325 -12.174   4.161  1.00  0.00           H
HETATM  115 H13S OL  a   1      -5.670 -12.857   4.208  1.00  0.00           H
HETATM  116 C114 OL  a   1      -5.995 -11.819   2.357  1.00  0.00           C
HETATM  117 H14R OL  a   1      -4.903 -11.898   2.168  1.00  0.00           H
HETATM  118 H14S OL  a   1      -6.282 -10.759   2.192  1.00  0.00           H
HETATM  119 C115 OL  a   1      -6.855 -12.761   1.560  1.00  0.00           C
HETATM  120 H15R OL  a   1      -7.917 -12.453   1.672  1.00  0.00           H
HETATM  121 H15S OL  a   1      -6.820 -13.811   1.919  1.00  0.00           H
HETATM  122 C116 OL  a   1      -6.528 -12.913   0.035  1.00  0.00           C
HETATM  123 H16R OL  a   1      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM  124 H16S OL  a   1      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM  125 C117 OL  a   1      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM  126 H17R OL  a   1      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM  127 H17S OL  a   1      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM  128 C118 OL  a   1      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM  129 H18R OL  a   1      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM  130 H18S OL  a   1      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM  131 H18T OL  a   1      -5.769 -10.553  -2.572  1.00  0.00           H
END}

	set f [open "DOPG_A.pdb" "w"]
	puts $f "$DOPG"
	close $f

	# DOPS

	set DOPS {HETATM    1  C12 OL  a   1      -6.858   0.132  16.646  1.00  0.00           C
HETATM    2  H2R OL  a   1      -6.888   0.769  17.556  1.00  0.00           H
HETATM    3  H2S OL  a   1      -5.825  -0.187  16.389  1.00  0.00           H
HETATM    4  C13 OL  a   1      -7.562   0.928  15.493  1.00  0.00           C
HETATM    5  H3R OL  a   1      -6.987   1.879  15.477  1.00  0.00           H
HETATM    6  H3S OL  a   1      -8.586   1.164  15.852  1.00  0.00           H
HETATM    7  C14 OL  a   1      -7.423   0.265  14.080  1.00  0.00           C
HETATM    8  H4R OL  a   1      -6.332   0.219  13.876  1.00  0.00           H
HETATM    9  H4S OL  a   1      -7.783  -0.784  14.015  1.00  0.00           H
HETATM   10  C15 OL  a   1      -8.065   1.114  12.979  1.00  0.00           C
HETATM   11  H5R OL  a   1      -7.564   2.105  12.962  1.00  0.00           H
HETATM   12  H5S OL  a   1      -9.143   1.224  13.226  1.00  0.00           H
HETATM   13  C16 OL  a   1      -7.893   0.479  11.608  1.00  0.00           C
HETATM   14  H6R OL  a   1      -8.269   1.162  10.817  1.00  0.00           H
HETATM   15  H6S OL  a   1      -6.824   0.333  11.346  1.00  0.00           H
HETATM   16  C17 OL  a   1      -8.758  -0.745  11.358  1.00  0.00           C
HETATM   17  H7R OL  a   1      -8.498  -1.466  12.163  1.00  0.00           H
HETATM   18  H7S OL  a   1      -9.796  -0.434  11.607  1.00  0.00           H
HETATM   19  C18 OL  a   1      -8.854  -1.363   9.941  1.00  0.00           C
HETATM   20  H8R OL  a   1      -9.651  -2.134   9.866  1.00  0.00           H
HETATM   21  H8S OL  a   1      -9.199  -0.577   9.237  1.00  0.00           H
HETATM   22  C19 OL  a   1      -7.480  -1.825   9.539  1.00  0.00           C
HETATM   23  H9R OL  a   1      -6.780  -2.286  10.252  1.00  0.00           H
HETATM   24 C110 OL  a   1      -7.026  -1.598   8.290  1.00  0.00           C
HETATM   25 H10R OL  a   1      -6.034  -1.930   7.948  1.00  0.00           H
HETATM   26 C111 OL  a   1      -7.682  -0.927   7.106  1.00  0.00           C
HETATM   27 H11R OL  a   1      -6.909  -0.455   6.463  1.00  0.00           H
HETATM   28 H11S OL  a   1      -8.341  -0.046   7.262  1.00  0.00           H
HETATM   29 C112 OL  a   1      -8.471  -1.871   6.179  1.00  0.00           C
HETATM   30 H12R OL  a   1      -8.997  -1.220   5.448  1.00  0.00           H
HETATM   31 H12S OL  a   1      -9.218  -2.408   6.802  1.00  0.00           H
HETATM   32 C113 OL  a   1      -7.601  -2.823   5.324  1.00  0.00           C
HETATM   33 H13R OL  a   1      -7.203  -3.596   6.016  1.00  0.00           H
HETATM   34 H13S OL  a   1      -6.716  -2.287   4.918  1.00  0.00           H
HETATM   35 C114 OL  a   1      -8.378  -3.708   4.310  1.00  0.00           C
HETATM   36 H14R OL  a   1      -8.716  -3.036   3.492  1.00  0.00           H
HETATM   37 H14S OL  a   1      -9.128  -4.359   4.808  1.00  0.00           H
HETATM   38 C115 OL  a   1      -7.512  -4.739   3.658  1.00  0.00           C
HETATM   39 H15R OL  a   1      -8.134  -5.302   2.930  1.00  0.00           H
HETATM   40 H15S OL  a   1      -7.213  -5.393   4.505  1.00  0.00           H
HETATM   41 C116 OL  a   1      -6.203  -4.200   3.130  1.00  0.00           C
HETATM   42 H16R OL  a   1      -5.562  -3.924   3.995  1.00  0.00           H
HETATM   43 H16S OL  a   1      -6.391  -3.367   2.419  1.00  0.00           H
HETATM   44 C117 OL  a   1      -5.336  -5.295   2.619  1.00  0.00           C
HETATM   45 H17R OL  a   1      -5.485  -6.195   3.254  1.00  0.00           H
HETATM   46 H17S OL  a   1      -4.280  -4.959   2.702  1.00  0.00           H
HETATM   47 C118 OL  a   1      -5.619  -5.708   1.154  1.00  0.00           C
HETATM   48 H18R OL  a   1      -5.251  -4.944   0.437  1.00  0.00           H
HETATM   49 H18S OL  a   1      -5.192  -6.728   1.045  1.00  0.00           H
HETATM   50 H18T OL  a   1      -6.684  -5.937   0.934  1.00  0.00           H
HETATM   51  N31 PS  a   1      -2.846  -0.084  21.117  1.00  0.00           N
HETATM   52 HN1A PS  a   1      -3.412   0.414  21.833  1.00  0.00           H
HETATM   53 HN1B PS  a   1      -2.638  -1.077  21.343  1.00  0.00           H
HETATM   54 HN1C PS  a   1      -1.939   0.423  21.168  1.00  0.00           H
HETATM   55  C32 PS  a   1      -3.466   0.177  19.756  1.00  0.00           C
HETATM   56  H2A PS  a   1      -2.673   0.286  19.031  1.00  0.00           H
HETATM   57  C33 PS  a   1      -4.278   1.497  19.725  1.00  0.00           C
HETATM   58  O35 PS  a   1      -4.627   2.132  20.745  1.00  0.00           O
HETATM   59  O36 PS  a   1      -4.579   1.862  18.532  1.00  0.00           O
HETATM   60  C31 PS  a   1      -4.329  -1.061  19.401  1.00  0.00           C
HETATM   61  H1A PS  a   1      -4.872  -0.723  18.493  1.00  0.00           H
HETATM   62  H1B PS  a   1      -5.104  -1.290  20.164  1.00  0.00           H
HETATM   63  P31 PS  a   1      -3.229  -2.850  17.912  1.00  0.00           P
HETATM   64  O33 PS  a   1      -2.422  -4.067  18.138  1.00  0.00           O
HETATM   65  O34 PS  a   1      -2.596  -1.847  17.032  1.00  0.00           O
HETATM   66  O32 PS  a   1      -3.463  -2.228  19.293  1.00  0.00           O
HETATM   67  O31 PS  a   1      -4.650  -3.288  17.349  1.00  0.00           O
HETATM   68  C3  PS  a   1      -5.422  -4.229  18.076  1.00  0.00           C
HETATM   69  HA  PS  a   1      -5.409  -3.963  19.155  1.00  0.00           H
HETATM   70  HB  PS  a   1      -5.134  -5.251  17.749  1.00  0.00           H
HETATM   71  C2  PS  a   1      -6.924  -4.137  17.656  1.00  0.00           C
HETATM   72  HX  PS  a   1      -7.371  -5.091  18.009  1.00  0.00           H
HETATM   73  O21 PS  a   1      -7.011  -3.982  16.268  1.00  0.00           O
HETATM   74  C21 PS  a   1      -6.512  -4.955  15.417  1.00  0.00           C
HETATM   75  O22 PS  a   1      -6.356  -6.114  15.733  1.00  0.00           O
HETATM   76  C1  PS  a   1      -7.704  -2.947  18.282  1.00  0.00           C
HETATM   77  HR  PS  a   1      -7.468  -3.080  19.359  1.00  0.00           H
HETATM   78  HS  PS  a   1      -8.802  -2.944  18.112  1.00  0.00           H
HETATM   79  O11 PS  a   1      -7.187  -1.725  17.828  1.00  0.00           O
HETATM   80  C11 PS  a   1      -7.762  -0.994  16.926  1.00  0.00           C
HETATM   81  O12 PS  a   1      -8.939  -1.158  16.438  1.00  0.00           O
HETATM   82  C12 OL  a   1      -6.308  -4.268  14.069  1.00  0.00           C
HETATM   83  H2R OL  a   1      -7.166  -3.643  13.741  1.00  0.00           H
HETATM   84  H2S OL  a   1      -6.124  -5.002  13.256  1.00  0.00           H
HETATM   85  C13 OL  a   1      -4.986  -3.524  14.137  1.00  0.00           C
HETATM   86  H3R OL  a   1      -4.276  -4.202  14.658  1.00  0.00           H
HETATM   87  H3S OL  a   1      -5.153  -2.602  14.733  1.00  0.00           H
HETATM   88  C14 OL  a   1      -4.475  -3.117  12.774  1.00  0.00           C
HETATM   89  H4R OL  a   1      -5.151  -2.310  12.417  1.00  0.00           H
HETATM   90  H4S OL  a   1      -4.363  -4.008  12.120  1.00  0.00           H
HETATM   91  C15 OL  a   1      -2.981  -2.524  12.862  1.00  0.00           C
HETATM   92  H5R OL  a   1      -2.247  -3.357  12.878  1.00  0.00           H
HETATM   93  H5S OL  a   1      -2.943  -2.061  13.871  1.00  0.00           H
HETATM   94  C16 OL  a   1      -2.749  -1.414  11.838  1.00  0.00           C
HETATM   95  H6R OL  a   1      -3.603  -0.708  11.928  1.00  0.00           H
HETATM   96  H6S OL  a   1      -2.804  -1.849  10.817  1.00  0.00           H
HETATM   97  C17 OL  a   1      -1.344  -0.722  12.115  1.00  0.00           C
HETATM   98  H7R OL  a   1      -0.627  -1.531  12.371  1.00  0.00           H
HETATM   99  H7S OL  a   1      -1.420   0.129  12.825  1.00  0.00           H
HETATM  100  C18 OL  a   1      -0.832   0.011  10.807  1.00  0.00           C
HETATM  101  H8R OL  a   1      -0.099   0.752  11.190  1.00  0.00           H
HETATM  102  H8S OL  a   1      -1.713   0.544  10.390  1.00  0.00           H
HETATM  103  C19 OL  a   1      -0.175  -0.878   9.700  1.00  0.00           C
HETATM  104  H9R OL  a   1       0.489  -1.675  10.065  1.00  0.00           H
HETATM  105 C110 OL  a   1      -0.412  -0.925   8.296  1.00  0.00           C
HETATM  106 H10R OL  a   1       0.211  -1.630   7.726  1.00  0.00           H
HETATM  107 C111 OL  a   1      -1.348  -0.149   7.462  1.00  0.00           C
HETATM  108 H11R OL  a   1      -0.884   0.205   6.516  1.00  0.00           H
HETATM  109 H11S OL  a   1      -1.562   0.772   8.044  1.00  0.00           H
HETATM  110 C112 OL  a   1      -2.709  -0.700   7.258  1.00  0.00           C
HETATM  111 H12R OL  a   1      -3.338  -0.056   6.606  1.00  0.00           H
HETATM  112 H12S OL  a   1      -3.042  -0.903   8.298  1.00  0.00           H
HETATM  113 C113 OL  a   1      -2.899  -2.127   6.683  1.00  0.00           C
HETATM  114 H13R OL  a   1      -3.977  -2.374   6.795  1.00  0.00           H
HETATM  115 H13S OL  a   1      -2.390  -2.881   7.320  1.00  0.00           H
HETATM  116 C114 OL  a   1      -2.533  -2.303   5.163  1.00  0.00           C
HETATM  117 H14R OL  a   1      -2.958  -1.553   4.463  1.00  0.00           H
HETATM  118 H14S OL  a   1      -3.045  -3.250   4.891  1.00  0.00           H
HETATM  119 C115 OL  a   1      -1.081  -2.448   4.817  1.00  0.00           C
HETATM  120 H15R OL  a   1      -0.526  -3.204   5.414  1.00  0.00           H
HETATM  121 H15S OL  a   1      -0.553  -1.513   5.099  1.00  0.00           H
HETATM  122 C116 OL  a   1      -0.813  -2.496   3.288  1.00  0.00           C
HETATM  123 H16R OL  a   1       0.282  -2.326   3.209  1.00  0.00           H
HETATM  124 H16S OL  a   1      -1.405  -1.695   2.797  1.00  0.00           H
HETATM  125 C117 OL  a   1      -1.124  -3.873   2.687  1.00  0.00           C
HETATM  126 H17R OL  a   1      -2.227  -3.984   2.621  1.00  0.00           H
HETATM  127 H17S OL  a   1      -0.730  -4.634   3.394  1.00  0.00           H
HETATM  128 C118 OL  a   1      -0.462  -4.127   1.251  1.00  0.00           C
HETATM  129 H18R OL  a   1      -1.126  -4.856   0.741  1.00  0.00           H
HETATM  130 H18S OL  a   1       0.539  -4.565   1.449  1.00  0.00           H
HETATM  131 H18T OL  a   1      -0.374  -3.179   0.679  1.00  0.00           H
END}

	set f [open "DOPS_A.pdb" "w"]
	puts $f "$DOPS"
	close $f

	# DPPC

	set DPPC {HETATM    1  C12 PA A    1      -3.196   1.629  11.662  1.00  0.00           C
HETATM    2  H2R PA A    1      -3.568   2.480  11.051  1.00  0.00           H
HETATM    3  H2S PA A    1      -4.098   1.199  12.147  1.00  0.00           H
HETATM    4  C13 PA A    1      -2.525   0.636  10.668  1.00  0.00           C
HETATM    5  H3R PA A    1      -3.371   0.135  10.152  1.00  0.00           H
HETATM    6  H3S PA A    1      -2.067  -0.172  11.278  1.00  0.00           H
HETATM    7  C14 PA A    1      -1.671   1.229   9.652  1.00  0.00           C
HETATM    8  H4R PA A    1      -0.925   1.980   9.989  1.00  0.00           H
HETATM    9  H4S PA A    1      -2.428   1.764   9.040  1.00  0.00           H
HETATM   10  C15 PA A    1      -0.945   0.254   8.813  1.00  0.00           C
HETATM   11  H5R PA A    1      -0.068  -0.206   9.317  1.00  0.00           H
HETATM   12  H5S PA A    1      -0.607   0.789   7.900  1.00  0.00           H
HETATM   13  C16 PA A    1      -1.808  -0.927   8.414  1.00  0.00           C
HETATM   14  H6R PA A    1      -2.714  -0.555   7.889  1.00  0.00           H
HETATM   15  H6S PA A    1      -2.161  -1.437   9.336  1.00  0.00           H
HETATM   16  C17 PA A    1      -1.042  -1.890   7.517  1.00  0.00           C
HETATM   17  H7R PA A    1      -0.334  -2.392   8.211  1.00  0.00           H
HETATM   18  H7S PA A    1      -0.504  -1.401   6.677  1.00  0.00           H
HETATM   19  C18 PA A    1      -1.884  -2.984   6.894  1.00  0.00           C
HETATM   20  H8R PA A    1      -2.631  -3.241   7.675  1.00  0.00           H
HETATM   21  H8S PA A    1      -1.288  -3.914   6.775  1.00  0.00           H
HETATM   22  C19 PA A    1      -2.743  -2.645   5.696  1.00  0.00           C
HETATM   23  H9R PA A    1      -3.210  -1.641   5.782  1.00  0.00           H
HETATM   24  H9S PA A    1      -3.529  -3.387   5.441  1.00  0.00           H
HETATM   25 C110 PA A    1      -1.891  -2.651   4.422  1.00  0.00           C
HETATM   26 H10R PA A    1      -1.738  -3.736   4.237  1.00  0.00           H
HETATM   27 H10S PA A    1      -0.937  -2.138   4.667  1.00  0.00           H
HETATM   28 C111 PA A    1      -2.488  -2.036   3.122  1.00  0.00           C
HETATM   29 H11R PA A    1      -2.384  -0.932   3.178  1.00  0.00           H
HETATM   30 H11S PA A    1      -3.516  -2.359   2.850  1.00  0.00           H
HETATM   31 C112 PA A    1      -1.620  -2.512   1.923  1.00  0.00           C
HETATM   32 H12R PA A    1      -2.373  -2.805   1.160  1.00  0.00           H
HETATM   33 H12S PA A    1      -0.993  -3.408   2.117  1.00  0.00           H
HETATM   34 C113 PA A    1      -0.763  -1.395   1.292  1.00  0.00           C
HETATM   35 H13R PA A    1      -0.098  -1.076   2.123  1.00  0.00           H
HETATM   36 H13S PA A    1      -1.551  -0.661   1.021  1.00  0.00           H
HETATM   37 C114 PA A    1      -0.014  -1.789  -0.012  1.00  0.00           C
HETATM   38 H14R PA A    1       0.707  -2.591   0.254  1.00  0.00           H
HETATM   39 H14S PA A    1       0.546  -0.918  -0.415  1.00  0.00           H
HETATM   40 C115 PA A    1      -0.990  -2.342  -1.120  1.00  0.00           C
HETATM   41 H15R PA A    1      -1.838  -1.642  -0.960  1.00  0.00           H
HETATM   42 H15S PA A    1      -1.179  -3.401  -0.845  1.00  0.00           H
HETATM   43 C116 PA A    1      -0.502  -2.324  -2.594  1.00  0.00           C
HETATM   44 H16R PA A    1      -1.391  -2.396  -3.256  1.00  0.00           H
HETATM   45 H16S PA A    1       0.193  -3.179  -2.735  1.00  0.00           H
HETATM   46 H16T PA A    1       0.081  -1.387  -2.718  1.00  0.00           H
HETATM   47  N31 PC A    1       0.025   0.157  20.702  1.00  0.00           N
HETATM   48  C33 PC A    1      -1.081   0.096  19.689  1.00  0.00           C
HETATM   49  H3A PC A    1      -1.513   1.080  19.583  1.00  0.00           H
HETATM   50  H3B PC A    1      -1.860  -0.580  20.010  1.00  0.00           H
HETATM   51  H3C PC A    1      -0.662  -0.267  18.762  1.00  0.00           H
HETATM   52  C34 PC A    1       0.830   1.428  20.386  1.00  0.00           C
HETATM   53  H4A PC A    1       0.144   2.228  20.149  1.00  0.00           H
HETATM   54  H4B PC A    1       1.427   1.260  19.502  1.00  0.00           H
HETATM   55  H4C PC A    1       1.401   1.613  21.285  1.00  0.00           H
HETATM   56  C35 PC A    1      -0.594   0.257  22.082  1.00  0.00           C
HETATM   57  H5A PC A    1      -1.104   1.207  22.140  1.00  0.00           H
HETATM   58  H5B PC A    1       0.197   0.333  22.814  1.00  0.00           H
HETATM   59  H5C PC A    1      -1.254  -0.562  22.328  1.00  0.00           H
HETATM   60  C32 PC A    1       0.937  -1.062  20.647  1.00  0.00           C
HETATM   61  H2A PC A    1       0.480  -1.913  21.130  1.00  0.00           H
HETATM   62  H2B PC A    1       1.800  -0.834  21.256  1.00  0.00           H
HETATM   63  C31 PC A    1       1.422  -1.670  19.299  1.00  0.00           C
HETATM   64  H1A PC A    1       0.554  -1.937  18.660  1.00  0.00           H
HETATM   65  H1B PC A    1       2.018  -2.576  19.537  1.00  0.00           H
HETATM   66  P31 PC A    1       2.730  -0.906  17.207  1.00  0.00           P
HETATM   67  O33 PC A    1       3.342  -2.234  17.104  1.00  0.00           O
HETATM   68  O34 PC A    1       3.531   0.284  16.845  1.00  0.00           O
HETATM   69  O31 PC A    1       1.406  -1.009  16.307  1.00  0.00           O
HETATM   70  O32 PC A    1       2.332  -0.693  18.731  1.00  0.00           O
HETATM   71  C3  PC A    1       0.574   0.107  16.200  1.00  0.00           C
HETATM   72  HA  PC A    1       1.069   1.037  15.848  1.00  0.00           H
HETATM   73  HB  PC A    1       0.069   0.273  17.175  1.00  0.00           H
HETATM   74  C2  PC A    1      -0.461  -0.104  15.129  1.00  0.00           C
HETATM   75  HX  PC A    1      -1.379  -0.508  15.607  1.00  0.00           H
HETATM   76  O21 PC A    1      -0.041  -1.024  14.085  1.00  0.00           O
HETATM   77  C21 PC A    1      -0.127  -2.342  14.264  1.00  0.00           C
HETATM   78  O22 PC A    1      -0.616  -2.931  15.228  1.00  0.00           O
HETATM   79  C1  PC A    1      -0.855   1.229  14.495  1.00  0.00           C
HETATM   80  HR  PC A    1      -0.024   1.758  13.983  1.00  0.00           H
HETATM   81  HS  PC A    1      -1.277   1.821  15.335  1.00  0.00           H
HETATM   82  O11 PC A    1      -1.857   0.986  13.482  1.00  0.00           O
HETATM   83  C11 PC A    1      -2.194   2.025  12.783  1.00  0.00           C
HETATM   84  O12 PC A    1      -1.820   3.209  12.898  1.00  0.00           O
HETATM   85  C12 PA A    1       0.500  -3.074  13.101  1.00  0.00           C
HETATM   86  H2R PA A    1       0.119  -2.615  12.163  1.00  0.00           H
HETATM   87  H2S PA A    1       0.290  -4.164  13.110  1.00  0.00           H
HETATM   88  C13 PA A    1       2.003  -2.927  13.164  1.00  0.00           C
HETATM   89  H3R PA A    1       2.498  -3.368  14.056  1.00  0.00           H
HETATM   90  H3S PA A    1       2.107  -1.821  13.183  1.00  0.00           H
HETATM   91  C14 PA A    1       2.865  -3.338  11.966  1.00  0.00           C
HETATM   92  H4R PA A    1       2.738  -4.438  11.873  1.00  0.00           H
HETATM   93  H4S PA A    1       3.932  -3.053  12.086  1.00  0.00           H
HETATM   94  C15 PA A    1       2.307  -2.718  10.644  1.00  0.00           C
HETATM   95  H5R PA A    1       2.207  -1.635  10.870  1.00  0.00           H
HETATM   96  H5S PA A    1       1.331  -3.167  10.361  1.00  0.00           H
HETATM   97  C16 PA A    1       3.252  -2.693   9.370  1.00  0.00           C
HETATM   98  H6R PA A    1       3.267  -3.694   8.887  1.00  0.00           H
HETATM   99  H6S PA A    1       4.298  -2.522   9.701  1.00  0.00           H
HETATM  100  C17 PA A    1       2.825  -1.584   8.367  1.00  0.00           C
HETATM  101  H7R PA A    1       1.719  -1.524   8.448  1.00  0.00           H
HETATM  102  H7S PA A    1       3.026  -1.842   7.305  1.00  0.00           H
HETATM  103  C18 PA A    1       3.448  -0.195   8.584  1.00  0.00           C
HETATM  104  H8R PA A    1       4.515  -0.206   8.273  1.00  0.00           H
HETATM  105  H8S PA A    1       3.281   0.006   9.664  1.00  0.00           H
HETATM  106  C19 PA A    1       2.602   0.872   7.766  1.00  0.00           C
HETATM  107  H9R PA A    1       3.208   1.801   7.829  1.00  0.00           H
HETATM  108  H9S PA A    1       1.660   0.970   8.346  1.00  0.00           H
HETATM  109 C110 PA A    1       2.324   0.538   6.294  1.00  0.00           C
HETATM  110 H10R PA A    1       1.684  -0.363   6.181  1.00  0.00           H
HETATM  111 H10S PA A    1       3.327   0.335   5.862  1.00  0.00           H
HETATM  112 C111 PA A    1       1.574   1.703   5.582  1.00  0.00           C
HETATM  113 H11R PA A    1       2.095   2.638   5.881  1.00  0.00           H
HETATM  114 H11S PA A    1       0.536   1.768   5.974  1.00  0.00           H
HETATM  115 C112 PA A    1       1.760   1.558   4.043  1.00  0.00           C
HETATM  116 H12R PA A    1       1.372   0.556   3.763  1.00  0.00           H
HETATM  117 H12S PA A    1       2.840   1.661   3.803  1.00  0.00           H
HETATM  118 C113 PA A    1       1.107   2.634   3.247  1.00  0.00           C
HETATM  119 H13R PA A    1       1.529   3.661   3.284  1.00  0.00           H
HETATM  120 H13S PA A    1       0.045   2.826   3.512  1.00  0.00           H
HETATM  121 C114 PA A    1       1.116   2.216   1.720  1.00  0.00           C
HETATM  122 H14R PA A    1       0.722   1.186   1.586  1.00  0.00           H
HETATM  123 H14S PA A    1       2.214   2.175   1.557  1.00  0.00           H
HETATM  124 C115 PA A    1       0.499   3.186   0.747  1.00  0.00           C
HETATM  125 H15R PA A    1       0.903   4.167   1.077  1.00  0.00           H
HETATM  126 H15S PA A    1      -0.603   3.045   0.771  1.00  0.00           H
HETATM  127 C116 PA A    1       0.947   2.791  -0.682  1.00  0.00           C
HETATM  128 H16R PA A    1       0.285   3.374  -1.356  1.00  0.00           H
HETATM  129 H16S PA A    1       0.770   1.714  -0.891  1.00  0.00           H
HETATM  130 H16T PA A    1       2.017   3.015  -0.877  1.00  0.00           H
END}
	set f [open "DPPC_A.pdb" "w"]
	puts $f "$DPPC"
	close $f

	# POPC

	set POPC {HETATM    1  C12 PA  a   1      -9.708 -12.351  13.081  1.00  0.00           C
HETATM    2  H2R PA  a   1      -8.730 -11.952  12.735  1.00  0.00           H
HETATM    3  H2S PA  a   1     -10.364 -11.455  13.078  1.00  0.00           H
HETATM    4  C13 PA  a   1     -10.194 -13.260  11.958  1.00  0.00           C
HETATM    5  H3R PA  a   1     -11.266 -13.521  12.088  1.00  0.00           H
HETATM    6  H3S PA  a   1      -9.599 -14.195  12.043  1.00  0.00           H
HETATM    7  C14 PA  a   1     -10.141 -12.774  10.514  1.00  0.00           C
HETATM    8  H4R PA  a   1     -10.724 -13.465   9.869  1.00  0.00           H
HETATM    9  H4S PA  a   1      -9.102 -12.816  10.124  1.00  0.00           H
HETATM   10  C15 PA  a   1     -10.692 -11.318  10.244  1.00  0.00           C
HETATM   11  H5R PA  a   1     -10.066 -10.636  10.858  1.00  0.00           H
HETATM   12  H5S PA  a   1     -11.724 -11.216  10.645  1.00  0.00           H
HETATM   13  C16 PA  a   1     -10.691 -10.938   8.740  1.00  0.00           C
HETATM   14  H6R PA  a   1     -11.156  -9.929   8.736  1.00  0.00           H
HETATM   15  H6S PA  a   1     -11.264 -11.659   8.121  1.00  0.00           H
HETATM   16  C17 PA  a   1      -9.269 -10.611   8.263  1.00  0.00           C
HETATM   17  H7R PA  a   1      -8.829 -11.578   7.937  1.00  0.00           H
HETATM   18  H7S PA  a   1      -8.676 -10.182   9.099  1.00  0.00           H
HETATM   19  C18 PA  a   1      -9.347  -9.603   7.073  1.00  0.00           C
HETATM   20  H8R PA  a   1      -8.281  -9.330   6.922  1.00  0.00           H
HETATM   21  H8S PA  a   1      -9.998  -8.725   7.270  1.00  0.00           H
HETATM   22  C19 PA  a   1      -9.957 -10.274   5.740  1.00  0.00           C
HETATM   23  H9R PA  a   1     -10.989 -10.669   5.859  1.00  0.00           H
HETATM   24  H9S PA  a   1      -9.239 -11.054   5.410  1.00  0.00           H
HETATM   25 C110 PA  a   1     -10.110  -9.338   4.493  1.00  0.00           C
HETATM   26 H10R PA  a   1     -10.322  -8.293   4.803  1.00  0.00           H
HETATM   27 H10S PA  a   1     -10.899  -9.701   3.800  1.00  0.00           H
HETATM   28 C111 PA  a   1      -8.800  -9.349   3.587  1.00  0.00           C
HETATM   29 H11R PA  a   1      -7.988  -9.036   4.278  1.00  0.00           H
HETATM   30 H11S PA  a   1      -8.970  -8.540   2.844  1.00  0.00           H
HETATM   31 C112 PA  a   1      -8.524 -10.729   2.868  1.00  0.00           C
HETATM   32 H12R PA  a   1      -9.551 -10.969   2.519  1.00  0.00           H
HETATM   33 H12S PA  a   1      -8.264 -11.532   3.591  1.00  0.00           H
HETATM   34 C113 PA  a   1      -7.595 -10.485   1.727  1.00  0.00           C
HETATM   35 H13R PA  a   1      -8.259  -9.889   1.065  1.00  0.00           H
HETATM   36 H13S PA  a   1      -7.489 -11.377   1.074  1.00  0.00           H
HETATM   37 C114 PA  a   1      -6.239  -9.845   2.069  1.00  0.00           C
HETATM   38 H14R PA  a   1      -5.640 -10.306   2.883  1.00  0.00           H
HETATM   39 H14S PA  a   1      -6.463  -8.779   2.286  1.00  0.00           H
HETATM   40 C115 PA  a   1      -5.445  -9.989   0.758  1.00  0.00           C
HETATM   41 H15R PA  a   1      -5.230 -11.078   0.716  1.00  0.00           H
HETATM   42 H15S PA  a   1      -4.439  -9.521   0.813  1.00  0.00           H
HETATM   43 C116 PA  a   1      -6.247  -9.412  -0.408  1.00  0.00           C
HETATM   44 H16R PA  a   1      -6.583  -8.394  -0.115  1.00  0.00           H
HETATM   45 H16S PA  a   1      -7.079 -10.086  -0.706  1.00  0.00           H
HETATM   46 H16T PA  a   1      -5.688  -9.368  -1.366  1.00  0.00           H
HETATM   47  N31 PC  a   1      -5.289 -14.415  19.437  1.00  0.00           N
HETATM   48  C32 PC  a   1      -4.530 -13.102  19.824  1.00  0.00           C
HETATM   49  H2A PC  a   1      -4.105 -12.733  18.901  1.00  0.00           H
HETATM   50  H2B PC  a   1      -3.771 -13.390  20.536  1.00  0.00           H
HETATM   51  C33 PC  a   1      -4.488 -15.224  18.488  1.00  0.00           C
HETATM   52  H3A PC  a   1      -3.798 -14.655  17.884  1.00  0.00           H
HETATM   53  H3B PC  a   1      -3.877 -15.934  19.025  1.00  0.00           H
HETATM   54  H3C PC  a   1      -5.095 -15.862  17.861  1.00  0.00           H
HETATM   55  C34 PC  a   1      -5.537 -15.175  20.698  1.00  0.00           C
HETATM   56  H4A PC  a   1      -5.968 -14.455  21.377  1.00  0.00           H
HETATM   57  H4B PC  a   1      -6.168 -16.030  20.508  1.00  0.00           H
HETATM   58  H4C PC  a   1      -4.678 -15.647  21.152  1.00  0.00           H
HETATM   59  C35 PC  a   1      -6.507 -14.076  18.667  1.00  0.00           C
HETATM   60  H5A PC  a   1      -7.138 -14.909  18.392  1.00  0.00           H
HETATM   61  H5B PC  a   1      -6.219 -13.490  17.806  1.00  0.00           H
HETATM   62  H5C PC  a   1      -7.058 -13.423  19.327  1.00  0.00           H
HETATM   63  C31 PC  a   1      -5.330 -11.912  20.412  1.00  0.00           C
HETATM   64  H1A PC  a   1      -5.673 -12.640  19.646  1.00  0.00           H
HETATM   65  H1B PC  a   1      -4.269 -11.605  20.292  1.00  0.00           H
HETATM   66  P31 PC  a   1      -6.117  -9.590  19.434  1.00  0.00           P
HETATM   67  O33 PC  a   1      -4.859  -8.888  19.504  1.00  0.00           O
HETATM   68  O34 PC  a   1      -7.255  -8.747  19.647  1.00  0.00           O
HETATM   69  O32 PC  a   1      -6.201 -10.733  20.498  1.00  0.00           O
HETATM   70  O31 PC  a   1      -6.310 -10.402  18.099  1.00  0.00           O
HETATM   71  C3  PC  a   1      -7.553 -10.767  17.495  1.00  0.00           C
HETATM   72  HA  PC  a   1      -7.686 -10.145  16.584  1.00  0.00           H
HETATM   73  HB  PC  a   1      -8.440 -10.752  18.164  1.00  0.00           H
HETATM   74  C2  PC  a   1      -7.434 -12.244  17.046  1.00  0.00           C
HETATM   75  HX  PC  a   1      -7.159 -12.670  18.034  1.00  0.00           H
HETATM   76  O21 PC  a   1      -6.374 -12.372  16.159  1.00  0.00           O
HETATM   77  C21 PC  a   1      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   78  O22 PC  a   1      -5.946 -14.499  16.694  1.00  0.00           O
HETATM   79  C1  PC  a   1      -8.669 -12.908  16.452  1.00  0.00           C
HETATM   80  HR  PC  a   1      -9.572 -12.745  17.077  1.00  0.00           H
HETATM   81  HS  PC  a   1      -8.470 -14.000  16.396  1.00  0.00           H
HETATM   82  O11 PC  a   1      -8.960 -12.322  15.217  1.00  0.00           O
HETATM   83  C11 PC  a   1      -9.722 -13.074  14.448  1.00  0.00           C
HETATM   84  O12 PC  a   1     -10.331 -14.103  14.689  1.00  0.00           O
HETATM   85  C12 OL  a   1      -5.038 -12.692  14.886  1.00  0.00           C
HETATM   86  H2R OL  a   1      -3.962 -12.929  15.034  1.00  0.00           H
HETATM   87  H2S OL  a   1      -5.223 -11.604  15.010  1.00  0.00           H
HETATM   88  C13 OL  a   1      -5.565 -13.291  13.565  1.00  0.00           C
HETATM   89  H3R OL  a   1      -6.674 -13.323  13.501  1.00  0.00           H
HETATM   90  H3S OL  a   1      -5.106 -14.290  13.402  1.00  0.00           H
HETATM   91  C14 OL  a   1      -5.066 -12.562  12.388  1.00  0.00           C
HETATM   92  H4R OL  a   1      -5.609 -12.965  11.507  1.00  0.00           H
HETATM   93  H4S OL  a   1      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   94  C15 OL  a   1      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   95  H5R OL  a   1      -4.722 -10.479  12.952  1.00  0.00           H
HETATM   96  H5S OL  a   1      -6.378 -10.838  12.716  1.00  0.00           H
HETATM   97  C16 OL  a   1      -5.328 -10.465  10.922  1.00  0.00           C
HETATM   98  H6R OL  a   1      -5.283  -9.355  10.939  1.00  0.00           H
HETATM   99  H6S OL  a   1      -6.083 -10.781  10.172  1.00  0.00           H
HETATM  100  C17 OL  a   1      -4.072 -10.915  10.088  1.00  0.00           C
HETATM  101  H7R OL  a   1      -3.803 -11.985  10.224  1.00  0.00           H
HETATM  102  H7S OL  a   1      -3.228 -10.288  10.447  1.00  0.00           H
HETATM  103  C18 OL  a   1      -4.096 -10.707   8.561  1.00  0.00           C
HETATM  104  H8R OL  a   1      -3.023 -10.549   8.317  1.00  0.00           H
HETATM  105  H8S OL  a   1      -4.539  -9.771   8.158  1.00  0.00           H
HETATM  106  C19 OL  a   1      -4.543 -11.872   7.693  1.00  0.00           C
HETATM  107  H9R OL  a   1      -4.056 -12.833   7.917  1.00  0.00           H
HETATM  108 C110 OL  a   1      -5.517 -11.984   6.700  1.00  0.00           C
HETATM  109 H10R OL  a   1      -5.665 -12.996   6.296  1.00  0.00           H
HETATM  110 C111 OL  a   1      -6.391 -10.916   6.211  1.00  0.00           C
HETATM  111 H11R OL  a   1      -6.147  -9.972   6.743  1.00  0.00           H
HETATM  112 H11S OL  a   1      -7.457 -11.176   6.387  1.00  0.00           H
HETATM  113 C112 OL  a   1      -6.043 -10.766   4.699  1.00  0.00           C
HETATM  114 H12R OL  a   1      -4.989 -10.419   4.634  1.00  0.00           H
HETATM  115 H12S OL  a   1      -6.774 -10.059   4.253  1.00  0.00           H
HETATM  116 C113 OL  a   1      -6.272 -11.989   3.860  1.00  0.00           C
HETATM  117 H13R OL  a   1      -7.325 -12.174   4.161  1.00  0.00           H
HETATM  118 H13S OL  a   1      -5.670 -12.857   4.208  1.00  0.00           H
HETATM  119 C114 OL  a   1      -5.995 -11.819   2.357  1.00  0.00           C
HETATM  120 H14R OL  a   1      -4.903 -11.898   2.168  1.00  0.00           H
HETATM  121 H14S OL  a   1      -6.282 -10.759   2.192  1.00  0.00           H
HETATM  122 C115 OL  a   1      -6.855 -12.761   1.560  1.00  0.00           C
HETATM  123 H15R OL  a   1      -7.917 -12.453   1.672  1.00  0.00           H
HETATM  124 H15S OL  a   1      -6.820 -13.811   1.919  1.00  0.00           H
HETATM  125 C116 OL  a   1      -6.528 -12.913   0.035  1.00  0.00           C
HETATM  126 H16R OL  a   1      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM  127 H16S OL  a   1      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM  128 C117 OL  a   1      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM  129 H17R OL  a   1      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM  130 H17S OL  a   1      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM  131 C118 OL  a   1      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM  132 H18R OL  a   1      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM  133 H18S OL  a   1      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM  134 H18T OL  a   1      -5.769 -10.553  -2.572  1.00  0.00           H
END}

	set f [open "POPC_A.pdb" "w"]
	puts $f "$POPC"
	close $f

	# POPE

	set POPE {HETATM    1  C12 OL  a   1      -5.038 -12.692  14.886  1.00  0.00           C
HETATM    2  H2R OL  a   1      -3.962 -12.929  15.034  1.00  0.00           H
HETATM    3  H2S OL  a   1      -5.223 -11.604  15.010  1.00  0.00           H
HETATM    4  C13 OL  a   1      -5.565 -13.291  13.565  1.00  0.00           C
HETATM    5  H3R OL  a   1      -6.674 -13.323  13.501  1.00  0.00           H
HETATM    6  H3S OL  a   1      -5.106 -14.290  13.402  1.00  0.00           H
HETATM    7  C14 OL  a   1      -5.066 -12.562  12.388  1.00  0.00           C
HETATM    8  H4R OL  a   1      -5.609 -12.965  11.507  1.00  0.00           H
HETATM    9  H4S OL  a   1      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   10  C15 OL  a   1      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   11  H5R OL  a   1      -4.722 -10.479  12.952  1.00  0.00           H
HETATM   12  H5S OL  a   1      -6.378 -10.838  12.716  1.00  0.00           H
HETATM   13  C16 OL  a   1      -5.328 -10.465  10.922  1.00  0.00           C
HETATM   14  H6R OL  a   1      -5.283  -9.355  10.939  1.00  0.00           H
HETATM   15  H6S OL  a   1      -6.083 -10.781  10.172  1.00  0.00           H
HETATM   16  C17 OL  a   1      -4.072 -10.915  10.088  1.00  0.00           C
HETATM   17  H7R OL  a   1      -3.803 -11.985  10.224  1.00  0.00           H
HETATM   18  H7S OL  a   1      -3.228 -10.288  10.447  1.00  0.00           H
HETATM   19  C18 OL  a   1      -4.096 -10.707   8.561  1.00  0.00           C
HETATM   20  H8R OL  a   1      -3.023 -10.549   8.317  1.00  0.00           H
HETATM   21  H8S OL  a   1      -4.539  -9.771   8.158  1.00  0.00           H
HETATM   22  C19 OL  a   1      -4.543 -11.872   7.693  1.00  0.00           C
HETATM   23  H9R OL  a   1      -4.056 -12.833   7.917  1.00  0.00           H
HETATM   24 C110 OL  a   1      -5.517 -11.984   6.700  1.00  0.00           C
HETATM   25 H10R OL  a   1      -5.665 -12.996   6.296  1.00  0.00           H
HETATM   26 C111 OL  a   1      -6.391 -10.916   6.211  1.00  0.00           C
HETATM   27 H11R OL  a   1      -6.147  -9.972   6.743  1.00  0.00           H
HETATM   28 H11S OL  a   1      -7.457 -11.176   6.387  1.00  0.00           H
HETATM   29 C112 OL  a   1      -6.043 -10.766   4.699  1.00  0.00           C
HETATM   30 H12R OL  a   1      -4.989 -10.419   4.634  1.00  0.00           H
HETATM   31 H12S OL  a   1      -6.774 -10.059   4.253  1.00  0.00           H
HETATM   32 C113 OL  a   1      -6.272 -11.989   3.860  1.00  0.00           C
HETATM   33 H13R OL  a   1      -7.325 -12.174   4.161  1.00  0.00           H
HETATM   34 H13S OL  a   1      -5.670 -12.857   4.208  1.00  0.00           H
HETATM   35 C114 OL  a   1      -5.995 -11.819   2.357  1.00  0.00           C
HETATM   36 H14R OL  a   1      -4.903 -11.898   2.168  1.00  0.00           H
HETATM   37 H14S OL  a   1      -6.282 -10.759   2.192  1.00  0.00           H
HETATM   38 C115 OL  a   1      -6.855 -12.761   1.560  1.00  0.00           C
HETATM   39 H15R OL  a   1      -7.917 -12.453   1.672  1.00  0.00           H
HETATM   40 H15S OL  a   1      -6.820 -13.811   1.919  1.00  0.00           H
HETATM   41 C116 OL  a   1      -6.528 -12.913   0.035  1.00  0.00           C
HETATM   42 H16R OL  a   1      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM   43 H16S OL  a   1      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM   44 C117 OL  a   1      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM   45 H17R OL  a   1      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM   46 H17S OL  a   1      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM   47 C118 OL  a   1      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM   48 H18R OL  a   1      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM   49 H18S OL  a   1      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM   50 H18T OL  a   1      -5.769 -10.553  -2.572  1.00  0.00           H
HETATM   51  N31 PE  a   1      -3.527 -12.627  18.830  1.00  0.00           N
HETATM   52 HN1A PE  a   1      -2.915 -11.859  19.172  1.00  0.00           H
HETATM   53 HN1B PE  a   1      -2.932 -13.427  18.533  1.00  0.00           H
HETATM   54 HN1C PE  a   1      -4.156 -12.281  18.078  1.00  0.00           H
HETATM   55  C32 PE  a   1      -4.558 -13.061  19.844  1.00  0.00           C
HETATM   56  H2A PE  a   1      -3.946 -13.622  20.583  1.00  0.00           H
HETATM   57  H2B PE  a   1      -5.206 -13.787  19.308  1.00  0.00           H
HETATM   58  C31 PE  a   1      -5.330 -11.912  20.412  1.00  0.00           C
HETATM   59  H1A PE  a   1      -5.673 -12.640  19.646  1.00  0.00           H
HETATM   60  H1B PE  a   1      -4.269 -11.605  20.292  1.00  0.00           H
HETATM   61  P31 PE  a   1      -6.117  -9.590  19.434  1.00  0.00           P
HETATM   62  O33 PE  a   1      -4.859  -8.888  19.504  1.00  0.00           O
HETATM   63  O34 PE  a   1      -7.255  -8.747  19.647  1.00  0.00           O
HETATM   64  O31 PE  a   1      -6.310 -10.402  18.099  1.00  0.00           O
HETATM   65  O32 PE  a   1      -6.201 -10.733  20.498  1.00  0.00           O
HETATM   66  C3  PE  a   1      -7.553 -10.767  17.495  1.00  0.00           C
HETATM   67  HA  PE  a   1      -7.686 -10.145  16.584  1.00  0.00           H
HETATM   68  HB  PE  a   1      -8.440 -10.752  18.164  1.00  0.00           H
HETATM   69  C2  PE  a   1      -7.434 -12.244  17.046  1.00  0.00           C
HETATM   70  HX  PE  a   1      -7.159 -12.670  18.034  1.00  0.00           H
HETATM   71  O21 PE  a   1      -6.374 -12.372  16.159  1.00  0.00           O
HETATM   72  C21 PE  a   1      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   73  O22 PE  a   1      -5.946 -14.499  16.694  1.00  0.00           O
HETATM   74  C1  PE  a   1      -8.669 -12.908  16.452  1.00  0.00           C
HETATM   75  HR  PE  a   1      -9.572 -12.745  17.077  1.00  0.00           H
HETATM   76  HS  PE  a   1      -8.470 -14.000  16.396  1.00  0.00           H
HETATM   77  O11 PE  a   1      -8.960 -12.322  15.217  1.00  0.00           O
HETATM   78  C11 PE  a   1      -9.722 -13.074  14.448  1.00  0.00           C
HETATM   79  O12 PE  a   1     -10.331 -14.103  14.689  1.00  0.00           O
HETATM   80  C12 PA  a   1      -9.708 -12.351  13.081  1.00  0.00           C
HETATM   81  H2R PA  a   1      -8.730 -11.952  12.735  1.00  0.00           H
HETATM   82  H2S PA  a   1     -10.364 -11.455  13.078  1.00  0.00           H
HETATM   83  C13 PA  a   1     -10.194 -13.260  11.958  1.00  0.00           C
HETATM   84  H3R PA  a   1     -11.266 -13.521  12.088  1.00  0.00           H
HETATM   85  H3S PA  a   1      -9.599 -14.195  12.043  1.00  0.00           H
HETATM   86  C14 PA  a   1     -10.141 -12.774  10.514  1.00  0.00           C
HETATM   87  H4R PA  a   1     -10.724 -13.465   9.869  1.00  0.00           H
HETATM   88  H4S PA  a   1      -9.102 -12.816  10.124  1.00  0.00           H
HETATM   89  C15 PA  a   1     -10.692 -11.318  10.244  1.00  0.00           C
HETATM   90  H5R PA  a   1     -10.066 -10.636  10.858  1.00  0.00           H
HETATM   91  H5S PA  a   1     -11.724 -11.216  10.645  1.00  0.00           H
HETATM   92  C16 PA  a   1     -10.691 -10.938   8.740  1.00  0.00           C
HETATM   93  H6R PA  a   1     -11.156  -9.929   8.736  1.00  0.00           H
HETATM   94  H6S PA  a   1     -11.264 -11.659   8.121  1.00  0.00           H
HETATM   95  C17 PA  a   1      -9.269 -10.611   8.263  1.00  0.00           C
HETATM   96  H7R PA  a   1      -8.829 -11.578   7.937  1.00  0.00           H
HETATM   97  H7S PA  a   1      -8.676 -10.182   9.099  1.00  0.00           H
HETATM   98  C18 PA  a   1      -9.347  -9.603   7.073  1.00  0.00           C
HETATM   99  H8R PA  a   1      -8.281  -9.330   6.922  1.00  0.00           H
HETATM  100  H8S PA  a   1      -9.998  -8.725   7.270  1.00  0.00           H
HETATM  101  C19 PA  a   1      -9.957 -10.274   5.740  1.00  0.00           C
HETATM  102  H9R PA  a   1     -10.989 -10.669   5.859  1.00  0.00           H
HETATM  103  H9S PA  a   1      -9.239 -11.054   5.410  1.00  0.00           H
HETATM  104 C110 PA  a   1     -10.110  -9.338   4.493  1.00  0.00           C
HETATM  105 H10R PA  a   1     -10.322  -8.293   4.803  1.00  0.00           H
HETATM  106 H10S PA  a   1     -10.899  -9.701   3.800  1.00  0.00           H
HETATM  107 C111 PA  a   1      -8.800  -9.349   3.587  1.00  0.00           C
HETATM  108 H11R PA  a   1      -7.988  -9.036   4.278  1.00  0.00           H
HETATM  109 H11S PA  a   1      -8.970  -8.540   2.844  1.00  0.00           H
HETATM  110 C112 PA  a   1      -8.524 -10.729   2.868  1.00  0.00           C
HETATM  111 H12R PA  a   1      -9.551 -10.969   2.519  1.00  0.00           H
HETATM  112 H12S PA  a   1      -8.264 -11.532   3.591  1.00  0.00           H
HETATM  113 C113 PA  a   1      -7.595 -10.485   1.727  1.00  0.00           C
HETATM  114 H13R PA  a   1      -8.259  -9.889   1.065  1.00  0.00           H
HETATM  115 H13S PA  a   1      -7.489 -11.377   1.074  1.00  0.00           H
HETATM  116 C114 PA  a   1      -6.239  -9.845   2.069  1.00  0.00           C
HETATM  117 H14R PA  a   1      -5.640 -10.306   2.883  1.00  0.00           H
HETATM  118 H14S PA  a   1      -6.463  -8.779   2.286  1.00  0.00           H
HETATM  119 C115 PA  a   1      -5.445  -9.989   0.758  1.00  0.00           C
HETATM  120 H15R PA  a   1      -5.230 -11.078   0.716  1.00  0.00           H
HETATM  121 H15S PA  a   1      -4.439  -9.521   0.813  1.00  0.00           H
HETATM  122 C116 PA  a   1      -6.247  -9.412  -0.408  1.00  0.00           C
HETATM  123 H16R PA  a   1      -6.583  -8.394  -0.115  1.00  0.00           H
HETATM  124 H16S PA  a   1      -7.079 -10.086  -0.706  1.00  0.00           H
HETATM  125 H16T PA  a   1      -5.688  -9.368  -1.366  1.00  0.00           H
END}

	set f [open "POPE_A.pdb" "w"]
	puts $f "$POPE"
	close $f

	# POPG

	set POPG {HETATM    1  C12 PA  a   2      -9.708 -12.351  13.081  1.00  0.00           C
HETATM    2  H2R PA  a   2      -8.730 -11.952  12.735  1.00  0.00           H
HETATM    3  H2S PA  a   2     -10.364 -11.455  13.078  1.00  0.00           H
HETATM    4  C13 PA  a   2     -10.194 -13.260  11.958  1.00  0.00           C
HETATM    5  H3R PA  a   2     -11.266 -13.521  12.088  1.00  0.00           H
HETATM    6  H3S PA  a   2      -9.599 -14.195  12.043  1.00  0.00           H
HETATM    7  C14 PA  a   2     -10.141 -12.774  10.514  1.00  0.00           C
HETATM    8  H4R PA  a   2     -10.724 -13.465   9.869  1.00  0.00           H
HETATM    9  H4S PA  a   2      -9.102 -12.816  10.124  1.00  0.00           H
HETATM   10  C15 PA  a   2     -10.692 -11.318  10.244  1.00  0.00           C
HETATM   11  H5R PA  a   2     -10.066 -10.636  10.858  1.00  0.00           H
HETATM   12  H5S PA  a   2     -11.724 -11.216  10.645  1.00  0.00           H
HETATM   13  C16 PA  a   2     -10.691 -10.938   8.740  1.00  0.00           C
HETATM   14  H6R PA  a   2     -11.156  -9.929   8.736  1.00  0.00           H
HETATM   15  H6S PA  a   2     -11.264 -11.659   8.121  1.00  0.00           H
HETATM   16  C17 PA  a   2      -9.269 -10.611   8.263  1.00  0.00           C
HETATM   17  H7R PA  a   2      -8.829 -11.578   7.937  1.00  0.00           H
HETATM   18  H7S PA  a   2      -8.676 -10.182   9.099  1.00  0.00           H
HETATM   19  C18 PA  a   2      -9.347  -9.603   7.073  1.00  0.00           C
HETATM   20  H8R PA  a   2      -8.281  -9.330   6.922  1.00  0.00           H
HETATM   21  H8S PA  a   2      -9.998  -8.725   7.270  1.00  0.00           H
HETATM   22  C19 PA  a   2      -9.957 -10.274   5.740  1.00  0.00           C
HETATM   23  H9R PA  a   2     -10.989 -10.669   5.859  1.00  0.00           H
HETATM   24  H9S PA  a   2      -9.239 -11.054   5.410  1.00  0.00           H
HETATM   25 C110 PA  a   2     -10.110  -9.338   4.493  1.00  0.00           C
HETATM   26 H10R PA  a   2     -10.322  -8.293   4.803  1.00  0.00           H
HETATM   27 H10S PA  a   2     -10.899  -9.701   3.800  1.00  0.00           H
HETATM   28 C111 PA  a   2      -8.800  -9.349   3.587  1.00  0.00           C
HETATM   29 H11R PA  a   2      -7.988  -9.036   4.278  1.00  0.00           H
HETATM   30 H11S PA  a   2      -8.970  -8.540   2.844  1.00  0.00           H
HETATM   31 C112 PA  a   2      -8.524 -10.729   2.868  1.00  0.00           C
HETATM   32 H12R PA  a   2      -9.551 -10.969   2.519  1.00  0.00           H
HETATM   33 H12S PA  a   2      -8.264 -11.532   3.591  1.00  0.00           H
HETATM   34 C113 PA  a   2      -7.595 -10.485   1.727  1.00  0.00           C
HETATM   35 H13R PA  a   2      -8.259  -9.889   1.065  1.00  0.00           H
HETATM   36 H13S PA  a   2      -7.489 -11.377   1.074  1.00  0.00           H
HETATM   37 C114 PA  a   2      -6.239  -9.845   2.069  1.00  0.00           C
HETATM   38 H14R PA  a   2      -5.640 -10.306   2.883  1.00  0.00           H
HETATM   39 H14S PA  a   2      -6.463  -8.779   2.286  1.00  0.00           H
HETATM   40 C115 PA  a   2      -5.445  -9.989   0.758  1.00  0.00           C
HETATM   41 H15R PA  a   2      -5.230 -11.078   0.716  1.00  0.00           H
HETATM   42 H15S PA  a   2      -4.439  -9.521   0.813  1.00  0.00           H
HETATM   43 C116 PA  a   2      -6.247  -9.412  -0.408  1.00  0.00           C
HETATM   44 H16R PA  a   2      -6.583  -8.394  -0.115  1.00  0.00           H
HETATM   45 H16S PA  a   2      -7.079 -10.086  -0.706  1.00  0.00           H
HETATM   46 H16T PA  a   2      -5.688  -9.368  -1.366  1.00  0.00           H
HETATM   47  C33 PGR a   2      -3.821 -12.598  18.543  1.00  0.00           C
HETATM   48  H3A PGR a   2      -3.126 -11.740  18.661  1.00  0.00           H
HETATM   49  H3B PGR a   2      -3.184 -13.397  18.105  1.00  0.00           H
HETATM   50  O36 PGR a   2      -4.780 -12.216  17.581  1.00  0.00           O
HETATM   51 HO6A PGR a   2      -4.795 -12.778  16.803  1.00  0.00           H
HETATM   52  C32 PGR a   2      -4.528 -13.106  19.822  1.00  0.00           C
HETATM   53  H2A PGR a   2      -5.082 -14.048  19.620  1.00  0.00           H
HETATM   54  O35 PGR a   2      -3.592 -13.540  20.829  1.00  0.00           O
HETATM   55 HO5A PGR a   2      -4.117 -13.381  21.616  1.00  0.00           H
HETATM   56  C31 PGR a   2      -5.330 -11.912  20.412  1.00  0.00           C
HETATM   57  H1A PGR a   2      -5.673 -12.640  19.646  1.00  0.00           H
HETATM   58  H1B PGR a   2      -4.269 -11.605  20.292  1.00  0.00           H
HETATM   59  P31 PGR a   2      -6.117  -9.590  19.434  1.00  0.00           P
HETATM   60  O33 PGR a   2      -4.859  -8.888  19.504  1.00  0.00           O
HETATM   61  O34 PGR a   2      -7.255  -8.747  19.647  1.00  0.00           O
HETATM   62  O32 PGR a   2      -6.201 -10.733  20.498  1.00  0.00           O
HETATM   63  O31 PGR a   2      -6.310 -10.402  18.099  1.00  0.00           O
HETATM   64  C3  PGR a   2      -7.553 -10.767  17.495  1.00  0.00           C
HETATM   65  HA  PGR a   2      -7.686 -10.145  16.584  1.00  0.00           H
HETATM   66  HB  PGR a   2      -8.440 -10.752  18.164  1.00  0.00           H
HETATM   67  C2  PGR a   2      -7.434 -12.244  17.046  1.00  0.00           C
HETATM   68  HX  PGR a   2      -7.159 -12.670  18.034  1.00  0.00           H
HETATM   69  O21 PGR a   2      -6.374 -12.372  16.159  1.00  0.00           O
HETATM   70  C21 PGR a   2      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   71  O22 PGR a   2      -5.946 -14.499  16.694  1.00  0.00           O
HETATM   72  C1  PGR a   2      -8.669 -12.908  16.452  1.00  0.00           C
HETATM   73  HR  PGR a   2      -9.572 -12.745  17.077  1.00  0.00           H
HETATM   74  HS  PGR a   2      -8.470 -14.000  16.396  1.00  0.00           H
HETATM   75  O11 PGR a   2      -8.960 -12.322  15.217  1.00  0.00           O
HETATM   76  C11 PGR a   2      -9.722 -13.074  14.448  1.00  0.00           C
HETATM   77  O12 PGR a   2     -10.331 -14.103  14.689  1.00  0.00           O
HETATM   78  C12 OL  a   2      -5.038 -12.692  14.886  1.00  0.00           C
HETATM   79  H2R OL  a   2      -3.962 -12.929  15.034  1.00  0.00           H
HETATM   80  H2S OL  a   2      -5.223 -11.604  15.010  1.00  0.00           H
HETATM   81  C13 OL  a   2      -5.565 -13.291  13.565  1.00  0.00           C
HETATM   82  H3R OL  a   2      -6.674 -13.323  13.501  1.00  0.00           H
HETATM   83  H3S OL  a   2      -5.106 -14.290  13.402  1.00  0.00           H
HETATM   84  C14 OL  a   2      -5.066 -12.562  12.388  1.00  0.00           C
HETATM   85  H4R OL  a   2      -5.609 -12.965  11.507  1.00  0.00           H
HETATM   86  H4S OL  a   2      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   87  C15 OL  a   2      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   88  H5R OL  a   2      -4.722 -10.479  12.952  1.00  0.00           H
HETATM   89  H5S OL  a   2      -6.378 -10.838  12.716  1.00  0.00           H
HETATM   90  C16 OL  a   2      -5.328 -10.465  10.922  1.00  0.00           C
HETATM   91  H6R OL  a   2      -5.283  -9.355  10.939  1.00  0.00           H
HETATM   92  H6S OL  a   2      -6.083 -10.781  10.172  1.00  0.00           H
HETATM   93  C17 OL  a   2      -4.072 -10.915  10.088  1.00  0.00           C
HETATM   94  H7R OL  a   2      -3.803 -11.985  10.224  1.00  0.00           H
HETATM   95  H7S OL  a   2      -3.228 -10.288  10.447  1.00  0.00           H
HETATM   96  C18 OL  a   2      -4.096 -10.707   8.561  1.00  0.00           C
HETATM   97  H8R OL  a   2      -3.023 -10.549   8.317  1.00  0.00           H
HETATM   98  H8S OL  a   2      -4.539  -9.771   8.158  1.00  0.00           H
HETATM   99  C19 OL  a   2      -4.543 -11.872   7.693  1.00  0.00           C
HETATM  100  H9R OL  a   2      -4.056 -12.833   7.917  1.00  0.00           H
HETATM  101 C110 OL  a   2      -5.517 -11.984   6.700  1.00  0.00           C
HETATM  102 H10R OL  a   2      -5.665 -12.996   6.296  1.00  0.00           H
HETATM  103 C111 OL  a   2      -6.391 -10.916   6.211  1.00  0.00           C
HETATM  104 H11R OL  a   2      -6.147  -9.972   6.743  1.00  0.00           H
HETATM  105 H11S OL  a   2      -7.457 -11.176   6.387  1.00  0.00           H
HETATM  106 C112 OL  a   2      -6.043 -10.766   4.699  1.00  0.00           C
HETATM  107 H12R OL  a   2      -4.989 -10.419   4.634  1.00  0.00           H
HETATM  108 H12S OL  a   2      -6.774 -10.059   4.253  1.00  0.00           H
HETATM  109 C113 OL  a   2      -6.272 -11.989   3.860  1.00  0.00           C
HETATM  110 H13R OL  a   2      -7.325 -12.174   4.161  1.00  0.00           H
HETATM  111 H13S OL  a   2      -5.670 -12.857   4.208  1.00  0.00           H
HETATM  112 C114 OL  a   2      -5.995 -11.819   2.357  1.00  0.00           C
HETATM  113 H14R OL  a   2      -4.903 -11.898   2.168  1.00  0.00           H
HETATM  114 H14S OL  a   2      -6.282 -10.759   2.192  1.00  0.00           H
HETATM  115 C115 OL  a   2      -6.855 -12.761   1.560  1.00  0.00           C
HETATM  116 H15R OL  a   2      -7.917 -12.453   1.672  1.00  0.00           H
HETATM  117 H15S OL  a   2      -6.820 -13.811   1.919  1.00  0.00           H
HETATM  118 C116 OL  a   2      -6.528 -12.913   0.035  1.00  0.00           C
HETATM  119 H16R OL  a   2      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM  120 H16S OL  a   2      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM  121 C117 OL  a   2      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM  122 H17R OL  a   2      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM  123 H17S OL  a   2      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM  124 C118 OL  a   2      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM  125 H18R OL  a   2      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM  126 H18S OL  a   2      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM  127 H18T OL  a   2      -5.769 -10.553  -2.572  1.00  0.00           H
END}

	set f [open "POPG_A.pdb" "w"]
	puts $f "$POPG"
	close $f

	# POPS

	set POPS {HETATM    1  C12 PA  a   1      -9.708 -12.351  13.081  1.00  0.00           C
HETATM    2  H2R PA  a   1      -8.730 -11.952  12.735  1.00  0.00           H
HETATM    3  H2S PA  a   1     -10.364 -11.455  13.078  1.00  0.00           H
HETATM    4  C13 PA  a   1     -10.194 -13.260  11.958  1.00  0.00           C
HETATM    5  H3R PA  a   1     -11.266 -13.521  12.088  1.00  0.00           H
HETATM    6  H3S PA  a   1      -9.599 -14.195  12.043  1.00  0.00           H
HETATM    7  C14 PA  a   1     -10.141 -12.774  10.514  1.00  0.00           C
HETATM    8  H4R PA  a   1     -10.724 -13.465   9.869  1.00  0.00           H
HETATM    9  H4S PA  a   1      -9.102 -12.816  10.124  1.00  0.00           H
HETATM   10  C15 PA  a   1     -10.692 -11.318  10.244  1.00  0.00           C
HETATM   11  H5R PA  a   1     -10.066 -10.636  10.858  1.00  0.00           H
HETATM   12  H5S PA  a   1     -11.724 -11.216  10.645  1.00  0.00           H
HETATM   13  C16 PA  a   1     -10.691 -10.938   8.740  1.00  0.00           C
HETATM   14  H6R PA  a   1     -11.156  -9.929   8.736  1.00  0.00           H
HETATM   15  H6S PA  a   1     -11.264 -11.659   8.121  1.00  0.00           H
HETATM   16  C17 PA  a   1      -9.269 -10.611   8.263  1.00  0.00           C
HETATM   17  H7R PA  a   1      -8.829 -11.578   7.937  1.00  0.00           H
HETATM   18  H7S PA  a   1      -8.676 -10.182   9.099  1.00  0.00           H
HETATM   19  C18 PA  a   1      -9.347  -9.603   7.073  1.00  0.00           C
HETATM   20  H8R PA  a   1      -8.281  -9.330   6.922  1.00  0.00           H
HETATM   21  H8S PA  a   1      -9.998  -8.725   7.270  1.00  0.00           H
HETATM   22  C19 PA  a   1      -9.957 -10.274   5.740  1.00  0.00           C
HETATM   23  H9R PA  a   1     -10.989 -10.669   5.859  1.00  0.00           H
HETATM   24  H9S PA  a   1      -9.239 -11.054   5.410  1.00  0.00           H
HETATM   25 C110 PA  a   1     -10.110  -9.338   4.493  1.00  0.00           C
HETATM   26 H10R PA  a   1     -10.322  -8.293   4.803  1.00  0.00           H
HETATM   27 H10S PA  a   1     -10.899  -9.701   3.800  1.00  0.00           H
HETATM   28 C111 PA  a   1      -8.800  -9.349   3.587  1.00  0.00           C
HETATM   29 H11R PA  a   1      -7.988  -9.036   4.278  1.00  0.00           H
HETATM   30 H11S PA  a   1      -8.970  -8.540   2.844  1.00  0.00           H
HETATM   31 C112 PA  a   1      -8.524 -10.729   2.868  1.00  0.00           C
HETATM   32 H12R PA  a   1      -9.551 -10.969   2.519  1.00  0.00           H
HETATM   33 H12S PA  a   1      -8.264 -11.532   3.591  1.00  0.00           H
HETATM   34 C113 PA  a   1      -7.595 -10.485   1.727  1.00  0.00           C
HETATM   35 H13R PA  a   1      -8.259  -9.889   1.065  1.00  0.00           H
HETATM   36 H13S PA  a   1      -7.489 -11.377   1.074  1.00  0.00           H
HETATM   37 C114 PA  a   1      -6.239  -9.845   2.069  1.00  0.00           C
HETATM   38 H14R PA  a   1      -5.640 -10.306   2.883  1.00  0.00           H
HETATM   39 H14S PA  a   1      -6.463  -8.779   2.286  1.00  0.00           H
HETATM   40 C115 PA  a   1      -5.445  -9.989   0.758  1.00  0.00           C
HETATM   41 H15R PA  a   1      -5.230 -11.078   0.716  1.00  0.00           H
HETATM   42 H15S PA  a   1      -4.439  -9.521   0.813  1.00  0.00           H
HETATM   43 C116 PA  a   1      -6.247  -9.412  -0.408  1.00  0.00           C
HETATM   44 H16R PA  a   1      -6.583  -8.394  -0.115  1.00  0.00           H
HETATM   45 H16S PA  a   1      -7.079 -10.086  -0.706  1.00  0.00           H
HETATM   46 H16T PA  a   1      -5.688  -9.368  -1.366  1.00  0.00           H
HETATM   47  N31 PS  a   1      -5.479 -14.200  19.423  1.00  0.00           N
HETATM   48 HN1A PS  a   1      -5.101 -15.168  19.396  1.00  0.00           H
HETATM   49 HN1B PS  a   1      -6.181 -14.004  20.165  1.00  0.00           H
HETATM   50 HN1C PS  a   1      -5.787 -13.841  18.497  1.00  0.00           H
HETATM   51  C32 PS  a   1      -4.518 -13.120  19.815  1.00  0.00           C
HETATM   52  H2A PS  a   1      -3.959 -12.880  18.923  1.00  0.00           H
HETATM   53  C33 PS  a   1      -3.509 -13.556  20.905  1.00  0.00           C
HETATM   54  O35 PS  a   1      -3.735 -14.675  21.449  1.00  0.00           O
HETATM   55  O36 PS  a   1      -2.563 -12.767  21.241  1.00  0.00           O
HETATM   56  C31 PS  a   1      -5.330 -11.912  20.412  1.00  0.00           C
HETATM   57  H1A PS  a   1      -5.673 -12.640  19.646  1.00  0.00           H
HETATM   58  H1B PS  a   1      -4.269 -11.605  20.292  1.00  0.00           H
HETATM   59  P31 PS  a   1      -6.117  -9.590  19.434  1.00  0.00           P
HETATM   60  O33 PS  a   1      -4.859  -8.888  19.504  1.00  0.00           O
HETATM   61  O34 PS  a   1      -7.255  -8.747  19.647  1.00  0.00           O
HETATM   62  O32 PS  a   1      -6.201 -10.733  20.498  1.00  0.00           O
HETATM   63  O31 PS  a   1      -6.310 -10.402  18.099  1.00  0.00           O
HETATM   64  C3  PS  a   1      -7.553 -10.767  17.495  1.00  0.00           C
HETATM   65  HA  PS  a   1      -7.686 -10.145  16.584  1.00  0.00           H
HETATM   66  HB  PS  a   1      -8.440 -10.752  18.164  1.00  0.00           H
HETATM   67  C2  PS  a   1      -7.434 -12.244  17.046  1.00  0.00           C
HETATM   68  HX  PS  a   1      -7.159 -12.670  18.034  1.00  0.00           H
HETATM   69  O21 PS  a   1      -6.374 -12.372  16.159  1.00  0.00           O
HETATM   70  C21 PS  a   1      -5.780 -13.468  16.066  1.00  0.00           C
HETATM   71  O22 PS  a   1      -5.946 -14.499  16.694  1.00  0.00           O
HETATM   72  C1  PS  a   1      -8.669 -12.908  16.452  1.00  0.00           C
HETATM   73  HR  PS  a   1      -9.572 -12.745  17.077  1.00  0.00           H
HETATM   74  HS  PS  a   1      -8.470 -14.000  16.396  1.00  0.00           H
HETATM   75  O11 PS  a   1      -8.960 -12.322  15.217  1.00  0.00           O
HETATM   76  C11 PS  a   1      -9.722 -13.074  14.448  1.00  0.00           C
HETATM   77  O12 PS  a   1     -10.331 -14.103  14.689  1.00  0.00           O
HETATM   78  C12 OL  a   1      -5.038 -12.692  14.886  1.00  0.00           C
HETATM   79  H2R OL  a   1      -3.962 -12.929  15.034  1.00  0.00           H
HETATM   80  H2S OL  a   1      -5.223 -11.604  15.010  1.00  0.00           H
HETATM   81  C13 OL  a   1      -5.565 -13.291  13.565  1.00  0.00           C
HETATM   82  H3R OL  a   1      -6.674 -13.323  13.501  1.00  0.00           H
HETATM   83  H3S OL  a   1      -5.106 -14.290  13.402  1.00  0.00           H
HETATM   84  C14 OL  a   1      -5.066 -12.562  12.388  1.00  0.00           C
HETATM   85  H4R OL  a   1      -5.609 -12.965  11.507  1.00  0.00           H
HETATM   86  H4S OL  a   1      -3.976 -12.722  12.243  1.00  0.00           H
HETATM   87  C15 OL  a   1      -5.379 -11.104  12.309  1.00  0.00           C
HETATM   88  H5R OL  a   1      -4.722 -10.479  12.952  1.00  0.00           H
HETATM   89  H5S OL  a   1      -6.378 -10.838  12.716  1.00  0.00           H
HETATM   90  C16 OL  a   1      -5.328 -10.465  10.922  1.00  0.00           C
HETATM   91  H6R OL  a   1      -5.283  -9.355  10.939  1.00  0.00           H
HETATM   92  H6S OL  a   1      -6.083 -10.781  10.172  1.00  0.00           H
HETATM   93  C17 OL  a   1      -4.072 -10.915  10.088  1.00  0.00           C
HETATM   94  H7R OL  a   1      -3.803 -11.985  10.224  1.00  0.00           H
HETATM   95  H7S OL  a   1      -3.228 -10.288  10.447  1.00  0.00           H
HETATM   96  C18 OL  a   1      -4.096 -10.707   8.561  1.00  0.00           C
HETATM   97  H8R OL  a   1      -3.023 -10.549   8.317  1.00  0.00           H
HETATM   98  H8S OL  a   1      -4.539  -9.771   8.158  1.00  0.00           H
HETATM   99  C19 OL  a   1      -4.543 -11.872   7.693  1.00  0.00           C
HETATM  100  H9R OL  a   1      -4.056 -12.833   7.917  1.00  0.00           H
HETATM  101 C110 OL  a   1      -5.517 -11.984   6.700  1.00  0.00           C
HETATM  102 H10R OL  a   1      -5.665 -12.996   6.296  1.00  0.00           H
HETATM  103 C111 OL  a   1      -6.391 -10.916   6.211  1.00  0.00           C
HETATM  104 H11R OL  a   1      -6.147  -9.972   6.743  1.00  0.00           H
HETATM  105 H11S OL  a   1      -7.457 -11.176   6.387  1.00  0.00           H
HETATM  106 C112 OL  a   1      -6.043 -10.766   4.699  1.00  0.00           C
HETATM  107 H12R OL  a   1      -4.989 -10.419   4.634  1.00  0.00           H
HETATM  108 H12S OL  a   1      -6.774 -10.059   4.253  1.00  0.00           H
HETATM  109 C113 OL  a   1      -6.272 -11.989   3.860  1.00  0.00           C
HETATM  110 H13R OL  a   1      -7.325 -12.174   4.161  1.00  0.00           H
HETATM  111 H13S OL  a   1      -5.670 -12.857   4.208  1.00  0.00           H
HETATM  112 C114 OL  a   1      -5.995 -11.819   2.357  1.00  0.00           C
HETATM  113 H14R OL  a   1      -4.903 -11.898   2.168  1.00  0.00           H
HETATM  114 H14S OL  a   1      -6.282 -10.759   2.192  1.00  0.00           H
HETATM  115 C115 OL  a   1      -6.855 -12.761   1.560  1.00  0.00           C
HETATM  116 H15R OL  a   1      -7.917 -12.453   1.672  1.00  0.00           H
HETATM  117 H15S OL  a   1      -6.820 -13.811   1.919  1.00  0.00           H
HETATM  118 C116 OL  a   1      -6.528 -12.913   0.035  1.00  0.00           C
HETATM  119 H16R OL  a   1      -7.453 -13.163  -0.528  1.00  0.00           H
HETATM  120 H16S OL  a   1      -5.855 -13.754  -0.236  1.00  0.00           H
HETATM  121 C117 OL  a   1      -5.877 -11.612  -0.614  1.00  0.00           C
HETATM  122 H17R OL  a   1      -4.837 -11.475  -0.247  1.00  0.00           H
HETATM  123 H17S OL  a   1      -6.432 -10.718  -0.257  1.00  0.00           H
HETATM  124 C118 OL  a   1      -5.907 -11.586  -2.188  1.00  0.00           C
HETATM  125 H18R OL  a   1      -6.911 -11.959  -2.484  1.00  0.00           H
HETATM  126 H18S OL  a   1      -5.200 -12.295  -2.668  1.00  0.00           H
HETATM  127 H18T OL  a   1      -5.769 -10.553  -2.572  1.00  0.00           H
END}

	set f [open "POPS_A.pdb" "w"]
	puts $f "$POPS"
	close $f

	# CHL

set CHL {HETATM    1  C1  CHL A   1      -0.403  -1.003   6.051  1.00  0.00           C
HETATM    2  H11 CHL a   1       0.658  -0.841   5.873  1.00  0.00           H
HETATM    3  H12 CHL a   1      -0.618  -2.015   5.729  1.00  0.00           H
HETATM    4  C2  CHL a   1      -0.659  -0.898   7.557  1.00  0.00           C
HETATM    5  H21 CHL a   1      -1.684  -1.159   7.802  1.00  0.00           H
HETATM    6  H22 CHL a   1      -0.023  -1.600   8.087  1.00  0.00           H
HETATM    7  C3  CHL a   1      -0.386   0.507   8.061  1.00  0.00           C
HETATM    8  H31 CHL a   1       0.668   0.738   7.906  1.00  0.00           H
HETATM    9  C4  CHL a   1      -1.232   1.516   7.287  1.00  0.00           C
HETATM   10  H41 CHL a   1      -2.272   1.342   7.555  1.00  0.00           H
HETATM   11  H42 CHL a   1      -0.994   2.527   7.608  1.00  0.00           H
HETATM   12  C5  CHL a   1      -1.034   1.401   5.788  1.00  0.00           C
HETATM   13  C6  CHL a   1      -0.727   2.468   5.071  1.00  0.00           C
HETATM   14  H61 CHL a   1      -0.604   3.414   5.575  1.00  0.00           H
HETATM   15  C7  CHL a   1      -0.531   2.487   3.582  1.00  0.00           C
HETATM   16  H71 CHL a   1      -1.033   3.357   3.168  1.00  0.00           H
HETATM   17  H72 CHL a   1       0.528   2.627   3.366  1.00  0.00           H
HETATM   18  C8  CHL a   1      -1.037   1.214   2.899  1.00  0.00           C
HETATM   19  H81 CHL a   1      -2.122   1.272   2.851  1.00  0.00           H
HETATM   20  C9  CHL a   1      -0.624  -0.022   3.729  1.00  0.00           C
HETATM   21  H91 CHL a   1       0.458   0.049   3.838  1.00  0.00           H
HETATM   22  C10 CHL a   1      -1.203   0.003   5.180  1.00  0.00           C
HETATM   23  C11 CHL a   1      -0.909  -1.336   2.975  1.00  0.00           C
HETATM   24 H111 CHL a   1      -0.476  -2.170   3.518  1.00  0.00           H
HETATM   25 H112 CHL a   1      -1.976  -1.524   2.943  1.00  0.00           H
HETATM   26  C12 CHL a   1      -0.363  -1.366   1.539  1.00  0.00           C
HETATM   27 H121 CHL a   1       0.725  -1.337   1.565  1.00  0.00           H
HETATM   28 H122 CHL a   1      -0.634  -2.318   1.095  1.00  0.00           H
HETATM   29  C13 CHL a   1      -0.870  -0.177   0.706  1.00  0.00           C
HETATM   30  C14 CHL a   1      -0.470   1.097   1.483  1.00  0.00           C
HETATM   31 H141 CHL a   1       0.613   1.032   1.601  1.00  0.00           H
HETATM   32  C15 CHL a   1      -0.723   2.236   0.495  1.00  0.00           C
HETATM   33 H151 CHL a   1      -1.765   2.545   0.517  1.00  0.00           H
HETATM   34 H152 CHL a   1      -0.132   3.118   0.718  1.00  0.00           H
HETATM   35  C16 CHL a   1      -0.352   1.620  -0.872  1.00  0.00           C
HETATM   36 H161 CHL a   1      -1.137   1.794  -1.601  1.00  0.00           H
HETATM   37 H162 CHL a   1       0.544   2.079  -1.273  1.00  0.00           H
HETATM   38  C17 CHL a   1      -0.139   0.091  -0.648  1.00  0.00           C
HETATM   39 H171 CHL a   1       0.922  -0.073  -0.457  1.00  0.00           H
HETATM   40  C18 CHL a   1      -2.392  -0.280   0.481  1.00  0.00           C
HETATM   41 H181 CHL a   1      -2.763   0.516  -0.155  1.00  0.00           H
HETATM   42 H182 CHL a   1      -2.952  -0.235   1.405  1.00  0.00           H
HETATM   43 H183 CHL a   1      -2.646  -1.219  -0.000  1.00  0.00           H
HETATM   44  C19 CHL a   1      -2.703  -0.371   5.206  1.00  0.00           C
HETATM   45 H191 CHL a   1      -3.272   0.232   4.507  1.00  0.00           H
HETATM   46 H192 CHL a   1      -3.136  -0.212   6.185  1.00  0.00           H
HETATM   47 H193 CHL a   1      -2.858  -1.414   4.958  1.00  0.00           H
HETATM   48  C20 CHL a   1      -0.499  -0.745  -1.900  1.00  0.00           C
HETATM   49 H201 CHL a   1      -1.539  -0.543  -2.154  1.00  0.00           H
HETATM   50  C21 CHL a   1      -0.354  -2.257  -1.677  1.00  0.00           C
HETATM   51 H211 CHL a   1      -1.068  -2.634  -0.957  1.00  0.00           H
HETATM   52 H212 CHL a   1       0.643  -2.504  -1.319  1.00  0.00           H
HETATM   53 H213 CHL a   1      -0.518  -2.808  -2.595  1.00  0.00           H
HETATM   54  C22 CHL a   1       0.365  -0.300  -3.100  1.00  0.00           C
HETATM   55 H221 CHL a   1       0.382   0.783  -3.151  1.00  0.00           H
HETATM   56 H222 CHL a   1       1.395  -0.609  -2.925  1.00  0.00           H
HETATM   57  C23 CHL a   1      -0.096  -0.818  -4.468  1.00  0.00           C
HETATM   58 H231 CHL a   1      -0.062  -1.901  -4.491  1.00  0.00           H
HETATM   59 H232 CHL a   1      -1.137  -0.539  -4.622  1.00  0.00           H
HETATM   60  C24 CHL a   1       0.752  -0.254  -5.612  1.00  0.00           C
HETATM   61 H241 CHL a   1       0.760   0.831  -5.534  1.00  0.00           H
HETATM   62 H242 CHL a   1       1.786  -0.574  -5.487  1.00  0.00           H
HETATM   63  C25 CHL a   1       0.280  -0.632  -7.026  1.00  0.00           C
HETATM   64 H251 CHL a   1      -0.762  -0.332  -7.120  1.00  0.00           H
HETATM   65  C26 CHL a   1       0.364  -2.138  -7.295  1.00  0.00           C
HETATM   66 H261 CHL a   1      -0.273  -2.709  -6.628  1.00  0.00           H
HETATM   67 H262 CHL a   1       1.382  -2.498  -7.174  1.00  0.00           H
HETATM   68 H263 CHL a   1       0.055  -2.365  -8.312  1.00  0.00           H
HETATM   69  C27 CHL a   1       1.084   0.140  -8.078  1.00  0.00           C
HETATM   70 H271 CHL a   1       1.001   1.212  -7.931  1.00  0.00           H
HETATM   71 H272 CHL a   1       0.737  -0.085  -9.082  1.00  0.00           H
HETATM   72 H273 CHL a   1       2.139  -0.120  -8.029  1.00  0.00           H
HETATM   73  O1  CHL a   1      -0.680   0.534   9.436  1.00  0.00           O
HETATM   74  HO1 CHL a   1      -0.471   1.387   9.789  1.00  0.00           H
END}

	set f [open "CHL_A.pdb" "w"]
	puts $f "$CHL"
	close $f

	set  k 0
	while { $k < [llength $p1] } {
		set count 0
		set k1 0

		while { [lindex $inp $k1] != "NA" } {
			if { [lindex $inp $k1] == [lindex $p1 $k] } {
				incr count
			}
			incr k1
		}
		if { $count == 0 } {
			#puts "[lindex $p1 $k]\_A.pdb"
			file delete [lindex $p1 $k]\_A.pdb
		}
		incr k
	} 
}
 
proc inp_grid {} {
	set g {POPC 
{ {A} { {0.0 0.0 0.50} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
{ {S} { {0.0 0.0 0.40} {0.0 0.0 0.70} } { {0.0 0.0 0.70} {0.0 0.0 0.40} } { {0.0 0.0 0.0} } } 
DLPC
{ {A} { {0.0 0.0 0.30} } { {0.0 0.0 0.30} } { {0.0 0.0 0.0} } }
{ {A} { {0.0 0.0 1.10} } { {0.0 0.0 1.10} } { {0.0 0.0 0.0} } }
DOPC
{ {S1} { {0.0 0.0 0.9} } { {0.0 0.0 1.00} } { {0.0 0.5 0.3} } { {0.0 0.0 5.9 } } } 
{ {A} { {0.0 0.0 0.20} } { {0.0 0.0 0.20} } { {0.0 0.0 0.20} } }
DOPG
{ {S} { {0.0 0.0 0.3} } { {0.0 0.0 3.70} } { {0.0 0.0 0.0} } }
{ {S2} { {0.0 0.0 0.1} } { {0.0 0.0 0.3} } { {0.0 0.0 0.2} } { {0.0 0.0 0.0} } }
DOPS
{ {A} { {0.0 0.0 0.4} } { {0.0 0.0 0.4} } { {0.0 0.0 0.4} } }
{ {A1} { { {0.0 0.0 0.9} } { {0.0 0.5 0.8} } } { { {0.0 0.0 0.9} } { {0.0 0.5 0.8} } } { { {0.0 0.0 0.7}} { {0.0 0.5 0.4} } } }
DPPC
{ {A} { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } }
{ {A} { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } }
POPG
{ {A} { {0.0 0.0 0.50} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
{ {S} { {0.0 0.0 0.40} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
POPG
{ {A} { {0.0 0.0 0.50} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
{ {S} { {0.0 0.0 0.40} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
POPE
{ {A} { {0.0 0.0 0.40} } { {0.0 0.0 0.40} } { {0.0 0.0 0.0} } }
{ {A} { {0.0 0.0 0.50} } { {0.0 0.0 0.50} } { {0.0 0.0 0.0} } }
DMPC
{ {A} { {0.0 0.0 0.10} } { {0.0 0.0 0.10} } { {0.0 0.0 0.0} } }
{ {S3} { {0.0 0.0 5.4} } { {0.0 0.0 5.4} } { {0.0 0.0 4.8} } { {0.0 0.0 4.8 } } { {0.0 0.0 0.0} } } 
DLPE
{ {A} { {0.0 0.0 0.3} } { {0.0 0.0 0.3} } { {0.0 0.0 0.0} } }
{ {A} { {0.0 0.0 1.1} } { {0.0 0.0 1.1} } { {0.0 0.0 0.0} } }
CHL
{ {A} { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } }
{ {A} { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } { {0.0 0.0 0.0} } }
}    

	set f [open "grid" "w"]
	
	puts $f "$g"
	close $f
}

proc solvation {opt} {
	set f [open "input_2" "w"]
	
	puts $f "\{ LEAP INPUT \}"

	puts "				#### ENTER THE NAME OF THE PARAMETER FILES (EACH SEPARATED BY A SPACE) ####"
	puts ""
	set inp1 [gets stdin]

	puts ""

	puts $f "\{ PARM: $inp1 \}"

	puts "				#### IS THERE ANY ADDITIONAL PARAMETER FILE? (frcmod file) (y/n) ####"
	puts ""
	set inp2 [gets stdin]

	if { $inp2 == "Y" || $inp2 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp3 [gets stdin]
	} else {
		set inp3 ""
	}

	puts $f "\{ frcmod: $inp3 \}"

	puts "				#### IS THERE ANY .LIB OR .OFF STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp4 [gets stdin]

	if { $inp4 == "Y" || $inp4 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp5 [gets stdin]
	} else {
		set inp5 ""
	}

	puts $f "\{ lib: $inp5 \}"

	puts "				#### IS THERE ANY .prepin STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp6 [gets stdin]

	if { $inp6 == "Y" || $inp6 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp7 [gets stdin]
	} else {
		set inp7 ""
	}

	puts $f "\{ prep: $inp7 \}"


	puts "				#### ENTER THE THICKNESS OF THE WATER LAYER IN ANGSTOMS YOU WANT TO ADD ####"
	puts ""
	set inp8 [gets stdin]

	puts $f "\{ water_layer: $inp8 \}"

	close $f

	# READING THE INPUT FILE

	set in [open "input_2" "r"]
	set inp [read $in]
	close $in

	set in1 [open "input" "r"]
	set inp1 [read $in1]
	close $in1

	set hm [lindex $inp1 10 1]

	if { [llength [lindex $inp1 7]] > 1 } {
		set pro [open "[lindex $inp1 7 1]" "r"]
		set pr [read $pro]
		close $pro

		set k 0

		while { [lindex $pr $k] != "ATOM" && [lindex $pr $k] != "HETATM"} {
			incr k
		}

		set resc [lindex $pr [expr { $k + 3 }]] 
	} else { 
		set resc 1
	}

	puts "STEP 2 :: LET LEAP PUT WATER IN"

	# FORMIING THE LEAP FILE

	set l [open "leap.in" "w"]

	set l1 [open "leap_1.in" "w"]

	puts $l "logfile leap.log"
	puts $l1 "logfile leap_1.log"

	# PARAMETER FILE
	set k 1
	while { $k < [llength [lindex $inp 1]] } {
		puts $l "source leaprc.[lindex $inp 1 $k]"
		puts $l1 "source leaprc.[lindex $inp 1 $k]"
		incr k
	}  

	# FRCMOD FILE
	set k 1
	while { $k < [llength [lindex $inp 2]] } {
		puts $l "loadamberparams [lindex $inp 2 $k]"
		puts $l1 "loadamberparams [lindex $inp 2 $k]"
		incr k
	}

	# LIB FILE
	set k 1
	while { $k < [llength [lindex $inp 3]] } {
		puts $l "loadoff [lindex $inp 3 $k]"
		puts $l1 "loadoff [lindex $inp 3 $k]"
		incr k
	}

	# PREPIN FILE
	set k 1
	while { $k < [llength [lindex $inp 4]] } {
		puts $l "loadamberprep [lindex $inp 4 $k]"
		puts $l1 "loadamberprep [lindex $inp 4 $k]"
		incr k
	}

	# PDB FILE
	if { $opt == 0 } {
		puts $l "mol=loadpdb \"lipids_no.pdb\""
	} else { 
		if { $hm == 1 } {
			puts $l "mol=loadpdb \"1_final_struct.pdb\""
		} else { 
			puts $l "mol=loadpdb \"lip_pro.pdb\""
		}
	}

	puts $l "solvatebox mol TIP3PBOX [lindex $inp 5 1]"

	puts $l "savepdb mol refm1.pdb"
	puts $l "quit"
	close $l

	exec tleap -f leap.in


	puts "STEP 3 :: FORMING THE FILE WITH JUST PROTEINS AND LIPIDS"

	set f [open "refm1.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "only_pro.pdb" "w"]

	set k 0
	while { [lindex $data $k] != "WAT" } {
		incr k
	}

	set strng [string first "WAT" $data]
	set strng [expr { $strng - 20 }]
	set data2 [string replace $data $strng end ""]

	puts $g "$data2"
	close $g	

	puts "STEP 4 :: SELECTING THE WATER ON THE TOP AND BOTTOM LAYER"
	
	wat_sel
	wat_sel_ps

	puts "STEP 5 :: REMOVING THE ERROR IN THE PDB FILE INTRODUCED BY LEAP"

	leap_error

	puts "STEP 6 :: REPLACING THE PDB FILE WITH THE MODIFIED LIPID FILE"

	set h1 [open "nd1_d2.pdb" "r"]
	set data1 [read $h1]
	close $h1

	if { $resc != 1 } {
		set strng [string first "$resc" $data1]
		set strng [expr { $strng - 18 }]
		set data2 [string replace $data1 0 $strng ""]
	} else {
		set strng [string first "WAT" $data1]
		set strng [expr { $strng - 18 }]
		set data2 [string replace $data1 0 $strng ""]
	}

	set h2 [open "mod_lipids.pdb" "r"]
	set data3 [read $h2]
	close $h2

	set data4 "$data3$data2"

	set h1 [open "nd1_d2.pdb" "w"]

	puts $h1 "$data4"
	close $h1

	puts "STEP 7 :: GETTING THE BOX DIMENSIONS"

	box_dimensions "nd1_d2.pdb"
		
	puts ""

	puts "STEP 8 :: FORMING THE FINAL LEAP FILE TO GET THE FINAL PARAMETERS AND THE COORDINATES"

	puts $l1 "mol=loadpdb \"nd1_d2.pdb\""

	set d [open "box_dimensions" "r"]
	set bd [read $d]
	close $d

	puts $l1 "set mol box { [lindex $bd 0] [lindex $bd 1] [lindex $bd 2] }"

	puts ""
	puts "				#### DO YOU WANT TO MANUALLY ENTER THE NUMBER OF IONS (Y/N) ####"
	puts ""
	set ionans [gets stdin]
	if { $ionans == "y" || $ionans == "Y" } {
		puts ""
		puts "			#### ENTER THE NUMBER OF CHLORIDE IONS ####"
		set cl [gets stdin]
		puts ""
		puts "			#### ENTER THE NUMBER OF POTASSIUM IONS ####"
		set K [gets stdin]
		puts ""
		puts $l1 "addions mol Cl- $cl"
		puts $l1 "addions mol K+ $K"
	} else {
		puts $l1 "addions mol Cl- 0.0"
		puts $l1 "addions mol K+ 0.0"
	}

	puts $l1 "saveamberparm mol mol.prmtop mol.inpcrd"
	puts $l1 "savepdb mol reff.pdb"
	puts $l1 "quit"

	close $l1

	exec tleap -f leap_1.in
	
}

proc solvation_100000 {opt} {

	set f [open "input_2" "w"]
	
	puts $f "\{ LEAP INPUT \}"

	puts "				#### ENTER THE NAME OF THE PARAMETER FILES (EACH SEPARATED BY A SPACE) ####"
	puts ""
	set inp1 [gets stdin]

	puts ""

	puts $f "\{ PARM: $inp1 \}"

	puts "				#### IS THERE ANY ADDITIONAL PARAMETER FILE? (frcmod file) (y/n) ####"
	puts ""
	set inp2 [gets stdin]

	if { $inp2 == "Y" || $inp2 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp3 [gets stdin]
	} else {
		set inp3 ""
	}

	puts $f "\{ frcmod: $inp3 \}"

	puts "				#### IS THERE ANY .LIB OR .OFF STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp4 [gets stdin]

	if { $inp4 == "Y" || $inp4 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp5 [gets stdin]
	} else {
		set inp5 ""
	}

	puts $f "\{ lib: $inp5 \}"

	puts "				#### IS THERE ANY .prepin STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp6 [gets stdin]

	if { $inp6 == "Y" || $inp6 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp7 [gets stdin]
	} else {
		set inp7 ""
	}

	puts $f "\{ prep: $inp7 \}"


	puts "				#### ENTER THE THICKNESS OF THE WATER LAYER IN ANGSTOMS YOU WANT TO ADD ####"
	puts ""
	set inp8 [gets stdin]

	puts $f "\{ water_layer: $inp8 \}"

	close $f

	# READING THE INPUT FILE

	set in [open "input_2" "r"]
	set inp [read $in]
	close $in

	set in1 [open "input" "r"]
	set inp1 [read $in1]
	close $in1

	set hm [lindex $inp1 10 1]

	if { [llength [lindex $inp1 7]] > 1 } {
		set pro [open "[lindex $inp1 7 1]" "r"]
		set pr [read $pro]
		close $pro

		set k 0

		while { [lindex $pr $k] != "ATOM" && [lindex $pr $k] != "HETATM"} {
			incr k
		}

		set resc [lindex $pr [expr { $k + 3 }]] 
	} else { 
		set resc 1
	}

	puts "STEP 1 :: LET LEAP PUT WATER IN"

	# FORMIING THE LEAP FILE

	set l [open "leap.in" "w"]

	set l1 [open "leap_1.in" "w"]

	puts $l "logfile leap.log"
	puts $l1 "logfile leap_1.log"

	# PARAMETER FILE
	set k 1
	while { $k < [llength [lindex $inp 1]] } {
		puts $l "source leaprc.[lindex $inp 1 $k]"
		puts $l1 "source leaprc.[lindex $inp 1 $k]"
		incr k
	}  

	# FRCMOD FILE
	set k 1
	while { $k < [llength [lindex $inp 2]] } {
		puts $l "loadamberparams [lindex $inp 2 $k]"
		puts $l1 "loadamberparams [lindex $inp 2 $k]"
		incr k
	}

	# LIB FILE
	set k 1
	while { $k < [llength [lindex $inp 3]] } {
		puts $l "loadoff [lindex $inp 3 $k]"
		puts $l1 "loadoff [lindex $inp 3 $k]"
		incr k
	}

	# PREPIN FILE
	set k 1
	while { $k < [llength [lindex $inp 4]] } {
		puts $l "loadamberprep [lindex $inp 4 $k]"
		puts $l1 "loadamberprep [lindex $inp 4 $k]"
		incr k
	}

	# PDB FILE
	if { $opt == 0 } {
		puts $l "mol=loadpdb \"lipids_no.pdb\""
	} else { 
		if { $hm == 1 } {
			puts $l "mol=loadpdb \"1_final_struct.pdb\""
		} else { 
			puts $l "mol=loadpdb \"lip_pro.pdb\""
		}
	}
	
	puts ""
	puts "				#### DO YOU WANT TO MANUALLY ENTER THE NUMBER OF IONS (Y/N) ####"
	puts ""
	set ionans [gets stdin]
	if { $ionans == "y" || $ionans == "Y" } {
		puts ""
		puts "			#### ENTER THE NUMBER OF CHLORIDE IONS ####"
		set cl [gets stdin]
		puts ""
		puts "			#### ENTER THE NUMBER OF POTASSIUM IONS ####"
		set K [gets stdin]
		puts ""
		puts $l1 "addions mol Cl- $cl"
		puts $l1 "addions mol K+ $K"
	} else {
		puts $l1 "addions mol Cl- 0.0"
		puts $l1 "addions mol K+ 0.0"
	}

	puts $l "solvatebox mol TIP3PBOX \{ 0.0 0.0 [lindex $inp 5 1] \}"
	puts $l "saveamberparm mol mol1.prmtop mol1.inpcrd"
	puts $l "savepdb mol refm1.pdb"
	puts $l "quit"
	close $l

	exec tleap -f leap.in

	puts "STEP 2 :: FORMING THE FILE WITH JUST PROTEINS AND LIPIDS"

	set f [open "refm1.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "only_pro.pdb" "w"]

	set k 0
	while { [lindex $data $k] != "WAT" } {
		incr k
	}

	set strng [string first "WAT" $data]
	set strng [expr { $strng - 20 }]
	set data2 [string replace $data $strng end ""]

	puts $g "$data2"
	close $g	

	puts "STEP 3 :: SELECTING THE UNWANTED WATER MOLECULES"
	
	wat_sel
	wat_sel_ps

	file_overlap
}

proc wat_sel {} {
	# THIS CODE READS THE INITIAL PDB FILE GENERATED BY LEAP AND SHIFT THE DEFAULT RADIUS TO THE ATOM OF THE FIRST RESIDUE AND THEN REMOVES THE WATER MOLECULES ALONG THE BILYER LEAFLETS

	set in [open "input" "r"]
	set inp [read $in]
	close $in

	set f [open "refm1.pdb" "r"]
	set data [read $f]
	close $f

	# SHIFTING THE ORIGIN TO THE FIRST ATOM OF THE FIRST RESIDUE, WHICH IS THE LIPID BELONGING TO THE EXTREME LEFT OF THE UPPER LEAFLET

	set k 0
	set xmin 1000.0
	set ymin 1000.0
	set zmax 0.000
	while { $k < [llength $data] } {
		if { [lindex $data $k] == "[lindex $inp 2 1]" } {
			set x1 [lindex $data [expr { $k + 3 }]]
			set sx1 [string length $x1]

			if { $sx1 > 8 } {
				set t 0
				while { [string range $x1 $t $t] != "." } {
					incr t
				}
				set xori [string range $x1 0 [expr { $t + 3 }]]
				set yori [string range $x1 [expr { $t + 4 }] end]
				set zori [lindex $data [expr { $k + 4 }]]
			} else { 
				set xori $x1
				set y1 [lindex $data [expr { $k + 4 }]]
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set yori [string range $y1 0 [expr { $t + 3 }]]
			 	 	set zori [string range $y1 [expr { $t + 4 }] end]
				}	else {
					set yori [lindex $data [expr { $k + 4 }]]
					set zori [lindex $data [expr { $k + 5 }]]
				}
			}
			if { $xori < $xmin } {
				set xmin $xori
			} 
			if { $yori < $ymin } {
				set ymin $yori
			}
			if { $zori > $zmax } {
				set zmax $zori
			}
		}
		incr k
	}
	set xori $xmin
	set yori $ymin
	set zori $zmax
	puts "$xori $yori $zori"

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	set p(5) ""
	set p(6) ""
	
	set p1(1) "    "
	set p1(2) "   "
	set p1(3) "  "
	set p1(4) " "
	set p1(5) ""
	set p1(6) ""

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

	# REMOVING THE WATER MOLECULE FROM THE LEFT HALF (X < 0.0, -t < Y < 0.0) AND RIGHT HALF (X > Xl, -t < Y < 0.0), HERE Xl = 150.0 AND t = 18.0

	set h [open "water.pdb" "w" ]
	set hh [open "water_resid" "w"]

	set k 0 
	set i 1
	set j 1
	set atype ""
	set count 0
	while { $k < [llength $data] } {
		set term [lindex $data $k]
		set t1 [string range $term 0 5]
		if { $atype == "H2" && $count != 0} {
			puts $h "TER"
			set count 0
			incr j
		}
		if { [lindex $data $k] == "ATOM" || $t1 == "HETATM" } {
			if { [lindex $data [expr { $k + 3 }]] == "WAT" } { 
				if { $t1 == "HETATM" } {
					set sterm [string length $term]
					if { $sterm > 6 } {
						set shift 1
						set rn1 $term
						set srn1 5
						set ft ""
					} else {
							set shift 0
							#set rn1 [lindex $data [expr { $k + 1 }]]
							set rn1 $i
							set srn1 [string length $rn1]
							set ft "HETATM"
					}
				} else { 
						set shift 0
						#set rn1 [lindex $data [expr { $k + 1 }]]
						set rn1 $i
						set srn1 [string length $rn1]
						if { $srn1 == 6 } {
							set ft "ATOM "
						} else { 
							set ft "ATOM  "
						}
				}
			
				set atype [lindex $data [expr { $k + 2 - $shift}]]
				set satype [string length $atype]

				if { $satype > 5} {
					set shift2 1
					set at1 [lindex $data [expr { $k + 2 - $shift }]]		
					set sat1 4
					set resn ""
					set sresn 4
				} else {
					set shift2 0
					set at1 [lindex $data [expr { $k + 2 - $shift }]]		
					set sat1 [string length $at1]
					set resn [lindex $data [expr { $k + 3 - $shift}]]
					set sresn [string length $resn]
				}			

				set an1 [lindex $data [expr { $k + 4 -$shift - $shift2 }]]
				#set an1 $j			
				set san1 [string length $an1]
				set shift1 0

				set x1 [lindex $data [expr { $k + 5 - $shift - $shift1 -$shift2}]]
				set sx1 [string length $x1]
				if { $sx1 > 8 } {
					set shift3 1
					set t 0
					while { [string range $x1 $t $t] != "." } {
						incr t
					}
					set corx [string range $x1 0 [expr { $t + 3 }]]
					set cory [string range $x1 [expr { $t + 4 }] end]
					set x1 [expr { $corx - $xori } ]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

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
	
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [expr { $corz - $zori }]
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					} else {
						set shift4 0
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [lindex $data [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]] 
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					}
				} else { 
					set shift3 0
					set x1 [format "%.3f" [expr { $x1 - $xori }]]
					set sx1 [string length $x1]
					set y1 [lindex $data [expr { $k + 6 - $shift -$shift1 - $shift2 - $shift3}]] 
					set sy1 [string length $y1]
					if { $sy1 > 8 } {
						set shift4 1
						set t 0
						while { [string range $y1 $t $t] != "." } {
							incr t
						}
						set cory [string range $y1 0 [expr { $t + 3 }]]
						set corz [string range $y1 [expr { $t + 4 }] end]
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [expr { $corz - $zori }]
						set z1 [format "%.3f" $z1]
						set sz1 [string length $z1]
					} else {
						set shift4 0
						set y1 [format "%.3f" [expr { $y1 - $yori }]]
						set sy1 [string length $y1]
						set z1 [lindex $data [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					}
				}
		
				# PARAMTERS HERE NEED TO BE TUNED TO GET THE CORRECT WATER LAYER AROUND THE BILAYER :: BASED ON BILYER CROSS SECTION AND THICKNESS

				set box_x [lindex $inp 5 1]
				set box_y [lindex $inp 5 2]
				set box_z [lindex $inp 5 3]
				set ext 5.0
				if { $z1 > $ext || $z1 < [expr { (-2*$box_z) + (-1*$ext) }] } {
					if { $x1 > [expr { -1*$ext }] && $x1 <= [expr { $box_x - $ext }] } {
						if { $y1 > [expr { -1*$ext }] && $y1 <= [expr { $box_y - $ext }]} {
							set x1 [format "%.3f" [expr { $x1 + $xori }]]
							set sx1 [string length $x1]
							set y1 [format "%.3f" [expr { $y1 + $yori }]]
							set sy1 [string length $y1]
							set z1 [format "%.3f" [expr { $z1 + $zori }]]
							set sz1 [string length $z1]
							incr count
							incr i
							if { $srn1 > 6 } {
								set ft "ATOM "
								set sat2 4
								set srn1 6
							} else {
								set sat2 $sat1
							}
							puts $h "$ft$p1($srn1)$rn1 $ic($sat2)$at1$sat($sat1)$p($sresn)$resn  $an1$p1($san1)   $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00"
						} else { 
							puts "			**** REMOVING WATER $an1 ****"
							puts $hh "$an1"
						}
					} else { 
						puts "			**** REMOVING WATER $an1 ****"
						puts $hh "$an1"
					}
				} else {
					puts "			**** REMOVING WATER $an1 ****"
					puts $hh "$an1"
				}
			}
		}
		incr k
	}
	puts $h "END"
	close $h
	close $hh

	# COMBINING THE PDB CONTATING THE PROTEIN AND LIPIDS AND THE WATER OF INTEREST

	set g1 [open "only_pro.pdb" "r"]
	set dat1 [read $g1]
	close $g1

	set g2 [open "water.pdb" "r"]
	set dat2 [read $g2]
	close $g2

	set g3 [open "nd1_d2.pdb" "w"]

	set dat3 "$dat1$dat2"

	puts $g3 "$dat3"

	close $g3
}

proc wat_sel_ps {} {

	set f [open "water.pdb" "r" ]
	set data [read $f]
	close $f

	set k 0
	set i 0
	while { $k < [llength $data] } {
		if { [lindex $data $k] == "ATOM" } {
			set elem1 [lindex $data [expr { $k + 2 }]]

			if { $elem1 == "O" } {
				set k1 [expr {$k + 1 }]

				while { [lindex $data $k1] != "ATOM" && [lindex $data $k1] != "TER" } {
					incr k1
				}
				set elem2 [lindex $data [expr { $k1 + 2 }]] 

				if { $elem2 == "H1" } {
					set k2 [expr { $k1 + 1}]
					set k $k2
				
					while { [lindex $data $k2] != "ATOM" && [lindex $data $k2] != "TER" } {
						incr k2
					}
					set elem3 [lindex $data [expr { $k2 + 2 }]] 
			
					if { $elem3 == "H2" } {
						set k [expr { $k2 + 1 }]
					} else {
						set bl($i) [lindex $data [expr { $k + 3 }]]
						incr i
					}
				} else {
						set bl($i) [lindex $data [expr { $k + 4 }]]
						incr i 
				}
			} else { 
					set bl($i) [lindex $data [expr { $k + 4 }]]
					incr i 
			}
		}
		incr k
	}

	set nbl $i
	puts "$nbl"
	after 1000
				
	set xori 0.0
	set yori 0.0
	set zori 0.0

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	set p(5) ""
	set p(6) ""
	
	set p1(1) "    "
	set p1(2) "   "
	set p1(3) "  "
	set p1(4) " "
	set p1(5) ""
	set p1(6) ""

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

	# REMOVING THE WATER MOLECULE FROM THE LEFT HALF (X < 0.0, -t < Y < 0.0) AND RIGHT HALF (X > Xl, -t < Y < 0.0), HERE Xl = 150.0 AND t = 18.0

	set f [open "water.pdb" "r" ]
	set data [read $f]
	close $f

	set h [open "mod_water.pdb" "w"]

	set h1 [open "mod_water_vis.pdb" "w"]

	set hh [open "water_resid1" "w"]

	set k 0 
	set i 1
	set j 1
	set atype ""
	set count 0
	while { $k < [llength $data] } {
		set term [lindex $data $k]
		set t1 [string range $term 0 5]
		if { $atype == "H2" && $count != 0} {
			puts $h "TER"
			set count 0
			incr j
			incr i
		}
		if { [lindex $data $k] == "ATOM" || $t1 == "HETATM" } {
			if { [lindex $data [expr { $k + 3 }]] == "WAT" } {
				if { $i > 99999 } {
					set i 1
				}
				if { $t1 == "HETATM" } {
					set sterm [string length $term]
					if { $sterm > 6 } {
						set shift 1
						set rn1 $term
						set srn1 5
						set ft ""
					} else {
							set shift 0
							set rn1 [lindex $data [expr { $k + 1 }]]
							#set rn1 $i
							set srn1 [string length $rn1]
							set ft "HETATM"
					}
				} else { 
						set shift 0
						set rn1 [lindex $data [expr { $k + 1 }]]
						#set rn1 $i
						set srn1 [string length $rn1]
						if { $srn1 == 6 } {
							set ft "ATOM "
						} else { 
							set ft "ATOM  "
						}
				}
			
				set atype [lindex $data [expr { $k + 2 - $shift}]]
				set satype [string length $atype]

				if { $satype > 5} {
					set shift2 1
					set at1 [lindex $data [expr { $k + 2 - $shift }]]		
					set sat1 4
					set resn ""
					set sresn 4
				} else {
					set shift2 0
					set at1 [lindex $data [expr { $k + 2 - $shift }]]		
					set sat1 [string length $at1]
					set resn [lindex $data [expr { $k + 3 - $shift}]]
					set sresn [string length $resn]
				}			

				set an1 [lindex $data [expr { $k + 4 -$shift - $shift2 }]]
				set an1n $i
				set san1n [string length $an1n]
				set san1 [string length $an1]
				set shift1 0

				set x1 [lindex $data [expr { $k + 5 - $shift - $shift1 -$shift2}]]
				set sx1 [string length $x1]
				if { $sx1 > 8 } {
					set shift3 1
					set t 0
					while { [string range $x1 $t $t] != "." } {
						incr t
					}
					set corx [string range $x1 0 [expr { $t + 3 }]]
					set cory [string range $x1 [expr { $t + 4 }] end]
					set x1 [expr { $corx - $xori } ]
					set x1 [format "%.3f" $x1]
					set sx1 [string length $x1]

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
	
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [expr { $corz - $zori }]
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					} else {
						set shift4 0
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [lindex $data [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3}]] 
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					}
				} else { 
					set shift3 0
					set x1 [format "%.3f" [expr { $x1 - $xori }]]
					set sx1 [string length $x1]
					set y1 [lindex $data [expr { $k + 6 - $shift -$shift1 - $shift2 - $shift3}]] 
					set sy1 [string length $y1]
					if { $sy1 > 8 } {
						set shift4 1
						set t 0
						while { [string range $y1 $t $t] != "." } {
							incr t
						}
						set cory [string range $y1 0 [expr { $t + 3 }]]
						set corz [string range $y1 [expr { $t + 4 }] end]
						set y1 [expr { $cory - $yori }]
						set y1 [format "%.3f" $y1]
						set sy1 [string length $y1]
						set z1 [expr { $corz - $zori }]
						set z1 [format "%.3f" $z1]
						set sz1 [string length $z1]
					} else {
						set shift4 0
						set y1 [format "%.3f" [expr { $y1 - $yori }]]
						set sy1 [string length $y1]
						set z1 [lindex $data [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3 - $shift4}]] 
						set z1 [format "%.3f" [expr { $z1 - $zori }]]
						set sz1 [string length $z1]
					}
				}
				set count1 0

				for {set ij 0} {$ij < $nbl} {incr ij} {
					if { $an1 == $bl($ij) } {
						incr count1
					}
				}

				if { $count1 == 0 } {
					if { $srn1 > 6 } {
						set ft "ATOM "
						set sat2 4
						set srn1 6
					} else {
						set sat2 $sat1
					}
					puts $h "$ft$p1($srn1)$rn1 $ic($sat2)$at1$sat($sat1)$p($sresn)$resn  $an1$p1($san1)   $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00"
					puts $h1 "$ft$p1($srn1)$rn1 $ic($sat2)$at1$sat($sat1)$p($sresn)$resn  $an1n$p1($san1n)   $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00"
					incr count
				} else { 
					puts "			**** REMOVING WATER $an1 ****"
					puts $hh "$an1"
				}
			}
		}
		incr k
	}
	puts $h "END"
	close $h
	close $hh

	puts $h1 "END"
	close $h1

	set g1 [open "only_pro.pdb" "r"]
	set dat1 [read $g1]
	close $g1

	set g2 [open "mod_water.pdb" "r"]
	set dat2 [read $g2]
	close $g2

	set g3 [open "nd1_d2.pdb" "w"]
	set g5 [open "nd1_d2_vis.pdb" "w"]

	set g4 [open "mod_water_vis.pdb" "r"]
	set dat4 [read $g4]
	close $g4

	set dat3 "$dat1$dat2"
	set dat5 "$dat1$dat4"

	puts $g3 "$dat3"
	puts $g5 "$dat5"

	close $g3
	close $g5
}

proc leap_error {} {
	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""


	set f [open "nd1_d2.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "mod_lipids.pdb" "w"]

	set in [open "input" "r"]
	set inp [read $in]
	close $in

	if { [llength [lindex $inp 7]] > 1 } {

		set pro [open "[lindex $inp 7 1]" "r"]
		set pr [read $pro]
		close $pro

		set k 0

		while { [lindex $pr $k] != "ATOM" && [lindex $pr $k] != "HETATM"} {
			incr k
		}

		set resc [lindex $pr [expr { $k + 3 }]] 
	} else { 
		
		set k [llength $data]

		while { [lindex $data $k] != "ATOM" && [lindex $data $k] != "HETATM" } {
			set k [expr { $k - 1 }]
		}

		set resc [lindex $data [expr { $k + 3 }]] 
	}
	puts "$resc"

	set k 0
	set i 1
	set j 1
	set num2 0
	set resname ""
	while { [lindex $data [expr { $k + 3 }]] != "$resc" } {
		if { $resname == "CHL" || $resname == "SM" } {
			if { [lindex $data $k] == "TER" } {
				puts $g "TER"
				incr i
				incr num2
			}
		} elseif { [lindex $data $k] == "TER" && [expr { ($j-$num2) % 3 }] == 0 } {
			puts $g "TER" 
			incr i
		}
	
		if { [lindex $data $k] == "ATOM" } {

			set j [lindex $data [expr { $k + 4 }]]

			set rn1 [lindex $data [expr { $k + 1 }]]
			set srn1 [string length $rn1]

			set at1 [lindex $data [expr { $k + 2 }]]
			set sat1 [string length $at1]

			set an1 $i
			set san1 [string length $an1]

			set x1 [lindex $data [expr { $k + 5 }]]
			set sx1 [string length $x1]
			if { $sx1 > 8 } {
				set shift3 1
				set t 0
				while { [string range $x1 $t $t] != "." } {
					incr t
				}
				set corx [string range $x1 0 [expr { $t + 3 }]]
				set cory [string range $x1 [expr { $t + 4 }] end]
				set x1 $corx
				set x1 [format "%.3f" $x1]
				set sx1 [string length $x1]

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
	
					set y1 $cory 
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]
					set z1 $corz
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
				} else {
					set shift4 0
					set y1 $cory
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]
					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
				}
			} else { 
				set shift3 0
				set x1 [format "%.3f" $x1]
				set sx1 [string length $x1]
				set y1 [lindex $data [expr { $k + 6 }]] 
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set shift4 1
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set cory [string range $y1 0 [expr { $t + 3 }]]
					set corz [string range $y1 [expr { $t + 4 }] end]
					set y1 $cory
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]
					set z1 $corz
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
				} else {
					set shift4 0
					set y1 [format "%.3f" $y1]
					set sy1 [string length $y1]
					set z1 [lindex $data [expr { $k + 7 }]] 
					set z1 [format "%.3f" $z1]
					set sz1 [string length $z1]
				}
			}

			set resname [lindex $data [expr { $k + 3 }]]
			set sresname [string length [lindex $data [expr { $k + 3 }]]]
	
			puts $g "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)[lindex $data [expr { $k + 3 }]]$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           "

		}
		incr k
	}
	close $g
}

proc solvation_nc {} {
	# THIS CODE PUTS THE WATER MOLECULE IN A BILAYER SYSTEM WITH NO PROTEIN INSIDE 

	# THE VARIABLES "ext" and "wat_lay" ARE THE IMPORTANT VARIABLE WHICH SHOULD BE CHANGED FIRST IN CASE THE WATER LAYER BEHAVE UNEXPECTEDLY 

	# ext = EXTENSION OF THE HEAD GROUP ABOVE THE PIVOT ATOM AND 
	# wat_lay = THE DISTANCE OF WATER LAYER FROM THE TOP OF THE BILAYER TO THE SURFACE OF THE BOX 

	sol_nc_inp

	set inp [open "water_input" "r"]
	set mn [read $inp]
	close $inp

	set nwat [lindex $mn 0 1]
	set ext [lindex $mn 2 1]
	set wat_lay [lindex $mn 1 1]

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""

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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set f [open "lipids_no.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "wat.pdb" "r"]
	set data1 [read $g]
	close $g

	set h [open "input" "r"]
	set data2 [read $h]
	close $h

	set w [open "nwat.pdb" "w"]

	set k 0
	set i 0

	while { $k < [llength $data1] } {
		if { [lindex $data1 $k] == "HETATM" || [lindex $data1 $k] == "ATOM" } {
			set xwat($i) [lindex $data1 [expr { $k + 5 }]]
			set ywat($i) [lindex $data1 [expr { $k + 6 }]]
			set zwat($i) [lindex $data1 [expr { $k + 7 }]]
			set at1($i) [lindex $data1 [expr { $k + 2 }]]
			incr i
		}
		incr k
	}

	set ty(0) "O"
	set ty(1) "H"
	set ty(2) "H"

	set lx [lindex $data2 5 1]
	set ly [lindex $data2 5 2]
	set lz [lindex $data2 5 3]

	set box_x [expr { $lx + 0.0 }]
	set box_y [expr { $ly + 0.0 }]
	set box_z [expr { (2*$ext)+ (2*$wat_lay)}]

	puts "				_____________________________________"
	puts "				|                                    |"
	puts "					WATER LAYER = $wat_lay Ang"
	puts "				|                                    |"
	puts "				....................................."
	puts "					EXTENSION = $ext Ang"
	puts "				|                                    |"
	puts "				-------------------------------------"
	puts "				|                                    |"
	puts "				|                                    |"
	puts "				   LIPID BILAYER = [expr { 2 * $lz }] Ang thick"
	puts "				|                                    |"
	puts "				|                                    |"
	puts "				-------------------------------------"
	puts "				|                                    |"
	puts "					EXTENSION = $ext Ang"
	puts "				-------------------------------------"
	puts "					WATER LAYER = $wat_lay Ang"
	puts "				|____________________________________|"


	set j 1
	set l 1
	set nwatin $j
	# UPPER WATER LAYER
	while { $nwatin < [expr { $nwat / 2 }] } {
		set shiftx($j) [expr { rand() * $box_x }]
		set shiftx1 [expr { rand() * $ext }]
		set shiftx($j) [expr { $shiftx($j) - $shiftx1 }]
		set shifty($j) [expr { rand() * $box_y }]
		set shifty1 [expr { rand() * $ext }]
		set shifty($j) [expr { $shifty($j) - $shifty1 }]
		set shiftz($j) [expr { rand() * $box_z }]
		set count 0
		for {set m 1} {$m < $j} {incr m} {
			set delx [expr { abs($shiftx($j) - $shiftx($m)) }]
			set delx2 [expr { $delx * $delx }]
			set dely [expr { abs($shifty($j) - $shifty($m)) }]
			set dely2 [expr { $dely * $dely }]
			set delz [expr { abs($shiftz($j) - $shiftz($m)) }]
			set delz2 [expr { $delz * $delz }]
			set del [expr { $delx2 + $dely2 + $delz2 }]
			if { $del < 0.5 } {
				incr count 
			}
		}

		if { $shiftz($j) > $ext && $count == 0} {
			puts "				**** PUTTING $j WATER OF $nwat ****"

			for {set i 0} {$i < 3} {incr i} {
				set newx($i) [expr { $xwat($i) + $shiftx($j) }]
				set newy($i) [expr { $ywat($i) + $shifty($j) }]
				set newz($i) [expr { $zwat($i) + $shiftz($j) }]
				set rn1 $l
				incr l
				set srn1 [string length $rn1]
				set sat1 [string length $at1($i)]

				if { $j > 99999 } {
					set j 1
				}

				set an1 $j
				set san1 [string length $an1]

				set x1 [format "%.3f" $newx($i)]
				set sx1 [string length $x1]

				set y1 [format "%.3f" $newy($i)]
				set sy1 [string length $y1]

				set z1 [format "%.3f" $newz($i)]
				set sz1 [string length $z1]
		
				set sresname 3

				if { $srn1 < 6 } {
					puts $w "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1($i)$sat($sat1)WAT$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           $ty($i)"
				} else {
					set srn1 5 
					puts $w "ATOM $p1($srn1)$rn1 $ic($sat1)$at1($i)$sat($sat1)WAT$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           $ty($i)"
				}
			}
			puts $w "TER"
			incr j
			incr nwatin
		}
	}

	set ul $j

	puts ""
	puts ""
	puts "			***** 50 % DONE *****"

	# LOWER WATER LAYER

	while { $nwatin < $nwat } {
		set shiftx($j) [expr { rand() * $box_x }]
		set shiftx1 [expr { rand() * $ext }]
		set shiftx($j) [expr { $shiftx($j) - $shiftx1 }]
		set shifty($j) [expr { rand() * $box_y }]
		set shifty1 [expr { rand() * $ext }]
		set shifty($j) [expr { $shifty($j) - $shifty1 }]
		set shiftz($j) [expr { rand() * ((2*$lz) + $box_z) * -1.0 }]
		set count 0
		for {set m 1} {$m < $j} {incr m} {
			set delx [expr { abs($shiftx($j) - $shiftx($m)) }]
			set delx2 [expr { $delx * $delx }]
			set dely [expr { abs($shifty($j) - $shifty($m)) }]
			set dely2 [expr { $dely * $dely }]
			set delz [expr { abs($shiftz($j) - $shiftz($m)) }]
			set delz2 [expr { $delz * $delz }]
			set del [expr { $delx2 + $dely2 + $delz2 }]
			if { $del < 0.5 } {
				incr count 
			}
		}

		if { [expr { abs($shiftz($j)) }] > [expr { (2*$lz) + $ext }] && $count == 0 } {
			puts "				**** PUTTING $j WATER OF $nwat ****"

			for {set i 0} {$i < 3} {incr i} {
				set newx($i) [expr { $xwat($i) + $shiftx($j) }]
				set newy($i) [expr { $ywat($i) + $shifty($j) }]
				set newz($i) [expr { $zwat($i) + $shiftz($j) }]
				set rn1 $l
				incr l
				set srn1 [string length $rn1]
				set sat1 [string length $at1($i)]

				if { $j > 99999 } {
					set j 1
				}

				set an1 $j
				set san1 [string length $an1]

				set x1 [format "%.3f" $newx($i)]
				set sx1 [string length $x1]

				set y1 [format "%.3f" $newy($i)]
				set sy1 [string length $y1]

				set z1 [format "%.3f" $newz($i)]
				set sz1 [string length $z1]
		
				set sresname 3

				if { $srn1 < 6 } {
					puts $w "ATOM  $p1($srn1)$rn1 $ic($sat1)$at1($i)$sat($sat1)WAT$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           $ty($i)"
				} else {
					set srn1 5 
					puts $w "ATOM $p1($srn1)$rn1 $ic($sat1)$at1($i)$sat($sat1)WAT$p($sresname)$p1($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           $ty($i)"
				}
			}
			puts $w "TER"
			incr j
			incr nwatin
		}
	}
	puts $w "END"
	close $w

	set n [open "nwat.pdb" "r"]
	set data3 [read $n]
	close $n

	set data4 "$data$data3"

	set fil [open "lip.pdb" "w"]
	puts $fil "$data4"
	close $fil

	# FORMING THE INPUT FOR THE EXECUTION OF LEAP

	set f [open "input_2" "w"]
	
	puts $f "\{ LEAP INPUT \}"
	puts ""
	puts ""
	puts "				#### ENTER THE NAME OF THE PARAMETER FILES (EACH SEPARATED BY A SPACE) ####"
	puts ""
	set inp1 [gets stdin]

	puts ""

	puts $f "\{ PARM: $inp1 \}"

	puts "				#### IS THERE ANY ADDITIONAL PARAMETER FILE? (frcmod file) (y/n) ####"
	puts ""
	set inp2 [gets stdin]

	if { $inp2 == "Y" || $inp2 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp3 [gets stdin]
	} else {
		set inp3 ""
	}

	puts $f "\{ frcmod: $inp3 \}"

	puts "				#### IS THERE ANY .LIB OR .OFF STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp4 [gets stdin]

	if { $inp4 == "Y" || $inp4 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp5 [gets stdin]
	} else {
		set inp5 ""
	}

	puts $f "\{ lib: $inp5 \}"

	puts "				#### IS THERE ANY .prepin STRUCTURE FILES YOU WANT TO ADD? (Y/N) ###"
	puts ""

	set inp6 [gets stdin]

	if { $inp6 == "Y" || $inp6 == "y" } {
		puts "			#### ENTER THE NAME OF THE FILE (EACH SEPARATED BY SPACE) ####"
		puts ""
		set inp7 [gets stdin]
	} else {
		set inp7 ""
	}

	puts $f "\{ prep: $inp7 \}"

	puts $f "\{ water_layer: NA \}"

	close $f

	# READING THE INPUT FILE

	set in [open "input_2" "r"]
	set inp [read $in]
	close $in

	# FORMIING THE LEAP FILE

	set l [open "leap.in" "w"]

	# PARAMETER FILE
	set k 1
	while { $k < [llength [lindex $inp 1]] } {
		puts $l "source leaprc.[lindex $inp 1 $k]"
		incr k
	}  

	# FRCMOD FILE
	set k 1
	while { $k < [llength [lindex $inp 2]] } {
		puts $l "loadamberparams [lindex $inp 2 $k]"
		incr k
	}

	# LIB FILE
	set k 1
	while { $k < [llength [lindex $inp 3]] } {
		puts $l "loadoff [lindex $inp 3 $k]"
		incr k
	}

	# PREPIN FILE
	set k 1
	while { $k < [llength [lindex $inp 4]] } {
		puts $l "loadamberprep [lindex $inp 4 $k]"
		incr k
	}
		
	# PDB FILE
	puts $l "mol=loadpdb \"lip.pdb\""

	# BOX_DIMENSIONS
	puts "			### GETTING THE BOX DIMENSIONS ###"

	box_dimensions "lip.pdb"
		
	set d [open "box_dimensions" "r"]
	set bd [read $d]
	close $d
	puts $l "set mol box { [lindex $bd 0] [lindex $bd 1] [lindex $bd 2] }"

	# IONS
	puts ""
	puts "				#### DO YOU WANT TO MANUALLY ENTER THE NUMBER OF IONS (Y/N) ####"
	puts ""
	set ionans [gets stdin]
	if { $ionans == "y" || $ionans == "Y" } {
		puts ""
		puts "			#### ENTER THE NUMBER OF CHLORIDE IONS ####"
		set cl [gets stdin]
		puts ""
		puts "			#### ENTER THE NUMBER OF POTASSIUM IONS ####"
		set K [gets stdin]
		puts ""
		puts $l "addions mol Cl- $cl"
		puts $l "addions mol K+ $K"
	} else {
		puts $l "addions mol Cl- 0.0"
		puts $l "addions mol K+ 0.0"
	}

	puts $l "saveamberparm mol mol.prmtop mol.inpcrd"
	puts $l "savepdb mol reff.pdb"
	puts $l "quit"

	close $l

	exec tleap -f leap.in
	puts ""
	puts "				**** DONE ****"
	puts ""
}

proc file_overlap {} {

	set f [open "water_resid" "r"]
	set data [read $f]
	close $f

	set g [open "wat_resid_sorted" "w"]

	set k 0

	set i 0
	while { $k < [llength $data] } {

		puts "$k of [llength $data]"
		set elem($i) [lindex $data $k]
	
		set count 0
		for {set j 0} {$j  < $i} {incr j} {
			if { $elem($j) == $elem($i) } {
				incr count
			}
		}
		if { $count == 0 } {
			puts $g "$elem($i),\\"
			incr i
		}
		incr k
	}
	close $g

	set f [open "water_resid1" "r"]
	set data [read $f]
	close $f

	set g [open "wat_resid_sorted1" "w"]

	set k 0

	set i 0
	while { $k < [llength $data] } {

		puts "$k of [llength $data]"
		set elem($i) [lindex $data $k]
	
		set count 0
		for {set j 0} {$j  < $i} {incr j} {
			if { $elem($j) == $elem($i) } {
				incr count
			}
		}
		if { $count == 0 } {
			puts $g "$elem($i),\\"
			incr i
		}
		incr k
	}
	close $g

	set g1 [open "wat_resid_sorted" "r"]
	set data [read $g1]
	close $g1

	set g2 [open "wat_resid_sorted1" "r"]
	set data1 [read $g2]
	close $g2

	set data2 $data$data1

	set g3 [open "1wat_resid_sorted" "w"]

	puts $g3 "$data2"

	close $g3
}

proc mod1 { inp1 inp2 inp3 inp4 inp5 } {
	# THIS PROCEDURE READS THE INPUT PDB FILE AND REMOVES THE CHAIN ID IF ITS MISSING AND ADD THE OCCUPANCY (1.00) AND TEMPERATURE FACTOR (0.00)AND SHIFTS THE ORIGIN TO THE PIVOT ATOM OF RESIDUE 1

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


	set f [open "$inp1" "r"]
	set data1 [read $f]
	close $f

	set k 1
	while { [lindex $data1 $k] != "ATOM" && [lindex $data1 $k] != "HETATM" } {
		incr k
	}

	set occ [lindex $data1 [expr { $k - 3 }]]
	set temp [lindex $data1 [expr { $k - 2 }]]

	set t 0
	while { [string range $occ $t $t] != "." } {
		incr t
	}
	set docc [string range $occ [expr { $t + 1 }] end]

	set t 0
	while { [string range $temp $t $t] != "." } {
		incr t
	}
	set dtemp [string range $temp [expr { $t + 1 }] end]

	if { [string length $docc] > 2 && [string length $dtemp] > 2 } {
		puts "			# PUTTING IN THE TEMPERATURE FACTOR VALUE AS 0.00 AND OCCUPANCY VALUE AS 1.00 IN THE PDB"
		set oshift1 2
	} elseif { [string length $docc] > 2 && [string length $dtemp] < 2 } { 
		puts "			# PUTTING IN THE OCCUPANCY VALUE AS 1.00 IN THE PDB"
		set oshift1 1
	} elseif { [string length $docc] < 2 && [string length $dtemp] > 2 } {
		puts "			# PUTTING IN THE TEMPERATURE FACTOR VALUE AS 0.00 IN THE PDB"
		set oshift1 1
	} elseif { [string length $docc] <= 2 && [string length $dtemp] <= 2 } {
		puts "			# NO PROBLEM WITH THE TEMPERATURE FACTOR VALUE AND OCCUPANCY VALUE IN THE PDB"
		set oshift1 0
	}


	set test [lindex $data1 5]
	set st [string length $test]

	if { $st < 5 } { 
		set oshift 0
		puts "			# REMOVING CHAIN ID"
	} else { 
		set oshift 1
		puts "				# NO CHAIN ID FOUND"
	} 


	set h [open "$inp5" "w"]
	set xori $inp2
	set yori $inp3
	set zori $inp3

	set k 0
	set rnum 1
	while { $k < [llength $data1] } {
		set term [lindex $data1 $k]
		set t1 [string range $term 0 5]
		if { $term == "TER" } {
			puts $h "TER"
			set rnum 1
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
			set rn1 $rnum
			incr rnum
			set srn1 [string length $rn1]

			set ft "HETATM"

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

			if { $oshift == 0 } {
				set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
				set schain_id [string length $chain_id]
			} else { 
				set chain_id " "
				set schain_id 1
			}

			if { $schain_id > 1 } {
				set shift1 1
				set an1 [lindex $data1 [expr { $k + 4 - $shift - $shift2 }]]
				set san1 4
				set cn ""
			} else { 
				#set cn [lindex $data1 [expr { $k + 4 - $shift - $shift2  }]]
				set cn " "
				set an1 [lindex $data1 [expr { $k + 5 -$shift - $shift2 - $oshift}]]
				set san1 [string length $an1]
				set shift1 0
			}


			set x1 [lindex $data1 [expr { $k + 6 - $shift - $shift1 -$shift2 - $oshift}]]
			set sx1 [string length $x1]
			if { $sx1 > 8 } {
				set shift3 1
				set t 0
				while { [string range $x1 $t $t] != "." } {
					incr t
				}
				set corx [string range $x1 0 [expr { $t + 3 }]]
				set corx [expr { $corx - $xori }]
				set corx [format "%.3f" $corx]
				set cory [string range $x1 [expr { $t + 4 }] end]
				set cory [expr { $cory - $yori }]
				set cory [format "%.3f" $cory]
				set x1 $corx$cory
				set sx1 [string length $corx]
				set y1 ""
				set sy1 8
				set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $oshift}]] 
				set z1 [format "%.3f" [expr { $z1 - $zori }]]
				set sz1 [string length $z1]
			} else { 
				set shift3 0
				set x1 [format "%.3f" [expr { $x1 - $xori }]]
				set sx1 [string length $x1]
				set y1 [lindex $data1 [expr { $k + 7 - $shift -$shift1 - $shift2 - $shift3 - $oshift}]] 
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set shift4 1
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set cory [string range $y1 0 [expr { $t + 3 }]]
					set cory [expr { $cory - $yori }]
					set cory [format "%.3f" $cory]
					set corz [string range $y1 [expr { $t + 4 }] end]
					set corz [expr { $corz - $zori }]
					set corz [format "%.3f" $corz]					
					set sy1 [string length $cory]
					set y1 $cory$corz
					set z1 ""
					set sz1 8
				} else {
					set shift4 0
					set y1 [format "%.3f" [expr { $y1 - $yori }]]
					set sy1 [string length $y1]
					set z1 [lindex $data1 [expr { $k + 8 - $shift -$shift1 - $shift2 - $shift3 - $shift4 - $oshift}]] 
					set z1 [format "%.3f" [expr { $z1 - $zori }]]
					set sz1 [string length $z1]
				}
			}
			if { [string length [lindex $data1 [expr { $k + 9 - $shift -$shift1 - $shift2 - $shift3 - $shift4 - $oshift }]]] > 5 } {
				set shift5 1
			} else {
				set shift5 0
			}
		
			puts $h "$ft$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$p($sresn)$resn $cn$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data1 [expr { $k + 11 - $shift - $shift1 - $shift2 - $shift3 - $shift4 - $shift5 - $oshift - $oshift1 }]]"
		}
		incr k
	}
	puts $h "END"
	close $h
}

proc sol_nc_inp {} {
	set f [open "water_input" "w"]

	puts "				#### ENTER THE NUMBER OF WATERS YOU WANT TO ADD ####"
	puts ""
	
	set inp1 [gets stdin]

	puts $f "\{ NUMBER_OF_WATER_MOLECULES: $inp1 \}"
	
	puts "				#### ENTER THE LAYER THICKNESS IN ANGSTOMS ####"
	puts ""
		
	set inp2 [gets stdin]

	puts $f "\{ LAYER_SIZE: $inp2 \}"

	puts $f "\{ EXTENSION: 5.0 \}"

	close $f

	set g [open "wat.pdb" "w"]

	set water {ATOM      1  O   WAT     1       0.000   0.000   0.000  1.00  0.00
ATOM      2  H1  WAT     1       0.957   0.000   0.000  1.00  0.00
ATOM      3  H2  WAT     1      -0.240   0.927   0.000  1.00  0.00
TER   
END}

	puts $g "$water"

	close $g   
}
proc delete {} {
	file delete dummy
	file delete dummy1
	file delete grid
	#file delete input
	file delete input_2
	file delete water_input
	file delete wat.pdb
	file delete water.pdb
	file delete lipids.pdb
	file delete nwat.pdb
	file delete lipid_order.pdb
	file delete box_dimensions
	file delete leap.in
	file delete leap_1.log
	file delete leap_1.in
	file delete mod_water.pdb
	file delete mod_water_vis.pdb
	file delete nd1_d2.pdb
	file delete nd1_d2_vis.pdb
	file delete neighbours
	file delete only_pro.pdb
	file delete refm1.pdb
	#file delete water_resid
	file delete mod_lipids.pdb
	#file delete only_mod_protein.pdb
	file delete lip_pro_final.pdb
	#file delete lip_pro.pdb
	file delete overlap_atom
	file delete spherical_coord
	file delete spherical_coord_ll
	file delete water_resid
	file delete water_resid1
	file delete wat_resid_sorted
	file delete wat_resid_sorted1
	file delete water_pore.pdb
	file delete sol_water.pdb
	file delete spherical_coord_wp
	file delete lip_pro.pdb
	file delete only_mod_protein.pdb
}

proc file_description {} {
	puts ""
	puts "******************************************************************************************"	
	puts "				#### DESCRIPTION OF THE FILES ####"
	puts "				(NOTE: SOME FILES MAY NOT BE GENERATED DEPENDING UPON THE OPTION YOU HAVE CHOSEN)"
	puts ""

	puts "lipids_no.pdb : UNSOLVATED PDB OF PURE LIPID BILAYER"
	puts "1_final_struct.pdb : UNSOLVATED PDB OF LIPID BILAYER WITH INSERTED COMPONENT"
	puts "lip.pdb: SOLVATED PURE LIPID BILAYER"
	puts "reff.pdb: FINAL SOLVATED PDB OF THE SYSTEM"	
	puts "mol.prmtop: FINAL AMBER PRMTOP FILE"
	puts "mol.inpcrd: FINAL AMBER COORDINATE FILE"
	puts "dir1: PROFILE OF THE PROTEIN INSIDE THE BILAYER"
	puts "lipids_wo.pdb: FINAL DRY VESSICLE OR DIBs"
	puts "sol_lipid.pdb: FINAL SOLVATED VESICLE OR DIBs"
	puts "*******************************************************************************************"	
}

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
	
	puts ""
	puts "			#### ENTER THE NAME OF THE PDB ####"
	puts "			    (MAKE SURE TO REMOVE THE REMARK SECTION OF THE PDB)"
	set inppdb [gets stdin]
	exec ls $inppdb

	puts ""
	puts "				#### DO YOU WANT TO MANUALLY INSERT THE PROTEIN BASED ON THE HYDROPHOBIC AND HYDROPHYLIC AMINO ACIDS? (Y/N)"
	puts ""

	set inp [gets stdin]

	if { $inp == "y" || $inp == "Y" } {
		
		puts "			***** NOTE : BASED ON THE INPUTS YOU GIVE BELOW THE INPUT PROTEIN PDB WILL BE TRANSLATED AND ROTATED ACCORDING TO THESE VALUES *****"
		puts ""

		set inp1 $tx
		set inp2 $ty

		puts "			#### THE CODE ALLOWS MAXIMUM OF 10 ATTEMPTS TO INSERT THE PROTEIN, CHECK 'pro_ori.pdb' AFTER EACH ATTEMPTS TO CHECK THE CORRECTNESS OF THE ORIENTATION ####"
		puts ""

		for {set attempt 0} {$attempt < 10} {incr attempt} {
			puts "			**** ATTEMPT $attempt *****"
			puts ""
			puts "			#### ENTER THE VALUE FOR Z TRANSLATION (HINT: MOST IMPORTANT VALUE AS IT CONTROLS THE MOVEMENT OF THE PROTEIN ALONG THE BILAYER THICKNESS) ####"
			puts "			                                 START WITH A WILD GUESS, PLOT 'ac_nature' FILE'S 2ND AND 3RD COLUMN TO GET THE BEST GUESS"
			set inp3 [gets stdin]

			puts "			#### ENTER THE VALUE FOR X ROTATION (NOTE: CONTROLS THE TILT OF THE MEMBRANE PROTEIN) ####"
			puts ""
			set inp4 [gets stdin]

			puts "			#### ENTER THE VALUE FOR Y ROTATION (NOTE: CONTROLS THE TILT OF THE MEMBRANE PROTEIN) ####"
			puts ""
			set inp5 [gets stdin]

			puts "			#### ENTER THE VALUE FOR Z ROTATION (NOTE: PROTEIN IS INVARIANT UNDER Z ROTATION, THIS VALUE ONLY CONTROLS THE SIMULATION SETUP) ####"
			puts ""
			set inp6 [gets stdin]
	
			temp_grid
			pi 0 $inppdb $inp1 $inp2 $inp3 $inp4 $inp5 $inp6
			protein_profile_pi ref1_no.pdb 2

			puts "				#### IMPORTANT - CHECK FILE pro_ori.pdb TO CHECK IF THE ORIENTATION IS CORRECT ####"
			puts "				    PLOTING SECOND (HYDROPHOBIC) AND THIRD (HYDROPHILIC) COLUMN AGAINST FIRST COLUMN IS ANOTHER WAY TO CHECK THE CORRECTNESS"
			puts ""		
			puts "				#### ARE YOU SATISFIED WITH THE ORIENTATION??? (Y/N)####"
			puts ""
			set ussat [gets stdin]
			if { $ussat == "Y" || $ussat == "y" } {
				set attempt 10
			}
		}
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
		protein_profile_pi ref1_no.pdb 2
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
		set ztrans [expr { $zori + $ztrans }]
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
				set tcoord [list $corx $cory $z1] 
		
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
				set tcoord [list $x1 $cory $corz] 
		
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
	puts "			*** THE SPREAD OF THE PROTEIN ALONG X AXIS IS [expr { $maxx - $minx }] ***"
	puts "			*** THE SPREAD OF THE PROTEIN ALONG Y AXIS IS [expr { $maxy - $miny }] ***"
	puts "			*** THE SPREAD OF THE PROTEIN ALONG Z AXIS IS [expr { $maxz - $minz }] ***"
	set g [open "dummy" "w"]
	puts $g " { $minx $maxx } "
	puts $g " { $miny $maxy } "
	puts $g " { $minz $maxz } "
	puts $g " { [expr { $maxx - $minx }] [expr { $maxy - $miny }] [expr { $maxz - $minz }] }"
	close $g
}

proc protein_profile_pi { inp dir} {

	protein_spread_pi $inp
	
	puts ""
	puts "				**** DETERIMING THE PROTEIN PROFILE ****"
	puts ""

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
	puts "						(x,y,z) = (0,1,2) "
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

		puts "				**** $ij ALONG DIRECTION $dir ****"
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

proc box_dimensions {PDB} {

	# THIS PROCEDURE WILL CALCULATE THE FINAL BOX DIMENSIONS SIMILAR TO THE VMD.sh SCRIPT

	set f [open "$PDB" "r"]
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
	set g [open "box_dimensions" "w"]
	puts $g " { [expr { ($maxx - $minx) + 2.0 }] [expr { ($maxy - $miny) + 2.0 }] [expr { ($maxz - $minz) + 2.0 }] }"
	close $g
}
proc file_delete {} {
	file delete dummy
	file delete temp_grid.pdb
	file delete res_num
}

proc files {} {
	puts ""
	puts "				#### DESCRIPTION OF THE VARIOUS FILES FORMED ####"
	puts ""

	puts "	ac_nature : DISTRICUTION OF VARIOUS TYPES OF AMINO ACIDS ALONG THE Z DIRECTION"
	puts "	amino_acid_dis : DISTRICUTION INDIVIDUAL AMINO ACIDS ALONG THE Z DIRECTION"
	puts "	dir1: SPREAD OF THE INSERTED PROTEIN ALONF THE XY PLANE"
	puts ""
	puts ""
	puts "	ref1_no.pdb: THE FINAL PROTEIN PDB WHICH GOES INTO THE MEMBRANE BUILDER"
	puts "	pro_ori.pdb: ORIENTATION OF THE ABOVE PDB WITH RESPECT TO REFERENCE BILAYER OF THICKNESS 38.0 Ang"
}

# INCLUSION OF LEARNING ALL ATOM GRID OF AMBAT INTO THE MAIN CODE

#####################################################################################################################
proc non_random_grid {} {
	# THIS FUCNTION HELPS BUILD NON RANDOM GRID BASED ON THE RULES DEFINED IN FILE 'RULE'

	# SCORING FUNCTIONS: conditions - aj+1/aj=b and (ai*Ri/ni)+(ai+1*Ri+1/ni+1)+(ai+2*Ri+2/ni+2) ...... (anRn/nn)=1

	# GLOBAL VARIABLES

	global lip1_rule
	global lip2_rule
	global lip_ul
	global lip_ll
	global lip_grid
	global rule
	global grid

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "	
	set p(4) ""
	
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

	set sat(1) "   "
	set sat(2) "  "
	set sat(3) " "
	set sat(4) " "

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""

	set f [open "input" "r"]
	set data [read $f]
	close $f

	# DIFFERENT LIPID TYPES

	set i 0
	set k 1
	while { $k < [llength [lindex $data 1]] } {
		set lip_type($i) [lindex $data 1 $k]
		set atom_type($i) [lindex $data 2 $k]
		incr i
		incr k 
	}
	set num_lip_type $i

	# COUNTING THE NUMBER OF LIPIDS AND THE LIPID PIVOT ATOM IN EACH LEAFLET

	set k 1
	set num_lip_ul 0
	set i 0
	while { $k < [llength [lindex $data 3]] } {
		set num_lip_ul [expr { $num_lip_ul + [lindex $data 3 $k] }]
		set lip_ul($i) [lindex $data 3 $k]
		incr i
		incr k
	}

	set k 1
	set num_lip_ll 0
	set i 0
	while { $k < [llength [lindex $data 4]] } {
		set num_lip_ll [expr { $num_lip_ll + [lindex $data 4 $k] }]
		set lip_ll($i) [lindex $data 4 $k]
		incr i
		incr k
	}

	# BOX INFORMATION

	set X [lindex $data 5 1]
	set Y [lindex $data 5 2]
	set Z [lindex $data 5 3]

	if { $num_lip_ul > $num_lip_ll } {
		set num_lip $num_lip_ul
	} else {
		set num_lip $num_lip_ll
	}

	set gridx [expr { $X / sqrt($num_lip) }]
	set gridy [expr { $Y / sqrt($num_lip) }]

	# DEFINATION OF RULES

	puts "				ENTER THE RULES AS SIMPLE EQUAL AND NON-EQUAL STATEMENTS"
	puts "				LIKE :::: \$lipi==\/!=\$lip2 \$grid"
	puts ""

	puts "				#### ENTER THE NUMBERS OF RULES YOU WANT TO DEFINE ####"
	set nrules [gets stdin]
	puts ""

	set k 0
	set i 0
	while { $k < $nrules } {
		puts "				#### ENTER RULE [expr { $k + 1 }] ####"
		set dum [gets stdin]
		
		set t 0
		while { [string range $dum $t $t] != "=" && [string range $dum $t $t] != "!" } {
			incr t
		}

		set lip1_rule($i) [string range $dum 0 [expr { $t - 1 }]]
		set rule($i) [string range $dum $t [expr { $t + 1 }]]
		set lip2_rule($i) [string range $dum [expr { $t + 2 }] [expr { [string length $dum] - 3 }]]
		set grid($i) [string range $dum [expr { [string length $dum] - 1 }] [string length $dum]]
		incr i
		incr k
	}

	# NEIGHBOURS DEFINATION (CUMMULATIVE)

	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }]} {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			set nei($i,$j) ""
		}
	}

	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }]} {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			for {set k1 0.00} {$k1 <= $i} {set k1 [expr { $k1 + $gridx }]} {
				for {set k2 0.00} {$k2 <= $j} {set k2 [expr { $k2 + $gridy }]} {
					if { $k1 != $i || $k2 != $j } {
						set xdis [expr { abs($i-$k1) }]
						set ydis [expr { abs($j-$k2) }]
						if { $xdis > [expr { $X / 2.0 }] } {
							set xdis [expr { $X - $xdis }]
						}
						if { $ydis > [expr { $Y / 2.0 }] } {
							set ydis [expr { $Y - $ydis }]
						}
						set xdis [format "%.3f" $xdis]
						set ydis [format "%.3f" $ydis]
						if { [expr { $xdis / $gridx }] <= 1 && [expr { $ydis /$gridy }] <= 1 } {
							lappend nei($i,$j) $k1
							lappend nei($i,$j) $k2
						}
					}
				}
			}
		}
	}

	# TESTING NEIGHBOUR DEFINATION

	#set net [open "nei_test" "w"]

	#for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }]} {
		#for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			#puts $net "FOR $i $j"
			#puts $net "$nei($i,$j)"
		#}
	#}
	#close $net

	# FORMING THE GRID :: UPPER LEAFLET

	#set te1 [open "test_UL" "w"]

	set cum_score 0.0
	set k 0
	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }] } {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			#puts "			### GRID $i $j ###"
			set lip1 [expr { rand() * [expr { $num_lip_type - 1 }] }]
			set lip1 [format "%.0f" $lip1]
			set type_lip(0) $lip1
			set lip(0) $lip_type($lip1)
			set k1 1
			for {set k 0} {$k < $num_lip_type} {incr k} {
				if {$k != $lip1 } {
					set type_lip($k1) $k
					set lip($k1) $lip_type($k)
					incr k1
				}
			}
			if { $lip_ul($type_lip(0)) != 0 } {

				# CHECKING IF ANY RULE DEFINED ON lip1

				set count 1
				#for {set k 0} {$k < $nrules} {incr k} {
					#if { $lip(0) == $lip1_rule($k) || $lip(0) == $lip2_rule($k) } {
						#incr count
					#}
				#}

				# SCORING FUNCTION

				if { $count != 0 } {
					set score(0) [scoring $nrules $nei($i,$j) $lip(0)]
					for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
						set score($k1) [scoring $nrules $nei($i,$j) $lip($k1)]
					}
				} else {
					set score(0) 1.0
					for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
						set score($k1) 0.0
					}
				}
			} else {
				set score(0) -1
				for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
					set score($k1) [scoring $nrules $nei($i,$j) $lip($k1)]
				}
			}	

			# CHOOSING THE LIPID FOR THIS GRID

			set gridlip 0
			set final_score $score(0)
			for {set k 1} {$k < $num_lip_type} {incr k} {
				if { $score($k) > $final_score && $lip_ul($type_lip($k)) > 0 } {
					set gridlip $k
					set final_score $score($k)
				}
			}
			set cum_score [expr { $final_score + $cum_score }]
			set lip_grid($i,$j) $lip($gridlip)
			set lip_grid_type($i,$j) $atom_type($type_lip($gridlip))
			#puts "			$lip_grid($i,$j) $lip_grid_type($i,$j)"
			set lip_ul($type_lip($gridlip)) [expr { $lip_ul($type_lip($gridlip)) - 1 }]		
		
			#puts $te1 "$lip(0) $lip(1) $lip(2) $score(0) $score(1) $score(2)"
			#puts $te1 "chosen : $gridlip"
		
		}
	}	
	#close $te1
	puts "		#### CUMULATIVE SCORE UPPER LEAFLET = $cum_score #### "

	# FORMING A PDB FILE

	set g [open "upper_leaflet.pdb" "w"]

	set k 1
	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }]} {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			set rn1 $k
			set srn1 [string length $rn1]
			
			set an1 $k
			set san1 [string length $an1]
			incr k

			set at1 $lip_grid_type($i,$j) 
			set sat1 [string length $at1]

			set resname $lip_grid($i,$j)
			set sresname [string length $resname]

			set x1 [format "%.3f" $i]
			set sx1 [string length $x1]

			set y1 [format "%.3f" $j]
			set sy1 [string length $y1]

			set z1 [format "%.3f" 0.000]
			set sz1 [string length $z1]

			puts $g "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname) $p($san1)$an1     $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
			puts $g "TER"
		}
	}
	close $g
	set cumk $k

	# FORMING THE GRID :: LOWER LEAFLET

	#set te [open "test_LL" "w"]

	set cum_score 0.0
	set k 0
	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }] } {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			#puts "			### GRID $i $j ###"
			set lip1 [expr { rand() * [expr { $num_lip_type - 1 }] }]
			set lip1 [format "%.0f" $lip1]
			set type_lip(0) $lip1
			set lip(0) $lip_type($lip1)
			set k1 1
			for {set k 0} {$k < $num_lip_type} {incr k} {
				if {$k != $lip1 } {
					set type_lip($k1) $k
					set lip($k1) $lip_type($k)
					incr k1
				}
			}
			if { $lip_ll($type_lip(0)) != 0 } {

				# CHECKING IF ANY RULE DEFINED ON lip1

				set count 1
				#for {set k 0} {$k < $nrules} {incr k} {
					#if { $lip(0) == $lip1_rule($k) || $lip(0) == $lip2_rule($k) } {
						#incr count
					#}
				#}

				# SCORING FUNCTION

				if { $count != 0 } {
					set score(0) [scoring $nrules $nei($i,$j) $lip(0)]
					for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
						set score($k1) [scoring $nrules $nei($i,$j) $lip($k1)]
					}
				} else {
					set score(0) 1.0
					for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
						set score($k1) 0.0
					}
				}
			} else {
				set score(0) -1
				for {set k1 1} {$k1 < $num_lip_type} {incr k1} {
					set score($k1) [scoring $nrules $nei($i,$j) $lip($k1)]
				}
			}	

			# CHOOSING THE LIPID FOR THIS GRID

			set gridlip 0
			set final_score $score(0)
			for {set k 1} {$k < $num_lip_type} {incr k} {
				if { $score($k) > $final_score && $lip_ll($type_lip($k)) > 0 } {
					set gridlip $k
					set final_score $score($k)
				}
			}
			set cum_score [expr { $final_score + $cum_score }]
			set lip_grid($i,$j) $lip($gridlip)
			set lip_grid_type($i,$j) $atom_type($type_lip($gridlip))
			#puts "			$lip_grid($i,$j) $lip_grid_type($i,$j)"
			set lip_ll($type_lip($gridlip)) [expr { $lip_ll($type_lip($gridlip)) - 1 }]		
		
			#puts $te "$lip(0) $lip(1) $lip(2) $score(0) $score(1) $score(2)"
			#puts $te "chosen : $gridlip"
		}
	}	
	#close $te
	puts "		#### CUMULATIVE SCORE LOWER LEAFLET $cum_score #### "

	# FORMING A PDB FILE

	set g1 [open "lower_leaflet.pdb" "w"]

	set k $cumk
	for {set i 0.00} {$i < $X} {set i [expr { $i + $gridx }]} {
		for {set j 0.00} {$j < $Y} {set j [expr { $j + $gridy }]} {
			set rn1 $k
			set srn1 [string length $rn1]
			
			set an1 $k
			set san1 [string length $an1]
			incr k

			set at1 $lip_grid_type($i,$j) 
			set sat1 [string length $at1]

			set resname $lip_grid($i,$j)
			set sresname [string length $resname]

			set x1 [format "%.3f" $i]
			set sx1 [string length $x1]

			set y1 [format "%.3f" $j]
			set sy1 [string length $y1]

			set z1 [format "%.3f" [expr { -1 * $Z }]]
			set sz1 [string length $z1]

			puts $g1 "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$sat($sat1)$resname$p($sresname) $p($san1)$an1     $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           O"
			puts $g1 "TER"
		}
	}
	puts $g1 "END"
	close $g1

	set g [open "upper_leaflet.pdb" "r"]
	set data1 [read $g]
	close $g

	set g1 [open "lower_leaflet.pdb" "r"]
	set data2 [read $g1]
	close $g1

	set h [open "lipids.pdb" "w"]
	puts $h "$data1$data2"
	close $h

	file delete lower_leaflet.pdb
	file delete upper_leaflet.pdb

	set method [lindex $data 9 1]
	set version [lindex $data 13 1]
	set num_dl [lindex $data 13 [expr { [llength [lindex $data 13]] - 1 } ]]
	set lip_name [lindex $data 13 2]
	set gd [lindex $data 9 2]

	set k 2

	if { $method == 1 } {
		if { $version == 3.0 } {
			puts "				#### USING VERSION 3.0 ####"
			if { $num_dl > 1 } {
				grid_mapping_assy
			} else {
				grid_mapping "[lindex $data1 13 $k]"
			}		
		} else {
			puts "				#### USING VERSION 2.0 ####" 
			lipid_growth
		}
	} else { 
		lipid_overlap
	}
}

proc scoring {nrules nei lip1} {

	global lip1_rule
	global lip2_rule
	global lip_ul
	global lip_ll
	global lip_grid
	global rule
	global grid

	# DETERMING THE COEFFICIENT FOR THE SCORING FUNCTION

	set beta 0.95
	set value 1.0
	set term 1.0

	for {set i 1} {$i < $nrules} {incr i} {
		set value [expr { $value*$beta }]
		set term [expr { $term + $value }]
	}
	set alpha(0) [expr { 1 / $term }]
	
	set value 1
	for {set i 1} {$i < $nrules} {incr i} {
		set value [expr { $value*$beta }]
		set alpha($i) [expr { $alpha(0)*$value }]
	}

	for {set k 0} {$k < $nrules} {incr k} {
		set broken($k) 0
	}
	set k 0
	while { $k < [llength $nei] } {
		set t1 [lindex $nei $k]
		set t2 [lindex $nei [expr { $k + 1 }]]
		set lip2 $lip_grid($t1,$t2)
		for {set k1 0} {$k1 < $nrules} {incr k1} {
			#puts "$lip2 $lip1 $lip1_rule($k1) $lip2_rule($k1) $rule($k1)"
			if { $rule($k1) == "==" } {
				if { $lip2 == $lip1_rule($k1) && $lip1 != $lip2_rule($k1) } {
					incr broken($k1)
				}	elseif { $lip2 == $lip2_rule($k1) && $lip1 != $lip1_rule($k1) } {
					incr broken($k1)
				} elseif { $lip1 == $lip1_rule($k1) && $lip2 != $lip2_rule($k1) } {
					incr broken($k1)
				}	elseif { $lip1 == $lip2_rule($k1) && $lip2 != $lip1_rule($k1) } {
					incr broken($k1)
				} 
			} elseif { $rule($k1) == "!=" } {
				if { $lip2 == $lip1_rule($k1) && $lip1 == $lip2_rule($k1) } {
					incr broken($k1)
				}	elseif { $lip2 == $lip2_rule($k1) && $lip1 == $lip1_rule($k1) } {
					incr broken($k1)
				}
			}
		}					
		incr k 2
	}
			
	set score 0
	set nvoilation 0
	for {set k 0} {$k < $nrules} {incr k} {
		set score [expr { $score + ($alpha($k) / (1+$broken($k))) }]
	}
	return $score
}
#####################################################################################################################

mem_builder
delete
file_description
