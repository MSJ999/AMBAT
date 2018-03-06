proc execution {} {
	# ******************************************************************************************************************************************************************************#

	puts "	************************************************"
	puts "	THIS IS A AMBER BASED MEMBRANE ANALYSIS TOOL"
	puts "	************************************************"
	puts ""
	puts ""
	puts "	DEVELOPED BY TARUN KHANNA, Dr. IAN GOULD"
	puts "	IMPERIAL COLLEGE LONDON, U.K."	
	puts ""
	puts ""
	puts "	FOR DOCUMENTATION OF THE CODE SEE :"
	puts ""
	puts ""
	puts "	THE CODE ALLOWS THE CALCULATION OF LOCAL MEMBRANE PROPERTIES,CONSTITUTING LIPID PROFILE AND CONSTITUTING PROTEIN PROFILE"
	puts ""
	puts "	IMP NOTE : THE CODE REQUIRES THE BASIC FUNCTIONALITY OF CPPTRAJ TO READ NETCDF FILE FOR ITS EXECUTION"
	puts "	CPPTRAJ CAN BE DOWNLOAD AS A PART OF AMBER TOOLS FROM http://ambermd.org/#AmberTools"
	puts ""

	# ******************************************************************************************************************************************************************************#

  # CODE EXECUTION OPTIONS

	puts ""
	puts "		### THE CODE CALCULATES THE FOLLOWING PROPERTIES ###"
	puts ""
	puts "	- LOCAL MEMBRANE PROPERTIES"
	puts "		- LOCAL AREA PER LIPID"
	puts "		- LOCAL BILAYER THICKNESS"
	puts "		- LOCAL MEMBRANE CURVATURE"
	puts "		- ION/WATER FLUX"
	puts ""
	puts "	- PROTEIN PROFILE"
	puts ""
	puts "	- LIPID PROFILE"
	puts "	- WATER PROFILE"
	puts ""

	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE LOCAL AREA PER LIPID AND LOCAL BILAYER THICKNESS? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		set MRES 1
	} else {
		set MRES 0
	}
	
	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE LOCAL MEMBRANE CURVATURE? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		set CRES 1
	} else {
		set CRES 0
	}

	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE ION FLUX THROUGH THE CHANNEL? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		puts ""
		puts "	### ENTER THE PDB NAME OF THE ION ###"
		set ionn [gets stdin]
		set IRES 1
	} else {
		set IRES 0
	}

	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE WATER FLUX THROUGH THE CHANNEL? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		set WRES 1
	} else {
		set WRES 0
	}

	puts ""
	puts "		### DO YOU WANT TO CALCULATE PROTEIN PROFILE? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		puts ""
		puts "		### ENTER THE NUMBER OF SLICES YOU WANT TO CREATE ###"
		puts ""
		puts "HINT : zlice == 1 AVERAGE OVER THE WHOLE BILAYER AND RESULTS IN 2D PROFILE (R,THETA)"
		puts "     AND zlice > 1 RESULTS IN 3D PROFILE WITH EACH SLICE AS THE ADDITIONAL DIMENSION"
		puts ""
		set zslice [gets stdin]
		set PRES 1
	} else {
		set PRES 0
	}
	
	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE LIPID PROFILE? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		set LRES 1
	} else {
		set LRES 0
	}

	puts ""
	puts "		### DO YOU WANT TO CALCULATE THE WATER PROFILE? (Y/N) ###"
	set res [gets stdin]

	if { $res == "y" || $res == "Y" } {
		puts ""
		puts "		### ENTER THE NUMBER OF SLICES YOU WANT TO CREATE ###"
		puts ""
		set zwpslice [gets stdin]
		set WPRES 1
	} else {
		set WPRES 0
	}

	# ******************************************************************************************************************************************************************************#

	# INPUT OPTIONS

	puts "		### ENTER THE NAME OF THE PRMTOP FILE ###"
	set prmtop [gets stdin]
	exec ls $prmtop
	
	puts ""
	puts "		### ENTER THE NAME OF THE COORDINATE FILE ###"
	set crd [gets stdin]
	exec ls $crd

	puts ""
	puts "		### ENTER THE STARTING FRAME ###"
	set start_frame [gets stdin]

	puts ""
	puts "		### ENTER THE END FRAME ###"
	set end_frame [gets stdin]

	puts ""
	puts "		### ENTER THE STEP SIZE BETWEEN THE TWO CONSECUTIVE FRAMES ###"
	set step [gets stdin]

	puts ""
	puts "	NOTE: TO MAKE THE CIRCULAR GRID, THE CODE NEED TO CALCULATE THE CENTRE OF THIS GRID"
	puts "	    - FOR PROTEIN-MEMBRANE OR MOLECULE-MEMBRANE SYSTEMS THE GEOMETRIC CENTRE OF THE INSERTED COMPONENT IS THE BEST CHOICE"
	puts "	    - AND FOR PURE BILAYER, THE CENTRE OF THE BOX (ENTER '-1' FOR NEXT TWO INPUTS)"
	puts ""
	puts "		### ENTER THE RESIDUE NUMBER OF THE STARTING RESIDUE ###"
	set start_res [gets stdin]
	puts ""
	puts "		### ENTER THE RESIDUE NUMBER OF THE LAST RESIDUE ###"
	set end_res [gets stdin]
	puts ""

	if { $start_res != "-1" || $end_res != "-1" } {
		puts "	THE CIRCULAR GRID WILL BE CENTRED AT THE AVERAGE GEOMETRIC CENTRE OF $start_res TO $end_res AVERAGED OVER $start_frame TO $end_frame"
	} else {
		puts "	THE CIRCULAR GRID WILL BE CENTRED AT THE AVERAGE CENTRE OF THE BOX"
	}

	puts ""
	puts "				#### ANALYSING THE INPUTS ####"

	# DETERMING THE LIPID TYPES

	set lipid ""
	set ind 0
	set pr [open "$prmtop" "r"]
	set dapr [read $pr]
	close $pr

	set k 0

	while { [lindex $dapr $k] != "RESIDUE_LABEL" } {
		incr k
	}
	incr k 2

	
	while { [string length [lindex $dapr $k]] == 2 && [lindex $dapr $k] != "K+" && [lindex $dapr $k] != "Cl-" && [lindex $dapr $k] != "Na+" } {
		set res1 [lindex $dapr $k]
		set res2 [lindex $dapr [expr { $k + 1 }]]
		set res3 [lindex $dapr [expr { $k + 2 }]]
	
		set k1 0
		set count 0
		while { $k1 < [llength $lipid] } {
			if { [lindex $lipid $k1] == $res1 && [lindex $lipid [expr { $k1 + 1 }]] == $res2 && [lindex $lipid [expr { $k1 + 2 }]] == $res3 } {
				incr count
			}
			incr k1 3
		}
		if { $count == 0 } {
			set lipid [linsert $lipid $ind $res1]
	 		incr ind
			set lipid [linsert $lipid $ind $res2]
	 		incr ind
			set lipid [linsert $lipid $ind $res3]
	 		incr ind
		}
		incr k 3
	}
		

	puts ""
	puts "		### DOES THE LIPID BILAYER CONTAIN CHOLESTROL? (y/n) ###"
	set ans [gets stdin]

	if { $ans == "y" || $ans == "Y" } {
		set nlipid [llength $lipid]
		set i 0
		set j 0
		while { $i < $nlipid } {
			set lt($j) [lindex $lipid [expr { $i + 1 }]]
			incr j
			incr i 3
		}
		set lt($j) "CHL"
		incr j
		set lipid [linsert $lipid $ind CHL]
	} else {
		set nlipid [llength $lipid]
		set i 0
		set j 0
		while { $i < $nlipid } {
			set lt($j) [lindex $lipid [expr { $i + 1 }]]
			incr j
			incr i 3
		}
	}						
	set ndlt $j

	set ulres [res_no $prmtop $crd]

	puts "				#### CHECK CAREFULLY :: LIPID 14 LIPID DIVISION IN PRMTOP = $lipid WITH THE RESIDUE 1 TO $ulres BELONG TO THE UPPER LEAFLET ####"
	puts "				#### IF CORRECT PRESS Y ####"
	
	set inp [gets stdin]

	if { $inp != "y" && $inp != "Y" } {
		puts "				#### ERROR IN AUTOMATIC DETERMINATION OF THE SYSTEM LIPIDS ####"
		puts "				#### MANUALLY ENTER THE LIPID 14 DIVISION OF THE LIPIDS (T H T ..... FORMAT] ####"
		set lipinp [gets stdin]
		set lipid $lipinp
		puts ""
		puts "					#### INPUT LIPID 14 DIVISION IS $lipid ####"
		puts ""
		set nlipid [llength $lipid]
		set i 0
		set j 0
		while { $i < $nlipid } {
			set lt($j) [lindex $lipid [expr { $i + 1 }]]
			incr j
			incr i 3
		}
		if { $ans == "Y" || $ans == "y" } {
			set nj [expr { $j - 1 }]
			set lt($nj) "CHL"
		}
		set ndlt $j
		puts "					#### IS THE $ulres THE POSITION OF UPPER LEAFLET LAST RESIDUE? (Y/N) ####"
		set inpp [gets stdin]
		if { $inpp == "n" || $inpp == "N" } {
			puts ""
			puts "				#### ENTER THE INPUT COORDINATE FILE (AT t =0) ####"
			set inpcrd [gets stdin]
			set ulres [res_no $prmtop $inpcrd]
			puts "				#### IS $ulres THE RIGHT VALUE? (Y/N) ####"
			set inpp [gets stdin]
			if { $inpp == "n" || $inpp == "N" } {
				puts "				#### ENTER THE VALUE MANUALLY ####"
				set ulres [gets stdin]
			}
		}	
	} 

	# ******************************************************************************************************************************************************************************#

	# GRID INPUTS 

	puts "				#### ENTER THE RADIUS OF THE CIRCULAR MESH ####"
	set rstep [gets stdin]
	puts ""
	
	puts "				#### ENTER THE STARTING POSITION OF THE CIRCULAR MESH ####"
	set sstep [gets stdin]
	puts ""	

	puts "				#### ENTER THE ANGLE IN DEGREE OF THE SMALLEST ARC ON THE CIRCULAR MESH ####"
	puts "				**** FOR 1D WITH RADIUS AS THE ONLY VARIABLE, INPUT 360.0 ****"
	set thestep [gets stdin]
	puts ""
	
	set ntstep [expr { 360.0 / $thestep }]
	set ntstep [format "%.0f" $ntstep]

	# ******************************************************************************************************************************************************************************#

	# GRID CONSTRUCTION

	set nres [expr { $end_res - $start_res + 1 }]

	set nframe [expr { ($end_frame - $start_frame) / $step }]

	# DETERMING THE AVERAGE CENTER OF MASS OF THE PROTEIN

	set xcom 0.0
	set ycom 0.0
	set zcom 0.0
	set xula 0.0
	set yula 0.0
	set zula 0.0
	set xlla 0.0
	set ylla 0.0
	set zlla 0.0	
	set avg_bx 0.0
	set avg_by 0.0
	set avg_bz 0.0

	puts "				#### DETERMING THE GEOMETRIC CENTRE OF THE RESIDUES $start_res to $end_res AND THE AVERAGE POSITION OF P31'S ON THE UPPER AND THE LOWER LEAFLET ####"
	puts ""
	
	set gc [open "geometric_centre.txt" "w"]

	for {set i $start_frame} {$i < $end_frame} {incr i $step} {

		puts "				#### FRAME $i ####"
	
		# EXECUTING CPPTRAJ

		set f [open "input" "w"]

		puts $f "trajin $crd $i $i"
		puts $f "trajout output.pdb pdb"
		puts $f "trajout output.rst rst"
		puts $f "go"

		close $f

		exec cpptraj -p $prmtop -i input

		set g [open "output.pdb" "r"]
		set data [read $g]
		close $g

		set g1 [open "output.rst" "r"]
		set data1 [read $g1]
		close $g1

		set k 0

		while { $k < [llength $data1] } {
			incr k
		}

		set box_x [lindex $data1 [expr { $k - 3 }]]
		set avg_bx [expr { $avg_bx + $box_x }] 
		set box_y [lindex $data1 [expr { $k - 2 }]]
		set avg_by [expr { $avg_by + $box_y }]
		set box_z [lindex $data1 [expr { $k - 1 }]]
		set avg_bz [expr { $avg_bz + $box_z }]

		set k 0
		set natom 0
		set natom1 0
		set natom2 0

		set xps 0.0
		set yps 0.0
		set zps 0.0
		set xul 0.0
		set xll 0.0
		set yul 0.0
		set yll 0.0
		set zul 0.0
		set zll 0.0

		while { $k < [llength $data] } {
			if { [lindex $data $k] == "ATOM" } {
				if { [lindex $data [expr { $k + 4 }]] >= $start_res && [lindex $data [expr { $k + 4 }]] <= $end_res } {
					set x1 [lindex $data [expr { $k + 5 }]]
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
							set y2 $yori
							set t 0
							while { [string range $y2 $t $t] != "." } {
								incr t
							}
							set yori [string range $y2 0 [expr { $t + 3 }]]
							set zori [string range $y2 [expr { $t + 4 }] end]
						} else {
							set zori [lindex $data [expr { $k + 6 }]] 
						}		
					} else { 
						set xori $x1
						set y1 [lindex $data [expr { $k + 6 }]]
						set sy1 [string length $y1]
						if { $sy1 > 8 } {
							set t 0
							while { [string range $y1 $t $t] != "." } {
								incr t
							}
							set yori [string range $y1 0 [expr { $t + 3 }]]
							set zori [string range $y1 [expr { $t + 4 }] end]
						} else {
							set yori [lindex $data [expr { $k + 6 }]]
							set zori [lindex $data [expr { $k + 7 }]]
						}
					}
					set xps [expr { $xps + $xori }]
					set yps [expr { $yps + $yori }]
					set zps [expr { $zps + $zori }]
					incr natom
				} elseif { [lindex $data [expr { $k + 4 }]] >= 1 && [lindex $data [expr { $k + 4 }]] <= $ulres } {
					if { [lindex $data [expr { $k + 2 }]] == "P31" } {
						set x1 [lindex $data [expr { $k + 5 }]]
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
								set y2 $yori
								set t 0
								while { [string range $y2 $t $t] != "." } {
									incr t
								}
								set yori [string range $y2 0 [expr { $t + 3 }]]
								set zori [string range $y2 [expr { $t + 4 }] end]
							} else {
								set zori [lindex $data [expr { $k + 6 }]] 
							}		
						} else { 
							set xori $x1
							set y1 [lindex $data [expr { $k + 6 }]]
							set sy1 [string length $y1]
							if { $sy1 > 8 } {
								set t 0
								while { [string range $y1 $t $t] != "." } {
									incr t
								}
								set yori [string range $y1 0 [expr { $t + 3 }]]
								set zori [string range $y1 [expr { $t + 4 }] end]
							} else {
								set yori [lindex $data [expr { $k + 6 }]]
								set zori [lindex $data [expr { $k + 7 }]]
							}
						}
						set xul [expr { $xul + $xori }]
						set yul [expr { $yul + $yori }]
						set zul [expr { $zul + $zori }]
						incr natom1
					}
				} elseif { [lindex $data [expr { $k + 4 }]] > $ulres && [lindex $data [expr { $k + 4 }]] < $start_res } {
					if { [lindex $data [expr { $k + 2 }]] == "P31" } {
						set x1 [lindex $data [expr { $k + 5 }]]
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
								set y2 $yori
								set t 0
								while { [string range $y2 $t $t] != "." } {
									incr t
								}
								set yori [string range $y2 0 [expr { $t + 3 }]]
								set zori [string range $y2 [expr { $t + 4 }] end]
							} else {
								set zori [lindex $data [expr { $k + 6 }]] 
							}		
						} else { 
							set xori $x1
							set y1 [lindex $data [expr { $k + 6 }]]
							set sy1 [string length $y1]
							if { $sy1 > 8 } {
								set t 0
								while { [string range $y1 $t $t] != "." } {
									incr t
								}
								set yori [string range $y1 0 [expr { $t + 3 }]]
								set zori [string range $y1 [expr { $t + 4 }] end]
							} else {
								set yori [lindex $data [expr { $k + 6 }]]
								set zori [lindex $data [expr { $k + 7 }]]
							}
						}
						set xll [expr { $xll + $xori }]
						set yll [expr { $yll + $yori }]
						set zll [expr { $zll + $zori }]
						incr natom2
					}
				}
			}
			incr k
		}

		if { $natom == 0 } {
			set natom 1
		}
		if { $natom1 == 0 } {
			set natom1 1
		}
		if { $natom2 == 0 } {
			set natom2 1
		}
		set xcom [expr { $xcom + ($xps / $natom) }]
		set ycom [expr { $ycom + ($yps / $natom) }]
		set zcom [expr { $zcom + ($zps / $natom) }]
		set xula [expr { $xula + ($xul / $natom1) }]
		set yula [expr { $yula + ($yul / $natom1) }]
		set zula [expr { $zula + ($zul / $natom1) }]
		set xlla [expr { $xlla + ($xll / $natom2) }]
		set ylla [expr { $ylla + ($yll / $natom2) }]
		set zlla [expr { $zlla + ($zll / $natom2) }]
	}

	set xcom [expr { $xcom / $nframe }]
	set ycom [expr { $ycom / $nframe }]
	set zcom [expr { $zcom / $nframe }]
	set xula [expr { $xula / $nframe }]
	set yula [expr { $yula / $nframe }]
	set zula [expr { $zula / $nframe }]
	set xlla [expr { $xlla / $nframe }]
	set ylla [expr { $ylla / $nframe }]
	set zlla [expr { $zlla / $nframe }]
	set avg_bx [expr { $avg_bx / $nframe }]
	set avg_by [expr { $avg_by / $nframe }]
	set avg_bz [expr { $avg_bz / $nframe }]

	if { $start_res == -1 && $end_res == -1 } {
		set xcom [expr { $avg_bx / 2.0 }]
		set ycom [expr { $avg_by / 2.0 }]
		set zcom [expr { $avg_bz / 2.0 }]
	}

	puts $gc "#xcoord_centre	ycoord_centre	zcoord_centre	avg_z_lower_leaflet	avg_z_upper_leaflet"
	puts $gc "$xcom	$ycom	$zcom	$zula	$zlla"
	close $gc

	if { $avg_bx < $avg_by } {
		set tsteps [expr { $avg_bx / 2.1 }]
		set tsteps [expr { round($tsteps/$rstep) }]
		set tsteps [expr { $rstep * $tsteps }]
	} else { 
		set tsteps [expr { $avg_by / 2.1 }]
		set tsteps [expr { round($tsteps/$rstep) }]
		set tsteps [expr { $rstep * $tsteps }]
	}

	set ngrid [expr { (($tsteps - $sstep) / $rstep) * $ntstep }] 
	
	set ngrid [format "%.0f" $ngrid]

	puts ""
	puts "				***************************************************************************************************************************************"
	puts "				#### IMP NOTE :: THE CODE WILL FORM A GRID OF $ngrid SIZE, HIGHER NUMBER WILL REDUCE THE SPEED APPROXIMATELY BY N^2 PER FRAME #####"
	puts "				***************************************************************************************************************************************"

	# ******************************************************************************************************************************************************************************#
	
	# LOCAL AREA PER LIPID AND BILAYER THICKNESS

	if { $MRES == 1 } {
		puts ""
		puts "			***** WELCOME TO THE LOCAL MEMBRANE PROPERTIES CALCULATIONS SUITE OF AMBAT *****"
		puts ""
		puts "" 
		puts "			**** CALCULATION OF AREA PER LIPID AND BILAYER THICKNESS WILL BE CENTERED AT AT $xcom , $ycom , $zcom ****"

		# DETERMINING THE AREA PER LIPID ALONG THE XY LEAFLET BY MAKING USE OF THE CONCENTRIC CIRCULAR MESH CENTERED ON THE GEOMETRIC CENTER OF THE PROTEIN

		for {set i 0} {$i < 1000} {incr i} {
			for {set ij 0} {$ij < 100} {incr ij} {	
				set nrul($i,$ij) 0.0
				set nrll($i,$ij) 0.0
				set blt($i,$ij) 0.0
				for {set j $start_frame} {$j < $end_frame} {incr j $step} {
					set nrulsd($j,$i,$ij) 0.0
					set nrllsd($j,$i,$ij) 0.0
					set bltsd($j,$i,$ij) 0.0
				}
			}
		}

		set th [open "bilayer_thickness.txt" "w"]
		puts $th "#radius	theta	bilayer_thickness	SD_bilayer_thickness"	

		puts ""
		puts ""
		puts "				#### FORMING A CIRCULAR MESH OF $rstep A radius DIVIDED INTO $ntstep EQUAL ARCS CENTERED AT $xcom,$ycom,$zcom STARTING FROM $sstep ####"
	
		set tes [open "test" "w"]

		set nframe1 0
		set nframe2 0
		set nframe3 0
		set countframe1 0
		set countframe2 0
		set countframe3 0

		for {set i $start_frame} {$i < $end_frame} {incr i $step} {

			puts "				#### FRAME $i ####"
	
			# EXECUTING CPPTRAJ

			set f [open "input" "w"]

			puts $f "trajin $crd $i $i"
			puts $f "trajout output.pdb pdb"
			puts $f "trajout output.rst rst"
			puts $f "go"

			close $f

			exec cpptraj -p $prmtop -i input

			# DETERMING THE CENTER OF MASS OF EACH RESIDUE 

			com $prmtop $lipid

			set g [open "center_of_mass" "r"]
			set data [read $g]
			close $g
	
			set g1 [open "output.rst" "r"]
			set data1 [read $g1]
			close $g1

			set k 0

			while { $k < [llength $data1] } {
				incr k
			}

			set box_x [lindex $data1 [expr { $k - 3 }]]
			set box_y [lindex $data1 [expr { $k - 2 }]]
			set box_z [lindex $data1 [expr { $k - 1 }]]

			set j 0

			#set tsteps 4.0

			set countframe1 0
			set countframe2 0
			set countframe3 0

			for { set r $sstep } { $r <= $tsteps } { set r [expr { $r + $rstep }] } {
				puts "				#### RADIUS OF $r ang FROM THE PROTEIN GEOMETRIC CENTRE ####"
				set ij 0
				for {set theta 0.0} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set reslist ""
					set reslist1 ""
					set reslist2 ""
					#set reslist ""
					set ind 0
					set ind1 0
					set ind2 0
					set rr($j) $r
					set the($j,$ij) $theta 
					set k 0
					while { $k < [llength $data] } {
						set x1 [lindex $data [expr { $k + 0 }]]
						set y1 [lindex $data [expr { $k + 1 }]]
						set z [lindex $data [expr { $k + 2 }]]
						set res1 [lindex $data [expr { $k + 3 }]]

						for { set imx [expr { -1 * $box_x }] } {$imx <= $box_x} {set imx [expr { $imx + $box_x }]} {
							for { set imy [expr { -1 * $box_y }] } {$imy <= $box_y} {set imy [expr { $imy + $box_y }]} {
								set x [expr { $x1 + $imx }]
								set y [expr { $y1 + $imy }]

								set xc [expr { $x - $xcom }]
								set yc [expr { $y - $ycom }]
								set tan [expr { $yc / $xc }]

								set tan [expr { atan($tan) }]
								set tan [expr { ($tan * 180.0) / 3.14 }]

								# 2ND QUARDRANT

								if { $yc > 0 && $xc < 0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 3RD QUARDRANT
	 
								if { $xc < 0.0 && $yc < 0.0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 4TH QUARDRANT

								if { $xc > 0.0 && $yc < 0.0 } {
									set tan [expr { $tan + 360.0 }]
								}
					
								set xc [expr { $xc * $xc }]	
								set yc [expr { $yc * $yc }]
								set zc [expr { $z - $zcom }]
								set zc [expr { $zc * $zc }]

								set rc [expr { sqrt($xc + $yc) }]	

								if { $rc >= [expr { $r - $rstep }] && $rc <= $r } {
									if { $tan >= [expr { $theta - $thestep }] && $tan <= $theta } {
										if { [lindex $data [expr { $k + 5 }]] <= $ulres } {
											set nrul($j,$ij) [expr { $nrul($j,$ij) + 1.0 }]
											set nrulsd($i,$j,$ij) [expr { $nrulsd($i,$j,$ij) + 1.0 }]
											incr countframe1 
											set reslist1 [linsert $reslist $ind1 [lindex $data [expr { $k + 4 }]]]
											incr ind1
										} else { 
											set nrll($j,$ij) [expr { $nrll($j,$ij) + 1.0 }]
											set nrllsd($i,$j,$ij) [expr { $nrllsd($i,$j,$ij) + 1.0 }]
											incr countframe2  
											set reslist2 [linsert $reslist $ind2 [lindex $data [expr { $k + 4 }]]]
											incr ind2
										}
									}
								}
							}
						}
						incr k 6
					}
					if { $reslist1 != "" && $reslist2 != "" } { 
						set thick1 [ed $prmtop output.pdb $reslist1]
						set thick2 [ed $prmtop output.pdb $reslist2]
						set max1 [lindex $thick1 0]
						set max2 [lindex $thick2 0]
						set thickness [expr { abs($max1 - $max2) }]
						set blt($j,$ij) [expr { $blt($j,$ij) + $thickness }]
						set bltsd($i,$j,$ij) $thickness
						puts $tes "$i $r $theta $nrulsd($i,$j,$ij) $nrllsd($i,$j,$ij)"
						incr countframe3
					}
					incr ij
				}
				incr j
			}
			if { $countframe1 != 0 } {
				incr nframe1
			}
			if { $countframe2 != 0 } {
				incr nframe2
			}
			if { $countframe3 != 0 } {
				incr nframe3
			}
		}
		close $tes

		set nmesh $j
		set h [open "mesh_area.txt" "w"]
		puts $h "#radius	theta	apl_UL	SD_apl_UL	apl_LL	SD_apl_ll	num_lipid_UL	SD_num_lipids_UL	num_lipid_LL	SD_num_lipids_LL	Avg_APL	SD_Avg_APL"
		set rr(-1) [expr { $sstep - $rstep }]
		set pi 3.14
		for {set j 0} {$j < $nmesh} {incr j} {
			for {set ij 1} {$ij <= $ntstep} {incr ij} { 
				set nrul($j,$ij) [expr { $nrul($j,$ij) / $nframe1 }]
				set nrll($j,$ij) [expr { $nrll($j,$ij) / $nframe2 }]
				set blt($j,$ij) [expr { $blt($j,$ij) / $nframe3 }]
				set oj [expr {$j - 1 }]
				set or $rr($oj)
				if { $nrul($j,$ij) > 0.0 } {
					set apl_ul [expr { ($pi * (($rr($j) * $rr($j)) - ($or * $or))) / ($nrul($j,$ij) * $ntstep) }]
				} else {
					set apl_ul -1
				}
				if { $nrll($j,$ij) > 0.0 } {
					set apl_ll [expr { ($pi * (($rr($j) * $rr($j)) - ($or * $or))) / ($nrll($j,$ij) * $ntstep) }]
				} else {
					set apl_ll -1
				}
				set sd1 0.0
				set sd2 0.0
				set sd3 0.0	
				set sd_nlul 0.0
				set sd_nlll 0.0
				set nterm1 0
				set nterm2 0
				for {set i $start_frame} {$i < $end_frame} {incr i $step} {
					set diff [expr { $nrulsd($i,$j,$ij) - $nrul($j,$ij) }]
					set diff [expr { $diff * $diff }]
					set sd_nlul [expr { $sd_nlul + $diff }]

					set diff [expr { $nrllsd($i,$j,$ij) - $nrll($j,$ij) }]
					set diff [expr { $diff * $diff }]
					set sd_nlll [expr { $sd_nlll + $diff }]

					if { $nrulsd($i,$j,$ij) > 0 } {
						set apl_ulsd [expr { ($pi * (($rr($j) * $rr($j)) - ($or * $or))) / ($nrulsd($i,$j,$ij)*$ntstep) }]
						set diff [expr { $apl_ul - $apl_ulsd }]
						set diff [expr { $diff * $diff }]
						set sd1 [expr { $sd1 + $diff }]
						incr nterm1 
					} 
					if { $nrllsd($i,$j,$ij) > 0 } {
						set apl_llsd [expr { ($pi * (($rr($j) * $rr($j)) - ($or * $or))) / ($nrllsd($i,$j,$ij)*$ntstep) }]
						set diff [expr { $apl_ll - $apl_llsd }]
						set diff [expr { $diff * $diff }]
						set sd2 [expr { $sd2 + $diff }]
						incr nterm2 
					}
	
					set diff [expr { $bltsd($i,$j,$ij) - $blt($j,$ij) }]
					set diff [expr { $diff * $diff }]
					set sd3 [expr { $sd3 + $diff }]
				}
				set sd_nlul [expr { $sd_nlul / $nframe1 }]
				set sd_nlul [expr { sqrt($sd_nlul) }]
				set sd_nlll [expr { $sd_nlll / $nframe2 }]
				set sd_nlll [expr { sqrt($sd_nlll) }]
				 
				if { $nterm1 > 0 } {
					set sd1 [expr { $sd1 / $nterm1 }]
				} else {
					set sd1 0.0
				}
				if { $nterm2 > 0 } {			
					set sd2 [expr { $sd2 / $nterm2 }]
				} else {
					set sd2 0.0
				}
				set sd3 [expr { $sd3 / $nframe3 }]
				set sd3 [expr { sqrt($sd3) }]
				set sd3 [format "%.2f" $sd3]	
				set sd [expr { sqrt($sd1 + $sd2) }]
				set sd [format "%.2f" $sd]	
				set sd1 [expr { sqrt($sd1) }]
				set sd1 [format "%.2f" $sd1]	
				set sd2 [expr { sqrt($sd2) }]
				set sd2 [format "%.2f" $sd2]

				set avg_apl [expr { ($apl_ul + $apl_ll) / 2.0 }]
				set avg_apl [format "%.2f" $avg_apl]		

				set nrul($j,$ij) [format "%.2f" $nrul($j,$ij)]
				set nrll($j,$ij) [format "%.2f" $nrll($j,$ij)]
				set blt($j,$ij) [format "%.2f" $blt($j,$ij)]
				set apl_ul [format "%.2f" $apl_ul]
				set apl_ll [format "%.2f" $apl_ll]  

				set sd_nlul [format "%.2f" $sd_nlul]
				set sd_nlll [format "%.2f" $sd_nlll]

				if { $sd1 < 30.0 && $sd2 < 30.0 } {
					if { $apl_ul > 0.0 && $apl_ll > 0.0 } {
						puts $h "$rr($j)	$the($j,$ij)	$apl_ul	$sd1	$apl_ll	$sd2	$nrul($j,$ij)  $sd_nlul	 $nrll($j,$ij)	$sd_nlll	$avg_apl	$sd"
					} elseif { $apl_ul <  0.0 && $apl_ll > 0.0 } {
						puts $h "$rr($j)	$the($j,$ij)	-1	-1	$apl_ll	$sd2	-1	-1	$nrll($j,$ij)  $sd_nlll	-1	-1	"
					} elseif { $apl_ul > 0.0 && $apl_ll < 0.0 } {
						puts $h "$rr($j)	$the($j,$ij)	$apl_ul	$sd1	-1	-1	$nrul($j,$ij)  $sd_nlul	-1	-1	-1	-1	"
					} elseif { $apl_ul < 0.0 && $apl_ll < 0.0 } {
						puts $h "$rr($j)	$the($j,$ij)	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	"
					}
				} elseif { $sd1 < 30.0 && $sd2 > 30.0 } {
					puts $h "$rr($j)	$the($j,$ij)	$apl_ul	$sd1	-1	-1	$nrul($j,$ij)  $sd_nlul	-1	-1	-1	-1	"
				} elseif { $sd1 > 30.0 && $sd2 < 30.0 } {
					puts $h "$rr($j)	$the($j,$ij)	-1	-1	$apl_ll	$sd2	-1	-1	$nrll($j,$ij)  $sd_nlll	-1	-1	"
				} elseif { $sd1 > 30.0 && $sd2 > 30.0 } {
					puts $h "$rr($j)	$the($j,$ij)	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	"
				}
	
				if { $sd3 < 15.0 && $blt($j,$ij) != 0.00 } {
					puts $th "$rr($j)	$the($j,$ij)	$blt($j,$ij)	$sd3"
				} else {
					puts $th "$rr($j)	$the($j,$ij)	-1	-1"
				}
			}
		}
		close $h
		close $th
		file delete test

		# FOR A VIDEO

		exec mkdir -p video_files_blt

		set vi 0
		for {set i $start_frame} {$i < $end_frame} {incr i $step} {
			set v [open "$vi.txt" "w"]
			for {set j 0} {$j < $nmesh} {incr j} {
				for {set ij 1} {$ij <= $ntstep} {incr ij} {
					puts $v "$rr($j)	$the($j,$ij)	$bltsd($i,$j,$ij)"
				}
			}
			close $v
			exec mv $vi.txt ./video_files_blt
			incr vi
		} 
	}

	# ******************************************************************************************************************************************************************************#

	# LOCAL CURVATURE

	if { $CRES == 1 } {
		puts "			***** WELCOME TO THE LOCAL MEMBRANE PROPERTIES CALCULATIONS SUITE OF AMBAT *****"
		puts ""
		puts "" 
		puts "			#### ANALYSING THE LOCAL MEMBRANE CURVATURE ####"
		puts ""
		curvature $crd $prmtop $xcom $ycom $zcom $avg_bx $avg_by $start_frame $end_frame $step $sstep $rstep $thestep $ulres
	}

	# ******************************************************************************************************************************************************************************#

	# ION FLUX
	
	if { $IRES == 1 } {
		puts ""
		puts "			***** WELCOME TO THE LOCAL MEMBRANE PROPERTIES CALCULATIONS SUITE OF AMBAT *****"
		puts ""
		puts "" 
		puts "			#### DETERMING THE ION FLUX THROUGH THE CHANNEL ####"
		puts ""
		flux $prmtop $crd $zlla $zula $start_frame $end_frame $step $ionn "number_$ionn"
	}

	# ******************************************************************************************************************************************************************************#

	# WATER FLUX
	
	if { $WRES == 1 } {
		puts ""
		puts "			***** WELCOME TO THE LOCAL MEMBRANE PROPERTIES CALCULATIONS SUITE OF AMBAT *****"
		puts ""
		puts "" 
		puts "			#### DETERMING THE WATER FLUX THROUGH THE CHANNEL ####"
		puts ""
		flux $prmtop $crd $zlla $zula $start_frame $end_frame $step O "number_water"
	}

	# ******************************************************************************************************************************************************************************#

	# PROTEIN PROFILE

	if { $PRES == 1 } {

		puts ""
		puts "			***** WELCOME TO THE PROTEIN PROFILE SUITE OF AMBAT *****"
		puts ""
		puts "" 

		# EXECUTING CPPTRAJ

		set f [open "input" "w"]

		puts $f "trajin $crd $end_frame $end_frame"
		puts $f "strip !:$start_res-$end_res"
		puts $f "trajout output.pdb pdb"
		puts $f "go"

		close $f

		exec cpptraj -p $prmtop -i input

		# 1D PROFILE
		
		puts ""
		puts "		### DETERMINING THE 1D PROFILE FROM THE LAST FRAME ###"
		puts ""

		protein_profile_1D output.pdb 2

		# 2D SPREAD

		puts ""
		puts "		### DETERMINING THE 2D SPREAD FROM THE LAST FRAME ###"
		puts ""

		protein_profile output.pdb 2

		# 2D/3D (DEPENDING UPON ZSLICE)
		puts ""
		puts "		### DETERMINING THE 3D PROFILE ###"
		puts ""

		puts "$prmtop $crd $xcom $ycom $zcom $start_frame $end_frame $step $start_res $end_res $sstep $tsteps $rstep $thestep $zula $zlla $zslice"

		pp $prmtop $crd $xcom $ycom $zcom $start_frame $end_frame $step $start_res $end_res $sstep $tsteps $rstep $thestep $zlla $zula $zslice
		average $sstep $tsteps $rstep $thestep
	}
	
	# ******************************************************************************************************************************************************************************#

	# LIPID PROFILE

	if { $LRES == 1 } {

		puts ""
		puts "			***** WELCOME TO THE LIPID PROFILE SUITE OF AMBAT *****"
		puts ""
		puts "" 
		
		puts "			**** CALCULATION OF LIPID NUMBER PROFILE WILL BE CENTERED AT AT $xcom , $ycom , $zcom ****"

		# DETERMINING THE AREA PER LIPID ALONG THE XY LEAFLET BY MAKING USE OF THE CONCENTRIC CIRCULAR MESH CENTERED ON THE GEOMETRIC CENTER OF THE PROTEIN

		for {set i 0} {$i < 1000} {incr i} {
			for {set ij 0} {$ij < 100} {incr ij} {	
				for {set j 0} {$j < $ndlt} {incr j} {
					set nltul($j,$i,$ij) 0.0
					set nltll($j,$i,$ij) 0.0
					for {set ji $start_frame} {$ji < $end_frame} {incr ji $step} {
						set nltulsd($ji,$j,$i,$ij) 0.0
						set nltllsd($ji,$j,$i,$ij) 0.0
					}
				}
			}
		}

		puts ""
		puts ""
		puts "				#### FORMING A CIRCULAR MESH OF $rstep A radius DIVIDED INTO $ntstep EQUAL ARCS CENTERED AT $xcom,$ycom,$zcom STARTING FROM $sstep ####"
	
		set tes [open "test" "w"]

		set nframe1 0
		set nframe2 0
		set nframe3 0
		set countframe1 0
		set countframe2 0
		set countframe3 0

		for {set i $start_frame} {$i < $end_frame} {incr i $step} {

			puts "				#### FRAME $i ####"
	
			# EXECUTING CPPTRAJ

			set f [open "input" "w"]

			puts $f "trajin $crd $i $i"
			puts $f "trajout output.pdb pdb"
			puts $f "trajout output.rst rst"
			puts $f "go"

			close $f

			exec cpptraj -p $prmtop -i input

			# DETERMING THE CENTER OF MASS OF EACH RESIDUE 

			com_lp $prmtop $lipid

			set g [open "center_of_mass" "r"]
			set data [read $g]
			close $g
	
			set g1 [open "output.rst" "r"]
			set data1 [read $g1]
			close $g1

			set k 0

			while { $k < [llength $data1] } {
				incr k
			}

			set box_x [lindex $data1 [expr { $k - 3 }]]
			set box_y [lindex $data1 [expr { $k - 2 }]]
			set box_z [lindex $data1 [expr { $k - 1 }]]

			set j 0

			set countframe1 0
			set countframe2 0
			set countframe3 0

			for { set r $sstep } { $r < $tsteps } { set r [expr { $r + $rstep }] } {
				puts "				#### RADIUS OF $r ang FROM THE PROTEIN GEOMETRIC CENTRE ####"
				set ij 0
				for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set reslist ""
					set reslist1 ""
					set reslist2 ""
					#set reslist ""
					set ind 0
					set ind1 0
					set ind2 0
					set rr($j) $r
					set the($j,$ij) $theta 
					set k 0
					while { $k < [llength $data] } {
						set x1 [lindex $data [expr { $k + 0 }]]
						set y1 [lindex $data [expr { $k + 1 }]]
						set z [lindex $data [expr { $k + 2 }]]
						set res1 [lindex $data [expr { $k + 3 }]]

						for { set imx [expr { -1 * $box_x }] } {$imx <= $box_x} {set imx [expr { $imx + $box_x }]} {
							for { set imy [expr { -1 * $box_y }] } {$imy <= $box_y} {set imy [expr { $imy + $box_y }]} {
								set x [expr { $x1 + $imx }]
								set y [expr { $y1 + $imy }]

								set xc [expr { $x - $xcom }]
								set yc [expr { $y - $ycom }]
								set tan [expr { $yc / $xc }]

								set tan [expr { atan($tan) }]
								set tan [expr { ($tan * 180.0) / 3.14 }]

								# 2ND QUARDRANT

								if { $yc > 0 && $xc < 0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 3RD QUARDRANT
	 
								if { $xc < 0.0 && $yc < 0.0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 4TH QUARDRANT

								if { $xc > 0.0 && $yc < 0.0 } {
									set tan [expr { $tan + 360.0 }]
								}
					
								set xc [expr { $xc * $xc }]	
								set yc [expr { $yc * $yc }]
								set zc [expr { $z - $zcom }]
								set zc [expr { $zc * $zc }]

								set rc [expr { sqrt($xc + $yc) }]	

								if { $rc >= [expr { $r - $rstep }] && $rc <= $r } {
									if { $tan >= [expr { $theta - $thestep }] && $tan <= $theta } {
										if { [lindex $data [expr { $k + 5 }]] <= $ulres } {
											set clt [lindex $data [expr { $k + 6 }]]
											for {set ci 0} {$ci < $ndlt} {incr ci} {
												if { $clt == $lt($ci) } {
													set cltn $ci
												}
											}
											set nltul($cltn,$j,$ij) [expr { $nltul($cltn,$j,$ij) + 1.0 }]
											set nltulsd($i,$cltn,$j,$ij) [expr { $nltulsd($i,$cltn,$j,$ij) + 1.0 }]
										} else { 
											set clt [lindex $data [expr { $k + 6 }]]
											for {set ci 0} {$ci < $ndlt} {incr ci} {
												if { $clt == $lt($ci) } {
													set cltn $ci
												}
											}
											set nltll($cltn,$j,$ij) [expr { $nltll($cltn,$j,$ij) + 1.0 }]
											set nltllsd($i,$cltn,$j,$ij) [expr { $nltllsd($i,$cltn,$j,$ij) + 1.0 }]
										}
									}
								}
							}
						}
						incr k 7
					}
					incr ij
				}
				incr j
			}
		}
		close $tes

		set nmesh $j
		for {set i 0} {$i < $ndlt} {incr i} {
			set h [open "3D_lipid_profile_$i.txt" "w"]
			puts $h "#r	theta	lipids_UL	SD_lipids_UL	lipids_LL	SD_lipids_LL"
			puts $h "# $lt($i)"
			for {set j 0} {$j < $nmesh} {incr j} {
				for {set ij 0} {$ij < $ntstep} {incr ij} {
					set nltul($i,$j,$ij) [expr { $nltul($i,$j,$ij) / $nframe }]
					set nltll($i,$j,$ij) [expr { $nltll($i,$j,$ij) / $nframe }]
					set sd_nlul 0.0
					set sd_nlll 0.0
					set nterm1 0
					set nterm2 0
					for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
						set diff [expr { $nltulsd($fr,$i,$j,$ij) - $nltul($i,$j,$ij) }]
						set diff [expr { $diff * $diff }]
						set sd_nlul [expr { $sd_nlul + $diff }]

						set diff [expr { $nltllsd($fr,$i,$j,$ij) - $nltll($i,$j,$ij) }]
						set diff [expr { $diff * $diff }]
						set sd_nlll [expr { $sd_nlll + $diff }]
					}
					set sd_nlul [expr { $sd_nlul / $nframe }]
					set sd_nlul [expr { sqrt($sd_nlul) }]
					set sd_nlul [format "%.2f" $sd_nlul]

					set sd_nlll [expr { $sd_nlll / $nframe }]
					set sd_nlll [expr { sqrt($sd_nlll) }]
					set sd_nlll [format "%.2f" $sd_nlll]

					set normfactor [expr { (($rr($j) * $rr($j)) - (($rr($j)-$rstep)*($rr($j)-$rstep))) / 1000.0 }]
			
					set var1 [expr { $nltul($i,$j,$ij) / $normfactor  }]
					set var1 [format "%.2f" $var1]
					set var1c [expr { $sd_nlul / $normfactor }]
					set var1c [format "%.2f" $var1c]
		
					set var2 [expr { $nltll($i,$j,$ij) / $normfactor }]
					set var2 [format "%.2f" $var2]
					set var2c [expr { $sd_nlll / $normfactor }]
					set var2c [format "%.2f" $var2c]

					puts $h "$rr($j)	$the($j,$ij)	$var1	$var1c	$var2	$var2c"
				}
			}
			close $h
		}

		# FOR A VIDEO

		exec mkdir -p video_files_lp


		for {set i 0} {$i < $ndlt} {incr i} {
			set vi 0
			for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
				set v [open "$i.$vi.txt" "w"]
				for {set j 0} {$j < $nmesh} {incr j} {
					for {set ij 0} {$ij < $ntstep} {incr ij} {
						set normfactor [expr { (($rr($j) * $rr($j)) - (($rr($j)-$rstep)*($rr($j)-$rstep))) / 1000.0 }]
						puts $v "$rr($j)	$the($j,$ij)	[expr { $nltulsd($fr,$i,$j,$ij) / $normfactor }]	 [expr { $nltllsd($fr,$i,$j,$ij) / $normfactor }]"
					}
				}
				close $v
				exec mv $i.$vi.txt ./video_files_lp
				incr vi
			}
		}
	}
	# ******************************************************************************************************************************************************************************#

	# WATER PROFILE

	if { $WPRES == 1 } {
		wat_prof $prmtop $crd $xcom $ycom $zcom $start_frame $end_frame $step $sstep $tsteps $rstep $thestep $zlla $zula $zwpslice
	}

	# ******************************************************************************************************************************************************************************#
}

proc res_no { prmtop crd } {
	# TO DETERMINE THE UPPER AND LOWER LAYER RESIDUE NUMBER

	# EXECUTING CPPTRAJ TO GET THE RESIDUE NUMBERS FOR THE UPPER AND LOWER RESIDUES

	set inp [open "input" "w"]

	puts $inp "trajin $crd 1 1"
	puts $inp "trajout output.pdb"
	puts $inp "go"

	close $inp

	exec cpptraj -p $prmtop -i input

	set f [open "output.pdb" "r"]
	set data [read $f]
	close $f

	set k 0
	while { [lindex $data [expr { $k + 2 }]] != "O21" } {
		incr k
	}

	set x1 [lindex $data [expr { $k + 5 }]]
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
			set y2 $yori
			set t 0
			while { [string range $y2 $t $t] != "." } {
				incr t
			}
			set yori [string range $y2 0 [expr { $t + 3 }]]
			set zori [string range $y2 [expr { $t + 4 }] end]
			} else {
			set zori [lindex $data [expr { $k + 6 }]] 
		}		
	} else { 
		set xori $x1
		set y1 [lindex $data [expr { $k + 6 }]]
		set sy1 [string length $y1]
		if { $sy1 > 8 } {
			set t 0
			while { [string range $y1 $t $t] != "." } {
				incr t
			}
			set yori [string range $y1 0 [expr { $t + 3 }]]
			set zori [string range $y1 [expr { $t + 4 }] end]
		} else {
			set yori [lindex $data [expr { $k + 6 }]]
			set zori [lindex $data [expr { $k + 7 }]]
		}
	}

	set zul $zori
		
	set k 0

	while { $k < [llength $data] } {
		if { [lindex $data $k] == "ATOM" && [lindex $data [expr { $k + 2 }]] == "O21" } {
			set x1 [lindex $data [expr { $k + 5 }]]
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
					set y2 $yori
					set t 0
					while { [string range $y2 $t $t] != "." } {
						incr t
					}
					set yori [string range $y2 0 [expr { $t + 3 }]]
					set zori [string range $y2 [expr { $t + 4 }] end]
					} else {
					set zori [lindex $data [expr { $k + 6 }]] 
				}		
			} else { 
				set xori $x1
				set y1 [lindex $data [expr { $k + 6 }]]
				set sy1 [string length $y1]
				if { $sy1 > 8 } {
					set t 0
					while { [string range $y1 $t $t] != "." } {
						incr t
					}
					set yori [string range $y1 0 [expr { $t + 3 }]]
					set zori [string range $y1 [expr { $t + 4 }] end]
				} else {
					set yori [lindex $data [expr { $k + 6 }]]
					set zori [lindex $data [expr { $k + 7 }]]
				}
			}
			if { $zori < [expr { $zul - 25.0 }] } {
				set resid [lindex $data [expr { $k + 4 }]]
				set k [llength $data]
				return [expr { $resid - 2 }]
			}
		}
		incr k
	}
}

proc com { prmtop lipid } { 

	# FOR CALCULATION OF CENTRE OF MASS

	set f [open "$prmtop" "r"]
	set data [read $f]
	close $f

	set g [open "output.pdb" "r"]
	set data1 [read $g]
	close $g

	set h [open "center_of_mass" "w"]

	set k 0

	while { [lindex $data $k] != "MASS" || [lindex $data [expr { $k - 1 }]] != "%FLAG" } {
		incr k
	}
	incr k 2
	set j 0
	while { [lindex $data $k] != "%FLAG" } {
		set mass($j) [lindex $data $k]
		incr j
		incr k
	}

	set k 0
	set j 0
	while { $k < [llength $data1] } {
		set total_M 0.0
		set xcom 0.0
		set ycom 0.0
		set zcom 0.0
		set c1 0
		set c2 0
		set c3 0
		while { [lindex $data1 $k] != "TER" && [lindex $data1 $k] != "END" } {
			if { [lindex $data1 $k] == "ATOM" } {
				set count 0
				set k1 0
				while { $k1 < [llength $lipid] } {
					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid $k1] } {
						incr count
					}
					incr k1
				}
				if { $count != 0 || [lindex $data1 [expr { $k + 3 }]] == "CHL" } {
					if { [lindex $data1 [expr { $k + 3} ]] == [lindex $lipid 0] && $c1 == 0 } {
						set res1 [lindex $data1 [expr { $k + 4 }]]
						incr c1
					}

					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid 1] && $c2 == 0 } {
						set res2 [lindex $data1 [expr { $k + 4 }]]
						incr c2
					}

					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid 2] && $c3 == 0 } {
						set res3 [lindex $data1 [expr { $k + 4 }]]
						incr c3
					}

					if { [lindex $data1 [expr { $k + 3 }]] == "CHL" } {
						set res1 [lindex $data1 [expr { $k + 4 }]]
						set res2 $res1
						set res3 $res1
					}

					set x1 [lindex $data1 [expr { $k + 5 }]]
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
							set y2 $yori
							set t 0
							while { [string range $y2 $t $t] != "." } {
								incr t
							}
							set yori [string range $y2 0 [expr { $t + 3 }]]
							set zori [string range $y2 [expr { $t + 4 }] end]
						} else {
							set zori [lindex $data1 [expr { $k + 6 }]] 
						}		
					} else { 
						set xori $x1
						set y1 [lindex $data1 [expr { $k + 6 }]]
						set sy1 [string length $y1]
						if { $sy1 > 8 } {
							set t 0
							while { [string range $y1 $t $t] != "." } {
								incr t
							}
							set yori [string range $y1 0 [expr { $t + 3 }]]
							set zori [string range $y1 [expr { $t + 4 }] end]
						} else {
							set yori [lindex $data1 [expr { $k + 6 }]]
							set zori [lindex $data1 [expr { $k + 7 }]]
						}
					}
					
					set x $xori
					set y $yori
					set z $zori

					set xcom [expr { $xcom + ($x*$mass($j)) }]
					set ycom [expr { $ycom + ($y*$mass($j)) }]
					set zcom [expr { $zcom + ($z*$mass($j)) }]
					set total_M [expr { $total_M + $mass($j) }]
	
				}
				incr j
			}
			incr k
		}
		if { $total_M > 0.0 } {
			puts $h "[expr { $xcom / $total_M }] [expr { $ycom / $total_M }] [expr { $zcom / $total_M }] $res1 $res2 $res3"
		}
		incr k
	}
	close $h
}

proc ed { prmtop crd reslist } {

	# FOR THE CALCULATION OF BILAYER THICKNESS (ELECTRON DENSITY CALCULATION)

	set f [open "$prmtop" "r"]
	set data [read $f]
	close $f

	set g [open "$crd" "r"]
	set data1 [read $g]
	close $g

	set box_z [lindex $data1 3]

	set bin 0.25
	
	# CHARGE LABEL

	set k2 0
	while { [lindex $data $k2] != "CHARGE" } {
		incr k2
	}
	incr k2
	set k3 $k2

	while { [lindex $data $k3] != "RESIDUE_POINTER" } {
		incr k3
	}
	incr k3
	set k 0
	set k1 0
	
	set re 0
	while { $k < [llength $reslist] } {
		set res [lindex $reslist $k]

		# COORDINATES OF THE RESIDUE

		set at 0
		while { [lindex $data1 $k1] != $res || [lindex $data1 [expr { $k1 - 4 }]] != "ATOM" } {
			incr k1
		}
		set k1 [expr { $k1 - 4 }]

		set fatom [lindex $data [expr { $k3 + $res }]]
		set latom [lindex $data [expr { $k3 + $res + 1}]]
		set natoms [expr { $latom - $fatom }]

		set scharge [expr { $k2 + $fatom }]
		set k4 0
		while { $k4 < $natoms } {
			if { [lindex $data1 $k1] == "ATOM" } {
				set x1 [lindex $data1 [expr { $k1 + 5 }]]
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
						set y2 $yori
						set t 0
						while { [string range $y2 $t $t] != "." } {
							incr t
						}
						set yori [string range $y2 0 [expr { $t + 3 }]]
						set zori [string range $y2 [expr { $t + 4 }] end]
						set atype [lindex $data1 [expr { $k1 + 8 }]] 
					} else {
						set zori [lindex $data1 [expr { $k1 + 6 }]]
						set atype [lindex $data1 [expr { $k1 + 9 }]] 
					}		
				} else { 
					set xori $x1
					set y1 [lindex $data1 [expr { $k1 + 6 }]]
					set sy1 [string length $y1]
					if { $sy1 > 8 } {
						set t 0
						while { [string range $y1 $t $t] != "." } {
							incr t
						}
						set yori [string range $y1 0 [expr { $t + 3 }]]
						set zori [string range $y1 [expr { $t + 4 }] end]
						set atype [lindex $data1 [expr { $k1 + 9 }]] 
					} else {
						set yori [lindex $data1 [expr { $k1 + 6 }]]
						set zori [lindex $data1 [expr { $k1 + 7 }]]
						set atype [lindex $data1 [expr { $k1 + 10 }]]
					}
				}
				if { $atype == "N" } {
					set Ne 7
				} 
				if { $atype == "C" } {
					set Ne 6
				} 
				if { $atype == "O" } {
					set Ne 8
				} 
				if { $atype == "P" } {
					set Ne 15
				} 
				if { $atype == "H" } {
					set Ne 1
				} 
						
				set x($re,$at) $xori
				set y($re,$at) $yori
				set z($re,$at) $zori
				set c($re,$at) [expr { $Ne - ([lindex $data $scharge] ) }]
				incr at	
				
				incr scharge 
				incr k4		
			}			
			incr k1
		}
		set nat($re) $at
		incr re
		incr k
	}
	# DETERMING THE BILAYER THICKNESS BASED ON ELECTRON DENSITY

	set max1 -1000
	set max2 -2000
	set zmax1 -1
	set zmax2 -1
	for {set zm 0} { $zm < $box_z } {set zm [expr { $zm + 0.25 }]} {
		set charge 0
		for {set i 0} {$i < [llength $reslist]} {incr i} {
			for {set j 0} {$j < $nat($i) } {incr j} {
				if { $z($i,$j) >= $zm && $z($i,$j) <= [expr { $zm + 0.25 }] } {
					set charge [expr { $charge + $c($i,$j) }]
				}
			}
		}
		if { $charge > $max1 } {
			set max2 $max1
			set max1 $charge	
			set zmax2 $zmax1
			set zmax1 $zm
		} elseif { $charge > $max2 } {
			set max2 $charge
			set zmax2 $zm
		}
	}
	return [list $zmax1 $zmax2]
}

proc curvature { crd prmtop xcom ycom zcom avg_bx avg_by start_frame end_frame step sstep rstep thestep ulres } {

	package require math::linearalgebra

	# THIS PROCEDURE DETERMINES THE LOCAL CURVATURE OF THE UNDERLYING LIPID BIILAYER CENTERED AROUND THE GEOMETRIC CENTRE OF THE PROTEIN

	set f1 [open "local_curvature.txt" "w"]

	puts $f1 "#radius	theta	curvature_UL	curvature_LL"

	set bln [list 0 0 1]
	
	for { set r $sstep } { $r < 200.0 } { set r [expr { $r + $rstep }] } {
		for {set theta 0.0} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
			set angleul($r,$theta) 0.0
			set anglell($r,$theta) 0.0
			set nterm1($r,$theta) 0
			set nterm2($r,$theta) 0
			for {set i $start_frame} {$i < $end_frame} {incr i $step} {
				set sdangleul($i,$r,$theta) -1.0
				set sdanglell($i,$r,$theta) -1.0
			}
		}
	}
	if { $avg_bx < $avg_by } {
		set tsteps [expr { $avg_bx / 2.1 }]
		set tsteps [expr { round($tsteps/$rstep) }]
		set tsteps [expr { $rstep * $tsteps }]
	} else { 
		set tsteps [expr { $avg_by / 2.1 }]
		set tsteps [expr { round($tsteps/$rstep) }]
		set tsteps [expr { $rstep * $tsteps }]
	}

	for {set i $start_frame} {$i < $end_frame} {incr i $step} {

		puts "				#### FRAME $i ####"
	
		# EXECUTING CPPTRAJ

		set f [open "input" "w"]

		puts $f "trajin $crd $i $i"
		puts $f "trajout output.pdb pdb"
		puts $f "trajout output.rst rst"
		puts $f "go"

		close $f

		exec cpptraj -p $prmtop -i input

		set g [open "output.pdb" "r"]
		set data [read $g]
		close $g

		set g1 [open "output.rst" "r"]
		set data1 [read $g1]
		close $g1	

		set k 0

		while { $k < [llength $data1] } {
			incr k
		}

		set box_x [lindex $data1 [expr { $k - 3 }]]
		set box_y [lindex $data1 [expr { $k - 2 }]]
		set box_z [lindex $data1 [expr { $k - 1 }]]

		set j 0

		for { set r $sstep } { $r <= $tsteps } { set r [expr { $r + $rstep }] } {
			puts "				#### RADIUS OF $r ang FROM THE PROTEIN GEOMETRIC CENTRE ####"
			set ij 0
			for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
				set xref [expr { cos(($theta * 3.14) / 180.0) }]
				set yref [expr { sin(($theta * 3.14) / 180.0) }]
				set refvec [list $xref $yref]
				set refvec [::math::linearalgebra::unitLengthVector $refvec]
				set coord1 ""
				set coord2 ""
				set reslist ""
				set reslist1 ""
				set reslist2 ""
				#set reslist ""
				set ind_ul 0
				set ind_ll 0
				set ind1 0
				set ind2 0
				set rr($j) $r
				set the($j,$ij) $theta 
				set k 0
				while { $k < [llength $data] } {
					if { [lindex $data $k] == "ATOM" && [lindex $data [expr { $k + 2 }]] == "P31" } {
						set x1 [lindex $data [expr { $k + 5 }]]
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
								set y2 $yori
								set t 0
								while { [string range $y2 $t $t] != "." } {
									incr t
								}
								set yori [string range $y2 0 [expr { $t + 3 }]]
								set zori [string range $y2 [expr { $t + 4 }] end]
							} else {
								set zori [lindex $data [expr { $k + 6 }]] 
							}		
						} else { 
							set xori $x1
							set y1 [lindex $data [expr { $k + 6 }]]
							set sy1 [string length $y1]
							if { $sy1 > 8 } {
								set t 0
								while { [string range $y1 $t $t] != "." } {
									incr t
								}
								set yori [string range $y1 0 [expr { $t + 3 }]]
								set zori [string range $y1 [expr { $t + 4 }] end]
							} else {
								set yori [lindex $data [expr { $k + 6 }]]
								set zori [lindex $data [expr { $k + 7 }]]
							}
						}
					
						set x1 $xori
						set y1 $yori
						set z1 $zori

						for { set imx [expr { -1 * $box_x }] } {$imx <= $box_x} {set imx [expr { $imx + $box_x }]} {
							for { set imy [expr { -1 * $box_y }] } {$imy <= $box_y} {set imy [expr { $imy + $box_y }]} {
								set x [expr { $x1 + $imx }]
								set y [expr { $y1 + $imy }]

								set xc [expr { $x - $xcom }]
								set yc [expr { $y - $ycom }]
								set tan [expr { $yc / $xc }]

								set tan [expr { atan($tan) }]
								set tan [expr { ($tan * 180.0) / 3.14 }]

								# 2ND QUARDRANT

								if { $yc > 0 && $xc < 0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 3RD QUARDRANT
	 
								if { $xc < 0.0 && $yc < 0.0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 4TH QUARDRANT

								if { $xc > 0.0 && $yc < 0.0 } {
									set tan [expr { $tan + 360.0 }]
								}
					
								set xc [expr { $xc * $xc }]	
								set yc [expr { $yc * $yc }]
								set zc [expr { $z1 - $zcom }]
								set zc [expr { $zc * $zc }]

								set rc [expr { sqrt($xc + $yc) }]

								if { $rc >= [expr { $r - $rstep }] && $rc <= $r } {
									if { $tan >= [expr { $theta - $thestep }] && $tan <= $theta } {	
										if { [lindex $data [expr { $k + 4 }]] <= $ulres } {
											set reslist1 [linsert $reslist1 $ind_ul [lindex $data [expr { $k + 4 }]]]
											set coord1 [linsert $coord1 $ind1 [expr { $x1 - $xcom }]]
											set coord1 [linsert $coord1 [expr { $ind1 + 1 }] [expr { $y1 - $ycom }]]
											set coord1 [linsert $coord1 [expr { $ind1 + 2 }] [expr { $z1 - $zcom }]]											
											incr ind1 3
											incr ind_ul
										} else {
											set reslist2 [linsert $reslist2 $ind_ll [lindex $data [expr { $k + 4 }]]]
											set coord2 [linsert $coord2 $ind2 [expr { $x1 - $xcom }]]
											set coord2 [linsert $coord2 [expr { $ind2 + 1 }] [expr { $y1 - $ycom }]]
											set coord2 [linsert $coord2 [expr { $ind2 + 2 }] [expr { $z1 - $zcom }]]																
											incr ind2 3
											incr ind_ll
										}
									}
								}
							}
						}
					}
					incr k
				}

				# DETERMINING THE CURVATURE OF THE LIPID BILAYER 

				# UPPPER LEAFLET

				if { [llength $reslist1] > 2 } {
					set lip1 [expr { rand() * [llength $reslist1] }]
					set lip1 [format "%.0f" $lip1]
					while { $lip1 == 0 } {
						set lip1 [expr { rand() * [llength $reslist1] }]
						set lip1 [format "%.0f" $lip1]
					}
					set lip2 [expr { rand() * [llength $reslist1] }]
					set lip2 [format "%.0f" $lip2]
					while { $lip2 == $lip1 || $lip2 == 0 } {	 
						set lip2 [expr { rand() * [llength $reslist1] }]
						set lip2 [format "%.0f" $lip2]
					}

					set lip3 [expr { rand() * [llength $reslist1] }]
					set lip3 [format "%.0f" $lip3]
					while { $lip3 == $lip1 || $lip3 == $lip2 || $lip3 == 0 } {	 
						set lip3 [expr { rand() * [llength $reslist1] }]
						set lip3 [format "%.0f" $lip3]
					}

					set xl1 [lindex $coord1 [expr { ($lip1 - 1) * 3 }]]
					set yl1 [lindex $coord1 [expr { (($lip1 - 1) * 3) + 1 }]]
					set zl1 [lindex $coord1 [expr { (($lip1 - 1) * 3) + 2 }]]
					set rl [expr { ($xl1*$xl1) + ($yl1*$yl1) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl1/$xl1) }]
					set vec12d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec12d [::math::linearalgebra::unitLengthVector $vec12d]

					set xl2 [lindex $coord1 [expr { ($lip2 - 1) * 3 }]]
					set yl2 [lindex $coord1 [expr { (($lip2 - 1) * 3) + 1 }]]
					set zl2 [lindex $coord1 [expr { (($lip2 - 1) * 3) + 2 }]]
					set rl [expr { ($xl2*$xl2) + ($yl2*$yl2) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl2/$xl2) }]
					set vec22d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec22d [::math::linearalgebra::unitLengthVector $vec22d]

					set xl3 [lindex $coord1 [expr { ($lip3 - 1) * 3 }]]
					set yl3 [lindex $coord1 [expr { (($lip3 - 1) * 3) + 1 }]]
					set zl3 [lindex $coord1 [expr { (($lip3 - 1) * 3) + 2 }]]
					set rl [expr { ($xl3*$xl3) + ($yl3*$yl3) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl3/$xl3) }]
					set vec32d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec32d [::math::linearalgebra::unitLengthVector $vec32d]

					set theta1 [::math::linearalgebra::dotproduct $refvec $uvec12d]
					set theta2 [::math::linearalgebra::dotproduct $refvec $uvec22d]
					set theta3 [::math::linearalgebra::dotproduct $refvec $uvec32d]
				
					if { $theta1 < $theta2 && $theta1 < $theta3 } {
						set coordl1 [list $xl1 $yl1 $zl1]
						if { $theta2 < $theta3 } {
							set coordl2 [list $xl2 $yl2 $zl2]
							set coordl3 [list $xl3 $yl3 $zl3]
						} else {
							set coordl3 [list $xl2 $yl2 $zl2]
							set coordl2 [list $xl3 $yl3 $zl3]	
						}
					}
					
					if { $theta2 < $theta1 && $theta2 < $theta3 } {
						set coordl1 [list $xl2 $yl2 $zl2]
						if { $theta1 < $theta3 } {
							set coordl2 [list $xl1 $yl1 $zl1]
							set coordl3 [list $xl3 $yl3 $zl3]
						} else {
							set coordl3 [list $xl1 $yl1 $zl1]
							set coordl2 [list $xl3 $yl3 $zl3]	
						}
					}
					
					if { $theta3 < $theta1 && $theta3 < $theta2 } {
						set coordl1 [list $xl3 $yl3 $zl3]
						if { $theta1 < $theta2 } {
							set coordl2 [list $xl1 $yl1 $zl1]
							set coordl3 [list $xl2 $yl2 $zl2]
						} else {
							set coordl3 [list $xl1 $yl1 $zl1]
							set coordl2 [list $xl2 $yl2 $zl2]	
						}
					}

					set vec1 [::math::linearalgebra::sub $coordl2 $coordl1]
					set vec2 [::math::linearalgebra::sub $coordl3 $coordl2]
	
					set vec1 [::math::linearalgebra::unitLengthVector $vec1]
					set vec2 [::math::linearalgebra::unitLengthVector $vec2]					

					set normal [::math::linearalgebra::crossproduct $vec2 $vec1]
					set normal [::math::linearalgebra::unitLengthVector $normal]

					set angle [::math::linearalgebra::dotproduct $bln $normal]
					set angle [expr { (acos($angle)) * 180.0 / 3.14 }]
				
					set sdangleul($i,$r,$theta) $angle

					set angleul($r,$theta) [expr { $angle + $angleul($r,$theta) }]
					set nterm1($r,$theta) [expr { $nterm1($r,$theta) + 1 }]
				}

				# LOWER LEAFLET

				if { [llength $reslist2] > 2 } {
					set lip1 [expr { rand() * [llength $reslist2] }]
					set lip1 [format "%.0f" $lip1]
					while { $lip1 == 0 } {
						set lip1 [expr { rand() * [llength $reslist2] }]
						set lip1 [format "%.0f" $lip1]
					}
					set lip2 [expr { rand() * [llength $reslist2] }]
					set lip2 [format "%.0f" $lip2]
					while { $lip2 == $lip1 || $lip2 == 0 } {	 
						set lip2 [expr { rand() * [llength $reslist2] }]
						set lip2 [format "%.0f" $lip2]
					}

					set lip3 [expr { rand() * [llength $reslist2] }]
					set lip3 [format "%.0f" $lip3]
					while { $lip3 == $lip1 || $lip3 == $lip2 || $lip3 == 0 } {	 
						set lip3 [expr { rand() * [llength $reslist2] }]
						set lip3 [format "%.0f" $lip3]
					}

					set xl1 [lindex $coord2 [expr { ($lip1 - 1) * 3 }]]
					set yl1 [lindex $coord2 [expr { (($lip1 - 1) * 3) + 1 }]]
					set zl1 [lindex $coord2 [expr { (($lip1 - 1) * 3) + 2 }]]
					set rl [expr { ($xl1*$xl1) + ($yl1*$yl1) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl1/$xl1) }]
					set vec12d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec12d [::math::linearalgebra::unitLengthVector $vec12d]

					set xl2 [lindex $coord2 [expr { ($lip2 - 1) * 3 }]]
					set yl2 [lindex $coord2 [expr { (($lip2 - 1) * 3) + 1 }]]
					set zl2 [lindex $coord2 [expr { (($lip2 - 1) * 3) + 2 }]]
					set rl [expr { ($xl2*$xl2) + ($yl2*$yl2) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl2/$xl2) }]
					set vec22d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec22d [::math::linearalgebra::unitLengthVector $vec22d]

					set xl3 [lindex $coord2 [expr { ($lip3 - 1) * 3 }]]
					set yl3 [lindex $coord2 [expr { (($lip3 - 1) * 3) + 1 }]]
					set zl3 [lindex $coord2 [expr { (($lip3 - 1) * 3) + 2 }]]
					set rl [expr { ($xl3*$xl3) + ($yl3*$yl3) }]	
					set rl [expr { sqrt($rl) }]
					set thetarl [expr { atan($yl3/$xl3) }]
					set vec32d [list [expr { $rl*cos($thetarl) }] [expr { $rl*sin($thetarl) }]]
					set uvec32d [::math::linearalgebra::unitLengthVector $vec32d]

					set theta1 [::math::linearalgebra::dotproduct $refvec $uvec12d]
					set theta2 [::math::linearalgebra::dotproduct $refvec $uvec22d]
					set theta3 [::math::linearalgebra::dotproduct $refvec $uvec32d]
				
					if { $theta1 < $theta2 && $theta1 < $theta3 } {
						set coordl1 [list $xl1 $yl1 $zl1]
						if { $theta2 < $theta3 } {
							set coordl2 [list $xl2 $yl2 $zl2]
							set coordl3 [list $xl3 $yl3 $zl3]
						} else {
							set coordl3 [list $xl2 $yl2 $zl2]
							set coordl2 [list $xl3 $yl3 $zl3]	
						}
					}
					
					if { $theta2 < $theta1 && $theta2 < $theta3 } {
						set coordl1 [list $xl2 $yl2 $zl2]
						if { $theta1 < $theta3 } {
							set coordl2 [list $xl1 $yl1 $zl1]
							set coordl3 [list $xl3 $yl3 $zl3]
						} else {
							set coordl3 [list $xl1 $yl1 $zl1]
							set coordl2 [list $xl3 $yl3 $zl3]	
						}
					}
					
					if { $theta3 < $theta1 && $theta3 < $theta2 } {
						set coordl1 [list $xl3 $yl3 $zl3]
						if { $theta1 < $theta2 } {
							set coordl2 [list $xl1 $yl1 $zl1]
							set coordl3 [list $xl2 $yl2 $zl2]
						} else {
							set coordl3 [list $xl1 $yl1 $zl1]
							set coordl2 [list $xl2 $yl2 $zl2]	
						}
					}

					set vec1 [::math::linearalgebra::sub $coordl2 $coordl1]
					set vec2 [::math::linearalgebra::sub $coordl3 $coordl2]
	
					set vec1 [::math::linearalgebra::unitLengthVector $vec1]
					set vec2 [::math::linearalgebra::unitLengthVector $vec2]					

					set normal [::math::linearalgebra::crossproduct $vec2 $vec1]
					set normal [::math::linearalgebra::unitLengthVector $normal]

					set angle [::math::linearalgebra::dotproduct $bln $normal]
					set angle [expr { (acos($angle)) * 180.0 / 3.14 }]

					set sdanglell($i,$r,$theta) $angle

					set anglell($r,$theta) [expr { $angle + $anglell($r,$theta) }]
					set nterm2($r,$theta) [expr { $nterm2($r,$theta) + 1 }]
				}
			}
		}
	}
	for { set r $sstep } { $r <= $tsteps } { set r [expr { $r + $rstep }] } {
		for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
			if { $nterm1($r,$theta) > 0 } {
				set angleul($r,$theta) [expr { $angleul($r,$theta) / $nterm1($r,$theta) }]
			} else {
				set angleul($r,$theta) "--"
			}
			if { $nterm2($r,$theta) > 0 } {
				set anglell($r,$theta) [expr { $anglell($r,$theta) / $nterm2($r,$theta) }]
			} else {
				set anglell($r,$theta) "-"
			}
			puts $f1 "$r	$theta	$angleul($r,$theta)	$anglell($r,$theta)"
		}
	}
	close $f1

	# FOR A VIDEO

	exec mkdir -p video_files_cur

	set vi 0
	for {set i $start_frame} {$i < $end_frame} {incr i $step} {
		set v [open "$vi.txt" "w"]
		for { set r $sstep } { $r <= $tsteps } { set r [expr { $r + $rstep }] } {
			for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
				puts $v "$r	$theta	$sdangleul($i,$r,$theta)	$sdanglell($i,$r,$theta)"
			}
		}
		close $v
		exec mv $vi.txt ./video_files_cur
		incr vi
	}
							
}

proc flux {prmtop crd ll ul start end step ionn name} {
	
	# FOR THE CACULATION OF FLUX THROUGH A PROTEIN CHANNEL

	# FOLLOWING THE TRAJECTORY OF ALL THE WATER MOLECULES BETWEEN THE ul AND ll

	set w [open "$name.txt" "w"]
	set w1 [open "$name-ent.txt" "w"]
	
	puts $w "#Frame number_of_water/ions"
	puts $w1 "#Frame number_of_water/ions"

	set nframes [expr { ($end - $start) / $step }]

	set n 0
	set om 0
	set nm 0
	set nt 0

	for {set i $start } {$i < $end} {incr i $step} {	

		set nwat 0
		set nwate 0

		puts "			**** FRAME $i ****"

		# EXECUTING CPPTRAJ TO GET THE PDB PER STEP 

		set inp [open "input" "w"]

		puts $inp "trajin $crd $i $i"
		puts $inp "strip !:$ionn"
		puts $inp "trajout output.pdb"
		puts $inp "trajout output.rst"
		puts $inp "go"

		close $inp

		exec cpptraj -p $prmtop -i input

		set f [open "output.pdb" "r"]
		set data [read $f]
		close $f

		set g1 [open "output.rst" "r"]
		set data1 [read $g1]
		close $g1

		set box_x [lindex $data1 [expr { [llength $data1] - 3 }]]
		set box_y [lindex $data1 [expr { [llength $data1] - 2 }]]
		set box_z [lindex $data1 [expr { [llength $data1] - 1 }]]

	# ******************************************************************************************************************************************************************************#

		# CONTROLLING PARAMETERS 

		set rad 10.0
		set zdim 8.0
		set shiftx -3.0
		set shifty 6.0

	# ******************************************************************************************************************************************************************************#

		set k 0
		while { $k < [llength $data] } {
			if { [lindex $data $k] == "ATOM" && [lindex $data [expr { $k + 2 }]] == "$ionn" } {
				set x1 [lindex $data [expr { $k + 5 }]]
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
						set y2 $yori
						set t 0
						while { [string range $y2 $t $t] != "." } {
							incr t
						}
						set yori [string range $y2 0 [expr { $t + 3 }]]
						set zori [string range $y2 [expr { $t + 4 }] end]
					} else {
						set zori [lindex $data [expr { $k + 6 }]] 
					}		
				} else { 
					set xori $x1
					set y1 [lindex $data [expr { $k + 6 }]]
					set sy1 [string length $y1]
					if { $sy1 > 8 } {
						set t 0
						while { [string range $y1 $t $t] != "." } {
							incr t
						}
						set yori [string range $y1 0 [expr { $t + 3 }]]
						set zori [string range $y1 [expr { $t + 4 }] end]
					} else {
						set yori [lindex $data [expr { $k + 6 }]]
						set zori [lindex $data [expr { $k + 7 }]]
					}
				}
				set x $xori
				set y $yori
				set z $zori

				if { $z > [expr { $ll + $zdim }] && $z < [expr { $ul - $zdim }] } {		

					#  A CUBOID OF DIMENSION ALONG X AND Y OF 2*$rad Ang

					if { $x > [expr { (($box_x / 2.0) + $shiftx) - $rad }] && $x < [expr { (($box_x / 2.0) + $shiftx) + $rad }] && $y > [expr { (($box_y / 2.0) + $shifty) - $rad }] && $y < [expr { (($box_y / 2.0) + $shifty) + (1*$rad) }]} { 
						incr nwat
						set resn [lindex $data [expr { $k + 4 }]]
						set atomid [lindex $data [expr { $k + 1}]]
						set aidt($nt) $atomid
						incr nt
						set count 0
						for {set j 0} {$j < $n} {incr j} {
							if { $atomid == $aid($j) } {
								incr count
							}
						}
						set count1 0
						for {set j $om} {$j < $nm} {incr j} {
							if { $atomid == $aidt($j) } {
								incr count1
							}
						}
						if { $count1 == 0 } {
							incr nwate
						}
						if { $count == 0 } {
							set res($n) $resn
							set aid($n) $atomid
							incr n
						}
					}
				}
			}
			incr k
		}
		set om $nm
		set nm $nt
		puts $w "$i $nwat"
		puts $w1 "$i $nwate"
	}
	close $w
	close $w1

	flux2 $ul $ll $name
}

proc flux2 { ul ll name } {

	# CALCULATING THE NUMBER OF WATER LEAVING THE AREA

	set f [open "$name.txt" "r"]
	set data [read $f]
	close $f

	set g [open "$name-ent.txt" "r"]
	set data1 [read $g]
	close $g

	set h [open "$name-lea.txt" "w"]
	puts $h "#Frame number_of_water/ions"

	set k 2

	set nwato 0

	while { $k < [llength $data] } {
		set t [lindex $data $k]
		set nwat [lindex $data [expr { $k + 1 }]]
		set nwate [lindex $data1 [expr { $k + 1 }]]

		set nwatl [expr { $nwato + $nwate - $nwat }]

		puts $h "$t $nwatl"

		set nwato $nwat

		incr k 2
	}
	close $h

	set t [format "%.1f" [expr { $ul - $ll - 20.0 }]]

	puts ""
	puts ""

	puts "		|"
	puts "		| nwate"
	puts "		v"
	puts "	   -----------"
	puts "	   | nwat     |"
	puts "	   | $t Ang |"
	puts "	   -----------"
	puts "		|"
	puts "		| nwatl"
	puts "		v"
}	

proc protein_spread_1D {inp} {

	# 1 DIMENSIONAL PROTEIN PROFILE (z)	

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

proc protein_profile_1D { inp dir} {

	# 1 DIMENSIONAL PROTEIN PROFILE (z)	

	protein_spread_1D $inp
	
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

	set f1 [open "protein_distribution.txt" "w"]
	puts $f1 "#Z	X	Y"
	
	set f [open "$inp" "r"]
	set data1 [read $f]
	close $f
	
	set step 1.0

	set f2 [open "amino_acid_dis.txt" "w"]
	puts $f2 "#Z	VAL	LEU	ILE	MET	PRO	PHE	ARG	LYS	ASP	ASN	GLU	GLN	SER	THR	HIS/HIE/HID	TYR	TRP	GLY	ALA	CYS	"
	
	set f3 [open "amino_acid_group_dis.txt" "w"]
	puts $f3 "#Z	gp1	gp2	gp3	gp4	gp5"

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
					set cory [string range $x1 [expr { $t + 3 }] end]
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
						set corz [string range $y1 [expr { $t + 3 }] end]
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
		set var1 [expr { $maxx - $minx }]
		set var1 [format "%.2f" $var1]

		set var2 [expr { $maxy - $miny }]
		set var2 [format "%.2f" $var2]

		puts $f1 " $ij	$var1	$var2"

		puts $f2 " [format %.3f $ij]	$ac(0)	$ac(1)	$ac(2)	$ac(3)	$ac(4)	$ac(5)	$ac(6)	$ac(7)	$ac(8)	$ac(9)	$ac(10)	$ac(11)	$ac(12)	$ac(13)	$ac(14)	$ac(15)	$ac(16)	$ac(17)	$ac(18)	$ac(19)"

		puts $f3 " [format %.3f $ij]	$gp(0)	$gp(1)	$gp(2)	$gp(3)	$gp(4)"
	}
	close $f1
	close $f2
	close $f3
	close $f4
}

proc protein_spread {inp} {	

	# 2D PROTEIN PROFILE (z,x);(z,y)

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

			#set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
			#set schain_id [string length $chain_id]
			#if { $schain_id > 1 } {
				set shift1 1
			#} else { 
				#set shift1 0
			#}

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

proc protein_profile { inp dir} {

	# 2D PROTEIN PROFILE

	protein_spread $inp
	
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
	set f1 [open "protein_distribution.txt" "w"]
	puts $f1 "#Z X	Y"
	
	set f [open "$inp" "r"]
	set data1 [read $f]
	close $f
	
	set step 1.0

	for {set ij $min} {$ij < $max} {set ij [expr { $ij + $step }]} {

		set oldresnum -1

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
				
				#set chain_id [lindex $data1 [expr { $k + 4 - $shift - $shift2}]]
				#set schain_id [string length $chain_id]
				#if { $schain_id > 1 } {
					set shift1 1
				#} else { 
					#set shift1 0
				#}


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

		set var1 [expr { $maxx - $minx }]
		set var1 [format "%.2f" $var1]

		set var2 [expr { $maxy - $miny }]
		set var2 [format "%.2f" $var2]

		puts $f1 " $ij	$var1	$var2"
	}
	close $f1
}

proc pp {prmtop crd xori yori zori start_frame end_frame step start_res end_res sstep estep rstep thestep lp31 up31 nzstep} {

	# 3D PROTEIN PROFILE

	set lp31 [format "%.4f" $lp31]
	set up31 [format "%.4f" $up31]

	set zstep [expr { ($up31 - $lp31) / $nzstep }]
	set zstep [format "%.4f" $zstep]

	if { $nzstep == 1 } {
		set zstart $up31
	} else {
		set zstart [expr { $lp31 + $zstep }]
	}

	# SLIGHT INCTREASE IN up31 VALUE TO TAKE OF ROUNDING OFF

	set up31 [expr { $up31 + 0.5 }]

	set nframe [expr { ($end_frame - $start_frame) / $step }]

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

	set f3 [open "amino_acid_group_dis_3D.txt" "w"]
	puts $f3 "#R theta	z	gp1	sdgp1	gp2	sdgp2	gp3	sdgp3	gp4	sdgp4	gp5	sdgp5"
	
	# AMINO ACID VARIABLES
	puts "				#### DEFINING THE VARIABLES #####"

	for {set i 0} {$i < 20} {incr i} {
		set per [expr { ($i * 100)/ 20.0 }]
		set per [format "%.2f" $per]
		puts "				#### $per PERCENT DONE ####"
		for {set j 0} {$j < 20} {incr j} {
			for {set o 0} {$o < 20} {incr o} {
				for {set k 0} {$k < 20} {incr k} {
					set ac($i,$j,$k,$o) 0.0
					set gp($i,$j,$k,$o) 0.0
					for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
						set gpsd($i,$j,$k,$o,$fr) 0.0
					}
				}
			}
		}
	}
	puts ""

	for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
		puts ""
		puts "				#### FRAME $fr ####"
	
		# EXECUTING CPPTRAJ

		set f [open "input" "w"]

		puts $f "trajin $crd $fr $fr"
		puts $f "strip !:$start_res-$end_res"
		puts $f "trajout output.pdb pdb"
		puts $f "trajout output.rst rst"
		puts $f "go"

		close $f

		exec cpptraj -p $prmtop -i input

		set g [open "output.pdb" "r"]
		set data1 [read $g]
		close $g

		set m 0
			
		for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
			set zd [format "%.4f" $zd]
			set i 0
			puts ""
			puts "				#### LAYER $zd OF THE LEAFLET ####"
			puts ""
			for {set r $sstep} {$r <= $estep } {set r [expr { $r + $rstep }]} {
				puts "				#### RADIUS $r FROM GEOMETRIC CENTRE ####"
				set j 0
				for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set oldresnum -1

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

							set xc [expr { $corx - $xori }]
							set yc [expr { $cory - $yori }]

							set tan [expr { $yc / $xc }]

							set tan [expr { atan($tan) }]
							set tan [expr { ($tan * 180.0) / 3.14 }]

							# 2ND QUARDRANT
							if { $yc > 0 && $xc < 0 } {
								set tan [expr { 180.0 + $tan }]
							}
	
							# 3RD QUARDRANT
	 
							if { $xc < 0.0 && $yc < 0.0 } {
								set tan [expr { 180.0 + $tan }]
							}

							# 4TH QUARDRANT
							if { $xc > 0.0 && $yc < 0.0 } {
								set tan [expr { $tan + 360.0 }]
							}
					
							set xc [expr { $xc * $xc }]	
							set yc [expr { $yc * $yc }]

							set rc [expr { sqrt($xc + $yc) }]	

							if { $corz >= [expr { $zd - $zstep }] && $corz <= $zd } {
								if { $rc >= [expr { $r - $rstep }] && $rc <= $r } {
									if { $tan >= [expr { $theta - $thestep }] && $tan <= $theta } {

										# DETERMING THE RESIDUE DISTRIBUTION

										if { $oldresnum != [lindex $data1 [expr { $k + 5 - $shift - $shift1 -$shift2}]] } { 

											# HYDROPHOBIC RESIDUES :: GROUP 1 :: MEMBRANE CORE

											if { $resn == "VAL" } {
												#set ac(0,$i,$j,$zd) [expr { $ac(0,$i,$j,$zd) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "LEU" } {
												#set ac(1,$i,$j) [expr { $ac(1,$i,$j) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "ILE" } {
												#set ac(2,$i,$j) [expr { $ac(2,$i,$j) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "MET" } {
												#set ac(3,$i,$j) [expr { $ac(3,$i,$j) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "PRO" } {
												#set ac(4,$i,$j) [expr { $ac(4,$i,$j) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "PHE" } {
												#set ac(5,$i,$j) [expr { $ac(5,$i,$j) + 1.0 }]
												set gp(0,$i,$j,$m) [expr { $gp(0,$i,$j,$m) + 1.0 }]
												set gpsd(0,$i,$j,$m,$fr) [expr { $gpsd(0,$i,$j,$m,$fr) + 1.0 }]
											}

											# CHARGED AND LARGE POLAR RESIDUES :: GROUP 2 :: MEMBRANE WATER INTERFACE

											if { $resn == "ARG" } {
												#set ac(6,$i,$j) [expr { $ac(6,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "LYS" } {
												#set ac(7,$i,$j) [expr { $ac(7,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "ASP" } {
												#set ac(8,$i,$j) [expr { $ac(8,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "ASN" } {
												#set ac(9,$i,$j) [expr { $ac(9,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "GLU" } {
												#set ac(10,$i,$j) [expr { $ac(10,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "GLN" } {
												#set ac(11,$i,$j) [expr { $ac(11,$i,$j) + 1.0 }]
												set gp(1,$i,$j,$m) [expr { $gp(1,$i,$j,$m) + 1.0 }]
												set gpsd(1,$i,$j,$m,$fr) [expr { $gpsd(1,$i,$j,$m,$fr) + 1.0 }]
											}

											# SMALL POLAR RESIDUES :: GROUP 3 :: INTERFACE AS WELL AS WITHIN BILAYER

											if { $resn == "SER" } {
												#set ac(12,$i,$j) [expr { $ac(12,$i,$j) + 1.0 }]
												set gp(2,$i,$j,$m) [expr { $gp(2,$i,$j,$m) + 1.0 }]
												set gpsd(2,$i,$j,$m,$fr) [expr { $gpsd(2,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "THR" } {
												#set ac(13,$i,$j) [expr { $ac(13,$i,$j) + 1.0 }]
												set gp(2,$i,$j,$m) [expr { $gp(2,$i,$j,$m) + 1.0 }]
												set gpsd(2,$i,$j,$m,$fr) [expr { $gpsd(2,$i,$j,$m,$fr) + 1.0 }]
											}
						
											# OTHER CHARGED RESIDUES :: GROUP 4  :: BOUNDARY BETWEEN HYDROPHOBIC CORE AND HEADGROUP REGION

											if { $resn == "HIS" || $resn == "HIE" || $resn == "HID" } {
												#set ac(14,$i,$j) [expr { $ac(14,$i,$j) + 1.0 }]
												set gp(3,$i,$j,$m) [expr { $gp(3,$i,$j,$m) + 1.0 }]
												set gpsd(3,$i,$j,$m,$fr) [expr { $gpsd(3,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "TYR" } {
												#set ac(15,$i,$j) [expr { $ac(15,$i,$j) + 1.0 }]
												set gp(3,$i,$j,$m) [expr { $gp(3,$i,$j,$m) + 1.0 }]
												set gpsd(3,$i,$j,$m,$fr) [expr { $gpsd(3,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "TRP" } {
												#set ac(16,$i,$j) [expr { $ac(16,$i,$j) + 1.0 }]
												set gp(3,$i,$j,$m) [expr { $gp(3,$i,$j,$m) + 1.0 }]
												set gpsd(3,$i,$j,$m,$fr) [expr { $gpsd(3,$i,$j,$m,$fr) + 1.0 }]
											}

											# OTHER SMALL RESIDUES :: GROUP 5 :: FOUND THROUGHOUT THE BILAYER

											if { $resn == "GLY" } {
												#set ac(17,$i,$j) [expr { $ac(17,$i,$j) + 1.0 }]
												set gp(4,$i,$j,$m) [expr { $gp(4,$i,$j,$m) + 1.0 }]
												set gpsd(4,$i,$j,$m,$fr) [expr { $gpsd(4,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "ALA" } {
												#set ac(18,$i,$j) [expr { $ac(18,$i,$j) + 1.0 }]
												set gp(4,$i,$j,$m) [expr { $gp(4,$i,$j,$m) + 1.0 }]
												set gpsd(4,$i,$j,$m,$fr) [expr { $gpsd(4,$i,$j,$m,$fr) + 1.0 }]
											}
											if { $resn == "CYS" } {
												#set ac(19,$i,$j) [expr { $ac(19,$i,$j) + 1.0 }]
												set gp(4,$i,$j,$m) [expr { $gp(4,$i,$j,$m) + 1.0 }]
												set gpsd(4,$i,$j,$m,$fr) [expr { $gpsd(4,$i,$j,$m,$fr) + 1.0 }]
											}					
										}
									}	
								} 
							}
							set oldresnum [lindex $data1 [expr { $k + 5 - $shift - $shift1 -$shift2}]] 
						}
						incr k
					}
				#puts $f2 " [format %.3f $ij] $ac(0) $ac(1) $ac(2) $ac(3) $ac(4) $ac(5) $ac(6) $ac(7) $ac(8) $ac(9) $ac(10) $ac(11) $ac(12) $ac(13) $ac(14) $ac(15) $ac(16) $ac(17) $ac(18) $ac(19)"
				incr j
				}	
			incr i
			}
		incr m
		}
	}
	set o 0
	for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
		set m 0
		for {set r $sstep} {$r <= $estep} {set r [expr { $r + $rstep }]} {
			set n 0
			for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
				set gp(0,$m,$n,$o) [expr { $gp(0,$m,$n,$o) / $nframe }]
				set gp(1,$m,$n,$o) [expr { $gp(1,$m,$n,$o) / $nframe }]
				set gp(2,$m,$n,$o) [expr { $gp(2,$m,$n,$o) / $nframe }]
				set gp(3,$m,$n,$o) [expr { $gp(3,$m,$n,$o) / $nframe }]
				set gp(4,$m,$n,$o) [expr { $gp(4,$m,$n,$o) / $nframe }]
				set gp(5,$m,$n,$o) [expr { $gp(5,$m,$n,$o) / $nframe }]
				set gp0sd 0.0
				set gp1sd 0.0
				set gp2sd 0.0
				set gp3sd 0.0
				set gp4sd 0.0
				for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
					set diff [expr { $gp(0,$m,$n,$o) - $gpsd(0,$m,$n,$o,$fr) }]
					set diff [expr { $diff * $diff }]
					set gp0sd [expr { $gp0sd + $diff }]
					set diff [expr { $gp(1,$m,$n,$o) - $gpsd(1,$m,$n,$o,$fr) }]
					set diff [expr { $diff * $diff }]
					set gp1sd [expr { $gp1sd + $diff }]
					set diff [expr { $gp(2,$m,$n,$o) - $gpsd(2,$m,$n,$o,$fr) }]
					set diff [expr { $diff * $diff }]
					set gp2sd [expr { $gp2sd + $diff }]
					set diff [expr { $gp(3,$m,$n,$o) - $gpsd(3,$m,$n,$o,$fr) }]
					set diff [expr { $diff * $diff }]
					set gp3sd [expr { $gp3sd + $diff }]
					set diff [expr { $gp(4,$m,$n,$o) - $gpsd(4,$m,$n,$o,$fr) }]
					set diff [expr { $diff * $diff }]
					set gp4sd [expr { $gp4sd + $diff }]
				}

				set gp0sd [expr { $gp0sd / $nframe }]
				set gp0sd [expr { sqrt($gp0sd) }]
		 		set gp1sd [expr { $gp1sd / $nframe }]
				set gp1sd [expr { sqrt($gp1sd) }]			
				set gp2sd [expr { $gp2sd / $nframe }]
				set gp2sd [expr { sqrt($gp2sd) }]
				set gp3sd [expr { $gp3sd / $nframe }]
				set gp3sd [expr { sqrt($gp3sd) }]
				set gp4sd [expr { $gp4sd / $nframe }]
				set gp4sd [expr { sqrt($gp4sd) }]

				set normfactor [expr { (($r * $r) - (($r-$rstep)*($r-$rstep))) / 100.0 }]

				set gp(0,$m,$n,$o) [format "%.2f" $gp(0,$m,$n,$o)]
				set gp0sd [format "%.2f" $gp0sd]
				set gp(1,$m,$n,$o) [format "%.2f" $gp(1,$m,$n,$o)]
				set gp1sd [format "%.2f" $gp1sd]
				set gp(2,$m,$n,$o) [format "%.2f" $gp(2,$m,$n,$o)]
				set gp2sd [format "%.2f" $gp2sd]
				set gp(3,$m,$n,$o) [format "%.2f" $gp(3,$m,$n,$o)]
				set gp3sd [format "%.2f" $gp3sd]
				set gp(4,$m,$n,$o) [format "%.2f" $gp(4,$m,$n,$o)]
				set gp4sd [format "%.2f" $gp4sd]

				set var1 [expr { $gp(0,$m,$n,$o) / $normfactor }]
				set var1 [format "%.2f" $var1]
				set var1c [expr { $gp0sd / $normfactor }]
				set var1c [format "%.2f" $var1c]

				set var2 [expr { $gp(1,$m,$n,$o) / $normfactor }]
				set var2 [format "%.2f" $var2]
				set var2c [expr { $gp1sd / $normfactor }]
				set var2c [format "%.2f" $var2c]

				set var3 [expr { $gp(2,$m,$n,$o) / $normfactor }]
				set var3 [format "%.2f" $var3]
				set var3c [expr { $gp2sd / $normfactor }]
				set var3c [format "%.2f" $var3c]

				set var4 [expr { $gp(3,$m,$n,$o) / $normfactor }]
				set var4 [format "%.2f" $var4]
				set var4c [expr { $gp3sd / $normfactor }]
				set var4c [format "%.2f" $var4c]

				set var5 [expr { $gp(4,$m,$n,$o) / $normfactor }]
				set var5 [format "%.2f" $var5]
				set var5c [expr { $gp4sd / $normfactor }]
				set var5c [format "%.2f" $var5c]

				set zd1 [format "%.3f" $zd]
				
				puts $f3 "$r	$theta	$zd1	$var1	$var1c	$var2	$var2c	$var3	$var3c	$var4	$var4c	$var5	$var5c"
				incr n
			}
			incr m
		}
		incr o 
	}
	close $f3

	if { $nzstep == 1 } {

		# FOR A VIDEO

		exec mkdir -p video_files

		set vi 0
		for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
			set v [open "$vi" "w"]
			set o 0
			for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
				set m 0
				for {set r $sstep} {$r <= $estep} {set r [expr { $r + $rstep }]} {
					set n 0
					for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
						set normfactor [expr { (($r * $r) - (($r-$rstep)*($r-$rstep))) / 100.0 }]
						puts $v "$r	$theta	$zd	[expr { $gpsd(0,$m,$n,$o,$fr) / $normfactor }]	[expr { $gpsd(1,$m,$n,$o,$fr) / $normfactor }]	[expr { $gpsd(2,$m,$n,$o,$fr) / $normfactor }]	[expr { $gpsd(3,$m,$n,$o,$fr) / $normfactor }]	[expr { $gpsd(4,$m,$n,$o,$fr) / $normfactor }]"
						incr n
					}
					incr m
				}
				incr o
			}
			close $v
			incr vi
		}

		# CONVERTING INTO STANDARD FORMAT

		set numfile $vi
		 
		for {set i 0} {$i < $numfile} {incr i} {
			set f [open "$i" "r"]
			set data [read $f]
			close $f

			set go [open "$i.txt" "w"]
			set ho [open "avg_$i.txt" "w"]

			set k 0

			while { $k < [llength $data] } {
				set t1 [lindex $data $k]
				set t2 [lindex $data [expr { $k + 1 }]]
				set zd [lindex $data [expr { $k + 2 }]]
				set t3 [lindex $data [expr { $k + 3 }]]
				set t4 [lindex $data [expr { $k + 4 }]]
				set t5 [lindex $data [expr { $k + 5 }]]
				set t6 [lindex $data [expr { $k + 6 }]]
				set t7 [lindex $data [expr { $k + 7 }]]
				set ttotal [expr { $t3 + $t4 + $t5 + $t6 + $t7 }]

				puts $go "$t1	$t2	$zd		$t3		$t4		$t5		$t6		$t7"
				puts $ho "$t1 	$t2	$ttotal"

				incr k 8
			}
			close $go
			close $ho
			exec mv $i ./video_files
			exec mv $i.txt ./video_files
			exec mv avg_$i.txt ./video_files
		}
	}	
}	

proc average {sstep estep rstep thestep} {

	# TO AVERAGE OVER THE TOTAL NUMBER OF AMINO ACID ALONG THE BILAYER THICKNESS

	set f [open "amino_acid_group_dis_3D.txt" "r"]
	set data [read $f]
	close $f

	set g [open "amino_acid_group_dis_3D_avg.txt" "w"]

	puts $g "#r	theta	number_AC	SD_number_AC"

	set rrstep $estep

	for {set i 0} {$i < 1000} {incr i} {
		set t($i) 0.0
		set st($i) 0.0
	}
	set i 0
	set k 13
	set zc [lindex $data 15]
	while { $k < [llength $data] } {
		set z [lindex $data [expr { $k + 2 }]]
		if { $z == $zc } {
			set t1 [lindex $data [expr { $k + 3 }]]
			set st1 [lindex $data [expr { $k + 4 }]]
			set t2 [lindex $data [expr { $k + 5 }]]
			set st2 [lindex $data [expr { $k + 6 }]]
			set t3 [lindex $data [expr { $k + 7 }]]
			set st3 [lindex $data [expr { $k + 8 }]]
			set t4 [lindex $data [expr { $k + 9 }]]
			set st4 [lindex $data [expr { $k + 10 }]]
			set t5 [lindex $data [expr { $k + 11 }]]
			set st5 [lindex $data [expr { $k + 12 }]]

			set t($i) [expr { $t($i) + $t1 + $t2 + $t3 + $t4 + $t5 }]
			set st($i) [expr  { $st($i) + $st1 + $st2 + $st3 + $st4 + $st5 }]	
			incr i			
			incr k 13
		} else {
			set zc [lindex $data [expr { $k + 2 }]]
			set i 0
		}
	}
	set i 0
	for {set r $sstep} {$r <= $rrstep} {set r [expr { $r + $rstep }] } {
		for {set the $thestep} {$the <= 360.0} {set the [expr { $the + $thestep }]} {
			set t($i) [format "%.0f" $t($i)]
			set st($i) [format "%.0f" $st($i)]
			puts $g "$r	$the	$t($i)	$st($i)"
			incr i
		}
	}
	close $g
}

proc com_lp { prmtop lipid } { 

	set f [open "$prmtop" "r"]
	set data [read $f]
	close $f

	set g [open "output.pdb" "r"]
	set data1 [read $g]
	close $g

	set h [open "center_of_mass" "w"]

	set k 0

	while { [lindex $data $k] != "MASS" || [lindex $data [expr { $k - 1 }]] != "%FLAG" } {
		incr k
	}
	incr k 2
	set j 0
	while { [lindex $data $k] != "%FLAG" } {
		set mass($j) [lindex $data $k]
		incr j
		incr k
	}

	set k 0
	set j 0
	while { $k < [llength $data1] } {
		set total_M 0.0
		set xcom 0.0
		set ycom 0.0
		set zcom 0.0
		set c1 0
		set c2 0
		set c3 0
		while { [lindex $data1 $k] != "TER" && [lindex $data1 $k] != "END" } {
			if { [lindex $data1 $k] == "ATOM" } {
				set count 0
				set k1 0
				while { $k1 < [llength $lipid] } {
					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid $k1] } {
						incr count
					}
					incr k1
				}
				if { $count != 0 || [lindex $data1 [expr { $k + 3 }]] == "CHL" } {
					if { [lindex $data1 [expr { $k + 3} ]] == [lindex $lipid 0] && $c1 == 0 } {
						set res1 [lindex $data1 [expr { $k + 4 }]]
						incr c1
					}

					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid 1] && $c2 == 0 } {
						set nk $k
						set res2 [lindex $data1 [expr { $k + 4 }]]
						incr c2
					}

					if { [lindex $data1 [expr { $k + 3 }]] == [lindex $lipid 2] && $c3 == 0 } {
						set res3 [lindex $data1 [expr { $k + 4 }]]
						incr c3
					}

					if { [lindex $data1 [expr { $k + 3 }]] == "CHL" } {
						set res1 [lindex $data1 [expr { $k + 4 }]]
						set res2 $res1
						set res3 $res1
					}

					set x1 [lindex $data1 [expr { $k + 5 }]]
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
							set y2 $yori
							set t 0
							while { [string range $y2 $t $t] != "." } {
								incr t
							}
							set yori [string range $y2 0 [expr { $t + 3 }]]
							set zori [string range $y2 [expr { $t + 4 }] end]
						} else {
							set zori [lindex $data1 [expr { $k + 6 }]] 
						}		
					} else { 
						set xori $x1
						set y1 [lindex $data1 [expr { $k + 6 }]]
						set sy1 [string length $y1]
						if { $sy1 > 8 } {
							set t 0
							while { [string range $y1 $t $t] != "." } {
								incr t
							}
							set yori [string range $y1 0 [expr { $t + 3 }]]
							set zori [string range $y1 [expr { $t + 4 }] end]
						} else {
							set yori [lindex $data1 [expr { $k + 6 }]]
							set zori [lindex $data1 [expr { $k + 7 }]]
						}
					}
					
					set x $xori
					set y $yori
					set z $zori

					set xcom [expr { $xcom + ($x*$mass($j)) }]
					set ycom [expr { $ycom + ($y*$mass($j)) }]
					set zcom [expr { $zcom + ($z*$mass($j)) }]
					set total_M [expr { $total_M + $mass($j) }]
	
				}
				incr j
			}
			incr k
		}
		if { $total_M > 0.0 } {
			if { [lindex $data1 [expr { $k + 2 }]] == "CHL" } { 
				puts $h "[expr { $xcom / $total_M }] [expr { $ycom / $total_M }] [expr { $zcom / $total_M }] $res1 $res2 $res3 CHL"
			} else {
				puts $h "[expr { $xcom / $total_M }] [expr { $ycom / $total_M }] [expr { $zcom / $total_M }] $res1 $res2 $res3 [lindex $data1 [expr { $nk + 3 }]]"
			}
		}
		incr k
	}
	close $h
}

proc wat_prof {prmtop crd xcom ycom zcom start_frame end_frame step sstep tsteps rstep thestep lp31 up31 zwpslice} {

	# THIS PROCEDURE WILL CALCULATE THE WATER PROFILE IN 3D CYLINDRICAL COORDINATES STARTING FROM 0 TO AVG_BZ DIVIDED INTO ZWPSLICE

	set nframe [expr { ($end_frame - $start_frame) / $step }]

	set lp31 [format "%.4f" $lp31]
	set up31 [format "%.4f" $up31]

	set zstep [expr { ($up31 - $lp31) / $zwpslice }]
	set zstep [format "%.4f" $zstep]

	if { $zwpslice == 1 } {
		set zstart $up31
	} else {
		set zstart [expr { $lp31 + $zstep }]
		set zstart [format "%.4f" $zstart]
	}
	
	# SLIGHT INCTREASE IN up31 VALUE TO TAKE OF ROUNDING OFF

	set up31 [expr { $up31 + 0.5 }]

	# VARIABLES

	puts ""
	puts "		### CREATING THE VARIABLES ###"

	set vi 0
	for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
		set o 0
		for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
			set m 0
			for {set r $sstep} {$r <= $tsteps} {set r [expr { $r + $rstep }]} {
				set n 0
				for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set nw($vi,$m,$n) 0
					set nwpf($o,$vi,$m,$n) 0
					incr n
				}
				incr m
			}
			incr o
		}
		incr vi
	}

	set i1 0
	for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
		puts ""
		puts "				#### FRAME $fr ####"
	
		# EXECUTING CPPTRAJ

		set f [open "input" "w"]

		puts $f "trajin $crd $fr $fr"
		puts $f "strip !:WAT"
		puts $f "trajout output.pdb pdb"
		puts $f "go"

		close $f

		exec cpptraj -p $prmtop -i input

		set g [open "output.pdb" "r"]
		set data1 [read $g]
		close $g

		set i2 0
			
		for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
			set zd [format "%.4f" $zd]
			puts ""
			puts "				#### LAYER $zd OF THE LEAFLET ####"
			puts ""
			set i3 0
			for {set r $sstep} {$r <= $tsteps } {set r [expr { $r + $rstep }]} {
				puts "				#### RADIUS $r FROM GEOMETRIC CENTRE ####"
				set i4 0
				for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set k 0
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
				
							set resname [lindex $data1 [expr { $k + 3 - $shift - $shift2}]]
							set sresname [string length $resname]

							# TO ACCOUNT FOR NO CHAIN ID

							set shift1 1

							if { $atype == "O" && $resname == "WAT" } {
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

								set xc [expr { $corx - $xcom }]
								set yc [expr { $cory - $ycom }]

								set tan [expr { $yc / $xc }]

								set tan [expr { $yc / $xc }]

								set tan [expr { atan($tan) }]
								set tan [expr { ($tan * 180.0) / 3.14 }]

								# 2ND QUARDRANT
								if { $yc > 0 && $xc < 0 } {
									set tan [expr { 180.0 + $tan }]
								}
	
								# 3RD QUARDRANT
		 
								if { $xc < 0.0 && $yc < 0.0 } {
									set tan [expr { 180.0 + $tan }]
								}

								# 4TH QUARDRANT
								if { $xc > 0.0 && $yc < 0.0 } {
									set tan [expr { $tan + 360.0 }]
								}

								set xc [expr { $xc * $xc }]	
								set yc [expr { $yc * $yc }]

								set rc [expr { sqrt($xc + $yc) }]	

								if { $corz >= [expr { $zd - $zstep }] && $corz <= $zd } {
									if { $rc >= [expr { $r - $rstep }] && $rc <= $r } {
										if { $tan >= [expr { $theta - $thestep }] && $tan <= $theta } {
											incr nw($i2,$i3,$i4)
											incr nwpf($i1,$i2,$i3,$i4)
										}
									}
								}
							}
						}
						incr k
					}
					incr i4
				}
				incr i3
			}
			incr i2
		}
		incr i1
	}

	# SAVING THE FILES PER SLICE AND PER FRAME

	exec mkdir -p video_files_wp

	set vi 0
	for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
		set zdn [format "%.0f" "$zd"]
		set v [open "z_$zdn.txt" "w"]
		set m 0
		for {set r $sstep} {$r <= $tsteps} {set r [expr { $r + $rstep }]} {
			set n 0
			for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
				set nw($vi,$m,$n) [expr { $nw($vi,$m,$n) / $nframe }]
				set normfactor [expr { (($r * $r) - (($r-$rstep)*($r-$rstep))) / 100.0 }]
				puts $v "$r	$theta	$zdn	 [expr { $nw($vi,$m,$n) / $normfactor }]"
				incr n
			}
			incr m
		}
		incr vi
		close $v
		exec mv z_$zdn.txt ./video_files_wp
	}

	# FOR A VIDEO

	set o 0
	for {set fr $start_frame} {$fr < $end_frame} {incr fr $step} {
		set vi 0
		for {set zd $zstart} {$zd <= $up31} {set zd [expr { $zd + $zstep }]} {
			set zdn [format "%.0f" "$zd"]
			set v1 [open "frame_$zdn.$o.txt" "w"]
			set m 0
			for {set r $sstep} {$r <= $tsteps} {set r [expr { $r + $rstep }]} {
				set n 0
				for {set theta $thestep} {$theta <= 360.0} {set theta [expr { $theta + $thestep }]} {
					set normfactor [expr { (($r * $r) - (($r-$rstep)*($r-$rstep))) / 100.0 }]
					puts $v1 "$r	$theta	$zdn	 [expr { $nwpf($o,$vi,$m,$n) / $normfactor }]"
					incr n
				}
				incr m
			}
			close $v1
			exec mv frame_$zdn.$o.txt ./video_files_wp
			incr vi
		}
		incr o
	}
} 	
	
proc file_delete {} {
	file delete dummy
	file delete res_num
	file delete test
	file delete center_of_mass
	file delete output.pdb
	file delete output.rst
	file delete input
}

execution
file_delete















	
	
