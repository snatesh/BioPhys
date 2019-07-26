proc evec dim {
	set vec {}
	for {set d 0} {$d < $dim} {incr d} {
		lappend vec {}
	}
	return $vec
}



source diffusion_new1.tcl
source align.tcl
package require pbctools

proc standardError {Li L} {
  set n [llength $Li]
  set stdv 0
  foreach x $Li {
    set var [expr pow($x-$L,2)]
    set stdv [expr $stdv + $var]
  }
  set stdv [expr sqrt($stdv/($n-1))/sqrt($n)] 
  return $stdv
}

mol new ../waterbox-protein-removed-neutralized.psf waitfor all
mol addfile ../1UBQ-wrap.dcd waitfor all
set molids {}
lappend molids [molinfo top]
align "protein and backbone" [lindex $molids 0]
mol new ../waterbox-protein-removed-neutralized.psf waitfor all
mol addfile ../1UBQ-mini-heatup-equil-concat-NPT-20ns-1ps.dcd waitfor all
lappend molids [molinfo top]
align "protein and backbone" [lindex $molids 1]

set abeg 0
set aend 19999
set wind_size_list "2"
set next_ref_list "10"
#set wind_size 500
##set next_ref 250

set wind_and_ref_list {}
for {set i 0} {$i < [llength $wind_size_list]} {incr i} {
	lappend wind_and_ref_list [lindex $wind_size_list $i]
	lappend wind_and_ref_list [lindex $next_ref_list $i]
	}


# pore radius
set rad 30
set rad2 [expr pow($rad,2)]
# shell thickness and number of shells
set thickness 2.0
set nshells [expr $rad/$thickness]
# remove existing data 
for {set i 1} {$i <= $nshells} {incr i} {
  file delete -force "shell${i}"
}
# loop over selected time origins with next_ref gap between them and for
# wind_size frames to get ensemble average of MSD for each lag time in wind_size_list for each shell
for {set nshell 1} {$nshell <= $nshells} {incr nshell} {
  foreach {j k} $wind_and_ref_list {
      	puts "$j $k"
      	set wind_size $j
     	set next_ref $k
        set gap [expr {$next_ref - $wind_size}]
    # initialize containers for MSD averaged over windows
    set MSDx [evec 2000]
    set MSDy [evec 2000]
    set MSDz [evec 2000]
    # initialize list for MSD at each window (for stdv calc)
    set preMSDx {}
    set preMSDy {}
    set preMSDz {}
    set num_wind 0
    # set radii for shell 1 A thick
    set rIn [expr ($nshell-1)*$thickness]
    set rOut [expr $nshell*$thickness]
    for {set i $abeg} {$i <= $aend - $wind_size}  {incr i $next_ref} {
		# set water wrap traj to top
		mol top [lindex $molids 0]
		# define central core at current frame
		#set selz [atomselect top "protein" frame $i]
		# measure center of core
		#set cent [measure center $selz]
		# delete core sel
		#$selz delete
		# get x and y components of core center
		#set cx [lindex $cent 0]
		#set cy [lindex $cent 1]
		#set cz [lindex $cent 2]
		# set last frame for this window
		set last [expr {$i + $wind_size}]
		# get waters in core based on location in water-wrapped traj and update
		set sel [atomselect top "name OH2 and (within $rOut of protein) and not (within $rIn of protein)" frame $i]
      set selnum [$sel num]
      if {$selnum > 1} { 
	puts "sel num [$sel num]" 
		  # get previously determined core waters in unwrapped traj and update
		  mol top [lindex $molids 1]
		  set sel1 [atomselect top "index [$sel get index]" frame $i]	
		  # write off dcd of core waters from unwrapped traj
		  animate write dcd $i-$last.dcd beg $i end $last waitfor all sel $sel1 [molinfo top] 
		  $sel1 writepsf $i-$last.psf 
		  # load the water traj
		  mol new $i-$last.psf waitfor all
		  mol addfile $i-$last.dcd waitfor all
		  # get diffusion in x, y and z
	lappend preMSDx [diffusion 0 0 $wind_size $wind_size 1]
	lappend preMSDy [diffusion 1 0 $wind_size $wind_size 1]
	lappend preMSDz [diffusion 2 0 $wind_size $wind_size 1]
		  lset MSDx $num_wind [lindex [lindex $preMSDx end] end]
		  lset MSDy $num_wind [lindex [lindex $preMSDy end] end]
		  lset MSDz $num_wind [lindex [lindex $preMSDz end] end]
	set num_wind [expr $num_wind + 1]
		  # cleanup
		  $sel1 delete
		  mol delete top
		  file delete $i-$last.dcd
	file delete $i-$last.psf
	  }
		$sel delete
      }
      #if {$num_wind > 1} {
#	set MSDx [vecscale $MSDx [expr {1.0/$num_wind}]]
#	set MSDy [vecscale $MSDy [expr {1.0/$num_wind}]]
#	set MSDz [vecscale $MSDz [expr {1.0/$num_wind}]]
	# put preMSD at each lag time into own list
	# and calculate standard error
#	set stdvX {} 
#	set stdvY {}
#	set stdvZ {}
#	for {set i 0} {$i < $wind_size} {incr i} {
#	  set Lx {}
#	  set Ly {}
#	  set Lz {}
#	  for {set j 0} {$j < $num_wind} {incr j} {
#	    lappend Lx [lindex [lindex $preMSDx $j] $i]
#	    lappend Ly [lindex [lindex $preMSDy $j] $i]
#	    lappend Lz [lindex [lindex $preMSDz $j] $i]
#	  }
#	  lappend stdvX [standardError $Lx [lindex $MSDx $i]]
#	  lappend stdvY [standardError $Ly [lindex $MSDy $i]]
#	  lappend stdvZ [standardError $Lz [lindex $MSDz $i]]
#	}
	# write the msd and standard errors for each axis to folder/file
#	file mkdir "shell${nshell}/tau${wind_size}f${gap}g"
#	set msdXfile "shell${nshell}/tau${wind_size}f${gap}g/unwrap-diff-x_out.txt"
#	set msdYfile "shell${nshell}/tau${wind_size}f${gap}g/unwrap-diff-y_out.txt"
#	set msdZfile "shell${nshell}/tau${wind_size}f${gap}g/unwrap-diff-z_out.txt"
#	set stdvXfile "shell${nshell}/tau${wind_size}f${gap}g/stdev-x_out.txt"
#	set stdvYfile "shell${nshell}/tau${wind_size}f${gap}g/stdev-y_out.txt"
#	set stdvZfile "shell${nshell}/tau${wind_size}f${gap}g/stdev-z_out.txt"
#	set outMSDx [open $msdXfile w]
#	set outMSDy [open $msdYfile w]
#	set outMSDz [open $msdZfile w]
#	set outStdvX [open $stdvXfile w]
#	set outStdvY [open $stdvYfile w]
#	set outStdvZ [open $stdvZfile w]
#	for {set i 0} {$i < $wind_size} {incr i} {
#	      puts $outMSDx [format "%.4f" [lindex $MSDx $i]]
#	      puts $outMSDy [format "%.4f" [lindex $MSDy $i]]
#	      puts $outMSDz [format "%.4f" [lindex $MSDz $i]]
#	      puts $outStdvX [format "%.4f" [lindex $stdvX $i]]
#	      puts $outStdvY [format "%.4f" [lindex $stdvY $i]]
#	      puts $outStdvZ [format "%.4f" [lindex $stdvZ $i]]
#	}
#	close $outMSDx
#	close $outMSDy
#	close $outMSDz
#	close $outStdvX
#	close $outStdvY
#	close $outStdvZ
 #     }
  #  }
#  }

set outMSDx [open "testMSDx.txt" w]
set outMSDy [open "testMSDy.txt" w]
set outMSDz [open "testMSDz.txt" w]

for {set i 0} {$i < 2000} {incr i} {
	puts $outMSDx "[expr 10 * $i],[lindex $MSDx $i]"
	puts $outMSDy "[expr 10 * $i],[lindex $MSDy $i]"
	puts $outMSDz "[expr 10 * $i],[lindex $MSDz $i]"
}

close $outMSDx
close $outMSDy
close $outMSDz
