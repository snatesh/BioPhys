source /home/sachin/Research/BioPhys/diffusion/codes/diffusion_new1.tcl
source /home/sachin/Research/BioPhys/diffusion/codes/align.tcl
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

mol new /home/sachin/Research/BioPhys/diffusion/data/2lmp/2lmp-inf-ion.psf waitfor all
mol addfile /home/sachin/Research/BioPhys/diffusion/data/2lmp/2lmp-inf-NPT5-21-stride100-water_around_prot.dcd waitfor all
set molids {}
lappend molids [molinfo top]
align "protein and backbone" [lindex $molids 0]
mol new /home/sachin/Research/BioPhys/diffusion/data/2lmp/2lmp-inf-ion.psf waitfor all
mol addfile /home/sachin/Research/BioPhys/diffusion/data/2lmp/2lmp-inf-NPT5-21-stride100.dcd waitfor all
lappend molids [molinfo top]
align "protein and backbone" [lindex $molids 1]

set abeg 1
set aend 2501
set wind_size_list "10 20 40 80 160 320"
set next_ref_list "15 25 45 85 165 325"
#set wind_size 500
#set next_ref 250

set wind_and_ref_list {}
for {set i 0} {$i < [llength $wind_size_list]} {incr i} {
	lappend wind_and_ref_list [lindex $wind_size_list $i]
	lappend wind_and_ref_list [lindex $next_ref_list $i]
}

# loop over selected time origins with next_ref gap between them and for
# wind_size frames to get ensemble average of MSD for each lag time in wind_size_list
foreach {j k} $wind_and_ref_list {
	puts "$j $k"
	set wind_size $j
	set next_ref $k
  # initialize containers for MSD averaged over windows
  set MSDx [veccreate $wind_size]
  set MSDy [veccreate $wind_size]
  set MSDz [veccreate $wind_size]
  # initialize list for MSD at each window (for stdv calc)
  set preMSDx {}
  set preMSDy {}
  set preMSDz {}
	set num_wind 0
  for {set i $abeg} {$i <= $aend - $wind_size}  {incr i $next_ref} {
		# set water wrap traj to top
		mol top [lindex $molids 0]
		# define central core at current frame
		set selz [atomselect top "protein and resid 33 35 37" frame $i]
		# measure center of core
		set cent [measure center $selz]
		# delete core sel
		$selz delete
		# get x and y components of core center
		set cx [lindex $cent 0]
		set cy [lindex $cent 1]
		# set last frame for this window
		set last [expr {$i + $wind_size}]
		# get waters in core based on location in water-wrapped traj and update
		set sel [atomselect top "name OH2 and sqr(x-$cx) + sqr(y-$cy) <= 100" frame $i]
		puts "[$sel num]" 
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
		set MSDx [vecadd $MSDx [lindex $preMSDx end]]
		set MSDy [vecadd $MSDy [lindex $preMSDy end]]
		set MSDz [vecadd $MSDz [lindex $preMSDz end]]
		# cleanup
		$sel delete
		$sel1 delete
		mol delete top
		file delete $i-$last.dcd
    file delete $i-$last.psf
    set num_wind [expr $num_wind + 1]
	}
  set MSDx [vecscale $MSDx [expr {1.0/$num_wind}]]
  set MSDy [vecscale $MSDy [expr {1.0/$num_wind}]]
  set MSDz [vecscale $MSDz [expr {1.0/$num_wind}]]
  # put preMSD at each lag time into own list
  # and calculate standard error
  set stdvX {} 
  set stdvY {}
  set stdvZ {}
  for {set i 0} {$i < $wind_size} {incr i} {
    set Lx {}
    set Ly {}
    set Lz {}
    for {set j 0} {$j < $num_wind} {incr j} {
      lappend Lx [lindex [lindex $preMSDx $j] $i]
      lappend Ly [lindex [lindex $preMSDy $j] $i]
      lappend Lz [lindex [lindex $preMSDz $j] $i]
    }
    lappend stdvX [standardError $Lx [lindex $MSDx $i]]
    lappend stdvY [standardError $Ly [lindex $MSDy $i]]
    lappend stdvZ [standardError $Lz [lindex $MSDz $i]]
  }
  # write the msd and standard errors for each axis to folder/file
	file mkdir "tau${wind_size}f"
  set msdXfile "tau${wind_size}f/unwrap-diff-x_out.txt"
  set msdYfile "tau${wind_size}f/unwrap-diff-y_out.txt"
  set msdZfile "tau${wind_size}f/unwrap-diff-z_out.txt"
  set stdvXfile "tau${wind_size}f/stdev-x_out.txt"
  set stdvYfile "tau${wind_size}f/stdev-y_out.txt"
  set stdvZfile "tau${wind_size}f/stdev-z_out.txt"
  set outMSDx [open $msdXfile w]
  set outMSDy [open $msdYfile w]
  set outMSDz [open $msdZfile w]
  set outStdvX [open $stdvXfile w]
  set outStdvY [open $stdvYfile w]
  set outStdvZ [open $stdvZfile w]
  for {set i 0} {$i < $wind_size} {incr i} {
	  puts $outMSDx [format "%.4f" [lindex $MSDx $i]]
	  puts $outMSDy [format "%.4f" [lindex $MSDy $i]]
	  puts $outMSDz [format "%.4f" [lindex $MSDz $i]]
	  puts $outStdvX [format "%.4f" [lindex $stdvX $i]]
	  puts $outStdvY [format "%.4f" [lindex $stdvY $i]]
	  puts $outStdvZ [format "%.4f" [lindex $stdvZ $i]]
  }
  close $outMSDx
  close $outMSDy
  close $outMSDz
  close $outStdvX
  close $outStdvY
  close $outStdvZ
}
