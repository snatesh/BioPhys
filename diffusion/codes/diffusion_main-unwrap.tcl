source /home/sachin/Research/BioPhys/diffusion/codes/diffusion_new.tcl
source /home/sachin/Research/BioPhys/diffusion/codes/align.tcl
package require pbctools

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

foreach {j k} $wind_and_ref_list {
	puts "$j $k"
	set wind_size $j
	set next_ref $k	
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
		file mkdir "tau${j}f"
		# get diffusion in x, y and z
		diffusion "tau${j}f/2lmp-unwrap-x-$i-$last-diffData" 0 0 $wind_size $wind_size .02 1
		diffusion "tau${j}f/2lmp-unwrap-y-$i-$last-diffData" 1 0 $wind_size $wind_size .02 1
		diffusion "tau${j}f/2lmp-unwrap-z-$i-$last-diffData" 2 0 $wind_size $wind_size .02 1
		# cleanup
		$sel delete
		$sel1 delete
		mol delete top
		file delete $i-$last.dcd
		file delete $i-$last.psf
	}
}
