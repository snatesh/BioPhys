# Procedure to compute the self diffusivity of an atomic species given 
# its trajectory through configuration space
#
# Sachin Natesh
# <sachin.natesh@gmail.com>
#
# Computed as an ensemble average; i.e. over all atoms and representative number
# of time origins
#
# Example usage:
# set sel [atomselect top all]
# source diffusion5.tcl
# diffusion "diff_out.txt" "resname SOD" 2 100 600 100 .001 
# selection must be single atomic species
# axis is that along which diffusion calculated (values= 0, 1, 2 : x, y, z)
# center of mass drift accounted for by default
# Make sure to take out center of mass drift if computing diffusion
# for single atom; i.e. if system diffusion is reduced to com diffusion 
# saved in a dummy atom
# abeg/end=first/last frame to be analysed  (frames start at 0)

# nT is number of frames in each time window (use appropriate multiple of abeg-aend)

proc diffusion {outfile axis abeg aend nT ts com_on} {

set out [open $outfile w]
set sel [atomselect top "name OH2"]
#"sqr(x + 3.23) + sqr(y - 9.617) <= sqr(10.3) and resname SOD"]
set sel0 [atomselect top "name OH2"]
# "sqr(x + 3.23) + sqr(y - 9.617) <= sqr(10.3) and resname SOD"]
set nummol [$sel num]

#----------------Procedure to convert axis input to axis string------------#
proc axname {axis} {
if {$axis==0} {
        set ax_name "x"
        return $ax_name
} elseif {$axis==1} {
        set ax_name "y"
        return $ax_name
} else {
        set ax_name "z"
        return $ax_name

}
}
set ax_name [axname $axis]

#----------------Procedure to create list of numbers between two inputs------#
proc num_list {first last} {
set num_list {}
set i $first
while {$i<=$last} {
	lappend num_list $i
	set i [expr {$i + 1}]
	}
return $num_list
}

#--------------Procedure to create numeric vector of length dim--------#
proc veccreate {dim} {
set v {}
for {set i 0} {$i < $dim} {incr i} {
	lappend v 0
	}
return $v
}
#-------------Generate list of lists of frames in each time window-------#

set nw [expr {($aend-$abeg)/$nT}] 
set wlist {}
for {set i 0} {$i<=$nw} {incr i} {
	lappend wlist [num_list [expr {$i*$nT}] [expr {($i + 1)*$nT}]]
	}
	
#--------------Compute MSD for each time window, append results and average--------#
set MSD [veccreate $nT]
set norm [expr {$nummol*$nw}]
if {$com_on == 1} {
for {set j $abeg} {$j<$aend} {incr j $nT} {
	$sel0 frame $j
	$sel0 update
	$sel frame $j
	$sel update
	set coms0 [measure center $sel0 weight mass]
	# set coms0_xy [list [lindex $coms0  0] [lindex $coms0 1]]	
	set com0 [lindex $coms0 $axis]
	set coords0 [$sel0 get $ax_name]
	# set coords0 [$sel0 get {x y}]
	set sdlist {}
	puts "frame $j \t[$sel0 num]"
	for {set i 0} {$i<$nT} {incr i} {
		$sel frame [expr {$j + $i}]
		# $sel update
		set coms [measure center $sel weight mass]
		# set coms_xy [list [lindex $coms 0] [lindex $coms 1]]
		set com [lindex $coms $axis]
		set coords [$sel get $ax_name]
		# set coords [$sel get {x y}]
		set sd 0
		for {set k 0} {$k < $nummol} {incr k} {
			set rt0t [expr {[lindex $coords $k] - $com}]
			# set rt0t [expr {[veclength [vecsub [lindex $coords $k] $coms_xy]]}]
			set rt0 [expr {[lindex $coords0 $k] - $com0}]
			# set rt0 [expr {[veclength [vecsub [lindex $coords0 $k] $coms0_xy]]}]
			set d [expr {$rt0t - $rt0}]
			set sd [expr {$sd + $d*$d}]
			}
		lappend sdlist $sd
		}
	set MSD [vecadd $MSD $sdlist]
	puts "computing: [expr {floor([expr {($j-$abeg)*100/($aend - $abeg)}])}]%"
	}
} else {
for {set j $abeg} {$j<=$aend} {incr j $nT} {
        $sel0 frame $j
        $sel0 update
        $sel frame $j
        $sel update
	set coords0 [$sel0 get $ax_name]
	set sdlist {}
  	puts "frame $j \t[$sel0 num]"
        for {set i 0} {$i<$nT} {incr i} {
                $sel frame [expr {$j + $i}]
		set coords [$sel get $ax_name]
		set sd 0
                for {set k 0} {$k < $nummol} {incr k} {
			set rt0t [expr {[lindex $coords $k]}]
			set rt0 [expr {[lindex $coords0 $k]}]
			set d [expr {$rt0t - $rt0}]
                	set sd [expr {$sd + $d*$d}]
			}
                lappend sdlist $sd
                }
        set MSD [vecadd $MSD $sdlist]
        puts "computing: [expr {floor([expr {($j-$abeg)*100/($aend - $abeg)}])}]%"
        }
}

set MSD [vecscale $MSD [expr {1.0/$norm}]]
for {set i 0} {$i<$nT} {incr i} {
	# puts $out "[expr {$ts*$i}] \t [lindex $MSD $i]"
	puts $out [format "%.4f" [lindex $MSD $i]]
        }
close $out
}


	
		
	
	

