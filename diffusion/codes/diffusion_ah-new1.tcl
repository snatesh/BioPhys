proc diffusion {axis abeg aend nT com_on} {

  set sel [atomselect top "name OH2"]
  set sel0 [atomselect top "name OH2"]
  set nummol [$sel num]
  	
  #--------------Compute MSD for each time window, append results and average--------#
  set MSD [veccreate $nT]
  set nw [expr {($aend-$abeg)/$nT}] 
  set norm [expr {$nummol*$nw}]
  set ax_name [axname $axis]
  if {$com_on == 1} {
      puts "com_on"
      for {set j $abeg} {$j<$aend} {incr j $nT} {
        $sel0 frame $j
        $sel0 update
        $sel frame $j
        $sel update
        set coms0 [measure center $sel0 weight mass]
        # set coms0_xy [list [lindex $coms0  0] [lindex $coms0 1]]	
        set com0 [lindex $coms0 $axis]
        set coords0 [$sel0 get $ax_name]
	puts "init center of mass, coords set"
        # set coords0 [$sel0 get {x y}]
        set sdlist {}
        for {set i 0} {$i<$nT} {incr i} {
  		$sel frame [expr {$j + $i}]
  		# $sel update
  		set coms [measure center $sel weight mass]
  		# set coms_xy [list [lindex $coms 0] [lindex $coms 1]]
  		set com [lindex $coms $axis]
  		set coords [$sel get $ax_name]
		puts "new center of mass, coords set"
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
		puts "sd appended"
	}
	puts $sdlist
	puts $MSD
	set MSD [vecadd $MSD $sdlist]
	puts "sd added to MSD"
      }
  } else {
      for {set j $abeg} {$j<=$aend} {incr j $nT} {
	$sel0 frame $j
	$sel0 update
	$sel frame $j
	$sel update
	set coords0 [$sel0 get $ax_name]
	set sdlist {}
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
      }
  }
  set MSD [vecscale $MSD [expr {1.0/$norm}]]
  puts "MSD scaled"
  return $MSD
}

#--------------Procedure to create numeric vector of length dim--------#
  proc veccreate {dim} {
    set v {}
    for {set i 0} {$i < $dim} {incr i} {
  	lappend v 0
    }
    return $v
  }
  

#-------------Procedure to convert ax_name input to ax_name string------------#
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
