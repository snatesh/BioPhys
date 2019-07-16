#-------------Procedure to create list of lists of length dim---------------------------------#
proc eveccreate {dim} {
	set v {}
	for {set i 0} {$i < $dim} {incr i} {
        	lappend v {}
        }
        return $v
}

#-------------Procedure to create list of lists of shells of increasing distance--------------#
proc shells {frame max_radius thick} {
	set prot [atomselect top "protein and noh and resid 32 33 35 40" frame $frame]
	set prot_ids [$prot get index]

	set cen [measure center $prot]
	set cx [lindex $cen 0]
	set cy [lindex $cen 1]

	set water [atomselect top "name OH2 and sqr(x-$cx) + sqr(y-$cy) <= 10" frame $frame]
        set ids [$water get index]
        set nummol [llength $ids]
        set dists {}

	# Generate distances from the center
	foreach {id} $ids {
		#set sel0 [atomselect top "index [lindex $ids $i]" frame $frame]
		#set coords [$sel0 get {x y}]
		set bond {}
		lappend bond $id
		set dist0 $max_radius
		foreach {prid} $prot_ids {
			lappend bond $prid
			set dist [measure bond $bond]
			if {$dist < $dist0} {set dist0 $dist}
			set bond [lreplace $bond end end]
		}
		#set x [lindex [lindex $coords 0] 0]
		#set y [lindex [lindex $coords 0] 1]
		#set dx [expr {$x - $cx}]
		#set dy [expr {$y - $cy}]
		#set dist0 [expr {sqrt(pow($dx, 2) + pow($dy, 2))}]
		lappend dists $dist0
	}

	# Create shells
	set numshells [expr {int(ceil($max_radius / $thick) + 1)}]
	set shellnums {}
	set shells [eveccreate $numshells]
	foreach {dist} $dists {
		set num [expr {int(floor($dist / $thick))}]
		lappend shellnums $num
	}
	for {set i 0} {$i < $nummol} {incr i} {
		set idx [lindex $shellnums $i]
		set shell_i [lindex $shells $idx]
		lappend shell_i [lindex $ids $i]
		lset shells $idx $shell_i
	}
	return $shells
}

