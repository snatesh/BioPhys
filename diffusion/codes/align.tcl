proc align {selections mol} {
	set sel [atomselect $mol $selections]
	set sel0 [atomselect $mol $selections frame 0]
	set sela [atomselect $mol all]
	set n [molinfo $mol get numframes]
	for {set i 0} {$i < $n} {incr i} {
		$sel frame $i
		$sel update
		set trans_mat [measure fit $sel $sel0]
		$sela frame $i
		$sela update
		$sela move $trans_mat
	}
}

