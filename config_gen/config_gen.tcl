# Sachin Natesh
# <sachin.natesh@gmail.com>


# procedure that takes a pdb as input and 
#     generates a psf for the pdb
#     solvates and neutralizes the system, generating corresponding psf and pdb
#       - with 10A padding and rotating protein to minimize volume
#     restrains and centers the system, generating corresponding pdbs
#     note psf uses top_all36_prot.rtf and toppar_water_ions_esmael.str (should be in working dir)
#       - if needed, change topology info in 'autopsf' command

# usage: config_gen input.pdb output_name

package require psfgen
package require solvate
package require autoionize

puts "usage: config_gen input.pdb output_name"
puts "output_name should be same as in namd .conf file"
puts "WARNING: ALL MOLECULES ARE DELETED AT START OF PROC"

proc config_gen {args} {
   
  mol delete all
 
  # error message for incorrect number of command line inputs
  if {![llength $args]} {
    puts "incorrect number of arguments"
    puts "usage: config_gen input.pdb output_name"
    error ""
  }

   
  # setting prefix name
  set prefix "[string trimright [lindex $args 0] ".pdb"]"

  mol new $prefix.pdb
  set id [molinfo top]

  ## begin writing auto_psf/pdb, solvating with 10A and neutralizing ##

  # writing psf for input pdb
  resetpsf
  					# LIST TOPOLOGIES BY RELATIVE PATH HERE 
  autopsf -mol $id -prefix $prefix -top {top_all36_prot.rtf toppar_water_ions_esmael.str}  
  mol delete $id
  set id [molinfo top]
  # solvating input pdb with 10A padding and rotating to minimize volume
  solvate "$prefix.psf" "$prefix.pdb" -o "$prefix-solvate" -rotate -t 10
  mol delete $id
  set id [molinfo top]
  autoionize -psf "$prefix-solvate.psf" -pdb "$prefix-solvate.pdb" -neutralize -o "$prefix-solvate-neutralized"
  mol delete $id
  
  ## begin restraining and centering ##
  
  mol new "$prefix-solvate-neutralized.psf"
  mol addfile "$prefix-solvate-neutralized.pdb"
   
  # centering
  set sel [atomselect top all]
  set v [measure center $sel]
  set nv [vecscale -1 $v] 
  $sel moveby $nv 
  $sel writepdb $prefix-solvate-neutralized-centered.pdb
  mol delete top
  
  # restraining
  mol load pdb $prefix-solvate-neutralized-centered.pdb
  set sel [atomselect top "protein and backbone"]
  $sel set beta 10
  set sel [atomselect top all]
  $sel writepdb $prefix-solvate-neutralized-restrain.pdb
  
  # estimating cell size
  set v [measure minmax $sel]
  set v1 [lindex $v 0]
  set v2 [lindex $v 1]
  set p1 [expr { abs([lindex $v1 0]) + abs([lindex $v2 0]) }]
  set p2 [expr { abs([lindex $v1 1]) + abs([lindex $v2 1]) }]
  set p3 [expr { abs([lindex $v1 2]) + abs([lindex $v2 2]) }]
  
  set x [format "%.2f" $p1]
  set y [format "%.2f" $p2]
  set z [format "%.2f" $p3]

  puts "x basis: $p1"
  puts "y basis: $p2"
  puts "z basis: $p3"
  mol delete [molinfo top] 
  puts "structure file generation complete"
  
  ## begin generating config file ##
  puts "generating config file ..."
   
  set outputname [lindex $args 1]
 
  exec cp template.conf tmp.conf 
  exec sed "s/coordinates/coordinates	$prefix-solvate-centered.pdb/;  
   	    s/structure/structure	$prefix-solvate-neutralized.psf/;  
   	    s/outputname/outputname	$outputname/;  
   	    s/restartname/restartname	$outputname/;   
   	    s/Consref/Consref		$prefix-solvate-neutralized-restrain.pdb;/;  	
   	    s/ConsKFile/ConsKFile	$prefix-solvate-neutralized-restrain.pdb;/;  
   	    s/cellBasisVector1/cellBasisVector1 	$x	0.00	0.00/;   
   	    s/cellBasisVector2/cellBasisVector2	0.00	$y	0.00/;  
   	    s/cellBasisVector3/cellBasisVector3	0.00	0.00	$z/" tmp.conf > $outputname.conf

   exec rm tmp.conf
   puts "all files generated. see working directory" 
  
}




 
   
