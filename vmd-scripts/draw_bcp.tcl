proc visualize {molid resolution} {

	draw delete all
  lassign [molinfo top get "{selection 0}"] select_text
  set select_atoms [atomselect top $select_text]
  set natoms  [$select_atoms num]

	if {$natoms <= 0} {
		error "No atoms selected"
	}

  set radius 0.3

	set natomsperNP 10
  set nA 3
  set nB [expr $natomsperNP - $nA]
  set nAm1 [expr $nA - 1]
  set nchains [expr $natoms/$natomsperNP]

  puts "Num atoms selected $natoms; num chains = $nchains"
  puts "nAm1 = $nAm1; nB = $nB"

	for {set i 0} {$i < $nchains} {incr i}  {

    puts "draw chain number $i of $nchains"

    # Block A

		draw color white

    for {set j 0} {$j < $nAm1} {incr j} {
  		set atom_index1 [expr $i * $natomsperNP + $j + 0]
	  	set atom_index2 [expr $i * $natomsperNP + $j + 1]
    
      set coord1 [lindex [$select_atoms get {x y z}] $atom_index1]
      set coord2 [lindex [$select_atoms get {x y z}] $atom_index2]

	  	draw cylinder $coord1 $coord2 radius $radius resolution $resolution filled yes
		  draw sphere $coord1 radius $radius resolution $resolution
  		draw sphere $coord2 radius $radius resolution $resolution
    }

    # Block B

    draw color blue

    for {set j 0} {$j < $nB} {incr j} {
  		set atom_index1 [expr $i * $natomsperNP + $j + 2]
     	set atom_index2 [expr $i * $natomsperNP + $j + 3]
    
      set coord1 [lindex [$select_atoms get {x y z}] $atom_index1]
      set coord2 [lindex [$select_atoms get {x y z}] $atom_index2]

	  	draw cylinder $coord1 $coord2 radius $radius resolution $resolution filled yes
	    draw sphere $coord1 radius $radius resolution $resolution
    	draw sphere $coord2 radius $radius resolution $resolution

#      puts "draw bonds between $atom_index1 and $atom_index2"
    }
	}
}

color Display Background white
pbc box -center bb -style tubes -width 1.0 -color silver
mol modstyle 0 0 Points
visualize 0 16

