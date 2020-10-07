proc arrows {scale_v scale_geo resolution color} {
	graphics top delete all
  lassign [molinfo top get "{selection 0}"] select_text

  set select_atoms [atomselect top $select_text]
  set natoms  [$select_atoms num]
  
  if {$natoms <= 0} {
    error "No atoms selected"
  }
  puts "Drawing $natoms arrows"
  
  draw color $color
  set scale1 [expr $scale_geo * 0.01 ]
  set scale2 [expr $scale_geo * 0.04 ]

	for {set i 0} {$i < $natoms} {incr i} {  
    # get the atom coordinates
	  set coords1 [lindex [$select_atoms get {x y z}] $i]

    # get the velocity components from file
    set coords3 [lindex [$select_atoms get {vx vy vz}] $i]
    
    # scale the velocity vector
    set coords3 [vecscale $scale_v $coords3]

    # base point of the arrow cone
    set coords2 [vecscale 0.7 $coords3]
    set coords2 [vecadd $coords1 $coords2]

    # end point of the arrow (tip)
    set coords3 [vecadd $coords1 $coords3]
			
    # draw arrow body and tip
    draw cylinder $coords1 $coords2 radius $scale1
    draw cone     $coords2 $coords3 radius $scale2
	}
}
