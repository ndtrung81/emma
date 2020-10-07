
proc visualize {molid res} {

	graphics $molid delete all
	set select_atoms [atomselect $molid "all"]
 
	set natoms [$select_atoms num]

	set natomspercopy 5

# set res 16

	set ncopies [expr $natoms/$natomspercopy]

  puts "Drawing $ncopies copies"

	if {$natoms <= 0} {
		error "No atoms selected"
	}

	for {set i 0} {$i < $ncopies} {incr i}  {

    if {$i % 2 == 0} {
      set color lime
    } else {
      set color mauve
    }

		set startidx [expr $i * $natomspercopy]
    draw_triangle $molid $startidx $color $res

		set startidx [expr $i * $natomspercopy + 2]
    draw_triangle $molid $startidx $color $res

		set startidx [expr $i * $natomspercopy]
		set idx1 0
		set idx2 2
		set idx3 4
    set drawing_edges 1
    draw_triangle3 $molid $startidx $idx1 $idx2 $idx3 $color $res $drawing_edges

	}
}

proc draw_square {molid startidx color res} {
	graphics $molid color $color	
  # edges
  set atom_index1 [expr $startidx + 0]
	set atom_index2 [expr $startidx + 1]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + 1]
	set atom_index2 [expr $startidx + 2]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + 2]
	set atom_index2 [expr $startidx + 3]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

  set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res

	set atom_index1 [expr $startidx + 3]
	set atom_index2 [expr $startidx + 0]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

  # faces
  set atom_index1 [expr $startidx + 0]
	set atom_index2 [expr $startidx + 1]
  set atom_index3 [expr $startidx + 2]

	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]
	set select_tip3 [atomselect $molid "index $atom_index3"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
  set x [lindex $coord1 0]
  set y [lindex $coord1 1]
  lappend newcoord1 $x
  lappend newcoord1 $y
  lappend newcoord1 0.5

	set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord2 $x
  lappend newcoord2 $y
  lappend newcoord2 0.5

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord3 $x
  lappend newcoord3 $y
  lappend newcoord3 0.5
	graphics $molid triangle $newcoord1 $newcoord2 $newcoord3
	
  set atom_index1 [expr $startidx + 2]
	set atom_index2 [expr $startidx + 3]
  set atom_index3 [expr $startidx + 0]

	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]
	set select_tip3 [atomselect $molid "index $atom_index3"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
  set x [lindex $coord1 0]
  set y [lindex $coord1 1]
  lappend newcoord11 $x
  lappend newcoord11 $y
  lappend newcoord11 0.5

  set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord21 $x
  lappend newcoord21 $y
  lappend newcoord21 0.5

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord31 $x
  lappend newcoord31 $y
  lappend newcoord31 0.5
	graphics $molid triangle $newcoord11 $newcoord21 $newcoord31
}

proc draw_triangle3 {molid startidx idx1 idx2 idx3 color res drawning_edges} {
	graphics $molid color $color
  if {$drawning_edges == 1} {

  # edges
  set atom_index1 [expr $startidx + $idx1]
	set atom_index2 [expr $startidx + $idx2]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + $idx2]
	set atom_index2 [expr $startidx + $idx3]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + $idx3]
	set atom_index2 [expr $startidx + $idx1]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	
  }

  # faces
  set atom_index1 [expr $startidx + $idx1]
	set atom_index2 [expr $startidx + $idx2]
  set atom_index3 [expr $startidx + $idx3]

	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]
	set select_tip3 [atomselect $molid "index $atom_index3"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
  set x [lindex $coord1 0]
  set y [lindex $coord1 1]
  lappend newcoord1 $x
  lappend newcoord1 $y
  lappend newcoord1 0.5

	set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord2 $x
  lappend newcoord2 $y
  lappend newcoord2 0.5

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord3 $x
  lappend newcoord3 $y
  lappend newcoord3 0.5
	graphics $molid triangle $newcoord1 $newcoord2 $newcoord3
}


proc draw_triangle {molid startidx color res} {
	graphics $molid color $color
  # edges
  set atom_index1 [expr $startidx + 0]
	set atom_index2 [expr $startidx + 1]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + 1]
	set atom_index2 [expr $startidx + 2]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

	set atom_index1 [expr $startidx + 2]
	set atom_index2 [expr $startidx + 0]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius 0.5 resolution $res filled yes
	graphics $molid sphere $coord1 radius 0.5 resolution $res
	graphics $molid sphere $coord2 radius 0.5 resolution $res	

  # faces
  set atom_index1 [expr $startidx + 0]
	set atom_index2 [expr $startidx + 1]
  set atom_index3 [expr $startidx + 2]

	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]
	set select_tip3 [atomselect $molid "index $atom_index3"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
  set x [lindex $coord1 0]
  set y [lindex $coord1 1]
  lappend newcoord1 $x
  lappend newcoord1 $y
  lappend newcoord1 0.5

	set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord2 $x
  lappend newcoord2 $y
  lappend newcoord2 0.5

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord3 $x
  lappend newcoord3 $y
  lappend newcoord3 0.5
	graphics $molid triangle $newcoord1 $newcoord2 $newcoord3
}

visualize 0 16
