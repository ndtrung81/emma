
proc visualize {molid res radius} {

	graphics $molid delete all
	set select_atoms [atomselect $molid "all"]
 
	set natoms [$select_atoms num]

	set natomspercopy 8
  set rad $radius
# set res 16

	set ncopies [expr $natoms/$natomspercopy]

  puts "Drawing $ncopies copies"

	if {$natoms <= 0} {
		error "No atoms selected"
	}

	for {set i 0} {$i < $ncopies} {incr i}  {

		set startidx [expr $i * $natomspercopy]
    draw_square $molid $startidx blue $res $radius

		set startidx [expr $i * $natomspercopy + 4]
    draw_square $molid $startidx orange $res $radius

	}
}

proc draw_square {molid startidx color res rad} {
	graphics $molid color $color	
  # edges
  set atom_index1 [expr $startidx + 0]
	set atom_index2 [expr $startidx + 1]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

	set atom_index1 [expr $startidx + 1]
	set atom_index2 [expr $startidx + 2]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

	set atom_index1 [expr $startidx + 2]
	set atom_index2 [expr $startidx + 3]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

  set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res

	set atom_index1 [expr $startidx + 3]
	set atom_index2 [expr $startidx + 0]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

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
  lappend newcoord1 $rad

	set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord2 $x
  lappend newcoord2 $y
  lappend newcoord2 $rad

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord3 $x
  lappend newcoord3 $y
  lappend newcoord3 $rad
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
  lappend newcoord11 $rad

  set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord21 $x
  lappend newcoord21 $y
  lappend newcoord21 $rad

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord31 $x
  lappend newcoord31 $y
  lappend newcoord31 $rad
	graphics $molid triangle $newcoord11 $newcoord21 $newcoord31
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
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

	set atom_index1 [expr $startidx + 1]
	set atom_index2 [expr $startidx + 2]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

	set atom_index1 [expr $startidx + 2]
	set atom_index2 [expr $startidx + 0]
	set select_tip1 [atomselect $molid "index $atom_index1"]
	set select_tip2 [atomselect $molid "index $atom_index2"]

	set coord1 [lindex [$select_tip1 get {x y z}] 0]
	set coord2 [lindex [$select_tip2 get {x y z}] 0]
	graphics $molid cylinder $coord1 $coord2 radius $rad resolution $res filled yes
	graphics $molid sphere $coord1 radius $rad resolution $res
	graphics $molid sphere $coord2 radius $rad resolution $res	

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
  lappend newcoord1 $rad

	set coord2 [lindex [$select_tip2 get {x y z}] 0]
  set x [lindex $coord2 0]
  set y [lindex $coord2 1]
  lappend newcoord2 $x
  lappend newcoord2 $y
  lappend newcoord2 $rad

	set coord3 [lindex [$select_tip3 get {x y z}] 0]
  set x [lindex $coord3 0]
  set y [lindex $coord3 1]
  lappend newcoord3 $x
  lappend newcoord3 $y
  lappend newcoord3 $rad
	graphics $molid triangle $newcoord1 $newcoord2 $newcoord3
}

set diameter [lindex $argv 0]
set radius [expr $diameter * 0.5 ]

mol modstyle 0 0 Points
visualize 0 32 $radius
