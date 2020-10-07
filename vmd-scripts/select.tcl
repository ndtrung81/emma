#   export LAMMPSREMAPFIELDS=q=c_q[1],vx=c_q[2],vy=c_q[3],vz=c_q[4]
# or
#   export LAMMPSREMAPFIELDS=vx=mux,vy=muy,vz=muz
#global env
#set env(LAMMPSREMAPFIELDS) vx=mux,vy=muy,vz=muz

proc write2file {select_text outputfile} {
  set select_atoms [atomselect top $select_text]
  set natoms  [$select_atoms num]
  
  if {$natoms <= 0} {
    error "No atoms selected"
  }

  set output [open $outputfile w]
  puts $output 

  for {set i 0} {$i < $natoms} {incr i} {  
    # get the atom coordinates
	  set coords [lindex [$select_atoms get {x y z}] $i]
    
    set x [lindex $coords 0]
    set y [lindex $coords 1]
    set z [lindex $coords 2]

    puts $output $coords

  }

  close $output
}


write2file "type AC" "AC.txt"

draw material AOEdgy

color Display Background white
#mol modstyle 0 0 points

axes location off
display cuedensity 0.20000
display cuemode Exp2 
display shadows on
display ambientocclusion on

#mol delrep 0

#mol selection {type 1 2 3 4}
#mol selection {type 1}
#mol addrep top
#mol modstyle   1 0 VDW 0.333 16.000

#mol selection {type 6 7}
#mol addrep top
#mol modstyle   2 0 VDW 1.0 16.000



