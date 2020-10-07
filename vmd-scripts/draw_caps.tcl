## 
## Modified from ellipse.tcl by chris forman <cjf41@cam.ac.uk>
## by Trung Nguyen (ORNL)
##

# Replace the header for 4 quaterion components with charge (q) and velocities
# sed "s/c_q\[1\]/q/g;s/c_q\[2\]/vx/g;s/c_q\[3\]/vy/g;s/c_q\[4\]/vz/g" < dump.txt > dump_modified.txt
# vmd -lammpstrj dump_modified.txt
# or set the environment variable LAMMPSREMAPFIELDS before opening the original dump file
# export LAMMPSREMAPFIELDS=q=c_q[1],vx=c_q[2],vy=c_q[3],vz=c_q[4]
# vmd -lammps dump.txt

# How to use:
# 1- Modify the first selection for the atoms of interest: 
# e.g. type 2 and z > 0 and z < 10 and x > 30 and x < 50 and y > 30 and y < 50
# 2- Open the Tcl Console
# 3- Load this script: source drawlc.tcl
# 4- Invoke the draw function: visualframe 3 1 1 2.0 12 green2
# where 1 1 3 are the lengths of the ellipsoids in x, y and z
#       2.0    = the scale factor for the arrow length
#       12     = the resolution
#       green2 = the arrow color

proc f {a u v} {
  expr {$a*sin($u) * cos($v)} 
}  

proc g {b u v} {
  expr {$b*cos($u) * cos($v)}
}

proc h {c u v} {
  expr {$c*sin($v)} 
}

proc mesh {a b c theta res} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv

  # to cover the whole surface: u runs from 0 to 2*pi and v runs from -pi to pi
  set PI 3.141592654
  set theta_rad [expr $theta * $PI / 180.0]
  set minu 0
  set maxu [expr $PI * 2]
#  set minu [expr $PI / 3]
#  set maxu [expr $PI / 2]
  set resu $res
  set stepu [expr ($maxu - $minu) / $resu]

  # v runs from -pi to pi
#  set minv [expr -1*$PI]
#  set maxv [expr $PI]
  set minv [expr $PI/2-$theta_rad]
  set maxv [expr $PI/2]
  set resv $res
  set stepv [expr ($maxv - $minv) / $resv]

  for {set iu 0} {$iu <= $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv <= $resv} {set iv [expr $iv + 1]} {
      set u [expr $minu + $iu * $stepu]
      set v [expr $minv + $iv * $stepv]
      set fdata($u,$v) [f $a $u $v]
      set gdata($u,$v) [g $b $u $v]
      set hdata($u,$v) [h $c $u $v]     
    }
  }
}

# a = semi-axis in x direction
# b = semi-axis in y direction
# c = semi-axis in z direction
# iquatw = real component of the quaternion
# iquati/j/k = vector components of the quaternion
# x = translation in x direction.
# y = translation in y direction
# z = translation in z direction
# orient = relative orientation of neighbors (read in)
# eColour is the colour of the ellipse in VMD colour table.
# res = resolution
# mode = 0: read-in relative orientaion
#      = 1: atom orientaion wrt z axis

proc draw_ellipsoid_quat { a b c iquatw iquati iquatj iquatk x y z color res} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv
  global minColor maxColor NrColors coloringMode

  # first normalize the quaternion
  set quatmag [expr sqrt($iquatw*$iquatw + $iquati*$iquati + $iquatj*$iquatj + $iquatk*$iquatk)]
  set quatw   [expr $iquatw/$quatmag]
  set quati   [expr $iquati/$quatmag]
  set quatj   [expr $iquatj/$quatmag]
  set quatk   [expr $iquatk/$quatmag]

  # compute rotation matrix to get desired orientation
  set a11 [expr $quatw*$quatw + $quati*$quati - $quatj*$quatj - $quatk*$quatk]
  set a21 [expr 2.0 * ($quati*$quatj + $quatw*$quatk)]
  set a31 [expr 2.0 * ($quati*$quatk - $quatw*$quatj)]
  set a12 [expr 2.0 * ($quati*$quatj - $quatw*$quatk)]
  set a22 [expr $quatw*$quatw - $quati*$quati + $quatj*$quatj - $quatk*$quatk]
  set a32 [expr 2.0 * ($quatj*$quatk + $quatw*$quati)]
  set a13 [expr 2.0 * ($quati*$quatk + $quatw*$quatj)]
  set a23 [expr 2.0 * ($quatj*$quatk - $quatw*$quati)]
  set a33 [expr $quatw*$quatw - $quati*$quati - $quatj*$quatj + $quatk*$quatk]

  set cnum [colorinfo num]
  set cmax [colorinfo max]

  # uniaxial order parameter S = (3 (ux dot ez)^2 - 1) / 2
  # S = -0.5 for ux perpendicular to ez
  # S = 1.0 for ux parallel or anti-parallel to ez

#  if {$coloringMode == 0} {
#    set normS $orient
#  } else {
#   set S [expr (3.0 * $a31 * $a31 - 1.0)/2.0]
#   set normS [expr ($S + 0.5)/(1.0 + 0.5)]
#  }

#  draw color [expr $minColor + int($normS * $NrColors)]

  draw color $color

  #perform rotation and translation
  for {set iu 0} {$iu <= $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv <= $resv} {set iv [expr $iv + 1]} {
      set u [expr $minu + $iu * $stepu]
      set v [expr $minv + $iv * $stepv]
      # note the rotation matrix is tranposed here
      set fdata_r($iu,$iv) [expr {($fdata($u,$v)*$a11+$gdata($u,$v)*$a12+$hdata($u,$v)*$a13) + $x}]
      set gdata_r($iu,$iv) [expr {($fdata($u,$v)*$a21+$gdata($u,$v)*$a22+$hdata($u,$v)*$a23) + $y}]
      set hdata_r($iu,$iv) [expr {($fdata($u,$v)*$a31+$gdata($u,$v)*$a32+$hdata($u,$v)*$a33) + $z}]
    }
  }

  for {set iu 0} {$iu < $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv < $resv} {set iv [expr $iv + 1]} {

      # get the next two corners 
      set iu2 [expr $iu + 1]
      set iv2 [expr $iv + 1]
      draw triangle "$fdata_r($iu,$iv)  $gdata_r($iu,$iv)  $hdata_r($iu,$iv)" \
                    "$fdata_r($iu2,$iv)  $gdata_r($iu2,$iv)  $hdata_r($iu2,$iv)" \
                    "$fdata_r($iu2,$iv2) $gdata_r($iu2,$iv2) $hdata_r($iu2,$iv2)"
      draw triangle "$fdata_r($iu2,$iv2) $gdata_r($iu2,$iv2) $hdata_r($iu2,$iv2)" \
                    "$fdata_r($iu,$iv2) $gdata_r($iu,$iv2) $hdata_r($iu,$iv2)" \
                    "$fdata_r($iu,$iv)  $gdata_r($iu,$iv)  $hdata_r($iu,$iv)"
    }
  }
}

proc draw_arrow_quat { x y z iquatw iquati iquatj iquatk color res scale} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv
  global minColor maxColor NrColors coloringMode

  # first normalize the quaternion
  set quatmag [expr sqrt($iquatw*$iquatw + $iquati*$iquati + $iquatj*$iquatj + $iquatk*$iquatk)]
  set quatw   [expr $iquatw/$quatmag]
  set quati   [expr $iquati/$quatmag]
  set quatj   [expr $iquatj/$quatmag]
  set quatk   [expr $iquatk/$quatmag]

  # compute rotation matrix to get desired orientation
  set a11 [expr $quatw*$quatw + $quati*$quati - $quatj*$quatj - $quatk*$quatk]
  set a21 [expr 2.0 * ($quati*$quatj + $quatw*$quatk)]
  set a31 [expr 2.0 * ($quati*$quatk - $quatw*$quatj)]
  set a12 [expr 2.0 * ($quati*$quatj - $quatw*$quatk)]
  set a22 [expr $quatw*$quatw - $quati*$quati + $quatj*$quatj - $quatk*$quatk]
  set a32 [expr 2.0 * ($quatj*$quatk + $quatw*$quati)]
  set a13 [expr 2.0 * ($quati*$quatk + $quatw*$quatj)]
  set a23 [expr 2.0 * ($quatj*$quatk - $quatw*$quati)]
  set a33 [expr $quatw*$quatw - $quati*$quati - $quatj*$quatj + $quatk*$quatk]

  draw color $color 

  # get the atom coordinates
	set coords1 ""
  lappend coords1 $x
  lappend coords1 $y
  lappend coords1 $z

  # get the velocity components from file
  set coords3 ""
  lappend coords3 $a13
  lappend coords3 $a23
  lappend coords3 $a33
    
  # scale the ez vector/ arrow length
  set arrow_len [expr $scale * 1.5]
  set coords3 [vecscale $arrow_len $coords3]

  # base point of the arrow cone
  set coords2 [vecscale 0.7 $coords3]
  set coords2 [vecadd $coords1 $coords2]

  # end point of the arrow (tip)
  set coords3 [vecadd $coords1 $coords3]
			
  # draw arrow body and tip
  set radius_cyl [expr $scale * 0.15]
  draw cylinder $coords1 $coords2 radius $radius_cyl
  set radius_cone [expr $scale * 0.3]
  draw cone     $coords2 $coords3 radius $radius_cone
}


# shapex = dimension in x direction
# shapey = dimension in y direction
# shapez = dimension in z direction
# resolution
# mode

proc visualframe {shapex shapey shapez theta scale resolution color} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv
  global minColor maxColor NrColors colorScheme coloringMode

  set colorScheme BGryR
  color scale method $colorScheme
  color scale min 0
  color scale max 1
  color scale midpoint 0.5
  
  set minColor [expr [colorinfo num]]
  set maxColor [expr [colorinfo max] - 1]
  set NrColors [expr $maxColor- $minColor- 1]

  graphics top delete all
  lassign [molinfo top get "{selection 0}"] select_text

  set select_atoms [atomselect top $select_text]
  set natoms  [$select_atoms num]

  if {$natoms <= 0} {
    error "No atoms selected"
  }

  # copy radius into user, set radius to 0.5
  $select_atoms set user [$select_atoms get radius]
# $select_atoms set radius 0.5

  puts "Drawing $natoms arrows"
    
  set a [expr $shapex / 2.0]
  set b [expr $shapey / 2.0]
  set c [expr $shapez / 2.0]
  
  mesh $a $b $c $theta $resolution

	for {set i 0} {$i < $natoms} {incr i} {  
    # get the atom coordinates
	  set coords [lindex [$select_atoms get {x y z}] $i]

    set x [lindex $coords 0]
    set y [lindex $coords 1]
    set z [lindex $coords 2]
    
    # get the quaternion components from file
    set qijk [lindex [$select_atoms get {vx vy vz}] $i]
    set qi [lindex $qijk 0]
    set qj [lindex $qijk 1]
    set qk [lindex $qijk 2]
    
    set qw [lindex [$select_atoms get {charge}] $i]
 
    set radius [lindex [$select_atoms get {radius}] $i]
    set orient [lindex [$select_atoms get {user}] $i]
    set normS [expr ($orient + 0.5)/(1.0 + 0.5)]
    
#    puts "Drawing $i / $natoms: quat = $qw $qi $qj $qk; radius = $radius; orient = $orient; normS = $normS"
    draw_ellipsoid_quat $a $b $c $qw $qi $qj $qk $x $y $z $color $resolution
#    draw_arrow_quat $x $y $z $qw $qi $qj $qk $color $resolution $scale
	}
}

# Create a movie from the generated ppm files:
# ffmpeg -i snap%4d.ppm -sameq out.mpg

proc visualtraj {shapex shapey shapez resolution color start end step} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv
  
  set nf [molinfo top get numframes]
  lassign [molinfo top get "{selection 0}"] select_text
  
  set a [expr $shapex / 2.0]
  set b [expr $shapey / 2.0]
  set c [expr $shapez / 2.0]
  
  mesh $a $b $c $resolution
  
  for {set f $start} {$f <= $end} {incr f $step} {
    set mol [molinfo top]
    molinfo $mol set frame $f
    set select_atoms [atomselect top $select_text]
    $select_atoms frame $f
    set natoms  [$select_atoms num]

    graphics top delete all
    
    for {set i 0} {$i < $natoms} {incr i} {  
      # get the atom coordinates
      set coords [lindex [$select_atoms get {x y z}] $i]

      set x [lindex $coords 0]
      set y [lindex $coords 1]
      set z [lindex $coords 2]
      
      # get the quaternion components from file
      set qijk [lindex [$select_atoms get {vx vy vz}] $i]
      set qi [lindex $qijk 0]
      set qj [lindex $qijk 1]
      set qk [lindex $qijk 2]
      
      set qw [lindex [$select_atoms get {charge}] $i]
      
      draw_ellipsoid_quat $a $b $c $qw $qi $qj $qk $x $y $z $color $resolution
    }
    
    set filename snap[format "%04d" $f].ppm
		render TachyonInternal $filename
    
#    rotate y by 1
    
  }
  
}

axes location off

mol delrep 0 top
mol representation VDW 0.33000 16.0000
#mol representation Points
mol color Name
mol material AOEdgy
#mol selection {type 1 8 and ((x-50)*(x-50)+(y-50)*(y-50)+(z-50)*(z-50)<=32*32)}
mol addrep top
color Display Background white
set theta [lindex $argv 0]
visualframe 1 1 1 $theta 1.0 32 green2

display cuedensity 0.20000
display cuemode Exp2 
display shadows on
display ambientocclusion on

pbc box -style tubes -color silver -width 0.8 -center bb

#render Tachyon scene.dat "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -res 1024 1024 -format TGA -res 1024 1024 -o %s.tga"

#set filename snap.tga
#render Tachyon "-aasamples 12 %s -format TARGA -res 1024 1024 -o $filename"
#"/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -o $filename

