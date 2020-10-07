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
# 4- Invoke the draw function: visualframe 3 1 1 12 1
# where 3 1 1 are the lengths of the ellipsoids in x, y and z
#       12 is the resolution
#       1 is the coloring mode:
#            0 = use read-in variables (c_orient) in the dump file
#            1 = use the calculated value for the orientation of the ellipsoid x axis wrt to the global z axis


proc f {a u v} {
  expr {$a*sin($u) * cos($v)} 
}  

proc g {b u v} {
  expr {$b*cos($u) * cos($v)}
}

proc h {c u v} {
  expr {$c*sin($v)} 
}

proc mesh {a b c res} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv

  set PI 3.1415926
  set minu 0
  set resu $res
  set maxu [expr $PI * 2]
  set stepu [expr $maxu / $resu]

  set minv [expr -1*$PI]
  set maxv [expr $PI]
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

proc draw_ellipsoid_quat { a b c iquatw iquati iquatj iquatk x y z orient eColour res} {
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

  if {$coloringMode == 0} {
    set normS $orient
  } else {
   set S [expr (3.0 * $a31 * $a31 - 1.0)/2.0]
   set normS [expr ($S + 0.5)/(1.0 + 0.5)]
  }

  draw color [expr $minColor + int($normS * $NrColors)]

#  draw color $eColour

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

# a = semi-axis in x direction
# b = semi-axis in y direction
# c = semi-axis in z direction
# phi = first rotation about z axis in xy plane  [see http://mathworld.wolfram.com/EulerAngles.html]
# theta = second rotation about x axis in zy plane.
# psi = third rotation about z' axis in x'y' plane.
# x = translation in x direction.
# y = translation in y direction
# z = translation in z direction
# eColour is the colour of the ellipse in VMD colour table.
# res = resolution

proc draw_ellipsoid_euler { a b c phi theta psi x y z eColour res} {
  set PI 3.1415926
  set phi [expr $phi*$PI/180]
  set theta [expr $theta*$PI/180]
  set psi [expr $psi*$PI/180]

  set resu $res
  set resv $res

  set minu 0
  set maxu [expr $PI * 2]
  set stepu [expr $maxu / $resu]

  set minv [expr -1*$PI]
  set maxv [expr $PI]
  set stepv [expr ($maxv - $minv) / $resv]

  # first, get the data (this isn't the most data efficient way of   
  # doing things)   
  for {set iu 0} {$iu <= $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv <= $resv} {set iv [expr $iv + 1]} {
      set u [expr $minu + $iu * $stepu]
      set v [expr $minv + $iv * $stepv]
      set fdata($u,$v) [f $a $u $v]
      set gdata($u,$v) [g $b $u $v]
      set hdata($u,$v) [h $c $u $v]     
    }
  }

  # compute rotation matrix to get desired orientation
  set a11 [expr {cos($psi)*cos($phi)-cos($theta)*sin($phi)*sin($psi)}]
  set a12 [expr {cos($psi)*sin($phi)+cos($theta)*cos($phi)*sin($psi)}]
  set a13 [expr {sin($psi)*sin($theta)}]
  set a21 [expr {-1*sin($psi)*cos($phi)-cos($theta)*sin($phi)*cos($psi)}]
  set a22 [expr {-1*sin($psi)*sin($phi)+cos($theta)*cos($phi)*cos($psi)}]
  set a23 [expr {cos($psi)*sin($theta)}]
  set a31 [expr {sin($theta)*sin($phi)}]
  set a32 [expr {-1*sin($theta)*cos($phi)}]
  set a33 [expr {cos($theta)}]

  #perform rotation and translation
  for {set iu 0} {$iu <= $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv <= $resv} {set iv [expr $iv + 1]} {
      set u [expr $minu + $iu * $stepu]
      set v [expr $minv + $iv * $stepv]
      set fdata_r($iu,$iv) [expr {($fdata($u,$v)*$a11+$gdata($u,$v)*$a21+$hdata($u,$v)*$a31) + $x}]       
      set gdata_r($iu,$iv) [expr {($fdata($u,$v)*$a12+$gdata($u,$v)*$a22+$hdata($u,$v)*$a32) + $y}]
      set hdata_r($iu,$iv) [expr {($fdata($u,$v)*$a13+$gdata($u,$v)*$a23+$hdata($u,$v)*$a33) + $z}]
    }
  }

  # make another pass through to plot it   
  set cnum [colorinfo num]
  set cmax [colorinfo max]

  for {set iu 0} {$iu < $resu} {set iu [expr $iu + 1]} {
    for {set iv 0} {$iv < $resv} {set iv [expr $iv + 1]} {

      # draw color [expr $cnum + (int($cmax * $u / $maxu) + int($cmax * $v / $maxv)) % $cmax]
      draw color $eColour

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

# shapex = dimension in x direction
# shapey = dimension in y direction
# shapez = dimension in z direction
# resolution
# mode

proc visualframe {shapex shapey shapez resolution mode} {
  global fdata gdata hdata minu maxu minv maxv stepu stepv resu resv
  global minColor maxColor NrColors colorScheme color coloringMode

  set colorScheme BGryR
  color scale method $colorScheme
  color scale min 0
  color scale max 1
  color scale midpoint 0.5
  
  set coloringMode $mode
  set color red
  
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

  puts "Drawing $natoms ellipsoids"
    
  set a [expr $shapex / 2.0]
  set b [expr $shapey / 2.0]
  set c [expr $shapez / 2.0]
  
  mesh $a $b $c $resolution

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
    
    puts "Drawing $i / $natoms: quat = $qw $qi $qj $qk; radius = $radius; orient = $orient; normS = $normS"
    draw_ellipsoid_quat $a $b $c $qw $qi $qj $qk $x $y $z $normS $color $resolution 
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
  
  set normS 1
  
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
      
      draw_ellipsoid_quat $a $b $c $qw $qi $qj $qk $x $y $z $normS $color $resolution
    }
    
    set filename snap[format "%04d" $f].ppm
		render TachyonInternal $filename
    
#    rotate y by 1
    
  }
  
}


