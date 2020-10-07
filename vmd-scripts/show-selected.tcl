mol delrep 0 top
mol representation VDW 1.0000 12.0000
mol color Name
mol selection {type 1 2 3 8 and ((x-50)*(x-50)+(y-50)*(y-50)+(z-50)*(z-50)<=40*40)}
mol addrep top

draw sphere {50 50 50} radius 16 resolution 32
