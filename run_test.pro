
g = read_neqdsk("efit/neqdsk")

R = REFORM(g.r[*,0])
Z = REFORM(g.z[0,*])

boundary=fltarr(2,4)
boundary[0,*] = [1.0, 1.0, 2.5, 2.5]
boundary[1,*] = [-1.4, 1.4, 1.4, -1.4]

;boundary = TRANSPOSE([[g.xlim], [g.ylim]])

mesh = create_grid(g.psi, R, Z, boundary=boundary, /strict)

; Need to get pressure and f on this mesh
;pressure = INTERPOL(g.pres, 
