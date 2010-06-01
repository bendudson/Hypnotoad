;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Line crossing detection
;  
; (r1, z1) and (r2, z2) are two lines
; period1 and period2 determine whether the lines are periodic
; 
; Returns a list of locations where they intersect

FUNCTION line_crossings, r1, z1, period1, r2, z2, period2, ncross=ncross, $
                         inds1=inds1, inds2=inds2
  n1 = N_ELEMENTS(r1)
  n2 = N_ELEMENTS(r2)
  
  result = 0
  ncross = 0

  FOR i=0, n1-1 DO BEGIN
    ip = i + 1
    IF i EQ n1-1 THEN BEGIN
      IF period1 THEN ip = 0 ELSE BREAK
    ENDIF
    
    FOR j=0, n2-1 DO BEGIN
      jp = j+1
      IF j EQ n2-1 THEN BEGIN
        IF period2 THEN jp = 0 ELSE BREAK
      ENDIF
      
      ; Test if line (i to ip) and (j to jp) intersects
      ; cast as a 2x2 matrix
      
      a = r1[ip] - r1[i]
      b = r2[j] - r2[jp]
      c = z1[ip] - z1[i]
      d = z2[j] - z2[jp]
      
      dr = r2[j] - r1[i]
      dz = z2[j] - z1[i]

      det = a*d - b*c
      
      ; Get location along the line segments
      IF ABS(det) GT 1.e-6 THEN BEGIN
        alpha = (d*dr - b*dz)/det
        beta =  (a*dz - c*dr)/det
      ENDIF ELSE BEGIN
        alpha = -1.
        beta = -1.
      ENDELSE
      
      IF (alpha GE 0.0) AND (alpha LE 1.0) AND (beta GE 0.0) AND (beta LE 1.0) THEN BEGIN
        ; Intersection
        
        r = r1[i] + alpha * a
        z = z1[i] + alpha * c
        
        IF ncross EQ 0 THEN BEGIN
          result = FLTARR(2,1)
          result[0,0] = r
          result[1,0] = z
          
          inds1 = [FLOAT(i)+alpha]
          inds2 = [FLOAT(j)+beta]
        ENDIF ELSE BEGIN
          rold = result
          result = FLTARR(2, ncross+1)
          result[*,0:(ncross-1)] = rold
          result[0,ncross] = r
          result[1,ncross] = z

          inds1 = [inds1, FLOAT(i)+alpha]
          inds2 = [inds2, FLOAT(j)+beta]
        ENDELSE
        ncross = ncross + 1
      ENDIF
    ENDFOR
  ENDFOR
  
  RETURN, result
END
