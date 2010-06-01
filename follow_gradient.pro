;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Follow the gradient from a given point to a target f
; 

; Calculate dR/df and dZ/df for use by LSODE
; Input: pos[0] = R, pos[1] = Z
; Output [0] = dR/df = -Bz/B^2 , [1] = dZ/df = Br/B^2
FUNCTION radial_differential, fcur, pos
  COMMON rd_com, dctf, lastgoodf, lastgoodpos, R, Z
  a = local_gradient(dctf, pos[0], pos[1], status=status)
  
  IF status EQ 0 THEN BEGIN
    ; No error in localgradient.
    lastgoodf = fcur 
    lastgoodpos = pos
    
    ; If status NE 0 then an error occurred.
    ; Allow a.dfdz to cause an error so escape LSODE
  ENDIF
  
  dRdi = INTERPOLATE(DERIV(R), pos[0])
  dZdi = INTERPOLATE(DERIV(Z), pos[1])
  
  ; Check mismatch between fcur and a.f ? 
  Br = a.dfdz/dZdi
  Bz = -a.dfdr/dRdi
  B2 = Br^2 + Bz^2

  RETURN, [-Bz/B2/dRdi, Br/B2/dZdi]
END

;
; F         (in)    2D psi data
; R, Z      (in)    1D arrays of position vs. index
; ri0,zi0   (in)    Starting indices
; ftarget   (in)    The f to aim for
; ri,zi     (out)   Final position
; status    (out)   Non-zero if hits a boundary. 1 if out of range at
;                   the start, 2 if before reaching ftarget
; boundary  (in)    Optional 2D array [2,n] of points on the boundary
; fbndry    (out)   If hits boundary, gives final f value
; ibndry    (out)   If hits boundary, index where hit
;
PRO follow_gradient, dctF, R, Z, ri0, zi0, ftarget, ri, zi, status=status, $
                     boundary=boundary, fbndry=fbndry, ibndry=ibndry
  COMMON rd_com, df, lastgoodf, lastgoodpos, Rpos, Zpos
  
  Rpos = R
  Zpos = Z
  
  ibndry = -1

  df = dctF
  
  ; Get starting f
  g = local_gradient(dctF, ri0, zi0, status=status)
  IF status EQ 1 THEN BEGIN
    ri = ri0
    zi = zi0
    status = 1
    RETURN
  ENDIF
  f0 = g.f

  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Call LSODE to follow gradient
    rzold = [ri0, zi0]
    rznew = LSODE(rzold,f0,ftarget - f0,'radial_differential')
    ri = rznew[0]
    zi = rznew[1]
  ENDIF ELSE BEGIN
    ; An error occurred in LSODE.
    ; lastgoodf contains the last known good f value
    CATCH, /cancel
    status = 2
    fbndry = lastgoodf
    ri = lastgoodpos[0]
    zi = lastgoodpos[1]
    PRINT, "Out of domain at f = ", fbndry
    
    ; Repeat to verify that this does work
    rzold = [ri0, zi0]
    CATCH, theError
    fbndry = lastgoodf - 0.1*(ftarget - f0)
    IF theError NE 0 THEN BEGIN
      PRINT, "   Error again at ", fbndry
    ENDIF
    rznew = LSODE(rzold,f0,fbndry - f0,'radial_differential')
    
    RETURN
  ENDELSE
  CATCH, /cancel
  
  IF KEYWORD_SET(boundary) THEN BEGIN
    ; Check if the line crossed a boundary
    ;PRINT, "Checking boundary ", boundary[*,1:2], [ri0, ri], [zi0, zi]
    cpos = line_crossings([ri0, ri], [zi0, zi], 0, $
                          boundary[0,*], boundary[1,*], 1, ncross=ncross, inds2=inds2)
    IF (ncross MOD 2) EQ 1 THEN BEGIN ; Odd number of boundary crossings
      ; Use divide-and-conquer to find crossing point
      
      ibndry = inds2[0] ; Index in boundary where hit
      
      fmax = ftarget ; Target (with boundary in the way)
      fcur = f0 ; Current known good position
      rzold = [ri0,zi0]
      rzcur = rzold
      REPEAT BEGIN
        fbndry = (fcur + fmax) / 2
        ; Try to go half-way to fmax
        rznew = LSODE(rzold,f0,fbndry-f0,'radial_differential')
        cpos = line_crossings([rzold[0], rznew[0]], [rzold[1], rznew[1]], 0, $
                              boundary[0,*], boundary[1,*], 1, ncross=ncross, inds2=inds2)
        IF (ncross MOD 2) EQ 0 THEN BEGIN
          ; Didn't cross. Make this the new current location
          fcur = fbndry
          rzcur = rznew
        ENDIF ELSE BEGIN
          ; Crossed boundary. Change fmax
          fmax = fbndry
          ibndry = inds2[0] ; refined boundary index
        ENDELSE
      ENDREP UNTIL ABS(fmax - fcur) LT 0.01*ABS(ftarget - f0)
      ri = rzcur[0]
      zi = rzcur[1]
      fbndry = fcur
      
      ;PRINT, "Hit boundary", ri, zi, " f =", fbndry
      status = 2
      RETURN
    ENDIF
  ENDIF

  status = 0
END
