;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; radial grid
;
; n             - number of grid points
; pin, pout     - range of psi
; seps          - locations of separatrices
; sep_factor    - separatrix peaking
; in_dp=in_dp   - Fix the dx on the lower side
; out_dp=out_dp - Fix the dx on the upper side

FUNCTION radial_grid, n, pin, pout, include_in, include_out, seps, sep_factor, $
         in_dp=in_dp, out_dp=out_dp

  IF n EQ 1 THEN BEGIN
    RETURN, [0.5*(pin+pout)]
  ENDIF

  x = FINDGEN(n)
  m = FLOAT(n-1)
  IF NOT include_in THEN BEGIN
    x = x + 0.5
    m = m + 0.5
  ENDIF
  
  IF NOT include_out THEN m = m + 0.5
  x = x / m
  
  IF (NOT KEYWORD_SET(in_dp)) AND (NOT KEYWORD_SET(out_dp)) THEN BEGIN
    ; Neither inner or outer gradients set. Just return equal spacing
    RETURN, pin + (pout - pin)*x
  ENDIF
  
  norm = (x[1] - x[0])*(pout - pin)

  IF KEYWORD_SET(in_dp) AND KEYWORD_SET(out_dp) THEN BEGIN
    ; Fit to dist = a*i^3 + b*i^2 + c*i
    c = in_dp/norm
    b = 3.*(1. - c) - out_dp/norm + c
    a = 1. - c - b
  ENDIF ELSE IF KEYWORD_SET(in_dp) THEN BEGIN
    ; Only inner set
    c = in_dp/norm
    a = 0.5*(c-1.)
    b = 1. - c - a

    ;a = 0
    ;c = in_dp/norm
    ;b = 1. - c
  ENDIF ELSE BEGIN
    ; Only outer set. Used in PF region
    ; Fit to (1-b)*x^a + bx for fixed b
    df = out_dp / norm
    b = 0.25 < df  ; Make sure a > 0
    a = (df - b) / (1. - b)
    vals = pin + (pout - pin)*( (1.-b)*x^a + b*x )
    RETURN, vals
  ENDELSE
  
  vals = pin + (pout - pin)*(c*x + b*x^2 + a*x^3)
  ;STOP
  RETURN, vals
END

