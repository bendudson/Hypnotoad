;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Separatrix line routine
; 
; Finds lines going through a given x-point. Uses this to
; return lines between the x-point and boundary
;
;
; boundary (in)  Optional 2D array [2,n] of points on the boundary

FUNCTION leg_separatrix, dctF, R, Z, xpt_ri, xpt_zi, $
                         opt_ri, opt_zi, $  ; Location of primary O-point
                         status=status, $
                         boundary=boundary
  COMMON td_com, fdata, lastgoodpos
  s = SIZE(dctF, /DIMENSION)
  nr = s[0]
  nz = s[1]

  fdata = dctF

  ; Get second derivatives at the x-point
  
  d = EvalCosP(dctF, x0=xpt_ri, y0=xpt_zi)
  ; extract derivatives (to make more readable)
  f0 = d[0]
  fxx = d[3]
  fyy = d[4]
  fxy = d[5]
  
  xr = INTERPOLATE(r,xpt_ri)
  xz = INTERPOLATE(z,xpt_zi)
  
  drdi = INTERPOLATE(DERIV(R), xpt_ri)
  dzdi = INTERPOLATE(DERIV(Z), xpt_zi)

  ; Use finite-differencing
  ;di = 2.
  ;axp = local_gradient(dctF, xpt_ri + di, xpt_zi)
  ;axm = local_gradient(dctF, xpt_ri - di, xpt_zi)
  ;ayp = local_gradient(dctF, xpt_ri, xpt_zi + di)
  ;aym = local_gradient(dctF, xpt_ri, xpt_zi - di)

  ;fxx = 0.5*(axp.dfdr - axm.dfdr)/di
  ;fyy = 0.5*(ayp.dfdz - aym.dfdz)/di
  ;fxy = 0.25*( (axp.dfdz - axm.dfdz) + (ayp.dfdr - aym.dfdr) ) / di
  
  IF ABS(fyy) GT 1e-4 THEN BEGIN
    ; Get gradients 1 and 2 (solutions y = g1 * x and y = g2 * x)
    
    g1 = ( -fxy + SQRT(fxy^2 - fxx*fyy) ) / fyy
    g2 = ( -fxy - SQRT(fxy^2 - fxx*fyy) ) / fyy
    
    ; Components of vectors

    v1 = [drdi, g1*dzdi]
    v2 = [drdi, g2*dzdi]
  ENDIF ELSE BEGIN
    ; One of the lines through the x-point is vertical (x = const)    
    v1 = [0, 1]
    v2 = [drdi, -fxx / (2.*fxy) * dzdi]
  ENDELSE
  
  ; For each line, work out which direction to go away from the
  ; primary O-point
  
  ; Simple (but possibly error-prone): Dot product with line from
  ; X-point to primary O-point

  v0 = [(opt_ri - xpt_ri)*drdi, (opt_zi - xpt_zi)*dzdi]
  
  ; Normalise to get unit vectors
  v1 = v1 / SQRT(TOTAL(v1^2))
  v2 = v2 / SQRT(TOTAL(v2^2))

  ;oplot, xr+[0,v1[0]], xz+[0,v1[1]], color=2
  ;oplot, xr+[0,v2[0]], xz+[0,v2[1]]
  
  vp = v1+v2
  vm = v1-v2
  vp = vp / SQRT(TOTAL(vp^2))
  vm = vm / SQRT(TOTAL(vm^2))
  
  di = 0.1
  
  dp = TOTAL(v0*vp)
  dm = TOTAL(v0*vm)
  IF ABS(dp) LT ABS(dm) THEN BEGIN
    ; One on core, one PF.
    v2 = -v2
    dp = dm
  ENDIF
  ; Either both pointing along core or both along pf
  IF dp GT 0. THEN BEGIN
    ; Both along core - reverse
    v1 = -v1
    v2 = -v2
  ENDIF
  
  ;oplot, xr+[0,vp[0]], xz+[0,vp[1]], color=2, thick=2
  ;oplot, xr+[0,vm[0]], xz+[0,vm[1]], color=1, thick=2  
  ;oplot, xr+[0,v0[0]], xz+[0,v0[1]], color=3, thick=2
  ;PRINT, v0
  ;PRINT, vp, dp
  ;PRINT, vm, dm
  ;STOP
  
  ; Need to decide which direction in theta this is
  
  dt = theta_differential(0., [xpt_ri + di*v1[0], xpt_zi + di*v1[1]])
  sign = 1.
  IF TOTAL(dt * v1) LT 0 THEN sign = -1.
  
  line1 = theta_line( dctF, $
                      xpt_ri + di*v1[0], xpt_zi + di*v1[1], $
                      sign*di, 1000, boundary=boundary)
  
  core1 = theta_line( dctF, $
                      xpt_ri - di*v1[0], xpt_zi - di*v1[1], $
                      sign*di, 100)

  dt = theta_differential(0., [xpt_ri + di*v2[0], xpt_zi + di*v2[1]])
  sign = 1.
  IF TOTAL(dt * v2) LT 0 THEN sign = -1.

  line2 = theta_line( dctF, $
                      xpt_ri + di*v2[0], xpt_zi + di*v2[1], $
                      sign*di, 1000, boundary=boundary)

  core2 = theta_line( dctF, $
                      xpt_ri - di*v2[0], xpt_zi - di*v2[1], $
                      sign*di, 100)

  OPLOT, INTERPOLATE(R, line1[*,0]), INTERPOLATE(Z, line1[*,1]), color=3, thick=2
  OPLOT, INTERPOLATE(R, line2[*,0]), INTERPOLATE(Z, line2[*,1]), color=4, thick=2
  
  OPLOT, INTERPOLATE(R, core1[*,0]), INTERPOLATE(Z, core1[*,1]), color=3, thick=2
  OPLOT, INTERPOLATE(R, core2[*,0]), INTERPOLATE(Z, core2[*,1]), color=4, thick=2
  
  RETURN, {leg1:line1, leg2:line2, core1:core1, core2:core2, ri:xpt_ri, zi:xpt_zi}
END
