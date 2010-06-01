; Takes the original R-Z data (from G-EQDSK), and the flux mesh
; from create_grid.pro
;
; o Derives additional quantities on the mesh
; o Enforces force-balance
; o Calculates integrated quantities for field-aligned codes
;
; Inputs
; ------
; 
; rz_grid - a structure containing
;    npsigrid - 1D normalised psi grid
;    fpol     - Poloidal current function
;    pres     - Plasma pressure in nt/m^2
;    qpsi     - q values
; 
; mesh - Structure produced by create_grid.pro
;

FUNCTION range, first, last
  IF first LT last THEN BEGIN
    RETURN, first + INDGEN(last - first + 1)
  ENDIF ELSE BEGIN
    RETURN, last + REVERSE(INDGEN(first - last + 1))
  ENDELSE
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Generator of continuous surfaces
;
; First call: 
;   status = gen_surface(mesh=mesh)     - Initialisation
; Subsequent calls
;   yi = gen_surface(period=period, last=last, xi=xi)
;
; period - Set to 1 if the surface is periodic, 0 otherwise
; last   - Set to 1 if this is the last surface
; xi     - The X index of this surface
;
FUNCTION gen_surface, mesh=mesh, period=period, last=last, xi=xi
  COMMON gen_surf_com, m, ys, xind, nd, domain, visited
  IF KEYWORD_SET(mesh) THEN BEGIN
    ; Starting
    m = mesh
    xind = 0 ; Radial surface
    nd = N_ELEMENTS(mesh.npol) ; Number of domains
    domain = 0 ; The domain to start in
    
    ; Running total of npol to get starting y index
    ys = LONARR(nd)
    FOR i=1, nd-1 DO ys[i] = ys[i-1] + mesh.npol[i-1]
    
    ; visited marks which domains have been used
    visited = INTARR(nd)
    
    RETURN, 0
  ENDIF

  IF xind GE TOTAL(m.nrad) THEN BEGIN
    last = 1
    RETURN, 1 ; Error
  ENDIF
  
  ; Get the next surface
  ny = 0
  period = 0 ; Mark as non-periodic
  last = 0 ; mark as not the last
  xi = xind
  REPEAT BEGIN
    IF visited[domain] EQ 1 THEN BEGIN
      ; Already visited this domain
      period = 1 ; Means this domain is periodic
      BREAK
    ENDIF
    
    ; Get the range of indices for this domain
    yi = [range(ys[domain], ys[domain]+m.npol[domain]-1)]
    IF ny EQ 0 THEN yinds = yi ELSE yinds = [yinds, yi]
    ny = ny + m.npol[domain]
    
    visited[domain] = 1 ; Mark domain as visited
    
    ; Find next domain
    IF xind LT m.yup_xsplit[domain] THEN BEGIN
      domain = m.yup_xin[domain]
    ENDIF ELSE BEGIN
      domain = m.yup_xout[domain]
    ENDELSE
  ENDREP UNTIL domain LT 0 ; Keep going until hit a boundary

  ; Find a domain which hasn't been visited
  w = WHERE(visited EQ 0, count)
  IF count NE 0 THEN BEGIN
    ; See if there are any regions with boundaries on lower side
    domain = -1
    FOR i=0, count-1 DO BEGIN
      IF xind LT m.ydown_xsplit[w[i]] THEN BEGIN
        d = m.ydown_xin[w[i]]
      ENDIF ELSE BEGIN
        d = m.ydown_xout[w[i]]
      END
      IF d LT 0 THEN BEGIN
        domain = w[i]
        BREAK
      ENDIF
    ENDFOR
    IF domain LT 0 THEN domain = w[0] ; Set the domain to the first one
    
  ENDIF ELSE BEGIN
    ; No domains left - increase x index (if possible)

    xind = xind + 1
    visited = INTARR(nd) ; Set all to zeros again
    domain = 0 ; Start again with the first domain
    IF xind EQ TOTAL(m.nrad) THEN last = 1 ; No more left
  ENDELSE
  
  IF ny EQ 0 THEN RETURN, 2 ; This shouldn't happen
  
  RETURN, yinds
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Average over flux-surfaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION surface_average, var, mesh
  f = var
  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi)
    f[xi,yi] = MEAN(var[xi,yi]) ; Average over this surface
  ENDREP UNTIL last
  RETURN, f
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Derivatives in X and Y
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; calculates x (psi) derivative for 2D variable
FUNCTION ddx, psi, var
  s = SIZE(var, /dimensions)
  nx = s[0]
  ny = s[1]

  dv = DBLARR(nx, ny)
  FOR i=0, ny-1 DO dv[*,i] = DERIV(psi[*,i], var[*,i])

  RETURN, dv
END

; Take derivative in y, taking into account branch-cuts
FUNCTION ddy, var, mesh
  f = var

  dtheta = 2.*!PI / FLOAT(TOTAL(mesh.npol))

  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    IF period THEN BEGIN
       f[xi,yi] = fft_deriv(var[xi,yi])
    ENDIF ELSE f[xi,yi] = DERIV(var[xi,yi])
  ENDREP UNTIL last
  RETURN, f / dtheta
END

; Integrate a function over y
FUNCTION int_y, var, mesh, loop=loop
  f = var
  
  s = SIZE(var, /dim)
  nx = s[0]
  loop = FLTARR(nx)
  
  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    
    IF period THEN BEGIN
      ; Periodic - use FFT
      f[xi,yi] = fft_integrate(var[xi,yi], loop=lo)
      loop[xi] = lo
    ENDIF ELSE BEGIN
      f[xi,yi] = int_func(var[xi,yi])
    ENDELSE
  ENDREP UNTIL last
  
  RETURN, f
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate f = R * Bt
; Using LSODE method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Function for use by LSODE to integrate Bt
FUNCTION bt_differential, x, Bt
  COMMON bt_com, psi, ax, bx
  a = INTERPOL(ax, psi, x)
  b = INTERPOL(bx, psi, x)
  
  RETURN, a*Bt + b/Bt
END

FUNCTION solve_f, Rxy, psixy, pxy, Bpxy, hthe
  
  MU = 4.e-7*!PI
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]
  
  a = -DDX(psixy, Rxy) / Rxy
  b = -MU*DDX(psixy, pxy) - Bpxy*DDX(Bpxy*hthe)/hthe
  
  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Call LSODE to follow gradient

  ENDIF ELSE BEGIN
     
  ENDELSE
END

FUNCTION force_balance, psixy, Rxy, Bpxy, Btxy, hthe, pxy
  MU =4.e-7*!PI
  
  a = DDX(psixy, Rxy) / Rxy
  b = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe
  
  RETURN, DDX(psixy, Btxy) + a*Btxy + b/Btxy
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate toroidal field
; Using NEWTON method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION Bt_func, Bt
  COMMON fnewt_com, psi, a, b
  
  RETURN, DERIV(psi, Bt) + a*Bt + b / Bt
END

FUNCTION newton_Bt, psixy, Rxy, Btxy, Bpxy, pxy, hthe, mesh
  COMMON fnewt_com, psi, a, b
  MU = 4.e-7*!PI
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]
  
  axy = DDX(psixy, Rxy) / Rxy
  bxy = MU*DDX(psixy, pxy) - Bpxy*DDX(psixy, Bpxy*hthe)/hthe
  
  Btxy2 = FLTARR(nx, ny)
  FOR i=0, ny-1 DO BEGIN
    psi = psixy[*,i]
    a = axy[*,i]
    b = bxy[*,i]
    PRINT, "Solving f for y=", i
    Btxy2[*,i] = NEWTON(Btxy[*,i], "Bt_func")
  ENDFOR
  
  ; Average f over flux surfaces
  fxy = surface_average(Btxy2*Rxy, mesh)

  RETURN, fxy / Rxy
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Solve for pressure and f=R*Bt using force balance
; Using CURVEFIT routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Calculate jxb force at every point, given f on flux-surfaces
; X = X  - dummy, not used
; profiles = [mu0*p(surface), f(surface)]
; force = force(nx, ny), reformed into single vector
PRO jxb_funct, X, profiles, force, pder
  COMMON jxb_com, nx, ny, indxy, psixy, Rxy, hthe, axy

  nsurf = N_ELEMENTS(profiles) / 2
  
  dpdx = profiles[0:(nsurf-1)]
  f = profiles[nsurf:*]

  ; Put profiles into 2D arrays
  fxy = DBLARR(nx, ny)
  dpxy = DBLARR(nx, ny)
  
  FOR x=0, nx-1 DO BEGIN
     FOR y=0, ny-1 DO BEGIN
        i = indxy[x,y]
        fxy[x,y] = f[i]
        dpxy[x,y] = dpdx[i]
     ENDFOR
  ENDFOR

  ; Components of B
  Btxy = fxy / Rxy

  force = Btxy*hthe*DDX(psixy, Btxy) $
    + fxy^2*axy $
    + hthe*dpxy
  
  pder = DBLARR(nx, ny, 2*nsurf)

  ; Diagonal dependencies (dF / dfi)
  dFdfi = (hthe/Rxy)*DDX(psixy, Btxy) $
         + 2.*fxy*axy
  
  ; Set the elements of pder
  FOR x=0, nx-1 DO BEGIN
    FOR y=0, ny-1 DO BEGIN
      ; Get indices into profiles
      xp = x+1 < nx-1
      xm = x-1 > 0
      i = indxy[x,y]
      ip = indxy[xp,y]
      im = indxy[xm,y]
      
      ; f components
      pder[x,y, nsurf+i]  = dFdfi[x,y]
      dx = psixy[xp,y] - psixy[xm,y]
      pder[x,y, nsurf+ip] = pder[x,y, nsurf+ip] + hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xp,y]*dx)
      pder[x,y, nsurf+im] = pder[x,y, nsurf+im] - hthe[x,y]*fxy[x,y]/ (Rxy[x,y]*Rxy[xm,y]*dx)
      
      ; p component
      pder[x,y, i] = hthe[x,y]
    ENDFOR
  ENDFOR
  
  force = REFORM(force, nx*ny)
  pder = REFORM(pder, nx*ny, 2*nsurf)
END

FUNCTION fit_profiles, mesh, psixy, Rxy, hthe, Bpxy, Btxy, dpdx
  COMMON jxb_com, nx, ny, indxy, psi, R, h, axy
  
  MU = 4.e-7*!PI
  
  psi = psixy
  r = Rxy
  h = hthe

  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  ; Map between location in xy and surface number
  indxy = INTARR(nx, ny)
  
  status = gen_surface(mesh=mesh) ; Start generator
  i = 0
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    indxy[xi,yi] = i
    
    IF i EQ 0 THEN BEGIN
      farr = [MEAN(Btxy[xi,yi]*Rxy[xi,yi])]
      parr = [MEAN(dpdx[xi,yi])]
    ENDIF ELSE BEGIN
      farr = [farr, MEAN(Btxy[xi,yi]*Rxy[xi,yi])]
      parr = [parr, MEAN(dpdx[xi,yi])]
    ENDELSE
    
    i = i + 1
  ENDREP UNTIL last
  nsurf = N_ELEMENTS(farr)
  
  profiles = [MU*parr, farr]

  ; Calculate useful quantities
  axy = hthe*DDX(psixy, Rxy)/(Rxy^3)
  
  fit = CURVEFIT(FINDGEN(nx*ny), $
                 REFORM(-Bpxy*DDX(psixy, Bpxy*hthe), nx*ny), $
                 weights, $
                 profiles, $
                 function_name="jxb_funct", /noder)
  
  Btxy2 = FLTARR(nx, ny)
  dpdx2 = FLTARR(nx, ny)
  
  status = gen_surface(mesh=mesh) ; Start generator
  i = 0
  REPEAT BEGIN
    yi = gen_surface(last=last, xi=xi, period=period)
    Btxy2[xi, yi] = profiles[nsurf+i] / Rxy[xi,yi]
    dpdx2[xi, yi] = profiles[i]
    i = i + 1
  ENDREP UNTIL last

  STOP
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Correct hthe using force balance

FUNCTION new_hfunc, h
  COMMON hvars, psi, fixpos, h0, a, b

  IF fixpos EQ 0 THEN BEGIN
    h2 = [h0, h]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    h2 = [h, h0]
  ENDIF ELSE BEGIN
    h2 = [h[0:(fixpos-1)], h0, h[fixpos:*]]
  ENDELSE

  f = a*h2 + b*DERIV(psi, h2)
  
  IF fixpos EQ 0 THEN BEGIN
    f = f[1:*]
  ENDIF ELSE IF fixpos EQ N_ELEMENTS(psi)-1 THEN BEGIN
    f = f[0:(N_ELEMENTS(f)-2)]
  ENDIF ELSE BEGIN
    f = [f[0:(fixpos-1)], f[(fixpos+1):*]]
  ENDELSE

  RETURN, f
END

FUNCTION correct_hthe, Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=fixhthe
  COMMON hvars, xarr, fixpos, h0, a, b
  
  s = SIZE(Rxy, /dim)
  nx = s[0]
  ny = s[1]

  MU = 4.e-7*!PI

  IF NOT KEYWORD_SET(fixhthe) THEN fixhthe = 0
  IF fixhthe LT 0 THEN fixhthe = 0
  IF fixhthe GT nx-1 THEN fixhthe = nx-1

  fixpos = fixhthe
  PRINT, "FIX = ", fixhthe
  
  axy = Btxy*DDX(psixy, Btxy) + Bpxy*DDX(psixy, Bpxy) $
    + Btxy^2*DDX(psixy, Rxy)/Rxy + MU*DDX(psixy, pressure)
  bxy = Bpxy^2

  nh = DBLARR(nx, ny)
  nh[fixhthe,*] = hthe[fixhthe,*]
  FOR i=0, ny-1 DO BEGIN
    PRINT, "Correcting y index ", i
    xarr = psixy[*,i]
    a = axy[*,i]
    b = bxy[*,i]
    h0 = hthe[fixhthe,i]
    
    IF fixhthe EQ 0 THEN BEGIN
      htmp = REFORM(hthe[1:*,i])
    ENDIF ELSE IF fixhthe GE nx-1 THEN BEGIN
      ; fix last point
      htmp = hthe[0:(nx-2),i]
      fixhthe = nx-1
    ENDIF ELSE BEGIN
      ; fix somewhere in the middle  
      htmp = [hthe[0:(fixhthe-1),i], hthe[(fixhthe+1):*,i]]
    ENDELSE
    
    htmp = NEWTON(htmp, "new_hfunc")
    
    IF fixhthe EQ 0 THEN BEGIN
      nh[1:*] = htmp
    ENDIF ELSE IF fixhthe GE nx-1 THEN BEGIN
      nh[0:(nx-2), i] = htmp
    ENDIF ELSE BEGIN
      nh[0:(fixhthe-1), i] = htmp[0:(fixhthe-1)]
      nh[(fixhthe+1):*, i]  = htmp[fixhthe:*]
    ENDELSE
    
    w = WHERE(nh[*,i] LT 0.0, count)
    IF count GT 0 THEN BEGIN
      PRINT, "Error in hthe solver: Negative solution at y = ", i
      ;STOP
    ENDIF
  ENDFOR

  RETURN, nh
END

PRO process_grid, rz_grid, mesh, output=output, poorquality=poorquality
  
  MU = 4.e-7*!PI

  poorquality = 0

  IF NOT KEYWORD_SET(output) THEN output="bout.grd.nc"
  
  ; Size of the mesh
  nx = TOTAL(mesh.nrad)
  ny = TOTAL(mesh.npol)


  Rxy = mesh.Rxy
  Zxy = mesh.Zxy
  psixy = mesh.psixy*mesh.fnorm + mesh.faxis ; Non-normalised psi

  pressure = INTERPOL(rz_grid.pres, rz_grid.npsigrid, mesh.psixy)
  w = WHERE(mesh.psixy GT MAX(rz_grid.npsigrid), count)
  IF count GT 0 THEN pressure[w] = 0.

  ; NOTE: WHAT ABOUT PF REGIONS?

  IF MIN(pressure) LT 0.0 THEN BEGIN
    PRINT, ""
    PRINT, "============= WARNING =============="
    PRINT, "Poor quality equilibrium: Pressure is negative"
    PRINT, ""
    poorquality = 1
  ENDIF
  
  w = WHERE(mesh.psixy GE 1.0, count)
  IF count GT 0 THEN pressure[w] = 0.0
  
  dpdpsi = DDX(psixy, pressure)

  IF MAX(dpdpsi)*mesh.fnorm GT 0.0 THEN BEGIN
    PRINT, ""
    PRINT, "============= WARNING =============="
    PRINT, "Poor quality equilibrium: Pressure is increasing radially"
    PRINT, ""
    poorquality = 1
  ENDIF

  ; Grid spacing
  dx = FLTARR(nx, ny)
  dx[*,0] = DERIV(psixy[*,0])
  FOR y=1, ny-1 DO dx[*,y] = dx[*,0]
  
  dtheta = 2.*!PI / FLOAT(ny)
  dy = FLTARR(nx, ny) + dtheta
  
  ; B field components
  
  Bpxy = SQRT(mesh.dpsidR^2 + mesh.dpsidZ^2) / Rxy
  
  ; Determine direction (dot B with grad y vector)
  
  ; Get toroidal field from poloidal current function fpol
  Btxy = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
    FOR j=0, ny-1 DO BEGIN
      IF mesh.psixy[i,j] GE 1.0 THEN BEGIN
        fpol = rz_grid.fpol[N_ELEMENTS(rz_grid.fpol)-1]
      ENDIF ELSE BEGIN
        fpol = INTERPOL(rz_grid.fpol, rz_grid.npsigrid, mesh.psixy[i,j])
      ENDELSE
      
      Btxy[i,j] = fpol / Rxy[i,j]
    ENDFOR
  ENDFOR

  ; Total B field
  Bxy = SQRT(Btxy^2 + Bpxy^2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Go through the domains to get a starting estimate
  ; of hthe
  hthe = FLTARR(nx, ny)

  m = MAX(Rxy[0, *], ymidplane) ; Pick a midplane index
  
  status = gen_surface(mesh=mesh) ; Start generator
  REPEAT BEGIN
    ; Get the next domain
    yi = gen_surface(period=period, last=last, xi=xi)
    
    ; Get distance along this line
    IF period THEN BEGIN
      ; Periodic, so can use FFT
      drdi = REAL_PART(fft_deriv(Rxy[xi, yi]))
      dzdi = REAL_PART(fft_deriv(Zxy[xi, yi]))
    ENDIF ELSE BEGIN
      ; Non-periodic
      drdi = DERIV(Rxy[xi, yi])
      dzdi = DERIV(Zxy[xi, yi])
    ENDELSE
    
    dldi = REFORM(SQRT(drdi^2 + dzdi^2))
    
    ; Need to smooth to get sensible results
    IF period THEN BEGIN
      n = N_ELEMENTS(dldi)
      dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
      dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
      dldi = (SMOOTH([ dldi[(n-2):*], dldi, dldi[0:1] ], 5))[2:(n+1)]
    ENDIF ELSE BEGIN
      dldi = SMOOTH(dldi, 5)
      dldi = SMOOTH(dldi, 5)
      dldi = SMOOTH(dldi, 5)
    ENDELSE
    
    hthe[xi, yi] = dldi / dtheta ; First estimate of hthe
    
    ; Get outboard midplane
    IF period AND xi EQ 0 THEN BEGIN
      m = MAX(Rxy[0,yi], ymidplane)
    ENDIF
  ENDREP UNTIL last

  PRINT, "Midplane index ", ymidplane

  fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pressure)
  PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Correct pressure using hthe
  
  PRINT, "Calculating pressure profile from force balance"

  ; Calculate force balance
  dpdx = ( -DDX(psixy, Bxy^2 * hthe) + hthe*Bpxy*DDX(psixy, Bpxy) + Btxy*Rxy*DDX(psixy, Btxy*hthe/Rxy) ) / (MU*hthe)
  
  ; Surface average
  dpdx2 = surface_average(dpdx, mesh)
  
  pres = FLTARR(nx, ny)
  ; Integrate to get pressure
  FOR i=0, ny-1 DO BEGIN
    pres[*,i] = int_func(psixy[*,i], dpdx2[*,i])
    pres[*,i] = pres[*,i] - pres[nx-1,i]
  ENDFOR
  
  w = WHERE(pres LT 0., count)
  IF count GT 0 THEN pres[w] = 0.
  
  ; Some sort of smoothing here?
  
  fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, hthe, pres)
  PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
  
  !P.MULTI=[0,0,2,0,0]
  SURFACE, pressure, xtitle="X", ytitle="Y", title="Input pressure", chars=2
  SURFACE, pres, xtitle="X", ytitle="Y", title="New pressure", chars=2
  
  IF get_yesno("Keep new pressure?") THEN BEGIN
    pressure = pres
    dpdpsi = dpdx2
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Correct f = RBt using force balance

  IF get_yesno("Correct f=RBt using force balance?") THEN BEGIN

    new_Btxy = newton_bt(psixy, Rxy, Btxy, Bpxy, pres, hthe, mesh)
    
    fb0 = force_balance(psixy, Rxy, Bpxy, new_Btxy, hthe, pressure)
    PRINT, "force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
    
    !P.MULTI=[0,0,2,0,0]
    SURFACE, Btxy, xtitle="X", ytitle="Y", title="Input Bt", chars=2
    SURFACE, new_Btxy, xtitle="X", ytitle="Y", title="New Bt", chars=2

    IF get_yesno("Keep new Bt?") THEN BEGIN
      Btxy = new_Btxy
      Bxy = SQRT(Btxy^2 + Bpxy^2)
    ENDIF
  ENDIF
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE HTHE
  ; Modify hthe to fit force balance using initial guess
  ; Does not depend on signs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF get_yesno("Adjust hthe using force balance?") THEN BEGIN
    ; This doesn't behave well close to the x-points
    fixhthe = FIX(nx / 2)
    nh = correct_hthe(Rxy, psixy, Btxy, Bpxy, hthe, pressure, fixhthe=fixhthe)
    
    fb0 = force_balance(psixy, Rxy, Bpxy, Btxy, nh, pressure)
    PRINT, "Force imbalance: ", MEAN(ABS(fb0)), MAX(ABS(fb0))
    
    PRINT, "Maximum difference in hthe: ", MAX(ABS(hthe - nh))
    PRINT, "Maximum percentage difference: ", 100.*MAX(ABS((hthe - nh)/hthe))

    !P.multi=[0,0,1,0,0]
    PLOT, hthe[*,0], title="Poloidal arc length at midplane. line is initial estimate", color=1
    OPLOT, nh[*,0], psym=1, color=2
    OPLOT, nh[*,0], color=2

    IF get_yesno("Keep new hthe?") THEN BEGIN
      hthe = nh
    ENDIF
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; CALCULATE PARALLEL CURRENT
  ; Provides a way to check if Btor should be reversed
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  PRINT, "Checking parallel current"

  jpar0 = Bxy * DDX(psixy, Btxy*Rxy) / MU + Rxy * dpdpsi

  j0 = ((Bpxy*Btxy*Rxy/(Bxy*hthe))*( DDX(psixy, Bxy^2*hthe/Bpxy) - Btxy*Rxy*DDX(psixy,Btxy*hthe/(Rxy*Bpxy)) ) $
        - Bxy*DDX(psixy, Btxy*Rxy)) / MU
  
  IF (MEAN(ABS(j0 + jpar0)) LT MEAN(ABS(j0 - jpar0))) THEN BEGIN
    PRINT, "****Equilibrium has -ve toroidal field"
    Btxy = -Btxy
    j0 = -j0
  ENDIF
  PRINT, "Maximum difference in jpar0: ", MAX(ABS(j0[2:*,*] - jpar0[2:*,*]))
  PRINT, "Maximum percentage difference: ", 200.*MAX(ABS((jpar0[2:*,*] - j0[2:*,*])/(jpar0[2:*,*]+j0[2:*,*])))

  !P.MULTI=[0,0,2,0,0]
  SURFACE, jpar0, xtitle="X", ytitle="Y", title="Jpar from f and p profiles", chars=2
  SURFACE, j0, xtitle="X", ytitle="Y", title="Jpar from B field", chars=2
  
  IF get_yesno("Use B field jpar?") THEN BEGIN
    jpar0 = j0
  ENDIF
  
  ; Calculate field-line pitch
  pitch = hthe * Btxy / (Bpxy * Rxy)
  
  ; derivative with psi
  dqdpsi = DDX(psixy, pitch)
  
  qinty = int_y(pitch, mesh, loop=qloop) * dtheta
  qloop = qloop * dtheta
  sinty = int_y(dqdpsi, mesh) * dtheta

  ; NOTE: This is only valid in the core
  pol_angle = FLTARR(nx,ny)
  FOR i=0, nx-1 DO pol_angle[i, *] = 2.0*!PI * qinty[i,*] / qloop[i]

  ;;;;;;;;;;;;;;;;;;;; THETA_ZERO ;;;;;;;;;;;;;;;;;;;;;;
  ; re-set zshift to be zero at the outboard midplane
  qm = qinty[*,ymidplane]
  sm = sinty[*,ymidplane]
  
  FOR i=0, ny-1 DO BEGIN
     qinty[*,i] = qinty[*,i] - qm
     sinty[*,i] = sinty[*,i] - sm
  ENDFOR
  
  ;;;;;;;;;;;;;;;;;;;; CURVATURE ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculating b x kappa
  
  dpb = DBLARR(nx, ny)      ; quantity used for y and z components
      
  FOR i=0, ny-1 DO BEGIN
    dpb[*,i] = MU*dpdpsi/Bxy[*,i]
  ENDFOR
  dpb = dpb + DDX(psixy, Bxy)
  
  bxcvx = Bpxy * DDY(Btxy*Rxy / Bxy, mesh) / hthe
  bxcvy = Bpxy*Btxy*Rxy*dpb / (hthe*Bxy^2)
  bxcvz = -dpb - sinty*bxcvx 
  
  ;;;;;;;;;;;;;;;;;;;; TOPOLOGY ;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate indices for backwards-compatibility
  
  IF (N_ELEMENTS(mesh.nrad) EQ 2) AND (N_ELEMENTS(mesh.npol) EQ 3) THEN BEGIN
    PRINT, "Single null equilibrium"
    
    ixseps1 = mesh.nrad[0]
    ixseps2 = nx
    
    jyseps1_1 = mesh.npol[0]-1
    jyseps1_2 = mesh.npol[0] + FIX(mesh.npol[1]/2)
    ny_inner = jyseps1_2
    jyseps2_1 = jyseps1_2
    jyseps2_2 = ny - mesh.npol[2]-1

  ENDIF ELSE IF (N_ELEMENTS(mesh.nrad) EQ 3) AND (N_ELEMENTS(mesh.npol) EQ 6) THEN BEGIN
    PRINT, "Double null equilibrium"
    
    ixseps1 = mesh.nrad[0]
    ixseps2 = ixseps1 + mesh.nrad[1]
    
    jyseps1_1 = mesh.npol[0]-1
    jyseps2_1 = jyseps1_1 + mesh.npol[1]
    
    ny_inner = jyseps2_1 + mesh.npol[2] + 1
    
    jyseps1_2 = ny_inner + mesh.npol[3] - 1
    jyseps2_2 = jyseps1_2 + mesh.npol[4]
    
  ENDIF ELSE BEGIN
    PRINT, "WARNING: Equilibrium not recognised."
    ixseps1 = -1
    ixseps2 = -1
    
    jyseps1_1 = -1
    jyseps1_2 = FIX(ny/2)
    jyseps2_1 = FIX(ny/2)
    ny_inner = FIX(ny/2)
    jyseps2_2 = ny
  ENDELSE

  PRINT, "Generating plasma profiles:"
          
  PRINT, "  1. Flat temperature profile"
  PRINT, "  2. Flat density profile"
  PRINT, "  3. Te proportional to density"
  REPEAT BEGIN
    opt = get_integer("Profile option:")
  ENDREP UNTIL (opt GE 1) AND (opt LE 3)
  
  IF opt EQ 1 THEN BEGIN
    ; flat temperature profile
    
    PRINT, "Setting flat temperature profile"
    REPEAT BEGIN
      Te_x = get_float("Temperature (eV):")
      
      
      ; get density
      Ni = pressure / (2.*Te_x* 1.602e-19*1.0e20)
      
      PRINT, "Maximum density (10^20 m^-3):", MAX(Ni)
      
      done = get_yesno("Is this ok?")
    ENDREP UNTIL done EQ 1
    
    Te = FLTARR(nx, ny)+Te_x
    Ti = Te
    Ni_x = MAX(Ni)
    Ti_x = Te_x
  ENDIF ELSE IF opt EQ 2 THEN BEGIN
    PRINT, "Setting flat density profile"
    
    REPEAT BEGIN
      ni_x = get_float("Density [10^20 m^-3]:")
      
      ; get temperature
      Te = pressure / (2.*ni_x* 1.602e-19*1.0e20)
      
      PRINT, "Maximum temperature (eV):", MAX(Te)
    ENDREP UNTIL get_yesno("Is this ok?") EQ 1
    
    Ti = Te
    Ni = FLTARR(nx, ny) + ni_x
    Te_x = MAX(Te)
    Ti_x = Te_x
  ENDIF ELSE BEGIN
    PRINT, "Setting te proportional to density"
    
    REPEAT BEGIN
      te_x = get_float("Maximum temperature [eV]:")
      
      ni_x = max(pressure) / (2.*Te_x* 1.602e-19*1.0e20)
      
      PRINT, "Maximum density [10^20 m^-3]:", ni_x
      
      Te = te_x * pressure / max(pressure)
      Ni = ni_x * pressure / max(pressure)
    ENDREP UNTIL get_yesno("Is this ok?") EQ 1
    Ti = Te
    Ti_x = Te_x
  ENDELSE
  
  rmag = MAX(ABS(Rxy))
  PRINT, "Setting rmag = ", rmag
  
  bmag = MAX(ABS(Bxy))
  PRINT, "Setting bmag = ", bmag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; save to file

  ;STOP

  PRINT, "Writing grid to file "+output

  handle = file_open(output, /CREATE)

  ; Size of the grid

  s = file_write(handle, "nx", nx)
  s = file_write(handle, "ny", ny)

  ; Topology for original scheme
  s = file_write(handle, "ixseps1", ixseps1)
  s = file_write(handle, "ixseps2", ixseps2)
  s = file_write(handle, "jyseps1_1", jyseps1_1)
  s = file_write(handle, "jyseps1_2", jyseps1_2)
  s = file_write(handle, "jyseps2_1", jyseps2_1)
  s = file_write(handle, "jyseps2_2", jyseps2_2)
  
  ; Topology for general configurations
  s = file_write(handle, "yup_xsplit", mesh.yup_xsplit)
  s = file_write(handle, "ydown_xsplit", mesh.ydown_xsplit)
  s = file_write(handle, "yup_xin", mesh.yup_xin)
  s = file_write(handle, "yup_xout", mesh.yup_xout)
  s = file_write(handle, "ydown_xin", mesh.ydown_xin)
  s = file_write(handle, "ydown_xout", mesh.ydown_xout)

  ; Grid spacing
  
  s = file_write(handle, "dx", dx)
  s = file_write(handle, "dy", dy)
  
  s = file_write(handle, "ShiftAngle", qloop)
  s = file_write(handle, "zShift", qinty)
  s = file_write(handle, "pol_angle", pol_angle)
  s = file_write(handle, "ShiftTorsion", dqdpsi)

  s = file_write(handle, "Rxy",  Rxy)
  s = file_write(handle, "Zxy",  Zxy)
  s = file_write(handle, "Bpxy", Bpxy)
  s = file_write(handle, "Btxy", Btxy)
  s = file_write(handle, "Bxy",  Bxy)
  s = file_write(handle, "hthe", hthe)
  s = file_write(handle, "sinty", sinty)
  s = file_write(handle, "psixy", psixy)

  ; plasma profiles

  s = file_write(handle, "pressure", pressure)
  s = file_write(handle, "Jpar0", Jpar0)
  s = file_write(handle, "Ni0", Ni)
  s = file_write(handle, "Te0", Te)
  s = file_write(handle, "Ti0", Ti)
  s = file_write(handle, "Ni_x", Ni_x)
  s = file_write(handle, "Te_x", Te_x)
  s = file_write(handle, "Ti_x", Ti_x)
  s = file_write(handle, "bmag", bmag)
  s = file_write(handle, "rmag", rmag)

  ; Curvature
  s = file_write(handle, "bxcvx", bxcvx)
  s = file_write(handle, "bxcvy", bxcvy)
  s = file_write(handle, "bxcvz", bxcvz)

  ; Psi range
  s = file_write(handle, "psi_axis", mesh.faxis)
  psi_bndry = mesh.faxis + mesh.fnorm
  s = file_write(handle, "psi_bndry", psi_bndry)

  file_close, handle
  PRINT, "DONE"
  
  !P.multi=[0,0,1,0,0]

  ;STOP
END
