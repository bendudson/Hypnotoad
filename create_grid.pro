; Tokamak grid generator
; ======================
; 
; Generates a flux-surface aligned grid from
; an R-Z mesh of psi values. 
;
; Features:
; --------
;  o An arbitrary number of X-points
;  o Automatic default settings when not
;    supplied. 
;
; Author: Ben Dudson, University of York, Nov 2009
; 
;
; NOTE: Throughout, "F" means un-normalised psi,
;       and "psi" means normalised psi
;
; Useage:
; ------
; 
; FUNCTION create_grid, F, R, Z, settings, critical=critical,
;                                boundary=boundary, fpsi=fpsi
;
; F - psi(nr, nz) 2D array
; R - R(nr)  1D array
; Z - Z(nz)  1D array
; 
; Settings is a structure containing:
;   
;   psi_inner, psi_outer  - Range of normalised psi
;   nrad                  - Number of radial grid points
;                           Scalar -> Total number. Distributed
;                                     automatically
;                           Array -> Specified for each section
;   rad_peaking           - Radial separatrix peaking factor
;                           Not supplied -> 0
;                           scalar -> same for all regions
;                           array -> Different for each region
;   npol                  - Number of poloidal points.
;                           Scalar -> Total number. Distributed
;                                     automatically
;                           Array -> Specified for each section
;   pol_peaking           - Poloidal peaking factor
;                           Not supplied -> 0
;                           Scalar -> Same for all poloidal sections
;                           Array -> Specified for each section
; 
; A minimal settings structure must contain 4 scalars:
;    psi_inner, psi_outer, nrad and npol.
;
; Critical is a structure containing critical points:
;  
;   n_opoint, n_xpoint   - Number of O- and X-points
;   primary_opt          - Index of plasma centre O-point
;   inner_sep            - X-point index of inner separatrix
;   opt_ri, opt_zi       - R and Z indices for each O-point
;   opt_f                - Psi value at each O-point
;   xpt_ri, xpt_zi       - R and Z indices for each X-point
;   xpt_f                - Psi value of each X-point
; 
; If critical is not supplied, it is calculated
; 
; Boundary is a 2D array of indices
;   boundary[0, *]       - R values
;   boundary[1, *]       - Z values
;
; fpsi is an (optional) current function as a 2D array
;   fpsi[0,*]            - Psi values
;   fpsi[1,*]            - f values
;
; Return structure contains:
; 
;   error                  - Non-zero if an error occurred
;   psi_inner, psi_outer   - Normalised ranges of psi used
;   nrad, npol             - Number of grid points used in each region
;   Rixy, Zixy             - 2D arrays of indices into R and Z array
;   Rxy, Zxy               - 2D arrays of (rad,pol) grid-point locations
;   psixy                  - 2D array of normalised psi at each point
;   faxis, fnorm           - Psi normalisation factors
;   settings               - Final settings used
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; radial grid
;
; n - number of grid points
; pin, pout - range of psi
; seps - locations of separatrices
; sep_factor - separatrix peaking

FUNCTION radial_grid, n, pin, pout, include_in, include_out, seps, sep_factor
  x = FINDGEN(n)
  m = FLOAT(n-1)
  IF NOT include_in THEN BEGIN
    x = x + 0.5
    m = m + 0.5
  ENDIF
  
  IF NOT include_out THEN m = m + 0.5
  x = x / m

  RETURN, pin + (pout - pin)*x
END

FUNCTION STR, val
  n = N_ELEMENTS(val)
  IF n GT 1 THEN BEGIN
    s = "(" + STRTRIM(STRING(val[0]),2)
    FOR i=1, n-1 DO s = s + ", " + STRTRIM(STRING(val[i]),2)
    s = s + ")"
    RETURN, s
  ENDIF ELSE RETURN, STRTRIM(STRING(val),2)
END

PRO swap, a, b
  v = a
  a = b
  b = v
END

FUNCTION diff, arr
  r = arr[1:*] - arr[0:(N_ELEMENTS(arr)-2)]
  RETURN, [r[0], r]
END

FUNCTION remove_ind, arr, ind
  n = N_ELEMENTS(arr)
  IF ind EQ 0 THEN BEGIN
    RETURN, arr[1:*]
  ENDIF ELSE IF ind EQ n-1 THEN BEGIN
    RETURN, arr[0:(ind-1)]
  ENDIF
  
  RETURN, [arr[0:(ind-1)], arr[(ind+1):*]]
END

; Return the contour lines of a given set of levels
; Just a small wrapper around the built-in CONTOUR
PRO contour_lines, z, x, y, levels=levels, $
                   path_info=path_info, path_xy=path_xy, $
                   _extra=_extra
  old_x = !x
  old_y = !y
  CONTOUR, z, x, y, levels=levels, $
    path_info=path_info, path_xy=path_xy, $
    closed=0, /PATH_DATA_COORDS, _extra=_extra
  !x = old_x
  !y = old_y
END

; Find the closest contour line to a given point
FUNCTION closest_line, info, xy, ri, zi, mind=mind
  mind = MIN( (xy[0,info[0].offset:(info[0].offset+info[0].n-1)] - ri)^2 + $
              (xy[1,info[0].offset:(info[0].offset+info[0].n-1)] - zi)^2 )
  ind = 0
  FOR i=1, N_ELEMENTS(info)-1 DO BEGIN
    d = MIN( (xy[0,info[i].offset:(info[i].offset+info[i].n-1)] - ri)^2 + $
             (xy[1,info[i].offset:(info[i].offset+info[i].n-1)] - zi)^2 )
    IF d LT mind THEN BEGIN
      mind = d
      ind = i
    ENDIF
  ENDFOR
  RETURN, ind
END

PRO oplot_contour, info, xy, R, Z, periodic=periodic, _extra=_extra
  ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
  zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
  
  IF KEYWORD_SET(periodic) THEN BEGIN
    ri = [ri, ri[0]]
    zi = [zi, zi[0]]
  ENDIF
  OPLOT, INTERPOLATE(R, ri), INTERPOLATE(Z, zi), _extra=_extra
END

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find the gradient at a given point by fitting
; 
; Uses DCT to interpolate and get gradients
; 
; dctF   - the DCT of F
; ri, zi - R and Z indices
;

FUNCTION local_gradient, dctF, ri, zi, status=status
  s = SIZE(dctF, /DIMENSION)
  nr = s[0]
  nz = s[1]
  
  IF (ri LT 0) OR (ri GT nr-1) OR (zi LT 0) OR (zi GT nz-1) THEN BEGIN
    status = 1
    RETURN, 0
  ENDIF
  
  res = EvalCosPfast(dctF, x0=ri, y0=zi)
  
  status = 0

  RETURN, {f:res[0], dfdr:res[1], dfdz:res[2]}
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; Poloidal grid
;
; Divide up a poloidal arc
;

FUNCTION poloidal_grid, dctF, R, Z, ri, zi, n, fpsi=fpsi, parweight=parweight

  IF NOT KEYWORD_SET(parweight) THEN parweight = 0.0
  
  np = N_ELEMENTS(ri)

  ; Calculate poloidal distance along starting line
  drdi = DERIV(INTERPOLATE(R, ri))
  dzdi = DERIV(INTERPOLATE(Z, zi))


  dldi = SQRT(drdi^2 + dzdi^2)
  poldist = int_func(findgen(np), dldi) ; Poloidal distance along line

  IF SIZE(fpsi, /dim) EQ 2 THEN BEGIN
    ; Parallel distance along line
    ; Need poloidal and toroidal field
    ni = N_ELEMENTS(ri)
    bp = FLTARR(ni)
    bt = FLTARR(ni)
    FOR i=0, ni-1 DO BEGIN
      g = local_gradient(dctF, ri[i], zi[i], status=status)
      bp[i] = SQRT(g.dfdr^2 + g.dfdz^2) / INTERPOLATE(R, ri[i])
      bt[i] = INTEPOL(REFORM(fpsi[1,*]), REFORM(fpsi[0,*]), g.f)
    ENDFOR
  ENDIF ELSE pardist = poldist ; Just use the same poloidal distance

  dist = parweight*pardist + (1. - parweight)*poldist

  ; Divide up distance. No points at the end (could be x-point)
  dloc = dist[np-1] * (FINDGEN(n)+0.5)/FLOAT(n)  ; Distance locations
  
  ; Get indices in ri, zi
  ind = INTERPOL(FINDGEN(np), dist, dloc)
  
  RETURN, ind
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Follow the gradient from a given point to a target f
; 

; Calculate dR/df and dZ/df for use by LSODE
; Input: pos[0] = R, pos[1] = Z
; Output [0] = dR/df = -Bz/B^2 , [1] = dZ/df = Br/B^2
FUNCTION radial_differential, fcur, pos
  COMMON rd_com, dctf, lastgoodf, lastgoodpos
  a = local_gradient(dctf, pos[0], pos[1], status=status)
  
  IF status EQ 0 THEN BEGIN
    ; No error in localgradient.
    lastgoodf = fcur 
    lastgoodpos = pos
    
    ; If status NE 0 then an error occurred.
    ; Allow a.dfdz to cause an error so escape LSODE
  ENDIF
  
  ; Check mismatch between fcur and a.f ? 
  Br = a.dfdz
  Bz = -a.dfdr
  B2 = Br^2 + Bz^2

  RETURN, [-Bz/B2, Br/B2]
END

;
; F         (in)    2D psi data
; ri0,zi0   (in)    Starting indices
; ftarget   (in)    The f to aim for
; ri,zi     (out)   Final position
; status    (out)   Non-zero if hits a boundary. 1 if out of range at
;                   the start, 2 if before reaching ftarget
; boundary  (in)    Optional 2D array [2,n] of points on the boundary
; fbndry    (out)   If hits boundary, gives final f value
; ibndry    (out)   If hits boundary, index where hit
;
PRO follow_gradient, dctF, ri0, zi0, ftarget, ri, zi, status=status, $
                     boundary=boundary, fbndry=fbndry, ibndry=ibndry
  COMMON rd_com, df, lastgoodf, lastgoodpos
  
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
    IF ncross GT 0 THEN BEGIN
      REPEAT BEGIN
        ibndry = inds2[0]
        g = local_gradient(dctF, cpos[0,0], cpos[1,0], status=status)
        fbndry = (g.f - f0)*0.99 + f0 ; Just inside boundary
        ; Re-follow to get location
        rzold = [ri0, zi0]
        rznew = LSODE(rzold,f0,fbndry - f0,'radial_differential')
        ri = rznew[0]
        zi = rznew[1]
        cpos = line_crossings([ri0, ri], [zi0, zi], 0, $
                              boundary[0,*], boundary[1,*], 1, ncross=ncross)
      ENDREP UNTIL ncross EQ 0
      ;PRINT, "Hit boundary", ri, zi, " f =", fbndry
      status = 2
      RETURN
    ENDIF
  ENDIF

  status = 0
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Create a grid around a given line
;
; Arguments:
;   dctF         DCT of 2D (R,Z) psi function to be gridded
;   R, Z         
;   ri, zi       1D indices into F for starting line
;   f0           Starting f for the line
;   fin, fout    Range of f to grid 
;   npar         Number of perpendicular lines to generate
;   nin, nout    Number of points inside and outside line
FUNCTION grid_region, dctF, R, Z, $
                      ri, zi, $       ; Starting line to grid.
                      fvals, $        ; Location of the surfaces
                      sind, $         ; Index in fvals of the starting line
                      npar, $         ; Number of points along the line
                      slast=slast, $  ; Index in fvals of last successful point
                      sfirst=sfirst, $
                      oplot=oplot, $
                      boundary=boundary, $
                      ffirst=ffirst, flast=flast, $
                      fpsi=fpsi ; f(psi) = R*Bt optional current function
  
  nsurf = N_ELEMENTS(fvals)
  
  IF sind GE 0 THEN BEGIN
    ; starting position is on one of the output surfaces
    f0 = fvals[sind]
    nin = sind
  ENDIF ELSE BEGIN
    ; Starting position between surfaces
    n = FIX(N_ELEMENTS(ri)/2)
    g = local_gradient(dctF, ri[n], zi[n], status=status)
    f0 = g.f
    
    IF fvals[0] LT fvals[nsurf-1] THEN BEGIN
      w = WHERE(fvals LE f0, nin)
    ENDIF ELSE BEGIN
      w = WHERE(fvals GE f0, nin)
    ENDELSE
  ENDELSE
  nout = nsurf - nin - 1

  sfirst = 0      ; Innermost successful index
  slast = nsurf-1 ; Last successful index 

  ffirst = fvals[sfirst]
  flast = fvals[slast]

  PRINT, "    => Gridding range: ", MIN(fvals), MAX(fvals)

  s = SIZE(dctF, /dimensions)
  nr = s[0]
  nz = s[1]

  ind = poloidal_grid(dctF, R, Z, ri, zi, npar, fpsi=fpsi)
  
  rii = INTERPOLATE(ri, ind)
  zii = INTERPOLATE(zi, ind)

  ; Refine the location of the starting point
  FOR i=0, npar-1 DO BEGIN
    follow_gradient, dctF, rii[i], zii[i], f0, ri1, zi1
    rii[i] = ri1
    zii[i] = zi1
  ENDFOR

  ; From each starting point, follow gradient in both directions
  
  rixy = FLTARR(nsurf, npar)
  zixy = FLTARR(nsurf, npar)
  FOR i=0, npar-1 DO BEGIN
    IF sind GE 0 THEN BEGIN
      rixy[nin, i] = rii[i]
      zixy[nin, i] = zii[i]
    ENDIF ELSE BEGIN
      ; fvals[nin] should be just outside the starting position
      ftarg = fvals[nin]
      follow_gradient, dctF, rii[i], zii[i], $
        ftarg, rinext, zinext, status=status
      rixy[nin, i] = rinext
      zixy[nin, i] = zinext
    ENDELSE
    FOR j=0, nout-1 DO BEGIN
      ftarg = fvals[nin+j+1]
      
      follow_gradient, dctF, rixy[nin+j, i], zixy[nin+j, i], $
        ftarg, rinext, zinext, status=status, $
        boundary=boundary, fbndry=fbndry

      IF status EQ 1 THEN BEGIN
        rixy[nin+j+1, i] = -1.0
        IF nin+j LT slast THEN slast = nin+j ; last good surface index
        fbndry = fvals[slast]
        IF (fvals[1] - fvals[0])*(flast - fbndry) GT 0 THEN flast = fbndry
        BREAK
      ENDIF ELSE IF status EQ 2 THEN BEGIN
        ; Hit a boundary 
        rixy[nin+j+1, i] = rinext
        zixy[nin+j+1, i] = zinext
        IF nin+j LT slast THEN slast = nin+j ; Set the last point
        IF (fvals[1] - fvals[0])*(flast - fbndry) GT 0 THEN flast = fbndry
        BREAK
      ENDIF ELSE BEGIN
        rixy[nin+j+1, i] = rinext
        zixy[nin+j+1, i] = zinext
      ENDELSE
    ENDFOR
    FOR j=0, nin-1 DO BEGIN
      ftarg = fvals[nin-j-1]
      
      follow_gradient, dctF, rixy[nin-j, i], zixy[nin-j, i], $
        ftarg, rinext, zinext, status=status, $
        boundary=boundary, fbndry=fbndry
      
      IF status EQ 1 THEN BEGIN
        rixy[nin-j-1, i] = -1.0
        IF nin-j GT sfirst THEN sfirst = nin-j
        fbndry = fvals[sfirst]
        IF (fvals[1] - fvals[0])*(ffirst - fbndry) LT 0 THEN ffirst = fbndry
        BREAK
      ENDIF

      rixy[nin-j-1, i] = rinext
      zixy[nin-j-1, i] = zinext

      IF status EQ 2 THEN BEGIN
        IF nin-j GT sfirst THEN sfirst = nin-j
        IF (fvals[1] - fvals[0])*(ffirst - fbndry) LT 0 THEN ffirst = fbndry
        BREAK
      ENDIF
    ENDFOR
  ENDFOR

  RETURN, {rixy:rixy, zixy:zixy, rxy:INTERPOLATE(R, rixy), zxy:INTERPOLATE(Z, zixy)}
END

PRO plot_grid_section, a, _extra=_extra
  s = SIZE(a.rxy, /dimension)
  nrad = s[0]
  npol = s[1]
  
  FOR j=0, nrad-1 DO BEGIN
    OPLOT, a.Rxy[j,*], a.Zxy[j,*], _extra=_extra
  ENDFOR
  
  FOR j=0, npol-1 DO BEGIN
    OPLOT, a.Rxy[*,j], a.Zxy[*,j], _extra=_extra
  ENDFOR
END

; Plot a line from a given starting point to a given f
PRO oplot_line, dctF, R, Z, ri0, zi0, fto, npt=npt, color=color, _extra=_extra
  IF NOT KEYWORD_SET(npt) THEN npt=10
  IF NOT KEYWORD_SET(color) THEN color=1
  
  ; Get starting f
  
  g = local_gradient(dctF, ri0, zi0)
  ffrom = g.f

  rixpt = FLTARR(npt+1)
  zixpt = rixpt
  rixpt[0] = ri0
  zixpt[0] = zi0
  FOR j=0, npt-1 DO BEGIN
    d = FLOAT(j+1)/FLOAT(npt)
    ftarg = d*fto + (1.0 - d)*ffrom
    follow_gradient, dctF, rixpt[j], zixpt[j], $
      ftarg, rinext, zinext
    rixpt[j+1] = rinext
    zixpt[j+1] = zinext
  ENDFOR
      
  OPLOT, INTERPOLATE(R, rixpt), INTERPOLATE(Z, zixpt), $
    color=color, _extra=_extra
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Locate an x-point 
; 
; F        F(R,Z) psi array
; ri, zi   1D array of indices into F describing a path
; xpt_ri, xpt_zi  Location index of x-point
; xpt_f    F at the x-point
;

FUNCTION loc_xpt_newt, mini
  COMMON loc_xpt_newt, dctF, fri, fzi, xpt_ri, xpt_zi, xpt_f
  
  ; Follow gradient towards the separatrix
  follow_gradient, dctF, fft_interp(fri, mini), fft_interp(fzi, mini), $
      xpt_f, eri, ezi
  
  ; Work out the distance to the x-point, including
  ; a sign to determine which direction to go in
  ; Cross product of the unit gradient and vector to the x-point
 
  ;PRINT, mini, fft_interp(fri, mini), fft_interp(fzi, mini)
  
  a = local_gradient(dctF, eri, ezi)

  gmag = SQRT(a.dfdr^2 + a.dfdz^2)
  gr = a.dfdr / gmag
  gz = a.dfdz / gmag
 
  ;PRINT, "  -> ", gmag, gr, gz, gr*(ezi - xpt_zi) - gz*(eri - xpt_ri)

  RETURN, gr*(ezi - xpt_zi) - gz*(eri - xpt_ri)
END

FUNCTION locate_xpoint, dctF, ri, zi, xpt_ri, xpt_zi, xpt_f, pos=pos
  COMMON loc_xpt_newt, fdata, fri, fzi, xr, xz, xf
  
  fdata = dctF
  fri = FFT(ri) ; for interpolating
  fzi = FFT(zi)
  xr = xpt_ri
  xz = xpt_zi
  xf = xpt_f
  
  ; Get starting index
  mind = MIN((ri - xpt_ri)^2 + (zi - xpt_zi)^2, mini)
  
  PRINT, "Starting index: ", mini

  CATCH, theError
  IF theError EQ 0 THEN BEGIN
    ; Use Newton method to find x-point
    mini = (NEWTON(mini, 'loc_xpt_newt'))[0]
    CATCH, /cancel
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; Just keep original mini
  ENDELSE
  

  ; interpolate to this position
  pos = [fft_interp(fri,mini), fft_interp(fzi,mini)]
  
  RETURN, mini MOD N_ELEMENTS(ri)
END


FUNCTION theta_differential, s, pos
;
; Expression of differentials
; dr/ds = Br/B
; dz/ds = Bz/B
;
; Inputs: x->s (distance along theta line)
; y[0:1] -> [r,z] (coords)
;------------------------------------------------;
  COMMON td_com, dctf, lastgoodpos
  a = local_gradient(dctf, pos[0], pos[1], status=status)
  
  IF status EQ 0 THEN BEGIN
    ; No error in localgradient.
    lastgoodpos = pos
    
    ; If status NE 0 then an error occurred.
    ; Allow a.dfdz to cause an error so escape LSODE
  ENDIF
  
  Br = a.dfdz
  Bz = -a.dfdr
  B = SQRT(Br^2 + Bz^2)

  RETURN, [Br/B, Bz/B]
END

FUNCTION theta_line, dctF, ri0, zi0, di, nstep, boundary=boundary
  COMMON td_com, fdata, lastgoodpos
  
  ri = [ri0]
  zi = [zi0]
  fdata = dctF
  
  pos = [ri0,zi0]
  FOR i=1, nstep DO BEGIN
    CATCH, theError
    IF theError EQ 0 THEN BEGIN
      pos = LSODE(pos, 0, di, 'theta_differential')
      ri = [ri, pos[0]]
      zi = [zi, pos[1]]
    ENDIF ELSE BEGIN
      ; Error occurred. Should have been caused by hitting a boundary
      ; and/or going off the end of the domain
      CATCH, /cancel
      ri = [ri, lastgoodpos[0]]
      zi = [zi, lastgoodpos[1]]
      RETURN, [[ri], [zi]]
    ENDELSE

    IF KEYWORD_SET(boundary) THEN BEGIN
      ; Check if crossed a boundary
      n = N_ELEMENTS(ri)
      cpos = line_crossings(ri[(n-2):*], zi[(n-2):*], 0, $
                            boundary[0,*], boundary[1,*], 1, ncross=ncross)
      IF ncross GT 0 THEN BEGIN
        ; Change the last point to the intersection location
        
        ri[n-1] = cpos[0,0]
        zi[n-1] = cpos[1,0]
        
        RETURN, [[ri], [zi]]
      ENDIF
    ENDIF
  ENDFOR
  RETURN, [[ri], [zi]]
END

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
  fxx = d[3]
  fyy = d[4]
  fxy = d[5]
  
  IF ABS(fyy) GT 1e-6 THEN BEGIN
    ; Get gradients 1 and 2 (solutions y = g1 * x and y = g2 * x)
    
    g1 = ( -fxy + SQRT(fxy^2 - fxx*fyy) ) / fyy
    g2 = ( -fxy - SQRT(fxy^2 - fxx*fyy) ) / fyy
    
    ; Components of vectors

    v1 = [1, g1]
    v2 = [1, g2]
  ENDIF ELSE BEGIN
    ; One of the lines through the x-point is vertical (x = const)    
    v1 = [0, 1]
    v2 = [1, -fxx / (2.*fxy)]
  ENDELSE

  ; For each line, work out which direction to go away from the
  ; primary O-point
  
  ; Simple (but possibly error-prone): Dot product with line from
  ; X-point to primary O-point

  v0 = [opt_ri - xpt_ri, opt_zi - xpt_zi]

  ; Dot-product
  d1 = TOTAL(v0 * v1)
  d2 = TOTAL(v0 * v2)
  
  ; For private flux region need to head away from primary O-point
  IF d1 GT 0.0 THEN v1 = -v1
  IF d2 GT 0.0 THEN v2 = -v2
  
  ; Normalise to get unit vectors
  v1 = v1 / SQRT(TOTAL(v1^2))
  v2 = v2 / SQRT(TOTAL(v2^2))
 
  di = 0.1
  
  ; Need to decide which direction in theta this is
  
  dt = theta_differential(0., [xpt_ri + di*v1[0], xpt_zi + di*v1[1]])
  sign = 1.
  IF TOTAL(dt * v1) LT 0 THEN sign = -1.

  line1 = theta_line( dctF, $
                      xpt_ri + di*v1[0], xpt_zi + di*v1[1], $
                      sign*di, 1000, boundary=boundary)
  
  dt = theta_differential(0., [xpt_ri + di*v2[0], xpt_zi + di*v2[1]])
  sign = 1.
  IF TOTAL(dt * v2) LT 0 THEN sign = -1.

  line2 = theta_line( dctF, $
                      xpt_ri + di*v2[0], xpt_zi + di*v2[1], $
                      sign*di, 1000, boundary=boundary)

  OPLOT, INTERPOLATE(R, line1[*,0]), INTERPOLATE(Z, line1[*,1]), color=3, thick=2
  OPLOT, INTERPOLATE(R, line2[*,0]), INTERPOLATE(Z, line2[*,1]), color=4, thick=2

  RETURN, {leg1:line1, leg2:line2}
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Equilibrium analysis routine
; 
; Takes a RZ psi grid, and finds x-points and o-points
; 
; F - F(nr, nz) 2D array of psi values
; R - R(nr) 1D array of major radii
; Z - Z(nz) 1D array of heights
;
; Returns a structure of critical points containing:
;
;   n_opoint, n_xpoint   - Number of O- and X-points
;   primary_opt          - Index of plasma centre O-point
;   inner_sep            - X-point index of inner separatrix
;   opt_ri, opt_zi       - R and Z indices for each O-point
;   opt_f                - Psi value at each O-point
;   xpt_ri, xpt_zi       - R and Z indices for each X-point
;   xpt_f                - Psi value of each X-point
; 
FUNCTION analyse_equil, F, R, Z
  s = SIZE(F, /DIMENSION)
  nx = s[0]
  ny = s[1]
  
  ;;;;;;;;;;;;;;;; Find critical points ;;;;;;;;;;;;;
  ;
  ; Need to find starting locations for O-points (minima/maxima)
  ; and X-points (saddle points)
  ;
  
  dfdr = FLTARR(nx, ny)
  FOR j=0, ny-1 DO BEGIN
    dfdr[*,j] = diff(f[*,j])
  ENDFOR
  
  dfdz = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
    dfdz[i,*] = diff(f[i,*])
  ENDFOR

  ; Use contour to get crossing-points where dfdr = dfdz = 0
  
  contour_lines, dfdr, findgen(nx), findgen(ny), levels=[0.0], $
    path_info=rinfo, path_xy=rxy
  
  contour_lines, dfdz, findgen(nx), findgen(ny), levels=[0.0], $
    path_info=zinfo, path_xy=zxy
  
  ; Find where these two cross
  nextrema = 0
  
  FOR i=0, N_ELEMENTS(rinfo)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(zinfo)-1 DO BEGIN
      cross = line_crossings(rxy[0,rinfo[i].offset:(rinfo[i].offset + rinfo[i].n - 1)], $
                             rxy[1,rinfo[i].offset:(rinfo[i].offset + rinfo[i].n - 1)], $
                             rinfo[i].type, $
                             zxy[0,zinfo[j].offset:(zinfo[j].offset + zinfo[j].n - 1)], $
                             zxy[1,zinfo[j].offset:(zinfo[j].offset + zinfo[j].n - 1)], $
                             zinfo[j].type, $
                             ncross=ncross)
      IF ncross GT 0 THEN BEGIN
        IF nextrema EQ 0 THEN BEGIN
          rex = REFORM(cross[0,*])
          zex = REFORM(cross[1,*])
        ENDIF ELSE BEGIN
          rex = [rex, REFORM(cross[0,*])]
          zex = [zex, REFORM(cross[1,*])]
        ENDELSE
        nextrema = nextrema + ncross
      ENDIF
    ENDFOR
  ENDFOR
  
  rex = ROUND(rex)
  zex = ROUND(zex)
  
  PRINT, ""
  PRINT, "Number of critical points:", nextrema
  
  ;;;;;;;;;;;;;;; Characterise extrema ;;;;;;;;;;;;;;;;;
  ; Fit a surface through local points using 6x6 matrix
  ; This is to determine the type of extrema, and to
  ; refine the location
  ; 
  
  n_opoint = 0
  n_xpoint = 0

  ; Calculate inverse matrix
  rio = [-1, 0, 0, 0, 1, 1] ; R index offsets
  zio = [ 0,-1, 0, 1, 0, 1] ; Z index offsets
  
  ; Fitting a + br + cz + drz + er^2 + fz^2
  A = TRANSPOSE([[FLTARR(6)+1], $
                 [rio], $
                 [zio], $
                 [rio*zio], $
                 [rio^2], $
                 [zio^2]])
  SVDC, A,W,U,V

  FOR e=0, nextrema-1 DO BEGIN
    ; Fit in index space so result is index number
    
    PRINT, "Critical point "+STRTRIM(STRING(e),2)
    
    valid = 1

    localf = FLTARR(6)
    FOR i=0, 5 DO localf[i] = F[rex[e]+rio[i], zex[e]+zio[i]]
    res=SVSOL(U,W,V,localf)
    
    ; Res now contains [a,b,c,d,e,f]
    ;                  [0,1,2,3,4,5]

    ; This determines whether saddle or extremum
    det = 4.*res[4]*res[5] - res[3]^2
    
    IF det LT 0.0 THEN BEGIN
      PRINT, "   X-point"
    ENDIF ELSE BEGIN
      PRINT, "   O-point"
    ENDELSE
    
    ; Get location (2x2 matrix of coefficients)
    
    rinew = (res[3]*res[2] - 2.*res[1]*res[5]) / det
    zinew = (res[3]*res[1] - 2.*res[4]*res[2]) / det

    IF (ABS(rinew) GT 1.) OR (ABS(zinew) GT 1.0) THEN BEGIN
      ; Method has gone slightly wrong. Try a different method.
      ; Get a contour line starting at this point. Should
      ; produce a circle around the real o-point. 
      PRINT, "   Fitted location deviates too much"
      IF det LT 0.0 THEN BEGIN
        PRINT, "   => X-point probably not valid"
        ;valid = 0
      ENDIF ELSE BEGIN
        
        contour_lines, F, findgen(nx), findgen(ny), levels=[F[rex[e], zex[e]]], $
          path_info=info, path_xy=xy
        
        IF N_ELEMENTS(info) GT 1 THEN BEGIN
          ; More than one contour. Select the one closest
          ind = closest_line(info, xy, rex[e], zex[e])
          info = info[ind]
        ENDIF ELSE info = info[0]
        
        rinew = 0.5*(MAX(xy[0, info.offset:(info.offset + info.n - 1)]) + $
                     MIN(xy[0, info.offset:(info.offset + info.n - 1)])) - rex[e]
        zinew = 0.5*(MAX(xy[1, info.offset:(info.offset + info.n - 1)]) + $
                     MIN(xy[1, info.offset:(info.offset + info.n - 1)])) - zex[e]
        
        IF (ABS(rinew) GT 2.) OR (ABS(zinew) GT 2.0) THEN BEGIN
          PRINT, "   Backup method also failed. Keeping initial guess"
          rinew = 0.
          zinew = 0.
        ENDIF
      ENDELSE
    ENDIF

    IF valid THEN BEGIN
      fnew = res[0] + res[1]*rinew + res[2]*zinew $
        + res[3]*rinew*zinew + res[4]*rinew^2 + res[5]*zinew^2
      
      rinew = rinew + rex[e]
      zinew = zinew + zex[e]
      
      PRINT, "   Starting index: " + STR(rex[e])+", "+STR(zex[e])
      PRINT, "   Refined  index: " + STR(rinew)+", "+STR(zinew)
      
      rnew = INTERPOLATE(R, rinew)
      znew = INTERPOLATE(Z, zinew)
      
      PRINT, "   Position: " + STR(rnew)+", "+STR(znew)
      PRINT, "   F = "+STR(fnew)
      
      IF det LT 0.0 THEN BEGIN
        
        IF n_xpoint EQ 0 THEN BEGIN
          xpt_ri = [rinew]
          xpt_zi = [zinew]
          xpt_f = [fnew]
          n_xpoint = n_xpoint + 1
        ENDIF ELSE BEGIN
          ; Check if this duplicates an existing point
          
          m = MIN((xpt_ri - rinew)^2 + (xpt_zi - zinew)^2, ind)
          IF m LT 2. THEN BEGIN
            PRINT, "   Duplicates existing X-point."
          ENDIF ELSE BEGIN
            xpt_ri = [xpt_ri, rinew]
            xpt_zi = [xpt_zi, zinew]
            xpt_f = [xpt_f, fnew]
            n_xpoint = n_xpoint + 1
          ENDELSE
        ENDELSE
        
      ENDIF ELSE BEGIN
        
        IF n_opoint EQ 0 THEN BEGIN
          opt_ri = [rinew]
          opt_zi = [zinew]
          opt_f = [fnew]
          n_opoint = n_opoint + 1
        ENDIF ELSE BEGIN
          ; Check if this duplicates an existing point
          
          m = MIN((opt_ri - rinew)^2 + (opt_zi - zinew)^2, ind)
          IF m LT 2. THEN BEGIN
            PRINT, "   Duplicates existing O-point"
          ENDIF ELSE BEGIN
            opt_ri = [opt_ri, rinew]
            opt_zi = [opt_zi, zinew]
            opt_f = [opt_f, fnew]
            n_opoint = n_opoint + 1
          ENDELSE
        ENDELSE
      ENDELSE
    ENDIF
  ENDFOR

  PRINT, "Number of O-points: "+STR(n_opoint)
  PRINT, "Number of X-points: "+STR(n_xpoint)
  
  ;;;;;;;;;;;;;;; Find plasma centre ;;;;;;;;;;;;;;;;;;;
  ; Find the O-point closest to the middle of the grid
  
  mind = (opt_ri[0] - (FLOAT(nx)/2.))^2 + (opt_zi[0] - (FLOAT(ny)/2.))
  ind = 0
  FOR i=1, n_opoint-1 DO BEGIN
    d = (opt_ri[i] - (FLOAT(nx)/2.))^2 + (opt_zi[i] - (FLOAT(ny)/2.))
    IF d LT mind THEN BEGIN
      ind = i
      mind = d
    ENDIF
  ENDFOR
  
  primary_opt = ind
  PRINT, "Primary O-point is at "+STR(INTERPOLATE(R, opt_ri[ind])) + $
    ", " + STR(INTERPOLATE(Z, opt_zi[ind]))
  PRINT, ""
  
  IF n_xpoint GT 0 THEN BEGIN
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Find separatrix closest to O-point

    mind = FLTARR(n_xpoint)
    FOR i=0, n_xpoint-1 DO BEGIN
      contour_lines, F, findgen(nx), findgen(ny), levels=[xpt_f[i]], $
        path_info=info, path_xy=xy
      ii = closest_line(info, xy, xpt_ri[i], xpt_zi[i])
      info=info[ii]
      xy = xy[*,info.offset:(info.offset+info.n-1)]
      mind[i] = MIN( (xy[0,*] - opt_ri[ind])^2 + (xy[1,*] - opt_zi[ind])^2 )
    ENDFOR

    si = SORT(mind) ; Sorted indices by distance
    xpt_f = xpt_f[si]
    xpt_ri = xpt_ri[si]
    xpt_zi = xpt_zi[si]
    
    ; Now make sure that psi is monotonic
    FOR i=1, n_xpoint-1 DO BEGIN
      IF ABS(xpt_f[i] - opt_f[ind]) LT ABS(xpt_f[i-1] - opt_f[ind]) THEN BEGIN
        ; Started to reverse. Remove all x-points after this
        xpt_f = xpt_f[0:(i-1)]
        xpt_ri = xpt_ri[0:(i-1)]
        xpt_zi = xpt_zi[0:(i-1)]
        BREAK
      ENDIF
    ENDFOR
    n_xpoint = N_ELEMENTS(xpt_f)

    inner_sep = 0
  ENDIF ELSE BEGIN
    ; No x-points. Pick mid-point in f
   
    xpt_f = 0.5*(MAX(F) + MIN(F))
    
    PRINT, "WARNING: No X-points. Setting separatrix to F = "+STR(xpt_f)

    xpt_ri = 0
    xpt_zi = 0
    inner_sep = 0
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Put results into a structure
  
  result = {n_opoint:n_opoint, n_xpoint:n_xpoint, $ ; Number of O- and X-points
            primary_opt:primary_opt, $ ; Which O-point is the plasma centre
            inner_sep:inner_sep, $ ; Innermost X-point separatrix
            opt_ri:opt_ri, opt_zi:opt_zi, opt_f:opt_f, $ ; O-point location (indices) and psi values
            xpt_ri:xpt_ri, xpt_zi:xpt_zi, xpt_f:xpt_f}  ; X-point locations and psi values
  
  RETURN, result
END

; Overplots the output of analyse_equil
PRO oplot_critical, F, R, Z, a
  
  ; Plot X-points and separatrices
  FOR i=0, a.n_xpoint-1 DO BEGIN
    ; plot the separatrix contour
    CONTOUR, F, R, Z, levels=[a.xpt_f[i]], c_colors=2, /overplot
    oplot, [INTERPOLATE(R, a.xpt_ri[i])], [INTERPOLATE(Z, a.xpt_zi[i])], psym=7, color=2
  ENDFOR

  ; Plot O-points
  FOR i=0, a.n_opoint-1 DO BEGIN
    oplot, [INTERPOLATE(R, a.opt_ri[i])], [INTERPOLATE(Z, a.opt_zi[i])], psym=7, color=3
  ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Select only those X- and O-points within the given range of
; normalised psi
FUNCTION restrict_psi_range, critical, psi_outer
  
  IF critical.n_xpoint EQ 0 THEN RETURN, critical

  xpt_psi = (critical.xpt_f - critical.opt_f[critical.primary_opt]) $
    / (critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt])
  
  PRINT, "X-point locations: "
  FOR i=0, critical.n_xpoint-1 DO BEGIN
    PRINT, "  "+STR(i)+": "+STR(xpt_psi[i])
  ENDFOR

  ; Select the x-points within this range
  w = WHERE(xpt_psi LT psi_outer, count)
  PRINT, "Number of x-points in range: "+STR(count)
  PRINT, ""
  
  IF count EQ 0 THEN BEGIN
    ; Keep the inner sep
    w = [critical.inner_sep]
    inner_sep = 0
  ENDIF ELSE BEGIN
    ; Need to update inner_sep
    w2 = WHERE(w EQ critical.inner_sep)
    inner_sep = w2[0]
  ENDELSE
  
  result = {n_opoint:critical.n_opoint, n_xpoint:count, $
            primary_opt:critical.primary_opt, $
            inner_sep:inner_sep, $
            opt_ri:critical.opt_ri, opt_zi:critical.opt_zi, opt_f:critical.opt_f, $
            xpt_ri:critical.xpt_ri[w], xpt_zi:critical.xpt_zi[w], xpt_f:critical.xpt_f[w]}
  RETURN, result
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make sure all critical points are inside the boundary
;
; bndryi[2,n]   r and z indices of the boundary
FUNCTION critical_bndry, critical, bndryi
  ; First the o-points
  
  primary_opt = critical.primary_opt
  opt_ri = critical.opt_ri
  opt_zi = critical.opt_zi
  opt_f  = critical.opt_f
  
  IF critical.n_opoint GT 1 THEN BEGIN
    w = [critical.primary_opt]
    primary_opt = 0
    
    FOR i=0,critical.n_opoint-1 DO BEGIN
      IF i NE critical.primary_opt THEN BEGIN
        ; Test if outside boundary by drawing a line to the primary o-point
        cp = line_crossings([opt_ri[i], opt_ri[primary_opt]], $
                            [opt_zi[i], opt_zi[primary_opt]], $
                            0, $
                            bndryi[0,*], bndryi[1,*], 1, $
                            ncross=ncross)
        IF ncross EQ 0 THEN w = [w,i]
      ENDIF
    ENDFOR
    
    n_opoint = N_ELEMENTS(w)
    opt_ri = opt_ri[w]
    opt_zi = opt_zi[w]
    opt_f  = opt_f[w]
  ENDIF

  ; Check the x-points
  n_xpoint = 0
  FOR i=0, critical.n_xpoint-1 DO BEGIN
    ; Test if outside boundary by drawing a line to the primary o-point
    cp = line_crossings([critical.xpt_ri[i], opt_ri[primary_opt]], $
                        [critical.xpt_zi[i], opt_zi[primary_opt]], $
                        0, $
                        bndryi[0,*], bndryi[1,*], 1, $
                        ncross=ncross)
    IF ncross EQ 0 THEN BEGIN
      ; hasn't crossed -> inside boundary. Add index to w
      IF n_xpoint EQ 0 THEN w = [i] ELSE w = [w,i]
      n_xpoint = n_xpoint + 1
    ENDIF
  ENDFOR
  
  IF n_xpoint EQ 0 THEN BEGIN
    ; Keep the inner sep (used for normalisation)
    w = [critical.inner_sep]
    inner_sep = 0
  ENDIF ELSE BEGIN
    ; Need to update inner_sep
    w2 = WHERE(w EQ critical.inner_sep)
    inner_sep = w2[0]
  ENDELSE

  ; Select x-points
  xpt_ri = critical.xpt_ri[w]
  xpt_zi = critical.xpt_zi[w]
  xpt_f  = critical.xpt_f[w]

  result = {n_opoint:n_opoint, n_xpoint:n_xpoint, $
            primary_opt:primary_opt, $
            inner_sep:inner_sep, $
            opt_ri:opt_ri, opt_zi:opt_zi, opt_f:opt_f, $
            xpt_ri:xpt_ri, xpt_zi:xpt_zi, xpt_f:xpt_f}
  RETURN, result
END

; Check if a given field is in a structure, and if not then add it
; with a default value
PRO str_check_present, str, name, default
  tag = TAG_NAMES(str)
  IF MAX(STRCMP(tag, name, /FOLD_CASE)) EQ 0 THEN BEGIN
    PRINT, "  "+name+" not specified. Setting to "+STR(default)
    str=CREATE_STRUCT(str, name, default)
  ENDIF
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; X-point refinement using DCT

FUNCTION ref_xpt_newt, X
  ; Return gradients
  COMMON rfn_xpt_com, dctf
  
  a = local_gradient(dctf, X[0], X[1], status=status)

  RETURN, [a.dfdr, a.dfdz]
END

PRO refine_xpoint, dctF, ri0, zi0, ri, zi
  ; Find where Grad(F) = 0
  
  COMMON rfn_xpt_com, fdata
  fdata = dctf
  
  CATCH, theError
  IF theError NE 0 THEN BEGIN
    CATCH, /cancel
    PRINT, "** Error occurred whilst refining x-point location"
    ri = ri0
    zi = zi0
    RETURN
  ENDIF
  
  rz0 = [ri0, zi0]
  rz = NEWTON(rz0, 'ref_xpt_newt')
  
  ri = rz[0]
  zi = rz[1]
  CATCH, /cancel
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main gridding function
;

FUNCTION create_grid, F, R, Z, in_settings, critical=critical, $
                      boundary=boundary, debug=debug, strictbndry=strictbndry, iter=iter, $
                      fpsi = fpsi ; f(psi) = R*Bt current function

  ; Create error handler
  err=0;CATCH, err
  IF err NE 0 THEN BEGIN
    PRINT, "CREATE_GRID failed"
    PRINT, "   Error message: "+!ERROR_STATE.MSG
    CATCH, /cancel
    RETURN, {error:1}
  ENDIF

  IF NOT KEYWORD_SET(iter) THEN iter = 0
  IF iter GT 3 THEN BEGIN
    PRINT, "ERROR: Too many iterations"
    RETURN, {error:1}
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Check the settings
  ; If a setting is missing, set a default value
  
  IF N_PARAMS() LT 3 THEN BEGIN
    PRINT, "ERROR: Need at least a 2D array of psi values, R and Z arrays"
    RETURN, {error:1}
  ENDIF ELSE IF N_PARAMS() LT 4 THEN BEGIN
    ; Settings omitted. Set defaults
    PRINT, "Settings not given -> using default values"
    settings = {psi_inner:0.9, $
                psi_outer:1.1, $
                nrad:36, $
                npol:64, $
                rad_peaking:0.0, $
                pol_peaking:0.0}
  ENDIF ELSE BEGIN
    PRINT, "Checking settings"
    settings = in_settings ; So the input isn't changed
    str_check_present, settings, 'psi_inner', 0.9
    str_check_present, settings, 'psi_outer', 1.1
    str_check_present, settings, 'nrad', 36
    str_check_present, settings, 'npol', 64
    str_check_present, settings, 'rad_peaking', 0.0
    str_check_present, settings, 'pol_peaking', 0.0
  ENDELSE

  s = SIZE(F, /DIMENSION)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
    PRINT, "ERROR: First argument must be 2D array of psi values"
    RETURN, {error:1}
  ENDIF
  nx = s[0]
  ny = s[1]

  s = SIZE(R, /DIMENSION)
  IF (N_ELEMENTS(s) NE 1) OR (s[0] NE nx) THEN BEGIN
    PRINT, "ERROR: Second argument must be 1D array of major radii"
    RETURN, {error:1}
  ENDIF
  
  s = SIZE(Z, /DIMENSION)
  IF (N_ELEMENTS(s) NE 1) OR (s[0] NE ny) THEN BEGIN
    PRINT, "ERROR: Second argument must be 1D array of heights"
    RETURN, {error:1}
  ENDIF

  IF KEYWORD_SET(boundary) THEN BEGIN
    s = SIZE(boundary, /DIMENSION)
    IF (N_ELEMENTS(s) NE 2) OR (s[0] NE 2) THEN BEGIN
      PRINT, "ERROR: boundary must be a 2D array: [2, n]"
      RETURN, {error:1}
    ENDIF
    
    ; Calculate indices
    bndryi = boundary
    bndryi[0,*] = INTERPOL(FINDGEN(nx), R, REFORM(boundary[0,*]))
    bndryi[1,*] = INTERPOL(FINDGEN(ny), Z, REFORM(boundary[1,*]))
  ENDIF
  
  ;;;;;;;;;;;;;;; Calculate DCT ;;;;;;;;;;;;;;
  
  PRINT, "Calculating DCT..."
  DCT2Dslow, F, dctF
  PRINT, "Finished DCT"

  ;;;;;;;;;;;;;;;; First plot ;;;;;;;;;;;;;;;;

  nlev = 100
  minf = MIN(f)
  maxf = MAX(f)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf

  safe_colors, /first
  CONTOUR, F, R, Z, levels=levels, color=1, /iso, xstyl=1, ysty=1
  
  IF KEYWORD_SET(boundary) THEN BEGIN
    OPLOT, [REFORM(boundary[0,*]), boundary[0,0]], [REFORM(boundary[1,*]), boundary[1,0]], $
           thick=2,color=2
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Analyse the equilibrium to find O- and X-points
  
  IF NOT KEYWORD_SET(critical) THEN critical = analyse_equil(F, R, Z)
  critical = restrict_psi_range(critical, MAX(settings.psi_outer))

  ; Check that the critical points are inside the boundary
  IF KEYWORD_SET(boundary) THEN critical = critical_bndry(critical, boundary)

  n_opoint = critical.n_opoint
  n_xpoint = critical.n_xpoint
  primary_opt = critical.primary_opt
  inner_sep   = critical.inner_sep
  opt_ri = critical.opt_ri
  opt_zi = critical.opt_zi
  opt_f  = critical.opt_f
  xpt_ri = critical.xpt_ri
  xpt_zi = critical.xpt_zi
  xpt_f  = critical.xpt_f

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  PRINT, "Refining x-point locations..."
  FOR i=0, n_xpoint-1 DO BEGIN
    PRINT, "  "+STR(i)+": "+STR([xpt_ri[i], xpt_zi[i]])+" "+STR(xpt_f[i])
    refine_xpoint, dctF, xpt_ri[i], xpt_zi[i], ri, zi
    a = local_gradient(dctF, ri, zi)
    PRINT, "   -> "+STR([ri, zi])+" "+STR(a.f)
    xpt_ri[i] = ri
    xpt_zi[i] = zi
    xpt_f[i] = a.f
  ENDFOR

  ; Overplot the separatrices, O-points
  oplot_critical, F, R, Z, critical

  ; Psi normalisation factors
  faxis = critical.opt_f[critical.primary_opt]
  fnorm = critical.xpt_f[critical.inner_sep] - critical.opt_f[critical.primary_opt]

  ; From normalised psi, get range of f
  f_inner = faxis + MIN(settings.psi_inner)*fnorm
  f_outer = faxis + MAX(settings.psi_outer)*fnorm
  
  ; Check the number of x-points
  IF critical.n_xpoint EQ 0 THEN BEGIN
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Grid entirely in the core

    PRINT, "Generating grid entirely in the core"

    nrad = LONG(TOTAL(settings.nrad)) ; Add up all points
    npol = LONG(TOTAL(settings.npol))
    rad_peaking = settings.rad_peaking[0] ; Just the first region
    pol_peaking = settings.pol_peaking[0]
    
    ; work out where to put the surfaces in the core
    fvals = radial_grid(nrad, f_inner, f_outer, 1, 1, [xpt_f[inner_sep]], rad_peaking)

    ; Create a starting surface
    sind = FIX(nrad / 2)
    start_f = fvals[sind]
    
    contour_lines, F, findgen(nx), findgen(ny), levels=[start_f], $
      path_info=info, path_xy=xy
  
    IF N_ELEMENTS(info) GT 0 THEN BEGIN
      ; Find the surface closest to the o-point
      
      ind = closest_line(info, xy, opt_ri[primary_opt], opt_zi[primary_opt])
      info = info[ind]
    ENDIF ELSE info = info[0]
    
    start_ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
    start_zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])

    oplot_contour, info, xy, R, Z, /periodic, color=2, thick=1.5
    
    a = grid_region(dctF, R, Z, $
                    start_ri, start_zi, $
                    fvals, $
                    sind, $
                    npol, fpsi=fpsi)
    
    OPLOT, [REFORM(a.rxy[0,*]), a.rxy[0,0]], [REFORM(a.zxy[0,*]), a.zxy[0,0]], color=4
      
    FOR i=1, nrad-1 DO BEGIN
      OPLOT, [REFORM(a.rxy[i,*]), a.rxy[i,0]], [REFORM(a.zxy[i,*]), a.zxy[i,0]], color=4
    ENDFOR

    FOR i=0, npol-1 DO BEGIN
      OPLOT, a.rxy[*,i], a.zxy[*,i], color=4
    ENDFOR

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Create result structure

    result = {error:0, $ ; Signal success
              psi_inner:settings.psi_inner, psi_outer:settings.psi_outer, $ ; Unchanged psi range
              nrad:nrad, npol:npol, $
              Rxy:a.Rxy, Zxy:a.Zxy}

    RETURN, result
    
  ENDIF ELSE BEGIN
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Grid contains at least one x-point

    ; Normalised psi value of each separatrix
    xpt_psi = (critical.xpt_f - faxis) / fnorm
    
    si = SORT(xpt_psi) ; Sort separatrices from inside out

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Choose primary x-point.
    ; Determines order of settings arrays
    primary_xpt = si[0]

    PRINT, "Primary X-point is number "+STR(primary_xpt)
    PRINT, "   at R = "+STR(INTERPOLATE(R, critical.xpt_ri[primary_xpt])) $
      +" Z = "+STR(INTERPOLATE(Z, critical.xpt_zi[primary_xpt]))
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; work out where to put the surfaces

    nrad = settings.nrad
    nnrad = N_ELEMENTS(nrad)

    nrad_flexible = 0
    IF nnrad NE (critical.n_xpoint + 1) THEN BEGIN
      ; Only given total number of radial points. Decide
      ; distribution automatically
      
      nrad_flexible = 1 ; Flag used later if nrad not fixed

      IF nnrad GT 1 THEN BEGIN
        PRINT, "WARNING: nrad has wrong number of elements ("+STR(nnrad)+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint + 1)+" elements"
      ENDIF

      PRINT, "Distributing radial points automatically"
      
      fvals = radial_grid(TOTAL(nrad), f_inner, f_outer, 1, 1, xpt_f, rad_peaking)
      psi_vals = (fvals - faxis) / fnorm ; Normalised psi
      
      ; Find out how many grid points were used
      nrad = LONARR(critical.n_xpoint + 1)
      tot = 0
      FOR i=0, critical.n_xpoint-1 DO BEGIN
        w = WHERE(psi_vals LT xpt_psi[si[i]], count)
        nrad[i] = count - tot
        tot = tot + nrad[i]
      ENDFOR
      w = WHERE(psi_vals GT MAX(xpt_psi), count)
      nrad[critical.n_xpoint] = count
      
    ENDIF ELSE BEGIN
      ; Specified number of points in each region
      
      ; Core first. Include first point
      fvals = radial_grid(nrad[0], f_inner, xpt_f[inner_sep], 1, 0, xpt_f, rad_peaking)
      ; Regions between separatrices. Don't include any end points
      FOR i=1, critical.n_xpoint-1 DO BEGIN
        fvals = [fvals, $
                 radial_grid(nrad[i], xpt_f[si[i-1]], xpt_f[si[i]], 0, 0, xpt_f, rad_peaking)]
      ENDFOR
      ; SOL outside all separatrices. End point is included
      fvals = [fvals, $
               radial_grid(nrad[critical.n_xpoint], $
                           xpt_f[si[critical.n_xpoint-1]], f_outer, 0, 1, xpt_f, rad_peaking)]
      psi_vals = (fvals - faxis) / fnorm ; Normalised psi
    ENDELSE
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Create arrays of psi values for each region
    sol_psi_vals = FLTARR(critical.n_xpoint, TOTAL(nrad))
    pf_psi_vals  = FLTARR(critical.n_xpoint, 2, TOTAL(nrad))
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      sol_psi_vals[i,*]  = psi_vals
      pf_psi_vals[i,0,*] = psi_vals
      pf_psi_vals[i,1,*] = psi_vals
    ENDFOR
    
    IF N_ELEMENTS(settings.psi_inner) EQ (critical.n_xpoint + 1) THEN BEGIN
      psi_inner = settings.psi_inner
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(settings.psi_inner) GT 1 THEN BEGIN
        PRINT, "WARNING: psi_inner has wrong number of elements ("+STR(N_ELEMENTS(settings.psi_inner))+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint+1)+" elements"
      ENDIF
      PRINT, "Keeping same inner psi for all regions"
      
      psi_inner = FLTARR(critical.n_xpoint+1) + MIN(settings.psi_inner)
    ENDELSE
    
    IF N_ELEMENTS(settings.psi_outer) EQ critical.n_xpoint THEN BEGIN
      psi_outer = settings.psi_outer
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(settings.psi_outer) GT 1 THEN BEGIN
        PRINT, "WARNING: psi_outer has wrong number of elements ("+STR(N_ELEMENTS(settings.psi_outer))+")"
        PRINT, "         Should have 1 or "+STR(critical.n_xpoint)+" elements"
      ENDIF 
      PRINT, "Keeping same outer psi for all regions"
      
      psi_outer = FLTARR(critical.n_xpoint) + MAX(settings.psi_outer)
    ENDELSE

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Need to work out where the x-points are
    ; and create lines just inside the separatrices
    
    ; take a contour surface just inside the innermost separatrix
    sind = FIX(nrad[0]/2) ;nrad[0]-1 ; the last point inside core
    f_cont = fvals[sind]
    
    contour_lines, F, findgen(nx), findgen(ny), levels=[f_cont], $
      path_info=info, path_xy=xy
    
    IF N_ELEMENTS(info) GT 0 THEN BEGIN
      ; Find the surface closest to the o-point

      info = info[closest_line(info, xy, $
                               critical.opt_ri[critical.primary_opt], critical.opt_zi[critical.primary_opt])]
    ENDIF ELSE info = info[0]
    
    oplot_contour, info, xy, R, Z, /periodic, color=3, thick=1.5

    start_ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
    start_zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])

    ; Make sure that the line goes clockwise
    
    m = MAX(INTERPOLATE(Z, start_zi), ind)
    IF (DERIV(INTERPOLATE(R, start_ri)))[ind] LT 0.0 THEN BEGIN
      ; R should be increasing at the top. Need to reverse
      start_ri = REVERSE(start_ri)
      start_zi = REVERSE(start_zi)
    ENDIF

    ; now have (start_ri, start_zi). For each x-point, find the radial
    ; line going through the x-point
    
    fri = FFT(start_ri) ; for interpolating periodic functions
    fzi = FFT(start_zi)

    xpt_ind = FLTARR(critical.n_xpoint)  ; index into start_*i
    
    pf_info = PTRARR(critical.n_xpoint)
    
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      PRINT, "Finding theta location of x-point "+STR(i)
      
      ; Follow a line from the x-point to this surface
      vdr = opt_ri[primary_opt] - xpt_ri[i]
      vdz = opt_zi[primary_opt] - xpt_zi[i]
      mag = SQRT(vdr^2 + vdz^2)
      vdr = vdr / mag
      vdz = vdz / mag
      ; Find index where it hits the surface (boundary)
      follow_gradient, dctF, $
        xpt_ri[i] + 0.5*vdr, xpt_zi[i] + 0.5*vdz, $
        0.95 * f_cont + 0.05*opt_f[primary_opt], $
        rhit, zhit, $
        boundary=TRANSPOSE([[start_ri], [start_zi]]), ibndry=hit_ind
      
      IF hit_ind GE 0 THEN BEGIN
        mini = hit_ind
      ENDIF ELSE BEGIN
        ; Didn't hit the surface. Use backup method
      
        PRINT, "WARNING: Using backup method"
        
        mini = locate_xpoint(dctF, start_ri, start_zi, $
                             xpt_ri[i], xpt_zi[i], xpt_f[i])
        
      ENDELSE
      
      xpt_ind[i] = mini  ; Record the index
        
      ; Plot the line to the x-point
      oplot_line, dctF, R, Z, $
        fft_interp(fri, mini), fft_interp(fzi, mini), critical.xpt_f[i]
      oplot_line, dctF, R, Z, $
        fft_interp(fri, mini), fft_interp(fzi, mini), f_inner

      ; Get tangent vector
      drdi = INTERPOLATE((DERIV(INTERPOLATE(R, start_ri))), mini)
      dzdi = INTERPOLATE((DERIV(INTERPOLATE(Z, start_zi))), mini)
      
      tmp = {core_ind:mini, drdi:drdi, dzdi:dzdi, $ ; Core index and tangent vector
             sol:LONARR(2)} ; Array to store SOL indices
      
      pf_info[i] = PTR_NEW(tmp)
    ENDFOR

    ; Sort x-points by index along this core surface
    ci = SORT(xpt_ind)
    
    ; Extract starting lines for each section
    sol_info = PTRARR(critical.n_xpoint)
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      IF i NE (critical.n_xpoint-1) THEN BEGIN
        ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
               start_ri[FIX(xpt_ind[ci[i]]+1.0):FIX(xpt_ind[ci[i+1]])], $
               fft_interp(fri,xpt_ind[ci[i+1]]) ]
        
        zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
               start_zi[FIX(xpt_ind[ci[i]]+1.0):FIX(xpt_ind[ci[i+1]])], $
               fft_interp(fzi,xpt_ind[ci[i+1]]) ]
      ENDIF ELSE BEGIN
        ; Index wraps around
        IF xpt_ind[ci[i]] GT N_ELEMENTS(start_ri)-2 THEN BEGIN
          ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
                 start_ri[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fri,xpt_ind[ci[0]]) ]
          
          zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
                 start_zi[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fzi,xpt_ind[ci[0]]) ]
          
        ENDIF ELSE BEGIN
          ri = [ fft_interp(fri,xpt_ind[ci[i]]), $
                 start_ri[FIX(xpt_ind[ci[i]]+1.0):*], $
                 start_ri[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fri,xpt_ind[ci[0]]) ]
          
          zi = [ fft_interp(fzi,xpt_ind[ci[i]]), $
                 start_zi[FIX(xpt_ind[ci[i]]+1.0):*], $
                 start_zi[0:FIX(xpt_ind[ci[0]])], $
                 fft_interp(fzi,xpt_ind[ci[0]]) ]
        ENDELSE
      ENDELSE
      
      ; Calculate length of the line
      drdi = DERIV(INTERPOLATE(R, ri))
      dzdi = DERIV(INTERPOLATE(Z, zi))
      dldi = SQRT(drdi^2 + dzdi^2)
      length = INT_TABULATED(findgen(N_ELEMENTS(dldi)), dldi)

      ; Change the grid in the far SOL
      w = WHERE(ci EQ primary_xpt)
      solid = (i - w[0] + critical.n_xpoint) MOD critical.n_xpoint
      xpt_psi_max = MAX([xpt_psi[ci[i]], xpt_psi[ci[(i+1) MOD critical.n_xpoint]]])
      w = WHERE(sol_psi_vals[i,*] GT xpt_psi_max, nsol)
      PRINT, "Generating sol psi: ", xpt_psi_max, psi_outer[solid]
      IF xpt_psi_max GT psi_outer[solid] THEN BEGIN
        PRINT, "**Far SOL cannot include both x-points"
        PRINT, " Probably caused by intersection with boundaries."
        
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "**Re-running, removing the strict keyword"
        ENDIF ELSE BEGIN
          
          IF N_ELEMENTS(boundary[0,*]) GT 4 THEN BEGIN
            PRINT, "**Re-running, simplifying boundary"
            boundary = TRANSPOSE([ [MIN(boundary[0,*]), MAX(boundary[0,*]), MAX(boundary[0,*]), MIN(boundary[0,*])], $
                                   [MIN(boundary[1,*]), MIN(boundary[1,*]), MAX(boundary[1,*]), MAX(boundary[1,*])] ])
          ENDIF ELSE BEGIN
            PRINT, "**Re-running, removing boundary"
            boundary = 0
          ENDELSE
        ENDELSE
          
        IF nrad_flexible THEN nrad = TOTAL(nrad) ; Allow nrad to change again

        new_settings = {psi_inner:psi_inner, psi_outer:(max(xpt_psi)+0.02), $
                        nrad:nrad, npol:settings.npol, $
                        rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}
        RETURN, create_grid(F, R, Z, new_settings, critical=critical, $
                            boundary=boundary, iter=iter+1)
      ENDIF
      sol_psi_vals[i,(TOTAL(nrad)-nsol):*] = radial_grid(nsol, $
                                                         xpt_psi_max, psi_outer[solid], $
                                                         0, 1, [xpt_psi_max], $
                                                         settings.rad_peaking)

      tmp = {xpt1:ci[i], xpt2:ci[(i+1) MOD critical.n_xpoint], $ ; X-point indices
             ri:ri, zi:zi, length:length, $  ; R and Z points and the line length
             nsol:nsol}
      
      sol_info[i] = PTR_NEW(tmp)

      ; Add this SOL to the PF structure
      (*(pf_info[ci[i]])).sol[i MOD 2] = i
      (*(pf_info[ci[(i+1) MOD critical.n_xpoint]])).sol[i MOD 2] = i
    ENDFOR
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Do the same for each PF region

    FOR i=0, critical.n_xpoint-1 DO BEGIN
      ; Change grid in PF regions (inner psi can vary)
      ; psi_inner sorted in separatrix psi from inside out
      xind = si[i]
      
      ; Get number of points in this PF region (npf)
      w = WHERE(psi_vals LT xpt_psi[xind], npf)
      w = WHERE(ci EQ xind)
      id = w[0]
      pf_psi_vals[xind,0,0:(npf-1)] = radial_grid(npf, psi_inner[id+1], xpt_psi[xind], 1, 0, $
                                                  [xpt_psi[xind]], settings.rad_peaking)
      pf_psi_vals[xind,1,0:(npf-1)] = pf_psi_vals[xind,0,0:(npf-1)]
      
      
      ; Get the lines from the x-point to the target plates
      legsep = leg_separatrix(dctF, R, Z, xpt_ri[i], xpt_zi[i], $
                              opt_ri[primary_opt], opt_zi[primary_opt], boundary=bndryi)
      
      pf_ri = [REVERSE(legsep.leg1[*,0]), xpt_ri[i], legsep.leg2[*,0]]
      pf_zi = [REVERSE(legsep.leg1[*,1]), xpt_zi[i], legsep.leg2[*,1]]
      mini = N_ELEMENTS(legsep.leg1[*,0])
      
      ; Use the tangent vector to determine direction
      ; relative to core and so get direction of positive theta
      
      drdi = INTERPOLATE((DERIV(INTERPOLATE(R, pf_ri))), mini)
      dzdi = INTERPOLATE((DERIV(INTERPOLATE(Z, pf_zi))), mini)
      
      IF drdi * (*pf_info[xind]).drdi + dzdi * (*pf_info[xind]).dzdi GT 0.0 THEN BEGIN
        ; Line is parallel to the core. Need to reverse
        pf_ri = REVERSE(pf_ri)
        pf_zi = REVERSE(pf_zi)
        mini = N_ELEMENTS(pf_ri) - 1. - mini
      ENDIF

      ; Make sure that the PF sections point to the correct SOLs
      w = WHERE(ci EQ xind)
      i1 = w[0]
      IF (*pf_info[xind]).sol[0] NE w[0] THEN BEGIN
        (*pf_info[xind]).sol = REVERSE((*pf_info[xind]).sol)
      ENDIF
      
      ; Copy radial grid in the SOL
      solind = (*pf_info[xind]).sol[0]
      nsol = (*sol_info[solind]).nsol
      pf_psi_vals[xind, 0, (TOTAL(nrad)-nsol):*] = sol_psi_vals[solind, (TOTAL(nrad)-nsol):*]
      
      solind = (*pf_info[xind]).sol[1]
      nsol = (*sol_info[solind]).nsol
      pf_psi_vals[xind, 1, (TOTAL(nrad)-nsol):*] = sol_psi_vals[solind, (TOTAL(nrad)-nsol):*]

      ; Put the starting line into the pf_info structure
      tmp = CREATE_STRUCT(*(pf_info[xind]), $
                          'npf', npf, $ ; Number of radial points in this PF region
                          'ri0', [pf_ri[0:mini], INTERPOLATE(pf_ri, mini)], $
                          'zi0', [pf_zi[0:mini], INTERPOLATE(pf_zi, mini)], $
                          'ri1', [INTERPOLATE(pf_ri, mini), pf_ri[(mini+1):*]], $
                          'zi1', [INTERPOLATE(pf_zi, mini), pf_zi[(mini+1):*]])

      ; Calculate length of each section
      dldi = SQRT(DERIV(INTERPOLATE(R, tmp.ri0))^2 + DERIV(INTERPOLATE(Z, tmp.zi0))^2)
      len0 = INT_TABULATED(FINDGEN(N_ELEMENTS(dldi)), dldi)
      dldi = SQRT(DERIV(INTERPOLATE(R, tmp.ri1))^2 + DERIV(INTERPOLATE(Z, tmp.zi1))^2)
      len1 = INT_TABULATED(FINDGEN(N_ELEMENTS(dldi)), dldi)
      tmp = CREATE_STRUCT(tmp, 'len0', len0, 'len1', len1)
      
      PTR_FREE, pf_info[xind]
      pf_info[xind] = PTR_NEW(tmp)
    ENDFOR
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Decide poloidal resolutions
    ; If automatic, divide up based on length, with a
    ; minimum of 2 points per region
    
    npol = settings.npol
    nnpol = N_ELEMENTS(npol)
    IF nnpol NE 3*critical.n_xpoint THEN BEGIN
      IF nnpol GT 1 THEN BEGIN
        PRINT, "WARNING: npol has wrong number of elements ("+STR(nnpol)+")"
        PRINT, "         Should have 1 or "+STR(3*critical.n_xpoint)+" elements"
        npol = TOTAL(npol)
      ENDIF
      
      IF npol LT 6*critical.n_xpoint THEN BEGIN
        PRINT, "WARNING: supplied npol ("+STR(npol)+") too small"
        npol = 6*critical.n_xpoint
        PRINT, "   => Increasing npol to "+ STR(npol)
      ENDIF
    
      nnpol = npol
      npol = LONARR(3*critical.n_xpoint) + 2
      nnpol = nnpol - 6*critical.n_xpoint ; Extra points to divide up
      
      ; Get lengths
      length = FLTARR(3*critical.n_xpoint)
      FOR i=0, critical.n_xpoint-1 DO BEGIN
        ; PF regions
        length[i]                     = (*pf_info[i]).len0
        length[critical.n_xpoint + i] = (*pf_info[i]).len1
        ; SOL
        length[2*critical.n_xpoint + i] = (*sol_info[i]).length
      ENDFOR
      
      FOR i=0, nnpol-1 DO BEGIN
        ; Add an extra point to the longest length
        
        dl = length / FLOAT(npol)
        dl[0:(2*critical.n_xpoint-1)] = length[0:(2*critical.n_xpoint-1)] $
          / (FLOAT(npol[0:(2*critical.n_xpoint-1)]) - 0.5)
        
        m = MAX(dl, ind)
        npol[ind] = npol[ind] + 1
      ENDFOR
      
      ; Now sort out the order of npol. Starts from innermost xpoint
      ; and goes clockwise around: PF, SOL, PF, PF, SOL, PF
      
      npol2 = npol
      xpt = si[0] ; X-point index to start with
      FOR i=0, critical.n_xpoint-1 DO BEGIN
        npol2[3*i] = npol[xpt]  ; PF part 0
        
        ; Get the SOL ID
        solid = (*pf_info[xpt]).sol[0]
        IF (*sol_info[solid]).xpt1 NE xpt THEN BEGIN
          PRINT, "ERROR: Indexing inconsistent => Bug!"
          STOP
        ENDIF
        npol2[3*i+1] = npol[2*critical.n_xpoint + solid]
        
        ; Get the next x-point
        xpt = (*sol_info[solid]).xpt2
        npol2[3*i+2] = npol[critical.n_xpoint + xpt] ; PF part 1
      ENDFOR
      
      npol = npol2
    ENDIF
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Now have all the radial locations for each region
    ; and number of poloidal points
    ; 
    ; => Grid each section
    ; 
    
    IF KEYWORD_SET(boundary) THEN BEGIN
      IF KEYWORD_SET(strictbndry) THEN BEGIN
        ; No part of the grid is allowed to be outside boundary
        PRINT, "Keeping the grid stictly within the boundary"
        gridbndry = bndryi
      ENDIF ELSE BEGIN
        ; Grid can leave boundary
        gridbndry = FLTARR(2,4)
        gridbndry[0,*] = [0, 0, nx-1, nx-1]
        gridbndry[1,*] = [0, ny-1, ny-1, 0]
      ENDELSE
    ENDIF
    
    ; Create 2D arrays for the grid
    Rxy = FLTARR(TOTAL(nrad), TOTAL(npol))
    Zxy = Rxy
    Rixy = Rxy
    Zixy = Rxy
    Psixy = Rxy
    
    xpt = si[0] ; Start with the innermost x-point
    ypos = 0
    rerun = 0   ; Flag. If 1 then have to re-run the grid generator
    FOR i=0, critical.n_xpoint-1 DO BEGIN
      ; Calculate maximum psi of x-point
      xpt_psi_max = MAX([xpt_psi[ci[i]], xpt_psi[ci[(i+1) MOD critical.n_xpoint]]])
      PRINT, "Gridding regions "+STR(3*i)+" to " +STR(3*i+2)
      ; Grid the lower PF region
      PRINT, "   x-point index ", xpt
      a = grid_region(dctF, R, Z, $
                      (*pf_info[xpt]).ri0, (*pf_info[xpt]).zi0, $
                      faxis + fnorm*pf_psi_vals[xpt,0,*], $
                      (*pf_info[xpt]).npf-1, $
                      npol[3*i], $
                      sfirst=sfirst1, $
                      slast=slast1, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast1, fpsi=fpsi)
      Rxy[*, ypos:(ypos+npol[3*i]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i]-1 DO Psixy[*, j] = pf_psi_vals[xpt,0,*]
      ypos = ypos + npol[3*i]
      
      IF (flast1 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          ; Only if haven't already marked as restarting
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast1 = faxis + fnorm*MAX(pf_psi_vals[xpt,0,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4

      IF sfirst1 GT 0 THEN BEGIN
        PRINT, "  PF region "+STR(3*i)+" incomplete"
        PRINT, "    => Inner only good to f = "+STR(ffirst)
        pfirst = (ffirst - faxis)/fnorm
        PRINT, "    => Inner normalised psi = "+STR(pfirst)
        
        ; Change the inner psi
        w = WHERE(ci EQ xpt)
        IF pfirst GT psi_inner[w[0]+1] THEN BEGIN
          psi_inner[w[0]+1] = pfirst
        ENDIF
        rerun = 1  ; Signal that the grid needs to be rerun
      ENDIF

      ; SOL region
      solid = (*pf_info[xpt]).sol[0]
      
      PRINT, "   SOL index ", solid

      a = grid_region(dctF, R, Z, $
                      (*sol_info[solid]).ri, (*sol_info[solid]).zi, $
                      faxis + fnorm*sol_psi_vals[solid,*], $
                      nrad[0]-1, $
                      npol[3*i+1], $
                      sfirst=sfirst2, $
                      slast=slast2, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast2, fpsi=fpsi)
      Rxy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i+1]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i+1]-1 DO Psixy[*, j] = pf_psi_vals[xpt,0,*]
      ypos = ypos + npol[3*i+1]

      IF (flast2 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast2 = faxis + fnorm*MAX(sol_psi_vals[solid,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4
      
      ; Second PF region
      xpt = (*sol_info[solid]).xpt2
      PRINT, "   x-point index ", xpt
      a = grid_region(dctF, R, Z, $
                      (*pf_info[xpt]).ri1, (*pf_info[xpt]).zi1, $
                      faxis + fnorm*pf_psi_vals[xpt,1,*], $
                      (*pf_info[xpt]).npf-1, $
                      npol[3*i+2], $
                      sfirst=sfirst3, $
                      slast=slast3, $
                      boundary=gridbndry, $
                      ffirst=ffirst, flast=flast3, fpsi=fpsi)
      Rxy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Rxy
      Zxy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Zxy
      Rixy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Rixy
      Zixy[*, ypos:(ypos+npol[3*i+2]-1)] = a.Zixy
      FOR j=ypos, ypos+npol[3*i+2]-1 DO Psixy[*, j] = pf_psi_vals[xpt,0,*]
      ypos = ypos + npol[3*i+2]

      IF (flast3 - faxis)/fnorm LT xpt_psi_max THEN BEGIN
        PRINT, "WARNING: Due to intersections with the boundary,"
        PRINT, "         the SOL can't cover both x-points"
        IF KEYWORD_SET(strictbndry) THEN BEGIN
          PRINT, "** Switching off strict boundary"
          strictbndry = 0
        ENDIF ELSE IF rerun EQ 0 THEN BEGIN
          PRINT, "** Removing the boundary"
          boundary = 0
        ENDIF
        
        flast2 = faxis + fnorm*MAX(pf_psi_vals[xpt,1,*]) ; don't restrict range
        rerun = 1
      ENDIF

      plot_grid_section, a, color=4
      
      IF sfirst3 GT 0 THEN BEGIN
        PRINT, "  PF region "+STR(3*i+2)+" incomplete"
        PRINT, "    => Inner only good to f = "+STR(ffirst)
        pfirst = (ffirst - faxis)/fnorm
        PRINT, "    => Inner normalised psi = "+STR(pfirst)
        
        ; Change the inner psi. Make sure hasn't already been reduced
        w = WHERE(ci EQ xpt)
        IF pfirst GT psi_inner[w[0]+1] THEN BEGIN
          psi_inner[w[0]+1] = pfirst
        ENDIF
        rerun = 1  ; Signal that the grid needs to be rerun
      ENDIF
      
      ; Check if any of the outer edges failed
      slast = MIN([slast1, slast2, slast3])
      IF slast LT (TOTAL(nrad)-1) THEN BEGIN
        plast = MIN(([flast1, flast2, flast3] - faxis) / fnorm)
        flast = faxis + fnorm*plast
        PRINT, "   => These regions only good until f = "+STR(flast)
        PRINT, "   => i.e. Normalised psi of "+STR( plast )
        
        ; Set the new outer psi for these regions
        w = WHERE(ci EQ primary_xpt)
        id = (solid - w[0] + critical.n_xpoint) MOD critical.n_xpoint
        psi_outer[id] = plast
        PRINT, "SETTING PSI_OUT", id, solid
        rerun = 1
      ENDIF
    ENDFOR

    IF rerun THEN BEGIN
      ; Create new settings
      
      IF nrad_flexible THEN nrad = TOTAL(nrad)  ; Allow nrad to change again

      new_settings = {psi_inner:psi_inner, psi_outer:psi_outer, $ ; New ranges
                      nrad:nrad, npol:npol, $
                      rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}
      
      PRINT, "** Re-running grid generator with changed settings"
      PRINT, "psi outer = ", psi_outer
      RETURN, create_grid(F, R, Z, new_settings, critical=critical, $
                          boundary=boundary, strictbndry=strictbndry, iter=iter+1)
      
    ENDIF
    
    ; Successfully created grid
    
    new_settings = {psi_inner:psi_inner, psi_outer:psi_outer, $ ; New ranges
                    nrad:nrad, npol:npol, $
                    rad_peaking:settings.rad_peaking, pol_peaking:settings.pol_peaking}

    result = {error:0, $ ; Signals success
              psi_inner:psi_inner, psi_outer:psi_outer, $ ; Range of psi
              nrad:nrad, npol:npol, $  ; Number of points in each domain
              Rixy:Rixy, Zixy:Zixy, $  ; Indices into R and Z of each point
              Rxy:Rxy, Zxy:Zxy, $ ; Locations (not really necessary)
              psixy:psixy, $ ; Normalised psi for each point
              faxis:faxis, fnorm:fnorm, $ ; Psi normalisation factors
              settings:new_settings, $ ; Settings used to create grid
              critical:critical}  ; Critical points
  ENDELSE

  CATCH, /cancel

  RETURN, result
END
