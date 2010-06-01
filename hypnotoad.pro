;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                  HYPNO-TOAD grid generator
;
; 
; This is a graphical interface to some grid generation routines.
; Aims to allow tokamak grids to be easily generated from
; a variety of input sources.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_region, R, Z, ymin, ymax, _extra=_extra
  s = SIZE(R, /dimen)
  nx = s[0]
  ny = s[1]
  
  FOR j=0, nx-1 DO BEGIN
    OPLOT, R[j,ymin:ymax], Z[j,ymin:ymax], _extra=_extra
  ENDFOR
  
  FOR j=ymin, ymax DO BEGIN
    OPLOT, R[*,j], Z[*,j], _extra=_extra
  ENDFOR
END

PRO plot_rz_equil, data
  CONTOUR, data.psi, data.r, data.z, /iso, nlev=40
END


PRO oplot_mesh, rz_mesh, flux_mesh
  ; Plot X-points and separatrices
  FOR i=0, flux_mesh.critical.n_xpoint-1 DO BEGIN
    ; plot the separatrix contour
    CONTOUR, rz_mesh.psi, rz_mesh.R, rz_mesh.Z, levels=[flux_mesh.critical.xpt_f[i]], c_colors=2, /overplot
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.xpt_ri[i])], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.xpt_zi[i])], psym=7, color=2
  ENDFOR

  ; Plot O-points
  FOR i=0, flux_mesh.critical.n_opoint-1 DO BEGIN
    oplot, [INTERPOLATE(rz_mesh.R, flux_mesh.critical.opt_ri[i])], [INTERPOLATE(rz_mesh.Z, flux_mesh.critical.opt_zi[i])], psym=7, color=3
  ENDFOR
  
  ypos = 0
  FOR i=0, N_ELEMENTS(flux_mesh.npol)-1 DO BEGIN
    plot_region, flux_mesh.rxy, flux_mesh.zxy, ypos, ypos+flux_mesh.npol[i]-1, color=i+1
    ypos = ypos + flux_mesh.npol[i]
  ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Event handling procedure

PRO event_handler, event
  ; Get the UVALUE
  widget_control, event.id, get_uvalue=uvalue
  
  ; Retrieve a copy of information stored in tlb
  widget_control, event.top, get_uvalue=info
  
  IF N_ELEMENTS(uvalue) EQ 0 THEN RETURN ; Undefined
  
  CASE uvalue OF
    'aandg': BEGIN
      PRINT, "Open G-eqdsk (neqdsk) file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="neqdsk", /read)
      PRINT, "Trying to read file "+filename
      g = read_neqdsk(filename)
      
      IF SIZE(g, /TYPE) EQ 8 THEN BEGIN
        ; Got a structure
        PRINT, "Successfully read equilibrium"
        WIDGET_CONTROL, info.status, set_value="Successfully read "+filename
        
        ; Extract needed data from g-file struct
        
        rz_grid = {nr:g.nx, nz:g.ny, $  ; Number of grid points
                   r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $  ; R and Z as 1D arrays
                   simagx:g.simagx, sibdry:g.sibdry, $ ; Range of psi
                   psi:g.psi, $  ; Poloidal flux in Weber/rad on grid points
                   npsigrid:(FINDGEN(N_ELEMENTS(g.pres))/N_ELEMENTS(g.pres)), $ ; Normalised psi grid for fpol, pres and qpsi
                   fpol:g.fpol, $ ; Poloidal current function on uniform flux grid
                   pres:g.pres, $ ; Plasma pressure in nt/m^2 on uniform flux grid
                   qpsi:g.qpsi, $ ; q values on uniform flux grid
                   nlim:g.nlim, rlim:g.xlim, zlim:g.ylim} ; Wall boundary
        

        IF info.rz_grid_valid GT 0 THEN BEGIN
          ; Need to free existing data
          PTR_FREE, info.rz_grid
        ENDIF

        ; Put pointer to data into info struct
        info.rz_grid = PTR_NEW(rz_grid)
        info.rz_grid_valid = 1
        
        ; Plot the equilibrium
        plot_rz_equil, rz_grid

        ; Set info to new values
        widget_control, event.top, set_UVALUE=info
      ENDIF ELSE BEGIN
        ; Couldn't read data
        PRINT, "ERROR: Failed to read grid file"
        WIDGET_CONTROL, info.status, set_value="   *** Failed to read G-EQDSK file "+filename+" ***"
      ENDELSE
    END
    'mesh': BEGIN
      ; Create a mesh
      IF info.rz_grid_valid EQ 0 THEN BEGIN
        PRINT, "ERROR: No valid equilibrium data. Read from file first"
        a = DIALOG_MESSAGE("No valid equilibrium data. Read from file first", /error)
        RETURN
      ENDIF
      ;boundary=fltarr(2,4)
      ;boundary[0,*] = [1.0, 1.0, 2.5, 2.5]
      ;boundary[1,*] = [-1.4, 1.4, 1.4, -1.4]
      
      boundary = TRANSPOSE([[(*info.rz_grid).rlim], [(*info.rz_grid).zlim]])
      
      IF info.detail_set THEN BEGIN
        settings = {dummy:0}
      ENDIF ELSE BEGIN
        ; Get settings
        widget_control, info.nrad_field, get_value=nrad
        widget_control, info.npol_field, get_value=npol

        widget_control, info.psi_inner_field, get_value=psi_inner
        widget_control, info.psi_outer_field, get_value=psi_outer
        settings = {nrad:nrad, npol:npol, psi_inner:psi_inner, psi_outer:psi_outer}
      ENDELSE
      
      WIDGET_CONTROL, info.status, set_value="Generating mesh ..."
      
      mesh = create_grid((*(info.rz_grid)).psi, (*(info.rz_grid)).r, (*(info.rz_grid)).z, settings, $
                         boundary=boundary, strict=info.strict_bndry)
      IF mesh.error EQ 0 THEN BEGIN
        PRINT, "Successfully generated mesh"
        WIDGET_CONTROL, info.status, set_value="Successfully generated mesh. All glory to the Hypnotoad!"
        oplot_mesh, *info.rz_grid, mesh
        
        info.flux_mesh_valid = 1
        info.flux_mesh = PTR_NEW(mesh)
        widget_control, event.top, set_UVALUE=info
      ENDIF ELSE BEGIN
        a = DIALOG_MESSAGE("Could not generate mesh", /error)
        WIDGET_CONTROL, info.status, set_value="  *** FAILED to generate mesh ***"
      ENDELSE
    END
    'process': BEGIN
      ; Process mesh to produce output
      PRINT, "Write output file"
      filename = DIALOG_PICKFILE(dialog_parent=event.top, file="bout.grd.nc", $
                                 /write, /overwrite_prompt)
      
      IF info.rz_grid_valid AND info.flux_mesh_valid THEN BEGIN
        process_grid, *(info.rz_grid), *(info.flux_mesh), $
                      output=filename, poorquality=poorquality
        
        IF poorquality THEN BEGIN
          r = DIALOG_MESSAGE("Poor quality equilibrium")
        ENDIF
      ENDIF ELSE BEGIN
        PRINT, "ERROR: Need to generate a mesh first"
      ENDELSE
    END
    'detail': BEGIN
      ; Checkbox with detail settings
      info.detail_set = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'strict': BEGIN
      ; Checkbox with boundary strictness
      info.strict_bndry = event.select
      widget_control, event.top, set_UVALUE=info
    END
    'draw': BEGIN
      IF info.flux_mesh_valid EQ 0 THEN RETURN
      
      pos = CONVERT_COORD(event.x, event.y, /device, /to_data)
      r = pos[0]
      z = pos[1]
      ; (r,z) position where clicked
      
      m = MIN( (REFORM((*info.flux_mesh).rxy - r))^2 + $
               (REFORM((*info.flux_mesh).zxy - r))^2 , ind)
      xi = FIX(ind / TOTAL((*info.flux_mesh).npol))
      yi = FIX(ind MOD TOTAL((*info.flux_mesh).npol))
      PRINT, xi, yi
    END
    ELSE: PRINT, "Unknown event", uvalue
  ENDCASE
END

PRO handle_resize, event
  WIDGET_CONTROL, event.top, get_uvalue=info, /No_Copy
  
  statusgeom = WIDGET_INFO(info.status, /geom) 

  WIDGET_CONTROL, info.draw, $
                  Draw_XSize=(event.x - info.leftbargeom.xsize) > statusgeom.xsize, $
                  Draw_YSize=(event.y - statusgeom.ysize) > info.leftbargeom.ysize
  

  IF info.rz_grid_valid THEN BEGIN
    ; Plot the equilibrium
    plot_rz_equil, *info.rz_grid

     IF info.flux_mesh_valid THEN BEGIN
       ; Overplot the mesh
       mesh = *info.flux_mesh
       oplot_mesh, *info.rz_grid, *info.flux_mesh
     ENDIF
  ENDIF

  Widget_Control, event.top, Set_UValue=info, /No_Copy
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main procedure

PRO hypnotoad
  
  ; Make IDL retain a backing store
  DEVICE, retain=2

  ; Create the main window
  base = WIDGET_BASE(title="Hypnotoad", /ROW, $ ; mbar=mbar
                     EVENT_PRO = 'event_handler', TLB_size_events=1)

  ; Put items in the menu
  ;input_menu = WIDGET_BUTTON(mbar, VALUE='Input', /MENU)
  ;input_bttn1=WIDGET_BUTTON(input_menu, VALUE='Open G-eqdsk (neqdsk)...',$
  ;                          UVALUE='aandg', EVENT_PRO = 'event_handler')
  ;input_bttn2=WIDGET_BUTTON(input_menu, VALUE='From IDAM...',$
  ;                          UVALUE='idam', EVENT_PRO = 'event_handler')
  ;input_bttn2=WIDGET_BUTTON(input_menu, VALUE='Test data...',$
  ;                          UVALUE='test', EVENT_PRO = 'event_handler')

  ; Create a bar down left side for buttons and settings
  bar = WIDGET_BASE(base, /COLUMN, EVENT_PRO = 'event_handler')
  
  read_button = WIDGET_BUTTON(bar, VALUE='Read G-EQDSK', $
                              uvalue='aandg', tooltip="Read RZ equilibrium from EFIT")
  
  nrad_field = CW_FIELD( bar,                            $
                         title  = 'Radial points:',      $ 
                         uvalue = 'nrad',                $ 
                         /long,                          $ 
                         value = 36,                     $
                         xsize=8                         $
                       )
  npol_field = CW_FIELD( bar,                            $
                         title  = 'Poloidal points:',    $ 
                         uvalue = 'npol',                $ 
                         /long,                          $ 
                         value = 64,                     $
                         xsize=8                         $
                       )
  
  psi_inner_field = CW_FIELD( bar,                            $
                              title  = 'Inner psi:',          $ 
                              uvalue = 'inner_psi',           $ 
                              /floating,                      $ 
                              value = 0.9,                    $
                              xsize=8                         $
                            )
  psi_outer_field = CW_FIELD( bar,                            $
                              title  = 'Outer psi:',          $ 
                              uvalue = 'outer_psi',           $ 
                              /floating,                      $ 
                              value = 1.1,                    $
                              xsize=8                         $
                            )
  

  mesh_button = WIDGET_BUTTON(bar, VALUE='Generate mesh', $
                              uvalue='mesh', tooltip="Generate a new mesh")

  checkboxbase = WIDGET_BASE(bar, /COLUMN, EVENT_PRO = 'event_handler', /NonExclusive)
  region_check = WIDGET_BUTTON(checkboxbase, VALUE="Detailed settings", uvalue='detail', $
                               tooltip="Set parameters for each region separately")
  strict_check = WIDGET_BUTTON(checkboxbase, VALUE="Strict boundaries", uvalue='strict', $
                               tooltip="Enforce boundaries strictly")
  Widget_Control, strict_check, Set_Button=1
  
  process_button = WIDGET_BUTTON(bar, VALUE='Output mesh', $
                                 uvalue='process', tooltip="Process mesh and output to file")

  leftbargeom = WIDGET_INFO(bar, /Geometry)

  rightbar = WIDGET_BASE(base, /COLUMN, EVENT_PRO = 'event_handler')
  status_box = WIDGET_TEXT(rightbar)
  ; Create an area for drawing
  draw = WIDGET_DRAW(rightbar, xsize=400, ysize=600, /button_events, $
                     uvalue="draw", EVENT_PRO = 'event_handler')

  widget_control, status_box, set_value="Hypnotoad flux grid generator. Read equilibrium G-EQDSK file to begin"

  ; Create a structure for storing the state
  ; This is shared 

  info = { nrad_field:nrad_field, $ ; nrad input box
           npol_field:npol_field, $ ; npol input box
           rz_grid:(PTRARR(1))[0], $ ; Pointer to R-Z grid data
           rz_grid_valid:0, $  ; Flag to indicate if rz_mesh is valid
           flux_mesh:(PTRARR(1))[0], $ ; Pointer to flux surface aligned mesh
           flux_mesh_valid:0, $
           boundary:(PTRARR(1))[0], $ ; Pointer to boundary array [2,*]
           boundary_valid: 0, $
           settings:(PTRARR(1))[0], $ ; Settings structure
           draw:draw, $ ; Drawing widget
           detail_set:0, $ ; 1 if using detailed settings
           strict_bndry:1, $ ; 1 if boundaries should be strict
           psi_inner_field:psi_inner_field, psi_outer_field:psi_outer_field, $
           status:status_box, $
           leftbargeom:leftbargeom $
         } 

  ; Store this in the base UVALUE
  WIDGET_CONTROL, base, set_uvalue=info 

  ; Draw everything
  WIDGET_CONTROL, base, /real

  XMANAGER, 'hypnotoad', base, /no_block, event_handler='handle_resize'
END
