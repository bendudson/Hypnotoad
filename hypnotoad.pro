;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                  HYPNO-TOAD grid generator
;
; 
; This is a graphical interface to some grid generation routines.
; Aims to allow tokamak grids to be easily generated from
; a variety of input sources.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


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
        
        ; Extract needed data from g-file struct
        
        rz_grid = {nr:g.nx, nz:g.ny, $  ; Number of grid points
                   r:REFORM(g.r[*,0]), z:REFORM(g.z[0,*]), $  ; R and Z as 1D arrays
                   simagx:g.simagx, sibdry:g.sibdry, $ ; Range of psi
                   psi:g.psi, $  ; Poloidal flux in Weber/rad on grid points
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
        
        ; Set info to new values
        widget_control, event.top, set_UVALUE=info
      ENDIF ELSE BEGIN
        ; Couldn't read data
        PRINT, "ERROR: Failed to read grid file"
      ENDELSE
    END
    'mesh': BEGIN
      ; Create a mesh
      IF info.rz_grid_valid EQ 0 THEN BEGIN
        PRINT, "ERROR: No valid equilibrium data. Read from file first"
        a = DIALOG_MESSAGE("No valid equilibrium data. Read from file first", /error)
        RETURN
      ENDIF
      boundary=fltarr(2,4)
      boundary[0,*] = [1.0, 1.0, 2.5, 2.5]
      boundary[1,*] = [-1.4, 1.4, 1.4, -1.4]
      
      
      mesh = create_grid((*(info.rz_grid)).psi, (*(info.rz_grid)).r, (*(info.rz_grid)).z, $
                         boundary=boundary, /strict)
      IF mesh.error EQ 0 THEN BEGIN
        PRINT, "Successfully generated mesh"
      ENDIF ELSE BEGIN
        a = DIALOG_MESSAGE("Could not generate mesh", /error)
      ENDELSE
    END
    'detail': BEGIN
      ; Checkbox with detail settings
      IF event.select THEN BEGIN
        ; Detailed settings
        
        ; Hide the global settings
        WIDGET_CONTROL, info.nrad_field, map = 0
        WIDGET_CONTROL, info.npol_field, map = 0
      ENDIF ELSE BEGIN
        ; Simple settings
        
        ; Show global settings
        WIDGET_CONTROL, info.nrad_field, map = 1
        WIDGET_CONTROL, info.npol_field, map = 1
      ENDELSE
    END
    ELSE: PRINT, "Unknown event", uvalue
  ENDCASE
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main procedure

PRO hypnotoad
  
  ; Make IDL retain a backing store
  DEVICE, retain=2

  ; Create the main window
  base = WIDGET_BASE(title="Hypnotoad", mbar=mbar, /ROW, $
                     EVENT_PRO = 'event_handler')

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
  
  read_button = WIDGET_BUTTON(bar, VALUE='Read G-eqdsk', $
                              uvalue='aandg', tooltip="Read RZ equilibrium from EFIT")

  nrad_field = CW_FIELD( bar,                           $
                         title  = 'Radial points:',      $ 
                         uvalue = 'nrad',                $ 
                         /long,                          $ 
                         value = 36,                     $
                         xsize=8,                        $
                         /return_events)
  npol_field = CW_FIELD( bar,                           $
                         title  = 'Poloidal points:',    $ 
                         uvalue = 'npol',                $ 
                         /long,                          $ 
                         value = 64,                     $
                         xsize=8,                        $
                         /return_events)
  
  mesh_button = WIDGET_BUTTON(bar, VALUE='Generate mesh', $
                              uvalue='mesh', tooltip="Generate a new mesh")

  checkboxbase = WIDGET_BASE(bar, /ROW, EVENT_PRO = 'event_handler', /NonExclusive)
  region_check = WIDGET_BUTTON(checkboxbase, VALUE="Detailed settings", uvalue='detail', $
                               tooltip="Set parameters for each region separately")

  ; Create an area for drawing
  draw = WIDGET_DRAW(base, xsize=400, ysize=400)

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
           draw:draw $ ; Drawing widget
         } 

  ; Store this in the base UVALUE
  WIDGET_CONTROL, base, set_uvalue=info 

  ; Draw everything
  WIDGET_CONTROL, base, /real

  XMANAGER, 'hypnotoad', base, /no_block
END
