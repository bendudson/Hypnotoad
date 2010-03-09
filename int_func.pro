; integrate a function, always using the maximum
; number of grid-points possible for highest accuracy

FUNCTION int_func, x, f
   n = N_ELEMENTS(f)

   n2 = FIX(n/2)

   g = FLTARR(n)

   g[0] = 0.0
   FOR i=n2, n-1 DO BEGIN
       g[i] = INT_TABULATED(x[0:i], f[0:i])
   ENDFOR

   FOR i=1, n2-1 DO BEGIN
       g[i] = g[n-1] - INT_TABULATED(x[i:*], f[i:*])
   ENDFOR

   RETURN, g
END
