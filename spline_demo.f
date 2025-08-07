      PROGRAM SPLINE_DEMO
      INTEGER X(100), Y(100), N
      INTEGER IX, IY

      CALL ps_init_window()

      N = 0
10    CONTINUE
        PRINT *, 'Click control point (or press Enter to quit)...'
        CALL ps_wait_click(IX, IY)
        N = N + 1
        X(N) = IX
        Y(N) = IY
	IF (IX .LT. 0 .OR. IY .LT. 0) THEN
	    PRINT *, 'Exiting...'
	    GOTO 20
	ENDIF

        CALL ps_draw_point(IX, IY)
        
        IF (N .GE. 4) THEN
            PRINT *, 'Drawing spline through', N, 'points'
            CALL ps_draw_spline(X, Y, N)
        ENDIF

        GOTO 10

20    END

