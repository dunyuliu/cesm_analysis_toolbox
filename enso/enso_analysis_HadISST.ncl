;********************************************************
;Need to load the following modules on Lonestar6 at TACC at UT Austin.
;module load ncl_ncarg nco
;********************************************************

;********************************************************
;Created on 09/21/2022. 
;The script has the following 4 functions:
; 1. Read in SST time series and calculate Nino 3.4 Index.
; 2. Calculate SST climatology and SST climotology anomolies. 
; 3. Perform EOF analysis on SST_Nino_anomolies, and produce EOF patterns, PC time series, and PC scatter plots.  
;*********************************************************

begin
; =================================================================
; Specify input data and a few parameters.
; =================================================================
	;yr = (/401, 500/) ; specify the year range to be computed.
	;f=addfile("b.e12.12ka.001.hose.pop_rgrd.SST.040101-050012.nc","r") ;# read in SST data.
    	yr = (/1870, 2010/); for HARD
	f  = addfile("HadISST_sst.nc","r")	
	
	neof = 2 ;number of EOFs

	latS = -10
	latN = 10
	lonL = 140 ;140E/
	lonR = 280 ;280E/80W

	mn = (yr(0:1)-yr(0))*12 + (/0,11/)
	sst  = f->sst(mn(0):mn(1),:,:) ; Load SST arrays to variable sst.
	sst  = lonFlip(sst) 
	time = sst&time
	lat  = sst&latitude
	lon  = sst&longitude
	
; =================================================================
; Compute climatology of SST, its anomolies, and the variables in the Nino 3.4 region.
; =================================================================
	sst_climatology      = clmMonLLT(sst(latitude|:,longitude|:,{time|:})) ;# Compute SST climatology
	sst_clim_anom        = calcMonAnomLLT(sst(latitude|:,longitude|:,time|:),sst_climatology) ;# compute the sst anomoly.
 	sst_nino34_clim_anom = dim_avg_Wrap(dim_avg_Wrap(sst_clim_anom(time|:,{latitude|latS:latN},{longitude|lonL:lonR}))) ;# compute average monthly SST anomalies within the Nino 3.4 region.
	sst_nino34           = dim_avg_Wrap(dim_avg_Wrap(sst(time|:,{latitude|latS:latN},{longitude|lonL:lonR}))) ;# compute monthly SST from Nino 3.4 region.
	sst_nino34_anom      = sst_clim_anom(time|:,{latitude|latS:latN},{longitude|lonL:lonR})
	
	nummonth = 12
	sst_nino34_annual_smoothed = runave_Wrap(sst_nino34,nummonth,0) ;# Running nummonth mean SST of Nino 3.4.
	
	rad      = atan(1.)/45.
    	clat     = lat
    	clat     = sqrt(cos( rad*clat)) ;# sqrt is better for EOF analysis.
    	clat!0   = "lat"
    	clat&lat = lat
	
	wgt_sst_nino34_clim_anom = wgt_areaave_Wrap(sst_clim_anom(time|:,{latitude|latS:latN},{longitude|lonL:lonR}), clat({latS:latN}), 1., 0)
	wgt_sst_nino34 = wgt_areaave_Wrap(sst(time|:,{latitude|latS:latN},{longitude|lonL:lonR}), clat({latS:latN}), 1., 0)
	
	printVarSummary(time)
	printVarSummary(sst_nino34_clim_anom)
	printVarSummary(wgt_sst_nino34_clim_anom)
	
	print("Finish calculating SST climotology and anomolies ... ...")

; =================================================================
; Figure 0. Plot Nino 3.4 Index.
; =================================================================
	wks = gsn_open_wks("pdf","sst_nino34_anom, smooth vs weighted")
	;plot_stack= new((/2,dimsizes(time)/),float)
	;plot_stack(0,:) = sst_nino34_clim_anom
	;plot_stack(1,:) = wgt_sst_nino34_clim_anom
	res = True
	;res@trYMinF = -10.0
	;res@trYMaxF = 500.0
	;res@trXMinF = yr(0)
	;res@trXMaxF = yr(1)
	;res@vpHeightF= 0.1;0.1
	;res@vpWidthF = 0.8
	plot = new(2, "graphic")
	res@tiYAxisString = "Degrees Celsius"
	res@gsnXYAboveFillColors = "red"
	res@gsnXYBelowFillColors = "blue"
	res@gsnYRefLine           = 0.4
	res@gsnYRefLine           = -0.4
	res@xyLineColor = "blue"
	res@xyLineThicknesses = (/  1.0,   2.0/)          ; make second line thicker
	res@xyLineColors      = (/"blue","red"/)          ; change line color
	data = new((/2,dimsizes(time)/),float)
	data(0,:) = sst_nino34_clim_anom
	data(1,:) = wgt_sst_nino34_clim_anom
	plot = gsn_csm_xy(wks,time,data,res)
	frame(wks)
	
	print("Finish plotting Figure 0 - Nino 3.4 Index ... ...")
	
; =================================================================
; EOF analysis.
; =================================================================
	;neof         = 4        ; number of EOFs
	optEOF       = True       
	optEOF@jopt  = 0   ; This is the default; most commonly used; no need to specify.
	optETS       = True
	optETS@jopt  = 1
	
	xw = sst_nino34_anom(latitude|:,longitude|:,time|:) ;# reorder to make 'time' the last dim.
	
	eof      = eofunc_Wrap(xw, neof, optEOF)      
	eof_ts   = eofunc_ts_Wrap (xw, eof, optETS)
	
	printVarSummary( eof )                         ; examine EOF variables
	printVarSummary( eof_ts )

	eof_ts   = dim_standardize( eof_ts, 0)         ; normalize by the rightmost dim, which is the highest order of EOF time series.
; =================================================================
; Regress
; =================================================================
	eof_regress = eof
	do ne = 0, neof-1
		eof_regress(ne,:,:) = (/ regCoef(eof_ts(ne,:), xw) /)
	end do

; =================================================================
; Normalize time series: Sum spatial weights over the area of used
; =================================================================
	;dimxw  = dimsizes( xw )
	;mln    = dimxw(1)
	;sumWgt = mln*sum( clat({lat|latS:latN}) )
	;print(sumWgt)
	;eof_ts1 = eof_ts/sumWgt
	;printVarSummary( eof_ts1 )
;============================================================
; PLOTS
;============================================================
	wks = gsn_open_wks("png","enso_eof_regress")         ; send graphics to PNG file
	plot_eof = new(neof,graphic)                ; create graphic array
                                          ; only needed if paneling
; EOF patterns

	res                      = True         
	res@gsnDraw              = False        ; don't draw yet
	res@gsnFrame             = False        ; don't advance frame yet
	res@gsnAddCyclic         = False        ; plotted data are not cyclic
	
	res@mpCenterLonF         = 180.
	;res@mpFillOn             = False        ; turn off map fill
	res@mpMinLatF            = latS         ; zoom in on map
	res@mpMaxLatF            = latN
	res@mpMinLonF            = lonL
	res@mpMaxLonF            = lonR
	;res@tiYAxisString        = " "
	res@cnFillOn             = True         ; turn on color fill
	res@cnLinesOn            = True        ; True is default
	res@cnLineLabelsOn       = True        ; True is default
	;res@cnFillPalette        = "BlWhRe"     ; set color map
	res@cnFillPalette        = "posneg_1"
	res@lbLabelBarOn         = False        ; turn off individual lb's
										  ; set symmetric plot min/max
	;symMinMaxPlt(eof, 16, False, res)       ; contributed.ncl

	; panel plot only resources
	resP                     = True         ; modify the panel plot
	resP@gsnMaximize         = True         ; large format
	resP@gsnPanelLabelBar    = True         ; add common colorbar
	
  	;yStrt                    = yyyymm(0)/100
  	;yLast                    = yyyymm(nyrs-1)/100
  	yStrt = yr(0)
  	yLast = yr(1)
	resP@gsnPanelMainString  = "EOF patterns"

;*******************************************
; Figure 1. EOF pattern plots.
;*******************************************
	do n=0,neof-1
		res@gsnLeftString  = "EOF "+(n+1)
		res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
		plot_eof(n)        = gsn_csm_contour_map(wks,eof_regress(n,:,:),res)
	end do
	gsn_panel(wks,plot_eof,(/neof,1/),resP)     ; now draw as one plot

;*******************************************
; Figure 2. EOF time series plot.
;*******************************************
	eof_ts@long_name        = "Amplitude"
	wks = gsn_open_wks("png","enso_eof_ts")         ; send graphics to PNG file
	plot_eof_ts = new(neof,graphic)                ; create graphic array

	; EOF time series  [bar form]
	rts           = True
	rts@gsnDraw   = False       ; don't draw yet
	rts@gsnFrame  = False       ; don't advance frame yet
	rts@gsnScale  = True        ; force text scaling               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

	rts@vpHeightF = 0.40        ; Changes the aspect ratio
	rts@vpWidthF  = 0.85
	rts@vpXF      = 0.10        ; change start locations
	rts@vpYF      = 0.75        ; the plot

	;rts@tiYAxisString = "Pa"                    ; y-axis label      

	rts@gsnYRefLine           = 0.              ; reference line   
	rts@gsnXYBarChart         = True            ; create bar chart 
	rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
	rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

	; panel plot only resources
	rtsP                      = True            ; modify the panel plot
	rtsP@gsnMaximize          = True            ; large format
	;rtsP@gsnPanelMainString   = "SLP: "+season+": "+yStrt+"-"+yLast

	;year = yyyymm/100

	; create individual plots
	do n=0,neof-1
		rts@gsnLeftString  = "EOF "+(n+1)
		rts@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
		plot_eof_ts(n) = gsn_csm_xy (wks,time,eof_ts(n,:),rts)
	end do
	gsn_panel(wks,plot_eof_ts,(/neof,1/),rtsP)     ; now draw as one plot
  
;*******************************************
; Figure 3. Scatter plot of PC1 and PC2.
;*******************************************  
	wks = gsn_open_wks("png","enso_pc_scatter")     ; send graphics to PNG file
	rcs                   = True                     ; plot mods desired
	rcs@tiMainString      = "PC1-PC2 scatter plot"   ; add title
	rcs@tiXAxisString        = "PC1"  
	rcs@tiYAxisString        = "PC2"
	rcs@xyMarkLineModes   = "Markers"                ; choose which have markers
	rcs@xyMarkers         =  16                      ; choose type of marker  
	rcs@xyMarkerColor     = "red"                    ; Marker color
	rcs@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
	rcs@trYMinF           = -6
	rcs@trYMaxF	      = 4
	rcs@trXMinF 	      = -4
	rcs@trXMaxF           = 6
	rcs@tmLabelAutoStride = True                     ; nice tick mark labels
	
	plot_pc_scatter = gsn_csm_xy (wks,eof_ts(0,:), eof_ts(1,:), rcs)
end
