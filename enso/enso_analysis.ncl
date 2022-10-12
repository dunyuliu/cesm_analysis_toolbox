;********************************************************
;Need to load the following modules on Lonestar6 at TACC at UT Austin.
;module load ncl_ncarg nco
;********************************************************

;********************************************************
;Created on 09/21/2022. 
;The script has the following 4 functions:
; 1. Read in SST time series and calculate Nino 3.4 Index.
; 2. Calculate SST climatology and SST anomolies. 
; 3. Perform EOF analysis on SST anomolies and generate EOF patterns.
; 4. Plot PC scatter plots based on EOF time series.  
;*********************************************************

begin
; =================================================================
; Specify input data and a few parameters.
; This block contains most adjustable parameters. 
; =================================================================
	case            = 1
	; 1: HadISST_sst.nc
	; 2: b.e12.0kaHose02SvNA.SST.030101-090012.nc
	; 3: b.e12.3ka.SST.010101-060012.nc
	; 4: b.e12.6ka.SST.010101-060012.nc
	; 5: b.e12.9ka.SST.040101-090012.nc
	; 6: b.e12.12ka.SST.030101-080012.nc
	; 7: b.e12.12kaHose01SvNA.SST.030101-090012.nc
	; 8: b.e12.12kaHose02SvNA.SST.030101-090012.nc
	; 9: b.e12.15ka.SST.020101-070012.nc
	; 10: b.e12.15kaHose02SvNA.SST.030101-090012.nc
	; 11: b.e12.18ka.SST.020101-070012.nc
	; 12: b.e12.21ka.SST.050101-100012.nc
	; 13: b.e12.21kaHose02SvNA.SST.030101-060012.nc
	; 14: b.e12.2XCO2.SST.030101-090012.nc
	; 15: b.e12.4XCO2.SST.030101-090012.nc
	;
	;b.e12.6ka.SST.010101-060012.nc
	;==============================================================
	if (case.eq.1) then 
		; case 1. Reproduction of Takahashi results on HadISST dataset.
		yr              = (/1870, 2010/)                            ; the whole year range for HadISST dataset.
		yrc             = (/1979, 2009/)                            ; year range to calculate climotology following Takahashi et al. (2011).
		fname           = "HadISST_sst.nc"
		f               = addfile(fname,"r")
	;==============================================================
	elseif (case.eq.2) then 
		; case 2. 0ka.hose.030101-090012
		yr              = (/301, 900/)
		yrc             = yr
		fname           = "b.e12.0kaHose02SvNA.SST.030101-090012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	elseif (case.eq.3) then 
		; case 3. b.e12.3ka.SST.010101-060012.nc
		yr              = (/101, 600/)
		yrc             = yr
		fname           = "b.e12.3ka.SST.010101-060012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	elseif (case.eq.4) then 
		; case 4. b.e12.6ka.SST.010101-060012.nc
		yr              = (/101, 600/)
		yrc             = yr
		fname           = "b.e12.6ka.SST.010101-060012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	elseif (case.eq.5) then 
		; case 5. b.e12.9ka.SST.040101-090012.nc
		yr              = (/401, 900/)
		yrc             = yr
		fname           = "b.e12.9ka.SST.040101-090012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	elseif (case.eq.6) then 
		; case 6. b.e12.12ka.SST.030101-080012.nc
		yr              = (/301, 800/)
		yrc             = yr
		fname           = "b.e12.12ka.SST.030101-080012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	elseif (case.eq.7) then 
		; case 7. b.e12.12kaHose01SvNA.SST.030101-090012.nc
		yr              = (/301, 900/)
		yrc             = yr
		fname           = "b.e12.12kaHose01SvNA.SST.030101-090012.nc"
		f               = addfile("/disk/cg6/yukoo/cesm1_timeslice/" + fname, "r")
	end if
	
	fontsize        = 0.01                                 ; font size used in the EOF pattern plots. 
	
	neof            = 2                                        ; number of EOFs

	latS            = -10                                      ; map range in lat and lon. 
	latN            = 10
	lonL            = 110                                      ;140E/. HadISST dataset lon in 0-360. 
	lonR            = 280                                      ;280E/80W

	mn              = (yr(0:1)-yr(0))*12 + (/0,11/)               ; time in month index for the whole HadISST data from yr 1870 to 2010.
	mnc             = (/(yrc(0)-yr(0))*12,(yrc(1)-yr(0))*12+11/) ; time in month index for the whole HadISST data from yr 1979 to 2009.
	
	nmonths_apr_sep = (yr(1)-yr(0)+1)*6 
	mn_apr_sep      = ispan(0,nmonths_apr_sep-1, 1)
	mn_oct_mar      = mn_apr_sep

	itag0           = 0
	itag1           = 0
	do i = 0,yr(1)-yr(0)                            ; loop over all the years and allocate month indexes for Arpil to September and October to March.
		do j = 0,11                                 ; loop over 12 months.
			if ((j.ge.3).and.(j.le.8)) then
				mn_apr_sep(itag0) = i*12 + j
				itag0             = itag0 + 1
			else
				mn_oct_mar(itag1) = i*12 + j 
				itag1             = itag1 +1 
			end if 
		end do 
	end do 
	
	if (case.eq.1) then
		sst             = f->sst(mn(0):mn(1),:,:) ; Load SST arrays to variable sst.
		lat             = sst&latitude
		sst             = lonFlip(sst) 
		sst             = where(sst.eq.-1000., -1.8, sst)
	else 
		sst             = f->SST(mn(0):mn(1),:,:) ; Load SST arrays to variable sst.
		sst!1="lat"
		sst!2="lon"
		lat             = sst&lat
	end if
	time            = fspan(yr(0)+0.5/12, yr(1)+11.5/12, (yr(1)-yr(0)+1)*12)
	
; =================================================================
; Compute climatology of SST, its anomolies, and the variables in the Nino 3.4 region.
; =================================================================
	sst_climatology      = clmMonTLL(sst(mnc(0):mnc(1),:,:)) ;# Compute SST climatology
	sst_clim_anom        = calcMonAnomTLL(sst,sst_climatology) ;# compute the sst anomoly.
	;sst_clim_anom        = dtrend_leftdim(sst_clim_anom, False) ; detrend the whole dataset.
	
	;nummonth = 12
	;sst_nino34_annual_smoothed = runave_Wrap(sst_nino34,nummonth,0) ;# Running nummonth mean SST of Nino 3.4.
	
	rad                  = atan(1.)/45.
	clat                 = lat
	clat                 = sqrt(cos( rad*clat)) ;# sqrt is better for EOF analysis.
	clat!0               = "lat"
	clat&lat             = lat

	print("Finish calculating SST climotology and anomolies ... ...")

; =================================================================
; Figure 0. Plot Nino 3.4 Index.
; =================================================================
	wks                  = gsn_open_wks("pdf","sst_nino34_anom, smooth vs weighted")
	res                  = True
	plot                 = new(2, "graphic")
	res@tiYAxisString    = "Degrees Celsius"
	res@gsnXYAboveFillColors = "red"
	res@gsnXYBelowFillColors = "blue"
	res@gsnYRefLine      = 0.4
	res@gsnYRefLine      = -0.4
	res@xyLineColor      = "blue"
	res@xyLineThicknesses = (/  1.0,   2.0/)          ; make second line thicker
	res@xyLineColors      = (/"blue","red"/)          ; change line color
	data                 = new((/2,dimsizes(time)/),float)
	;data(0,:) = sst_nino34_clim_anom
	;data(1,:) = wgt_sst_nino34_clim_anom
	;plot = gsn_csm_xy(wks,time,data,res)
	frame(wks)
	
	print("Finish plotting Figure 0 - Nino 3.4 Index ... ...")
	
; =================================================================
; EOF analysis.
; =================================================================
	optEOF       = True       
	optEOF@jopt  = 0   ; This is the default; most commonly used; no need to specify.
	optETS       = True
	optETS@jopt  = 0   ; if we calculate EOFs using a covariance matrix.
	
	xw           = sst_clim_anom(mnc(0):mnc(1),{latS:latN},{lonL:lonR})
	xw_all       = sst_clim_anom(:,{latS:latN},{lonL:lonR})     ; xw and xw_all are only differed for the HadISST set.
	
	;xw           = dtrend_leftdim(xw0, False)
	;xw_all       = dtrend_leftdim(xw_all0, False)
	
	;printVarSummary(xw0)
	printVarSummary(xw)
	
	eof          = eofunc_n_Wrap(xw, neof, optEOF, 0)                     ; eof is computed with selected time periods for the HadISST set.

	if (ismissing(eof(0,{0.},{270.}) .lt. 0.) ) then; reverse eof if eof1 is negative in the eastern equatorial Pacific.
		print("EOF1 (0,270)<0 missing ... ... Try (0,255) ... ...")
		if (eof(0,{0.},{255.}) .lt. 0.) then
			eof(0,:,:) = -eof(0,:,:)
		end if
	else 
		if(eof(0,{0.},{270.}) .lt. 0.) then
			eof(0,:,:) = -eof(0,:,:)
		end if
	end if 
	if (ismissing(eof(1,{0.},{270.}) .gt. 0.) ) then; reverse eof if eof1 is negative in the eastern equatorial Pacific.
		print("EOF2 (0,270)>0 missing ... ... Try (0,255) ... ...")
		if (eof(1,{0.},{255.}) .gt. 0.) then
			eof(1,:,:) = -eof(1,:,:)
		end if
	else 
		if(eof(1,{0.},{270.}) .gt. 0.) then
			eof(1,:,:) = -eof(1,:,:)
		end if
	end if 
	eof_ts       = eofunc_ts_n_Wrap (xw_all, eof, optETS, 0)              ; then, the whole dataset is projected to the calculated eof patterns to get PC time series.
	printVarSummary( eof )                         ; examine EOF variables
	printVarSummary( eof_ts )

	eof_ts       = dim_standardize( eof_ts, 0)         ; normalize by the rightmost dim, which is the highest order of EOF time series.
; =================================================================
; Regress
; =================================================================
	;eof_regress  = eof
	;eof_regress  = (/ regCoef_n(eof_ts, xw_all, 1, 0) /)
	
	eof_regress   = sst_clim_anom(0:neof-1,:,:)
	eof_regress  = (/ regCoef_n(eof_ts, sst_clim_anom, 1, 0) /)

;*******************************************
; Figure 1. EOF pattern plots.
;*******************************************
	wks          = gsn_open_wks("pdf", fname + ".enso_eof_regress")         ; send graphics to PNG file
	plot_eof     = new(neof,graphic)                ; create graphic array
                                          ; only needed if paneling
; EOF patterns
	res                      = True         
	res@gsnDraw              = False        ; don't draw yet
	res@gsnFrame             = False        ; don't advance frame yet
	res@gsnAddCyclic         = True        ; plotted data are not cyclic
	res@gsnStringFontHeightF = fontsize
	res@gsnPolarLableFontHeightF = fontsize
	res@tmXBLabelFontHeightF = fontsize
	res@tmYLLabelFontHeightF = fontsize
	
	res@mpCenterLonF         = 210.
	;res@mpFillOn            = False        ; turn off map fill
	;res@mpMinLatF            = latS         ; zoom in on map
	;res@mpMaxLatF            = latN
	;res@mpMinLonF            = lonL
	;res@mpMaxLonF            = lonR
	res@mpLandFillColor      = "gray70"
	
	res@tiYAxisString        = " "          ; no y label.
	
	res@cnFillOn             = True         ; turn on color fill
	res@cnLinesOn            = True        ; True is default
	res@cnLineLabelsOn       = True        ; True is default
	;res@cnFillPalette       = "BlWhRe"     ; set color map
	;res@cnFillPalette       = "amwg_blueyellowred"
	res@cnFillPalette        = "temp_19lev"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels             = fspan(-0.9,0.9,19)
	res@lbLabelBarOn         = False       ; turn off individual lb's
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
	resP@gsnPanelMainString   = "EOF patterns"
	
	;============================================================
	; add custom label bar
	;============================================================
	cmap                       = read_colormap_file( "temp_19lev")   ;20 colors
	col                        = cmap 
	ncol                       = dimsizes( col(:,0))
	lev                        = fspan( -0.1*9, 0.1*9, ncol-1)
	lab = sprintf( "%3.1f", lev)
	lab(0::2) = ""
	
	lbres                      = True
	lbres@lbPerimOn            = False               ; no label bar box
	lbres@lbOrientation        = "Horizontal"        ; orientation
	lbres@vpWidthF             = 0.5                 ; size
	lbres@vpHeightF            = 0.042
	lbres@lbLabelFontHeightF   = 0.01                ; label font height
	lbres@lbLabelAlignment     = "InteriorEdges"     ; where to label
	lbres@lbMonoFillPattern    = True                ; fill sold
	lbres@lbFillColors         = col      ; must be RGB triplets
	
	; To customize the color levels and labels.
	;levels = (/-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7/)
	;labels = (/"-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7"/)
	;nlevels = dimsizes(levels) 
	
	;posx = (/0.02, 0.02/)
	;posy = (/0.95, 0.75/)
	posx = (/0.05, 0.05/)
	posy = (/0.95, 0.48/)
	res@vpWidthF              = 0.90
	res@vpHeightF             = 0.35
	;res@vpKeepAspect         = False
	do n = 0, neof-1
		res@vpXF   = posx(n)
		res@vpYF   = posy(n)
		res@gsnLeftString     = "EOF "+(n+1)
		res@gsnRightString    = sprintf("%5.1f", eof@pcvar(n)) +"%"
		plot_eof(n)           = gsn_csm_contour_map(wks,eof_regress(n,:,:),res)
	end do
	; plot customized colobar.
	;gsn_labelbar_ndc( wks, ncol, lab, 0.5-lbres@vpWidthF/2, 0.57, lbres)
	gsn_labelbar_ndc( wks, ncol, lab, 0.5-lbres@vpWidthF/2, 0.06, lbres)
	draw(wks)
	frame(wks)

;*******************************************
; Figure 2. EOF time series plot.
;*******************************************
	eof_ts@long_name        = "Amplitude"
	wks                     = gsn_open_wks("pdf","enso_eof_ts")         ; send graphics to PNG file
	plot_eof_ts             = new(neof,graphic)                ; create graphic array

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
	wks                   = gsn_open_wks("pdf", fname + ".enso_pc_scatter")     ; send graphics to PNG file
	rcs                   = True                     ; plot mods desired
	rcs@tiMainString      = "PC1-PC2 scatter plot"   ; add title
	rcs@tiXAxisString     = "PC1"  
	rcs@tiYAxisString     = "PC2"
	rcs@xyMarkLineModes   = "Markers"                ; choose which have markers
	rcs@xyMarkers         =  16                      ; choose type of marker  
	rcs@xyMarkerColor     = "red"                    ; Marker color
	rcs@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
	rcs@trYMinF           = -6
	rcs@trYMaxF	          = 4
	rcs@trXMinF 	      = -4
	rcs@trXMaxF           = 6
	rcs@tmLabelAutoStride = True                     ; nice tick mark labels
	
	;plot_pc_scatter      = gsn_csm_xy (wks,eof_ts(0,:), eof_ts(1,:), rcs)
	plot_pc_scatter_apr_sep = gsn_csm_xy (wks,eof_ts(0,mn_apr_sep), eof_ts(1,mn_apr_sep), rcs)
	
	rcs@xyMarkerColor     = "blue"
	plot_pc_scatter_oct_mar = gsn_csm_xy (wks,eof_ts(0,mn_oct_mar), eof_ts(1,mn_oct_mar), rcs)
	
	overlay(plot_pc_scatter_apr_sep, plot_pc_scatter_oct_mar)
	
	draw( wks)
	frame( wks)
end
