FUNCTION dst_bkg_fourier_interp, bg_, bg_info, info, str_indx

; NAME:
;	dst_bkg_fourier_interp
; PURPOSE:
;	compute the background estimation for a given set of images starting from some not 
;	saturated ones.
; EXPLANATION:
;	This technique starts from a set of not saturated images during the time interval under
;	analysis. The method starts computing the FFT for the not saturated mages and 
;	interpolating their filtered components for the saturated time frames.
;	Inverting the FFT we obtain an estimation of the background intensity over the image FOV. 
; CALLING SEQUENCE:
;	result = dst_bkg_fourier_interp(bg_, bg_info, info, str_indx, psf)
; INPUTS:
;	bg_		 = data array containing not saturated images corresponding to *str.bg_ind	
;	bg_info	 = index structure of the *str.bg_ind images
;	info	 = index for the saturated images *str.sat_ind
;	str_indx = dst_str global structure for the desaturation routine.	
;
;	OUTPUT:
;	the corresponding background estimation for saturated images.
; CALLS:
;	CALLED BY:
;		data_restore_2
;	CALLS:
;		interpol
; PROCEDURE:
;	

	ct_bg  = size(/n_e, *str_indx.bg_ind)
	ct_sat = size(/n_e, *str_indx.sat_ind)
	
	dim = [bg_info[0].naxis1, bg_info[0].naxis2]

	fbg_ = make_array(dim[0], dim[1], ct_bg, /complex)

	for i=0, ct_bg-1 do fbg_[*, *, i] = fft(bg_[*, *, i], -1)

	filter = BUTTERWORTH( dim[0], dim[1], ct_bg , CUTOFF=100, order = 4 )

	ffbg_  = filter * fbg_

	x_orig = anytim(bg_info.t_obs)-anytim(bg_info[0].t_obs)
	x_intr = anytim(info.t_obs)-anytim(bg_info[0].t_obs)

	ij = array_indices(filter[*, *, 0],where(filter[*, *, 0] gt 0.1, n_indx))

	interp = make_array(dim[0], dim[1], ct_sat, /complex)

	for i=0L, n_indx-1 do begin 

		r_orig = real_part(reform(ffbg_[ij[0, i], ij[1, i], *]))
		i_orig = imaginary(reform(ffbg_[ij[0, i], ij[1, i], *]))

		rres = interpol(r_orig, x_orig, x_intr, /quadratic)
		ires = interpol(i_orig, x_orig, x_intr, /quadratic)

		interp[ij[0, i], ij[1, i], *] = complex(rres, ires)

	endfor
	
	bg = make_array(dim[0], dim[1], ct_sat)
		
	for i=0, ct_sat-1 do begin 
		interp[*, *, i] = (interp[*, *, i]*filter[*, *, 0]) + (fft(bg_[*, *, 0],-1)*(1-filter[*, *, 0]))
		bg[*, *, i]= real_part(fft(interp[*, *, i],1))>0
	endfor
	
	return, bg

end
