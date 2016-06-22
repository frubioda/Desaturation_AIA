function path_folder_concat, path
	concat_path = path[0]
	for i = 1 , size(/n_e, path)-1 do concat_path = concat_dir(concat_path,path[i])
return, concat_path 
end

pro flare_peak_wavelength, passband, peaklam = peaklam, printpk = printpk, plotall = plotall 
;Project     : Diffraction  
;                   
; Name        : flare_peak_wavelength
;               
; Purpose     : To identify the brightest line in a passband during a flare 
;               
; Use         : IDL> flare_peak_wavelength, 171, peaklam = peaklam    
;                       This will identify the brightest line
;						in the 131 A passband and return it as peaklam
;    
; Inputs      : passband = AIA passband wavelength (integer)
;               
; Outputs	  : peaklam = Wavelength of brightest line in Angstroms
;                              
; Keywords    : printpk = if set then print the peak wavelength to the screen 
;  			  : plotall = plot the spectrum with the brightest line marked 
;
; Written     : Claire L. Raftery, SSL, Nov 2015      


;;retrieve the AIA response functions 
aia_responses=aia_get_response(/eff,/dn)

;;Define which abundance, ionization fraction and dem file to use
;;I've been using coronal abundances and Mazotta ionizations

path = [get_logenv('SSW_CHIANTI') , 'dbase' , 'ioneq' , 'chianti.ioneq']
ioneq_name = path_folder_concat(path)

path = [get_logenv('SSW_CHIANTI') , 'dbase' , 'abundance' , 'sun_coronal_1999_fludra_ext.abund']
abund_name = path_folder_concat(path)

path = [get_logenv('SSW_CHIANTI') , 'dbase' , 'dem' , 'flare.dem']
dem_name = path_folder_concat(path)

;;Define denisty at which spectrum is calcuated 
den = 1e10

;Identify the passband (and the effective area while we're at it)
case passband of
94: begin
	lambda = aia_responses.a94.wave
	effarea = aia_responses.a94.ea
	end
131:begin
	lambda = aia_responses.a131.wave
	effarea = aia_responses.a131.ea
	end
171:begin
	lambda = aia_responses.a171.wave
	effarea = aia_responses.a171.ea
	end
193:begin
	lambda = aia_responses.a193.wave
	effarea = aia_responses.a193.ea
	end
211:begin
	lambda = aia_responses.a211.wave
	effarea = aia_responses.a211.ea
	end
304:begin
	lambda = aia_responses.a304.wave
	effarea = aia_responses.a304.ea
	end
335:begin
	lambda = aia_responses.a335.wave
	effarea = aia_responses.a335.ea
	end
else: begin
	lambda = aia_responses.wave
	effarea = aia_responses.all 
	end
endcase


help, lambda
pmm, lambda
help, effarea 

;;Identify the spectral components for a given DEM function for the 
;;wavelength range we're analyzing.
ch_synthetic, min(lambda), max(lambda), output=str, density=den, $
/all,/photons, ioneq_name = ioneq_name, dem_name = dem_name, err_msg=ms

make_chianti_spec, str, lambda, spec, /cont, abund = abund_name,/photons;, bin_size = bin;, err_msg = em
scale_spec = spec.spectrum * effarea

peaklam = lambda[where(scale_spec eq max(scale_spec))]

if keyword_set(printpk) then print, 'Wavelength of brightest line = '+peaklam+' Angstroms'
if keyword_set(plotall) then begin
	set_line_color
	plot, lambda, scale_spec, xtitle = 'Wavelength ['+spec.wvl_units+']', ytitle = 'Intensity ['+spec.int_units+']'
	verline, peaklam, col = 3   
	oplot, lambda, scale_spec
endif 

done:
end




















