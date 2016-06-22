function dst_psf_gen, wavelength, npix, dwavelength, core_dim

; NAME:
;	dst_psf_gen
; PURPOSE:
;	returns the psf structure that contains the diffraction and dispersion and the complete
;	aia psf.
; EXPLANATION:
;	The dispersion component of the psf is considered cutting out a circular portion on the
;	central component of the psf computed by aia_calc_psf_mod.pro. The radius of this circle is 
;	fixed to 5 pixels but it can be modified by users.
; CALLING SEQUENCE:
;	psf = dst_psf_gen( info , dwavelength , core_dim )
; INPUTS: 
;	info	  =	index structure for the aia data, necessary to retrieve infos about the psf
;				to generate
;	dwavelength	  = correction parameter on the wavelength of the psf (look aia_calc_psf_mod.pro)
;	core_dim  =	radius of the central core of the psf.
; OUTPUT:
;	psf.cpsf = diffraction component of the psf
;	psf.opsf = dispersion component of the psf
;	psf.psf  = complete psf
; CALLS:
;	CALLED BY:
;		DESATURATION
;	CALLS:
;		aia_calc_psf_mod.pro
;	

	default, core_dim, 5

	psf_in = aia_calc_psf_mod( wavelength, $
	  npix = npix, dwavelength = dwavelength)

	dim = size(/dim,psf_in)

	xgrid	= (fltarr(dim[1])+1)##indgen(dim[0])
	ygrid	= indgen(dim[1])##(fltarr(dim[0])+1)

	center	= [fix(dim[0]/2.),fix(dim[1]/2.)]
	w		= where((xgrid-center[0])^2+(ygrid-center[1])^2 le core_dim^2)

	cpsf	= psf_in*0.
	opsf	= psf_in*0.

	opsf[w] = psf_in[w]
	cpsf    = psf_in - opsf

	psf		= {wavelength: wavelength, dwavelength: dwavelength, npix: npix, cpsf:cpsf, opsf:opsf, psf:psf_in}
	
	return, psf
end	