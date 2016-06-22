pro dst_sqm , str , cpsf 
; NAME:
;	dst_sqm
; PURPOSE:
;	It computes the sqm matrix for dst_em routine
; EXPLANATION:
;	select the components of the convolution matrix that links the saturated pixels str.s 
;	with the diffraction fringes in str.g
; CALLING SEQUENCE:
;	dst_sqm , str , cpsf
; INPUTS: 
;	str	= dst_str global structure for the desaturation routine.	\
;	cps	= diffraction component of the PSF 
; OUTPUT:
;	str.sqm = sqm matrix for dst_em routine
;	str.y	= vector containing the diffraction fringes intensities
;	str.bg	= vector containing the estimated background intensities for the same pixels 
;			  considered in y
; CALLS:
;	CALLED BY:
;		DESATURATION
;	

s     = *str.s
g     = *str.g
n_sat = str.ns

dim		= size(/dim , *str.im)
dim_psf = size(/dim , cpsf)
n_g		= size(/n_e , g)

sqm = fltarr(n_g , n_sat)

ij_sat = array_indices(*str.im , s)

im1 = *str.im * 0.0 & im1[dim[0]/2, dim[1]/2] = 1 & c0 = convolve(/corr,im1, cpsf) > 0.0
for i=0, n_sat-1 do sqm[*,i] = (shift(c0 , ij_sat[0,i] - dim[0]/2 , ij_sat[1,i] - dim[1]/2 ))[g] 

str.sqm = ptr_new(sqm)
str.y	= ptr_new((*str.im)[g] > 0.)
str.bg	= ptr_new((*str.bg)[g] > 0.)

end

