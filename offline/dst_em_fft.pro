pro dst_em_fft , str , cpsf , it=it , level=level 

; NAME:
;	dst_em_fft
; PURPOSE:
;	Deconvolution/desaturation routine based on EM method implemented by means the FFT.
; EXPLANATION:
;	This routine has a multiple usage in the desaturation method:
;		1) deconvolution method routine for the not saturated image;
;		2) compute the correlation product for determining primary saturated pixels; 
;		3) desaturation routine in the case of large saturated areas ( > 1000 );
; CALLING SEQUENCE:
;	dst_em_fft , str , cpsf , it=it , level=level , pow=pow
; INPUTS:
;	str		= dst_str global structure for the desaturation routine.	
;	cpsf	= diffraction component of the PSF
; OPTIONAL:
;	it		= maximum iteration number for EM method (default 300)
;	level	= EM stopping rule parameter (default 0.01)
; OUTPUT:	
;	str.x 		= retrieved intensities for saturated pixels
;	str.y 		= array corresponding to the diffraction fringes intensities
;	str.c_exp	= expectation computed by the method (cpsf*x + bg at the last iteration)
;	str.c_stat 	= c_statistic computed between y and c_exp					
; CALLS:
;	CALLED BY:
;		DST_X_Y
;		DECONV
;		DESATURATION
;	CALLS:
;		CONVOLVE, C_SATATISTIC, AVG
;
	
default , it , str.it
default , level , str.lev
default , pow , 1

s   = *str.s
g	= *str.g

;;;ZEROPADING THE DATA ARRAY TO AN N*N ARRAY WHER N IS A POWER OF 2 (INCREASE FFT SPEED) 
dim = size(/dim,*str.im)
check = alog(dim)/alog(2)
if round(check[0]) ne check[0] or round(check[1]) ne check[1] then dim = [2^(floor(check[0])+1),2^(floor(check[1])+1)]

im 		= fltarr(dim)
bkg 	= fltarr(dim)
mask_g  = fltarr(dim)
x		= fltarr(dim)
;;;

emp_kkt = fltarr(it+1)
exp_kkt = fltarr(it+1)

mask = im * 0 & mask[dim[0]/2,dim[1]/2] = 1

im[0,0]  = *str.im
bkg[0,0] = *str.bg

cpsf = convolve(mask , cpsf)>0
cpsf2 = cpsf * cpsf

x_tmp 	 = *str.im*0. 
x_tmp[s] = total(*str.bg_bkup) gt 0 ? (*str.bg_bkup)[s] : 1

x[0,0] 	 = x_tmp

sat_pos  = where(f_div(x , x) gt 0)

mask_g_4_zeropad = *str.im * 0. & mask_g_4_zeropad[g] = 1.
mask_g[0,0] = mask_g_4_zeropad

y 	= ( im  * mask_g )>0.
bkg = ( bkg * mask_g )>0.				 ;;; GT modification

V = convolve(mask_g , cpsf , /corr )>0.

for I = 0, it-1 do begin

	y_i = (convolve( x , cpsf ) + bkg)>0
	
	Z = (f_div(y , y_i)) * mask_g

	U = convolve( z , cpsf , /corr )>0.
	x[sat_pos] *= ( f_div( U , V ) )[sat_pos]

	emp_kkt[i] = total( ( x * (V - U))^2.)
	exp_kkt[i] = total( x^2. * (convolve( f_div(1. , y_i * mask_g)>0 , cpsf2 , /corr )>0.) )

	if emp_kkt[i] le level * exp_kkt[i] then break
	
	print, i ,' --- ', emp_kkt[i]/exp_kkt[i] , level
	
endfor

dim = size(/dim , *str.im)

x   = x[0:dim[0]-1 , 0:dim[1]-1] 
y   = y[0:dim[0]-1 , 0:dim[1]-1] 
y_i = y_i[0:dim[0]-1 , 0:dim[1]-1] 
 
;;; C_STATISTIC COMPUTED ON A SUBSET OF THE DATA
c_stat = fltarr(10)
n_samp = size(/n_e,g)

for iii = 0, size(/n_e,c_stat)-1 do begin 
	ind = randomu(seed,100)*n_samp
	ind = g[ind[UNIQ(ind, SORT(ind))]]
	C_stat[iii] = c_statistic( y[ind] , y_i[ind] ) 
endfor

ind = 220+findgen(50) 
pmm,c_statistic( y[ind] , y_i[ind] ) 

C_stat = AVG(C_stat)
;;;

print , i , '  R_flx:' , total( y_i[g]) , '  Exp_flx:' , total( y[g] )  , '  C_stat:' , c_stat , ' ' 

str.x = ptr_new(x[s])
str.y = ptr_new(y)
str.c_exp  = ptr_new(y_i[g])
str.c_stat = c_stat					

end
