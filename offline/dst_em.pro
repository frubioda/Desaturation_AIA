pro dst_em, str ,level=level, it=it 
;+
; NAME:
;	dst_em
; PURPOSE:
;	Computes the intensities of saturated pixels by mean EM algorithm.
; EXPLANATION:
;	It computes the intensities for the primary saturated pixels by means EM algorithm 
;	implemented using the convolution matrix (sqm). This routine is used for small saturated
;	region (x < 1000).
; CALLING SEQUENCE:
;	dst_em, str ,level=level, it=it 
; INPUTS:
;	str		= dst_str global structure for the desaturation routine.	
;	cpsf	= diffraction component of the psf
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
;		DESATURATION
;-
	
y   = *str.y 
sqm = *str.sqm 
bg  = *str.bg 

default , it , str.it
default , level , str.lev

dim_sqm = size(/dim , sqm)

emp_kkt = fltarr(it)
exp_kkt = fltarr(it)
sqm2	= sqm * sqm

;x = fltarr(str.ns) + 1.
;x = (*str.bg_bkup)[*str.s]
x = total((*str.bg_bkup)[*str.s]) gt 0 ? (*str.bg_bkup)[*str.s] : fltarr(str.ns) + 1.  ;;; GT modification

V = ((y * 0.) + 1.) # sqm 

for i = 0 , it-1 do begin 

	y_i = isarray(x) eq 1 ? sqm # x + bg : sqm * x + bg
	z = f_div( y , y_i )

	U = reform( z # sqm )

	x = x * f_div( U , V )

	emp_kkt[i] = total( ( x * (V - U) )^2. )
	exp_kkt[i] = total(  x^2. * ( f_div(1.,y_i) # sqm2 ) ) 
	
	if emp_kkt[i] le level * exp_kkt[i] then break
	
	print, i ,' --- ', emp_kkt[i]/exp_kkt[i] , level
	
endfor

c_stat = fltarr(10)
n_samp = size(/n_e,y)

for iii = 0, size(/n_e,c_stat)-1 do begin 
	ind = randomu(seed,100)*n_samp
	ind = UNIQ(ind, SORT(ind))
	C_stat[iii] = c_statistic( y[ind] , y_i[ind] ) 
endfor

C_stat = avg(c_stat)

print , i , '  R_flx:' , total( y_i ) , '  Exp_flx:' , total( y )  , '  C_stat:' , c_stat , ' ' 
 
str.x = ptr_new(x)
str.c_exp = ptr_new(y_i)
str.c_stat = C_stat

end

