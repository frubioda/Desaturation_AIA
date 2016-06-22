pro demo_for_desaturation, obj, result, tstart = tstart, tend=tend, wav = wav, path_save=path_save ,	$
						   sat_lev=sat_lev , lev=lev , it=it , npix=npix, core_dim=core_dim ,			$
						   save_fts=save_fts , aec = aec , loud=loud

; NAME:
;	demo_for_desaturation
; PURPOSE:
;	Calls the desaturation routine 
; CALLING SEQUENCE:
;	demo_for_desaturation, obj, results
; INPUTS:
;	tstart		= user defined analysis starting time
;	tend		= user defined analysis end time
;	wav 		= string array containing the wavelength to be process 
; OUTPUT:
;	obj 		= 'desat' object for desaturated images.
;	result 		= Structure containing the desaturated images
;   				DATA            POINTER   <PtrHeapVar4635> 
;   				INFO            POINTER   <PtrHeapVar4636>
; OPTIONAL:
;	sat_level	= saturation intensity level (default = 16300) 
;	lev			= EM stopping rule parameter (default = 0.01)
;	it			= maximum iteration number (default = 300)
;	core_dim	= radius for the PSF core (default = 5)
;	npix		= setting the resulting image to be npix*npix (default = 499)
; KEYWORDS:
;   path_save	= path to folder where .fts file will be saved
;	save_fts	= if set save fts file at the end of the run (default = 1)
;	aec 		= if set will desaturate also sat. images with short exp. time (default = 1)
;	loud		= if 1 intermiediate and final results are shown in windows (default = 1)
; CALLS:
;	CALLED BY: - 
;	CALLS:	
;		readfts
;		indices_analysis
;		savefts
;		data_restore_2
;		aia_psf_gen
;		dst_x_y
;		dst_sqm
;		dst_EM
;		dst_em_FFT
;		dst_disp_factor_2
; PROCEDURE:
;

;;; Demo for the desaturation procedure for AIA images 
cd , c = curdir
default , tstart 	, '06-Sep-2011 22:16:00'		; start time selection
default , tend  	, '06-Sep-2011 22:20:00'		; end time selection
default , path		, curdir						; data storage folder path
default , path_save , curdir						; data storage folder path
default , sat_lev	, 15000
default , lev		, 0.01
default , it		, 300
default , n_pix		, 499
default , core_dim	, 5
default , save_fts	, 1
default , aec		, 1
default , loud		, 1
default , wav		, ['94','131','171','193','211','304','335']

;;; Cutout data request line.
;ssw_cutout_service,'2011/09/06 22:00:00','2011/09/06 22:30:00',ref_helio='N15W18',fovx=350,fovy=350,email='torre@dima.unige.it',waves=[94,131,171,193,211,304,335],max_frames=1000,instrument='aia',aec=1 

tstart = '06-Sep-2011 22:16:00'					
tend   = '06-Sep-2011 22:20:00'					

path = '~/SDO_data/20110906_2'					; data storage folder path

obj = obj_new('desat')

result = obj -> desaturation( wav , tstart , tend , path , path_save = path_save , sat_lev = sat_lev ,$
							  lev = lev , it = it , npix = npix, core_dim = core_dim ,$
							  save_fts = save_fts , aec = aec , loud = loud )
end