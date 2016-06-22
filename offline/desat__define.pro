
function desat::find_file, path, wav, range = valid_range

  ; NAME:
  ;	find_file
  ; PURPOSE:
  ;	search all the .fts files in the path folder associated to wav wavelength.
  ; EXPLANATION:
  ;	This routine will find all the .fts files inside the path folder. If will return
  ;	all the files associated to the wav wavelength.
  ; CALLING SEQUENCE:
  ;	result = obj -> find_file(path, '131')
  ; INPUTS:
  ;	pix = string containing the path to the folder containing the aia .fts files
  ;	wav = string array or value for the data wavelength to search

  file_name	= file_search(concat_dir(path,'*.f*ts'), COUNT=ct)
  if ct eq 0 then file_name	= file_search(concat_dir(path,'*.fits'), COUNT=ct)

  mreadfits, file_name, info ;, data ;ras, remove data because you only want to read the headers here
  ind 	= where(info.WAVELNTH eq float(wav))
  info  = info[ ind ]
  file_name = file_name[ ind ]
  ind  = where_within( anytim( info.date_obs ), valid_range )
  info  = info[ ind ]
  file_name = file_name[ ind ]

  return,file_name

end

pro desat::readfts, file

  ; NAME:
  ;	readfts
  ; PURPOSE:
  ;	read the .fts files and fill the object.
  ; EXPLANATION:
  ;	Read the data contained in file filling the object with data and infos
  ; CALLING SEQUENCE:
  ;	obj -> readfts, file
  ; INPUTS:
  ;	file = string or array of string containing the file names to be read.
  ; CALLS:
  ;	deconv
  ; PROCEDURE:
  ;	it is used to read the data from .fts files, by means mreadfits routine, and fill
  ;	the 'desat' object.
  ;

  mreadfits, file, info, data
  index2map, info, data, map
  self-> set, index=info, map=map, /no_copy

end

pro desat::data_restore_2, str, data, bg, info, dwavelength = dwavelength

  ; NAME:
  ;	data_restore_2
  ; PURPOSE:
  ;	returns the data, background and info array for the selected time interval for
  ;	saturated frames.
  ; EXPLANATION:
  ;	it will return the data array for the saturated images inside the selected time
  ;	interval.
  ;	At the same time this routine will return a computed estimation for data
  ;	background by means the dst_bkg_fourier_interp.pro routine. Moreover, the unsaturated
  ;	images inside the selected time interval are deconvolved with the
  ;	Expectation Maximization algorithm.
  ; CALLING SEQUENCE:
  ;	obj -> data_restore_2, str, data, bg, info
  ; INPUTS:
  ;	str  = dst_str global structure for the desaturation routine.
  ; OUTPUT:
  ; 	data = data 3D corresponding to saturated images into the selected time interval
  ;	bg	 = computed background estimation obtained by dst_bkg_fourier_interp routine
  ;	info = index structure associated to data array
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;	CALLS:
  ;		DECONV
  ;		dst_bkg_fourier_interp
  ; PROCEDURE:
  ;

  ;;; cut data array to be npix*npix pixel
  self -> ctrl_data, str.npix

  ;;; data/info array storage
  info	= (self ->get(/index))[*str.sat_ind]
  data    = ((self ->get(/data))>0)[*,*,*str.sat_ind]

  ;;; deconvolution step on un-saturated data
  bg_real = self -> deconv( *str.bg_ind, it = str.it, lev = str.lev, dwavelength=dwavelength,  /fill )

  ;;; background estimation routine
  bg_real_info	= (self -> get(/index))[*str.bg_ind]
  bg = dst_bkg_fourier_interp(bg_real, bg_real_info, info, str)

end

pro desat::indices_analysis, str, aec = aec

  ; NAME:
  ;	indices_analysis
  ; PURPOSE:
  ;	It will return indices for the considering event for saturated, background in the
  ;	considered time interval.
  ; EXPLANATION:
  ;	Analyzing the infos contained into the index structure it is possible to retrive
  ;	the indices for saturated frames, the time frames that will be used for the background
  ;	estimation step.
  ;	Moreover, are defined the indices the selected time interval.
  ; CALLING SEQUENCE:
  ;	obj -> indices_analysis, str, aec = aec
  ; INPUTS:
  ;	str  = dst_str global structure for the desaturation routine.
  ; OUTPUT:
  ;	str.time_ind = selected time interval indices
  ;	str.sat_ind  = saturated frames indices
  ;	str.bg_ind	 = indices for the frame that will be used into the background estimation step
  ; KEYWORDS:
  ;	aec	= if set also saturated images with short exp. time were desaturated
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;

  default, aec, 1

  info = self ->get(/index) ;done after this, we get info from reading all the headers, not the data
  n_el = size(/n_e, info)

  ;;; file indices for the selected time interval
  t_ind = where(anytim(info.date_obs) gt anytim(str.ts) and $
    anytim(info.date_obs) lt anytim(str.te), ct )

  ;;; saturated frame indices
  sat_t = where(info.DATAMAX ge str.sat_lev)

  ;;; short exposure time file indices
  s_exp = where(info.exptime lt 0.8*median(info.exptime))

  ;;; indices for saturated frames in the selected time interval
  if keyword_set(aec) eq 0 then begin
    sat_ind = where( histogram(t_ind,min=0,max=n_el) + histogram(sat_t,min=0,max=n_el) - $
      histogram(s_exp,min=0,max=n_el) eq 2, ct_sat )
  endif else begin
    sat_ind = where( histogram(t_ind,min=0,max=n_file) + histogram(sat_t,min=0,max=n_file) eq 2, ct_sat )
  endelse

  ;;; background indices estimation over the whole dataset
  if isarray(s_exp) eq 1 then begin
    bg_ind_dataset = where( histogram(findgen(n_el)) 						-$
      histogram(sat_t, min=0, max=n_el) 				+$
      histogram(s_exp, min=0, max=n_el) ge 1, ct_bg )
  endif else begin
    bg_ind_dataset = where( histogram(findgen(n_el)) 						-$
      histogram(sat_t, min=0, max=n_el) ge 1, ct_bg )
  endelse

  ;;; background indices estimation in the select time interval
  bg_ind_t_ind = where( histogram(t_ind, min=0, max=n_el) 		+$
    histogram(bg_ind_dataset, min=0, max=n_el) ge 2, ct_bg )

  ;;; background indices for background estimation routine (bg_ind_t_ind-1 < bg_ind < bg_ind_t_ind+1)
  bg_ind = [bg_ind_dataset[where( bg_ind_dataset eq min(bg_ind_t_ind)) -1]	,$
    bg_ind_t_ind														,$
    bg_ind_dataset[where( bg_ind_dataset eq max(bg_ind_t_ind)) +1]	]
  if stregex( /boolean, /fold, str.bkg_method, 'quad' ) then begin ;limit indices to no more than two before sat_ind
    sel = where( min( sat_ind ) - bg_ind le 2, nsel)
    if nsel ge 1 then bg_ind = bg_ind[ sel ] else bg_ind = -1
  endif
  ;limit indices to only two before  sat_ind

  str.time_ind	= ptr_new(t_ind)
  str.sat_ind 	= ptr_new(sat_ind)
  str.bg_ind	= ptr_new(bg_ind)
end

function desat::deconv, index, fill = fill, it = it, lev = lev, dwavelength = dwavelength

  ; NAME:
  ;	deconv
  ; PURPOSE:
  ;	Deconvolution procedure EM based for SDO/AIA images
  ; EXPLANATION:
  ;	Deconvolution were computed with the whole psf over a reduced FOV for SDO/AIA data
  ; CALLING SEQUENCE:
  ;	dst_result = obj -> desat( index, fill=fill )
  ; INPUTS:
  ;	index  = indices of data images in the obj on wich deconvolution has to be performed.
  ; OUTPUT:
  ;	dst_result = deconvolved images
  ; KEYWORDS:
  ;	fill = Deconvolved data were substituted into the obj in index position
  ;	it   = Maximum number of iterations
  ;	lev  = Stopping rule parameter for EM
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;	CALLS:
  ;		dst_em_fft
  ;

  default, it , 1000
  default, lev, 0.1

  info = (self-> get(/index))[0]
  data = (self-> get(/data))[*,*,index]

  data_real = data * 0

  psf	= Self->psf( info, dwavelength, core_dim ) ; dst_psf_gen( info, dwavelength, core_dim )
  dim	= [[info.naxis1],[info.naxis2]]

  ;;; DECONVOLUTION FOR BACKGROUND ESTIMATION
  for i=0,size(/n_e, index)-1 do begin

    print, strtrim(i+1), '	of ', strtrim(size(/n_e, index))

    res	= fltarr(info.naxis1,info.naxis2) + 1

    in  = {s   		: ptr_new(findgen(dim[0]*dim[1])),$
      g   		: ptr_new(findgen(dim[0]*dim[1])),$
      im  		: ptr_new(data[*,*,i])			 ,$
      bg  		: ptr_new(data[*,*,i]*0)		 ,$
      bg_bkup 	: ptr_new(data[*,*,i]*0)		 ,$
      x   		: ptr_new(/allocate_heap)		 ,$
      y   		: ptr_new(/allocate_heap)		 ,$
      c_exp	: ptr_new(/allocate_heap)		 ,$
      c_stat	: 0								 ,$
      it		: it							 ,$
      lev		: lev							 }

    dst_em_fft, in, psf.psf, it=it, level=level

    res = reform(*in.x,dim[0],dim[1])

    data_real[*,*,i] = res						;;;real image (de-convolved one)
    data[*,*,i] = convolve(res, psf.opsf)>0		;;;real image convolution with central core -> un-sat. image restoration

  endfor

  if keyword_set(fill) eq 1 then begin
    dec_info = self-> get(/index)
    dec		 = self-> get(/data)
    dec[*,*,index] = data

    index2map, dec_info, dec, map
    self-> set, index=dec_info, map=map
  endif

  return, data_real

end
pro desat::_bld_filenames,  index, data_path
  ; NAME:
  ; _bld_filenames
  ; PURPOSE:
  ; save data array into an .fts file into defined data_path folder
  ; EXPLANATION:
  ; if original is set it will save original data
  ; if desat is set it will save original data
  ; CALLING SEQUENCE:
  ; obj -> savefts, index, data_path, original=original, desat=desat
  ; INPUTS:
  ; index  = indices of data images in the obj that has to be saved.\
  ; data_path = path of the save folder
  ; info_str = informational structure with details of deconvolution
  ; KEYWORDS:
  ; original = set if original data has to be saved (before the desaturation procedure)
  ; desat = set if desat. data ha to be saved (after the desaturation procedure)
  ; CALLS:
  ; CALLED BY:
  ;   DESATURATION
  ; CALLS:
  ;   mwritefts, mkdir, file_mkdir
  ; PROCEDURE:
  ;

  n_el = size(/n_e,index)

  info = (self -> get(/index))[index]

  stringtime = anytim(info.date_obs,/vms,/date_only)
  stringwave = strtrim(string(info.wavelnth),1)

  path = concat_dir(stringtime,stringwave)
  path = concat_dir(data_path, path)

  path_or  = concat_dir(path,'data')
  path_sat = concat_dir(path,'desat')

  if ~dir_exist( path_or[0] ) then file_mkdir, path_or
  if ~dir_exist( path_sat[0] ) then file_mkdir, path_sat
  ;use time2file in constructing filenames!
  orig_file = concat_dir(path_or, 'aia_orig_'  + time2file( info.date_obs, /sec )+'_'+stringwave+'.fts')
  desat_file = concat_dir(path_sat,'aia_desat_' + time2file( info.date_obs, /sec )+'_'+stringwave+'.fts')
  filenames = replicate( {date_obs:'', orig: '', desat: ''}, n_el)
  filenames.date_obs = info.date_obs
  filenames.orig     = orig_file
  filenames.desat    = desat_file

  self.filenames = ptr_new( filenames )
end

function desat::_Get_Filename, date_obs, original = original,  index = index
  filenames = *Self.filenames
  default, original, 0
  desat = 1 - original
  info = (self -> get(/index))
  data_class = desat ? 'desat' : 'orig'
  tg_id = stregex(/fold, /boo, 'desat', data_class) ? tag_index( filenames, 'DESAT' ) : tag_index( filenames, 'ORIG' )


  if exist( date_obs ) then begin
    select = where(  anytim( filenames.date_obs ) eq anytim( date_obs ), nsel )
    if nsel eq 0 then message,'date_obs does not match any of the filenames!'
    filename = filenames[select].(tg_id)
    filenames[select].(tg_id) = '' ;remove so we don't overwrite the desat image with the original at the last fits write
    *Self.filenames = filenames
  endif else begin
    ;We must want all remaining desat filenames
    filename = filenames.(tg_id)
    ;Only retrieve the remaining valid filenames
    select = where( filename ne '', nsel)

    filename = nsel ge 1 ? filename[ select ] : -1
  endelse
  date_obs_4_index = filenames[select].date_obs
  index = where_arr( info.date_obs, date_obs_4_index )

  return, filename
end

pro desat::savefts, date_obs, original=original, desat = desat, use_prep=use_prep, info_str = info_str

  ; NAME:
  ;	savefts
  ; PURPOSE:
  ;	save data array into an .fts file into defined data_path folder
  ; EXPLANATION:
  ;	if original is set it will save original data
  ;	if desat is set it will save original data
  ; CALLING SEQUENCE:
  ;	obj -> savefts, index, data_path, original=original, desat=desat
  ; INPUTS:
  ;	index  = indices of data images in the obj that has to be saved.\
  ;	data_path = path of the save folder
  ;	info_str = informational structure with details of deconvolution
  ; KEYWORDS:
  ;	original = set if original data has to be saved (before the desaturation procedure)
  ;	desat = set if desat. data ha to be saved (after the desaturation procedure)
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;	CALLS:
  ;		mwritefts, mkdir, file_mkdir
  ; PROCEDURE:
  ;

  default, original, 0
  default, use_prep, 1
  desat = 1 - original

  tmp_file = Self->_Get_Filename( date_obs, original=original, index = index )
  data = (self -> get(/data))[*,*, index]
  info = (self -> get(/index))[index]
  if keyword_set( use_prep ) then begin
    aia_prep, info, data, oindex, odata, /cutout
    info = oindex
    data = odata
    Self->update_history, info, 'desat alg'

  endif
  MWRITEFITS, info, data, outfile=tmp_file
  if exist( info_str ) && is_struct( info_str ) then begin
    info = { sat_core: *info_str.s, $
      sat_fringe: *info_str.g, $
      sat_bloom: *info_str.b, $
      c_stat: info_str.c_stat, $
      ex_max_target: info_str.lev, $
      sat_level: info_str.sat_lev }
    ;mwrfits, info, ssw_strsplit( tmp_file, '.fts' ) + '_desat_info.fts'
    mwrfits, info, tmp_file, /silent
    fxhmodify, tmp_file,'EXTEND',1
    fxhmodify, tmp_file, 'EXTNAME  ','DESAT_INFO', /extension

  endif
  data_path = file_dirname( file_search(tmp_file[0], /full) )
  if keyword_set(original) then print, 'Original .fts file are saved in	' + data_path + '	folder'
  if keyword_set(desat) then print, 'De-saturated .fts file are saved in	' + data_path + '	folder'

end
pro desat::identify_saturation_regions, info, str, psf_str, loud

  ; NAME:
  ; identify_saturation_regions
  ; PURPOSE:
  ; By means of the correlation product computed with dst_em_fft routine return the
  ; position for primary saturated pixels, bloomed pixels, diffraction fringes.
  ; EXPLANATION:
  ; Selecting two regions of the image (I_1 with intensities > 15000 and the I_2 with
  ; intensities < 10000) the correlation into I_1 is computed using I_2 as data. Note that
  ; in our formulation the correlation product correspond to the first iteration of
  ; dst_em_fft procedure. Convolving the result of the cross-correlation with the dispersion
  ; core of the PSF we recognize the primary saturation site S (the one with intensity > saturation
  ; level) and as a consequence the I_2 pixels that are not in S are the bloomed ones B.
  ; Moreover, considering the support of the convolution product between mask_s, defined as
  ; mask_s[s] = 1. and 0 elsewhere, and the psf diffraction component is possible to find out
  ; the diffraction fringes directly generated by the primary saturation site.
  ; CALLING SEQUENCE:
  ; desat::identify_saturation_regions, str, psf_str
  ; INPUTS:
  ; str     = dst_str global structure for the desaturation routine.
  ; psf_str   = structure containing all the psf components, diffraction, dispersion and
  ;       the complete one. This structure is returned as output from dst_psf_gen.pro
  ;       routine.
  ; OUTPUT:
  ; str.s   = primary saturated pixels position
  ; str.z   = flaring pixels (not saturated pixels, with high intensity avoided by the method)
  ; str.g   = diffraction fringes pixels position
  ; str.b   = bloomed pixels positions
  ; str.ns    = number of saturated pixel
  ; str.sat_lev = saturation intensity level
  ; CALLS:
  ; CALLED BY:
  ;   DESATURATION
  ;

  cpsf  = psf_str.cpsf
  im    = *str.im
  bkg   = *str.bg
  bg_bkup = *str.bg_bkup

  ;;; Classify the map intensity in three groups (saturation / fringes / pure background emission)
  ;;; Possible other way to do that???
  ;ind = where(im ge 1)

  Q_2 = im lt 10000
  Q_1 = im gt 15000
  I_1 = where( Q_1 )
  I_2 = where( Q_2 )

  ; data selection mask for correlation
  mask = im * 0. & mask[I_1] = 1
  corr_fring = where(convolve( mask, psf_str.cpsf ) ge 1.e-3 * max(cpsf))
  corr_g = where((histogram(corr_fring, min = 0, max=size(/n_e,im)) + $
    histogram(I_2, min = 0, max=size(/n_e,im))) eq 2)

  in  = {s      : ptr_new(I_1)      ,$
    g      : ptr_new(corr_g)     ,$
    im     : ptr_new(im)       ,$
    bg     : ptr_new(bkg)        ,$
    bg_bkup  : ptr_new(bg_bkup)    ,$
    x      : ptr_new(/allocate_heap)   ,$
    y      : ptr_new(/allocate_heap)   ,$
    c_exp  : ptr_new(/allocate_heap)   ,$
    c_stat : 0             ,$
    it   : 1             ,$
    lev    : 1             }

  dst_em_fft, in, psf_str.cpsf, it=it, level=level
  c = im *0. & c[I_1] = *in.x

  ;;; SATURATED PIXEL POSITION
  c1 = sigma_filter( c, 5, n_sigma = 1, /iterate)     ;;;Necessary to avoid artifacts in the
  cc_bg = bkg * 0 & cc_bg[I_1] = bkg[I_1]
  cc = (convolve(c1,psf_str.opsf)) + cc_bg        ;;;correlation map


  print, 'Maximum correlation map:', max(cc)
  q_s = max(cc) gt str.sat_lev ? cc ge str.sat_lev : cc ge 0.25*max(cc)
  ;s = max(cc) gt str.sat_lev ? where(cc ge str.sat_lev, ns) : where(cc ge 0.25*max(cc), ns)
  s = where( q_s, ns )
  help, s

  ;;; BLOOMED PIXELS
  ;  b = where(histogram(I_1 , min = 0, max=size(/n_e, im)) - $
  ;    histogram(s , min = 0, max=size(/n_e, im)) eq 1.)
  b = where( Q_1 and ~q_s, nbloom)

  ;s = where(cc ge 0.2*max(cc), ns)

  mask_s = im * 0. & mask_s[s] = 1.

  qfwd_ind = convolve( mask_s, psf_str.cpsf ) ge 1.e-3 * max(cpsf)
  fwd_ind = where( qfwd_ind )
  ;;; ...flaring intense part of the image that is unsaturated
  ;mask_z = im * 0. & mask_z[where( im gt 10000 )] = 1.


  ;;mask_z = im * 0. & mask_z[ I_2 ] = 1.
  mask_z = im * 0.+1 & mask_z[ I_2 ] = 0. ;;;; GT modification


  mask_z = convolve(mask_z, psf_str.opsf) > 0
  z = where(mask_z gt min(psf_str.opsf[where(psf_str.opsf gt 0)]))

  ;;; DATA (Fringe) PIXELS
  g = where(histogram(I_2, min = 0, max=size(/n_e, im)) + $
    histogram(fwd_ind, min = 0, max=size(/n_e, im)) - $
    histogram(z, min=0,max=size(/n_e, im)) eq 2.)
  mg = byte( im * 0 ) & mg[z] = 1b & gg = where( Q_2 and Qfwd_ind and ~mg )
  if ~same_data2( gg, g) then stop

  mask_g = im * 0. & mask_g[g] = 1.


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if loud eq 1 then begin
    tt = im*0 & tt[I_2] = 1
    index2map, info ,((im-bkg)*tt)^0.3, map
    plot_map, map, /positive, /square ,title = 'Core Masked Bkg-subtracted data', thick=1.5

    index2map, info ,cc^0.3, map
    plot_map, map, /positive, /square ,title = 'Correlation', thick=1.5
  endif
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  str.s = ptr_new(s)
  str.z = ptr_new(z)
  str.g = ptr_new(g)
  str.b = ptr_new(b)
  str.ns = ns

end
pro desat::update_history, info, records, _extra=extra
  update_history, info, 'desat alg', _extra = extra ;initial implementation of history record, ras, 21-feb-2015
  ;  update_history, index,  records,   mode=mode , debug=debug, $
  ;    caller=caller, routine=routine, noroutine=noroutine, version=version
  ;+
  ;    Name: update_history
  ;
  ;    Purpose: add history record(s) to input structure(s)
  ;
  ;    Input Parameters:
  ;      index   - structure vector or FITs header array
  ;      records - info to add - string/string array
  ;
  ;    Keyword Parameters:
  ;      routine - if set , routine name to prepend to each history record -
  ;                default is via: 'get_caller' unless /noroutine set
  ;      caller - synonym for routine
  ;      noroutine - if set, dont prepend 'Caller' to records
  ;      version   - if set, verion number, include 'VERSION:' string
  ;
  ;      mode - if set and #records=#index, add record(i)->index(i)
  ;             (default mode = record(*)->index(*) (all records->all structure)
end
function desat::psf, info, dwavelength, core_dim

  ; NAME:
  ; dst_psf_gen
  ; PURPOSE:
  ; returns the psf structure that contains the diffraction and dispersion and the complete
  ; aia psf.
  ; EXPLANATION:
  ; The dispersion component of the psf is considered cutting out a circular portion on the
  ; central component of the psf computed by aia_calc_psf_mod.pro. The radius of this circle is
  ; fixed to 5 pixels but it can be modified by users.
  ; CALLING SEQUENCE:
  ; psf = dst_psf_gen( info , dwavelength , core_dim )
  ; INPUTS:
  ; info    = index structure for the aia data, necessary to retrieve infos about the psf
  ;       to generate
  ; dwavelength   = correction parameter on the wavelength of the psf (look aia_calc_psf_mod.pro)
  ; core_dim  = radius of the central core of the psf.
  ; OUTPUT:
  ; psf.cpsf = diffraction component of the psf
  ; psf.opsf = dispersion component of the psf
  ; psf.psf  = complete psf
  ; CALLS:
  ; CALLED BY:
  ;   DESATURATION
  ; CALLS:
  ;   aia_calc_psf_mod.pro
  ;

  npix = min([info.naxis1,info.naxis2])
  wavelength = strtrim( string(info[0].WAVELNTH),1)
  have_psf = 0
  if ptr_valid( Self.psf ) then begin
    psf = *Self.psf
    have_psf =  npix eq psf.npix && wavelength eq psf.wavelength && dwavelength eq psf.dwavelength
  endif
  if ~have_psf then begin
    default, core_dim, 5
    default, dwavelength, 0.0
    psf   = dst_psf_gen( wavelength, npix, dwavelength, core_dim )
    self.psf = ptr_new( psf )
  endif
  return, psf

end

;+
; :Description:
;    Describe the procedure.
;
; :Params:
;    wav
;    ts
;    te
;    path
;
; :Keywords:
;    path_save
;    sat_lev
;    lev
;    it
;    npix
;    core_dim
;    dwavelength
;    use_prep
;    save_fts
;    aec
;    loud
;    psplot
;    onewindow
; :History: 28-apr-2015, pass USE_PREP to savefts
; 2-may-2015, G Torre, fix xcen and ycen when cutting a data subset
; :Author: gabriele torre and richard schwartz
;-
function desat::desaturation, wav, ts, te, path, path_save=path_save, sat_lev=sat_lev, lev=lev, it=it, npix=npix,$
  core_dim=core_dim ,$
  dwavelength = dwavelength, $
  use_prep = use_prep, $
  save_fts=save_fts, aec = aec, loud=loud, psplot = psplot, onewindow = onewindow

  ; NAME:
  ;	desaturation
  ; PURPOSE:
  ;	Desaturation routine for saturated SDO/AIA images
  ; EXPLANATION:
  ;	Recover the pixel intensity inside the primary saturation region by means a EM
  ;	based routine for the saturated images inside the user defined time interval.
  ;	Moreover, the fringes generated by the recovered components are subtracted and
  ;	bloomed intensity are substituted by the background intensities.
  ; CALLING SEQUENCE:
  ;	result = obj -> desaturation( wav, ts, te, path, sat_lev=sat_lev, lev=lev, it=it, core_dim=core_dim ,$
  ;					  save_fts=save_fts, aec = aec, DWAVELENGTH = dwavelength)
  ; INPUTS:
  ;	file	= string array containing the file names for the event taken into account
  ;	ts		= user defined analysis starting time
  ;	te		= user defined analysis end time
  ;	wav 	= string array containing the wavelength to be process
  ; OUTPUT:
  ;	** Structure RESULTS, 2 tags, length=8, data length=8:
  ;   DATA            POINTER   <PtrHeapVar4635>
  ;   INFO            POINTER   <PtrHeapVar4636>
  ; OPTIONAL:
  ;	sat_level	= saturation intensity level (default = 16300)
  ;	lev			= EM stopping rule parameter (default = 0.01)
  ;	it			= maximum iteration number (default = 300)
  ;	core_dim	= radius for the PSF core (default = 5)
  ;	npix		= setting the resulting image to be npix*npix (default = 499)
  ; KEYWORDS:
  ; DWAVELENGTH - dwavelength is the fractional change in the wavelength over the nominal wavelength
  ;               so for 131, 131.5 dwavelength would be 0.5/131 = 0.00381679
  ; USE_PREP - default 1, if set then use aia_prep to register the data and desat cutout outputs
  ; path_save	= path to folder where .fts file will be saved
  ;	save_fts	= if set save fts file at the end of the run (default = 1)
  ;	aec 		= if set will desaturate also sat. images with short exp. time (default = 1)
  ;	loud		= if 1 intermediate and final results are shown in windows (default = 1)
  ;	ONEWINDOW - if set and LOUD is set then plots are written and then overwritten in the same graphic window
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

  ;;;DESAT STRUCTURE DEFINITION
  default, save_fts	, 1
  default, loud 	 	, 1
  default, lev	 	, 0.01
  default, it 	 	, 300
  default, npix 	 	, 499
  default, core_dim	, 5
  default, sat_lev 	, 15000
  default, aec	 	, 1
  default, path_save	, './'
  default, bkg_method, 'quadratic' ;requires two observations prior to ts
  default, psplot, 0
  default, onewindow, 0
  default, dwavelength, 0.0
  results 	 = replicate({dst_result},size(/n_e, wav))
  info_results = PTR_NEW(/ALLOCATE_HEAP)

  for I_WAV=0, size(/n_e, wav)-1 do begin
    range = anytim( [ts,te] ) + [-180.0, 180]
    file = self->find_file(path, wav[I_WAV], range = range)
    str = {dst_str}

    self -> readfts, file

    str.file = ptr_new(file)
    str.ts	 = ts
    str.te	 = te
    str.bkg_method = bkg_method ;need this for time selection of background

    if exist(sat_lev)  then str.sat_lev	= sat_lev
    if exist(lev)	   then str.lev		= lev
    if exist(npix)	   then str.npix	= npix
    if exist(it)	   then str.it		= it

    self -> indices_analysis, str, aec = aec
    self -> readfts, file
    self->_bld_filenames, *str.time_ind , path_save ;Prep the filenames. Use them as needed. Write out unused Desat filenames at the end
    self -> ctrl_data, str.npix ;set all the maps to npix x npix pixels
    if keyword_set(save_fts) then begin

      self -> savefts, /original, use_prep = use_prep
    endif
    self -> data_restore_2, str, data, bg, info, dwavelength = dwavelength
    ;;; POINT SPREAD FUNCTION GENERATION
    psf = Self->psf( info, dwavelength, core_dim )

    n_rec = size(/n_e,*str.sat_ind)

    for i=0, n_rec-1 do begin

      print, 'Desaturation ' + string(info[i].wavelnth) + ':' + string(i+1) + ' of ' + $
        string(n_rec), ' Time:', anytim(/atim,info[i].date_obs)

      str.im  	= ptr_new(data[*, *, i])
      str.bg  	= ptr_new(convolve(bg[*,*,i],psf.opsf)>0)
      str.bg_bkup = ptr_new(bg[*,*,i])

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      if loud eq 1 then begin
        numwin = onewindow ? 31 : i
        title= 'AIA Desaturation Routine '
        ;text = info[i].TELESCOP+' '+info[i].INSTRUME+' '+strtrim(info[i].WAVELNTH,1)+' '+anytim(info[i].DATE_OBS, /vms)+' UT'
        text = info[i].INSTRUME+' '+strtrim(info[i].WAVELNTH,1)+' '+anytim(info[i].DATE_OBS, /yoh)
        if ~psplot then window, numwin, xsize = 900, ysize = 600, title= title + text else $
          sps, /colo, /landscape
        !p.multi=[0,3,2]
        !p.charsize = 1.5
        index2map, info[i], *str.im^0.3, map
        plot_map, map, /positive, /square ,title = 'Original', thick=1.5
        ;if psplot then  ssw_legend, [title, text], /bottom, /left

        index2map, info[i], *str.bg^0.3, map
        plot_map, map, /positive, /square ,title = 'Background'
      endif
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      self->identify_saturation_regions, info[i], str, psf, loud				;;;FIND SATURATED PIXEL/DATA POSITION

      mask_im  = lonarr(n_elements(*str.im)) + 1
      ;tmp_g    = lindgen( n_elements(*str.im) )
      ; s    = Primary saturation pixels position
      ; z    = Not saturated pixel of the flare (pixels close to the s saturation region)
      ;        The use of these intensities is avoid to reduce the dispersion effects
      ;        from s to g not taken into account by the method
      ; g    = Diffraction fringes pixel position
      ; b    = Bloomed pixels position
      mask_im[ [*str.z, *str.z, *str.g ] ] = 0
      outer = where( mask_im )
      ;		tmp_g = where( histogram(indx,  min = 0, max = max(indx)) - $
      ;					   histogram(*str.s, min = 0, max = max(indx)) - $
      ;					   histogram(*str.z, min = 0, max = max(indx)) - $
      ;					   histogram(*str.g, min = 0, max = max(indx)) gt 0)

      ;mult_fact = total((*str.im)[tmp_g])/total((*str.bg)[tmp_g])
      ;get a scaling factor between the non-flaring regions of the flare image and the background image
      mult_fact = total((*str.im)[outer])/total((*str.bg)[outer])

      *str.bg *= mult_fact

      ;;; EM FOR SMALL SATURATED REGION
      Fourier = str.ns lt 1000 ? 0 : 1
      print,'*****	BEGIN EM	*****'
      case Fourier of
        0:begin
          print, '*** SQM expectation maximization ****'
          dst_sqm, str, psf.cpsf
          dst_em , str, level=level, it=it
        end
        1:begin
          print, '*** FFT expectation maximization ****'
          dst_em_FFT, str, psf.cpsf, level=level, it=it
        end
      endcase
      print,'*****	END EM		*****'

      ;;;CENTRAL PSF CORE CONVOLUTION
      dst_disp_factor_2, str, mult_fact, psf


      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      if loud eq 1 then begin
        mask_g = *str.im*0 & mask_g[*str.g] = *str.c_exp
        D = (2./n_elements(*str.g))*(*str.im * alog(f_div(*str.im, mask_g)) + mask_g - *str.im )

        index2map, info[i], d, map
        plot_map, map, /positive, /square ,title = 'Residuals', thick=1.5

        index2map, info[i] ,(*str.x>0)^0.3, map
        plot_map, map, /positive, /square ,title = 'Reconstruction', thick=1.5
        !p.charsize = 1

        !p.multi = 0
        !p.thick = 1

      endif
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      data_fill = (self->get(/data))
      data_fill[*,*,(*str.sat_ind)[i]] = *str.x

      index2map, (self->get(/index)), data_fill, map

      self -> set, index=(self->get(/index)), map=map
      self -> set, grid=30, /limb

      if keyword_set(save_fts) then begin
        ;self -> savefts, reform( info[i].date_obs ), /original
        self -> savefts, reform( info[i].date_obs ), info_str = str, original = 0
      endif

      print, i
    endfor

    results[i_wav].data = ptr_new((self->get(/data))[*,*,*str.time_ind])
    results[i_wav].info = ptr_new((self->get(/index))[*str.time_ind])

  endfor
  if psplot then begin & device, /close & x & endif
  ;Write the remaining desat fits files,
  ;self -> savefts, original = 0

  return, results

end

pro desat::ctrl_data, pix, data = data_, info = info_

  ; NAME:
  ;	ctrl_data
  ; PURPOSE:
  ;	Reduce the FOV of a given data array.
  ; EXPLANATION:
  ;	Reduce the dimension of the FOV of a given data array to [pixel,pixel].
  ; 	If data and info are given, the routine vill redice them, otherwise it will act on
  ;	the object data and info
  ; CALLING SEQUENCE:
  ;	ctrl_data, pix 										acting on the object
  ;	ctrl_data, pix, data = data_in, info = info_in	acting on data_in and info_in
  ; INPUTS:
  ;	pix = number of pixel for the new FOV
  ; OPTIONAL:
  ;	data_ = 3D data array on wich we what to act
  ;	info_ = index array structure for data_ array
  ; KEYWORDS:
  ; CALLS:
  ;	data_restore_2
  ; PROCEDURE:
  ;	uses the extract_slice function to extract the new FOV around the center of the
  ;	of the image. The new resized data are store in the previous object refreshing
  ;	the infos array.

  var = exist(data_)+exist(info_) eq 2 ? 1 : 0

  default, data, self -> get(/data)
  default, info, self -> get(/index)

  if var eq 1 then data = data_
  if var eq 1 then info = info_

  n_f = size(/n_e, info)

  tmp = max(data, pos);[(size(/dim, data))[0:1]]/2
  ij = array_indices(data, pos)
  img_cent =[avg(ij[0,*]) , avg(ij[1,*])]

  data_rescaled = fltarr(pix, pix, n_f)

  for i = 0, n_f - 1 do begin

    zi = i
    data_rescaled[*,*,i] = EXTRACT_SLICE( data, pix, pix, img_cent[0], img_cent[1], zi, 0, 0, 0 )

    xp=(mk_map_xp(info[i].xcen,info[i].cdelt1,info[i].naxis1,1))[img_cent[0]]
    yp=(mk_map_yp(info[i].ycen,info[i].cdelt2,1,info[i].naxis1))[img_cent[1]]

    info[i].xcen = xp
    info[i].ycen = yp

    info[i].naxis1 = pix
    info[i].naxis2 = pix

  endfor

  if var eq 0 then begin
    index2map, info, data_rescaled, map
    self-> set, index=info, map=map
  endif else begin
    data_ = data_rescaled
    info_ = info
  endelse

end


pro desat__define, void

  void={desat, $
    filenames: ptr_new(), $
    psf: ptr_new(), $
    inherits sdo}
  return

end


