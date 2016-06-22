function desat::find_file, path, wav, range = valid_range
  ;+
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
  ;-
  file_name	= file_search(concat_dir(path,'*.f*ts'), COUNT=ct)
  if ct eq 0 then file_name	= file_search(concat_dir(path,'*.fits'), COUNT=ct)

  mreadfits, file_name, info ;, data ;ras, remove data because you only want to read the headers here
  ind 	= where(info.WAVELNTH eq float(wav))
  info  = info[ ind ]
  file_name = file_name[ ind ]

  ;ind  = where_within( anytim( info.date_obs ), valid_range )
  ind = where(anytim( info.date_obs ) ge valid_range[0] and anytim( info.date_obs ) le valid_range[1])

  info  = info[ ind ]
  file_name = file_name[ ind ]

  return,file_name

end


pro desat::readfts, file
  ;+
  ; NAME:
  ;	desat::readfts
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
  ;-

  mreadfits, file, info, data
  index2map, info, data, map
  self-> set, index=info, map=map, /no_copy

end


pro desat::get_data_and_info, str, data, info
  ;+
  ; NAME:
  ;	desat::get_data_and_info
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
  ;	obj -> get_data_and_info, str, data, bg, info
  ; INPUTS:
  ;	str  = dst_str global structure for the desaturation routine.
  ; OUTPUT:
  ; data = data 3D corresponding to saturated images into the selected time interval
  ;	bg	 = computed background estimation obtained by dst_bkg_fourier_interp routine
  ;	info = index structure associated to data array
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;	CALLS:
  ;		DECONV
  ;		dst_bkg_fourier_interp
  ; PROCEDURE:
  ;-

  ;;; cut data array to be npix*npix pixel
  self -> ctrl_data, str.npix

  ;;; data/info array storage
  if size(/n_e, *str.sat_ind) eq 1 and max(*str.sat_ind) eq -1 then begin
    print, 'Not saturated frames within the selected time range'
  endif else begin
    info	= (self ->get(/index))[*str.sat_ind]
    data    = ((self ->get(/data))>0)[*,*,*str.sat_ind]
  endelse
end


function desat::get_background_estimation, str, info

  ;;; deconvolution step on un-saturated data
  bg_real = self -> deconv( *str.bg_ind, str, it = str.it, lev = str.lev, dwavelength=dwavelength,  /fill )

  ;;; background estimation routine
  bg_real_info = (self -> get(/index))[*str.bg_ind]

  bg = dst_bkg_fourier_interp(bg_real, bg_real_info, info, str)

  return, bg
end


pro desat::indices_analysis, str, aec = aec
  ;+
  ; NAME:
  ;	desat::indices_analysis
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
  ;-

  default, aec, 1

  info = self ->get(/index) ;done after this, we get info from reading all the headers, not the data

  ;;; file indices for the selected time interval
  q_time_range = anytim(info.date_obs) gt anytim(str.ts) and anytim(info.date_obs) lt anytim(str.te)

  ;;; saturated frame indices
  q_sat_frame = info.DATAMAX ge str.sat_lev

  ;;; short exposure time file indices

  q_short_exp = size(/n_e,info) gt 1 ? info.exptime lt 0.8*median(info.exptime) : info.exptime lt 0.7*max(info.exptime)

  ;;; indices for saturated frames in the selected time interval
  sat_ind = ~keyword_set(aec) ? where(q_time_range and q_sat_frame and ~q_short_exp) : where(q_time_range and q_sat_frame)

  ;;; background indices estimation over the whole dataset
  q_bg_dataset = total(q_short_exp) gt 0 ? ~q_sat_frame or q_short_exp : ~q_sat_frame
  bg_dataset = where(q_bg_dataset)

  ;;; background indices estimation in the select time interval
  bg_ind = where(q_time_range and q_bg_dataset)

  ;;; background indices for background estimation routine (bg_ind_t_ind-2 < bg_ind < bg_ind_t_ind+2)
  for i = 0 , 1 do begin

    prew_ind = where(bg_dataset lt min(bg_ind), ct)
    if ct gt 0 then bg_ind = [ max(bg_dataset[prew_ind]) , bg_ind ]

    post_ind = where(bg_dataset gt max(bg_ind), ct)
    if ct gt 0 then bg_ind = [ bg_ind , min(bg_dataset[post_ind]) ]

  endfor

  str.time_ind	= ptr_new(where(q_time_range))
  str.sat_ind 	= ptr_new(sat_ind)
  str.bg_ind	= ptr_new(bg_ind)
end


function desat::deconv, index, str, fill = fill, it = it, lev = lev, dwavelength = dwavelength
  ;+
  ; NAME:
  ;	desat::deconv
  ; PURPOSE:
  ;	Deconvolution procedure EM based for SDO/AIA images
  ; EXPLANATION:
  ;	Deconvolution were computed with the whole psf over a reduced FOV for SDO/AIA data
  ; CALLING SEQUENCE:
  ;	dst_result = obj -> desat( index, fill=fill )
  ; INPUTS:
  ;	index  = indices of data images in the obj on wich deconvolution has to be performed.
  ; OUTPUT:
  ;	flux = deconvolved images
  ; KEYWORDS:
  ;	fill = Deconvolved data were substituted into the obj in index position
  ;	it   = Maximum number of iterations
  ;	lev  = Stopping rule parameter for EM
  ; CALLS:
  ;	CALLED BY:
  ;		DESATURATION
  ;	CALLS:
  ;		dst_em_fft
  ;-

  default, it , 1000
  default, lev, 0.1

  info = (self-> get(/index))[0]

  imag = (self-> get(/data))[*,*,index]
  flux = imag * 0

  dim	= [[info.naxis1],[info.naxis2]]

  ;;; DECONVOLUTION FOR BACKGROUND ESTIMATION
  for i=0, n_elements(index)-1 do begin

    print, 'Deconvolution: ' + strtrim(i+1,1), ' of ', strtrim(size(/n_e, index),1)

    in = {	s: ptr_new(findgen(product(dim))),$
      g: ptr_new(findgen(product(dim))),$
      c_exp: ptr_new(/allocate_heap),$
      x: ptr_new(/allocate_heap),$
      y: ptr_new(/allocate_heap),$
      im: ptr_new(imag[*,*,i]),$
      bg_bkup: ptr_new(flux[*,*,i]),$
      bg: ptr_new(flux[*,*,i]),$
      ns: product(dim),$
      npix: str.npix,$
      c_stat: 0,$
      lev: lev,$
      it: it}

    if size((*self.psf).psf,/n_d) eq 3 then begin
      message, 'Not Enabled'
;      dst_em_fft_mw, in, (*self.psf).psf, it=in.it, level=level
;      flux[*,*,i] = reform( total(*in.x, 2), dim )
    endif else begin
      self->dst_em_fft, in, (*self.psf).psf, it=in.it, level=level
      flux[*,*,i] = reform( *in.x, dim )
    endelse

    imag[*,*,i] = convolve(flux[*,*,i], ((*self.psf).opsf)[*,*,0])>0
  endfor

  if keyword_set(fill) eq 1 then begin
    dec_info = self-> get(/index)
    dec		 = self-> get(/data)
    dec[*,*,index] = imag

    index2map, dec_info, dec, map
    self-> set, index=dec_info, map=map
  endif

  return, flux

end


pro desat::_bld_filenames, index, data_path, peaklam = peaklam
  ;+
  ; NAME:
  ; desat::_bld_filenames
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
  ;-

  n_el = size(/n_e,index)

  info = (self -> get(/index))[index]

  stringtime = anytim(info.date_obs,/vms,/date_only)
  stringwave = strtrim(string(peaklam),1)

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
  ;+
  ; NAME:
  ;	desat::savefts
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
  ;-

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
    info = { desat_flux: *info_str.sat_flux, $
      background: *info_str.bg, $
      sat_core: *info_str.s, $
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


pro desat::identify_saturation_regions, info, str, loud
  ;+
  ; NAME:
  ; desat::identify_saturation_regions
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
  ;-

  cpsf = (*self.psf).cpsf
  opsf = ((*self.psf).opsf)[*,*,0]

  im    = *str.im
  bkg   = *str.bg
  bg_bkup = *str.bg_bkup

  Q_2 = im lt 10000
  Q_1 = im gt 15000
  I_1 = where( Q_1 )
  I_2 = where( Q_2 )

  ; data selection mask for correlation
  mask = im * 0. & mask[I_1] = 1

  if size(cpsf, /n_d) eq 3 then begin
    corr_fring = im * 0
    for iw = 0, (size(cpsf, /dim))[2] - 1 do corr_fring += convolve( mask, cpsf[*,*,iw] )
    corr_fring = corr_fring ge (1.e-3 * max(cpsf))
  endif else begin
    corr_fring = convolve( mask, cpsf ) ge (1.e-3 * max(cpsf))
  endelse

  corr_g = where(corr_fring and Q_2)

  in  = {s      : ptr_new(I_1)				,$
    g      : ptr_new(corr_g)			,$
    im		: ptr_new(im)				,$
    bg		: ptr_new(bkg)				,$
    bg_bkup	: ptr_new(bg_bkup)			,$
    x		: ptr_new(/allocate_heap)	,$
    y		: ptr_new(/allocate_heap)	,$
    c_exp 	: ptr_new(/allocate_heap)	,$
    ns		: total(Q_1)				,$
    npix	: str.npix					,$
    c_stat	: 0							,$
    it   	: 1							,$
    lev		: 1							}

  if size(cpsf,/n_d) eq 3 then begin
    message, 'Not Enabled'
;    dst_em_fft_mw, in, cpsf, it=it, level=level
;    c = im *0. & c[I_1] = total(*in.x, 2)
  endif else begin
    self->dst_em_fft, in, cpsf, it=it, level=level
    c = im *0. & c[I_1] = *in.x
  endelse

  ;;; SATURATED PIXEL POSITION
  c1 = sigma_filter( c, 5, n_sigma = 1, /iterate)     ;;;Necessary to avoid artifacts in the
  cc_bg = bkg * 0 & cc_bg[I_1] = bkg[I_1]
  cc = (convolve(c1, opsf)) + cc_bg        ;;;correlation map

  print, 'Maximum correlation map:', max(cc)
  q_s = max(cc) gt str.sat_lev ? cc ge 14000 : cc ge 0.1*max(cc)
  s = where(q_s)

  ;;; BLOOMED PIXELS
  b = where( Q_1 and ~q_s, nbloom)

  mask_s = im * 0. & mask_s[s] = 1.

  if size(cpsf, /n_d) eq 3 then begin
    tmp = im * 0
    for iw = 0, (size(cpsf, /dim))[2] - 1 do tmp += convolve( mask_s, cpsf[*,*,iw] )
    qfwd_ind = tmp ge (1.e-3 * max(cpsf))
  endif else begin
    qfwd_ind = convolve( mask_s, cpsf ) ge 1.e-3 * max(cpsf)
  endelse
  fwd_ind = where( qfwd_ind )

  mask_z = im * 0 + 1 & mask_z[ I_2 ] = 0.
  q_z = convolve(mask_z, opsf) gt min(opsf[where(opsf gt 0)])
  z= where(q_z)

  ;;; DATA (Fringe) PIXELS
  g = where(Q_2 and qfwd_ind and ~q_z)
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
  str.ns = total(q_s)

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
  ;-
end


function desat::psf, info, dwavelength, core_dim
  ;+
  ; NAME:
  ; desat::psf
  ; PURPOSE:
  ; returns the psf structure that contains the diffraction and dispersion and the complete
  ; aia psf.
  ; EXPLANATION:
  ; The dispersion component of the psf is considered cutting out a circular portion on the
  ; central component of the psf computed by aia_calc_psf_mod.pro. The radius of this circle is
  ; fixed to 5 pixels but it can be modified by users.
  ; CALLING SEQUENCE:
  ; psf = self->psf( info , dwavelength , core_dim )
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
  ;-

  npix = min([info.naxis1,info.naxis2])
  wavelength = strtrim( string(info[0].WAVELNTH),1)
  have_psf = 0

  if ptr_valid( Self.psf ) then begin
    psf = *Self.psf
    have_psf = npix eq psf.npix && wavelength eq psf.wavelength && ARRAY_EQUAL(dwavelength, psf.dwavelength) eq 1;dwavelength eq psf.dwavelength
  endif

  if ~have_psf then begin
    default, core_dim, 5
    default, dwavelength, 0.0
    psf = self->dst_psf_gen( wavelength, npix, dwavelength, core_dim )
    self.psf = ptr_new( psf )
  endif


  return, psf

end


function desat::dst_psf_gen, wavelength, npix, dwavelength, core_dim
  ;+
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
  ;-

  default, core_dim, 5

  psf		= fltarr(npix, npix, n_elements(dwavelength))
  cpsf	= fltarr(npix, npix, n_elements(dwavelength))
  opsf	= fltarr(npix, npix, n_elements(dwavelength))

  ; build a circular mask
  xgrid	= (fltarr(npix)+1)##indgen(npix)
  ygrid	= indgen(npix)##(fltarr(npix)+1)
  center	= [fix(npix/2.),fix(npix/2.)]
  w		= where((xgrid-center[0])^2+(ygrid-center[1])^2 le core_dim^2)
  mask_c	= fltarr(npix,npix) & mask_c[w] = 1

  for iw = 0 , n_elements(dwavelength)-1 do begin
    print, iw
    wav = (dwavelength[iw] * float(wavelength)) + float(wavelength)
    file = file_search(concat_dir('AIA_PSF' ,'AIA_PSF_'+ strtrim(wav,1) +'_'+strtrim(npix,1)+'.fts'), count=ct)

    if ct gt 0 then begin
      mreadfits, file, info, psf_IN
    endif else begin
      FILE_MKDIR , 'AIA_PSF'
      psf_in = aia_calc_psf_mod( wavelength, npix = npix, dwavelength = dwavelength[iw])
      m1		= fltarr(npix,npix) & m1[npix/2 , npix/2 ] = 1
      psf_in  = convolve(/corr, m1, psf_in) > 0.0
      mwrfits, psf_in , concat_dir('AIA_PSF' ,'AIA_PSF_'+ strtrim(wav,1) +'_'+strtrim(npix,1)+'.fts')
    endelse

    psf[*,*,iw]  = psf_in
    opsf[*,*,iw] = REFORM(psf[*,*,iw] * mask_c)
    cpsf[*,*,iw] = REFORM(psf[*,*,iw] * (1 - mask_c))

  endfor
  psf		= {wavelength: wavelength, dwavelength: dwavelength, npix: npix, cpsf:cpsf, opsf:opsf, psf:psf}
  return, psf
end


pro desat::dst_em_fft , str , cpsf , it=it , level=level
  ;+
  ; NAME:
  ;	dst_em_fft
  ; PURPOSE:
  ;	Deconvolution/desaturation routine based on EM method implemented by using FFT convolution.
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
  ;-

  default , it , str.it
  default , level , str.lev
  default , pow , 1

  s   = *str.s
  g	= *str.g
  psf = cpsf
  ;;;ZEROPADDING THE DATA ARRAY TO AN N*N ARRAY WHERE N IS A POWER OF 2 (INCREASES FFT SPEED)
  dim = size(/dim,*str.im)
  check = alog(dim)/alog(2)
  if round(check[0]) ne check[0] or round(check[1]) ne check[1] then dim = [2^(floor(check[0])+1),2^(floor(check[1])+1)]

  im 		  = fltarr(dim)
  bkg 	  = fltarr(dim)
  mask_g  = fltarr(dim)
  x		    = fltarr(dim)
  ;;;

  emp_kkt = fltarr(it+1)
  exp_kkt = fltarr(it+1)

  mask = im * 0 & mask[dim[0]/2,dim[1]/2] = 1

  im[0,0]  = *str.im
  bkg[0,0] = *str.bg

  psf  = convolve(mask , psf)>0
  psf2 = psf * psf

  x_tmp 	 = *str.im*0.
  x_tmp[s] = total(*str.bg_bkup) gt 0 ? (*str.bg_bkup)[s] : 1

  x[0,0] 	 = x_tmp

  sat_pos  = where(f_div(x , x) gt 0)

  mask_g_4_zeropad = *str.im * 0. & mask_g_4_zeropad[g] = 1.
  mask_g[0,0] = mask_g_4_zeropad

  y 	= ( im  * mask_g )>0.
  bkg = ( bkg * mask_g )>0.				 ;;; GT modification

  V = convolve(mask_g , psf , /corr )>0.

  for I = 0, it-1 do begin

    y_i = (convolve( x , psf ) + bkg)>0

    Z   = (f_div(y , y_i)) * mask_g

    U = convolve( z , psf , /corr )>0.
    x[sat_pos] *= ( f_div( U , V ) )[sat_pos]

    emp_kkt[i] = total( ( x * (V - U))^2.)
    exp_kkt[i] = total( x^2. * (convolve( f_div(1. , y_i * mask_g)>0 , psf2 , /corr )>0.) )

    if emp_kkt[i] le level * exp_kkt[i] then break

    print, i ,' --- ', emp_kkt[i]/exp_kkt[i] , level

  endfor

  dim = size(/dim , *str.im)

  x   = x[0:dim[0]-1 , 0:dim[1]-1]
  y   = y[0:dim[0]-1 , 0:dim[1]-1]
  y_i = y_i[0:dim[0]-1 , 0:dim[1]-1]

  C_stat = c_statistic( y , y_i )

  print , i , '  R_flx:' , total( y_i[g]) , '  Exp_flx:' , total( y[g] )  , '  C_stat:' , c_stat , ' '

  str.x = ptr_new(x[s])
  str.y = ptr_new(y)
  str.c_exp  = ptr_new(y_i[g])
  str.c_stat = c_stat

end
pro desat::dst_sqm , str , cpsf
  ; NAME:
  ; dst_sqm
  ; PURPOSE:
  ; It computes the sqm matrix for dst_em routine
  ; EXPLANATION:
  ; select the components of the convolution matrix that links the saturated pixels str.s
  ; with the diffraction fringes in str.g
  ; CALLING SEQUENCE:
  ; dst_sqm , str , cpsf
  ; INPUTS:
  ; str = dst_str global structure for the desaturation routine.  \
  ; cps = diffraction component of the PSF
  ; OUTPUT:
  ; str.sqm = sqm matrix for dst_em routine
  ; str.y = vector containing the diffraction fringes intensities
  ; str.bg  = vector containing the estimated background intensities for the same pixels
  ;       considered in y
  ; CALLS:
  ; CALLED BY:
  ;   DESATURATION
  ;

  s     = *str.s
  g     = *str.g
  n_sat = str.ns

  dim   = size(/dim , *str.im)
  dim_psf = size(/dim , cpsf)
  n_g   = size(/n_e , g)

  sqm = fltarr(n_g , n_sat)

  ij_sat = array_indices(*str.im , s)

  im1 = *str.im * 0.0 & im1[dim[0]/2, dim[1]/2] = 1 & c0 = convolve(/corr,im1, cpsf) > 0.0
  for i=0, n_sat-1 do sqm[*,i] = (shift(c0 , ij_sat[0,i] - dim[0]/2 , ij_sat[1,i] - dim[1]/2 ))[g]

  str.sqm = ptr_new(sqm)
  str.y = ptr_new((*str.im)[g] > 0.)
  str.bg  = ptr_new((*str.bg)[g] > 0.)

end



pro desat::dst_em, str ,level=level, it=it
  ;+
  ; NAME:
  ;	dst_em
  ; PURPOSE:
  ;	Computes the intensities of saturated pixels by using an EM algorithm.
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

  C_stat = c_statistic( y , y_i )

  print , i , '  R_flx:' , total( y_i ) , '  Exp_flx:' , total( y )  , '  C_stat:' , c_stat , ' '

  str.x = ptr_new(x)
  str.c_exp = ptr_new(y_i)
  str.c_stat = C_stat

end


pro desat::dst_image_synthesis, str , mult_fact , psf
;+
; NAME:
;	dst_image_synthesis
; PURPOSE:
;	compute the diffraction component to subtract from the image, the convolution between 
;	the retrieved intensties, inside the saturation region, and the central core of the 
;	PSF and the intensities for the bloomed pixels.
; EXPLANATION:
;	1) Computes the convolution product between the intensities obtained from the EM method 	
;		and the central core of the psf. 
;	2) The diffraction fringes generated by the retrieved intensities are computed and 
;		subtracted from the original image
;	3) The bloomed intensities are substituted by the estimated background values.
; CALLING SEQUENCE:
;	dst_image_synthesis, str , mult_fact , psf
; INPUTS:
;	str		  =	dst_str global structure for the desaturation routine.	
;	mult_fact = multiplicative factor obtained in desaturation routine necessary to 
;				renormalize the estimated background, minimizing the interpolation error 
;				made during the interpolation of the zero-th frequency component in 
;				dst_bkg_fourier_interp.pro
;	psf 	  = structure containing all the psf components, diffraction, dispersion and 
;				the complete one. This structure is returned as output from dst_psf_gen.pro
;				routine.	
;	OUTPUT:
;	It fills the x tag of str
; CALLS:
;	CALLED BY:
;		DESATURATION
;	CALLS:
;		CONVOLVE
; PROCEDURE:
;-

	im = *str.im > 0
	s  = *str.s 
	b  = *str.b

	opsf = psf.opsf
	cpsf = psf.cpsf
		
	mwav = size(/n_dim, *str.x) - 1

	nwav = mwav ? (size(/dim, *str.x))[1] : 1
	fin_res = mwav ? fltarr(str.npix,str.npix,nwav) : fltarr(str.npix,str.npix)

	x_im = mwav ? fltarr(str.npix,str.npix,nwav) : fltarr(str.npix,str.npix)
	f_im = *str.im
	
	check = where(b gt 0 , ct)
	fq = intarr(long(str.npix)*long(str.npix))+1
	fq[where(im gt 15000)] = 0
	;if ct gt 1 then fq[[s,b,*str.z]] = 0 else fq[[s,*str.z]] = 0 
	f = where(fq eq 1)
	ijs = array_indices(im , s)
	
	for iw = 0, nwav - 1 do begin 

		xtmp = fltarr(str.npix,str.npix)
		b_im = fltarr(str.npix,str.npix)
		
		xtmp[s] = (*str.x)[*,iw]
		
		f_im[f] -= (convolve(xtmp, cpsf[*,*,iw])>0)[f]
		
		if ct gt 0 then b_im[b] = mult_fact*(convolve(*str.bg_bkup,opsf[*,*,iw])>0)[b]
		
		xtmp[s] = (convolve(xtmp, opsf[*,*,iw])>0)[s]
		
		x_im[*,*,iw] = xtmp + b_im
		
	endfor
	
	for iw = 0, nwav - 1 do x_im[*,*,iw] += f_im * reform(fq, str.npix,str.npix)
	
	str.x = ptr_new(x_im)	
	
end



function desat::desaturation, wav, ts, te, path, path_save=path_save, sat_lev=sat_lev, lev=lev, it=it, npix=npix,$
  core_dim=core_dim , dwavelength = dwavelength, wavstrategy = wavstrategy, peaklam = peaklam, $
  use_prep = use_prep, save_fts=save_fts, aec = aec, loud=loud, psplot = psplot, onewindow = onewindow, $
  bkg_filename = bkg_filename
  ;+
  ; NAME:
  ;	desat::desaturation
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
  ; peaklam	    = set the wavelength value to use for the generation of the PSF if wavstrategy = 1
  ;	wavstrategy = wavelength strategy definition
  ;					 0: generate the PSF using the nominal wavelength value of the passband
  ;					 1: generate the PSF using the wavelength value defined by peaklam
  ;					 2: generate the PSF using the wavelength associated to the BRIGHTEST EMISSIVITY
  ;						value computed by flare_peak_wavelength.pro
  ;	bkg_filename= filename of the image to use as background
  ; KEYWORDS:
  ; DWAVELENGTH - dwavelength is the fractional change in the wavelength over the nominal wavelength
  ;               so for 131, 131.5 dwavelength would be 0.5/131 = 0.00381679
  ; PASSBANDWAV -
  ; USE_PREP - default 1, if set then use aia_prep to register the data and desat cutout outputs
  ; path_save	= path to folder where .fts file will be saved
  ;	save_fts	= if set save fts file at the end of the run (default = 1)
  ;	aec 		= if set will desaturate also sat. images with short exp. time (default = 1)
  ;	loud		= if 1 intermediate and final results are shown in windows (default = 1)
  ;	ONEWINDOW - if set and LOUD is set then plots are written and then overwritten in the same graphic window\
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
  ;		dst_image_synthesis
  ; :History: 28-apr-2015, pass USE_PREP to savefts
  ;  2-may-2015, G Torre, fix xcen and ycen when cutting a data subset
  ; :Author: gabriele torre and richard schwartz
  ;-
  ;

  ;;;DESAT STRUCTURE DEFINITION
  default, save_fts, 	1
  default, loud,		1
  default, lev,			1
  default, it,			300
  default, npix,		499
  default, core_dim,	5
  default, sat_lev,		15000
  default, aec,			1
  default, path_save,	'./'
  default, bkg_method,	'quadratic' ;requires two observations prior to ts
  default, psplot,		0
  default, onewindow,	0
  default, wavstrategy, 0
  default, dwavelength, 0.0
  default, peaklam,		float(wav)
  default, bkg_filename,''

  results 	 = replicate({dst_result},size(/n_e, wav))
  info_results = PTR_NEW(/ALLOCATE_HEAP)

  for i_wav = 0, size(/n_e, wav)-1 do begin

    range = anytim( [ts,te] ) + [-180.0, 180]
    file = self->find_file(path, wav[i_wav], range = range)

    str = {dst_str}

    self -> readfts, file

    str.file = ptr_new(file)
    str.ts	 = ts
    str.te	 = te
    str.bkg_method = bkg_method ;need this for time selection of background

    if exist(sat_lev) 		then str.sat_lev	= sat_lev
    if exist(lev)	  		then str.lev		= lev
    if exist(npix)	  		then str.npix		= npix
    if exist(it)	  		then str.it			= it
    if exist(bkg_filename)	then str.bkg_filename= bkg_filename

    self -> indices_analysis, str, aec = aec
    self -> readfts, file

	;;;desaturation wavelength strategy definition    
    if wavstrategy eq 1 and ~exist(peaklam) then begin
      MESSAGE, 'Error: peaklam not defined for strategy ' +strtrim(wavstrategy,1)
    endif
    if wavstrategy eq 2 then begin
      flare_peak_wavelength, float(wav[I_WAV]), peaklam = peaklam;, /plotall
    endif
	if wavstrategy ne 0 then begin
      print, "Desaturation performed at lambda = " + string(peaklam) +" Angstrom"
      dwavelength =fltarr(size(/n_e,peaklam))
      for pw = 0 , size(/n_e,peaklam)-1 do dwavelength[pw] = (peaklam[pw] - float(wav))/float(wav)
    endif

    self ->_bld_filenames, *str.time_ind , path_save, peaklam = peaklam

    self -> ctrl_data, str.npix 					 ;set all the maps to npix x npix pixels

    if keyword_set(save_fts) then self -> savefts, /original, use_prep = use_prep

    self -> get_data_and_info, str, data, info

    psf_ptr = self -> psf( info, dwavelength, core_dim )

    if strmatch(str.bkg_filename, '') eq 0 then begin
      mreadfits,concat_dir(path,str.bkg_filename), info_, bg_
      bg_   = [[[bg_]],[[bg_]]]
      info_ = [[info_],[info_]]
      self -> ctrl_data, str.npix, data = bg_, info = info_
      bg = bg_[*,*,0]
    endif else begin
      bg = self -> get_background_estimation(str, info)
    endelse

    n_rec = size(/n_e,*str.sat_ind)

    for i=0, n_rec-1 do begin

      print, 'Desaturation ' + string(info[i].wavelnth) + ':' + string(i+1) + ' of ' + $
        string(n_rec), ' Time:', anytim(/atim,info[i].date_obs)

      str.im = ptr_new(data[*, *, i])
      if strmatch(str.bkg_filename, '') eq 0 then begin
        str.bg = ptr_new(bg)
        str.bg_bkup = ptr_new(bg)
      endif else begin
        str.bg = ptr_new(convolve(bg[*,*,i],((*self.psf).opsf)[*,*,0])>0)
        str.bg_bkup = ptr_new(bg[*,*,i])
      endelse

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      if loud eq 1 then begin
        numwin = onewindow ? 31 : i
        title= 'AIA Desaturation Routine '
        text = info[i].INSTRUME+' '+strtrim(info[i].WAVELNTH,1)+' '+anytim(info[i].DATE_OBS, /yoh)
        if ~psplot then window, numwin, xsize = 900, ysize = 600, title= title + text else $
          sps, /colo, /landscape
        !p.multi=[0,3,2]
        !p.charsize = 1.5
        index2map, info[i], *str.im^0.3, map
        plot_map, map, /positive, /square ,title = 'Original', thick=1.5

        index2map, info[i], *str.bg^0.3, map
        plot_map, map, /positive, /square ,title = 'Background'
      endif

      ;psf_ptr = self -> psf( info, 0, core_dim )
      self->identify_saturation_regions, info[i], str, loud				;;;FIND SATURATED PIXEL/DATA POSITION
      ;psf_ptr = self -> psf( info, dwavelength, core_dim )

      mask_im  = lonarr(n_elements(*str.im)) + 1
      mask_im[ [*str.z, *str.z, *str.g ] ] = 0
      outer = where( mask_im )

      mult_fact = total((*str.im)[outer])/total((*str.bg)[outer])
      *str.bg *= mult_fact
	  
      ;;; EM FOR SMALL SATURATED REGION
      Fourier = str.ns lt 1000 ? 0 : 1
      print,'*****	BEGIN EM	*****'
      if size((*self.psf).cpsf,/n_d) eq 2 then begin
        case Fourier of
          0:begin
            print, '*** SQM expectation maximization ****'
            self->dst_sqm, str, (*self.psf).cpsf
            self->dst_em , str, level=level, it=it
          end
          1:begin
            print, '*** FFT expectation maximization ****'
            self->dst_em_FFT, str, (*self.psf).cpsf, level=level, it=it
          end
        endcase
      endif
      print,'*****	END EM		*****'

      str.sat_flux = str.x

      ;;;CENTRAL PSF CORE CONVOLUTION
      self->dst_image_synthesis, str, mult_fact, *self.psf

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;PLOT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      if loud eq 1 and size((*self.psf).cpsf,/n_d) eq 2 then begin

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
        self -> savefts, reform( info[i].date_obs ), info_str = str, original = 0
      endif
      
    endfor

    results[i_wav].data = ptr_new((self->get(/data))[*,*,*str.time_ind])
    results[i_wav].info = ptr_new((self->get(/index))[*str.time_ind])

  endfor
  if psplot then begin & device, /close & x & endif

  return, results

end


pro desat::ctrl_data, pix, data = data_, info = info_
  ;+
  ; NAME:
  ;	desat::ctrl_data
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
  ;-
  var = exist(data_)+exist(info_) eq 2 ? 1 : 0

  default, data, self -> get(/data)
  default, info, self -> get(/index)

  if var eq 1 then data = data_
  if var eq 1 then info = info_

  n_f = size(/n_e, info)

  img_cent = [(size(/dim, data))[0:1]]/2

  data_rescaled = n_f eq 1 ? fltarr(pix, pix) : fltarr(pix, pix, n_f)

  for i = 0, n_f - 1 do begin

    zi = i
    if n_f eq 1 then begin
      data_rescaled[*,*] = EXTRACT_SLICE( [[[data]],[[data]]], pix, pix, img_cent[0], img_cent[1], zi, 0, 0, 0 )
    endif else begin
      data_rescaled[*,*,i] = EXTRACT_SLICE( data, pix, pix, img_cent[0], img_cent[1], zi, 0, 0, 0 )
    endelse

    xp=(mk_map_xp(info[i].xcen,info[i].cdelt1,info[i].naxis1,1))[img_cent[0]]
    yp=(mk_map_yp(info[i].ycen,info[i].cdelt2,1,info[i].naxis1))[img_cent[1]]

    info[i].xcen = xp
    info[i].ycen = yp

    info[i].naxis1 = pix
    info[i].naxis2 = pix

    info[i].crpix1 = info[i].crpix1 - img_cent[0] + (pix/2.)
    info[i].crpix2 = info[i].crpix2 - img_cent[1] + (pix/2.)

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


