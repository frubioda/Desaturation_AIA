pro dst_str__define
  ;+
  ; NAME:
  ;	dst_str__define
  ; PURPOSE:
  ;	Definition of the global structure for the desaturation procedure
  ; CALLING SEQUENCE:
  ;	str = {dst_str}
  ;	OUTPUT:
  ;	s		 = Primary saturation pixels position
  ;	z		 = unsaturated pixel of the flare (pixels close to the s saturation region)
  ;				 The use of these intensities is avoid to reduce the dispersion effects
  ;				 from s to g not taken into account by the method
  ;	g		 = Diffraction fringes pixel position
  ;	b		 = Bloomed pixels position
  ;	bg		 = Estimated background in the image space
  ;	bg_bkup	 = Estimated background in the object space
  ;	y		 = Data vector for EM method
  ;	sqm		 = Convolution matrix used for EM method (y = sqm#x + bg)
  ;	x		 = Vector of primary saturated intensities
  ;	c_exp	 = Expectation vector computed at the end of EM
  ;	file 	 = File names string array on wich the analysis is performed
  ;	bg_ind   = Temporal indices for images used into the background estimation procedure
  ;	time_ind = Temporal indices for frames within the time range selected by the user
  ;	sat_ind  = Temporal indices for saturated images into the dataset
  ;	ts 		 = Analysis starting time selected by the user
  ;	te 		 = Analysis end time selected by the user
  ;	sat_lev	 = Saturation level (default = 16300)
  ;	c_stat 	 = C_statistics computed between data (y) and the computed expectation (c_exp)
  ;	lev		 = EM stopping rule level (default = 0.01)
  ;	it		 = Maximum number of EM iterations
  ;	ns		 = Number of saturated pixels
  ;	npix	 = Number of pixel for the result images (npix*npix)
  ;-
  str = {	dst_str							,$
    im		:PTR_NEW(/ALLOCATE_HEAP),$
    s		:PTR_NEW(/ALLOCATE_HEAP),$
    z		:PTR_NEW(/ALLOCATE_HEAP),$
    g		:PTR_NEW(/ALLOCATE_HEAP),$
    b		:PTR_NEW(/ALLOCATE_HEAP),$
    bg		:PTR_NEW(/ALLOCATE_HEAP),$
    bg_bkup	:PTR_NEW(/ALLOCATE_HEAP),$
    y		:PTR_NEW(/ALLOCATE_HEAP),$
    sqm		:PTR_NEW(/ALLOCATE_HEAP),$
    x		:PTR_NEW(/ALLOCATE_HEAP),$
    c_exp	:PTR_NEW(/ALLOCATE_HEAP),$
    file 	:PTR_NEW(/ALLOCATE_HEAP),$
    bg_ind  :PTR_NEW(/ALLOCATE_HEAP),$
    time_ind:PTR_NEW(/ALLOCATE_HEAP),$
    sat_ind :PTR_NEW(/ALLOCATE_HEAP),$
    sat_flux:PTR_NEW(/ALLOCATE_HEAP),$
    ts 		:''						,$
    te 		:''						,$
    sat_lev	:0.						,$
    c_stat 	:0.						,$
    lev		:0.						,$
    it		:0	 					,$
    ns		:0						,$
    npix	:0            			,$
    bkg_method:''					,$
    bkg_filename:''					}
end