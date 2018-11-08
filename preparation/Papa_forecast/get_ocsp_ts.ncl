; get_ocsp_ts
;==========================================================================
;
; USAGE:
;  ncl get_ocsp_ts.ncl
;
; DESCRIPTION:
;  Function to retrieve latest subsurface temperature and salinity data from
;  Ocean Climate Station Papa mooring and export as ASCII files
;
; INPUT:
;
; OUTPUT:
;  ASCII files containing subsurface temperature and salinity data
;
; AUTHOR:
;  September 29 2018. Zhihua Zheng                     [ zhihua@uw.edu ]
;==========================================================================

begin

  ;--- LOAD netCDF file served by OPeNDAP server at NDBC -----------------------
  url = "https://dods.ndbc.noaa.gov/thredds/dodsC/oceansites/DATA_GRIDDED/PAPA/"
  f_name = "OS_PAPA_200706_M_TSVMBP_50N145W_hr.nc"

  exists = isfilepresent(url+f_name)
    if(.not.exists) then
      print("OPeNDAP test unsuccessful.")
      print("Either the file doesn't exist, or NCL does")
      print("not have OPeNDAP cabilities on this system.")
    else
      f = addfile(url + f_name,"r")
      ;vnames = getfilevarnames(f)
      ;print(vnames)
    end if

  ;---- READ VARIABLES ---------------------------------------------------------

  time = f->TIME
  t_length = dimsizes(time)
  t = time(t_length-1) ; only grab most recent data

  ; temperature profile
  depth_t = f->DEPTH
  tprof = f->TEMP(t_length-1,:,0,0)

  ; salinity profile
  depth_s = f->DEPPSAL
  sprof = f->PSAL(t_length-1,:,0,0)

  ; get the date string for mooring data
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
  mooring_date = cd_string(t,"%Y-%N-%D %H:%M:%S")

  ; get date for yesterday
  yest_date = systemfunc("date -v-1d +%F")
  format = "%Y-%N-%D"

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"
  yest_t = cd_inv_string(yest_date,format)
  format@units = t@units

  ; if the mooring data is delayed by more than one day, then use the mooring
  ; data as from yesterday, hence the GFS forecast data can update as usual.
  if ((yest_t-t) .gt. 0) then
    date = cd_string(yest_t,"%Y-%N-%D %H:%M:%S")
    ; mooring data usually delays one day
  else
    date = mooring_date
  end if

  ;---- EXPORTATION TO ASCII FILES ----------------------------------------------

  ; linear interpolate the missing values
  tprof@_FillValue = 1e+35
  t_prof = linmsg(tprof,-1)  ; endPoint= -1 missing end point(s)
                             ;              set to nearest non-missing value
  sprof@_FillValue = 1e+35
  s_prof = linmsg(sprof,-1)  ; endPoint= -1 missing end point(s)
                             ;              set to nearest non-missing value

  ; the directory storing forcing data
  basecase = "../../data/OCSPapa_forecast/"

  file_I = basecase+"/tprof_file.dat"
  system("/bin/rm -f "+file_I)   ; remove any pre-existing file
  write_table(file_I, "w", [/date,dimsizes(t_prof),2/], "%s  %i  %i")
  write_table(file_I, "a", [/-depth_t,t_prof/], "%6.1f   %9.6f") ; append

  file_I = basecase+"/sprof_file.dat"
  system("/bin/rm -f "+file_I)   ; remove any pre-existing file
  write_table(file_I, "w", [/date,dimsizes(s_prof),2/], "%s  %i  %i")
  write_table(file_I, "a", [/-depth_s,s_prof/], "%6.1f   %9.6f") ; append

end
