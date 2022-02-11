
; This mk_input_ls_forcing.ncl 
;  reads input_ls_forcing.txt and creates input_ls_forcing.nc
;  in the format defined in input_ls_forcing.cdl
;  for idealized crm simulation. 
;
;  sendo at bnl.gov

; 12/12/2011 sendo  included u,v and tau


begin

  infile   =  "input_ls_forcing.txt"

  outfile = "input_ls_forcing.nc"
  cdlfile = "input_ls_forcing.cdl"

  nz = 4 ;

  dt = 10800; sec. for the calculation of tendency and time
 
  tlast = 54000 ; sec at last line

;#-------END USER MODIFICATIONS-----------------

 
  nt = tlast / dt + 1  ; num of time steps 


  ;# check for existence of cdl file
  ffile = systemfunc("ls "+cdlfile) 

  if ( ismissing(ffile) ) then

    print("Please supply a template "+cdlfile+" that is consistent with "+infile)
    exit

  end if


  ; create forcing file

  ierr = systemfunc("/bin/rm -f "+outfile)
  ierr = systemfunc("ncgen -o "+outfile+" "+cdlfile)


  ; open output file 

  oFl = addfile(outfile,"rw")


  ; read input file

  strs = asciiread(infile, -1, "string")

  delim = " "
  wrftime1 = str_get_field(strs, 1, delim)
  z1       = str_get_field(strs, 2, delim)
  w1       = str_get_field(strs, 3, delim)
  thadv1   = str_get_field(strs, 4, delim)
  qvadv1   = str_get_field(strs, 5, delim)
  thrlx1   = str_get_field(strs, 6, delim)
  qvrlx1   = str_get_field(strs, 7, delim)

  delete(strs)

  wrftime = onedtond(wrftime1, (/nt, nz/))
  z       = onedtond(z1,       (/nt, nz/))
  w       = onedtond(w1,       (/nt, nz/))
  thadv   = onedtond(thadv1,   (/nt, nz/))
  qvadv   = onedtond(qvadv1,   (/nt, nz/))
  thrlx   = onedtond(thrlx1,   (/nt, nz/))
  qvrlx   = onedtond(qvrlx1,   (/nt, nz/))

  
  ; * Get the length of the time strings

  dims = filevardimsizes(oFl,"Times")
  DateLen = dims(1)
  delete(dims)


  ; string to char

  DateStr = stringtochar(wrftime(:,0))
  ;print (DateStr)

  ; string to float

  zf      = stringtofloat(z)
  wf      = stringtofloat(w)
  thadvf  = stringtofloat(thadv)
  qvadvf  = stringtofloat(qvadv)
  thrlxf  = stringtofloat(thrlx)
  qvrlxf  = stringtofloat(qvrlx)

  ;print (zf)


  ; substitute  

  do itime = 0, nt-1

    oFl->Times(itime,:) = (/DateStr(itime,0:DateLen-1)/)

    ; large scale
    oFl->Z_LS(itime,:) = (/zf(itime,:)/)
    oFl->W_LS(itime,:) = (/wf(itime,:)/)
    oFl->TH_ADV(itime,:)  = (/thadvf(itime,:)/)
    oFl->TH_RLX(itime,:)  = (/thrlxf(itime,:)/)
    oFl->QV_ADV(itime,:)  = (/qvadvf(itime,:)/)
    oFl->QV_RLX(itime,:)  = (/qvrlxf(itime,:)/)

    if (itime .lt. nt - 1) then
      tenzf = (zf(itime + 1,:) - zf(itime,:)) / dt
      tenwf = (wf(itime + 1,:) - wf(itime,:)) / dt
 
      tenthadvf = (thadvf(itime + 1,:) - thadvf(itime,:)) / dt
      tenqvadvf = (qvadvf(itime + 1,:) - qvadvf(itime,:)) / dt
      tenthrlxf = (thrlxf(itime + 1,:) - thrlxf(itime,:)) / dt
      tenqvrlxf = (qvrlxf(itime + 1,:) - qvrlxf(itime,:)) / dt
    else 
      tenzf = (zf(itime,:) - zf(itime - 1,:)) / dt
      tenwf = (wf(itime,:) - wf(itime - 1,:)) / dt
 
      tenthadvf = (thadvf(itime,:) - thadvf(itime - 1,:)) / dt
      tenqvadvf = (qvadvf(itime,:) - qvadvf(itime - 1,:)) / dt
      tenthrlxf = (thrlxf(itime,:) - thrlxf(itime - 1,:)) / dt
      tenqvrlxf = (qvrlxf(itime,:) - qvrlxf(itime - 1,:)) / dt
    end if

    oFl->Z_LS_TEND(itime,:)    = (/tenzf/)
    oFl->W_LS_TEND(itime,:)    = (/tenwf/)

    oFl->TH_ADV_TEND(itime,:)  = (/tenthadvf/)
    oFl->TH_RLX_TEND(itime,:)  = (/tenthrlxf/)

    oFl->QV_ADV_TEND(itime,:)  = (/tenqvadvf/)
    oFl->QV_RLX_TEND(itime,:)  = (/tenqvrlxf/)

    oFl->U_LS(itime,:)    = (/0.0/)
    oFl->V_LS(itime,:)    = (/0.0/)

    oFl->U_LS_TEND(itime,:)    = (/0.0/)
    oFl->V_LS_TEND(itime,:)    = (/0.0/)

    oFl->INV_TAU_S(itime,:) = (/9999.9/)
    oFl->INV_TAU_M(itime,:) = (/9999.9/)



  end do

end

