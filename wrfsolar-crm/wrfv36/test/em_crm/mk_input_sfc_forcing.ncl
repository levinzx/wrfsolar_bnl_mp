
; This mk_input_sfc_forcing.ncl 
;  reads input_sfc_forcing.txt and creates input_sfc_forcing.nc
;  in the format defined in input_forcing.cdl
;  for idealized crm simulation. 
;
;  sendo at bnl

begin

  infile  = "input_sfc_forcing.txt"

  outfile = "input_sfc_forcing.nc"
  cdlfile = "input_sfc_forcing.cdl"

  dt    = 1800  ; time step sec 
  tlast = 54000 ; sec at last line


;#-------END USER MODIFICATIONS-----------------

  nt = tlast / dt + 1  ; num of time steps 


  ; check for existence of cdl file

  ffile = systemfunc("ls "+cdlfile) 

  if ( ismissing(ffile) ) then

    print("Please supply a template "+cdlfile+" that is consistent with "+sfc_infile)
    exit

  end if


  ; create forcing file

  ierr = systemfunc("/bin/rm -f "+outfile)
  ierr = systemfunc("ncgen -o "+outfile+" "+cdlfile)


  ; open output file 

  oFl = addfile(outfile,"rw")


  ;  read input file

  strs = asciiread(infile, -1, "string")

  delim = " "
  wrftime = str_get_field(strs, 1, delim)
  shflx   = str_get_field(strs, 2, delim)
  lhflx   = str_get_field(strs, 3, delim)
  albd    = str_get_field(strs, 4, delim)
  tsk     = str_get_field(strs, 5, delim)

  ;print (shflx)
  delete(strs)


  ; * Get the length of the time strings

  dims = filevardimsizes(oFl,"Times")
  DateLen = dims(1)
  delete(dims)


  ; string to char

  DateStr = stringtochar(wrftime)
  

  ; string to float

  shflxf = stringtofloat(shflx)
  lhflxf = stringtofloat(lhflx)
  albdf  = stringtofloat(albd)
  tskf   = stringtofloat(tsk)
  ;print(fshflx)


  ; substitute
 
  do itime = 0, nt-1
 
    ;print (nt)

    oFl->Times(itime,:) = (/DateStr(itime,0:DateLen-1)/)

    ; surface
    oFl->PRE_SH_FLX(itime)  = (/shflxf(itime)/)
    oFl->PRE_LH_FLX(itime)  = (/lhflxf(itime)/)
    oFl->PRE_ALBEDO(itime)  = (/albdf(itime)/)
    oFl->PRE_TSK(itime)     = (/tskf(itime)/)

   if (itime .lt. nt - 1) then
      tenshflxf = (shflxf(itime + 1) - shflxf(itime)) / dt
      tenlhflxf = (lhflxf(itime + 1) - lhflxf(itime)) / dt
      tenalbdf  = (albdf(itime + 1)  - albdf(itime)) / dt
      tentskf   = (tskf(itime  + 1)  - tskf(itime))  / dt
    else 
      tenshflxf = (shflxf(itime) - shflxf(itime - 1)) / dt
      tenlhflxf = (lhflxf(itime) - lhflxf(itime - 1)) / dt
      tenalbdf  = (albdf(itime)  - albdf(itime - 1)) / dt
      tentskf   = (tskf(itime)   - tskf(itime - 1))  / dt
    end if

    oFl->PRE_SH_FLX_TEND(itime)  = (/tenshflxf/)
    oFl->PRE_LH_FLX_TEND(itime)  = (/tenlhflxf/)

    oFl->PRE_ALBEDO_TEND(itime)  = (/tenalbdf/)
    oFl->PRE_TSK_TEND(itime)     = (/tentskf/)

  end do

end

