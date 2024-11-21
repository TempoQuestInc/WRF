subroutine da_read_obs_netcdf_gems2 (obstype,iv, infile)

   !---------------------------------------------------------------------------
   !  Purpose: read in OMS GEMS2 data to innovation structure
   !
   !   METHOD: use F90 sequential data structure to avoid reading file twice  
   !            so that da_scan_bufrtovs is not necessary any more.
   !            1. read file radiance data in sequential data structure
   !            2. do gross QC check
   !            3. assign sequential data structure to innovation structure
   !               and deallocate sequential data structure
   !   File History:
   !      Peiming Dong Added NPP atms, 2012/4/17
   !---------------------------------------------------------------------------

   implicit none

   character(5)      ,  intent (in)    :: obstype
   character(20)     ,  intent (in)    :: infile
   type (iv_type)    ,  intent (inout) :: iv

#ifdef BUFR

   integer          :: iost
   integer(i_kind), allocatable :: nread(:)

   logical outside, outside_all, iuse
   integer :: inst

   character(10) date
   integer(i_kind) i,j,k,ifov
   integer(i_kind) nchan
   integer :: num_netcdf(7), numnetcdf, inetcdf
   character(20) :: filename

   ! thinning variables
   integer(i_kind) itt,itx,iobs,iout
   real(r_kind) terrain,timedif,crit,dist
   real(r_kind) dlon_earth,dlat_earth

   real(r_kind) tbmin,tbmax, tbbad

   ! Instrument triplet, follow the convension of RTTOV 
   integer   :: platform_id, satellite_id, sensor_id

   ! pixel information
   integer   ::  year,month,day,hour,minute,second  ! observation time
   real*8    ::  obs_time
   ! real      ::  rlat, rlon                         !  lat/lon in degrees   for Anfovs
   real      ::  satzen, satazi, solzen ,solazi       !  scan angles for Anfovs
   integer   ::  landsea_mask
   real      ::  srf_height
   ! channels bright temperature
   real , allocatable ::   tb_inv(:)                    !  bright temperatures
   !  end type bright_temperature

   type (datalink_type), pointer    :: head, p, current, prev

   integer                        ::  ifgat
   type(info_type)                ::  info
   type(model_loc_type)           ::  loc


   data tbmin,tbmax,tbbad / 50.0_r_kind, 550.0_r_kind, -9.99e11_r_kind /
   integer :: num_tovs_local, num_tovs_file, num_tovs_global, num_tovs_selected
   integer :: num_tovs_thinned, num_tovs_used, num_tovs_used_tmp
   integer :: lnnetcdf
   integer :: n
   integer(i_kind), allocatable :: ptotal(:)
   real , allocatable :: in(:), out(:)
   logical :: found, head_found
   integer :: ncid, varid, ephemeris_data_gid, geolocation_gid, radiometric_data_gid, metadata_gid
   integer*8 :: nrecords
   integer*1, allocatable, dimension(:) :: on_earth
   real*8, allocatable, dimension(:) :: latitude,longitude, unix_timestamps
   real*8, allocatable, dimension(:) :: sat_az, sat_zen, sun_az, sun_zen
   real*8, allocatable, dimension(:) :: altitude
   real*8, allocatable, dimension(:,:) :: main_beam_antenna_temperatures
   integer :: rid
   real :: lat, lon

   call da_trace_entry("da_read_obs_netcdf_gems2")

   ! Initialize variables

   nchan = 24
   allocate(nread(1:rtminit_nsensor))
   allocate(ptotal(0:num_fgat_time))
   nread(1:rtminit_nsensor) = 0
   ptotal(0:num_fgat_time) = 0

   ! Set various variables depending on type of data to be read

   ! platform_id  = 1                 !! for NOAA series
   ! platform_id  = 10                !! for METOP series

   ! else if (mhs) then 
   !    sensor_id    =  15
   !    step   = 10.0_r_kind/9.0_r_kind
   !    start  = -445.0_r_kind/9.0_r_kind
   !    nchan=nchan_mhs
   !    subfgn='NC021027'
  sensor_id    =  -1
  nchan=nchan_gems2

   allocate (tb_inv(nchan))

   num_tovs_file     = 0    ! number of obs in file
   num_tovs_global   = 0    ! number of obs within whole domain
   num_tovs_local    = 0    ! number of obs within tile
   num_tovs_thinned  = 0    ! number of obs rejected by thinning
   num_tovs_used     = 0    ! number of obs entered into innovation computation
   num_tovs_selected = 0    ! number of obs limited for debuging
   iobs = 0                 ! for thinning, argument is inout

   ! 0.0  Open unit to satellite netcdf file and read file header
   !--------------------------------------------------------------

   num_netcdf(:)=0
   numnetcdf=0
   if (num_fgat_time>1) then
      do i=1,7
         call da_get_unit(lnnetcdf)
         write(filename,fmt='(A,2I1,A)') trim(infile),0,i,'.nc'
         open(unit   = lnnetcdf, FILE   = trim(filename),iostat =  iost, form = 'unformatted', STATUS = 'OLD')
         if (iost == 0) then
            numnetcdf=numnetcdf+1
	        num_netcdf(numnetcdf)=i
         else
            close (lnnetcdf)
         end if
         call da_free_unit(lnnetcdf)
      end do
   else
     numnetcdf=1
   end if
  
   if (numnetcdf==0) numnetcdf=1

netcdffile:  do inetcdf=1,numnetcdf   
   if (num_fgat_time==1) then
      filename=trim(infile)//'.nc'
   else
      if ((numnetcdf ==1) .and. (num_netcdf(inetcdf) == 0)) then
         filename=trim(infile)//'.nc'
      else
         write(filename,fmt='(A,2I1,A)') trim(infile),0,num_netcdf(inetcdf),'.nc'   
      end if
   end if


   ! SE - open file
   ierr=nf_open(trim(filename),nf_nowrite,ncid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error opening file "//trim(filename))

   ! SE - get group ids in netcdf file
   ierr=nf_inq_grp_ncid(ncid,'ephemeris_data',ephemeris_data_gid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing ephemeris_data group")
   ierr=nf_inq_grp_ncid(ncid,'geolocation_data',geolocation_gid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing geolocation_data group")
   ierr=nf_inq_grp_ncid(ncid,'radiometric_data',radiometric_data_gid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing radiometric_data group")
   ierr=nf_inq_grp_ncid(ncid,'metadata',metadata_gid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing metadata group")

   ! SE - read the number of records from the file
   ierr=nf_inq_varid(metadata_gid,'number_records',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing number_records variable")
   ierr=nf_get_var(metadata_gid,varid,nrecords)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing number_records variable")
   num_tovs_file=nrecords

   ! Loop to read netcdf file and assign information to a sequential structure
   !-------------------------------------------------------------------------

   if ( inetcdf == 1 ) then
      allocate (head)
      !  allocate ( head % tb_inv (1:nchan) )
      nullify  ( head % next )
      p => head
   endif

   if (tovs_start > 1) then
      write (unit=stdout,fmt='(A,I6)') "   Skipping tovs obs before", tovs_start
   end if

   ! SE - read variables from netcdf file
   ! allocate(VAR(nrecords))
   ! ierr=nf_inq_varid(ephemeris_data_gid,'VAR',varid)
   ! ierr=nf_get_var(ephemeris_data_gid,varid,VAR)
   allocate(main_beam_antenna_temperatures(nchan,nrecords))
   ierr=nf_inq_varid(radiometric_data_gid,'main_beam_antenna_temperatures',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing main_beam_antenna_temperatures variable")
   ierr=nf_get_var(radiometric_data_gid,varid,main_beam_antenna_temperatures)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing main_beam_antenna_temperatures variable")
   allocate(altitude(nrecords))
   ierr=nf_inq_varid(geolocation_gid,'altitude',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing altitude variable")
   ierr=nf_get_var(geolocation_gid,varid,altitude)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing altitude variable")
   allocate(sat_zen(nrecords))
   ierr=nf_inq_varid(ephemeris_data_gid,'sat_zen',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sat_zen variable")
   ierr=nf_get_var(ephemeris_data_gid,varid,sat_zen)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sat_zen variable")
   allocate(sat_az(nrecords))
   ierr=nf_inq_varid(ephemeris_data_gid,'sat_az',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sat_az variable")
   ierr=nf_get_var(ephemeris_data_gid,varid,sat_az)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sat_az variable")
   allocate(sun_zen(nrecords))
   ierr=nf_inq_varid(ephemeris_data_gid,'sun_zen',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sun_zen variable")
   ierr=nf_get_var(ephemeris_data_gid,varid,sun_zen)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sun_zen variable")
   allocate(sun_az(nrecords))
   ierr=nf_inq_varid(ephemeris_data_gid,'sun_az',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sun_az variable")
   ierr=nf_get_var(ephemeris_data_gid,varid,sun_az)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing sun_az variable")
   allocate(on_earth(nrecords))
   ierr=nf_inq_varid(geolocation_gid,'on_earth',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing on_earth variable")
   ierr=nf_get_var(geolocation_gid,varid,on_earth)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing on_earth variable")
   allocate(latitude(nrecords))
   ierr=nf_inq_varid(geolocation_gid,'latitude',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing latitude variable")
   ierr=nf_get_var(geolocation_gid,varid,latitude)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing latitude variable")
   allocate(longitude(nrecords))
   ierr=nf_inq_varid(geolocation_gid,'longitude',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing longitude variable")
   ierr=nf_get_var(geolocation_gid,varid,longitude)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing longitude variable")
   allocate(unix_timestamps(nrecords))
   ierr=nf_inq_varid(radiometric_data_gid,'unix_timestamps',varid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing unix_timestamps variable")
   ierr=nf_get_var(radiometric_data_gid,varid,unix_timestamps)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error accessing unix_timestamps variable")

   ! SE - close file
   ierr=nf_close(ncid)
   call check_netcdf_status(__FILE__, __LINE__, ierr, "error closing netcdf file")

   obs: do rid=1,nrecords
         ! write (unit=stdout,fmt=*) "Tb(ch=1)", main_beam_antenna_temperatures(1,rid), " for obs #", rid
        ! SE - 2.0 - read lat/lon of obs
        lat = real(latitude(rid))
        lon = real(longitude(rid))
        if (on_earth(rid) == 0) then
            ! write (unit=stdout,fmt=*) "   Obs point ", rid, "is not on Earth :: on_earth(rid) = ", on_earth(rid), " :: <lat,lon> = <", lat, ",", lon, ">"
            cycle
        endif
        info%lat = lat
        info%lon = lon
        call da_llxy (info, loc, outside, outside_all)
        if (outside_all) then
            ! write (unit=stdout,fmt=*) "   Obs point ", rid, "is outside of domain :: on_earth(rid) = ", on_earth(rid), " :: <lat,lon> = <", lat, ",", lon, ">"
            cycle
        endif

         !  3.0     Extract other information
         !------------------------------------------------------
         !  3.1     Extract satellite id and scan position. 
         ! SE - skip this for now. I think that these variables may not be necessary since they are used to calculate products that
         ! already exist in the GEMS dataset

        ! SE - 3.2 - extract date info
        call unix_to_date(unix_timestamps(rid), year, month, day, hour, minute, second)
        call da_get_julian_time(year,month,day,hour,minute,obs_time)

        if (obs_time < time_slots(0) .or.  &
            obs_time >= time_slots(num_fgat_time)) then
            ! write (unit=stdout,fmt='(A,I6,A,I4,I2,I2,I2,I2,I2)') "   Obs point ", rid, "is outside of time range :: ", year, month, day, hour, minute, second
            cycle
        ! else
        !     write (unit=stdout,fmt='(A,I6,A,I4,I2,I2,I2,I2,I2)') "   Obs point ", rid, "is inside of time range :: ", year, month, day, hour, minute, second
        endif

         ! 3.2.1 determine FGAT index ifgat
   
         do ifgat=1,num_fgat_time
            if (obs_time >= time_slots(ifgat-1) .and.  &
                obs_time  < time_slots(ifgat)) exit
         end do

         ! 3.3 Find wrfvar instrument index from RTTOV instrument triplet
         !     go to next data if id is not in the lists

         ! SE - These are dummy values for now that shouldn't overlap with existing platform/sat/sensor options. I'll come back to
         ! this once I have a better understanding of the conventions around these
         platform_id = 100
         satellite_id = 100
         sensor_id = 100

         inst = 0
         do i = 1, rtminit_nsensor
            if (platform_id  == rtminit_platform(i) &
               .and. satellite_id == rtminit_satid(i)    &
               .and. sensor_id    == rtminit_sensor(i)) then
               inst = i
               exit
            end if
         end do
         if (inst == 0) cycle

         ! 3.4 extract satellite and solar angle
         satzen = real(sat_zen(rid))*rad2deg
         satazi = real(sat_az(rid))*rad2deg
         solzen = real(sun_zen(rid))*rad2deg
         solazi = real(sun_az(rid))*rad2deg

         if ( satzen > 65.0 ) then   ! CRTM has a satzen > 65.0 check
            ! write (unit=stdout,fmt=*) "   Obs point ", rid, "obs angle is too shallow :: <satzen, satazi, solzen, solazi> = <", satzen, satazi, solzen, solazi, ">"
            cycle
        endif

         num_tovs_global = num_tovs_global + 1
         ptotal(ifgat) = ptotal(ifgat) + 1

         if (num_tovs_global < tovs_start) then
            cycle
         end if

         if (num_tovs_global > tovs_end) then
            write (unit=stdout,fmt='(A,I6)') "   Skipping radiance obs after", tovs_end
            exit obs
         end if

         num_tovs_selected = num_tovs_selected + 1
 
         if (num_tovs_selected > max_tovs_input) then
            write(unit=message(1),fmt='(A,I10,A)') &
               "Max number of tovs",max_tovs_input," reached"
            call da_warning(__FILE__,__LINE__,message(1:1))
            num_tovs_selected = num_tovs_selected-1
            num_tovs_global   = num_tovs_global-1
            ptotal(ifgat) = ptotal(ifgat) - 1
            exit obs
         end if

         if (outside) cycle ! No good for this PE
         num_tovs_local = num_tovs_local + 1


         !  Make Thinning
         !  Map obs to thinning grid
         !-------------------------------------------------------------------
         if (thinning) then
            dlat_earth = info%lat
            dlon_earth = info%lon
            if (dlon_earth<zero) dlon_earth = dlon_earth+r360
            if (dlon_earth>=r360) dlon_earth = dlon_earth-r360
            dlat_earth = dlat_earth*deg2rad
            dlon_earth = dlon_earth*deg2rad           
            timedif = 0.0 !2.0_r_kind*abs(tdiff)        ! range:  0 to 6
            terrain = 0.01_r_kind*abs(altitude(rid))
            crit = 1.0 !0.01_r_kind+terrain + timedif !+ 10.0_r_kind*float(iskip)
            call map2grids(inst,ifgat,dlat_earth,dlon_earth,crit,iobs,itx,1,itt,iout,iuse)
            if (.not. iuse) then
               num_tovs_thinned=num_tovs_thinned+1
               cycle
            end if
         end if

         num_tovs_used = num_tovs_used + 1
         nread(inst) = nread(inst) + 1

         ! 3.5 extract surface information
         srf_height = real(altitude(rid))          ! station height
         if (srf_height < 8888.0 .AND. srf_height > -416.0) then
         else
            srf_height = 0.0
         endif  

         ! SE - need to fix this once input file has land/sea mask
         landsea_mask = 1  ! 0:land ; 1:sea (same as RTTOV)
         ! rmask=one                          ! land
         ! if (nint(bfr1bhdr(netcdf_landsea_mask))==1) rmask=0.0_r_kind   ! reverse the land/sea mask in netcdf
         ! landsea_mask = rmask+.001_r_kind             ! land sea mask

         info%elv = srf_height

         ! 3.6 extract channel bright temperature
   
         do k = 1, nchan
            tb_inv(k) = real(main_beam_antenna_temperatures(k,rid))
         enddo
         do k = 1, nchan
            if ( tb_inv(k) < tbmin .or. tb_inv(k) > tbmax) &
               tb_inv(k) = missing_r
         end do
         if ( all(tb_inv<0.0) ) then
            num_tovs_local = num_tovs_local -1
            num_tovs_used = num_tovs_used - 1
            nread(inst) = nread(inst) - 1
            cycle
         end if

         !  4.0   assign information to sequential radiance structure
         !--------------------------------------------------------------------------
         allocate (p % tb_inv (1:nchan))
         p%info             = info
         p%loc              = loc
         p%landsea_mask     = landsea_mask
         p%scanpos          = ifov
         p%satzen           = satzen
         p%satazi           = satazi
         p%solzen           = solzen
         p%tb_inv(1:nchan)  = tb_inv(1:nchan)
         p%sensor_index     = inst
         p%ifgat            = ifgat
!RTTOV9_3
         p%solazi           = solazi
 !end of RTTOV9_3
         allocate (p%next)   ! add next data

         p => p%next
         nullify (p%next)


   enddo obs


   ! SE - deallocate vaiables
   deallocate(main_beam_antenna_temperatures)
   deallocate(altitude)
   deallocate(sat_zen)
   deallocate(sat_az)
   deallocate(sun_zen)
   deallocate(sun_az)
   deallocate(on_earth)
   deallocate(latitude)
   deallocate(longitude)
   deallocate(unix_timestamps)


end do netcdffile

   if (thinning .and. num_tovs_global > 0 ) then

#ifdef DM_PARALLEL

      ! Get minimum crit and associated processor index.
      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            j = j + thinning_grid(n,ifgat)%itxmax
         end do
      end do

      allocate ( in  (j) )
      allocate ( out (j) )

      j = 0
      do ifgat = 1, num_fgat_time
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               in(j) = thinning_grid(n,ifgat)%score_crit(i) 
            end do
         end do
      end do
      call mpi_reduce(in, out, j, true_mpi_real, mpi_min, root, comm, ierr)

      call wrf_dm_bcast_real (out, j)

      j = 0
      do ifgat = 1, num_fgat_time 
         do n = 1, iv%num_inst
            do i = 1, thinning_grid(n,ifgat)%itxmax
               j = j + 1
               if ( ABS(out(j)-thinning_grid(n,ifgat)%score_crit(i)) > 1.0E-10 ) thinning_grid(n,ifgat)%ibest_obs(i)  = 0
            end do
         end do
      end do

      deallocate( in  )
      deallocate( out )

#endif

      ! Delete the nodes which being thinning out
      p => head
      prev => head
      head_found = .false.
      num_tovs_used_tmp = num_tovs_used
      do j = 1, num_tovs_used_tmp
         n = p%sensor_index
         ifgat = p%ifgat 
         found = .false.

         do i = 1, thinning_grid(n,ifgat)%itxmax
            if ( thinning_grid(n,ifgat)%ibest_obs(i) == j .and. thinning_grid(n,ifgat)%score_crit(i) < 9.99e6_r_kind ) then
               found = .true.
               exit
            endif
         end do 
        
         ! free current data
         if ( .not. found ) then
            current => p
            p => p%next
            if ( head_found ) then
               prev%next => p
            else
               head => p
               prev => p
            endif
            deallocate ( current % tb_inv )
            deallocate ( current )
            num_tovs_thinned = num_tovs_thinned + 1
            num_tovs_used = num_tovs_used - 1
            nread(n) = nread(n) - 1
            continue
         endif

         if ( found .and. head_found ) then
            prev => p
            p => p%next
            continue
         endif

         if ( found .and. .not. head_found ) then
            head_found = .true.
            head => p
            prev => p
            p => p%next
         endif
        
      end do

   endif  ! End of thinning

   iv%total_rad_pixel   = iv%total_rad_pixel   + num_tovs_used
   iv%total_rad_channel = iv%total_rad_channel + num_tovs_used*nchan

   iv%info(radiance)%nlocal = iv%info(radiance)%nlocal + num_tovs_used
   iv%info(radiance)%ntotal = iv%info(radiance)%ntotal + num_tovs_global

   do i = 1, num_fgat_time
      ptotal(i) = ptotal(i) + ptotal(i-1) 
      iv%info(radiance)%ptotal(i) = iv%info(radiance)%ptotal(i) + ptotal(i) 
   end do
   if ( iv%info(radiance)%ptotal(num_fgat_time) /= iv%info(radiance)%ntotal ) then
      write(unit=message(1),fmt='(A,I10,A,I10)') &
          "Number of ntotal:",iv%info(radiance)%ntotal," is different from the sum of ptotal:", iv%info(radiance)%ptotal(num_fgat_time)
      call da_warning(__FILE__,__LINE__,message(1:1))
   endif

   write(unit=stdout,fmt='(a)') 'num_tovs_file num_tovs_global num_tovs_local num_tovs_used num_tovs_thinned'
   write(unit=stdout,fmt='(5i10)') num_tovs_file,num_tovs_global, num_tovs_local,num_tovs_used,num_tovs_thinned

   deallocate(tb_inv)  

   !  5.0 allocate innovation radiance structure
   !----------------------------------------------------------------  
   
   do i = 1, iv%num_inst
      if (nread(i) < 1) cycle
      iv%instid(i)%num_rad = nread(i)
      iv%instid(i)%info%nlocal = nread(i)
      write(UNIT=stdout,FMT='(a,i3,2x,a,3x,i10)') &
        'Allocating space for radiance innov structure', &
         i, iv%instid(i)%rttovid_string, iv%instid(i)%num_rad

      call da_allocate_rad_iv(i,nchan,iv)

   end do
   
   !  6.0 assign sequential structure to innovation structure
   !-------------------------------------------------------------
   nread(1:rtminit_nsensor) = 0 
   p => head
   ! do while ( associated(p) )

   do n = 1, num_tovs_used
      i = p%sensor_index
      nread(i) = nread(i) + 1

      call da_initialize_rad_iv (i, nread(i), iv, p)

      current => p
      p => p%next

      ! free current data
      deallocate ( current % tb_inv )
      deallocate ( current )
   end do

   deallocate ( p )

   deallocate (nread)
   deallocate (ptotal)

   ! check if sequential structure has been freed
   !
   ! p => head
   ! do i = 1, num_rad_selected
   !    write (unit=stdout,fmt=*)  i, p%tb_inv(1:nchan)
   !    p => p%next
   ! end do

   call da_trace_exit("da_read_obs_netcdf_gems2")
#else
   call da_error(__FILE__,__LINE__,(/"Needs to be compiled with a BUFR library"/)) 
#endif
  
contains

  ! Subroutine for error checking NetCDF status
  subroutine check_netcdf_status(file, line, status, action)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(in) :: status
    character(len=*), intent(in) :: action

    if (status /= nf_noerr) then
       call da_trace_exit("da_read_obs_netcdf_gems2")
       call da_error(file,line,(/"Error "//action//":"//trim(nf_strerror(status))/))
    end if
  end subroutine check_netcdf_status

    subroutine unix_to_date(unix_time, year, month, day, hour, minute, second)
        implicit none
        real*8, intent(in) :: unix_time
        integer, intent(out) :: year, month, day, hour, minute, second

        integer :: days, seconds_in_day, leap_year, days_in_month
        integer, parameter :: SECONDS_PER_DAY = 86400
        integer, dimension(12) :: month_days_normal = (/31, 28, 31, 30, 31, 30, &
                                                        31, 31, 30, 31, 30, 31/)
        integer, dimension(12) :: month_days_leap = (/31, 29, 31, 30, 31, 30, &
                                                      31, 31, 30, 31, 30, 31/)
        integer, dimension(12) :: days_in_months

        ! Calculate total days and remaining seconds in the day
        days = int(unix_time / SECONDS_PER_DAY)
        seconds_in_day = mod(unix_time, SECONDS_PER_DAY)

        ! Calculate hour, minute, and second
        hour = seconds_in_day / 3600
        minute = mod(seconds_in_day, 3600) / 60
        second = mod(seconds_in_day, 60)

        ! Start from the Unix epoch year
        year = 1970

        ! Loop to find the correct year
        do while (days >= 365)
            leap_year = (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0))
            if (leap_year) then
                if (days >= 366) then
                    days = days - 366
                    year = year + 1
                else
                    exit
                endif
            else
                days = days - 365
                year = year + 1
            endif
        end do

        ! Determine the month and day
        leap_year = (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0))

        ! Set days_in_months based on whether it's a leap year
        if (leap_year) then
            days_in_months = month_days_leap
        else
            days_in_months = month_days_normal
        endif

        month = 1
        do while (days >= days_in_months(month))
            days = days - days_in_months(month)
            month = month + 1
        end do

        day = days + 1

    end subroutine unix_to_date


end subroutine da_read_obs_netcdf_gems2


