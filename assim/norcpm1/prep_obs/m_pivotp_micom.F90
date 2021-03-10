module m_pivotp_micom
  use netcdf
  use nfw_mod
  use m_spherdist
  implicit none

contains
  ! F. Counillon (adapted from an algorithm of Mats Bentsen) 
  ! This subroutine search the pivot point for a given observations
  ! The search is linear moving toward the neigboring grid cell that minimize
  ! the distance to the obs. This search is in the worst case in O(n). The input
  ! ipiv, jpiv corresponds to the pivot point from the previous search. If there
  ! is a kind of order in the way the observation are given the search will be
  ! very fast.
  !
  subroutine pivotp_micom(lon, lat,modlon,modlat, ipiv, jpiv, nx, ny, min_r, max_r, &
  itw, jtw, its, jts, itn, jtn, ite, jte)
   real, intent(in) ::  lon, lat
   integer, intent(in) :: nx, ny
   real, intent(in), dimension(nx,ny) :: modlon,modlat, min_r, max_r
   integer, intent(in), dimension(nx,ny) :: itw,jtw, &
        its, jts, itn, jtn, ite, jte
   integer, intent(inout) :: ipiv, jpiv
   real*8 :: min_d, d
   integer :: i, j, ito, jto

   min_d = spherdist(modlon(ipiv,jpiv), modlat(ipiv,jpiv), lon, lat)
   do while (min_d > min_r(ipiv,jpiv))

      ito = ipiv
      jto = jpiv

      i = itw(ito,jto)
      j = jtw(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = ite(ito,jto)
      j = jte(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = its(ito,jto)
      j = jts(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = itn(ito,jto)
      j = jtn(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif

      if (ipiv == ito .and. jpiv == jto) exit

   enddo
end subroutine pivotp_micom


! Y. WANG (adapted from an algorithm of F. COUNILLON)                                             
! This subroutine search the pivot point for a given observation                                                                             
! The search is linear moving toward the neigboring grid cell that minimize
!  
subroutine pivotp_micom_new(lon, lat,modlon,modlat, ipiv, jpiv, nx, ny, min_r, max_r, &
     itw, jtw, itn, jtn, its, jts, ite, jte)
  real, intent(in) ::  lon, lat
  integer, intent(in) :: nx, ny
  real, intent(in), dimension(nx,ny) :: modlon,modlat, min_r, max_r
  integer, intent(in), dimension(nx,ny) :: itw,jtw, &
       its, jts, itn, jtn, ite, jte
  integer, intent(inout) :: ipiv, jpiv
  real*8 :: min_d, d
  integer :: i, j, ito, jto

  real :: xx(4), yy(4)
  logical :: inside

  min_d = spherdist(modlon(ipiv,jpiv), modlat(ipiv,jpiv), lon, lat)
  do while (min_d > min_r(ipiv,jpiv))

     ito = ipiv
     jto = jpiv
     
     i = itw(ito,jto)
     j = jtw(ito,jto)
     d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
     if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
     endif

     i = ite(ito,jto)
     j = jte(ito,jto)
     d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
     if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
     endif

     i = its(ito,jto)
     j = jts(ito,jto)
     d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
     if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
     endif

     i = itn(ito,jto)
     j = jtn(ito,jto)
     d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
     if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
     endif

     if (ipiv == ito .and. jpiv == jto) exit

  enddo

  ito = ipiv
  jto = jpiv
  if (jpiv == ny) then
     ! check two grids
     i = ito
     j = jto
     xx(1) = modlon(i, j)
     xx(2) = modlon(itw(i,j), jtw(i,j))
     xx(3) = modlon(its(itw(i,j), jtw(i,j)), jts(itw(i,j), jtw(i,j)))
     xx(4) = modlon(its(i,j), jts(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(itw(i,j), jtw(i,j))
     yy(3) = modlat(its(itw(i,j), jtw(i,j)), jts(itw(i,j), jtw(i,j)))
     yy(4) = modlat(its(i,j), jts(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = its(itw(i,j), jtw(i,j))
        jpiv = jts(itw(i,j), jtw(i,j))
        return
     end if
     xx(1) = modlon(i, j)
     xx(2) = modlon(ite(i,j), jte(i,j))
     xx(3) = modlon(its(ite(i,j), jte(i,j)), jts(ite(i,j), jte(i,j)))
     xx(4) = modlon(its(i,j), jts(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(ite(i,j), jte(i,j))
     yy(3) = modlat(its(ite(i,j), jte(i,j)), jts(ite(i,j), jte(i,j)))
     yy(4) = modlat(its(i,j), jts(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = its(i,j)
        jpiv = jts(i,j)
        return
     else
        jpiv = ny
     end if
  else if (jpiv > 1) then 
     ! check four gird
     i = ito
     j = jto
     
     ! check north-east                                                                                                                          
     xx(1) = modlon(i, j)
     xx(2) = modlon(ite(i,j), jte(i,j))
     xx(3) = modlon(itn(ite(i,j), jte(i,j)), jtn(ite(i,j), jte(i,j)))
     xx(4) = modlon(itn(i,j), jtn(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(ite(i,j), jte(i,j))
     yy(3) = modlat(itn(ite(i,j), jte(i,j)), jtn(ite(i,j), jte(i,j)))
     yy(4) = modlat(itn(i,j), jtn(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = i
        jpiv = j
        return
     end if
     ! check south-west
     xx(1) = modlon(i, j)
     xx(2) = modlon(itw(i,j), jtw(i,j))
     xx(3) = modlon(its(itw(i,j), jtw(i,j)), jts(itw(i,j), jtw(i,j)))
     xx(4) = modlon(its(i,j), jts(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(itw(i,j), jtw(i,j))
     yy(3) = modlat(its(itw(i,j), jtw(i,j)), jts(itw(i,j), jtw(i,j)))
     yy(4) = modlat(its(i,j), jts(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = its(itw(i,j), jtw(i,j))
        jpiv = jts(itw(i,j), jtw(i,j))
        return
     end if
     ! check south-east
     xx(1) = modlon(i, j)
     xx(2) = modlon(ite(i,j), jte(i,j))
     xx(3) = modlon(its(ite(i,j), jte(i,j)), jts(ite(i,j), jte(i,j)))
     xx(4) = modlon(its(i,j), jts(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(ite(i,j), jte(i,j))
     yy(3) = modlat(its(ite(i,j), jte(i,j)), jts(ite(i,j), jte(i,j)))
     yy(4) = modlat(its(i,j), jts(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = its(i,j)
        jpiv = jts(i,j)
        return
     end if
     ! check north-west
     xx(1) = modlon(i, j)
     xx(2) = modlon(itw(i,j), jtw(i,j))
     xx(3) = modlon(itn(itw(i,j), jtw(i,j)), jtn(itw(i,j), jtw(i,j)))
     xx(4) = modlon(itn(i,j), jtn(i,j))
     yy(1) = modlat(i, j)
     yy(2) = modlat(itw(i,j), jtw(i,j))
     yy(3) = modlat(itn(itw(i,j), jtw(i,j)), jtn(itw(i,j), jtw(i,j)))
     yy(4) = modlat(itn(i,j), jtn(i,j))
     call check_point(lon,lat,xx,yy,inside)
     if (inside) then
        ipiv = itw(i,j)
        jpiv = jtw(i,j)
        return
     else
        print *, 'ERROR: pivotp_micom_new()', ipiv, jpiv
        print *, lon, lat
        stop
     end if
  else
     return
  end if
end subroutine pivotp_micom_new
   
subroutine ini_pivotp(modlon,modlat, nx, ny, min_r, max_r, itw, jtw, itn, jtn, its, jts, ite, jte)
   integer, intent(in) :: nx, ny
   real, intent(in), dimension(nx,ny) :: modlon,modlat
   real, intent(out), dimension(nx,ny):: min_r, max_r
   integer, intent(out), dimension(nx,ny) :: itw, jtw, &
        itn, jtn, its, jts, ite, jte
   integer :: ncid,vITW_ID, vJTW_ID, vITE_ID ,vJTE_ID ,vITS_ID,  &
      vJTS_ID, vITN_ID, vJTN_ID, vVCLON_ID, vVCLAT_ID, vUCLON_ID &
      , vUCLAT_ID, vTCLON_ID, vTCLAT_ID ,ii,jj

   real, dimension(nx,ny,4) :: vclon, vclat, uclon, uclat, tclon, tclat 
   call nfw_open('grid.nc', nf_nowrite, ncid)

   call nfw_inq_varid('grid.nc', ncid,'inw' ,vITW_ID)
   call nfw_inq_varid('grid.nc', ncid,'jnw' ,vJTW_ID)
   call nfw_inq_varid('grid.nc', ncid,'ine' ,vITE_ID)
   call nfw_inq_varid('grid.nc', ncid,'jne' ,vJTE_ID)
   call nfw_inq_varid('grid.nc', ncid,'ins' ,vITS_ID)
   call nfw_inq_varid('grid.nc', ncid,'jns' ,vJTS_ID)
   call nfw_inq_varid('grid.nc', ncid,'inn' ,vITN_ID)
   call nfw_inq_varid('grid.nc', ncid,'jnn' ,vJTN_ID)
   call nfw_inq_varid('grid.nc', ncid,'vclon' ,vVCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'vclat' ,vVCLAT_ID)
   call nfw_inq_varid('grid.nc', ncid,'uclon' ,vUCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'uclat' ,vUCLAT_ID)
   call nfw_inq_varid('grid.nc', ncid,'pclon' ,vTCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'pclat' ,vTCLAT_ID)
   call nfw_get_var_int('grid.nc', ncid, vITW_ID, itw)
   call nfw_get_var_int('grid.nc', ncid, vJTW_ID, jtw)
   call nfw_get_var_int('grid.nc', ncid, vITS_ID, its)
   call nfw_get_var_int('grid.nc', ncid, vJTS_ID, jts)
   call nfw_get_var_int('grid.nc', ncid, vITN_ID, itn)
   call nfw_get_var_int('grid.nc', ncid, vJTN_ID, jtn)
   call nfw_get_var_int('grid.nc', ncid, vITE_ID, ite)
   call nfw_get_var_int('grid.nc', ncid, vJTE_ID, jte)
   call nfw_get_var_double('grid.nc', ncid, vVCLON_ID, vclon)
   call nfw_get_var_double('grid.nc', ncid, vVCLAT_ID, vclat)
   call nfw_get_var_double('grid.nc', ncid, vUCLON_ID, uclon)
   call nfw_get_var_double('grid.nc', ncid, vUCLAT_ID, uclat)
   call nfw_get_var_double('grid.nc', ncid, vTCLON_ID, tclon)
   call nfw_get_var_double('grid.nc', ncid, vTCLAT_ID, tclat)
   do jj = 1, ny
         do ii = 1, nx
            min_r(ii,jj) =                                   &
               min(spherdist(vclon(ii,jj,4), vclat(ii,jj,4), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(vclon(ii,jj,3), vclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(uclon(ii,jj,2), uclat(ii,jj,2), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(uclon(ii,jj,3), uclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)))
            max_r(ii,jj) =                                   &
               max(spherdist(tclon(ii,jj,1), tclat(ii,jj,1), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,2), tclat(ii,jj,2), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,3), tclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,4), tclat(ii,jj,4), &
                              modlon(ii,jj), modlat(ii,jj)))
         enddo
   enddo
end subroutine ini_pivotp

! Y. WANG: check if point in model grid (irregulaer four points)
subroutine check_point(x, y, xx, yy, inside)
  real, intent(in) :: x, y
  real, intent(in) :: xx(4), yy(4)
  logical, intent(inout) :: inside
  ! distance between (xx(i), yy(i)) and (x,y)
  real :: a, b, c
  ! semi-parameter
  real :: s
  ! area
  real :: sum_triangle_area
  real :: rectangle_area

  ! initialasation 
  inside = .false.
  sum_triangle_area = 0.
  rectangle_area = 0.
  
  ! first trangle area with (xx(1), yy(1)), (xx(2), yy(2)) and (x,y) 
  a = periodic_sqrt(xx(1), yy(1), x, y)
  b = periodic_sqrt(xx(2), yy(2), x, y)
  c = periodic_sqrt(xx(1), yy(1), xx(2), yy(2))
  s = (a + b + c) / 2
  sum_triangle_area = sum_triangle_area + sqrt(s*(s-a)*(s-b)*(s-c))
  
  ! second trangle area with (xx(2), yy(2)), (xx(3), yy(3)) and (x,y) 
  a = periodic_sqrt(xx(3), yy(3), x, y)
  b = periodic_sqrt(xx(2), yy(2), x, y)
  c = periodic_sqrt(xx(3), yy(3), xx(2), yy(2))
  s = (a + b + c) / 2
  sum_triangle_area = sum_triangle_area + sqrt(s*(s-a)*(s-b)*(s-c))
  
  ! third trangle area with (xx(3), yy(3)), (xx(4), yy(4)) and (x,y)                                                                           
  a = periodic_sqrt(xx(3), yy(3), x, y)
  b = periodic_sqrt(xx(4), yy(4), x, y)
  c = periodic_sqrt(xx(3), yy(3), xx(4), yy(4))
  s = (a + b + c) / 2
  sum_triangle_area = sum_triangle_area + sqrt(s*(s-a)*(s-b)*(s-c))

 ! fourth trangle area with (xx(1), yy(1)), (xx(4), yy(4)) and (x,y)
  a = periodic_sqrt(xx(1), yy(1), x, y)
  b = periodic_sqrt(xx(4), yy(4), x, y)
  c = periodic_sqrt(xx(1), yy(1), xx(4), yy(4))
  s = (a + b + c) / 2
  sum_triangle_area = sum_triangle_area + sqrt(s*(s-a)*(s-b)*(s-c))

  ! rectangle area = triangle_area((xx(1),yy(1)), (xx(2),yy(2)), (xx(3),yy(3))) 
  !                + triangle_area((xx(1),yy(1)), (xx(4),yy(4)), (xx(3),yy(3)))
  a = periodic_sqrt(xx(1), yy(1), xx(2), yy(2))
  b = periodic_sqrt(xx(2), yy(2), xx(3), yy(3))
  c = periodic_sqrt(xx(1), yy(1), xx(3), yy(3))
  s = (a + b + c) / 2
  rectangle_area = rectangle_area + sqrt(s*(s-a)*(s-b)*(s-c))

  a = periodic_sqrt(xx(1), yy(1), xx(4), yy(4))
  b = periodic_sqrt(xx(4), yy(4), xx(3), yy(3))
  c = periodic_sqrt(xx(1), yy(1), xx(3), yy(3))
  s = (a + b + c) / 2
  rectangle_area = rectangle_area + sqrt(s*(s-a)*(s-b)*(s-c))

  if (abs(rectangle_area - sum_triangle_area) < 1e-6) then
     inside = .true.
  end if
end subroutine check_point

! Y. WANG: calculate the distance of two points
! e.g. for two points (179 W, 0 N) and (-179 W, 0N)
!      sqrt = 358**(1/2)
!      periodic_sqrt = 2**(1/2)
!
function periodic_sqrt(lon1,lat1,lon2,lat2)

  real periodic_sqrt
  real, intent (in) :: lon1,lat1,lon2,lat2

  ! local 
  real :: tmp
  
  ! longitude is a periodic function 
  tmp = abs(lon1 - lon2)
  ! smaller difference
  tmp = min(tmp, 360 - tmp)

  periodic_sqrt = sqrt(tmp**2 + (lat1 - lat2)**2)

end function periodic_sqrt

end module m_pivotp_micom
