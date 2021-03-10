      module dimensions
c
c --- itdm  = total grid dimension in i direction
c --- jtdm  = total grid dimension in j direction
c --- kdm   =       grid dimension in k direction
      integer    itdm,jtdm,kdm
      parameter (itdm=320,jtdm=384,kdm=53)
c
c --- ijqr  = total number of active tiles
c --- iqr   = maximum number of tiles in i direction
c --- jqr   = maximum number of tiles in j direction
c --- idm   = maximum single tile grid dimension in i direction
c --- jdm   = maximum single tile grid dimension in j direction
      integer    ijqr,iqr,jqr,idm,jdm
      parameter (ijqr=16,iqr=4,jqr=4,idm=110,jdm=120)
c
c --- region type (nreg=-1: unknown,
c ---              nreg= 0:   closed/closed,
c ---              nreg= 1: periodic/closed,
c ---              nreg= 2: periodic/arctic,
c ---              nreg= 3: periodic/periodic
c ---              nreg= 4:   closed/periodic)
      integer    nreg
      data nreg/1/
c
c --- i0_pe = 1st dimension tile offsets
c --- ii_pe = 1st dimension tile extents
c --- j0_pe = 2nd dimension tile offsets
c --- jj_pe = 2nd dimension tile extents
      integer, dimension(iqr,jqr) :: i0_pe,ii_pe,j0_pe,jj_pe
      data i0_pe/
     .  0,79,170,241,0,85,158,228,0,110,180,243,0,68,176,232/
      data ii_pe/
     .  79,91,71,79,85,73,70,92,110,70,63,77,68,108,56,88/
      data j0_pe/
     .  0,0,0,0,86,86,86,86,174,174,174,174,264,264,264,264/
      data jj_pe/
     .  86,86,86,86,88,88,88,88,90,90,90,90,120,120,120,120/
c
      end module dimensions
