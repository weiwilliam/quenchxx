!----------------------------------------------------------------------
! Module: fieldsio_arome_fa_mod
!> AROME FA reader module
! Author: Benjamin Menetrier
! Copyright 2025 Meteorologisk Institutt
!----------------------------------------------------------------------
module fieldsio_arome_fa_mod

use atlas_module, only: atlas_functionspace_structuredcolumns,atlas_fieldset,atlas_field,atlas_structuredgrid,atlas_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use fckit_log_module, only: fckit_log
use kinds, only: kind_int, kind_real
use mpl_module, only: mpl_init

implicit none

! Includes
#include "setup_trans0.h"
#include "esetup_trans.h"
#include "etrans_inq.h"
#include "edist_spec.h"
#include "einv_trans.h"
#include "egath_grid.h"

! Constants
integer(kind_int),parameter :: trans_max_handles = 100

! Trans object
type type_trans
  integer(kind_int) :: handle
  integer(kind_int) :: ndgl
  integer(kind_int) :: nlon
  integer(kind_int) :: nmsmax
  integer(kind_int) :: nsmax
  integer(kind_int) :: nspec2
  integer(kind_int) :: nspec2g
  integer(kind_int) :: ngptot
  integer(kind_int) :: ngptotg
  integer(kind_int) :: nproma
  integer(kind_int) :: ngpblks
  real(kind_real) :: dx
  real(kind_real) :: dy
end type type_trans

! Handles counter
integer(kind_int) :: trans_count_handles = 0

! Hangles
type(type_trans),dimension(trans_max_handles) :: trans

private
public :: fieldsio_arome_fa

contains

!----------------------------------------------------------------------

subroutine fieldsio_arome_fa(conf,comm,fspace,akbk,fset)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_mpi_comm),intent(in) :: comm
type(atlas_functionspace_structuredcolumns),intent(in) :: fspace
type(atlas_fieldset),intent(inout) :: akbk
type(atlas_fieldset),intent(inout) :: fset

! Local variables
integer(kind_int),parameter :: ifile = 11
integer(kind_int) :: irep,imaxlev,imaxtrunc,imaxgl,imaxlon,inbari,ityptr,itronc,kflev
integer(kind_int) :: nvar2d,ivar2d,nfield,ifield,nlev,ilev,ingrib,inbits,istron,ipuila
integer(kind_int) :: nlon,ndgl,nmsmax,nsmax,itrans,from(1)
integer(kind_int) :: nprgpew,nprtrv,nprtrw,nprgpns,n_regions_ns,n_regions_ew
integer(kind_int) :: igpg,ix,iy,inode
integer(kind_int),allocatable :: inlopa(:),inozpa(:),nloen(:),levvec(:)
integer(kind_int),allocatable :: i_regions(:)
real(kind_real) :: dx,dy,zslapo,zclopo,zslopo,zcodil,zref,zeps,zundf,req
real(kind_real),allocatable :: zsinla(:),zvalh(:),zvbh(:)
real(kind_real),allocatable :: zgpg(:,:),zspg(:,:),zgp(:,:,:),zsp(:,:)
real(kind_real),pointer :: ak_ptr(:),bk_ptr(:),ptr(:,:)
character(len=256) :: clfile
character(len=16) :: clframe
character(len=1024) :: message,varname
character(len=1024),allocatable :: prevec(:),varvec(:)
character(len=:),allocatable :: str,str_array(:)
logical :: lgard,found,lexist,lcosp,lundf
!type(atlas_structuredgrid) :: grid
type(atlas_field) :: ak,bk,field

if (trans_count_handles == 0) then
  ! Setup parallelization
  nprgpew = max(1,int(sqrt(real(comm%size(),kind_real)),kind_int))
  call mpl_init(koutput=0,kunit=6,ldinfo=.false.)
  allocate(i_regions(comm%size()))
  nprtrv = 1
  nprgpns = comm%size()/nprgpew
  nprtrw = comm%size()/nprtrv;
  call conf%get_or_die("earth radius",req)
  call setup_trans0(kout = 99, &
                  & kerr = 99, &
                  & kprintlev = 0, &
                  & kmax_resol = trans_max_handles, &
                  & kprtrw = nprtrw, &
                  & ldeq_regions = .false., &
                  & kprgpns = nprgpns, &
                  & kprgpew = nprgpew, &
                  & prad = req, &
                  & k_regions_ns = n_regions_ns, &
                  & k_regions_ew = n_regions_ew, &
                  & k_regions = i_regions, &
                  & ldmpoff = .false. )
  deallocate(i_regions)
  if (comm%rank() == 0) then
    write(message,'(a)') "Info     : AROME FA reader: parallelization setup done (only once per execution)"
    call fckit_log%info(message)
  end if
end if

if (comm%rank() == 0) then
  ! Open file
  call conf%get_or_die("filepath",str)
  clfile = str
  WRITE(clframe,'(''CADRE_LECTURE'',I3.3)') 1
  call faitou(irep,ifile,.true.,clfile,'OLD',.true.,.true.,0,1,inbari,clframe)
  call lfimst(irep,ifile,.false.)

  ! Get file limits
  call falimu(imaxlev,imaxtrunc,imaxgl,imaxlon)

  ! Allocation
  allocate(inlopa(imaxgl))
  allocate(inozpa(imaxgl))
  allocate(zsinla(imaxgl))
  allocate(zvalh(0:imaxlev))
  allocate(zvbh(0:imaxlev))

  ! Read file characteristics
  kflev=imaxlev
  lgard=.false.
  call facies(clframe,ityptr,zslapo,zclopo,zslopo,zcodil,itronc,&
   & ndgl,nlon,inlopa,inozpa,zsinla,kflev,zref,zvalh,&
   & zvbh,lgard)

  ! Check number of levels
  if (kflev > imaxlev) then
    call abor1_ftn('max. number of level in model too small')
  endif

  ! Test file format
  if(zsinla(1) >= 0.0_kind_real) then
    call abor1_ftn("old eggx frame format")
  end if

  ! Copy into trans object
  nmsmax = inozpa(1)
  nsmax = inozpa(2)
  dx = zsinla(7)
  dy = zsinla(8)

  ! Copy ak/bk
  ak = atlas_field("ak",atlas_real(kind_real),(/kflev+1/))
  bk = atlas_field("bk",atlas_real(kind_real),(/kflev+1/))
  call akbk%add(ak)
  call akbk%add(bk)
  call ak%data(ak_ptr)
  call bk%data(bk_ptr)
  do ilev=0,kflev
    ak_ptr(ilev+1) = zvalh(ilev)
    bk_ptr(ilev+1) = zvbh(ilev)
  end do
end if

! Broadcast sizes
call comm%broadcast(nlon,0)
call comm%broadcast(ndgl,0)
call comm%broadcast(nmsmax,0)
call comm%broadcast(nsmax,0)
call comm%broadcast(dx,0)
call comm%broadcast(dy,0)

! Find if the transform already exists
found = .false.
itrans = 0
do while (.not.found)
  ! Update index
  itrans = itrans + 1

  ! Check index
  if (itrans > trans_max_handles) call abor1_ftn("too many trans handles")
  if (itrans > trans_count_handles) then
    ! Copy sizes in handle
    trans(itrans)%nlon = nlon
    trans(itrans)%ndgl = ndgl
    trans(itrans)%nmsmax = nmsmax
    trans(itrans)%nsmax = nsmax
    trans(itrans)%dx = dx
    trans(itrans)%dy = dy

    ! Setup new transform
    trans(itrans)%handle = itrans
    allocate(nloen(trans(itrans)%ndgl))
    nloen = trans(itrans)%nlon
    call esetup_trans(kmsmax = trans(itrans)%nmsmax, &
                    & ksmax = trans(itrans)%nsmax, &
                    & kdgl = trans(itrans)%ndgl, &
                    & kdgux = trans(itrans)%ndgl, &
                    & kloen = nloen, &
                    & ldsplit = .true., &
                    & kresol = trans(itrans)%handle, &
                    & pexwn = trans(itrans)%dx, &
                    & peywn = trans(itrans)%dy, &
                    & ldgridonly = .false.) ! could be true if dist_grid only, no transform
    deallocate(nloen)

    ! Get new transform info
    call etrans_inq(kresol = trans(itrans)%handle, &
                  & kspec2 = trans(itrans)%nspec2, &
                  & kspec2g = trans(itrans)%nspec2g, &
                  & kgptot = trans(itrans)%ngptot, &
                  & kgptotg = trans(itrans)%ngptotg)

    ! Set nproma/ngpblks
    trans(itrans)%nproma = trans(itrans)%ngptot
    trans(itrans)%ngpblks = 1

    ! New handle done
    trans_count_handles = trans_count_handles + 1
    found = .true.
  else
    ! Check nlon/nlat/nmsmax/nsmax/dx/dy
    if ((nlon == trans(itrans)%nlon).and.(ndgl == trans(itrans)%ndgl) &
   .and.(nmsmax == trans(itrans)%nmsmax).and.(nsmax == trans(itrans)%nsmax) &
 & .and.(.not.(abs(dx - trans(itrans)%dx) > 0.0)).and.(.not.(abs(dy - trans(itrans)%dy) > 0.0))) then
      ! Found a valid handle
      found = .true.
    end if
  end if
end do

if (comm%rank() == 0) then
  ! Get variables to read
  call conf%get_or_die("nvar2d",nvar2d)
  allocate(prevec(nvar2d))
  allocate(levvec(nvar2d))
  allocate(varvec(nvar2d))
  call conf%get_or_die("prefix vector",str_array)
  prevec = str_array
  call conf%get_or_die("level vector",levvec)
  call conf%get_or_die("variable vector",str_array)
  varvec = str_array

  ! Allocation
  allocate(zgpg(trans(itrans)%ngptotg,1))
  allocate(zspg(1,trans(itrans)%nspec2g))
end if
allocate(zgp(trans(itrans)%ngptot,1,1))
allocate(zsp(1,trans(itrans)%nspec2))
from = 1

! Get ATLAS grid
!grid = fspace%grid()

! Get number of fields
if (comm%rank() == 0) nfield = fset%size()
call comm%broadcast(nfield,0)

! Loop over fields
ivar2d = 0
do ifield=1,nfield
  if (comm%rank() == 0) then
    ! Get field
    field = fset%field(ifield)

    ! Check horizontal dimension
    if (field%shape(2) /= trans(itrans)%ngptotg) call abor1_ftn("wrong horizontal dimension")

    ! Get number of levels
    nlev = field%levels()
  end if
  call comm%broadcast(nlev,0)

  ! Loop over levels
  do ilev=1,nlev
    ! Update variable/level index
    ivar2d = ivar2d+1

    ! Get 2D field info
    if (comm%rank() == 0) call  fanion(irep,ifile,trim(prevec(ivar2d)),levvec(ivar2d),trim(varvec(ivar2d)), &
     & lexist,lcosp,ingrib,inbits,istron,ipuila)

    ! Broadcast flags
    call comm%broadcast(lexist,0)
    call comm%broadcast(lcosp,0)

    ! Check 2D field existence
    if (.not.lexist) then
      call abor1_ftn("field does not exist")
    end if

    ! Check 2D field storage (grid-point or spectral
    if (lcosp) then
      ! Spectral field
      if (comm%rank() == 0) then
        ! Read field
        call facilo(irep,ifile,trim(prevec(ivar2d)),levvec(ivar2d),trim(varvec(ivar2d)),zspg(:,1),lcosp,lundf,zundf)

        ! Scatter spectral field (send)
        call edist_spec(pspecg = zspg, &
                      & kfdistg = 1, &
                      & kfrom = from, &
                      & kresol = trans(itrans)%handle, &
                      & pspec = zsp)
      else
        ! Scatter spectral field (receive)
        call edist_spec(kfdistg = 1, &
                      & kfrom = from, &
                      & kresol = trans(itrans)%handle, &
                      & pspec = zsp)
      end if

      ! Inverse spectral transform
      call einv_trans(kresol = trans(itrans)%handle, &
                    & kproma = trans(itrans)%nproma, &
                    & ldscders = .false., &
                    & pspscalar = zsp, &
                    & pgp = zgp)

      if (comm%rank() == 0) then
        ! Gather grid-point field (receive)
        call egath_grid(kresol = trans(itrans)%handle, &
                      & kfgathg = 1, &
                      & kto = from, &
                      & kproma = trans(itrans)%nproma, &
                      & pgp = zgp, &
                      & pgpg = zgpg)
      else
        ! Gather grid-point field (send)
        call egath_grid(kresol = trans(itrans)%handle, &
                      & kfgathg = 1, &
                      & kto = from, &
                      & kproma = trans(itrans)%nproma, &
                      & pgp = zgp)
      end if
    else
      ! Read grid-point field
      if (comm%rank() == 0) call facilo(irep,ifile,trim(prevec(ivar2d)),levvec(ivar2d),trim(varvec(ivar2d)),zgpg(:,1), &
       & lcosp,lundf,zundf)
    end if

    ! Copy data
    if (comm%rank() == 0) then
      call field%data(ptr)
      do igpg=1,trans(itrans)%ngptotg
!        iy = (igpg-1)/trans(itrans)%nlon+1
!        ix = igpg-(iy-1)*trans(itrans)%nlon
!        inode = grid%index(ix,iy)
!        ptr(ilev,inode) = zgpg(igpg,1)
        ptr(ilev,igpg) = zgpg(igpg,1)
      end do
    end if
  end do
end do

if (comm%rank() == 0) then
  ! Release memory
  deallocate(inlopa)
  deallocate(inozpa)
  deallocate(zsinla)
  deallocate(zvalh)
  deallocate(zvbh)
  deallocate(prevec)
  deallocate(levvec)
  deallocate(varvec)
  deallocate(zgpg)
  deallocate(zspg)

  ! Close file
  call fairme(irep,ifile,'UNKNOWN')
end if

! Release memory
deallocate(zgp)
deallocate(zsp)

call comm%barrier()

end subroutine fieldsio_arome_fa

!----------------------------------------------------------------------

end module fieldsio_arome_fa_mod
