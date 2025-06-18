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
use trans_module, only: trans_t

implicit none

! Includes
#include "setup_trans0.h"
#include "esetup_trans.h"
#include "etrans_inq.h"
#include "edist_spec.h"
#include "einv_trans.h"
#include "egath_grid.h"

private
public :: fieldsio_arome_fa

contains

!----------------------------------------------------------------------

subroutine fieldsio_arome_fa(conf,comm,fspace,trans,akbk,fset)

implicit none

! Passed variables
type(fckit_configuration),intent(in) :: conf
type(fckit_mpi_comm),intent(in) :: comm
type(atlas_functionspace_structuredcolumns),intent(in) :: fspace
type(trans_t),intent(in) :: trans
type(atlas_fieldset),intent(inout) :: akbk
type(atlas_fieldset),intent(inout) :: fset

! Local variables
integer(kind_int),parameter :: ifile = 11
integer(kind_int) :: irep,imaxlev,imaxtrunc,imaxgl,imaxlon,inbari,ityptr,itronc,kflev
integer(kind_int) :: nvar2d,ivar2d,nfield,ifield,nlev,ilev,ingrib,inbits,istron,ipuila
integer(kind_int) :: nlon,ndgl,nmsmax,nsmax,from(1),nproma,ngpblks
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

  ! Get sizes
  nsmax = inozpa(1)
  nmsmax = inozpa(2)
  dx = zsinla(7)
  dy = zsinla(8)

  ! Compare file and transform sizes
  if (nlon /= trans%nlon) call abor1_ftn("inconsistent nlon")
  if (ndgl /= trans%ndgl) call abor1_ftn("inconsistent ndgl")
  if (nsmax /= trans%nsmax) call abor1_ftn("inconsistent nsmax")
  if (nmsmax /= trans%nmsmax) call abor1_ftn("inconsistent nmsmax")

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

! Set nproma/ngpblks
nproma = trans%ngptot
ngpblks = 1

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
  allocate(zgpg(trans%ngptotg,1))
  allocate(zspg(1,trans%nspec2g))
end if
allocate(zgp(trans%ngptot,1,1))
allocate(zsp(1,trans%nspec2))
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
    if (field%shape(2) /= trans%ngptotg) call abor1_ftn("wrong horizontal dimension")

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
                      & kresol = trans%handle, &
                      & pspec = zsp)
      else
        ! Scatter spectral field (receive)
        call edist_spec(kfdistg = 1, &
                      & kfrom = from, &
                      & kresol = trans%handle, &
                      & pspec = zsp)
      end if

      ! Inverse spectral transform
      call einv_trans(kresol = trans%handle, &
                    & kproma = nproma, &
                    & ldscders = .false., &
                    & pspscalar = zsp, &
                    & pgp = zgp)

      if (comm%rank() == 0) then
        ! Gather grid-point field (receive)
        call egath_grid(kresol = trans%handle, &
                      & kfgathg = 1, &
                      & kto = from, &
                      & kproma = nproma, &
                      & pgp = zgp, &
                      & pgpg = zgpg)
      else
        ! Gather grid-point field (send)
        call egath_grid(kresol = trans%handle, &
                      & kfgathg = 1, &
                      & kto = from, &
                      & kproma = nproma, &
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
      do igpg=1,trans%ngptotg
!        iy = (igpg-1)/trans%nlon+1
!        ix = igpg-(iy-1)*trans%nlon
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
