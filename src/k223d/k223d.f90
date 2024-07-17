program k223d
use utils
use typedef
use makepdf
use forparse
implicit none

  type(mesh) :: amesh

  real, parameter :: kilo = 1000.0000, time = 0.
  real, parameter :: pi=acos(-1.), mega = 10.e6
  integer :: m,n,i,na,sct,j,sws,inx,seed,swtape, distance, pdf_variable, seed_stoc, rigidity
  integer :: merc_zone, ii
  real :: dx,mu,rmax,rmin,sd, slip_aux ! moment,mag,

  character (len=20) :: fnamef, string_aux
  character (len=3) :: geo_zone
  character*200 :: index_local
  real :: nuc_ptx,nuc_pty,nuc_ptz  !nucleation location (x,y,z)
  real ::asp_radius,cslp,dummy    !radius of asperity,Esheleby const. , dummy vari able for distance
  integer :: nuc_idx                            ! nodal index for nucleation locat ion
  integer :: num_cells                ! number of nodes per cell
  real :: a,b,Mw,target_area          ! Strasser consts, est. of fault area
  real :: moment
  integer :: findcentre,findgausscentre,numb_gaussians

  integer :: iway,UTM_PROJECTION_ZONE, nuc_mode, coord_type
  logical :: SUPPRESS_UTM_PROJECTION, logical_scenario
  real :: lon,lat
  real :: ran
  type(pdfinputs) :: gausspar ! gaussian parameters
  real, dimension(:), allocatable :: pdf,slipout
  real :: zthr

  type(borderlist), pointer :: aborderlist,surfacelist

  integer,parameter :: nit=24
  integer, parameter :: nit_Mw=25
  integer, parameter :: nit_zone=26

  character*200 :: fname, file_out

  if (parse_arg('input',fname) /= PARSE_OK) stop 'missing input file name'
  !if (parse_arg('index',index_string) /= PARSE_OK) stop 'missing index string'
  !if (parse_arg('numb_gauss',numb_gaussians) /= PARSE_OK) stop 'missing gauss'
  open(nit,file=trim(fname))
  open(nit_Mw,file='input_magnitude')
  open(nit_zone,file='param_zone.dat')
  if (parse_arg('mu',mu,nit) /= PARSE_OK) stop 'rigidity missing/syntax error'
  if (parse_arg('rigidity',rigidity,nit) /= PARSE_OK) stop 'rigidity_logic missing/syntax error'
  if (parse_arg('nb cracks',na,nit) /= PARSE_OK) stop 'crak number missing/syntax error'
  if (parse_arg('rmin',rmin,nit) /= PARSE_OK) stop 'minimum radius missing/syntax error'
  if (parse_arg('rmax',rmax,nit) /= PARSE_OK) stop 'maximum radius missing/syntax error'
  if (parse_arg('seed',seed,nit) /= PARSE_OK) stop 'seed number missing/syntax error'
  if (parse_arg('magnitude',Mw,nit_Mw) /= PARSE_OK) stop 'magnitude missing/syntax error'
  if (parse_arg('z threshold',zthr,nit) /= PARSE_OK) stop 'z threshold missing/syntax error'
  if (parse_arg('nucleation mode',nuc_mode,nit) /= PARSE_OK) stop 'nucleation mode missing/syntax error'
  if (parse_arg('slip taper',swtape,nit) /= PARSE_OK) stop 'slip taper missing/syntax error'
  if (parse_arg('distance',distance,nit) /= PARSE_OK) stop 'distance missing/syntax error'
  if (parse_arg('geo zone',geo_zone,nit_zone) /= PARSE_OK) stop 'Geo zone missing/syntax error'
  if (parse_arg('pdf variable',pdf_variable,nit) /= PARSE_OK) stop ' PDF_variable yes/no missing/syntax error'
  if (parse_arg('seed stochastic',seed_stoc,nit) /= PARSE_OK) stop ' No information about seed /syntax error'
  if (parse_arg('coordinates',coord_type,nit) /= PARSE_OK) stop &
 ' No information about coordinates /syntax error'
  if (parse_arg('mercator',merc_zone,nit_zone) /= PARSE_OK) stop &
 ' No information about mercator zone/syntax error'
  close(nit); close(nit_Mw); close(nit_zone)
  
  amesh%Mw = Mw
  amesh%geo_zone=geo_zone
  amesh%nuc_mode=nuc_mode
  amesh%distance=distance
  amesh%pdf_variable=pdf_variable
  amesh%seed=seed_stoc
  amesh%rigidity=rigidity
  amesh%mu=mu
  amesh%rmax=rmax
  amesh%rmin=rmin
  amesh%coord_type=coord_type
  amesh%merc_zone=merc_zone

sct = 0       ! = 0 => moment scaling, =1 => stress drop scaling (dic)
! relationship based on Strasser et al. 2010  (interface events)
a = -3.476
b = 0.952
target_area = 10**(a+b*Mw)*(kilo**2.)
print*, 'Target area=',target_area, 'Magnitude=',Mw
moment=10**(1.5*Mw+9.1)


write(0,*) "starting ..."

! read input data files
call readfile(amesh) ! reading nodal points and ameshents
call faultborder(amesh,aborderlist)

! set random number
!seed=224
!call srand(seed)
do ii=1,amesh%numb_scen
if (mod(ii,250)==0) write(0,*) "scenario # ", ii
write(*,*) "scenario # ", ii
call set_seed(seed,amesh%seed,amesh%index_string(ii))
!amesh%seed=1
!print*, 'ciao'

if (nuc_mode==1) then
!print*, 'ciao0'
   call select_fault_zone_fromfile(amesh,target_area,ii)
else
!print*, 'ciao1'
   call select_fault_zone(amesh,target_area)
endif

if (amesh%QuakeElemNo<13*amesh%rmin) then    ! generate non-stochastic distributions for small enough ruptures
    index_local=amesh%index_string(ii)
    index_local=index_local(1:5)
    file_out='mu_Slip_aux_' // trim(index_local) // '.dat'
    file_out=trim(file_out)
    open(110,file=file_out)
    read(110,*) slip_aux
    slip_aux=slip_aux/amesh%mu
    close(110)
    amesh%index_string(ii)(7:9)='000'
    !!! WRITE Slip4cells files for linear combinations
    ! file_out='Slip4cells_' // trim(amesh%index_string(ii)) // '.dat'
     !file_out=trim(file_out)
    ! open(12,file=file_out,form='formatted')
    ! if (amesh%Ncells<1000) write(string_aux,"(I3)") amesh%Ncells
    ! if (amesh%Ncells>=1000 .and. amesh%Ncells<10000) write(string_aux,"(I4)") amesh%Ncells
    ! if (amesh%Ncells>=10000 .and. amesh%Ncells<100000) write(string_aux,"(I5)") amesh%Ncells
    ! string_aux='"(' // trim(string_aux) // 'F10.6)"'
    ! do i=1,amesh%Ncells
    ! write(12,"(F10.6)") (slipout(i))
    ! end do
     !1000 format (<amesh%Ncells>F10.6)
    ! call flush(12)
    ! close(12)
    !!! END OF Slip4Cells files writing 
    file_out='Slip4HySea' // trim(amesh%index_string(ii)) // '.dat'
    file_out=trim(file_out)
    open(15,file=file_out,form='formatted')
    write(15,*) "LON1     LAT1    DEPTH1(km)      LON2    LAT2    DEPTH2(km)      LON3    LAT3    DEPTH3(km)      RAKE  SLIP(m)"
    do i=1,amesh%QuakeElemNo
         if (amesh%rigidity==1) then
            allocate(slipout(amesh%Ncells))
            slipout=slip_aux
            call renormalize_slip(amesh,slipout)
            write(15,"(11F12.6)") (amesh%HySea(amesh%QuakeElem(i),j), j=1,10), slipout(amesh%QuakeElem(i))
            deallocate(slipout)
         else
            write(15,"(11F12.6)") (amesh%HySea(amesh%QuakeElem(i),j), j=1,10), slip_aux
         endif
    enddo
    close(15)
    deallocate(amesh%QuakeElem,amesh%QuakeNodes,amesh%QuakeBorder_elem,amesh%QuakeBorder_Nodes)
    cycle
endif

!print*, 'ciao2'
if (swtape == 1) then
   call faultsurface(amesh,aborderlist,surfacelist,zthr)
   call quakesurface(amesh,surfacelist)
endif

call random_number(ran)
! set gaussians that define pdf
i = amesh%numb_gauss(ii) !1+int(ran*3.) ! randomly choose between 1-4 gaussian functions
print*,'number of gaussians.....',i
gausspar%ng=i
allocate(gausspar%vertexcenter(gausspar%ng),gausspar%sizeandhigh(gausspar%ng,2))
!gausspar%vertexcenter(1)=findcentre(amesh)
!gausspar%vertexcenter(1)=findgausscentre(amesh,1.e3)
!gausspar%sizeandhigh(1,1)=sqrt(target_area)/5.
!gausspar%sizeandhigh(1,2)=1. ! relative weight
!print*, amesh%QuakeBorder_Nodes
do i = 1,gausspar%ng
  gausspar%vertexcenter(i)=findgausscentre(amesh,10.e3)   ! gaussian centre can be anywhere bar 10km from edge
  gausspar%sizeandhigh(i,1)=sqrt(target_area)/3.
  gausspar%sizeandhigh(i,2)=1. ! relative weight
enddo
allocate(pdf(amesh%QuakeElemNo))
call faultpdf(gausspar,pdf,amesh,ii)
rmax=amesh%rmax
rmin=amesh%rmin
logical_scenario=.true.
call pdftoslip(pdf,amesh,moment,sd,mu,na,rmin,rmax,sct,logical_scenario)
!if (.not. logical_scenario) then
!        deallocate(amesh%slip,pdf,amesh%QuakeElem,amesh%QuakeNodes,amesh%QuakeBorder_elem,amesh%QuakeBorder_Nodes)
!        deallocate(gausspar%vertexcenter,gausspar%sizeandhigh)
!        cycle
!endif
if (swtape == 2) call trimslip(amesh,zthr)
mu=amesh%mu
call normslip(amesh,moment,mu)

! cut slip as if at surface
!call cutslip(amesh,moment,surf_depth,mu)

write(*,*) minval(amesh%slip),maxval(amesh%slip)

allocate(slipout(amesh%Ncells))
slipout=0.
do i=1,amesh%QuakeElemNo
   slipout(amesh%QuakeElem(i))=amesh%slip(i)
enddo

! output for gmt script
open(10,file="slip.out")
do i=1,amesh%Ncells
  if (slipout(i) < 10) then  ! if slip is between 0-9
   write(10,'(a4,f9.7)') '> -Z',slipout(i)
 else      ! if slip is between 10 - 99
   write(10,'(a4,f10.7)') '> -Z',slipout(i)
 endif
   do j=1,3
      write(10,*) amesh%px(amesh%cell(i,j)),amesh%py(amesh%cell(i,j)),amesh%pz(amesh%cell(i,j))
   enddo
enddo
close(10)

! generate a slip distribution that has a value for all the whole mesh
!deallocate(amesh%slip)
!allocate(amesh%slip(amesh%Ncells))
!amesh%slip = slipout

if (amesh%rigidity==1) then
   print*, "Rigidity variable"
   call renormalize_slip(amesh,slipout)
end if

!open(11,file='Areas.dat',form='formatted')
!do i=1,amesh%Ncells
!write(11,*) amesh%area(i)
!enddo
!call flush(11)
!close(11)

call normarea(amesh,moment,slipout)

deallocate(amesh%slip)
allocate(amesh%slip(amesh%Ncells))
amesh%slip = slipout
!print*, "somma=",sum(amesh%slip)
  !file_out='Slip_' // trim(amesh%index_string(ii)) // '.vtk'
  !file_out=trim(file_out)
  !open(11,file=file_out,form = 'formatted')
  !call dumpmeshvtk(11,amesh)
  !call dumpcellattributevtk(11,amesh,slipout,'slip',.true.)
  !call flush(11)
  !close(11)
 
 ! file_out='Slip4cells_' // trim(amesh%index_string(ii)) // '.dat'
 ! file_out=trim(file_out)
 ! open(12,file=file_out,form='formatted')
 ! if (amesh%Ncells<1000) write(string_aux,"(I3)") amesh%Ncells
 ! if (amesh%Ncells>=1000 .and. amesh%Ncells<10000) write(string_aux,"(I4)") amesh%Ncells
 ! if (amesh%Ncells>=10000 .and. amesh%Ncells<100000) write(string_aux,"(I5)") amesh%Ncells    
 !  string_aux='"(' // trim(string_aux) // 'F10.6)"'
 ! write(12,"(F10.6)") (slipout(i), i=0,amesh%Ncells-1)
  !1000 format (<amesh%Ncells>F10.6)
 ! call flush(12)
 ! close(12)
  file_out='Slip4HySea' // trim(amesh%index_string(ii)) // '.dat'
  file_out=trim(file_out)
  open(15,file=file_out,form='formatted')
  write(15,*) "LON1     LAT1    DEPTH1(km)      LON2    LAT2    DEPTH2(km)      LON3    LAT3    DEPTH3(km)      RAKE    SLIP(m)"
  do i=1,amesh%QuakeElemNo
     write(15,"(11F12.6)") (amesh%HySea(amesh%QuakeElem(i),j), j=1,10), amesh%slip(amesh%QuakeElem(i))
  enddo
 close(15) 
deallocate(amesh%slip,slipout,pdf,amesh%QuakeElem,amesh%QuakeNodes,amesh%QuakeBorder_elem,amesh%QuakeBorder_Nodes)
deallocate(gausspar%vertexcenter,gausspar%sizeandhigh)
enddo
!!! TO WRITE matrix distance binary file
if (amesh%distance==0) then
 file_out= trim(amesh%geo_zone) // '_matrix_distance.bin'
  file_out=trim(file_out)
  open(13,file=file_out,form='unformatted',access='direct',recl=amesh%Nnodes*amesh%Nnodes*4)
!     do i=1,amesh%Nnodes
        write(13,rec=1) amesh%dist   !(i,1:amesh%Nnodes)
!     enddo
  close(13)
endif

!  open(11,file='TohokuLatLon.dat',form = 'formatted')
!      call tsunamiout(11,amesh)
!  close(11)
!   do i=1,amesh%Ncells
!      write(11,*) slipout(i)
!   enddo
!   close(11)

stop
end program k223d

!##########################################################################
subroutine  cutslip(amesh,moment,surf_depth,mu)
use typedef
implicit none

type(mesh) :: amesh
real :: mf,surf_depth,moment,mu
integer :: i,j,k

write(*,*) '-------account for surface effect------------'

do j=1,amesh%QuakeElemNo
   k = amesh%QuakeElem(j)
        if (minval(abs( amesh%pz(amesh%cell(k,:)) )) < surf_depth) then
              amesh%slip(j) = 0.
        endif
enddo

! renormalize
mf = 0.
do i=1,amesh%QuakeElemNo
   mf=mf+amesh%slip(i)*amesh%area(amesh%QuakeElem(i))
enddo
mf=mf*mu
print*, "mf=",mf
write(*,*) "difference to target : ",100*(mf-moment)/moment,"%"
do i=1,amesh%QuakeElemNo
      amesh%slip(i)=amesh%slip(i)/mf*moment
enddo

end subroutine cutslip
!##########################################################################
function findgausscentre(amesh,dist)
! pick gaussian centre with the condition that it cannot be within distance 'dist'
! from boundary edge
  use typedef
    type(mesh) :: amesh
    integer :: findgausscentre
    real :: dist
    real :: x,vmin,vmax
    integer :: i,j,inx,id,idx
    integer,allocatable,dimension(:) :: keep_id
    logical :: good2use
    allocate(keep_id(amesh%QuakeNodesNo))
    inx = 0
     do i=1,amesh%QuakeNodesNo
       good2use = .true.
       do j=1,amesh%QuakeBorder_NodesNo
        if (amesh%QuakeBorder_Nodes(j)>0 .and. amesh%QuakeBorder_Nodes(j)<=amesh%Nnodes) then
         if (amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j)) < dist) then
              good2use = .false.
              exit
         endif
        endif
       enddo
        if (good2use) then
          inx = inx+1
          keep_id(inx) = i
        endif
    enddo

    if (inx == 0) then
       inx=amesh%QuakeNodesNo
       keep_id=(/ (k, k=1,amesh%QuakeNodesNo)/)
    endif

    call random_number(x)
    id = floor(x*float(inx))+1
    findgausscentre = keep_id(id)
!
   ! do i=1,amesh%QuakeNodesNo
    !   if (amesh%QuakeNodes(i) == 66) then
     !     findgausscentre = i
      !    exit
      ! endif
   ! enddo
    write(*,*) "findgausscentre :",findgausscentre,amesh%QuakeNodes(findgausscentre)

end function findgausscentre
!##########################################################################
function findcentre(amesh)
use typedef
  type(mesh) :: amesh
  integer :: findcentre

  integer :: i,j
  real :: sd,mean,sdmin
  real, parameter :: infini=1.e32

  findcentre=0
  sdmin=infini
  do i=1,amesh%QuakeNodesNo
     mean=0.
     do j=1,amesh%QuakeBorder_NodesNo
        mean=mean+amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j))
     enddo
     mean=mean/float(amesh%QuakeBorder_NodesNo)
     sd=0.
     do j=1,amesh%QuakeBorder_NodesNo
        sd=sd+(amesh%dist(amesh%QuakeNodes(i),amesh%QuakeBorder_Nodes(j))-mean)**2.
     enddo
     if (sd < sdmin) then
        sdmin=sd
        findcentre=i
     endif
  enddo
  return
end function findcentre
!##########################################################################
