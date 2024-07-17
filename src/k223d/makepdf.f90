module makepdf
use lateration
implicit none

type pdfinputs
   integer :: ng
! pointer to a Node number of the Quake
   integer, dimension(:), allocatable :: vertexcenter
   real, dimension(:,:), allocatable :: sizeandhigh
end type pdfinputs

contains

!#######################################################################################
subroutine faultpdf(gausspar,pdf,amesh,ii)
implicit none

  integer, intent(in) :: ii 
! gaussian parameters
  type(pdfinputs) :: gausspar
! main structure
  type(mesh) :: amesh
! pdf on cells
  real, dimension(amesh%QuakeElemNo) :: pdf
  real, dimension(amesh%Ncells) :: pdfout

  integer :: i,j,k
  real, dimension(3) :: pdfvertexofcurrentcell
  real :: d,ug,pdfint,totarea,pdf_aux
  character*200 :: file_out, index_local


  totarea=0.
  do i=1,amesh%QuakeElemNo ! loop on all the cells
     do j=1,3
        pdfvertexofcurrentcell(j)=0.
     enddo
     do j=1,gausspar%ng ! loop on the gaussian
        do k=1,3 ! loop on the vertex of the cell
          d=amesh%dist(amesh%QuakeNodes(gausspar%vertexcenter(j)),amesh%cell(amesh%QuakeElem(i),k))
!          d=amesh%dist(amesh%QuakeNodes(gausspar%vertexcenter(j)),amesh%cell(amesh%QuakeElem(i),k))
!          d=amesh%dist(amesh%QuakeElem(gausspar%vertexcenter(j)),amesh%cell(amesh%QuakeElem(i),k))
!          d=amesh%dist(gausspar%vertexcenter(j),amesh%cell(amesh%QuakeElem(i),k))
           ug=exp(-d**2./2./gausspar%sizeandhigh(j,1)**2.)
           pdfvertexofcurrentcell(k)=pdfvertexofcurrentcell(k)+ug*gausspar%sizeandhigh(j,2)
!           pdfvertexofcurrentcell(k)=pdfvertexofcurrentcell(k)+d
        enddo
     enddo

     pdf(i)=0.
     do j=1,3
        pdf(i)=pdf(i)+pdfvertexofcurrentcell(j)
     enddo
!     pdf(i)=pdf(i)/3.*amesh%area(i)
     pdf(i)=pdf(i)/3.!*amesh%area(amesh%QuakeElem(i))
!     pdfint=pdfint+pdf(i)*elem%area(i)
     totarea=totarea+amesh%area(amesh%QuakeElem(i))
!     totarea=totarea+amesh%area(i)
  enddo
  pdfint=0.
   if (amesh%pdf_variable==1) then
      index_local=amesh%index_string(ii)
      index_local=index_local(1:5)
      file_out='Slip_PDF_' // trim(index_local) // '.dat'
      file_out=trim(file_out)
      open(23,file=file_out,action='read')
      do i=1,amesh%QuakeElemNo
         read(23,*) pdf_aux
         pdf(i)=(pdf(i)*pdf_aux)/totarea
         pdfint=pdfint+pdf(i)
      enddo
      close(23)
   else
     do i=1,amesh%QuakeElemNo
         pdf(i)=pdf(i)/totarea
         pdfint=pdfint+pdf(i)
      enddo
   endif

  ug=0.
  do i=1,amesh%QuakeElemNo
     pdf(i)=pdf(i)/pdfint
     ug=ug+pdf(i)
  enddo
  pdfout=0.
  do i=1,amesh%QuakeElemNo
     pdfout(amesh%QuakeElem(i))=pdf(i)
  enddo
  file_out='PDF' // trim(amesh%index_string(ii)) // '.vtk'
  file_out=trim(file_out)

  !open(22,file=file_out)
  !call dumpmeshvtk(22,amesh)
  !call dumpcellattributevtk(22,amesh,pdfout,'pdf',.true.)
  !close(22)
end subroutine faultpdf
!#######################################################################################
subroutine pdftoslip(pdf,amesh,moment,sd,mu,na,rmin,rmax,sct,logic_scenario)
  use typedef
! main structure
  type(mesh) :: amesh
  logical, intent (inout) :: logic_scenario
! pdf on cells
  real, dimension(amesh%QuakeElemNo) :: pdf
  integer :: na,sct
  real :: moment,sd,mu,rmin,rmax
  real :: random
  real :: cslp,p,length,width,area,dx,dmean,mf,db,dummy,numrnd,curcum
  real, dimension(:), allocatable :: r
  real, dimension(3) :: smean
  real, parameter :: pi=acos(-1.)
  integer :: i,j,k,idxmax,cellidx,casp,cnt

! scaling operation between moment and stress drop
  call wl(amesh,length,width,area,dx)
  print*, moment, length, area
  if (sct==1) then
     moment=sd*length*area
  else
     sd=moment/(length*area)
  endif
! Eshelby's constant
  cslp=24/7/pi*sd/mu
! min/max of the radius distribution
! TB BE DEFINED MORE CLEARLY <=====================================
  rmax=width*rmax
  write(*,*) "width and rmax :",width,rmax
  rmin=dx*rmin
  print*, "rmin=", rmin
! fractal parameter
  p=2.*7./16.*moment/sd/(rmax-rmin)
! memory allocation
  allocate(r(na))
! fractal distribution of the radii (Zeng el al. 1994)
  do i=1,na
! D=2
    call random_number(random)
    r(i)=(2.*random*float(na)/p+rmax**(-2.))**(-.5)
!    r(i)=(2.*rand(0)*float(na)/p+rmax**(-2.))**(-.5)
  enddo
! Sort of the distribution by size
  call reorder(r,na)
! locating the maximum of the pdf
  idxmax=maxloc(pdf,1)
! initializing the slip array
  allocate(amesh%slip(amesh%QuakeElemNo))
  amesh%slip=0.
! asperity loop
  do i=1,na
!  write(0,*) "working on asperity #",i,"/",na
! select the cell of the asperity center
! and one of its nodes
     if (i == 1) then
        cellidx=idxmax
        call random_number(random)
        casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(random*3+1),3))
        write(*,*) 'depth of max pdf..',amesh%pz(amesh%cell(amesh%QuakeElem(cellidx),:))
!        casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(rand(0)*3+1),3))
     else
        cnt=0
        db=r(i)-1.
        do while (r(i) > db)
           cnt=cnt+1
           if (cnt == 100 ) then
              write(*,*) "reach 100"
              logic_scenario=.false.
              goto 1000
           endif
           curcum=0.
           call random_number(numrnd)
!           numrnd=rand(0)
           cellidx=0
           do while (curcum <= numrnd)
              cellidx=cellidx+1
              if (cellidx == amesh%QuakeElemNo) exit
              curcum=curcum+pdf(cellidx)
           enddo
           call random_number(random)
           casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(random*3+1),3))
!           casp=amesh%cell(amesh%QuakeElem(cellidx),min(int(rand(0)*3+1),3))
           db=mindistbd(amesh,casp) !
        enddo
     endif
! adding the asperity to the slip distribution
     do j=1,amesh%QuakeElemNo
        dmean=0.
        do k=1,3
           dmean=dmean+amesh%dist(casp,amesh%cell(amesh%QuakeElem(j),k))
        enddo
        dmean=dmean/3.
        if (dmean > r(i) ) cycle
        amesh%slip(j)=amesh%slip(j)+cslp*sqrt(r(i)**2.-dmean**2.)

!        smean = 0.
!        do k=1,3
!          if ( amesh%dist(casp,amesh%cell(amesh%QuakeElem(j),k)) < r(i) ) then
!            smean(k) = amesh%slip(j)+cslp*sqrt(r(i)**2.-amesh%dist(casp,amesh%cell(amesh%QuakeElem(j),k))**2.)
!          endif
!        enddo
!        if ( minval(smean(:)) == 0. ) cycle
!        amesh%slip(j)=sum(smean)/3.


     enddo
  enddo
  1000 continue
  if (.not.logic_scenario) amesh%slip=pdf
  !return
end subroutine pdftoslip
!#######################################################################################
subroutine normslip(amesh,moment,mu)
! slip renormalization
use typedef
! main structure
  type(mesh) :: amesh
  real :: moment,mu

  real :: mf
  integer :: i

  mf=0.
  do i=1,amesh%QuakeElemNo
     mf=mf+amesh%slip(i)*amesh%area(amesh%QuakeElem(i))
  enddo
  mf=mf*mu
  write(*,*) "difference to target : ",100*(mf-moment)/moment,"%"
  do i=1,amesh%QuakeElemNo
     amesh%slip(i)=amesh%slip(i)/mf*moment
  enddo
  return
end subroutine normslip
!#######################################################################################
function mindistbd(amesh,anode)

  real :: mindistbd
! main structure
  type(mesh) :: amesh
  integer :: anode

  real, parameter :: infini=1.e32
  integer :: i

  mindistbd=infini
  do i=1,amesh%QuakeBorder_NodesNo
     if (amesh%QuakeBorder_Nodes(i)>0 .and. amesh%QuakeBorder_Nodes(i)<=amesh%Nnodes) &
     mindistbd=min(mindistbd,amesh%dist(amesh%QuakeBorder_Nodes(i),anode))
  enddo

  return
end function mindistbd
!#######################################################################################
subroutine wl(amesh,length,width,area,dx)

! main structure
  type(mesh) :: amesh
  real :: length,width,area,dx

  integer :: i,j
  real, parameter :: kilo=1.e3

! fault length approximation
  length=0.
  do i=1,amesh%QuakeBorder_NodesNo
   if (amesh%QuakeBorder_Nodes(i)>0 .and. amesh%QuakeBorder_Nodes(i)<=amesh%Nnodes) then
     do j=i,amesh%QuakeBorder_NodesNo
       if (amesh%QuakeBorder_Nodes(j)>0 .and. amesh%QuakeBorder_Nodes(j)<=amesh%Nnodes) then
        length=max(length,amesh%dist(amesh%QuakeBorder_Nodes(i),amesh%QuakeBorder_Nodes(j)))
       endif
     enddo
   endif
  enddo
  write(*,*) "length :",length
! area of the fault
  area=0.
  do i=1,amesh%QuakeElemNo
     area=area+amesh%area(amesh%QuakeElem(i))
  enddo
  write(*,*) "area :",area
! fault width approximation
  width=area/length
  write(*,*) "width :",width
! dx approximation
  dx=sqrt(area/float(amesh%QuakeElemNo))
  return
end subroutine wl
!#######################################################################################
subroutine trimslip(amesh,zthr)
use typedef
! main structure
  type(mesh) :: amesh
  real :: zthr
 
  integer :: i

  do i=1,amesh%QuakeElemNo
     if (minval(abs(amesh%pz(amesh%cell(amesh%QuakeElem(i),:)))) < zthr) then
        amesh%slip(i)=0.
     endif
  enddo
return
end subroutine trimslip
!#######################################################################################
end module makepdf
