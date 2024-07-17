module lateration
implicit none
type mesh
       integer :: Nnodes,Ncells, nuc_mode, distance, pdf_variable, seed, rigidity, coord_type ! Number of Nodes and number of ameshents
       integer :: merc_zone, numb_scen
       integer :: QuakeElemNo,QuakeNodesNo ! Number of ameshents and nodes for an earthquake
       real,allocatable,dimension(:,:) :: dist, HySea ! distance and HySea matrix 
       integer,allocatable,dimension(:,:) ::  cell
       integer,allocatable,dimension(:,:) ::  EToE,EToF  !EtoV ! elements to : nodes, neighbours, faces
       integer,allocatable,dimension(:) :: QuakeElem,QuakeNodes ! array of ameshents and nodes for an earthquake
       integer :: NoNucLocs, NoPropLocs  ! number of nucleation locations and propagation
       integer,allocatable,dimension(:) :: NucLocs, PropLocs, numb_gauss ! arry of nodal ids where earthquake can start and propagate
       integer :: QuakeBorder_NodesNo, QuakeBorder_ElemNo ! Number of ameshents and nodes forming the earthquake border
       integer,allocatable,dimension(:) :: QuakeBorder_Nodes, QuakeBorder_Elem ! ameshents and nodes forming the earthquake border
       real,allocatable,dimension(:) :: px,py,pz,lon,lat,area,slip ! Coordinates of the nodes, area array for the ameshents and slip arra for the ameshents
       real :: Mw, mu, rmax, rmin
       character (len=3) :: geo_zone
       character (len=9),allocatable, dimension(:) :: index_string
end type mesh
!  type mesh
!    integer :: Nnodes,Ncells
!    real, dimension(:), allocatable :: px,py,pz
!    integer, dimension(:,:), allocatable :: cell
!  end type mesh

  type listn
    integer :: idnode
    type(listn), pointer :: previous
    type(listn), pointer :: next
  end type listn

  type container
    type(listn), pointer :: ptr
  end type container

  real, parameter :: infinity=1.e32
  logical, parameter :: verbose=.false.

contains

!###############################################################################
subroutine allvsall2d(amesh)
  type(mesh) :: amesh
!  real, dimension(:,:), allocatable :: distarray

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last
  real :: d01,d02,d12,r1,r2,dface,dedge,dtest
  integer, dimension(2) :: otherninc
  integer :: i,j,k,toggle
  logical :: begin,hasbeenupdated

! intialisation of dist array
!     allocation
  allocate(amesh%dist(amesh%Nnodes,amesh%Nnodes))
!     set to infinity
  amesh%dist=infinity
!     trace to zero
  do i=1,amesh%Nnodes
     amesh%dist(i,i)=0.
  enddo
! computing distance to neighbor nodes cell by cell
  do i=1,amesh%Ncells
!     distance 1-2
      call updclosedist(amesh,amesh%dist,i,1,2)
!     distance 2-3
      call updclosedist(amesh,amesh%dist,i,2,3)
!     distance 3-1
      call updclosedist(amesh,amesh%dist,i,3,1)
  enddo
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
! main loop on nodes - nodes is a starting point
  if (verbose) write(*,*) 'starting main loop on nodes'
  do k=1,amesh%Nnodes
  if (verbose) write(*,*) '##############################################################'
  if (verbose) write(*,*) 'distance from node #',k,'/',amesh%Nnodes
     call printperc(k,amesh%Nnodes)
! initializing  node todo list
     if (verbose) write(*,*) 'entering vicinode'
     call vicinode(amesh,k,ntoc(k)%ptr,ntodo)
! loop on the to do list
     if (verbose) call printlist(ntodo)
     if (verbose) write(*,*) 'entering in the ntodo list management loop'
     do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
        if (verbose) write(*,*) 'propagating distance from node #',ntodo%idnode
        if (verbose) write(*,*) 'its own distance is :',amesh%dist(k,ntodo%idnode)
        hasbeenupdated=.true.
        toggle=2
        do while(hasbeenupdated)
        hasbeenupdated=.false.
! switch the toggle
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'sweep forward on cellist'
        if (toggle==2) write(*,*) 'sweep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(ntodo%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,ntodo)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,4(i2.2))') 'code',k,ntodo%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!      symmetry check
                if (verbose) write(*,*) '       symmetry check:'
                if (verbose) write(*,*) otherninc(i),'<',k,' ???'
                if (otherninc(i) < k) then
                   if (verbose) write(*,*) '       yes'
                   if (verbose) write(*,*) 'is dist(k,otherninc(i)) > dist(otherninc(i),k) ?'
                   if (verbose) write(*,*) amesh%dist(k,otherninc(i)),' > ',amesh%dist(otherninc(i),k)
                   if (amesh%dist(k,otherninc(i)) > amesh%dist(otherninc(i),k)) then
                      if (verbose) write(*,*) '       yes then take the complimentary'
                      amesh%dist(k,otherninc(i))=amesh%dist(otherninc(i),k)
                      call updatelist(otherninc(i),ntodo)
                   endif
                   cycle
                endif
!     edge propagation
               dedge=amesh%dist(k,ntodo%idnode)+amesh%dist(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',amesh%dist(k,ntodo%idnode),'+',amesh%dist(ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dedge
!     face propagation
               if (amesh%dist(k,otherninc(3-i)) /= infinity) then
               if (verbose) write(*,*) 'computing dface because ',amesh%dist(k,otherninc(3-i)),' is not infinity'
                  d12=amesh%dist(ntodo%idnode,otherninc(3-i))
                  d01=amesh%dist(otherninc(i),ntodo%idnode)
                  d02=amesh%dist(otherninc(i),otherninc(3-i))
               if (verbose) write(*,*) 'd12,d01,d02 :'
               if (verbose) write(*,*) d12,d01,d02
                  r1=amesh%dist(k,ntodo%idnode)
                  r2=amesh%dist(k,otherninc(3-i))
               if (verbose) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d01,d02,r1,r2)
               if (verbose) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
!               dtest=dedge
               if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose) write(*,*) 'dtest < distance between',k,' and ',otherninc(i)
               if (verbose) write(*,*) dtest,amesh%dist(k,otherninc(i))
               if (dtest < amesh%dist(k,otherninc(i))) then
                   if (verbose) write(*,*) 'better !'
! distance is better : if not in the list, add it
                  amesh%dist(k,otherninc(i))=dtest
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
! cancelling the cell vicinity list
! removing the first element of the to do list
        if (associated(ntodo%next)) then
           ntodo=>ntodo%next
           deallocate(ntodo%previous)
           nullify(ntodo%previous)
        else
           deallocate(ntodo)
           nullify(ntodo)
        endif
      enddo
  enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
end subroutine allvsall2d
!###############################################################################
subroutine onevsall2d(amesh,k,dist)

  type(mesh) :: amesh
  real, dimension(:), allocatable :: dist
  integer :: k

  type(container), dimension(amesh%Nnodes) :: ntoc
  type(listn), pointer :: ntodo,cellcur,pcur,last
  real :: d13,d23,d12,r1,r2,dface,dedge,dtest
  integer, dimension(2) :: otherninc
  integer :: i,j,toggle
  logical :: begin,hasbeenupdated
  integer :: nswp,mxswp

  mxswp=0
! intialisation of dist array
!     allocation
  allocate(dist(amesh%Nnodes))
!     set to infinity
  dist=infinity
!     k nodes set to zero
  dist(k)=0.
! computing the node-to-cell array ntoc
  call compntoc(amesh,ntoc)
  if (verbose) write(*,*) 'ntoc(500)%ptr%idnode',ntoc(500)%ptr%idnode
  if (verbose) write(*,*) '##############################################################'
  if (verbose) write(*,*) 'distance from node #',k,'/',amesh%Nnodes
! initializing  node todo list
  if (verbose) write(*,*) 'entering vicinode'
  if (verbose) write(*,*) 'ntoc(k)%ptr%idnode',ntoc(k)%ptr%idnode
  if (verbose) call printlist(ntoc(k)%ptr)
  call vicinode(amesh,k,ntoc(k)%ptr,ntodo)
  call startdist(amesh,k,ntodo,dist)
! loop on the to do list
  if (verbose) call printlist(ntodo)
  if (verbose) write(*,*) 'entering in the ntodo list management loop'
  do while (associated(ntodo))
! searching the cells attached to the first element of the to do list
        if (verbose) write(*,*) 'state of ntodo :'
        if (verbose) call printlist(ntodo)
        if (verbose) write(*,*) 'propagating distance from node #',ntodo%idnode
        if (verbose) write(*,*) 'its own distance is :',dist(ntodo%idnode)
        hasbeenupdated=.true.
        toggle=2
        nswp=0
        do while (hasbeenupdated)
        nswp=nswp+1
        hasbeenupdated=.false.
        toggle=3-toggle
        if (verbose) then
        if (toggle==1) write(*,*) 'seep forward on cellist'
        if (toggle==2) write(*,*) 'seep backward on cellist'
        endif
        if (toggle==1) then
           pcur=>ntoc(ntodo%idnode)%ptr
        else
           pcur=>last
        endif
        do while (associated(pcur))
!     find the two complemantory nodes in the current cell
            if (verbose) write(*,*) '    working on cell #',pcur%idnode
            call givencomp(amesh,otherninc,pcur,ntodo)
            if (verbose) write(*,*) '    complemantory nodes : ',otherninc
!     loop on the two complemantory nodes
            do i=1,2
               if (verbose) write(*,'(a4,4(i4.4))') 'code',k,ntodo%idnode,pcur%idnode,otherninc(i)
               if (verbose) write(*,*) '       working on complemantory nodes :',otherninc(i)
!     edge propagation
               dedge=dist(ntodo%idnode)+donedge(amesh,ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dist(ntodo%idnode),'+',donedge(amesh,ntodo%idnode,otherninc(i))
               if (verbose) write(*,*) 'dedge : ',dedge
!     face propagation
               if (dist(otherninc(3-i)) /= infinity) then
               if (verbose) write(*,*) 'computing dface because ',dist(otherninc(3-i)),' is not infinity'
                  d12=donedge(amesh,ntodo%idnode,otherninc(3-i))
                  d13=donedge(amesh,otherninc(i),ntodo%idnode)
                  d23=donedge(amesh,otherninc(i),otherninc(3-i))
               if (verbose) write(*,*) 'd12,d13,d23 :'
               if (verbose) write(*,*) d12,d13,d23
                  r1=dist(ntodo%idnode)
                  r2=dist(otherninc(3-i))
               if (verbose) write(*,*) 'r1,r2 : ',r1,r2
                  dface=dcircle(d12,d13,d23,r1,r2)
               if (verbose) write(*,*) 'dface : ',dface
               else
                  dface=infinity
               endif
               dtest=min(dedge,dface)
               if (verbose) write(*,*) 'dtest=min(dedge,dface) : ',dtest
               if (verbose) write(*,*) 'dtest < distance between',k,' and ',otherninc(i)
               if (verbose) write(*,*) dtest,dist(otherninc(i))
               if (dtest < dist(otherninc(i))) then
                   hasbeenupdated=.true.
                   if (verbose) write(*,*) 'better !'
! distance is better : if not in the list, add it
                  dist(otherninc(i))=dtest
                  if (verbose) write(*,*) 'updatelist with :',otherninc(i)
                  call updatelist(otherninc(i),ntodo)
                  if (verbose) write(*,*) 'state of ntodo :'
                  if (verbose) call printlist(ntodo)
               endif
            enddo
! moving on the vicinity list
            if (toggle==1) then
               last=>pcur
               pcur=>pcur%next
            else
               pcur=>pcur%previous
            endif
        enddo
        enddo
        mxswp=max(mxswp,nswp)
! cancelling the cell vicinity list
! removing the first element of the to do list
        if (associated(ntodo%next)) then
           ntodo=>ntodo%next
           deallocate(ntodo%previous)
           nullify(ntodo%previous)
        else
           nullify(ntodo)
        endif
      enddo
! deallocating ntoc
  call deallocntoc(ntoc,amesh%Nnodes)
  write(0,*) 'maximun sweep numner : ',mxswp
end subroutine onevsall2d
!###############################################################################
subroutine printperc(a,b)

  integer :: a,b

!  write(*,*) a,b
  write(*,'(a1,$)') char(8)
  write(*,'(a1,$)') char(8)
  write(*,'(a1,$)') char(8)
!  write(*,*) a,b
  if (a == b) then
     write(*,'(a4)') '100%'
  else
     write(*,'(i2.2,a1,$)') int(100.*(float(a)/float(b))),'%'
  endif
  return
end subroutine printperc
!###############################################################################
subroutine printlist(alist)

  type(listn), pointer :: alist

  type(listn), pointer :: pcur

  pcur=>alist
  do while (associated(pcur))
     write(*,'(i5,$)') pcur%idnode
     pcur=>pcur%next
  enddo
  write(*,*)
end subroutine printlist
!###############################################################################
subroutine updatelist(anode,alist)

  type(listn), pointer :: alist
  integer :: anode

  type(listn), pointer :: pcur,last
  logical :: isinit

! is anode in alist ?
  isinit=.false.
  pcur=>alist
  do while (associated(pcur))
     isinit=(pcur%idnode == anode)
     last=>pcur
     pcur=>pcur%next
     if (isinit) return
  enddo
  allocate(last%next)
  last%next%previous=>last
  nullify(last%next%next)
  last=>last%next
  last%idnode=anode
  return
end subroutine updatelist
!###############################################################################
function isinit(anode,alist)

  type(listn), pointer :: alist
  integer :: anode
  logical :: isinit

  type(listn), pointer :: pcur

  isinit=.false.
  pcur=>alist
  do while (associated(pcur))
     isinit=(pcur%idnode == anode)
     if (isinit) return
     pcur=>pcur%next
  enddo
  return
end function isinit
!###############################################################################
function donedge(amesh,i,j)

  type(mesh) :: amesh
  integer :: i,j
  real :: donedge

  donedge=sqrt((amesh%px(i)-amesh%px(j))**2+&
               (amesh%py(i)-amesh%py(j))**2+&
               (amesh%pz(i)-amesh%pz(j))**2)
  return
end function donedge
!###############################################################################
subroutine updclosedist(amesh,dist,k,i,j)

   type(mesh) :: amesh
   real, dimension(:,:), allocatable :: dist
   integer :: i,j,k

   if (dist(amesh%cell(k,i),amesh%cell(k,j)) == infinity) then
      dist(amesh%cell(k,i),amesh%cell(k,j))=&
           sqrt((amesh%px(amesh%cell(k,i))-amesh%px(amesh%cell(k,j)))**2+&
                (amesh%py(amesh%cell(k,i))-amesh%py(amesh%cell(k,j)))**2+&
                (amesh%pz(amesh%cell(k,i))-amesh%pz(amesh%cell(k,j)))**2)
! reciprocity
      dist(amesh%cell(k,j),amesh%cell(k,i))=dist(amesh%cell(k,i),amesh%cell(k,j))
   endif
end subroutine updclosedist
!###############################################################################
subroutine givencomp(amesh,otherninc,pc,pa)

  type(mesh) :: amesh
  integer, dimension(2) :: otherninc
  type(listn), pointer :: pc,pa

  integer :: i,k

  k=1
  do i=1,3
     if (amesh%cell(pc%idnode,i) /= pa%idnode) then
        otherninc(k)=amesh%cell(pc%idnode,i)
        k=k+1
     endif
  enddo
  return
end subroutine givencomp
!###############################################################################
function dcircle(d12,d13,d23,r1,r2)
! computes the coordinates in 2D plane of three vertices V1,V2,V3 making the
! assumptions that V1 is the origin (0,0), V2 is on the x axis (0,d12). The
! first part gives the coordinates x,y of V3.
! Same set of equations are used to estimate the origin of a point distant by
! r1 and r2 from V1 and V2 respectively (xc,+/-yc). The function returns the
! distance between (xc,-yc) and V0 (x,y)

  real :: dcircle,d12,d13,d23,r1,r2

  real :: x,y,xc,yc,a

  x=(d12**2-d23**2+d13**2)/(2*d12)
  y=sqrt(max(d13**2-x**2,0.))
  if (verbose) write(*,*) 'dcircle - x,y :',x,y
  xc=(d12**2-r2**2+r1**2)/(2*d12)
  yc=sqrt(max(r1**2-xc**2,0.))
  if (verbose) write(*,*) 'dcircle - xc,yc :',xc,yc
  a=abs(xc-x)/(1+(yc/y))
  if (x < xc) then
     a=x+a
  else
     a=x-a
  endif
  if (a < 0. .or. a > d12) then
     dcircle=infinity
     if (verbose) write(*,*) 'dcircle - a : ',a
     if (verbose) write(*,*) 'dcircle - a outside range'
     return
  endif
  dcircle=sqrt((x-xc)**2+(y+yc)**2.)  ! remember ... distance to (xc,-yc)
  return
end function dcircle
!###############################################################################
subroutine startdist(amesh,anode,nodelist,dist)

! Compute the distance at the neightbor nodes list (nodelist) of a given node (anode)
! This routine is needed by onevsall routine after the initial vicinode call.
! The routine allvsall uses the updclosedist routine instead

  type(mesh) :: amesh
  type(listn), pointer :: nodelist
  integer :: anode
  real, dimension(amesh%Nnodes) :: dist

  type(listn), pointer :: pcur

  pcur=>nodelist
  do while (associated(pcur))
     dist(pcur%idnode)=donedge(amesh,anode,pcur%idnode)
     pcur=>pcur%next
  enddo
  return
end subroutine startdist
!###############################################################################
subroutine vicinode(amesh,anode,celllist,nodelist)

! Compute the neightbor nodes list (nodelist) of a given node (anode)
! and its given neighbor cell list (celllist).

  type(mesh) :: amesh
  type(listn), pointer :: celllist,nodelist
  integer :: anode

  type(listn), pointer :: curcell,curnode
  logical :: begin
  integer :: i,j

  if (verbose) write(*,*) 'celllist%idnode',celllist%idnode
  begin=.true.
  curcell=>celllist
  do while (associated(curcell))
     do i=1,3
        if (amesh%cell(curcell%idnode,i).ne.anode) then
           if (begin) then
              allocate(nodelist)
              nullify(nodelist%previous)
              nullify(nodelist%next)
              curnode=>nodelist
              curnode%idnode=amesh%cell(curcell%idnode,i)
              begin=.false.
           else
              if (.not.isinit(amesh%cell(curcell%idnode,i),nodelist)) then
                 allocate(curnode%next)
                 curnode%next%previous=>curnode
                 nullify(curnode%next%next)
                 curnode=>curnode%next
                 curnode%idnode=amesh%cell(curcell%idnode,i)
              endif
           endif
        endif
     enddo
     curcell=>curcell%next
  enddo
end subroutine vicinode
!###############################################################################
subroutine compntoc(amesh,ntoc)
! compute a node to cell array. ntoc is a list array. it is defined by an array
! pointing at the first element of the cell list.
! warning : but fortran does not understand pointer arrays, so I use a container.
! ntoc is not tecnically a pointer array but an array of container containing a
! pointer.
! "last" is a pointer array which, for each node, points at the end of the cell list.

  type(mesh) :: amesh
  type(container), dimension(amesh%Nnodes) :: ntoc

  type(container), dimension(amesh%Nnodes) :: last
  integer :: i,j
  logical :: begin

  if (verbose) write(*,*) 'entering in compntoc'
! initialisation at null
  do i=1,amesh%Nnodes
     nullify(ntoc(i)%ptr)
     nullify(last(i)%ptr)
  enddo
! loop on the cell number
  do i=1,amesh%Ncells
! for each cell, loop on the nodes forming the cell
     do j=1,3
! if the cell list associated to the point amesh%cell(i,j) is not started, do it
! and point "last" to it
        if (.not.associated(ntoc(amesh%cell(i,j))%ptr)) then
           allocate(ntoc(amesh%cell(i,j))%ptr)
           nullify(ntoc(amesh%cell(i,j))%ptr%next)
           nullify(ntoc(amesh%cell(i,j))%ptr%previous)
           last(amesh%cell(i,j))%ptr=>ntoc(amesh%cell(i,j))%ptr
        else
! else add an element to the list after the last element of the list
           allocate(last(amesh%cell(i,j))%ptr%next)
           nullify(last(amesh%cell(i,j))%ptr%next%next)
           last(amesh%cell(i,j))%ptr%next%previous=>last(amesh%cell(i,j))%ptr
           last(amesh%cell(i,j))%ptr=>last(amesh%cell(i,j))%ptr%next
        endif
! add the cell i to the cell list for the node amesh%cell(i,j)
        last(amesh%cell(i,j))%ptr%idnode=i
!       if (verbose) write(*,*) 'last(amesh%cell(i,j))%ptr%idnode',last(amesh%cell(i,j))%ptr%idnode
     enddo
  enddo
  if (verbose) write(*,*) 'check on cell list on node :',amesh%Nnodes/2
  if (verbose) call printlist(ntoc(amesh%Nnodes/2)%ptr)
  if (verbose) write(*,*) 'ntoc(500)%ptr%idnode',ntoc(500)%ptr%idnode
end subroutine compntoc
!###############################################################################
subroutine deallocntoc(ntoc,n)
   integer :: n
   type(container), dimension(n) :: ntoc

   type(listn), pointer :: p1,p2
   integer :: i

   do i=1,n
      p2=>ntoc(i)%ptr
      do while (associated(p2%next))
         p1=>p2
         p2=>p2%next
         deallocate(p1)
      enddo
   enddo
end subroutine deallocntoc
!###############################################################################
subroutine dumpmeshvtk(dev,amesh)
  type(mesh) :: amesh
  integer :: dev

  integer :: i,j

  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a8)') 'distance'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) amesh%Ncells,amesh%Ncells*4
  do i=1,amesh%Ncells
     write(dev,'(a1,3i7)') '3',(amesh%cell(i,j)-1,j=1,3)
  enddo
  call flush(dev)
  return
end subroutine dumpmeshvtk
!###############################################################################
subroutine dumpnodeattributevtk(dev,amesh,field,attname,init)
 type(mesh) :: amesh
  real, dimension(amesh%Nnodes) :: field
  integer :: dev
  character*(*) :: attname
  logical :: init

  integer :: i

  if (init) write(dev,'(a11,i5)') 'POINT_DATA ',amesh%Nnodes
  write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpnodeattributevtk
!###############################################################################
subroutine dumpcellattributevtk(dev,amesh,field,attname,init)
 type(mesh) :: amesh
  real, dimension(amesh%Ncells) :: field
  integer :: dev
  character*(*) :: attname
  logical :: init

  integer :: i

  if (init) write(dev,*)
  if (init) write(dev,'(a10,i5)') 'CELL_DATA ',amesh%Ncells
  write(dev,'(a8,a,a8)') 'SCALARS ',trim(attname),' float 1'
  write(dev,'(a20)') 'LOOKUP_TABLE default'
  do i=1,amesh%Ncells
     write(dev,*) field(i)
  enddo
  return
end subroutine dumpcellattributevtk
!###############################################################################
end module lateration
