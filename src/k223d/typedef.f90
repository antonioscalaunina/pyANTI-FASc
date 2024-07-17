module typedef
use lateration
use utils
implicit none


type borderlist
     integer :: cell
     type(borderlist), pointer :: next
end type borderlist


contains

!==============================================================================
subroutine faultborder(amesh,aborderlist)
! search the cells which form the border of the mesh amesh.
! The cells are returned in a list form aborderlist
type(mesh) :: amesh
type(borderlist), pointer :: aborderlist,curr

integer :: i,j

  nullify(curr)
  do i=1,amesh%Ncells
     do j=1,3
        if (amesh%EToE(i,j) == i) then
           if (.not.associated(curr)) then
              allocate(aborderlist)
              curr=>aborderlist
              nullify(aborderlist%next)
           else
              allocate(curr%next)
              nullify(curr%next%next)
              curr=>curr%next
           endif
           curr%cell=i
           exit
        endif
     enddo
  enddo
  return
end subroutine faultborder
!==============================================================================
subroutine faultsurface(amesh,aborderlist,surfacelist,z)
! define a subset of the cell border list (aborderlist) of the
! mesh amesh which are close to the surface (depth condition lesser
! than a given threshold z). The subset is returned in a list
! format (surfacelist)
type(mesh) :: amesh
type(borderlist), pointer :: aborderlist,bcurr,scurr,surfacelist
real :: z
integer :: dev,cnt,i,j

   nullify(scurr)
   bcurr=>aborderlist
   cnt=0
   do while (associated(bcurr))
      if (minval(abs(amesh%pz(amesh%cell(bcurr%cell,:)))) < z) then
         if (.not.associated(scurr)) then
             allocate(surfacelist)
             scurr=>surfacelist
             nullify(surfacelist%next)
          else
             allocate(scurr%next)
             nullify(scurr%next%next)
             scurr=>scurr%next
          endif
          scurr%cell=bcurr%cell
          cnt=cnt+1
      endif
      bcurr=>bcurr%next
   enddo
!
  dev=23
  open(dev,file='surface.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) cnt,cnt*4
  scurr=>surfacelist
  do while (associated(scurr))
     write(dev,'(a1,3i5)') '3',(amesh%cell(scurr%cell,j)-1,j=1,3)
     scurr=>scurr%next
  enddo
 close(dev)
!
   return
end subroutine faultsurface
!==============================================================================
subroutine quakesurface(amesh,surfacelist)
! correct the array of nodes on the border of the quake removing the nodes
! which belong to the mesh surface list.
type(mesh) :: amesh
type(borderlist), pointer :: surfacelist,scurr

integer :: k,i,j
integer, dimension(:), allocatable :: dummy

integer :: dev
logical :: found
   
!
   k=1
   do while (k <= amesh%quakeborder_nodesno)
      scurr=>surfacelist
      found=.false.
      do while (associated(scurr))
         do i=1,3
            if (amesh%cell(scurr%cell,i) == amesh%quakeborder_nodes(k)) then
               do j=k,amesh%quakeborder_nodesno-1
                  amesh%quakeborder_nodes(j)=amesh%quakeborder_nodes(j+1)
               enddo
               amesh%quakeborder_nodesno=amesh%quakeborder_nodesno-1
               found=.true.
               exit
            endif
         enddo
         if (found) exit
         scurr=>scurr%next
      enddo
      if (.not.found) k=k+1
   enddo
   allocate(dummy(amesh%quakeborder_nodesno))
   dummy=amesh%quakeborder_nodes(1:amesh%quakeborder_nodesno)
   deallocate(amesh%quakeborder_nodes)
   allocate (amesh%quakeborder_nodes(amesh%quakeborder_nodesno))
   amesh%quakeborder_nodes=dummy
   deallocate(dummy)

  dev=23
  open(dev,file='quakeborderlineunsurf.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'VERTICES '
  write(dev,*) amesh%quakeborder_nodesno,2*amesh%quakeborder_nodesno
  do i=1,amesh%quakeborder_nodesno
     write(dev,'(a1,3i5)') '1',amesh%quakeborder_nodes(i)-1
  enddo
  close(dev)
   return
end subroutine  quakesurface
!==============================================================================
subroutine readfile(amesh)
type(mesh) :: amesh
integer :: i, j, numb_area, numb_slip 
character(len=200) :: fnamef, file_matrix_string
character :: dummy_char
! read in mesh from file
call read_mesh(amesh)
!call read_mesh_vtk(amesh)

if (amesh%nuc_mode==0) then
 call read_nuc(amesh)
 call read_prop(amesh)
endif

! area is in S.I.
call calc_area(amesh)



if (amesh%distance==0) then
 ! calculate distances
 call allvsall2d(amesh)
else
 open(1111,file='matrix_string.txt',form='formatted',action='read')
 read(1111,"(A)") file_matrix_string
 close(1111) 
 fnamef = trim(file_matrix_string) ! trim(amesh%geo_zone) // '_matrix_distance.dat'
 !write(0,*) fnamef
 open(11,file=fnamef,form='unformatted',access='direct',recl=amesh%Nnodes*4)
 allocate(amesh%dist(amesh%Nnodes,amesh%Nnodes))
 do i=1,amesh%Nnodes
   read(11,rec=i) (amesh%dist(i,j), j=1,amesh%Nnodes) 
 enddo
 close(11)
endif

open(11,file='index_file.dat',form='formatted',action='read')
read(11,*) numb_area, numb_slip
amesh%numb_scen=numb_area*numb_slip
write(0,*) "Number of scenarios is", amesh%numb_scen
write(*,*) "Number of scenarios is", amesh%numb_scen
allocate(amesh%index_string(amesh%numb_scen),amesh%numb_gauss(amesh%numb_scen))
do i=1,amesh%numb_scen
   read(11,"(A,A,I1)") amesh%index_string(i),dummy_char, amesh%numb_gauss(i)
   !print*, amesh%index_string(i), amesh%numb_gauss(i)
enddo
close(11)



! build connectivity between elements
call tiConnect2d(amesh)


end subroutine readfile
!==============================================================================
subroutine read_nuc(amesh)
type(mesh) :: amesh
integer :: stat,i,j
real :: lat, lon, x, y
integer :: iway,UTM_PROJECTION_ZONE
logical :: SUPPRESS_UTM_PROJECTION
integer :: NoNucLocs,INOUT
real,allocatable,dimension(:) :: nuc_x,nuc_y
!open(13,file='Nucleation_CyA.txt',form = 'formatted', action = 'read')
!open(13,file='Nucleation_HeAe.txt',form = 'formatted', action = 'read')
!open(13,file='Nucleation_HeAw.txt',form = 'formatted', action = 'read')
open(13,file='SoA_propagation.txt',form = 'formatted', action = 'read')
iway = 0 ! convert from Lat/Long to UTM

SUPPRESS_UTM_PROJECTION = .false.
UTM_PROJECTION_ZONE = amesh%merc_zone !27 !34 36 33 ! calabrian

! work out how long file is
i = 0
stat = 0
do while (stat == 0)
      read(13, fmt=*, iostat=stat)
      i = i+1
enddo
print*, i
NoNucLocs = i-2
rewind(13)
read(13,*) ! header
allocate(nuc_x(NoNucLocs),nuc_y(NoNucLocs))

do i = 1,NoNucLocs
    read(13,*) lon, lat
! convert coordinates to utm
    call utm_geo(lon,lat,x,y,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)
    nuc_x(i) = x
    nuc_y(i) = y
enddo
close(13)


! write out utm vtk version of nucleation coordinates
open(10,file="utm_nuc.vtk")
write(10,'(a26)') '# vtk DataFile Version 2.0'
write(10,'(a8)') 'distance'
write(10,'(a5)') 'ASCII'
write(10,'(a16)') 'DATASET POLYDATA'
write(10,'(a7,$)') 'POINTS '
write(10,*) NoNucLocs,' float'
do i=1,NoNucLocs
     write(10,*) nuc_x(i),nuc_y(i),0.
enddo
write(10,'(a9,$)') 'POLYGONS '
write(10,*) 1, NoNucLocs+1
write(10,'(i4,$)') NoNucLocs
write(10,*) (j,j=0,NoNucLocs-1)
close(10)


print*,NoNucLocs
! calculate how many nucleation points there are
i = 0
do j = 1,amesh%Nnodes
          call pnpoly(amesh%px(j),amesh%py(j),nuc_x,nuc_y,NoNucLocs,INOUT)
          if (INOUT >= 0) then
              i = i+1
          endif
enddo
print*,'No of Nucleation Locations...',i
amesh%NoNucLocs = i
allocate(amesh%NucLocs(amesh%NoNucLocs))
if (amesh%NoNucLocs==0) then
    write(*,*)'ERROR there are no nucleation points!!!!!'
endif

i = 0
do j = 1,amesh%Nnodes
          call pnpoly(amesh%px(j),amesh%py(j),nuc_x,nuc_y,NoNucLocs,INOUT)
          if (INOUT>=0) then
               i = i+1
               amesh%NucLocs(i) = j
          endif
enddo

return
end subroutine
!==============================================================================
subroutine read_prop(amesh)
type(mesh) :: amesh
integer :: stat,i,j
real :: lat, lon, x, y
integer :: iway,UTM_PROJECTION_ZONE
logical :: SUPPRESS_UTM_PROJECTION
integer :: NoPropLocs,INOUT
real,allocatable,dimension(:) :: prop_x,prop_y
!open(13,file='Nucleation_CyA.txt',form = 'formatted', action = 'read')
!open(13,file='Nucleation_HeAe.txt',form = 'formatted', action = 'read')
!open(13,file='Nucleation_HeAw.txt',form = 'formatted', action = 'read')
open(13,file='CaA_propagation.txt',form = 'formatted', action = 'read')
iway = 0 ! convert from Lat/Long to UTM

SUPPRESS_UTM_PROJECTION = .false.
UTM_PROJECTION_ZONE =  amesh%merc_zone !34 36 33 ! calabrian

! work out how long file is
i = 0
stat = 0
do while (stat == 0)
      read(13, fmt=*, iostat=stat)
      i = i+1
enddo
NoPropLocs = i-2
rewind(13)
read(13,*) ! header
allocate(prop_x(NoPropLocs),prop_y(NoPropLocs))

do i = 1,NoPropLocs
    read(13,*) lon, lat
! convert coordinates to utm
    call utm_geo(lon,lat,x,y,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)
    prop_x(i) = x
    prop_y(i) = y
enddo
close(13)


! write out utm vtk version of propagation coordinates
open(10,file="utm_prop.vtk")
write(10,'(a26)') '# vtk DataFile Version 2.0'
write(10,'(a8)') 'distance'
write(10,'(a5)') 'ASCII'
write(10,'(a16)') 'DATASET POLYDATA'
write(10,'(a7,$)') 'POINTS '
write(10,*) NoPropLocs,' float'
do i=1,NoPropLocs
     write(10,*) prop_x(i),prop_y(i),0.
enddo
write(10,'(a9,$)') 'POLYGONS '
write(10,*) 1, NoPropLocs+1
write(10,'(i4,$)') NoPropLocs
write(10,*) (j,j=0,NoPropLocs-1)
close(10)


print*,NoPropLocs
! calculate how many nucleation points there are
i = 0
do j = 1,amesh%Nnodes
          call pnpoly(amesh%px(j),amesh%py(j),prop_x,prop_y,NoPropLocs,INOUT)
          if (INOUT >= 0) then
              i = i+1
          endif
enddo
print*,'No of Propagation Locations...',i
amesh%NoPropLocs = i
allocate(amesh%propLocs(amesh%NoPropLocs))
if (amesh%NoPropLocs==0) then
    write(*,*)'ERROR there are no nucleation points!!!!!'
endif

i = 0
do j = 1,amesh%Nnodes
          call pnpoly(amesh%px(j),amesh%py(j),prop_x,prop_y,NoPropLocs,INOUT)
          if (INOUT>=0) then
               i = i+1
               amesh%PropLocs(i) = j
          endif
enddo

return
end subroutine
!==============================================================================
subroutine read_mesh_vtk(amesh)
type(mesh) :: amesh
real :: x(3),y(3),z(3)
integer :: i,j,id,inx
character(len=60) :: str
integer :: a,b,c
open(13,file='SmoothEdgeMeshV3.vtk',form = 'formatted', action = 'read')
! skip header
do i = 1,4
       read(13,*)
enddo
! find end of nodal section & element section
read(13,'(a6,i5,4a)') str,amesh%Nnodes,str
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))

inx = 1
! read in vertices
do i = 1,amesh%Nnodes/3
      read(13,*) x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3)

      do j = 1,3
       amesh%px(inx) = x(j)
       amesh%py(inx) = y(j)
       amesh%pz(inx) = z(j)
       inx = inx+1
     enddo
enddo
! skip blank
read(13,*)
! read in number of elements
read(13,'(a8,i5,i5)') str,i,j
amesh%Ncells = i
allocate(amesh%cell(amesh%Ncells,3))
! read element to vertices
do i = 1,amesh%Ncells
       read(13,*) id,a,b,c
       amesh%cell(i,1) = a+1
       amesh%cell(i,2) = b+1
       amesh%cell(i,3) = c+1
!       print*,amesh%cell(i,1),amesh%cell(i,2),amesh%cell(i,3)
!       pause
enddo

close(13)
end subroutine read_mesh_vtk
!==============================================================================
subroutine read_mesh(amesh)
type(mesh) :: amesh
integer :: i,j,id
character(len=60) :: str, fnamef
integer :: iway,UTM_PROJECTION_ZONE
logical :: SUPPRESS_UTM_PROJECTION
real :: lon,lat,rake


!open(13,file='CA_mesh_15km.inp',form = 'formatted', action = 'read')
!open(13,file='CY_mesh_15km.inp',form = 'formatted', action = 'read')
fnamef = trim(amesh%geo_zone) // '_mesh_15km.inp'
open(13,file=fnamef,form = 'formatted', action = 'read')

iway = 0 ! convert from Lat/Long to UTM
SUPPRESS_UTM_PROJECTION = .false.
!if (amesh%geo_zone=='CaA') then
!   UTM_PROJECTION_ZONE = 33 !34 36 33 ! calabrian
!elseif (amesh%geo_zone=='HeA') then
!   UTM_PROJECTION_ZONE = 34
!elseif (amesh%geo_zone=='CyA') then
!   UTM_PROJECTION_ZONE = 36
!elseif (amesh%geo_zone=='CaR') then
!   UTM_PROJECTION_ZONE = 20
!elseif (amesh%geo_zone=='AlR'.or.amesh%geo_zone=='AlY') then
!   UTM_PROJECTION_ZONE = 30
!elseif (amesh%geo_zone=='HoF') then
!   UTM_PROJECTION_ZONE = 27
!else
!   stop 'Geo zone unknown/syntax error'
!endif
UTM_PROJECTION_ZONE = amesh%merc_zone
! skip header
do i = 1,9
       read(13,*)
enddo
! find end of nodal section & element section
i = 0
id = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Nnodes = i-1
print*,'No. of Nodes ....',amesh%Nnodes
allocate(amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes))
allocate(amesh%lon(amesh%Nnodes),amesh%lat(amesh%Nnodes))

do i = 1,2
       read(13,*) str
enddo
id = 0
i = 0
do while ( id == 0 )
    read(13,*) str
    if (trim(str) == '**') then
        id = 1
    endif
    i = i+1
enddo
amesh%Ncells = i-1
print*,'No of cells......',amesh%Ncells
allocate(amesh%cell(amesh%Ncells,3))
allocate(amesh%HySea(amesh%Ncells,10))

!====== now read data in
rewind(13)
! skip header
do i = 1,9
       read(13,*)
enddo


do i = 1,amesh%Nnodes
       read(13,*) id,lon,lat,amesh%pz(i)
       if (amesh%coord_type==0) then
       ! convert coordinates to utm
         call utm_geo(lon,lat,amesh%px(i),amesh%py(i),UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)
              amesh%lon(i)=lon
              amesh%lat(i)=lat
       else
         amesh%px(i)=lon; amesh%py(i)=lat
       endif
enddo
print*,amesh%px(amesh%Nnodes),amesh%py(amesh%Nnodes),amesh%pz(amesh%Nnodes)

! skip header
do i = 1,3
       read(13,*) str
enddo
! read element to vertices and ! build a matrix for the HySea inputs
rake=90.
do i = 1,amesh%Ncells
       read(13,*) id,amesh%cell(i,1),amesh%cell(i,2),amesh%cell(i,3)
       amesh%HySea(i,1)=amesh%lon(amesh%cell(i,1))
       amesh%HySea(i,2)=amesh%lat(amesh%cell(i,1))
       amesh%HySea(i,3)=-amesh%pz(amesh%cell(i,1))/1000
       amesh%HySea(i,4)=amesh%lon(amesh%cell(i,2))
       amesh%HySea(i,5)=amesh%lat(amesh%cell(i,2))
       amesh%HySea(i,6)=-amesh%pz(amesh%cell(i,2))/1000
       amesh%HySea(i,7)=amesh%lon(amesh%cell(i,3))
       amesh%HySea(i,8)=amesh%lat(amesh%cell(i,3))
       amesh%HySea(i,9)=-amesh%pz(amesh%cell(i,3))/1000
       amesh%HySea(i,10)=rake
enddo
print*,id,amesh%cell(amesh%Ncells,1),amesh%cell(amesh%Ncells,2),amesh%cell(amesh%Ncells,3)
close(13)

end subroutine read_mesh
!==============================================================================
subroutine calc_area(amesh)
implicit none
type(mesh) :: amesh
integer :: i
real,dimension(3) :: x,y,z,vec_ab,vec_ac
real :: ab,ac,area,res,theta

allocate(amesh%area(amesh%Ncells))

do i = 1,amesh%Ncells
	x = amesh%px(amesh%cell(i,:))
	y = amesh%py(amesh%cell(i,:))
	z = amesh%pz(amesh%cell(i,:))
       vec_ab = (/ x(2)-x(1), y(2)-y(1), z(2)-z(1) /)
       vec_ac = (/ x(3)-x(1), y(3)-y(1), z(3)-z(1) /)
	ac = norm2(vec_ac)
	ab = norm2(vec_ab)

	res = dot_product(vec_ab, vec_ac)/ab/ac
	theta = acos(res)
	amesh%area(i) = ab*ac*sin(theta)/float(2)
enddo

end subroutine calc_area
!==============================================================================
!==============================================================================
subroutine tiConnect2d(amesh)
type(mesh) :: amesh
integer,allocatable,dimension(:,:) :: fnodes,spNodeToNode,sorted,matchL,matchR
integer,allocatable,dimension(:) :: id,buf_EToE,buf_EToF,indices_d,indices
!integer,allocatable,dimension(:) :: id,buf_EToE,buf_EToF,imult,indices_d,indices

integer :: Nfaces
integer :: irow,krow,nsize,inx,i,j

Nfaces = 3
allocate (fnodes(Nfaces*amesh%Ncells,2))
! create list of all faces 1, then 2, & 3
fnodes(1:amesh%Ncells,:) = amesh%cell(:,1:2)
fnodes(amesh%Ncells+1:2*amesh%Ncells,:) = amesh%cell(:,2:3)
fnodes(2*amesh%Ncells+1:3*amesh%Ncells,:) = amesh%cell(:,(/ 3, 1 /))


!fnodes = sort(fnodes,2)-1;
nsize = Nfaces*amesh%Ncells
do irow = 1, nsize
         krow = minloc( fnodes( irow, : ),dim=1)
         if (krow == 2) then
              fnodes( irow, : ) = fnodes( irow, (/ 2, 1 /) )
         endif
enddo

fnodes = fnodes-1

! set up default element to element and element to faces connectivity
allocate(amesh%EToE(amesh%Ncells,3),amesh%EToF(amesh%Ncells,3))
do i = 1,amesh%Ncells
       amesh%EToE(i,:) = i
       amesh%EToF(i,:) = (/ 1 , 2, 3 /)
enddo

! uniquely number each set of three faces by their node numbers
allocate(id(Nfaces*amesh%Ncells),spNodeToNode(Nfaces*amesh%Ncells,Nfaces+1))
id = fnodes(:,1)*amesh%Nnodes + fnodes(:,2)+1
deallocate(fnodes)
allocate(buf_EToE(Nfaces*amesh%Ncells),buf_EToF(Nfaces*amesh%Ncells))
buf_EToE = reshape(amesh%EToE,(/ Nfaces*amesh%Ncells/))
buf_EToF = reshape(amesh%EToF,(/ Nfaces*amesh%Ncells/))

do i = 1,Nfaces*amesh%Ncells
       spNodeToNode(i,:) = (/ id(i), i, buf_EToE(i), buf_EToF(i) /)
enddo
deallocate(buf_EToE,buf_EToF,id)

! Now we sort by global face number.
!allocate(imult(Nfaces*K),sorted(Nfaces*K,Nfaces+1))
allocate(sorted(Nfaces*amesh%Ncells,Nfaces+1))

!call sort_rows(spNodeToNode,Nfaces*K,Nfaces+1,imult)
call sort_rows(spNodeToNode,Nfaces*amesh%Ncells,Nfaces+1)
sorted = spNodeToNode;!(imult,:)

! find matches in the sorted face list
allocate(indices_d(Nfaces*amesh%Ncells-1))

inx = 0;
do i = 1,Nfaces*amesh%Ncells-1
       if (sorted(i,1) == sorted(i+1,1)) then
           inx = inx+1
           indices_d(inx) = i
       endif
enddo


allocate(indices(inx))
indices = indices_d(1:inx)
deallocate(indices_d)

! make links reflexive

allocate(matchL(inx*2,4),matchR(inx*2,4))
matchL(1:inx,:) = sorted(indices,:)
matchL(inx+1:2*inx,:) = sorted(indices+1,:)

matchR(1:inx,:) = sorted(indices+1,:)
matchR(inx+1:2*inx,:) = sorted(indices,:)

! insert matches
do i = 1,inx*2

		if (matchL(i,2) <=  amesh%Ncells) then
	       amesh%EToE(matchL(i,2),1) = matchR(i,3);
    	   amesh%EToF(matchL(i,2),1) = matchR(i,4);
		elseif (amesh%Ncells <  matchL(i,2).and. matchL(i,2) <=  2*amesh%Ncells) then
	       amesh%EToE(matchL(i,2)-amesh%Ncells,2) = matchR(i,3);
    	   amesh%EToF(matchL(i,2)-amesh%Ncells,2) = matchR(i,4);
		elseif (amesh%Ncells*2 <  matchL(i,2).and. matchL(i,2) <=  3*amesh%Ncells) then
    	   amesh%EToE(matchL(i,2)-amesh%Ncells*2,3) = matchR(i,3);
    	   amesh%EToF(matchL(i,2)-amesh%Ncells*2,3) = matchR(i,4);
    	elseif (amesh%Ncells*2 <  matchL(i,2).and. matchL(i,2) <=  3*amesh%Ncells) then
    	   amesh%EToE(matchL(i,2)-amesh%Ncells*3,4) = matchR(i,3);
    	   amesh%EToF(matchL(i,2)-amesh%Ncells*3,4) = matchR(i,4);
		endif
enddo

!open(1000,file='EToE.dat',status='unknown')
!do i=1,amesh%Ncells
!   write(1000,"(3I7)") amesh%EToE(i,:)
!end do
!close(1000)



deallocate(matchR,matchL)
return
end subroutine tiConnect2d
!==============================================================================
!==============================================================================
subroutine  select_fault_zone(amesh,target_area)
type(mesh),intent (inout) :: amesh
real,parameter :: kilo = 1000.
integer,allocatable,dimension(:) :: quake_nodes,quake_ameshents,irow,icol,haveit
real, intent(in) :: target_area
integer :: choice_amesh(3),ameshs(3)
!integer,allocatable,dimension(:) :: ameshs
real :: nuc_x,nuc_y,nuc_z
real :: quake_area
real :: random

integer :: bc_int,int
integer :: inx,get_out,nuc_id,node_idx,nxt_pt,ele_inx
integer :: kk,i,k,j,ii,jj,nc,before,pt,np,idxx,dev,nline,idxxx,nv,pf
logical :: found
integer, dimension(2) :: cj
integer, dimension(3) :: p

integer, allocatable,dimension(:) :: id,idx
integer :: i_nodes(3)!,j_nodes(3)
integer :: n,m,nlength
integer,allocatable,dimension(:,:) :: j_nodes,dummy2
integer,allocatable,dimension(:) :: int_amesh,bc_amesh,bc_nodes,common_node,id_keep,dummy

allocate(quake_nodes(amesh%Nnodes),quake_ameshents(amesh%Ncells))


quake_nodes = 0         ! NOTE: ameshent or node can have an index of 0 !IMP
quake_ameshents = 0      ! NOTE: ameshent or node can have an index of 0 !IMP

quake_area = 0.   ! initial area for earthquake zone on fault


inx = 1
node_idx = 1
ele_inx = 1
get_out = 0
do while (get_out < 1)
       if (inx == 1) then
! randomly pick starting cella
              call random_seed
              call random_number(random)
              write(*,*) "random=", random
            !  nuc_id = ceiling(random*amesh%Nnodes)   !   ! NOTE: setting fault location anywhere on fault!!!!
!               nuc_x = amesh%px(nuc_id)
!               nuc_y = amesh%py(nuc_id);
!               nuc_z = amesh%pz(nuc_id);
             nuc_id = ceiling(random*amesh%NoNucLocs)   !   ! NOTE: setting fault location based on nucleation zone!!!!
             nuc_id = amesh%NucLocs(nuc_id)
             nuc_x = amesh%px(nuc_id)
             nuc_y = amesh%py(nuc_id);
             nuc_z = amesh%pz(nuc_id);
              write(*,*) 'Nucleation id....',nuc_id
              write(*,*) nuc_x,nuc_y,nuc_z
              nxt_pt = nuc_id
       else
! take next nodal point
              node_idx = node_idx+1
              if (node_idx > size(quake_nodes)) then
                  write(*,*) '***WARNING: requested Mag. poss. too big ***'
                  get_out = 1
                  cycle
              endif
              nxt_pt = quake_nodes(node_idx)
       endif
! find ameshents connected to node
       call find(amesh%cell,'==',nxt_pt,irow,icol)       ! find the index at where amesh%dist = 0
       do kk = 1,size(irow)
              choice_amesh(:) = amesh%EToE(irow(kk),:)   !take ameshent with a joining face
              do i = 1,3
                     if (inx == 1) then
                            quake_nodes(inx) = amesh%cell(choice_amesh(i),1);
                            quake_nodes(inx+1) = amesh%cell(choice_amesh(i),2);
                            quake_nodes(inx+2) = amesh%cell(choice_amesh(i),3);
                            inx = inx+3;
                     else               ! check to see if node is already in list
                         do k = 1,3
                            call find(quake_nodes,'==',amesh%cell(choice_amesh(i),k),haveit)
                            if (haveit(1) == -1) then      ! don't have node
                                   quake_nodes(inx) =  amesh%cell(choice_amesh(i),k);
                                   inx = inx+1;
                            endif
                         enddo
                     endif
! check to see if element is already in quake_ameshents
                    call find(quake_ameshents,'==',choice_amesh(i),haveit)
                    if (haveit(1) == -1) then      ! don't have ameshent
                           quake_ameshents(ele_inx) = choice_amesh(i);
                           quake_area = quake_area + amesh%area(quake_ameshents(ele_inx)) ! assume area is in m^2
                           if (quake_area >= target_area) then
                                    get_out = 1;
                           endif
                           ele_inx = ele_inx+1;
                     endif

              enddo

       enddo

enddo


! Save section of fault used for earthquake
allocate(amesh%QuakeElem(ele_inx-1),amesh%QuakeNodes(inx-1))
amesh%QuakeElem = quake_ameshents(1:ele_inx-1)
amesh%QuakeNodes = quake_nodes(1:inx-1)
amesh%QuakeElemNo = ele_inx-1
amesh%QuakeNodesNo = inx-1
deallocate(quake_nodes,quake_ameshents)

  dev=23
  open(dev,file='Quake_Area.dat')
  do i=1,amesh%QuakeElemNo
     write(dev,*) amesh%QuakeElem(i)
  enddo

!
  dev=23
  open(dev,file='quake2.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) amesh%QuakeElemNo,amesh%QuakeElemNo*4
  do i=1,amesh%QuakeElemNo
     write(dev,'(a1,3i5)') '3',(amesh%cell(amesh%QuakeElem(i),j)-1,j=1,3)
  enddo
  close(dev)
!
  
! border elements
  amesh%QuakeBorder_elemNo=0
  do i=1,amesh%QuakeElemNo
     nc=0
     do j=1,3
        do k=1,amesh%QuakeElemNo
           if (amesh%EToE(amesh%QuakeElem(i),j) == amesh%QuakeElem(k).and.k /= i) nc=nc+1
        enddo
     enddo
     if (nc < 3) amesh%QuakeBorder_elemNo=amesh%QuakeBorder_elemNo+1
  enddo
  allocate(amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo))
  allocate(dummy2(amesh%QuakeelemNo,7))
  dummy2=0
  amesh%QuakeBorder_elemNo=0
  do i=1,amesh%QuakeElemNo
     nc=0
     do j=1,3
        before=nc
        do k=1,amesh%QuakeElemNo
           if (amesh%EToE(amesh%QuakeElem(i),j) == amesh%QuakeElem(k).and.k /= i) nc=nc+1
        enddo
        if (before == nc) then
           dummy2(i,1)=dummy2(i,1)+1
           idxx=dummy2(i,1)*2
           if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),j)) then
           do ii=1,3
              do jj=1,3
                 if (amesh%cell(amesh%QuakeElem(i),ii) == amesh%cell(amesh%EToE(amesh%QuakeElem(i),j),jj)) then
                    dummy2(i,idxx)=amesh%cell(amesh%QuakeElem(i),ii)
                    idxx=idxx+1
                 endif
              enddo
           enddo
           else
! case where the elements is on the mesh border
             nv=0
             do ii=1,3
                if (amesh%EToE(amesh%QuakeElem(i),ii) == amesh%QuakeElem(i)) nv=nv+1
             enddo
!      one element mesh (should not happened) : three edges facing the void
             if (nv==3) then
                dummy2(i,1)=3
                dummy2(i,2)=amesh%cell(amesh%QuakeElem(i),1)
                dummy2(i,3)=amesh%cell(amesh%QuakeElem(i),2)
                dummy2(i,4)=amesh%cell(amesh%QuakeElem(i),2)
                dummy2(i,5)=amesh%cell(amesh%QuakeElem(i),3)
                dummy2(i,6)=amesh%cell(amesh%QuakeElem(i),1)
                dummy2(i,7)=amesh%cell(amesh%QuakeElem(i),3)
                exit
             endif
!      two edges facing the void
             if (nv==2) then
                idxxx=1
!           finding the unique neighbour cell
                do ii=1,3
                   if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),ii)) then
                      cj(1)=amesh%EToE(amesh%QuakeElem(i),ii)
                      exit
                   endif
                enddo
!           finding the two common vertices between the reference cell and its unique neighbour
                do ii=1,3
                   do jj=1,3
                      if (amesh%cell(amesh%QuakeElem(i),ii) == amesh%cell(cj(1),jj)) then
                         p(idxxx)=amesh%cell(amesh%QuakeElem(i),ii)
                         idxxx=idxxx+1
                      endif
                   enddo
                enddo
!           finding the complemantary vertice of the reference cell
                do ii=1,3
                   pf=0
                   do jj=1,2
                      if (amesh%cell(amesh%QuakeElem(i),ii) /= p(jj)) pf=pf+1
                   enddo
                   if (pf==0) p(3)=amesh%cell(amesh%QuakeElem(i),ii)
                enddo
!           saving the two edges facing the void
                dummy2(i,1)=2
                dummy2(i,2)=p(1)
                dummy2(i,3)=p(3)
                dummy2(i,4)=p(2)
                dummy2(i,5)=p(3)
                exit
           endif
!      one edges facing the void
           if (nv==1) then
!           finding the two neighbor cells
              idxxx=1
              do ii=1,3
                 if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),ii)) then
                    cj(idxxx)=amesh%EToE(amesh%QuakeElem(i),ii)
                    idxxx=idxxx+1
                 endif
              enddo
!          finding the common vertice between the two neighbor cells
              do ii=1,3
                 do jj=1,3
                    if (amesh%cell(cj(1),ii) == amesh%cell(cj(2),jj)) then
                       p(1)=amesh%cell(cj(1),ii)
                    endif
                 enddo
              enddo
!          finding the two complementary vertices which by contruction formed the edge facing the void
              idxxx=2
              do ii=1,3
                 if (amesh%cell(amesh%QuakeElem(i),ii) /= p(1)) then
                    p(idxxx)=amesh%cell(amesh%QuakeElem(i),ii)
                    idxxx=idxxx+1
                 endif
              enddo
              dummy2(i,1)=1
              dummy2(i,2)=p(2)
              dummy2(i,3)=p(3)
           endif
          endif
        endif
     enddo
     if (nc < 3) then
        amesh%QuakeBorder_elemNo=amesh%QuakeBorder_elemNo+1
        amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo)=amesh%QuakeElem(i)
     endif
  enddo

  dev=23
  open(dev,file='quakeborder.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) amesh%QuakeBorder_elemNo,amesh%QuakeBorder_elemNo*4
  do i=1,amesh%QuakeBorder_elemNo
     write(dev,'(a1,3i5)') '3',(amesh%cell(amesh%QuakeBorder_elem(i),j)-1,j=1,3)
  enddo
  close(dev)
!

! border nodes

  allocate(dummy(amesh%QuakeBorder_elemNo*3))
  np=0
  do i=1,amesh%QuakeelemNo
     if (dummy2(i,1) > 0) then
     do j=1,dummy2(i,1)
        if (np==0) then
           np=2
           dummy(1)=dummy2(i,2)
           dummy(2)=dummy2(i,3)
           cycle
        endif
        do k=1,2
           pt=dummy2(i,j*2+k-1)
           found=.false.
           do ii=1,np
              if (pt==dummy(ii)) then
                 found=.true.
                 exit
              endif
           enddo
           if (.not.found) then
              np=np+1
              dummy(np)=pt
           endif
        enddo
     enddo
     endif
  enddo
  amesh%QuakeBorder_NodesNo=np
  allocate(amesh%QuakeBorder_Nodes(np))
  amesh%QuakeBorder_Nodes=dummy(1:np)


  nline=0
  do i=1,amesh%QuakeelemNo
     nline=nline+dummy2(i,1)
  enddo
  dev=23
  open(dev,file='quakeborderline.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a6,$)') 'LINES '
  write(dev,*) nline,3*nline
  do i=1,amesh%QuakeelemNo
     if (dummy2(i,1) /= 0) then
         do j=1,dummy2(i,1)
            write(dev,'(a1,3i5)') '2',(dummy2(i,j*2+k-1)-1,k=1,2)
         enddo
     endif
  enddo
  close(dev)
!
  deallocate(dummy,dummy2)

end subroutine select_fault_zone
!==============================================================================
subroutine  select_fault_zone_fromfile(amesh,target_area,sc_index)
        
type(mesh),intent (inout) :: amesh
integer, intent(in) :: sc_index
real,parameter :: kilo = 1000.
integer,allocatable,dimension(:) :: quake_nodes,quake_ameshents,irow,icol,haveit
real, intent(in) :: target_area
integer :: choice_amesh(3),ameshs(3)
!integer,allocatable,dimension(:) :: ameshs
real :: nuc_x,nuc_y,nuc_z
real :: quake_area
real :: random
integer :: bc_int,int
integer :: inx,get_out,nuc_id,node_idx,nxt_pt,ele_inx,INOUT
integer :: kk,i,k,j,ii,jj,nc,before,pt,np,idxx,dev,nline,idxxx,nv,pf
logical :: found
integer, dimension(2) :: cj
integer, dimension(3) :: p

integer, allocatable,dimension(:) ::idx
integer :: i_nodes(3)!,j_nodes(3)
integer :: n,m,nlength
integer,allocatable,dimension(:,:) :: j_nodes,dummy2
integer,allocatable,dimension(:) :: int_amesh,bc_amesh,bc_nodes,common_node,id_keep,dummy
character*200 :: file_out, index_local

allocate(quake_nodes(amesh%Nnodes))


quake_nodes = 0         ! NOTE: ameshent or node can have an index of 0 !IMP
!quake_ameshents = 0      ! NOTE: ameshent or node can have an index of 0 !IMP

quake_area = 0.   ! initial area for earthquake zone on fault

index_local=amesh%index_string(sc_index)
index_local=index_local(1:5)
file_out='QuakeArea_' // trim(index_local) // '.dat'
   file_out=trim(file_out)
open(11,file=file_out)

allocate(dummy(amesh%Ncells))
dummy=0
do i=1,amesh%Ncells
   read(11,*,end=80) dummy(i)
enddo

80 continue
amesh%QuakeElemNo = i-1
allocate (amesh%QuakeElem(i-1))
amesh%QuakeElem = dummy(1:i-1)
close(11)
deallocate(dummy)

inx = 1
quake_nodes(inx)=amesh%cell(amesh%QuakeElem(1),1)
inx=inx+1
quake_nodes(inx)=amesh%cell(amesh%QuakeElem(1),2)
inx=inx+1
quake_nodes(inx)=amesh%cell(amesh%QuakeElem(1),3)
inx=inx+1

do i=2,amesh%QuakeElemNo
   do j=1,3
      call find(quake_nodes,'==',amesh%cell(amesh%QuakeElem(i),j),haveit)
      if (haveit(1) == -1) then
      quake_nodes(inx) = amesh%cell(amesh%QuakeElem(i),j)
      inx = inx + 1
      endif
   enddo
enddo
   


! Save section of fault used for earthquake
allocate(amesh%QuakeNodes(inx-1))
amesh%QuakeNodes = quake_nodes(1:inx-1)
amesh%QuakeNodesNo = inx-1
deallocate(quake_nodes) !,quake_ameshents)



dev=23
  open(dev,file='quake2.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) amesh%QuakeElemNo,amesh%QuakeElemNo*4
  do i=1,amesh%QuakeElemNo
     write(dev,'(a1,3i5)') '3',(amesh%cell(amesh%QuakeElem(i),j)-1,j=1,3)
  enddo
  close(dev)
!
  
! border elements
  amesh%QuakeBorder_elemNo=0
  do i=1,amesh%QuakeElemNo
     nc=0
     do j=1,3
        do k=1,amesh%QuakeElemNo
           if (amesh%EToE(amesh%QuakeElem(i),j) == amesh%QuakeElem(k).and.k /= i) nc=nc+1
        enddo
     enddo
     if (nc < 3) amesh%QuakeBorder_elemNo=amesh%QuakeBorder_elemNo+1
  enddo
  allocate(amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo))
  allocate(dummy2(amesh%QuakeelemNo,7))
  dummy2=0
  amesh%QuakeBorder_elemNo=0
  do i=1,amesh%QuakeElemNo
     nc=0
     do j=1,3
        before=nc
        do k=1,amesh%QuakeElemNo
           if (amesh%EToE(amesh%QuakeElem(i),j) == amesh%QuakeElem(k).and.k /= i) nc=nc+1
        enddo
        if (before == nc) then
           dummy2(i,1)=dummy2(i,1)+1
           idxx=dummy2(i,1)*2
           if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),j)) then
           do ii=1,3
              do jj=1,3
                 if (amesh%cell(amesh%QuakeElem(i),ii) == amesh%cell(amesh%EToE(amesh%QuakeElem(i),j),jj)) then
                    dummy2(i,idxx)=amesh%cell(amesh%QuakeElem(i),ii)
                    idxx=idxx+1
                 endif
              enddo
           enddo
           else
! case where the elements is on the mesh border
             nv=0
             do ii=1,3
                if (amesh%EToE(amesh%QuakeElem(i),ii) == amesh%QuakeElem(i)) nv=nv+1
             enddo
!      one element mesh (should not happened) : three edges facing the void
             if (nv==3) then
                dummy2(i,1)=3
                dummy2(i,2)=amesh%cell(amesh%QuakeElem(i),1)
                dummy2(i,3)=amesh%cell(amesh%QuakeElem(i),2)
                dummy2(i,4)=amesh%cell(amesh%QuakeElem(i),2)
                dummy2(i,5)=amesh%cell(amesh%QuakeElem(i),3)
                dummy2(i,6)=amesh%cell(amesh%QuakeElem(i),1)
                dummy2(i,7)=amesh%cell(amesh%QuakeElem(i),3)
                exit
             endif
!      two edges facing the void
             if (nv==2) then
                idxxx=1
!           finding the unique neighbour cell
                do ii=1,3
                   if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),ii)) then
                      cj(1)=amesh%EToE(amesh%QuakeElem(i),ii)
                      exit
                   endif
                enddo
!           finding the two common vertices between the reference cell and its unique neighbour
                do ii=1,3
                   do jj=1,3
                      if (amesh%cell(amesh%QuakeElem(i),ii) == amesh%cell(cj(1),jj)) then
                         p(idxxx)=amesh%cell(amesh%QuakeElem(i),ii)
                         idxxx=idxxx+1
                      endif
                   enddo
                enddo
!           finding the complemantary vertice of the reference cell
                do ii=1,3
                   pf=0
                   do jj=1,2
                      if (amesh%cell(amesh%QuakeElem(i),ii) /= p(jj)) pf=pf+1
                   enddo
                   if (pf==0) p(3)=amesh%cell(amesh%QuakeElem(i),ii)
                enddo
!           saving the two edges facing the void
                dummy2(i,1)=2
                dummy2(i,2)=p(1)
                dummy2(i,3)=p(3)
                dummy2(i,4)=p(2)
                dummy2(i,5)=p(3)
                exit
           endif
!      one edges facing the void
           if (nv==1) then
!           finding the two neighbor cells
              idxxx=1
              do ii=1,3
                 if (amesh%QuakeElem(i) /= amesh%EToE(amesh%QuakeElem(i),ii)) then
                    cj(idxxx)=amesh%EToE(amesh%QuakeElem(i),ii)
                    idxxx=idxxx+1
                 endif
              enddo
!          finding the common vertice between the two neighbor cells
              do ii=1,3
                 do jj=1,3
                    if (amesh%cell(cj(1),ii) == amesh%cell(cj(2),jj)) then
                       p(1)=amesh%cell(cj(1),ii)
                    endif
                 enddo
              enddo
!          finding the two complementary vertices which by contruction formed the edge facing the void
              idxxx=2
              do ii=1,3
                 if (amesh%cell(amesh%QuakeElem(i),ii) /= p(1)) then
                    p(idxxx)=amesh%cell(amesh%QuakeElem(i),ii)
                    idxxx=idxxx+1
                 endif
              enddo
              dummy2(i,1)=1
              dummy2(i,2)=p(2)
              dummy2(i,3)=p(3)
           endif
          endif
        endif
     enddo
     if (nc < 3) then
        amesh%QuakeBorder_elemNo=amesh%QuakeBorder_elemNo+1
        amesh%QuakeBorder_elem(amesh%QuakeBorder_elemNo)=amesh%QuakeElem(i)
     endif
  enddo

  dev=23
  open(dev,file='quakeborder.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a9,$)') 'POLYGONS '
  write(dev,*) amesh%QuakeBorder_elemNo,amesh%QuakeBorder_elemNo*4
  do i=1,amesh%QuakeBorder_elemNo
     write(dev,'(a1,3i5)') '3',(amesh%cell(amesh%QuakeBorder_elem(i),j)-1,j=1,3)
  enddo
  close(dev)
!

! border nodes
!print*, dummy2
  allocate(dummy(amesh%QuakeBorder_elemNo*3))
  np=0
  do i=1,amesh%QuakeelemNo
     if (dummy2(i,1) > 0) then
     do j=1,dummy2(i,1)
        if (np==0) then
           np=2
           dummy(1)=dummy2(i,2)
           dummy(2)=dummy2(i,3)
           cycle
        endif
        do k=1,2
           pt=dummy2(i,j*2+k-1)
           found=.false.
           do ii=1,np
              if (pt==dummy(ii)) then
                 found=.true.
                 exit
              endif
           enddo
           if (.not.found) then
              np=np+1
              dummy(np)=pt
           endif
        enddo
     enddo
     endif
  enddo
  amesh%QuakeBorder_NodesNo=np
  allocate(amesh%QuakeBorder_Nodes(np))
  amesh%QuakeBorder_Nodes=dummy(1:np)


  nline=0
  do i=1,amesh%QuakeelemNo
     nline=nline+dummy2(i,1)
  enddo
  dev=23
  open(dev,file='quakeborderline.vtk')
  write(dev,'(a26)') '# vtk DataFile Version 2.0'
  write(dev,'(a7)') 'surface'
  write(dev,'(a5)') 'ASCII'
  write(dev,'(a16)') 'DATASET POLYDATA'
  write(dev,'(a7,$)') 'POINTS '
  write(dev,*) amesh%Nnodes,' float'
  do i=1,amesh%Nnodes
     write(dev,*) amesh%px(i),amesh%py(i),amesh%pz(i)
  enddo
  write(dev,'(a6,$)') 'LINES '
  write(dev,*) nline,3*nline
  do i=1,amesh%QuakeelemNo
     if (dummy2(i,1) /= 0) then
         do j=1,dummy2(i,1)
            write(dev,'(a1,3i5)') '2',(dummy2(i,j*2+k-1)-1,k=1,2)
         enddo
     endif
  enddo
  close(dev)
!
  deallocate(dummy,dummy2)



 
return
end subroutine select_fault_zone_fromfile
!##########################################################################
subroutine  renormalize_slip(amesh,slipout)


!!! To change the name of the subroutine to avoid misleading!!!
real,dimension(:),intent(inout):: slipout
type(mesh) :: amesh
integer :: i
real :: mu_local
character (len=30) :: name_file ! dimension(10) :: name_file

name_file='mu_' // trim(amesh%geo_zone) // '.dat'
name_file=trim(name_file)
open(11,file=name_file)

do i=1,amesh%Ncells
   read(11,*) mu_local
  ! print*, "Before",slipout(i)
   slipout(i)=slipout(i)*amesh%mu/(mu_local*1.E9)
  ! print*, "After",slipout(i)
enddo

close(11)

end subroutine renormalize_slip
!=============================================================================
subroutine normarea(amesh,moment,slipout)
! renormalization by area (to be adapted for different areas)
! main structure
real,dimension(:),intent(inout):: slipout
type(mesh) :: amesh
real :: moment,mu_local
real :: avg_area
!real, allocatable, dimension (:) :: mu_full, mu_quake
real :: mf
integer :: i
character (len=30) :: name_file ! dimension(10) :: name_file

!allocate(mu_full(amesh%Ncells))
!allocate(mu_quake(amesh%QuakeElemNo))

avg_area=0.
do i=1,amesh%QuakeElemNo
   avg_area=avg_area+amesh%area(amesh%QuakeElem(i))
enddo
avg_area=avg_area/amesh%QuakeElemNo

!if (amesh%rigidity==1) then
!   name_file='mu_' // trim(amesh%geo_zone) // '.dat'
!   name_file=trim(name_file)
!   open(11,file=name_file)
!   do i=1,amesh%Ncells
!      read(11,*) mu_full(i)
!      mu_full(i)=mu_full(i)*1.E9
!   enddo
!   close(11)
!   do i=1,amesh%QuakeElemNo
!      mu_quake(i)=mu_full(amesh%QuakeElem(i))
!   enddo
!else
!   mu_full=amesh%mu
!   mu_quake=amesh%mu
!endif

name_file='mu_' // trim(amesh%geo_zone) // '.dat'
name_file=trim(name_file)
open(11,file=name_file)
mu_local=amesh%mu/1.E9

mf=0
do i=1,amesh%Ncells
   if(amesh%rigidity==1) read(11,*) mu_local
   slipout(i)=slipout(i)*amesh%area(i)/avg_area
   mf=mf+mu_local*1.E9*slipout(i)*amesh%area(i)
enddo
close(11)

print*, "momento=",moment
print*, "mf=", mf
do i=1,amesh%Ncells
   slipout(i)=slipout(i)/mf*moment
enddo

!mf=0
!do i=1,amesh%QuakeElemNo
   !slipout(amesh%QuakeElem(i))=slipout(amesh%QuakeElem(i))*amesh%area(amesh%QuakeElem(i))/avg_area
!      mf=mf+mu_full(amesh%QuakeElem(i))*slipout(amesh%QuakeElem(i))*amesh%area(amesh%QuakeElem(i))
!      enddo
!print*, "mf_update=", mf

!deallocate(mu_full,mu_quake)
end subroutine normarea
!==============================================================================
subroutine intersect(A,B,v)
!  Find the common values to both A and B
!
!
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: v
integer, allocatable,dimension(:) :: dummy
integer,allocatable,dimension(:) :: haveit

integer :: i, inx

allocate(dummy(size(A)))
inx = 0
do  i = 1,size(A)
       if (allocated(haveit)) deallocate(haveit)

       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              inx = inx +1
              dummy(inx) = A(i)
       endif
enddo

allocate(v(inx))
v = dummy(1:inx)
deallocate(dummy,haveit)

return
end subroutine intersect
!==============================================================================
!==============================================================================
subroutine isamember(A,B,idx)
!  Check if ameshents of A are in B
!      idx(i) = 1 ameshent of A is present in B
!      idx(i) = 0 ameshent of A is not present in B
integer,dimension(:),intent(in) :: A, B
integer, allocatable,dimension(:),intent(out) :: idx
integer,allocatable,dimension(:) :: haveit
integer :: i

allocate(idx(size(A)))

do  i = 1,size(A)
       call find(B,'==',A(i),haveit)
       if ( haveit(1) > 0 ) then
              idx(i) = 1
       else
              idx(i) = 0
       endif
enddo

return
end subroutine isamember
!==============================================================================
!==============================================================================
! locate returns an array mask based on what points are inside a set an given ameshent
!
FUNCTION locate2d(npts,ptx,pty,amesh) result(arr)
integer,intent(in) ::npts
REAL, INTENT(IN) :: ptx(npts),pty(npts)
REAL, INTENT(IN) :: amesh(4)
LOGICAL :: arr(npts)
real :: ex_min,ex_max,ey_min,ey_max
integer :: i


ex_min = amesh(1)
ex_max = amesh(2)
ey_min = amesh(3)
ey_max = amesh(4)

arr = .false.
do i = 1,npts
       if ((ptx(i) <= ex_max.and.ptx(i) >= ex_min).and. &
            (pty(i) <= ey_max.and.pty(i) >= ey_min)) then

             arr(i) = .true.
       endif
enddo

END FUNCTION locate2d
!!*******************************************************
!*******************************************************
subroutine reorder(a,n)
!######################################################
! Author : André Herrero
! Contact : andherit@gmail.com, andre.herrero@ingv.it
! Public Domain (CC0 1.0 Universal)
!######################################################
implicit none

integer n
real a(n)
integer i
real mem
logical done

done=.false.
do while (.not.done)
       done=.true.
       do i=1,n-1
              if (a(i).lt.a(i+1)) then
                     mem=a(i)
                     a(i)=a(i+1)
                     a(i+1)=mem
                     done=.false.
              endif
       enddo
enddo

return
end subroutine reorder
!=============================================================================
!########################################################################################
!  ########################################################################################
subroutine interp_lin(x,y,xf,yf)
!function interp_lin(x,y,xf) result(yf)
implicit none
real,dimension(2) :: x,y
real :: c,m,xf,yf

m = (y(2)-y(1))/(x(2)-x(1))
c = y(1)-m*x(1)

yf = m*xf+c

!end function interp_lin
end subroutine interp_lin
!  ########################################################################################
subroutine locate(n,xx,x,ans)
IMPLICIT NONE
integer,intent(in) ::n
REAL, INTENT(IN) :: xx(n)
REAL, INTENT(IN) :: x
INTEGER :: ans
INTEGER :: jl,jm,ju
LOGICAL :: ascnd

ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
              jl=jm
       else
              ju=jm
       end if
end do
if (x == xx(1)) then
       !locate=1
       ans=1
else if (x == xx(n)) then
       !locate=n-1
       ans=n-1
else
       !locate=jl
       ans=jl
end if
END subroutine locate

!***********************************************************
subroutine set_seed(na,seed_stoc,index_string)
implicit none
integer,intent(in) :: na, seed_stoc
character*9,intent(in) :: index_string
integer :: jobid,tot_rnos,iseed
integer,dimension(:),allocatable :: a_seed
character (len=30) :: en_var, file_out
integer :: i,s
logical :: file_exists

! initialization of the random generator
!call srand(iseed)
!call sleepqq(1000) ! pause program for 1 second (cause different random seed to be taken)
!call random_seed(put=a_seed)     ! generate random seed


open(10, file='/dev/urandom', access='stream', form='UNFORMATTED')
read(10) i
close(10)

! another way to make sure random numbers don't overlap between runs
call getenv("PBS_ARRAYID",en_var) ! get environmental variable "test"
!call getenv("test",en_var)
read(en_var,'(I5)') jobid
call random_seed(size=iseed)
allocate(a_seed(1:iseed))
call random_seed(get=a_seed)
tot_rnos = na*2  ! this is an estimate of the total number of random numbers used in the programme
call system_clock(s)
a_seed = abs( mod((s*181)*((i-83)*359), 104729) )
if (seed_stoc==1) then
   file_out='Seed_' // trim(index_string) // '.dat'
   file_out=trim(file_out)
   inquire(FILE=file_out, EXIST=file_exists)
   if (file_exists) then
      open(1000,file=file_out)
      read(1000,*) a_seed
      close(1000)
   endif
endif
!a_seed = 123+(jobid-1)*tot_rnos
!!a_seed(iseed) = 123+(jobid-1)*tot_rnos
call random_seed(put=a_seed)     ! generate random seed
file_out='Seed_' // trim(index_string) // '.dat'
   file_out=trim(file_out)
open(1000,file=file_out)
write(1000,*) a_seed
close(1000)
deallocate(a_seed)

! random seed is choosen based on clock time
!call random_seed     ! generate random seed

end subroutine set_seed
!***********************************************************
!=================================================
subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)

  ! convert geodetic longitude and latitude to UTM, and back
  ! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to
  ! lat/long
  ! a list of UTM zones of the world is available at
  ! www.dmap.co.uk/utmworld.htm
    implicit none
  !!  include "constants.h"

  !
  !-----CAMx v2.03
  !
  !     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
  !
  !     This is a Fortran version of the BASIC program "Transverse
  !     Mercator
  !     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra,
  !     2/94)
  !     Based on algorithm taken from "Map Projections Used by the USGS"
  !     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
  !
  !     Input/Output arguments:
  !
  !        rlon                  Longitude (deg, negative for West)
  !        rlat                  Latitude (deg)
  !        rx                    UTM easting (m)
  !        ry                    UTM northing (m)
  !        UTM_PROJECTION_ZONE   UTM zone
  !        iway                  Conversion type
  !                              ILONGLAT2UTM = geodetic to UTM
  !                              IUTM2LONGLAT = UTM to geodetic
  !
    integer UTM_PROJECTION_ZONE,iway
    logical SUPPRESS_UTM_PROJECTION
!    double precision rx,ry,rlon,rlat
!    double precision,parameter :: pi = 3.14159265358979323846264338327950
!    double precision, parameter :: degrad=pi/180.d0, raddeg=180.d0/pi
!    double precision, parameter :: semimaj=6378206.4d0,semimin=6356583.8d0
!    double precision, parameter :: scfa=0.9996d0
    real rx,ry,rlon,rlat
    real,parameter :: pi = 3.14159265358979323846264338327950
    real, parameter :: degrad=pi/180., raddeg=180./pi
    real, parameter :: semimaj=6378206.4,semimin=6356583.8
    real, parameter :: scfa=0.9996

  ! some extracts about UTM:
  !
  ! There are 60 longitudinal projection zones numbered 1 to 60 starting
  ! at 180Â°W.
  ! Each of these zones is 6 degrees wide, apart from a few exceptions
  ! around Norway and Svalbard.
  ! There are 20 latitudinal zones spanning the latitudes 80Â°S to 84Â°N
  ! and denoted
  ! by the letters C to X, ommitting the letter O.
  ! Each of these is 8 degrees south-north, apart from zone X which is 12
  ! degrees south-north.
  !
  ! To change the UTM zone and the hemisphere in which the
  ! calculations are carried out, need to change the fortran code and
  ! recompile. The UTM zone is described
  ! actually by the central meridian of that zone, i.e. the longitude at
  ! the midpoint of the zone, 3 degrees
  ! from either zone boundary.
  ! To change hemisphere need to change the "north" variable:
  !  - north=0 for northern hemisphere and
  !  - north=10000000 (10000km) for southern hemisphere. values must be in
  !  metres i.e. north=10000000.
  !
  ! Note that the UTM grids are actually Mercators which
  ! employ the standard UTM scale factor 0.9996 and set the
  ! Easting Origin to 500,000;
  ! the Northing origin in the southern
  ! hemisphere is kept at 0 rather than set to 10,000,000
  ! and this gives a uniform scale across the equator if the
  ! normal convention of selecting the Base Latitude (origin)
  ! at the equator (0 deg.) is followed.  Northings are
  ! positive in the northern hemisphere and negative in the
  ! southern hemisphere.
     real, parameter :: north=10000000 ! in southern hemisphere   !north=0.d0
     real, parameter :: east=500000.
     real e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
     real f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
     real rx_save,ry_save,rlon_save,rlat_save
!    double precision, parameter :: north=10000.d3 ! in southern hemisphere   !north=0.d0
!    double precision, parameter :: east=500000.d0
!    double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
!    double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
!    double precision rx_save,ry_save,rlon_save,rlat_save

  ! flag for projection from latitude/longitude to UTM, and back
    integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1


    ! checks if conversion to utm has to be done
    if(SUPPRESS_UTM_PROJECTION) then
      if (iway == ILONGLAT2UTM) then
        rx = rlon
        ry = rlat
      else
        rlon = rx
        rlat = ry
      endif
      return
    endif

  ! save original parameters
    rlon_save = rlon
    rlat_save = rlat
    rx_save = rx
    ry_save = ry

    xx = 0.0
    yy = 0.0
    dlat = 0.0
    dlon = 0.0
!    xx = 0.d0
!    yy = 0.d0
!    dlat = 0.d0
!    dlon = 0.d0

  ! define parameters of reference ellipsoid
    e2=1.0-(semimin/semimaj)**2.0
    e4=e2*e2
    e6=e2*e4
    ep2=e2/(1.-e2)

    if (iway == IUTM2LONGLAT) then
      xx = rx
      yy = ry
    else
      dlon = rlon
      dlat = rlat
    endif
  !
  !----- Set Zone parameters
  !
    zone = dble(UTM_PROJECTION_ZONE)
    ! sets central meridian for this zone
    cm = zone*6.0 - 183.0
    cmr = cm*degrad
  !
  !---- Lat/Lon to UTM conversion
  !
    if (iway == ILONGLAT2UTM) then

    rlon = degrad*dlon
    rlat = degrad*dlat

    delam = dlon - cm
    if (delam < -180.) delam = delam + 360.
    if (delam > 180.) delam = delam - 360.
    delam = delam*degrad

    f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
    f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
    f2 = f2*sin(2.*rlat)
    f3 = 15.*e4/256.*45.*e6/1024.
    f3 = f3*sin(4.*rlat)
    f4 = 35.*e6/3072.
    f4 = f4*sin(6.*rlat)
    rm = semimaj*(f1 - f2 + f3 - f4)
    if (dlat == 90. .or. dlat == -90.) then
      xx = 0.
      yy = scfa*rm
    else
      rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
      t = tan(rlat)**2
      c = ep2*cos(rlat)**2
      a = cos(rlat)*delam

      f1 = (1. - t + c)*a**3/6.
      f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
      f2 = f2*a**5/120.
      xx = scfa*rn*(a + f1 + f2)
      f1 = a**2/2.
      f2 = 5. - t + 9.*c + 4.*c**2
      f2 = f2*a**4/24.
      f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
      f3 = f3*a**6/720.
      yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
    endif
    xx = xx + east
    yy = yy + north

  !
  !---- UTM to Lat/Lon conversion
  !
    else

    xx = xx - east
    yy = yy - north
    e1 = sqrt(1. - e2)
    e1 = (1. - e1)/(1. + e1)
    rm = yy/scfa
    u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
    u = rm/(semimaj*u)

    f1 = 3.*e1/2. - 27.*e1**3./32.
    f1 = f1*sin(2.*u)
    f2 = 21.*e1**2/16. - 55.*e1**4/32.
    f2 = f2*sin(4.*u)
    f3 = 151.*e1**3./96.
    f3 = f3*sin(6.*u)
    rlat1 = u + f1 + f2 + f3
    dlat1 = rlat1*raddeg
    if (dlat1 >= 90. .or. dlat1 <= -90.) then
      dlat1 = dmin1(dlat1,dble(90.) )
      dlat1 = dmax1(dlat1,dble(-90.) )
      dlon = cm
    else
      c1 = ep2*cos(rlat1)**2.
      t1 = tan(rlat1)**2.
      f1 = 1. - e2*sin(rlat1)**2.
      rn1 = semimaj/sqrt(f1)
      r1 = semimaj*(1. - e2)/sqrt(f1**3)
      d = xx/(rn1*scfa)

      f1 = rn1*tan(rlat1)/r1
      f2 = d**2/2.
      f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
      f3 = f3*d**2*d**2/24.
      f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
      f4 = f4*(d**2)**3./720.
      rlat = rlat1 - f1*(f2 - f3 + f4)
      dlat = rlat*raddeg

      f1 = 1. + 2.*t1 + c1
      f1 = f1*d**2*d/6.
      f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
      f2 = f2*(d**2)**2*d/120.
      rlon = cmr + (d - f1 + f2)/cos(rlat1)
      dlon = rlon*raddeg
      if (dlon < -180.) dlon = dlon + 360.
      if (dlon > 180.) dlon = dlon - 360.
    endif
    endif

    if (iway == IUTM2LONGLAT) then
      rlon = dlon
      rlat = dlat
      rx = rx_save
      ry = ry_save
    else
      rx = xx
      ry = yy
      rlon = rlon_save
      rlat = rlat_save
    endif

  end subroutine utm_geo
!=================================================
!  ########################################################################################
subroutine pnpoly(PX,PY,XX,YY,N,INOUT)
!
!        SUBROUTINE PNPOLY
!
!        PURPOSE
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
!
!        USAGE
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
!
!        DESCRIPTION OF THE PARAMETERS
!           PX      - X-COORDINATE OF POINT IN QUESTION.
!           PY      - Y-COORDINATE OF POINT IN QUESTION.
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF
!                     VERTICES OF POLYGON.
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF
!                     VERTICES OF POLYGON.
!           N       - NUMBER OF VERTICES IN THE POLYGON.
!           INOUT   - THE SIGNAL RETURNED:
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.
!
!        REMARKS
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY
!           OPTIONALLY BE INCREASED BY 1.
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
!
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE
!
!        METHOD
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE
!           POINT IS INSIDE OF THE POLYGON.
!
!     ..................................................................
	  implicit none
      real,intent(in) :: px,py
      integer,intent(in) :: n
      real :: XX(N),YY(N) !,X(200),Y(200)
      real,allocatable,dimension(:) :: X,Y
      LOGICAL MX,MY,NX,NY
      integer,intent(out) :: inout
      integer :: i,j!,maxdim
allocate(x(n))
allocate(y(n))

do  i = 1,n
     x(i)=xx(i)-px
     y(i)=yy(i)-py
enddo
inout = -1
do i = 1,n
      j = 1+MOD(I,N)
      mx = x(i).GE.0.0
      nx = x(j).GE.0.0
      my = y(i).GE.0.0
      ny = y(j).GE.0.0
!       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2
!       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3
!       INOUT=-INOUT
!       GO TO 2
!3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
!4     INOUT=0
!       RETURN
!5     INOUT=-INOUT
!2     enddo
!       RETURN



    if(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) then
! !     	continue
              cycle
    else
      	if(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) then
		    if (((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))).lt.0.0) then
!!      		 	continue
                            cycle
    		    elseif (((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))).eq.0.0) then
            	             inout = 0
             	             exit
                  elseif(((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))).gt.0.0) then
      			      inout = -inout
!      			continue
                            cycle
                  endif
       else
              inout = -inout
    	endif

     endif
enddo

return
end subroutine pnpoly
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end module typedef
