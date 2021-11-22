program born
     implicit none
     real,allocatable::b(:,:),polar(:,:),non_polar(:,:),diff(:,:),u(:,:),Res(:,:)
     real(kind=8)::cart_coord(3,3),S_T(3,3),px=0.0,py=0.0,pz=0.0
     integer::i,j,k,n=0,natom,ios
     real::dummyv,vol=0.0
     
     open (99,file='born')
     do 
      read(99,*,iostat=ios)dummyv
      if (ios/=0) then
         exit
      else
         n=n+1
      end if
     end do
     close(99) 
     
     natom=n/3
     allocate(b(n,3))
     allocate(polar(natom,3))
     allocate(non_polar(natom,3))
     allocate(diff(natom,3))
     allocate(u(3,natom))
     allocate(Res(3,natom))

! read the input files, born, coordinates, polar, non-polar

     open (99,file='born',action='read')
     do i=1,n
        read(99,*)(b(i,j),j=1,3)
     end do

     open(100,file='coord',action='read')
     do i=1,3
        read(100,*)(cart_coord(i,j),j=1,3)
     end do

     open(101,file='polar',action='read')
     do i=1,natom      
        read(101,*)(polar(i,j),j=1,3)
     end do

     open(102,file='non-polar',action='read')
     do i=1,natom
         read(102,*)(non_polar(i,j),j=1,3)
     end do

! Calculate volume of cell
     
     vol=cart_coord(1,2)*cart_coord(2,3)*cart_coord(3,1)-cart_coord(1,3)*cart_coord(2,2)*cart_coord(3,1)
     vol=vol-cart_coord(1,1)*cart_coord(2,3)*cart_coord(3,2)+cart_coord(1,3)*cart_coord(2,1)*cart_coord(3,2)
     vol=vol+cart_coord(1,1)*cart_coord(2,2)*cart_coord(3,3)-cart_coord(1,2)*cart_coord(2,1)*cart_coord(3,3)
   
     S_T=TRANSPOSE(cart_coord)
     diff=polar-non_polar
     u=TRANSPOSE(diff)

     do k=1,natom
        do i=1,3
           Res(i,k)=S_T(i,1)*u(1,k)+S_T(i,2)*u(2,k)+S_T(i,3)*u(3,k)
        end do
     end do
     
     i=1
     k=1
     do while (k<(natom+1))
        px=px+b(i,1)*Res(1,k)+b(i,2)*Res(2,k)+b(i,3)*Res(3,k)
        py=py+b(i+1,1)*Res(1,k)+b(i+1,2)*Res(2,k)+b(i+1,3)*Res(3,k)
        pz=pz+b(i+2,1)*Res(1,k)+b(i+2,2)*Res(2,k)+b(i+2,3)*Res(3,k)
        i=i+3
        k=k+1
     end do

! Polarization vector

     write(*,*)'Polarization vector',px,py,pz

! Magnitude of Polarization (in microC/cm^2, multiplying by e/V)

     write(*,*)'Magnitude',SQRT(px*px+py*py+pz*pz)*16.0*100/vol,'microC/cm^2'
     
end program born
