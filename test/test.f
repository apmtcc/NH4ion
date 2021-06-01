      program main
      implicit none
      integer,parameter::natom=5
      integer ios,i,j
      real*8 vab,ref,pes,V,RMSE
      real*8 cart(3,5)
      character*1 atom
      character*25 aa,bb,cc
      open(7,file="input.tt")
      open(20,file="output.tt")

      call PREPOT 
      write(20,*)"    NO.   Ref      pes "
      do j=1,6
         read(7,*)
         read(7,"(a19,f8.3,a24,f8.3,a9)")aa,vab,bb,ref,cc
         read(7,*)
         do i=1,natom
            read(7,*)atom,cart(:,i)
         end do
        call POT(cart,V,RMSE) !--->eV
         pes=V*23.0605d0-(-0.010715d0)     !eV--->kcal/mol shift
         write(20,"(i6,2f9.3)")j,ref,pes
         read(7,*)
      end do
      end
