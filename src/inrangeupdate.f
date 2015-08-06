      subroutine inrangeupdate(p,numlt,dmat,index,cf)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------

      implicit none

      integer    p,numlt
	  real       dmat(p,numlt)
	  logical    index(numlt)
	  real       cf(p)

      integer i,j

      do i = 1,numlt
	  
	    do j = 1,p

            if ( abs(dmat(j,i)) .ge. cf(j) ) then
              goto 100 
            endif

         enddo
		 
		 index(i) = .TRUE.

 100  continue
      enddo

      return
      end
