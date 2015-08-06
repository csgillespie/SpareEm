      subroutine distintri(nrc,np,x,cf,nlt,nm,rc,rowinds,colinds,dmat)
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------

      implicit none

	  integer   nrc,np
      real      x(nrc,np)
      real      cf
	  integer   nlt,nm,rc
      integer   rowinds(nm)
      integer   colinds(nm)
	  real      dmat(np,nm)      

      integer i,j,k
	  real    sm
C --- ------------------------------------------------------------------
C --- Do Something
C --- ------------------------------------------------------------------

C --- The diagonal is always true - so don't check it.
C --- The matrices are symmetric, so we're only checking half.

      nlt = 0
	  rc = 0

      do i = 1,nrc-1
	  do j = i+1,nrc

         if ( nlt > nm ) then
		   rc = -i
		   goto 100
	     endif

C     write(*,*)'checking ',i,j
		 
		 sm = 0
         do k = 1,np

            dmat(k,nlt+1) = abs(x(i,k) - x(j,k))
			sm = sm + dmat(k,nlt+1)
            if ( sm > cf ) then
C             write(*,*)'getting out at ',i,j,k
              goto 100 
            endif

         enddo

         nlt          = nlt + 1
         rowinds(nlt) = i
         colinds(nlt) = j

 100  continue
      enddo
      enddo

C     write(*,*)'The rows are',rowinds(1:nlt)
C     write(*,*)'The cols are',colinds(1:nlt)

      return
      end
