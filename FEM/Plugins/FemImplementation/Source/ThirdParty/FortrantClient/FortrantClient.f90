!  FortrantClient.f90 
!
!  FUNCTIONS:
!  FortrantClient - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: FortrantClient
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    module mathLibFEM  
        
    contains      
       SUBROUTINE beemat(bee, deriv)
            !
            ! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
             REAL(iwp),INTENT(IN)::deriv(:,:)
             REAL(iwp),INTENT(OUT)::bee(:,:)
             INTEGER::k,l,m,n,ih,nod
             REAL::x,y,z
             bee=0.0_iwp
             ih=UBOUND(bee,1)
             nod=UBOUND(deriv,2)
             SELECT CASE (ih)
             CASE(3,4)
               DO m=1,nod
                 k=2*m
                 l=k-1
                 x=deriv(1,m)
                 y=deriv(2,m)
                 bee(1,l)=x
                 bee(3,k)=x
                 bee(2,k)=y
                 bee(3,l)=y
               END DO
             CASE(6)
               DO m = 1, nod
                 n = 3 * m
                 k = n - 1
                 l = k - 1
                 x = deriv(1, m)
                 y = deriv(2, m)
                 z = deriv(3, m)
                 bee(1, l) = x
                 bee(4, k) = x
                 bee(6, n) = x
                 bee(2, k) = y
                 bee(4, l) = y
                 bee(5, n) = y
                 bee(3, n) = z
                 bee(5, k) = z
                 bee(6, l) = z
               END DO
             CASE DEFAULT
               WRITE(*,*)'wrong dimension for nst in bee matrix'        
             END SELECT   
            RETURN
            END SUBROUTINE beemat
   
    end module mathLibFEM 
    
program FortrantClient

use mathLibFEM  
implicit none
    
    
    INTEGER, PARAMETER::iwp=SELECTED_REAL_KIND(15)
    
    INTEGER, PARAMETER::dims = 3
    INTEGER, PARAMETER::nbPoints = 4
        
    INTEGER res, i, j
    REAL(iwp) deriv(dims, nbPoints)
    REAL(iwp) bee(6, nbPoints * 3)
    !REAL ARRAY_B (2:8, -3:20)
    
    !float der[dim * nbPoints] = {
    !    0, 1, 0, -1,
    !    0,-1, 1, 0,
    !    1, 0, 0, -1
    !};
    
    deriv(1, : ) = (/ 0, 1, 0, -1 /)
    deriv(2, : ) = (/ 0,-1, 1, 0 /)
    deriv(3, : ) = (/ 1, 0, 0, -1 /)
    
    call beemat(bee, deriv)
        
    do i = LBOUND (bee, 1), UBOUND (bee, 1)
      do j = LBOUND (bee, 2), UBOUND (bee, 2)
         Print *, bee(i, j)
         !write(*, '(f9.1)') bee(i, j)
      end do
    end do
    !
    !
    !res2 = LBOUND (ARRAY_B, 2)
    !
    !Print *, res2
    
    !call show_consts(res2)

    pause
    
end program FortrantClient

