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
       
       SUBROUTINE deemat(dee, e, v)
        !
        ! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
        ! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
        ! (three dimensions).
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
         REAL(iwp),INTENT(IN)::e,v
         REAL(iwp),INTENT(OUT)::dee(:,:)
         REAL(iwp)::v1, v2, c, vv, zero = 0.0_iwp, pt5 = 0.5_iwp, one = 1.0_iwp, two = 2.0_iwp
         INTEGER::i,ih
         dee = zero  
         ih = UBOUND(dee, 1)
         v1= one - v
         c = e / ((one + v)*(one - two * v))
         SELECT CASE(ih)
         CASE(6)
           v2=v/(one-v)
           vv=(one-two*v)/(one-v)*pt5
           DO i=1,3
             dee(i,i)=one
           END DO
           DO i=4,6
             dee(i,i)=vv
           END DO
           dee(1,2)=v2
           dee(2,1)=v2
           dee(1,3)=v2
           dee(3,1)=v2
           dee(2,3)=v2
           dee(3,2)=v2
           dee=dee*e/(two*(one+v)*vv)
         CASE DEFAULT
           WRITE(*,*)'wrong size for dee matrix'
         END SELECT
        RETURN
       END SUBROUTINE deemat
       
   
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
    REAL(iwp) dee(6, 6)
    REAL(iwp) btdb( nbPoints * 3, nbPoints * 3)
    REAL(iwp)::v = 0.3_iwp, e = 10000.0_iwp
        
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
    call deemat(dee, e, v)
    
    btdb =  MATMUL(MATMUL(TRANSPOSE(bee), dee), bee)
    
    do i = LBOUND (btdb, 1), UBOUND (btdb, 1)
      do j = LBOUND (btdb, 2), UBOUND (btdb, 2)
         Print *, btdb(i, j)
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

