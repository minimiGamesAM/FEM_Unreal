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
module main

    contains
    
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
    
    SUBROUTINE invert(matrix)
            !
            ! This subroutine inverts a small square matrix onto itself.
            !
             IMPLICIT NONE
             INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
             REAL(iwp),INTENT(IN OUT)::matrix(:,:)
             REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
             INTEGER::ndim,i,k
             ndim=UBOUND(matrix,1)
             IF(ndim==2)THEN
               det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
               j11=matrix(1,1)
               matrix(1,1)=matrix(2,2)
               matrix(2,2)=j11
               matrix(1,2)=-matrix(1,2)
               matrix(2,1)=-matrix(2,1)
               matrix=matrix/det
             ELSE IF(ndim==3)THEN
               det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
               det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
               det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
               j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
               j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
               j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
               j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
               j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
               j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
               j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
               j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
               j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
               matrix(1,1)=j11
               matrix(1,2)=j12
               matrix(1,3)=j13
               matrix(2,1)=j21
               matrix(2,2)=j22
               matrix(2,3)=j23
               matrix(3,1)=j31
               matrix(3,2)=j32
               matrix(3,3)=j33
               matrix=matrix/det
             ELSE
               DO k=1,ndim
                 con=matrix(k,k)
                 matrix(k,k)=1.0_iwp
                 matrix(k,:)=matrix(k,:)/con
                 DO i=1,ndim
                   IF(i/=k)THEN
                     con=matrix(i,k)
                     matrix(i,k)=0.0_iwp
                     matrix(i,:)=matrix(i,:)-matrix(k,:)*con
                   END IF
                 END DO
               END DO
             END IF
            RETURN
    END SUBROUTINE invert
    
    SUBROUTINE getname(argv,nlen)
        !
        ! This subroutine reads the base name of data file.
        !
             IMPLICIT NONE
             INTEGER::narg
             INTEGER,INTENT(OUT)::nlen
             INTEGER::lnblnk,iargc
             CHARACTER(*),INTENT(OUT)::argv
             LOGICAL found
             narg=IARGC()
             IF(narg<1)THEN
               WRITE(*,*)'Please enter the base name of data file: '
               READ(*,*)argv
              ELSE
               CALL getarg(1,argv)
             ENDIF
             nlen=LNBLNK(argv)
             INQUIRE(FILE=argv(1:nlen)//'.dat',EXIST=found)
             IF(.NOT.found)THEN
              WRITE(*,*)'Data file not found: ',argv(1:nlen)//'.dat'
              WRITE(*,*)'Please create or check spelling.'
              STOP
             ENDIF
            RETURN
    END SUBROUTINE getname
    
    SUBROUTINE mesh(g_coord,g_num,argv,nlen,ips)
        !
        ! This subroutine produces a PostScript output file "*.msh" displaying
        ! the undeformed finite element mesh.
        !
         IMPLICIT NONE
         INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
         REAL(iwp),INTENT(IN)::g_coord(:,:)
         INTEGER,INTENT(IN)::g_num(:,:),ips,nlen
         CHARACTER(*),INTENT(IN)::argv
         REAL(iwp)::xmin,xmax,ymin,ymax,width,height,scale=72,sxy,xo,yo,x,y,      &
           pt5=0.5_iwp,opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,         &
           d11=11.0_iwp
         INTEGER::i,ii,j,jj,nn,nod,nel
         OPEN(ips,FILE=argv(1:nlen)//'.msh')
        !
        !                       compute size of mesh
        !
         nn=UBOUND(g_coord,2)
         xmin=g_coord(1,1)
         xmax=g_coord(1,1)
         ymin=g_coord(2,1)
         ymax=g_coord(2,1)
         DO i=2,nn
           IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)      
           IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)      
           IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)      
           IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)      
         END DO
         width =xmax-xmin
         height=ymax-ymin
        !
        !                       allow 1.5" margin minimum on each side of figure
        !
         IF(height.GE.d11/ept5*width)THEN
        !
        !                       height governs the scale
        !
           sxy=scale*d8/height
           xo=scale*pt5*(ept5-d8*width/height)
           yo=scale*opt5
         ELSE
        !
        !                       width governs the scale
        !
           sxy=scale*fpt5/width
           xo=scale*opt5
           yo=scale*pt5*(d11-fpt5*height/width)
         END IF
        !
        !                       start PostScript output
        !
         WRITE(ips,'(a)')'%!PS-Adobe-1.0'
         WRITE(ips,'(a)')'%%DocumentFonts: none'
         WRITE(ips,'(a)')'%%Pages: 1'
         WRITE(ips,'(a)')'%%EndComments'
         WRITE(ips,'(a)')'/m {moveto} def'
         WRITE(ips,'(a)')'/l {lineto} def'
         WRITE(ips,'(a)')'/s {stroke} def'
         WRITE(ips,'(a)')'/c {closepath} def'
         WRITE(ips,'(a)')'%%EndProlog'
         WRITE(ips,'(a)')'%%Page: 0 1'
         WRITE(ips,'(a)')'gsave'
         WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
         WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
        !
        !                       draw the mesh
        !
         nod=UBOUND(g_num,1)
         nel=UBOUND(g_num,2)
         IF(nod==5)nod=4
         IF(nod==9)nod=8
         IF(nod==10)nod=9
         IF(nod==15)nod=12
         DO i=1,nel
           ii=g_num(1,i)
           IF(ii==0)CYCLE
           x=sxy*(g_coord(1,ii)-xmin)
           y=sxy*(g_coord(2,ii)-ymin)
           WRITE(ips,'(2f9.2,a)')x,y,' m'
           DO j=2,nod
             jj=g_num(j,i)
             x=sxy*(g_coord(1,jj)-xmin)
             y=sxy*(g_coord(2,jj)-ymin)
             WRITE(ips,'(2f9.2,a)') x, y,' l'
           END DO
           WRITE(ips,'(a)')'c s'
         END DO
        !
        !                       close output file
        !
         WRITE(ips,'(a)')'grestore'
         WRITE(ips,'(a)')'showpage'
         CLOSE(ips)
        !
        RETURN
    END SUBROUTINE mesh
    
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
    
    SUBROUTINE num_to_g(num,nf,g)
        !
        ! This subroutine finds the g vector from num and nf.
        !
         IMPLICIT NONE
         INTEGER,INTENT(IN)::num(:),nf(:,:)  
         INTEGER,INTENT(OUT)::g(:)
         INTEGER::i,k,nod,nodof 
         nod=UBOUND(num,1) 
         nodof=UBOUND(nf,1)
         DO i=1,nod
           k=i*nodof
           g(k-nodof+1:k)=nf(:,num(i))
         END DO
        RETURN
    END SUBROUTINE num_to_g  
    
    SUBROUTINE fkdiag(kdiag,g)
        !
        ! This subroutine computes the skyline profile.
        !
         IMPLICIT NONE
         INTEGER,INTENT(IN)::g(:)
         INTEGER,INTENT(OUT)::kdiag(:)
         INTEGER::idof,i,iwp1,j,im,k
         idof=SIZE(g)
         DO i=1,idof
           iwp1=1
           IF(g(i)/=0)THEN
             DO j=1,idof
               IF(g(j)/=0)THEN
                 im=g(i)-g(j)+1
                 IF(im>iwp1)iwp1=im
               END IF
             END DO
             k=g(i)
             IF(iwp1>kdiag(k))kdiag(k)=iwp1
           END IF
         END DO
        RETURN
    END SUBROUTINE fkdiag
    
end module main   

module geom

    contains
    
    SUBROUTINE formnf(nf)
        !
        ! This subroutine forms the nf matrix.
        !
         IMPLICIT NONE
         INTEGER,INTENT(IN OUT)::nf(:,:)
         INTEGER::i,j,m
         m=0
         DO j=1,UBOUND(nf,2)
           DO i=1,UBOUND(nf,1)
             IF(nf(i,j)/=0)THEN
               m=m+1
               nf(i,j)=m
             END IF
           END DO
         END DO
        RETURN
    END SUBROUTINE formnf
    
end module geom
      
    FUNCTION determinant(jac)RESULT(det)
           !
           ! This function returns the determinant of a 1x1, 2x2 or 3x3
           ! Jacobian matrix.
           !
            IMPLICIT NONE    
            INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
            REAL(iwp),INTENT(IN)::jac(:,:)
            REAL(iwp)::det
            INTEGER::it 
            it=UBOUND(jac,1)  
            SELECT CASE(it)
            CASE(1)
              det=1.0_iwp
            CASE(2)
              det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
            CASE(3)
              det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
              det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
              det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
            CASE DEFAULT
              WRITE(*,*)' wrong dimension for Jacobian matrix'
            END SELECT
           RETURN
    END FUNCTION determinant
    
program FortrantClient

use main
use geom

implicit none
      
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndim,ndof,nels,neq,nip,nlen,&
      nn,nod,nodof,nprops=3,np_types,nr,nst 
    REAL(iwp)::det,penalty=1.0e20_iwp,zero=0.0_iwp
    CHARACTER(len=15)::argv,element
    !-----------------------dynamic arrays------------------------------------
    INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),kdiag(:),nf(:,:), &
      no(:),node(:),num(:),sense(:)    
    REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),dee(:,:),der(:,:),deriv(:,:), &
      eld(:),fun(:),gc(:),gravlo(:),g_coord(:,:),jac(:,:),km(:,:),kv(:),     &
      loads(:),points(:,:),prop(:,:),sigma(:),value(:),weights(:)  
    !-----------------------input and initialisation--------------------------
    CALL getname(argv,nlen)
    OPEN(10,FILE=argv(1:nlen)//'.dat') 
    OPEN(11,FILE=argv(1:nlen)//'.res')
    READ(10,*)element,nod,nels,nn,nip,nodof,nst,ndim,np_types 
    ndof=nod*nodof
    ALLOCATE(nf(nodof,nn),points(nip,ndim),dee(nst,nst),g_coord(ndim,nn),    &
      coord(nod,ndim),jac(ndim,ndim),weights(nip),num(nod),g_num(nod,nels),  &
      der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),km(ndof,ndof),eld(ndof),   &
      sigma(nst),g(ndof),g_g(ndof,nels),gc(ndim),fun(nod),etype(nels),       &
      prop(nprops,np_types))
    READ(10,*)prop 
    etype=1 
    IF(np_types>1)READ(10,*)etype
    READ(10,*)g_coord 
    READ(10,*)g_num
    IF(ndim==2)CALL mesh(g_coord,g_num,argv,nlen,12)
    nf=1 
    READ(10,*)nr,(k,nf(:,k),i=1,nr) 
    CALL formnf(nf) 
    neq=MAXVAL(nf) 
    ALLOCATE(kdiag(neq),loads(0:neq),gravlo(0:neq)) 
    kdiag=0
    !-----------------------loop the elements to find global arrays sizes-----
    elements_1: DO iel=1,nels
      num=g_num(:,iel) 
      CALL num_to_g(num,nf,g) 
      g_g(:,iel)=g
      CALL fkdiag(kdiag,g)
    END DO elements_1
    DO i=2,neq 
      kdiag(i)=kdiag(i)+kdiag(i-1) 
    END DO 
    ALLOCATE(kv(kdiag(neq)))
    WRITE(11,'(2(A,I5))')                                                    &
      " There are",neq," equations and the skyline storage is",kdiag(neq)
    !-----------------------element stiffness integration and assembly--------
 
 
    !INTEGER, PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !
    !INTEGER, PARAMETER::dims = 3
    !INTEGER, PARAMETER::nbPoints = 4
    !    
    !INTEGER res, i, j
    !REAL(iwp)::det
    !REAL(iwp) determinant ! para llamar la funcion.
    !REAL(iwp) coord(nbPoints, dims)
    !REAL(iwp) der(dims, nbPoints)
    !REAL(iwp) jac(dims, dims)
    !REAL(iwp) deriv(dims, nbPoints)
    !REAL(iwp) bee(6, nbPoints * 3)
    !REAL(iwp) dee(6, 6)
    !REAL(iwp) btdb( nbPoints * 3, nbPoints * 3)
    !REAL(iwp)::v = 0.3_iwp, e = 10000.0_iwp
    !  
    !coord(1, : ) = (/ 0, 0, 1 /)
    !coord(2, : ) = (/ 1, 0, 0 /)
    !coord(3, : ) = (/ 1, 1, 0 /)
    !coord(4, : ) = (/ 0, 0, 0 /)
    !
    !der(1, : ) = (/ 0, 1, 0, -1 /)
    !der(2, : ) = (/ 0,-1, 1, 0  /)
    !der(3, : ) = (/ 1, 0, 0, -1 /)
    !
    !jac = MATMUL(der, coord)
    !
    !det = determinant(jac)
    !
    !call invert(jac)
    !  
    !deriv = MATMUL(jac, der)
    ! 
    !call beemat(bee, deriv)
    !call deemat(dee, e, v)
    !
    !btdb =  MATMUL(MATMUL(TRANSPOSE(bee), dee), bee)
    !
    !!integration
    !btdb = btdb * det / 6
    !
    !do i = LBOUND (btdb, 1), UBOUND (btdb, 1)
    !  do j = LBOUND (btdb, 2), UBOUND (btdb, 2)
    !     Print *, btdb(i, j)
    !     !write(*, '(f9.1)') bee(i, j)
    !  end do
    !end do
    
    pause
    
end program FortrantClient