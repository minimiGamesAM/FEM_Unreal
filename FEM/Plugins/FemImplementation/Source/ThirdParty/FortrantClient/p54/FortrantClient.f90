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

program FortrantClient

USE main
USE geom

implicit none
    
    !ndim       = number of dimensions
    !ndof       = number of degree of freedom per element
    !nels       = total number of elements
    !neq        = number of equations (total number of non-zero freedoms)
    !nip        = number of intregation points per element
    !neq        = number of degree of freedom in the mesh
    !nn         = total number of nodes in the problem
    !nod        = number of node per element
    !nodof      = number of freedoms per node (x, y, z, q1, q2, q3 etc)
    !nprops     = number of material properties
    !nr         = number of restrained nodes (puede ser por lo menos en uno de los freedoms)
    !nst        = number of stress / strain terms
    !gc         = Integrating point coordinates
    !sigma      = stress terms
    !x0         = old displacements
    !x1         = new displacements
    !d1x0       = old velocity
    !d1x1       = new velocity
    !d2x0       = old acceleration
    !d2x1       = new acceleration
    !theta      = time integration weighting parameter
    !fk         = Rayleigh damping parameter on stiffness
    !fm         = Rayleigh damping parameter on mass
    !nres       = node number at witch time history is to be printed
    !val        = applied nodal load weightings
    !npri       = output printed every npri time step
    
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER::fixed_freedoms,i,j,iel,k,loaded_nodes,ndim,ndof,nels,neq,nip,nlen,&
      nn,nod,nodof,npri,nprops=3,np_types,nr,nst,nstep,nres 
    REAL(iwp)::det,penalty=1.0e20_iwp,zero=0.0_iwp,c1,c2,c3,c4,theta,time,dtim,fm,fk 
    CHARACTER(len=15)::argv,element
    !-----------------------dynamic arrays------------------------------------
    !num        = element node number vector
    INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),kdiag(:),nf(:,:), &
      no(:),node(:),num(:),sense(:)    
    REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),dee(:,:),der(:,:),deriv(:,:),ecm(:,:), &
      d1x0(:),d1x1(:),d2x0(:),d2x1(:),eld(:),fun(:),f1(:),gc(:),gravlo(:),g_coord(:,:),jac(:,:),km(:,:),kv(:),     &
      loads(:),mm(:,:),mv(:),points(:,:),prop(:,:),sigma(:),value(:),val(:,:),weights(:), &
      x0(:), x1(:)  
    
   !! REAL(iwp) determinant ! para llamar la funcion. ?????? no estaba en el codigo
    
    !-----------------------input and initialisation--------------------------
    
    !CALL getname(argv,nlen)
    argv = 'p54_2'
    nlen = 5
    
    OPEN(10,FILE=argv(1:nlen)//'.dat') 
    OPEN(11,FILE=argv(1:nlen)//'.res')
    
    !prop       = material property (e, v, gamma)
    !nodef      = nodes to be fixed
    
    READ(10,*) element, nod, nels, nn, nip, nodof, nst, ndim, np_types 
    ndof=nod*nodof
    
    !nf         = nodal freedom array (nodof rows and nn colums)
    !km         = element stiffness matrix
    !eld        = element displacement vector
    !g          = element steering vector
    !g_g        = global element steering matrix (ndof "number of deg of freedom per element", nels "number of elements")
    ALLOCATE( nf(nodof, nn), points(nip, ndim), dee(nst, nst), g_coord(ndim, nn),    &
      coord(nod, ndim), jac(ndim, ndim), weights(nip), num(nod), g_num(nod, nels),  &
      der(ndim, nod), deriv(ndim, nod), bee(nst, ndof), km(ndof, ndof), eld(ndof),   &
      sigma(nst), g(ndof), g_g(ndof, nels), mm(ndof,ndof), ecm(ndof,ndof), gc(ndim), fun(nod), etype(nels),       &
      prop(nprops, np_types))
    
    READ(10,*)prop
    
    etype=1
    
    IF(np_types>1)READ(10,*)etype
    
    READ(10,*)g_coord 
    READ(10,*)g_num
    
    IF(ndim==2)CALL mesh(g_coord,g_num,argv,nlen,12)
    
    nf=1 
    
    READ(10,*)nr,(k,nf(:,k),i=1,nr) 
    CALL formnf(nf) 
    
    neq = MAXVAL(nf) 
    ALLOCATE(x0(0:neq),d1x0(0:neq),x1(0:neq),d2x0(0:neq), &
        d1x1(0:neq),d2x1(0:neq),kdiag(neq), loads(0 : neq), gravlo(0 : neq)) 
    
    loaded_nodes = 2
    nres = 6
    ALLOCATE(node(loaded_nodes),val(loaded_nodes,ndim))
    val(1, :) = 0.25
    val(2, :) = 0.25
       
    node(1) = nres
    node(2) = 2
    
    kdiag=0
    
    !-----------------------loop the elements to find global arrays sizes-----
    elements_1: DO iel=1,nels
      !get element node number vector
      num=g_num(:,iel) 
      CALL num_to_g(num,nf,g) 
      g_g(:, iel) = g
      CALL fkdiag(kdiag,g)
    END DO elements_1
    DO i=2,neq 
      kdiag(i)=kdiag(i)+kdiag(i-1) 
    END DO 
    ALLOCATE(kv(kdiag(neq)), mv(kdiag(neq)),f1(kdiag(neq)))
    WRITE(11,'(2(A,I5))')                                                    &
      " There are",neq," equations and the skyline storage is",kdiag(neq)
    
    !-----------------------element stiffness integration and assembly--------
    CALL sample(element,points,weights)
    kv=zero 
    gravlo=zero
    elements_2: DO  iel=1,nels
      CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel))) 
      num=g_num(:,iel)
      coord=TRANSPOSE(g_coord(:,num)) 
      g=g_g(:,iel) 
      km=zero
      mm=zero
      eld=zero
      int_pts_1: DO i=1,nip
        CALL shape_fun(fun,points,i) 
        CALL shape_der(der,points,i)
        jac=MATMUL(der,coord) 
        det=determinant(jac) 
        CALL invert(jac)
        deriv=MATMUL(jac,der) 
        CALL beemat(bee,deriv)
        km=km+MATMUL(MATMUL(transpose(bee),dee),bee)*det*weights(i)
        CALL ecmat(ecm,fun,ndof,nodof)
        mm=mm+ecm*det*weights(i)*prop(3,etype(iel))
        eld(nodof:ndof:nodof)=eld(nodof:ndof:nodof)+fun(:)*det*weights(i)
      END DO int_pts_1
      CALL fsparv(kv,km,g,kdiag)
      CALL fsparv(mv,mm,g,kdiag)
      gravlo(g)=gravlo(g)-eld*prop(3,etype(iel))
    END DO elements_2
    
    !-----------------------initial conditions and factorise equations--------
    
    dtim = 1.0
    theta = 0.5
    fm = 0.005
    fk = 0.272
        
    x0=0.0
    d1x0=0.0
    d2x0=0.0
    c1=(1-theta)*dtim
    c2=fk-c1
    c3=fm+1.0/(theta*dtim)
    c4=fk+theta*dtim
    f1=c3*mv+c4*kv
    CALL sparin(f1,kdiag)
    time=zero
    
    !-----------------------time stepping loop--------------------------------
        
    WRITE(11,'(/A,I5)')" Result at node",nres
    WRITE(11,'(A)')"    time        load        x-disp      y-disp      z-disp"
    WRITE(11,'(5E12.4)')time,load(time),x0(nf(:,nres))
    
    nstep = 20
    npri = 1
    !timesteps: DO j=1,nstep
        
    timesteps: DO j=1,nstep
        time=time+dtim
        loads=zero
        x1=c3*x0+d1x0/theta
        DO i=1,loaded_nodes
          loads(nf(:,node(i)))=                                                &
            val(i,:)*(theta*dtim*load(time)+c1*load(time-dtim))
        END DO
        CALL linmul_sky(mv,x1,d1x1,kdiag)
        d1x1=loads+d1x1
        loads=c2*x0
        CALL linmul_sky(kv,loads,x1,kdiag)
        x1=x1+d1x1
        CALL spabac(f1,x1,kdiag)
        d1x1=(x1-x0)/(theta*dtim)-d1x0*(1.0-theta)/theta
        d2x1=(d1x1-d1x0)/(theta*dtim)-d2x0*(1.0-theta)/theta
        x0=x1
        d1x0=d1x1
        d2x0=d2x1
        IF(j/npri*npri==j)WRITE(11,'(5E12.4)')time,load(time),x0(nf(:,nres))
    END DO timesteps
    
    
    
    
    
    
    
    
    
    
    !loads=zero 
    !READ(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
    !loads=loads+gravlo 
    !READ(10,*)fixed_freedoms
    !IF(fixed_freedoms/=0)THEN
    !  ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),                   &
    !    value(fixed_freedoms),no(fixed_freedoms))
    !  READ(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
    !  DO  i=1,fixed_freedoms 
    !    no(i)=nf(sense(i),node(i)) 
    !  END DO 
    !  kv(kdiag(no))=kv(kdiag(no))+penalty 
    !  loads(no)=kv(kdiag(no))*value
    !END IF
    !
    !!-----------------------equation solution---------------------------------
    !CALL sparin(kv,kdiag) 
    !CALL spabac(kv,loads,kdiag) 
    !loads(0)=zero
    !
    !IF(ndim==3)THEN 
    !  WRITE(11,'(/A)')"  Node   x-disp      y-disp      z-disp"
    !ELSE 
    !  WRITE(11,'(/A)')"  Node   x-disp      y-disp"
    !END IF
    !DO k=1,nn 
    !  WRITE(11,'(I5,3E12.4)')k,loads(nf(:,k)) 
    !END DO
    !
    !!-----------------------recover stresses at element Gauss-points----------
    !
    !elements_3: DO iel=1,nels
    !    CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel))) 
    !    num=g_num(:,iel)
    !    coord=TRANSPOSE(g_coord(:,num)) 
    !    g=g_g(:,iel) 
    !    eld=loads(g)
    !    int_pts_2: DO i=1,nip
    !      CALL shape_der(der,points,i) 
    !      CALL shape_fun(fun,points,i)
    !      gc=MATMUL(fun,coord) 
    !      jac=MATMUL(der,coord) 
    !      CALL invert(jac)
    !      deriv=MATMUL(jac,der) 
    !      CALL beemat(bee,deriv)
    !      sigma=MATMUL(dee,MATMUL(bee,eld))
    !      IF(ndim==3)THEN 
    !        WRITE(11,'(I8,4X,3E12.4)')iel,gc
    !        WRITE(11,'(6E12.4)')sigma
    !      ELSE 
    !        WRITE(11,'(I8,2E12.4,5X,3E12.4)')iel,gc,sigma
    !      END IF
    !    END DO int_pts_2
    !END DO elements_3
    
    !-----------------------Load-time function--------------------------------
 CONTAINS
FUNCTION load(t) RESULT(load_result)
    IMPLICIT NONE     
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::t
    REAL(iwp)::load_result
    load_result=COS(0.3_iwp*t)
    RETURN
 END FUNCTION load
    
    !----------------------------------------------------------------------- 
    ! Probando la notacion eld(1:12:4) = 2 --> significa llenar de numero 2 desde 1 a 12 saltando 4 por ejemplo
    !INTEGER, PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !INTEGER i
    !
    !REAL(iwp)::zero = 0.0_iwp
    !REAL(iwp),ALLOCATABLE::eld(:)
    !
    !ALLOCATE(eld(12))
    !
    !eld = zero
    !
    !eld(1:12:4) = 2
    !
    !do i = LBOUND (eld, 1), UBOUND (eld, 1)
    !  Print *, eld(i)
    !end do
    !-------------------------------------------------------------------------

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
    
end program FortrantClient