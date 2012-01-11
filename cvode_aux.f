      subroutine add_fcvfun_usr(ydot,j)
c
c     plug-in specific contributions to rhs
c
      include 'SIZE'
      include 'TOTAL'
      common /scrns/ sum (lx1,ly1,lz1,lelt)
     $              ,dtmp(lx1,ly1,lz1,lelt)
      real ydot(1)

      if(.not.ifvarp(1)) return
      ntotv = nx1*ny1*nz1*nelv

      call make_p0th(ydot,j)

      dd = gamma0 ! CVref/CPref ! Note CVref denotes the inverse CPref
      dd = (dd - 1.)/dd

      xfacr= dd * dp0thdt
      call rzero (dtmp,ntotv)
      call cadd(dtmp,xfacr,ntotv)
      call invcol2(dtmp,vtrans(1,1,1,1,2),ntotv)
      call add2(ydot,dtmp,ntotv)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_unpack_sol(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy(t(1,1,1,1,i-1),y(j),ntot)
         j = j + ntot
      enddo

      p0th = y(j)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_pack_sol(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy (y(j),t(1,1,1,1,i-1),ntot)
         j = j + ntot
      enddo

      y(j) = p0th 

      return
      end
c----------------------------------------------------------------------
      subroutine make_p0th(ydot,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real ydot(1) 
      common /scrns/ sum (lx1,ly1,lz1,lelt)
     $              ,dtmp(lx1,ly1,lz1,lelt)

      real dcon(lx1,ly1,lz1,lelv)
      real div (lx1,ly1,lz1,lelv)
      real dpt1(lx1,ly1,lz1,lelv)
      real dpt2(lx1,ly1,lz1,lelv)
      real dpt3(lx1,ly1,lz1,lelv)

      nxyz = nx1*ny1*nz1
      ntotv = nxyz*nelv

      if (.not. ifvcor) then
         ydot(j) = 0.0
         dp0thdt = ydot(j)
         return
      endif

c      call invcol3(dtmp,vtrans,wave,ntotv)
c      call invcol2(dtmp,vtrans(1,1,1,1,2),ntotv)
      call rone(dtmp,ntotv)

      dd = gamma0 ! CVref/CPref ! Note CVref denotes the inverse CPref
      dd = -1.0*(dd - 1.)/dd
      call cmult(dtmp,dd,ntotv)
      call cadd(dtmp,1.0,ntotv)
      call invcol1(dtmp,ntotv)
      call cmult(dtmp,p0th,ntotv)

      call copy   (sum,ydot,  ntotv)

      call convop (dcon,T    ,ntotv)
      call col2   (dcon, bm1,ntotv)
c      call divws  (dpt3,T,dpt1,nelv,1)		! mesh moving terms should 
c      call sub2   (dcon,dpt3 ,ntotv)		! not be included with cvode

      call dssum  (dcon,nx1,ny1,nz1)
      call col2   (dcon,binvm1,ntotv)
      call add2   (sum ,dcon ,ntotv)
      call invcol2(sum,     T,ntotv)

c add the divergence term here!

      call opdiv   (div,vx,vy,vz)
      call dssum   (div,nx1,ny1,nz1)
      call col2    (div,binvm1,ntotv)
      call sub2    (sum,div,ntotv)

      call col2(sum,dtmp,ntotv)      
      call col2(sum,bm1,ntotv)

      ydot(j) = glsum(sum,ntotv) / volvm1 
      dp0thdt = ydot(j)

      return
      end
 
