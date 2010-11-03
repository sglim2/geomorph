c -----------------------------------------------------------------------------
*dk midpt
      subroutine midpt(x,x1,x2,nn)
 
c...  This routine finds the midpoint x along the shorter great circle
c...  arc between points x1 and x2 on the unit sphere.
 
      real x(nn,3), x1(nn,3), x2(nn,3)
 
      do j=1,3
         x(1,j) = x1(1,j) + x2(1,j)
      end do
 
      xnorm = 1./sqrt(x(1,1)**2 + x(1,2)**2 + x(1,3)**2)
 
      do j=1,3
         x(1,j) = xnorm*x(1,j)
      end do
 
      end

c -----------------------------------------------------------------------------
*dk grdgen
      subroutine grdgen(xn,nt)
      implicit none
c...  This routine generates the nodal coordinates xn for an
c...  icosahedral grid on the unit sphere.  The grid resolution
c...  corresponds to a subdivision of the edges of the original
c...  icosahedral triangles into nt equal parts.

c...  Note: nt is input as mt

c      real xn(nt+1,nt+1,10,3)
      integer, intent(inout) :: nt
c      real, allocatable :: xn(:,:,:,:) 
      real, intent(inout) :: xn((nt+1),(nt+1),10,3)

      real     :: fifthpi, w, cosw, sinw, lvt, nn, phi
      integer  :: id, k, i1, i2, j1, j2, l, l2, m, sgn
 
      fifthpi = 0.4*asin(1.)
       w       = 2.0*acos(1./(2.*sin(fifthpi)))
      cosw    = cos(w)
      sinw    = sin(w)
      lvt     = 1.45*log(real(nt))
      nn      = (nt+1)**2*10

c      allocate ( xn(nt+1,nt+1,10,3)   )
 
      do id=1,10
 
         sgn = 1.
         if(id .ge. 6) sgn = -1.
         phi = (2*mod(id, 5) - 3 + (id - 1)/5)*fifthpi
 
         xn(   1,   1,id,1) =  0.
         xn(   1,   1,id,2) =  0.
         xn(   1,   1,id,3) =  sgn
         xn(nt+1,   1,id,1) =  sinw*cos(phi)
         xn(nt+1,   1,id,2) =  sinw*sin(phi)
         xn(nt+1,   1,id,3) =  cosw*sgn
         xn(   1,nt+1,id,1) =  sinw*cos(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,2) =  sinw*sin(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,3) =  cosw*sgn
         xn(nt+1,nt+1,id,1) =  sinw*cos(phi + fifthpi)
         xn(nt+1,nt+1,id,2) =  sinw*sin(phi + fifthpi)
         xn(nt+1,nt+1,id,3) = -cosw*sgn
 
c         do k=0,int(lvt+.05)-1
         do k=0,lvt-1
 
            m  = 2**k
            l  = nt/m
            l2 = l/2
 
c...        rows of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j1-1)*l + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1,i2-l2,id,1),
     &                          xn(i1,i2+l2,id,1),nn)
               end do
            end do
 
c...        columns of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j2-1)*l + l2 + 1
                     i2 = (j1-1)*l + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2,id,1),
     &                          xn(i1+l2,i2,id,1),nn)
               end do
            end do
 
c...        diagonals of diamond--
 
            do j1=1,m
               do j2=1,m
                     i1 = (j1-1)*l + l2 + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2+l2,id,1),
     &                          xn(i1+l2,i2-l2,id,1),nn)
               end do
            end do
 
         end do
 
      end do

c      deallocate ( xn )

      end
