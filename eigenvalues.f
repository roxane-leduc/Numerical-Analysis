       subroutine eigenvalues(T,n)
               implicit none
               integer i,j,n,k
               double precision T(n,n)
               double precision, dimension(:,:), allocatable :: R,Q
               double precision, dimension(:), allocatable ::prec,act
               double precision seuil

               allocate(prec(n))
               allocate(act(n))
               allocate(R(n,n)) 
               allocate(Q(n,n))

               prec=0
               seuil=1e-6
               do i=1,500
                  call decompose_QR(T,n,Q,R)
                  T = matmul(R,Q)
                  do k=1,n
                        act(k)=T(k,k)
                  enddo
                  
                  if (norm2(prec-act)<seuil) then
                        write(*,*) "Nombres d iterations:",i
                        exit
                  endif
                  prec=act
              enddo

               write(*,*) "Valeurs propres:"
               do j=1,n
                write(*,*) T(j,j)
               enddo

               deallocate(prec)
               deallocate(act)
               deallocate(R) 
               deallocate(Q)

      end
