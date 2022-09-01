      subroutine decompose_QR(T,n,Q,R)
              implicit none

              double precision, dimension(:,:), allocatable :: H,ID
              double precision,dimension(:), allocatable :: u
              double precision, dimension(:,:), allocatable :: Ucol
              double precision, dimension(:,:), allocatable :: Ulig
              integer c,i,j,alpha,n
              double precision T(n,n),Q(n,n),R(n,n)
              
              allocate(H(n,n))
              allocate(ID(n,n))
              allocate(Ucol(n,1))
              allocate(Ulig(1,n))
              allocate(u(n))


              ID = 0
              do i=1,n
                ID(i,i) = 1
              enddo

              Q = ID
              R = T

              do c=1,(n-1)
                u = R(:,c)

                alpha = 0
                if (u(1) < 0) then
                        alpha = -1
                else
                        alpha = 1
                endif

                do i=1,c-1
                        u(i) = 0.d0
                enddo

                u(c) = u(c) - alpha * norm2(u)

                ! Reformater le vect. col. en une matrice n lig. x 1 col
                Ucol = reshape(u, (/ n, 1 /)) 
                ! Reformater le vect. col. en une matrice 1 lig. x n col
                Ulig = reshape(u, (/ 1, n /))
                
                H = ID - 2 * matmul( Ucol, Ulig ) / (norm2(u)**2)
                
                Q = matmul(Q,H)
                R = matmul(H,R)

              enddo

              R = matmul(transpose(Q),T)

              deallocate(ID)
              deallocate(u)
              deallocate(H)
              deallocate(Ucol)
              deallocate(Ulig)

      end
