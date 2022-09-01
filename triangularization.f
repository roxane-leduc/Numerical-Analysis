      subroutine triangularization(A,n,T)
              implicit none
              double precision, dimension(:,:), allocatable :: Q,ID,H
              double precision,dimension(:), allocatable :: u
              double precision, dimension(:,:), allocatable :: Ucol
              double precision, dimension(:,:), allocatable :: Ulig
              integer c,n,i,j,k
              double precision alpha,r
              double precision T(n,n),A(n,n)

              allocate(Q(n,n))
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
              T = A

              do c=1,(n-1)
                u = T(:,c)

                do i=1,c
                        u(i) = 0.d0
                enddo

                alpha = 0
                if (u(c+1) < 0) then
                        alpha = norm2(u)
                else
                        alpha = -norm2(u)
                endif

                r=sqrt(0.5*alpha*(alpha-u(c+1)))
                u(c+1)=(u(c+1)-alpha)/(2*r)

                !u(c+2:n)=u(c+2:n)/(2*r)
                do k=(c+2),n
                      u(k)=u(k)/(2*r)
                enddo  

                ! Reformater le vect. col. en une matrice n lig. x 1 col
                Ucol = reshape(u, (/ n, 1 /))
                ! Reformater le vect. col. en une matrice 1 lig. x n col
                Ulig = reshape(u, (/ 1, n /))

                H = ID - 2 * matmul( Ucol, Ulig ) / (norm2(u)**2)

                Q = matmul(Q,H)
                T = matmul(matmul(transpose(H),T),H)

              enddo

              deallocate(Q)
              deallocate(ID)
              deallocate(H)
              deallocate(Ucol)
              deallocate(Ulig)
              deallocate(u)

      end

