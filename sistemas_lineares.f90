!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PROGRAMA SISTEMAS LINEARES - VERSÃO 2008							 !
!	1 - MÉTODO ITERATIVO DE JACOBI								 !
!	2 - MÉTODO ITERATIVO DE GAUSS-SEIDEL							 !
!	3 - MÉTODO ITERATIVO DE GAUSS-SEIDEL COM RELAXAÇÃO					 !
!	4 - TDMA (PARA MATRIZES TRIDIAGONAIS)							 !
!												 !
!	PROGRAMADOR: RAFAEL GABLER GONTIJO							 !	
!												 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program sistemas_lineares
implicit none
integer i,j,n,k,metodo, entrada, saida
real, dimension(:,:),allocatable::A
real, dimension(:),allocatable::x,b
real w
character*16 arquivo 

print *,'******************************************************************************'
print *,'*                  Programa Sistemas Lineares - v.1.0 - 2009                 *'
print *,'*                                                                            *'
print *,'* Autor: Rafael Gabler Gontijo                                               *'
print *,'*                                                                            *'
print *,'* Programa utilizado para resolucao de sistemas lineares                     *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''
print *,'  Primeiramente escolha o metodo utilizado para resolucao do sistema linear'
print *,''
print *,'1 - Metodo iterativo de Jacobi'
print *,'2 - Metodo iterativo de Gauss-Seidel'
print *,'3 - Metodo iterativo de Gauss-Seidel com relaxacao'
print *,'4 - TDMA (para matrizes tridiagonais)'


read *,metodo

print *,'Digite a ordem do sistema linear'
read *,n

allocate(A(n,n))
allocate(x(n))
allocate(b(n))

print *,'Escolha a forma de entrada dos dados'
print *,'Tela (1) ou arquivo (2)'
read *,entrada
print *,''
print *,'Digite agora a forma de saida dos dados'
print *,'Tela (1) ou arquivo (2)'
read*,saida

if(entrada.eq.1) then
print *,'Digite cada termo da matriz A'
do i=1,n
do j=1,n
print *,'A',i,j
read*,A(i,j)
end do 
end do
print *,''
print *,'Digite agora cada termo do vetor b'
do i=1,n
print *,'b',i
read*,b(i)
end do
end if

if(entrada.eq.2) then
print *,'Digite o nome do arquivo de entrada:'
read(*,'(A)') arquivo
open(1,file=arquivo)

read(1,'(F9.3)') A
read(1,'(F9.3)') b

do i=1,n
do j=1,n
A(i,j)=A(j,i)
end do 
end do

end if


! Basicamente o programa direciona a solucao pelo metodo escolhido pelo usuario chamando as subro-
! tinas definidas posteriormente

if (metodo.eq.1) then
print *,'Digite o numero de iteracoes:'
read *,k
call jacobi(A,x,b,n,k,saida)
pause
end if


if (metodo.eq.2) then
print *,'Digite o numero de iteracoes:'
read *,k
call seidel(A,x,b,n,k,saida)
pause
end if

if (metodo.eq.3) then
print *,'Digite o numero de iteracoes:'
read *,k
print *,'Digite o valor de w:'
print *,'Se 0<w<1 (sub-relaxado) se w>1 super relaxado (SOR)'
read *,w
call seidel_r(A,x,b,n,k,w,saida)
pause
end if

if (metodo.eq.4) then
call tdma(A,x,b,n,saida)
pause
end if


end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TDMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tdma(A,x,b,n,saida)

implicit none
integer i,j,n,saida
real A(n,n),C(n,n+1),L(n,n),U(n,n)
real x(n),b(n),z(n)


do i=1,n
do j=1,n
 C(i,j)=0.0
end do
end do

do i=1,n
do j=1,n
 C(i,j)=A(i,j)
end do
end do


do i=1,n
 C(i,n+1)=b(i)
end do



L(1,1)=C(1,1)
U(1,2)=C(1,2)/L(1,1)

do i=2,n-1
L(i,i-1)=C(i,i-1)
L(i,i)=C(i,i)-L(i,i-1)*U(i-1,i)
U(i,i+1)=C(i,i+1)/L(i,i)
end do


L(n,n-1)=C(n,n-1)
L(n,n)=C(n,n)-(L(n,n-1)*U(n-1,n))

z(1)=C(1,n+1)/L(1,1)

do i=2,n
z(i)=(C(i,n+1)-L(i,i-1)*z(i-1))/L(i,i)
end do

x(n)=z(n)



do i=1,n-1
x(n-i)=z(n-i)-(U(n-i,n-i+1)*x(n-i+1))
end do

if(saida.eq.1)then
write(*,*) 'A solucao e:'
do i=1,n
write(*,*) x(i)
end do
end if
pause

if(saida.eq.2)then
open (2,file='saida.dat')
write(2,*) 'A solucao e:'
do i=1,n
write(2,'(F8.4)') x(i)
end do
end if

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GAUSS-SEIDEL COM RELAXAÇÃO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine seidel_r(A,x,b,n,k,w,saida)

implicit none
integer i,j,n,iter,k,saida
real A(n,n), C(n,n), h(n), l(n)
real x(n),b(n),f(n),g,p,d(n)
real w

iter=0

do i=1,n
x(i)=b(i)/A(i,i)
end do

100 do i=1,n-2
	do j=i+1,n
	f(i)=(-A(i,j)*x(j))
	end do
	end do


	g=sum(f)

	do i=1,n
	do j=1,i-1
	d(i)=(-A(i,j)*x(j))
	p=sum(d)
	x(i)=(1-w)*x(i)+w*(b(i)+g+p)/A(i,i)
	end do
	end do

	iter=iter+1
	
	do i=1,n
	f(i)=0
	end do

	g=0

	if(iter.ne.k) then
	go to 100
	end if

	iter=0

	do i=1,n
	h(i)=b(i)/A(i,i)
	C(i,i)=0
	end do

	do i=1,n
	do j=1,n
	if(i.ne.j) then
	C(i,j)=-A(i,j)/A(i,i)
	end if
	end do
	end do

200	l=matmul(C,x)

	do i=1,n
	x(i)=l(i)+h(i)
	end do

	iter=iter+1

	if(iter.ne.k) then
	go to 200
	end if

if(saida.eq.1)then
write(*,*) 'A solucao pelo metodo de Gauss-Seidel com relaxacao e:'
do i=1,n
write(*,*) x(i)
end do
end if

if(saida.eq.2)then
open (2,file='saida.dat')
write(2,*) 'A solucao e:'
do i=1,n
write(2,'(F8.4)') x(i)
end do
end if


end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GAUSS-SEIDEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine seidel(A,x,b,n,k,saida)

implicit none
integer i,j,n,iter,k,saida
real A(n,n), C(n,n), D(n,n)
real x(n),b(n),f(n),g(n),l(n)

iter=0

do i=1,n
x(i)=b(i)/a(i,i)
g(i)=0
f(i)=0
end do

500 do i=2,n
    do j=1,i-1
       C(i,j)=A(i,j)*x(j)
    end do
    end do

f=sum(C,dim=2)


do i=1,n
x(i)=(b(i)-f(i)-g(i))/A(i,i)
end do


do i=1,n
do j=1,n
f(i)=0
g(i)=0
 C(i,j)=0
 D(i,j)=0
end do
end do


do i=1,n-1
do j=i+1,n
       D(i,j)=A(i,j)*x(j)
end do
end do

g=sum(D,dim=2)



iter=iter+1

if(iter.ne.k) then
go to 500
end if

l=matmul(x,A)

if(saida.eq.1)then
write(*,*) 'A solucao pelo metodo de Gauss-Seidel e:'
do i=1,n
write(*,*) x(i)
end do
end if

if(saida.eq.2)then
open (2,file='saida.dat')
write(2,*) 'A solucao e:'
do i=1,n
write(2,'(F8.4)') x(i)
end do
end if

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JACOBI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine jacobi(A,x,b,n,k,saida)

implicit none
integer i,j,n,iter,k,saida
real A(n,n), C(n,n), h(n), l(n)
real x(n),b(n),f(n),g

iter=0

do i=1,n
! x(i)=b(i)/A(i,i)
x(i)=1.0
end do

100 do j=1,n
	do i=1,n
	if(j.ne.i) then
	f(i)=(-A(i,j)*x(j))
	end if
	end do
	end do

	g=sum(f)

	do i=1,n
	x(i)=(b(i)+g)/A(i,i)
	end do

	iter=iter+1
	
	do i=1,n
	f(i)=0
	end do

	g=0

	if(iter.ne.k) then
	go to 100
	end if

	iter=0

	do i=1,n
	h(i)=b(i)/A(i,i)
	C(i,i)=0
	end do

	do i=1,n
	do j=1,n
	if(i.ne.j) then
	C(i,j)=-A(i,j)/A(i,i)
	end if
	end do
	end do

200	l=matmul(C,x)

	do i=1,n
	x(i)=l(i)+h(i)
	end do

	iter=iter+1

	if(iter.ne.k) then
	go to 200
	end if

if(saida.eq.1)then
write(*,*) 'A solucao pelo metodo de Jacobi e:'
do i=1,n
write(*,*) x(i)
end do
end if

if(saida.eq.2)then
open (2,file='saida.dat')
write(2,*) 'A solucao e:'
do i=1,n
write(2,'(F8.4)') x(i)
end do
end if

end
