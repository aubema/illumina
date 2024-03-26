! subroutine for linear regression
subroutine linearfit(x,y,a,b)
implicit none
integer,parameter::n=4
integer i
real*8, dimension(n)::x,y
real*8::s1,s2,s3,s4,a,b
do i=1,n
  s1=s1+x(i)
  s2=s2+x(i)**2
  s3=s3+y(i)
  s4=s4+x(i)*y(i)
enddo
a=((n*s4)-(s1*s3))/((n*s2)-(s1**2)) ! slope
b=(s3-(s1*a))/n ! intercept
return
end
