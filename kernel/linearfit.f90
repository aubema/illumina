! subroutine for linear regression
subroutine linearfit(x,y,nfit,a,b)
implicit none
integer,parameter::n=6
integer i,nfit
real*8, dimension(n)::x,y
real*8 s1,s2,s3,s4,a,b
s1=0.D0
s2=0.D0
s3=0.D0
s4=0.D0
do i=1,nfit
  s1=s1+x(i)
  s2=s2+x(i)**2.
  s3=s3+y(i)
  s4=s4+x(i)*y(i)
enddo
a=((nfit*s4)-(s1*s3))/((nfit*s2)-(s1**2)) ! slope
b=(s3-(s1*a))/nfit ! intercept
return
end
