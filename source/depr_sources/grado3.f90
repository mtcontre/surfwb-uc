SUBROUTINE grado3(a,b,c,d,x)
!Subrutina que resuelve una ec. de 3 orden
implicit none
real (kind=8):: a,b,c,d,x,a1,a2,a3,p,Q,R,De,S,T, theta,x1,x2,x3,pi,y1,y2,y3, cero
a1=b/a
a2=c/a
a3=d/a
pi=3.141592654D0

! p=(3.0D0*a2-a1**2.0D0)/3.0D0
! !q=(9.0D0*a1*a2-27.0D0*a3-2.0D0*a1**3.0D0)/54.0D0
! q=(2.0D0*a1**3.0D0-9.0D0*a1*a2-27.0D0*a3 )/27.0D0
! De=(p/3.0D0)**3.0D0+(q/2.0D0)**2.0D0
! 
! if (De>=0.0D0) then
! S=(-q/2.0D0+sqrt(De))**(1.0D0/3.0D0)
! T=(-q/2.0D0-sqrt(De))**(1.0D0/3.0D0)
! y1=S+T
! x=y1-a1/3.0D0
! else
! theta=acos(-q/2.0D0/sqrt(abs(p)**(3.0D0)/27.0D0))
! 
! y1=2.0D0*sqrt(abs(p)/3.0D0)*cos(1.0D0/3.0D0*theta)
! y2=-2.0D0*sqrt(abs(p)/3.0D0)*cos(1.0D0/3.0D0*(theta+pi))
! y3=-2.0D0*sqrt(abs(p)/3.0D0)*cos(1.0D0/3.0D0*(theta-pi))
! 
! x1=y1-a1/3.0D0
! x2=y2-a1/3.0D0
! x3=y3-a1/3.0D0
! x=max(x1,x2,x3)
! 
! end if

Q=(3.0D0*a2-a1**2.0D0)/9.0D0
R=(9.0D0*a1*a2-27.0D0*a3-2.0D0*a1**3.0D0)/54.0D0
De=Q**3.0D0+R**2.0D0
if (De>0.0D0) then
S=(R+sqrt(De))**(1.0D0/3.0D0)
T=(R-sqrt(De))**(1.0D0/3.0D0)
x=S+T-1.0D0/3.0D0*a1

else if (De==0.0D0) then
x1=S+T-1.0D0/3.0D0*a1
x2=-0.5D0*(S+T)-1.0D0/3.0D0*a1
x=max(x1,x2)

else
theta=acos(-R/sqrt(-1.0D0*Q**3.0D0))
x1=2.0D0*sqrt(-Q)*cos(1.0D0/3.0D0*theta)
x2=2.0D0*sqrt(-Q)*cos(1.0D0/3.0D0*theta+2.0D0/3.0D0*pi)
x3=2.0D0*sqrt(-Q)*cos(1.0D0/3.0D0*theta+4.0D0/3.0D0*pi)
x=max(x1,x2,x3)

!Verificación

end if

 cero=a*x**3.0D0+b*x**2.0D0+c*x+d
 !print*,cero
 
! x = -b/(3.0D0*a) - (2.0D0**(1.0D0/3.0D0)*(-b**2.0D0 + 3.0D0*a*c))/(3.0D0*a*(-2.0D0*b**3.0D0 + 9.0D0*a*b*c - 27.0D0*a**2.0D0*d + Sqrt(4.0D0*(-b**2.0D0 + 3*a*c)**3.0D0 + (-2.0D0*b**3.0D0 + 9.0D0*a*b*c - 27.0D0*a**2.0D0*d)**2.0D0))**(1.0D0/3.0D0)) + (-2.0D0*b**3.0D0 + 9.0D0*a*b*c - 27.0D0*a**2.0D0*d + Sqrt(4.0D0*(-b**2.0D0 + 3.0D0*a*c)**3.0D0 + (-2.0D0*b**3.0D0 + 9.0D0*a*b*c - 27.0D0*a**2.0D0*d)**2.0D0))**(1.0D0/3.0D0)/(3.0D0*2.0D0**(1.0D0/3.0D0)*a)

!print*,'x=',Q,R,De,S,T,x



END SUBROUTINE grado3

! Solution: To find roots to the following cubic equation where a, b, c, and d are real.
! 
!    a*x3 + b*x2 + c*x + d = 0
! 
! Formula:
! 
!   Step 1: Calculate p and q
!           p = ( 3*c/a - (b/a)2 ) / 3
!           q = ( 2*(b/a)3 - 9*b*c/a/a + 27*d/a ) / 27
!   Step 2: Calculate discriminant D
!           D = (p/3)3 + (q/2)2
!   Step 3: Depending on the sign of D, you follow different strategy.
!           If D<0, three distinct real roots.
!           If D=0, three real roots of which at least two are equal.
!           If D>0, one real and two complex roots.
!   Step 3a: For D>0 and D=0
!           Calculate u and v
!           u = cubic_root(-q/2 + sqrt(D))
!           v = cubic_root(-q/2 - sqrt(D))
!           Find the three transformed roots
!           y1 = u + v
!           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
!           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
!   Step 3b: Alternately, for D<0, a trigonometric formulation is more convenient
!           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
!           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
!           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
!           where phi = acos(-q/2/sqrt(|p|3/27))
!                 pi  = 3.141592654...
!   Step 4  Finally, find the three roots
!           x = y - b/a/3
! 
!   Things to watch out for:
!     1. Make sure you know what is integer, real, and complex.
!     2. FORTRAN's SQRT (square root) function takes only non-negative arguments.
!     3. FORTRAN's exponentiation is "**", not "^".
!     4. There is no "cubic_root" function in FORTRAN; you do it by raising a
!        number to a factional power, e.g., 0.333... for cubic root.
!     5. Your can only raise a non-negative number to a fractional power.
!        e.g., (-2.)**(1./3.) will not work; you need to issue -(2.)**(1./3.)
!        And do not forget the decimal point.  For example, in computing 2.**(1/3),
!        FORTRAN evaluates 1/3 first, which is 0, not 0.33... as you might expect.
!     6  There is no "|" in FORTRAN.  The absolute value function is "ABS".
!     7. FORTRAN is case-insensitive; do not simply use "d" (for coefficient)
!        and "D" (for discriminant), as FORTRAN thinks both variables are equivalent.
