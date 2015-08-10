SUBROUTINE PropSecNat(seccion,N, h,Area,Pmojado,Tsup,Rh,Dh)

! Funcion que calcula las propiedades geometricas de cualquier seccion
! [A Pm T Rh Dh]=PropSecNat(seccion, h)
! 
! "seccion" debe ser una matriz de coordenadas de dos columnas que contiene los puntos que definen la matriz
! la unica condicion es que los puntos se den en el orden de recorrido
! de la forma de la seccion
! 
! "h" es la altura de agua sobre el punto mas bajo de la seccion
real (kind=8),dimension(N,2)::seccion
real (kind=8),dimension(N)::Y,Z, hi, A, Pm, T
real (kind=8):: Z0,h,Zh, Yizq, Yder, Area, Pmojado,Tsup,Rh,Dh
integer,dimension(N)::Indi 
integer:: N, i

Y=seccion(:,1)
Z=seccion(:,2)

Z0=minval(Z)
Zh=Z0+h

T(:)=0.0D0
Pm(:)=0.0D0
Yizq=0.0D0
Yder=0.0D0

!Nodos secos y mojados

DO i=1,N

hi(i)=Zh-Z(i)
    
if (Z(i).le.Zh) then !Nodo Mojado
   
      Indi(i)=1
            
else !Nodo Seco
      Indi(i)=0
      hi(i)=0.0D0
end if

END DO

!vamos por cada par de nodos
do i=1,(N-1)
    if (Indi(i)==0) then    !%primer nodo seco
        if (Indi(i+1)==0) then  !%y segundo nodo seco (ambos fuera)
            A(i)=0.0D0
            Pm(i)=0.0D0
            T(i)=0.0D0
        else        !%y segundo nodo mojado (CC:BB izq)
            Yizq=Y(i)-(Y(i+1)-Y(i))/(hi(i+1)-hi(i))*hi(i)
            T(i)=Y(i+1)-Yizq
            A(i)=hi(i+1)*T(i)/2.0D0
            Pm(i)=sqrt(hi(i+1)**2.0D0+T(i)**2.0D0)
        end if
    else        !%primer nodo mojado
        if (Indi(i+1)==0) then  !%y segundo nodo seco (CC:BB der)
            Yder=Y(i)-(Y(i+1)-Y(i))/(hi(i+1)-hi(i))*hi(i)
            T(i)=Yder-Y(i)
            A(i)=hi(i)*T(i)/2.0D0
            Pm(i)=sqrt(hi(i)**2.0D0+T(i)**2.0D0)
        else        !%y segundo nodo mojado (nodo central)
            T(i)=Y(i+1)-Y(i)
            A(i)=(hi(i)+hi(i+1))/2.0D0*T(i)
            Pm(i)=sqrt(T(i)**2.0D0+(hi(i)-hi(i+1))**2.0D0)
        end if
    end if

end do

Area=0.0D0
Do i=1,N-1
Area=Area+A(i)
end do
!Area=sum(A)

Pmojado=sum(Pm)
if (Indi(1)==1) then
    Pmojado=Pmojado+hi(1)
end if
if (Indi(N)==1) then
    Pmojado=Pmojado+hi(N)
end if

Rh=Area/Pmojado
Tsup=sum(T)
Dh=Area/Tsup
! print*,Area
! pause
! print*,hi
! pause
! print*,T
! pause
! print*,'A=',A
! pause

END SUBROUTINE PropSecNat