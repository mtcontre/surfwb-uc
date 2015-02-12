MODULE senales
  integer:: GA1,GA2,GA3,GA4,Nsenal1,Nsenal2,Nsenal3,Nsenal4, IO1,IO2,IO3,IO4
  real (kind=8):: h01,h02,h03,h04
  real (kind=8),dimension(:,:),save,allocatable:: etaL1,qs1,us1,hs1, &
						  etaR2,qs2,us2,hs2,&
						  etaL3,qs3,us3,hs3, &
						  etaR4,qs4,us4,hs4, &
						  qA1,qA2,qA3,qA4, &
						  qsx1,qsy1, qsx2, qsy2, &
						  qsx3,qsy3, qsx4, qsy4, &
						  etaL9, hs9, &
						  etas1,etas2,etas3,etas4
						  
  real (kind=8),dimension(:),save,allocatable:: zA1,zA2,zA3,zA4,timeS1,timeS2,timeS3,timeS4,timeS9
  
 
END MODULE senales
