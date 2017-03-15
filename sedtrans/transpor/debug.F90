Program transpor_debug
    real(4)::HD,UG,HS,TP,PHI,D50,D90,WS,RC,RW,DS,A,RNU,RHOW
    real(4)::RHOS,BETA,SSR,SSF,SBF
    HD = 2.
    UG = 0.5
    HS = 0.0
    TP = 0.0
    PHI = 0.0
    D50 = 0.00025
    D90 = 0.00050
    WS = 0.010
    RC = 0.100
    RW = 0.100
    DS = 0.100
    A = 0.100
    RNU = 1.0e-6
    RHOW = 1025.
    RHOS = 2650.
    BETA = 1.0
    CALL TRANSP(HD,UG,HS,TP,PHI,D50,D90,WS,RC,RW,DS,A,RNU,RHOW,&
           RHOS,BETA,SSR,SSF,SBF)
    PRINT*,'SSR = ',SSR
    PRINT*,'SSF = ',SSF
    PRINT*,'SBF = ',SBF
    PAUSE
    
end


