
module AMSU_constants

    real(8),parameter        :: A_SMALL_NUMBER = 1.0d-10
    real(8),parameter        :: TWO_PI = 6.283185307d0
    real(8),parameter        :: PI = 3.141592654d0

    integer(4),parameter            :: ILU                = 20
    integer(4),parameter            :: AMSU_SAA_LU        = 21
    integer(4),parameter            :: log_unit            = 24
    integer(4),parameter            :: gap_unit            = 23
    integer(4),parameter            :: AMSU_HEADER_LIST_UNIT = 25


    integer(4),parameter            :: NUM_SATS = 9+5+1+3

    integer(4),parameter            :: NOAA_15 = 10
    integer(4),parameter            :: NOAA_16 = 11
    integer(4),parameter            :: NOAA_17 = 12
    integer(4),parameter            :: NOAA_18 = 13
    integer(4),parameter            :: METOP_A = 14
    integer(4),parameter            :: AQUA    = 15
    integer(4),parameter            :: NOAA_19 = 16
    integer(4),parameter            :: METOP_B = 17
    integer(4),parameter            :: METOP_C = 18



    real(4),parameter,dimension(15) :: AMSU_A_Freq = &
                (/23.800,    & ! Channel 1
                  31.400,   & ! Channel 2
                  50.300,    & ! Channel 3
                  52.800,    & ! Channel 4
                  53.596,    & ! Channel 5
                  54.400,    & ! Channel 6
                  54.940,    & ! Channel 7
                  55.500,   & ! Channel 8
                  57.290334,& ! Channel 9
                  57.290334,& ! Channel 10
                  57.290334,& ! Channel 11
                  57.290334,& ! Channel 12
                  57.290334,& ! Channel 13
                  57.290334,& ! Channel 14
                  89.000/)      ! Channel 15

    real(4),parameter,dimension(15) :: AMSU_A_Stopband = &
                (/0.018,    & ! Channel 1
                  0.018,    & ! Channel 2
                  0.018,    & ! Channel 3
                  0.018,    & ! Channel 4
                  0.0,        & ! Channel 5
                  0.018,    & ! Channel 6
                  0.018,    & ! Channel 7
                  0.018,    & ! Channel 8
                  0.018,    & ! Channel 9
                  0.0,        & ! Channel 10
                  0.0,        & ! Channel 11
                  0.0,        & ! Channel 12
                  0.0,        & ! Channel 13
                  0.0,        & ! Channel 14
                  0.0/)      ! Channel 15

    real(4),parameter,dimension(15) :: AMSU_A_Freq_Split_1 = &
                (/0.0,        & ! Channel 1
                  0.0,        & ! Channel 2
                  0.0,        & ! Channel 3
                  0.0,        & ! Channel 4
                  0.115,    & ! Channel 5
                  0.0,        & ! Channel 6
                  0.0,        & ! Channel 7
                  0.0,        & ! Channel 8
                  0.0,        & ! Channel 9
                  0.217,    & ! Channel 10
                  0.3222,    & ! Channel 11
                  0.3222,    & ! Channel 12
                  0.3222,    & ! Channel 13
                  0.3222,    & ! Channel 14
                  0.00/)      ! Channel 15

    real(4),parameter,dimension(15) :: AMSU_A_Freq_Split_2 = &
                (/0.0,        & ! Channel 1
                  0.0,        & ! Channel 2
                  0.0,        & ! Channel 3
                  0.0,        & ! Channel 4
                  0.0,        & ! Channel 5
                  0.0,        & ! Channel 6
                  0.0,        & ! Channel 7
                  0.0,        & ! Channel 8
                  0.0,        & ! Channel 9
                  0.0,        & ! Channel 10
                  0.048,    & ! Channel 11
                  0.022,    & ! Channel 12
                  0.010,    & ! Channel 13
                  0.0045,    & ! Channel 14
                  0.00/)      ! Channel 15

    real(4),parameter,dimension(15) :: AMSU_A_BANDWIDTH = &
                (/0.251,        & ! Channel 1
                  0.161,        & ! Channel 2
                  0.161,        & ! Channel 3
                  0.3805,        & ! Channel 4
                  0.170,        & ! Channel 5
                  0.3805,        & ! Channel 6
                  0.3805,        & ! Channel 7
                  0.3103,        & ! Channel 8
                  0.3300,        & ! Channel 9
                  0.07658,        & ! Channel 10
                  0.03511,    & ! Channel 11
                  0.01529,    & ! Channel 12
                  0.00793,    & ! Channel 13
                  0.00294,    & ! Channel 14
                  1.9989/)      ! Channel 15

    integer(4),parameter                        :: NUM_AMSU_FOVS=30

    real(4),dimension(15),parameter :: AMSU_VIEW_ANGLES = (/1.6666666,5.0000000,8.3333333,11.666667,15.000000,    &
                                                            18.333333,21.666667,25.000000,28.333333,31.666667,    &
                                                            35.000000,38.333333,41.666667,45.000000,48.333333/)

    real(4),dimension(15),parameter :: AMSU_NOM_EIAS    = (/1.875947, 5.629541, 9.388301,13.155880,16.936250,        &
                                                           20.733890,24.554020,28.402830,32.287970,36.219100,        &
                                                           40.208800,44.274040,48.438570,52.737200,57.224260/)


    integer(4),parameter,dimension(15) :: AMSU_A_Polarization = & ! 1 = V, 2 = H
                (/1,1,1,1,2,2,1,2,2,2,2,2,2,2,1/)



    real(8),parameter                :: AMSU_SCAN_TIME    = 9.25925926d-5  ! in days: typical scan time 8 s
    integer(4),parameter            :: AMSU_MAX_SCANS_PER_ORBIT = 790
    integer(4),parameter            :: FILE_LEN = 120

    character(7),dimension(NUM_SATS),parameter :: amsu_sat_names_2 =(/'TIROS-N','NOAA-06', 'NOAA-07', 'NOAA-08','NOAA-09', &
                                                                 'NOAA-10','NOAA-11','NOAA-12','NOAA-14','NOAA-15','NOAA-16', &
                                                                 'NOAA-17','NOAA-18','METOP-A','AQUA   ','NOAA-19','METOP-B','METOP-C'/)

    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_PENTAD2 =        (/59, 110,255, 389, 509, 651, 785, 977, 1243,1504,1685,1800,2000, 2140, 1600,2240,2480,2480/)  !first pentad should be
    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_PENTAD =         (/59, 110,255, 389, 509, 651, 785, 977, 1243,1504,1685,1300,1685, 1504, 1800,1800,2480,2480/)
    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_DAY_NUM_1978   = (/294,546,1272,1941,2541,3251,3923,4885,6211,7520,8423,6500,8423, 6500, 8000,8000,12400,12400/)
    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_DAY_NUM_1978_2 = (/294,546,1272,1941,2541,3251,3923,4885,6211,7520,8423,9000,10000,10700,8000,11200,12400,12400/)  !fisrt day should be
                                                                            !  TN   N6  N7   N8   N9   N10  N11  N12  N14  N15  N16   N17  N18  MO-A  AQUA  N19
    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_MONTH          = (/0,  12,  36,  60,  72,  96,  120, 156, 204, 240, 276,  204, 324, 204,  276, 276, 408, 408/)
    integer(4),parameter,dimension(NUM_SATS)  :: AMSU_FIRST_MONTH2         = (/0,  12,  36,  60,  72,  96,  120, 156, 204, 240, 276,  204, 324, 348,  276, 372, 408, 408/)  !first month should be....

    integer(4),parameter                :: MAX_ORBITS_TO_PROCESS = 120000


end module AMSU_constants