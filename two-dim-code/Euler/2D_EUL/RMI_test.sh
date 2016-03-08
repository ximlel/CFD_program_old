#!/bin/bash


#=============================RMI KA test=========================
:<<KA
nohup ./EUL_source.out RMI/RMI_K1A10 RMI/RMI_K1A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A20 RMI/RMI_K1A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A40 RMI/RMI_K1A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A80 RMI/RMI_K1A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A10 RMI/RMI_K2A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A20 RMI/RMI_K2A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A40 RMI/RMI_K2A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A80 RMI/RMI_K2A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A10 RMI/RMI_K4A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A20 RMI/RMI_K4A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A40 RMI/RMI_K4A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A80 RMI/RMI_K4A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A10 RMI/RMI_K8A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A20 RMI/RMI_K8A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A40 RMI/RMI_K8A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A80 RMI/RMI_K8A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A10 RMI/RMI_K16A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A20 RMI/RMI_K16A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A40 RMI/RMI_K16A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A80 RMI/RMI_K16A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A80.out 2>&1 &
KA
#=============================RMI KA test=========================

#=============================RMI lnKA more test==========================
#:<<lnKA
nohup ./EUL_source.out RMI/RMI_lnKA_-3 RMI/RMI_lnKA_-3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.9 RMI/RMI_lnKA_-2.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.8 RMI/RMI_lnKA_-2.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.7 RMI/RMI_lnKA_-2.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.6 RMI/RMI_lnKA_-2.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.5 RMI/RMI_lnKA_-2.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.4 RMI/RMI_lnKA_-2.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.3 RMI/RMI_lnKA_-2.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.2 RMI/RMI_lnKA_-2.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.1 RMI/RMI_lnKA_-2.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2 RMI/RMI_lnKA_-2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.9 RMI/RMI_lnKA_-1.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.8 RMI/RMI_lnKA_-1.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.7 RMI/RMI_lnKA_-1.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.6 RMI/RMI_lnKA_-1.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.5 RMI/RMI_lnKA_-1.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.4 RMI/RMI_lnKA_-1.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.3 RMI/RMI_lnKA_-1.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.2 RMI/RMI_lnKA_-1.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.1 RMI/RMI_lnKA_-1.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1 RMI/RMI_lnKA_-1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.9 RMI/RMI_lnKA_-0.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.8 RMI/RMI_lnKA_-0.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.7 RMI/RMI_lnKA_-0.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.6 RMI/RMI_lnKA_-0.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.5 RMI/RMI_lnKA_-0.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.4 RMI/RMI_lnKA_-0.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.3 RMI/RMI_lnKA_-0.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.2 RMI/RMI_lnKA_-0.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.1 RMI/RMI_lnKA_-0.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0 RMI/RMI_lnKA_0 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.1 RMI/RMI_lnKA_0.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.2 RMI/RMI_lnKA_0.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.3 RMI/RMI_lnKA_0.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.4 RMI/RMI_lnKA_0.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.5 RMI/RMI_lnKA_0.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.6 RMI/RMI_lnKA_0.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.7 RMI/RMI_lnKA_0.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.8 RMI/RMI_lnKA_0.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.9 RMI/RMI_lnKA_0.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1 RMI/RMI_lnKA_1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.1 RMI/RMI_lnKA_1.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.2 RMI/RMI_lnKA_1.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.3 RMI/RMI_lnKA_1.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.4 RMI/RMI_lnKA_1.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.5 RMI/RMI_lnKA_1.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.6 RMI/RMI_lnKA_1.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.7 RMI/RMI_lnKA_1.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.8 RMI/RMI_lnKA_1.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.9 RMI/RMI_lnKA_1.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2 RMI/RMI_lnKA_2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.1 RMI/RMI_lnKA_2.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.2 RMI/RMI_lnKA_2.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.3 RMI/RMI_lnKA_2.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.4 RMI/RMI_lnKA_2.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.5 RMI/RMI_lnKA_2.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.6 RMI/RMI_lnKA_2.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.7 RMI/RMI_lnKA_2.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.8 RMI/RMI_lnKA_2.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.9 RMI/RMI_lnKA_2.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_3 RMI/RMI_lnKA_3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_3.out 2>&1 &
#lnKA
#=============================RMI lnKA more test==========================

#=============================RMI At more test==========================
:<<At
nohup ./EUL_source.out RMI_At/RMI_At_1 RMI_At/RMI_At_1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_1.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_2 RMI_At/RMI_At_2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_2.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_3 RMI_At/RMI_At_3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_3.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_4 RMI_At/RMI_At_4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_4.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_5 RMI_At/RMI_At_5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_5.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_6 RMI_At/RMI_At_6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_6.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_7 RMI_At/RMI_At_7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_7.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_8 RMI_At/RMI_At_8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_8.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_9 RMI_At/RMI_At_9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_9.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_10 RMI_At/RMI_At_10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_10.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_11 RMI_At/RMI_At_11 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_11.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_12 RMI_At/RMI_At_12 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_12.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_13 RMI_At/RMI_At_13 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_13.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_14 RMI_At/RMI_At_14 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_14.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_15 RMI_At/RMI_At_15 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_15.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_16 RMI_At/RMI_At_16 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_16.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_17 RMI_At/RMI_At_17 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_17.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_18 RMI_At/RMI_At_18 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_18.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_19 RMI_At/RMI_At_19 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_19.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_20 RMI_At/RMI_At_20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_20.out 2>&1 &
At
#=============================RMI At more test==========================
