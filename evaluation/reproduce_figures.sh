python3 spe11a_visualize_spatial_maps.py -g TetraTech OPM4 Pau-Inria IFPEN -f /media/bernd/bernd/spe11/data -t 120
python3 spe11a_assemble_time_series.py -g CSIRO1 CSIRO2 GEOS1 GEOS2 IFPEN OPM1 OPM2 OPM3 OPM4 TetraTech UT-CSEE1 UT-CSEE2 UT-CSEE3 -f /media/bernd/bernd/spe11/data -cAB Calgary CTC-CNE OpenGoSim Pau-Inria PFLOTRAN -cC Calgary CAU-Kiel CTC-CNE PFLOTRAN Pau-Inria SLB1 SLB2 -t /media/bernd/bernd/spe11/evaluation/spe11a/sparse
python3 spe11a_assemble_performance_time_series.py -g CAU-Kiel CSIRO1 CSIRO2 CTC-CNE IFPEN OPM1 OPM2 OPM3 OPM4 Pau-Inria UT-CSEE1 UT-CSEE2 UT-CSEE3 -f /media/bernd/bernd/spe11/data

python3 spe11b_visualize_spatial_maps.py -g IFPEN1 OPM4 SINTEF2 Rice2 -f /media/bernd/bernd/spe11/data -t 50
python3 spe11b_assemble_time_series.py -g CSIRO DARTS GEOS1 GEOS2 IFPEN1 IFPEN2 KFUPM OPM1 OPM2 OPM3 OPM4 Pau-Inria Rice1 Rice2 TetraTech1 TetraTech2 UT-CSEE1 UT-CSEE2 -f /media/bernd/bernd/spe11/data -cAB Calgary CTC-CNE OpenGoSim1 OpenGoSim2 OpenGoSim3 PFLOTRAN -cC Calgary CAU-Kiel CTC-CNE PFLOTRAN SINTEF1 SINTEF2 SINTEF3 SINTEF4 SLB Stuttgart1 Stuttgart2 Stuttgart3 Stuttgart4 -t /media/bernd/bernd/spe11/evaluation/spe11b/sparse
python3 spe11b_assemble_performance_time_series.py -g CAU-Kiel CSIRO CTC-CNE IFPEN1 IFPEN2 KFUPM OPM1 OPM2 OPM3 OPM4 Pau-Inria TetraTech2 UT-CSEE1 UT-CSEE2 -d IFPEN1 IFPEN2 SINTEF1 SINTEF2 SINTEF3 SINTEF4 -f /media/bernd/bernd/spe11/data

python3 spe11c_visualize_spatial_maps.py -g SINTEF3 OPM4 GEOS2 OpenGoSim1 -f /media/bernd/bernd/spe11/data -t 1000
python3 spe11c_assemble_time_series.py -g CSIRO GEOS1 GEOS2 IFPEN OPM1 OPM2 OPM3 OPM4 TetraTech1 TetraTech2 UT-CSEE -f /media/bernd/bernd/spe11/data -cAB Calgary CTC-CNE OpenGoSim1 OpenGoSim2 PFLOTRAN -cC Calgary CAU-Kiel CTC-CNE Pau-Inria PFLOTRAN SINTEF1 SINTEF2 SINTEF3 SLB -t /media/bernd/bernd/spe11/evaluation/spe11c/sparse
python3 spe11c_assemble_performance_time_series.py -g CAU-Kiel CSIRO CTC-CNE IFPEN OPM1 OPM2 OPM3 OPM4 Pau-Inria TetraTech2 UT-CSEE -d SINTEF1 SINTEF2 SINTEF3 -f /media/bernd/bernd/spe11/data
