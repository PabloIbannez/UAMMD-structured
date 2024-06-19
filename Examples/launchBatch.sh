N=5
for i in `seq 1 $N`
do
        UAMMDlauncher simulation/simulation.json 2> mps_$i.log &
done
