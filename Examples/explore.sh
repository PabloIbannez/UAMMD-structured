N=10000
filename="scaling"$N".dat"
max=30
rm -f $filename
for i in `seq 1 $max`
do
        echo "Creating simulation $i/$max"
        python LennardJonesBatches.py $i $N
        cd simulation
        echo "Running simulation $i/$max"
        UAMMDlauncher simulation.json 2> $i.log
        cd ..
        (echo -n $i"  ";grep "FPS" simulation/$i.log | awk '{print $5}' | tr '\n' ' '; grep "elapsed" simulation/$i.log | awk '{print $6}' | sed 's/.$//') | cat >> $filename
done
