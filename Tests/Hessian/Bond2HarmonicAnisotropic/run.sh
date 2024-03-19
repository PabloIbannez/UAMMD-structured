python scripts/generateTest.py

cd results

UAMMDlauncherDouble test.json 2>>stderr.log

cd ..

python ../scripts/analyzeResults.py 1e-5
