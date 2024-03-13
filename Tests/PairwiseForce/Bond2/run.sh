python scripts/generateTest.py

cd results

UAMMDlauncher test.json 2>>stderr.log

cd ..

python scripts/analyze.py 1e-5
