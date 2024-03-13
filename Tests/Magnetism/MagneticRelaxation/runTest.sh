#Generate the JSON files
python scripts/generateTestRelaxationFixed.py
python scripts/generateTestRelaxationRigidDipole.py
python scripts/generateTestRelaxationCombined.py

cd results

#Run the simulations

UAMMDlauncher simulationFixedParticles.json
UAMMDlauncher simulationRigidDipole.json
UAMMDlauncher simulationCombined.json

cd ..

#Analyze the results
python scripts/analysis.py
